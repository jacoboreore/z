import os
import json
import math
import time
import uuid
from typing import List, Dict, Any, Tuple

import boto3


# ----------------------------
# Config from environment
# ----------------------------
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

SQS_QUEUE_URL = os.environ["SQS_QUEUE_URL"]  # required
DDB_TABLE = os.environ["DDB_TABLE"]          # required (job metadata + chunk status)
JOB_PREFIX = os.environ.get("JOB_PREFIX", "valley")

# Chunking controls
MAX_CHUNKS = int(os.environ.get("MAX_CHUNKS", "1024"))  # total messages to create (upper bound)
TARGET_SECONDS_PER_CHUNK = int(os.environ.get("TARGET_SECONDS_PER_CHUNK", "600"))  # ~10min
MIN_CHUNKS = int(os.environ.get("MIN_CHUNKS", "64"))
MAX_WORKERS_HINT = int(os.environ.get("MAX_WORKERS_HINT", "256"))

# Precision defaults (workers can override if they want)
COARSE_PREC_BITS = int(os.environ.get("COARSE_PREC_BITS", "128"))

# N(t) strategy
N_STRATEGY = os.environ.get("N_STRATEGY", "sqrt_t_over_2pi")  # placeholder
N_MULTIPLIER = float(os.environ.get("N_MULTIPLIER", "1.0"))   # to tune empirically


sqs = boto3.client("sqs", region_name=AWS_REGION)
ddb = boto3.resource("dynamodb", region_name=AWS_REGION)
table = ddb.Table(DDB_TABLE)


# ----------------------------
# Helpers
# ----------------------------
def compute_n_max(t: float) -> int:
    """
    IMPORTANT:
    Replace this with your real Hardy+Gabcke cutoff N(t).
    This placeholder uses N ~ sqrt(t / (2π)) * multiplier.
    """
    if t <= 0:
        raise ValueError("t must be positive")

    if N_STRATEGY == "sqrt_t_over_2pi":
        n = math.sqrt(t / (2.0 * math.pi))
    else:
        # Default fallback: same as above
        n = math.sqrt(t / (2.0 * math.pi))

    n *= N_MULTIPLIER
    n_int = max(1, int(n))
    return n_int


def choose_num_chunks(n_max: int, max_workers_hint: int) -> int:
    """
    Choose number of chunks.
    We want:
      - enough chunks to keep heterogeneous fleet busy
      - avoid tiny chunks overhead
    Strategy:
      - start with >= max_workers_hint * 4
      - clamp between MIN_CHUNKS and MAX_CHUNKS
      - never exceed n_max
    """
    desired = max(MIN_CHUNKS, max_workers_hint * 4)
    desired = min(desired, MAX_CHUNKS)
    return min(desired, n_max)


def split_range_balanced(n_max: int, num_chunks: int) -> List[Tuple[int, int]]:
    """
    Split [1..n_max] into num_chunks contiguous balanced segments.
    """
    base = n_max // num_chunks
    rem = n_max % num_chunks

    ranges = []
    cur = 1
    for i in range(num_chunks):
        size = base + (1 if i < rem else 0)
        start = cur
        end = cur + size - 1
        ranges.append((start, end))
        cur = end + 1
    return ranges


def sqs_send_batches(queue_url: str, messages: List[Dict[str, Any]]) -> int:
    """
    Send SQS messages in batches of 10.
    Returns how many messages were successfully sent.
    """
    sent = 0
    batch = []

    for idx, msg in enumerate(messages):
        batch.append({
            "Id": str(idx % 10),
            "MessageBody": json.dumps(msg, separators=(",", ":")),
        })

        if len(batch) == 10:
            resp = sqs.send_message_batch(QueueUrl=queue_url, Entries=batch)
            sent += len(resp.get("Successful", []))
            # If failures exist, raise: better to fail fast
            failed = resp.get("Failed", [])
            if failed:
                raise RuntimeError(f"SQS batch send failed: {failed}")
            batch = []

    if batch:
        resp = sqs.send_message_batch(QueueUrl=queue_url, Entries=batch)
        sent += len(resp.get("Successful", []))
        failed = resp.get("Failed", [])
        if failed:
            raise RuntimeError(f"SQS batch send failed: {failed}")

    return sent


def write_job_header(job_id: str, t_str: str, n_max: int, num_chunks: int, precision_bits: int) -> None:
    """
    Writes a job header item to DynamoDB so reducers can find job metadata.
    Key design (simple single-table pattern):
      PK = JOB#<job_id>
      SK = META
    """
    now = int(time.time())
    table.put_item(
        Item={
            "PK": f"JOB#{job_id}",
            "SK": "META",
            "job_id": job_id,
            "t": t_str,
            "n_max": n_max,
            "total_chunks": num_chunks,
            "precision_bits": precision_bits,
            "status": "DISPATCHED",
            "created_at": now,
        }
    )


# ----------------------------
# Lambda handler
# ----------------------------
def lambda_handler(event, context):
    """
    Expected event formats (either is fine):

    1) Direct:
       {"t": "1e20"}  or {"t": 1e20}

    2) Batch:
       {"t_list": ["1e20", "2e20", ...]}

    Returns job_id, n_max, chunks created, etc.
    """

    # Parse t(s)
    t_list = []
    if "t_list" in event:
        t_list = event["t_list"]
    elif "t" in event:
        t_list = [event["t"]]
    else:
        raise ValueError("event must include 't' or 't_list'")

    results = []

    for t_in in t_list:
        # keep the original string representation as canonical
        if isinstance(t_in, (int, float)):
            t_val = float(t_in)
            t_str = f"{t_val:.0f}" if t_val >= 1e6 else str(t_val)
        else:
            t_str = str(t_in).strip()
            # try float for n_max computation
            t_val = float(t_str)

        n_max = compute_n_max(t_val)
        num_chunks = choose_num_chunks(n_max, MAX_WORKERS_HINT)

        # create job
        job_id = f"{JOB_PREFIX}-{uuid.uuid4().hex[:12]}"

        # store metadata for reducer
        write_job_header(job_id, t_str, n_max, num_chunks, COARSE_PREC_BITS)

        # chunk messages
        ranges = split_range_balanced(n_max, num_chunks)

        messages = []
        for chunk_id, (n_start, n_end) in enumerate(ranges):
            messages.append({
                "job_id": job_id,
                "chunk_id": chunk_id,
                "total_chunks": num_chunks,
                "t": t_str,  # string to avoid float issues
                "n_start": n_start,
                "n_end": n_end,
                "precision_bits": COARSE_PREC_BITS,
                "schema": "valley_chunk_v1",
            })

        sent = sqs_send_batches(SQS_QUEUE_URL, messages)

        results.append({
            "job_id": job_id,
            "t": t_str,
            "n_max": n_max,
            "chunks_planned": num_chunks,
            "chunks_sent": sent,
            "precision_bits": COARSE_PREC_BITS,
        })

    return {
        "ok": True,
        "results": results,
    }
