import json
import os
import time
import logging
from decimal import Decimal, getcontext

import boto3
from boto3.dynamodb.conditions import Attr, Key
from botocore.exceptions import ClientError

getcontext().prec = 100

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

JOBS_TABLE = os.environ.get("JOBS_TABLE", "valley-scanner-dev-jobs")
CHUNKS_TABLE = os.environ.get("CHUNKS_TABLE", "valley-scanner-dev-chunks")

RESULTS_BUCKET = os.environ.get("RESULTS_BUCKET", "jacob-z")

PRECALC_PREFIX = os.environ.get("PRECALC_PREFIX", "jobs")
CHUNK_RESULTS_PREFIX = os.environ.get("CHUNK_RESULTS_PREFIX", "z-dist-results")
FINAL_RESULT_FILENAME = os.environ.get("FINAL_RESULT_FILENAME", "final_result.json")

MAX_JOBS_PER_RUN = int(os.environ.get("MAX_JOBS_PER_RUN", "25"))

# -----------------------------------------------------------------------------
# LOGGING
# -----------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
log = logging.getLogger("monitor-assembler")

# -----------------------------------------------------------------------------
# AWS CLIENTS
# -----------------------------------------------------------------------------
dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
s3 = boto3.client("s3", region_name=AWS_REGION)

jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
def parse_decimal(s: str) -> Decimal:
    return Decimal(str(s).replace("E", "e"))


def parse_precalc(content: str) -> dict:
    data = {}
    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        data[k.strip()] = v.strip()

    if "N" in data:
        data["N"] = int(data["N"])
    if "gabcke" in data:
        data["gabcke_decimal"] = parse_decimal(data["gabcke"])

    return data


def parse_chunk_results(content: str) -> list[dict]:
    results = []
    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        row = {}
        parts = line.split(",")
        for p in parts:
            k, v = p.split("=", 1)
            row[k.strip()] = v.strip()

        results.append({
            "chunk_id": int(row["chunk_id"]),
            "n0": int(row["n0"]),
            "n1": int(row["n1"]),
            "sum": parse_decimal(row["sum"])
        })

    results.sort(key=lambda x: x["chunk_id"])
    return results


def read_s3_text(bucket: str, key: str) -> str:
    resp = s3.get_object(Bucket=bucket, Key=key)
    return resp["Body"].read().decode("utf-8")


def write_s3_json(bucket: str, key: str, payload: dict) -> None:
    s3.put_object(
        Bucket=bucket,
        Key=key,
        Body=json.dumps(payload, indent=2).encode("utf-8"),
        ContentType="application/json"
    )


def update_job_status(job_id: str, status: str, extra: dict | None = None) -> None:
    expr_names = {"#s": "status"}
    expr_values = {
        ":status": status,
        ":updated_at": int(time.time())
    }
    update_expr = "SET #s = :status, updated_at = :updated_at"

    if extra:
        for k, v in extra.items():
            if k == "status":
                continue
            expr_names[f"#{k}"] = k
            expr_values[f":{k}"] = v
            update_expr += f", #{k} = :{k}"

    jobs_table.update_item(
        Key={"job_id": str(job_id)},
        UpdateExpression=update_expr,
        ExpressionAttributeNames=expr_names,
        ExpressionAttributeValues=expr_values
    )


def try_mark_job_assembling(job_id: str) -> bool:
    """
    Move job from CHUNKS_CREATED to ASSEMBLING only if it is still CHUNKS_CREATED.
    Prevents duplicate assembly across overlapping Lambda invocations.
    """
    try:
        jobs_table.update_item(
            Key={"job_id": str(job_id)},
            UpdateExpression="SET #s = :assembling, updated_at = :updated_at",
            ConditionExpression="#s = :chunks_created",
            ExpressionAttributeNames={"#s": "status"},
            ExpressionAttributeValues={
                ":assembling": "ASSEMBLING",
                ":chunks_created": "CHUNKS_CREATED",
                ":updated_at": int(time.time())
            }
        )
        return True
    except ClientError as e:
        if e.response["Error"]["Code"] == "ConditionalCheckFailedException":
            return False
        raise


def scan_jobs_in_chunks_created(limit: int) -> list[dict]:
    items = []
    scan_kwargs = {
        "FilterExpression": Attr("status").eq("CHUNKS_CREATED")
    }

    resp = jobs_table.scan(**scan_kwargs)
    items.extend(resp.get("Items", []))

    while "LastEvaluatedKey" in resp and len(items) < limit:
        resp = jobs_table.scan(
            ExclusiveStartKey=resp["LastEvaluatedKey"],
            **scan_kwargs
        )
        items.extend(resp.get("Items", []))

    return items[:limit]


def query_chunks_for_job(job_id: str) -> list[dict]:
    """
    Assumes chunks table PK is job_id and SK is chunk_id.
    """
    items = []
    resp = chunks_table.query(
        KeyConditionExpression=Key("job_id").eq(str(job_id))
    )
    items.extend(resp.get("Items", []))

    while "LastEvaluatedKey" in resp:
        resp = chunks_table.query(
            KeyConditionExpression=Key("job_id").eq(str(job_id)),
            ExclusiveStartKey=resp["LastEvaluatedKey"]
        )
        items.extend(resp.get("Items", []))

    return items


def build_chunk_result_key(job_id: str, chunk_id: str) -> str:
    return f"{CHUNK_RESULTS_PREFIX}/{job_id}/{chunk_id}/chunk_results.txt"


def build_final_result_key(job_id: str) -> str:
    return f"{CHUNK_RESULTS_PREFIX}/{job_id}/{FINAL_RESULT_FILENAME}"


def build_precalc_key(job_id: str) -> str:
    return f"{PRECALC_PREFIX}/{job_id}/pre_calc.txt"


def extract_expected_chunks(job_item: dict) -> int:
    for key in ("number_of_chunks_requested", "num_chunks", "chunk_count", "expected_chunks"):
        if key in job_item:
            return int(job_item[key])
    raise KeyError("Job item does not contain number_of_chunks_requested/num_chunks/chunk_count/expected_chunks")


def assemble_job(job_item: dict, chunk_items: list[dict]) -> dict:
    job_id = str(job_item["job_id"])
    expected_chunks = extract_expected_chunks(job_item)

    done_chunks = [c for c in chunk_items if c.get("status") == "DONE"]
    failed_chunks = [c for c in chunk_items if c.get("status") == "FAILED"]

    if failed_chunks:
        raise RuntimeError(f"Job {job_id} has FAILED chunks: {len(failed_chunks)}")

    if len(done_chunks) != expected_chunks:
        raise RuntimeError(
            f"Job {job_id} not ready: done_chunks={len(done_chunks)} expected_chunks={expected_chunks}"
        )

    precalc_key = build_precalc_key(job_id)
    precalc_content = read_s3_text(RESULTS_BUCKET, precalc_key)
    precalc = parse_precalc(precalc_content)

    if "gabcke_decimal" not in precalc:
        raise ValueError(f"gabcke not found in pre_calc.txt for job {job_id}")

    gabcke = precalc["gabcke_decimal"]
    precalc_N = precalc.get("N")

    ordered_chunks = sorted(done_chunks, key=lambda x: int(x["chunk_index"]))

    chunk_details = []
    total = Decimal(0)

    for chunk in ordered_chunks:
        chunk_id = str(chunk["chunk_id"])
        chunk_index = int(chunk["chunk_index"])
        key = build_chunk_result_key(job_id, chunk_id)

        content = read_s3_text(RESULTS_BUCKET, key)
        rows = parse_chunk_results(content)

        row = next((r for r in rows if r["chunk_id"] == chunk_index), None)
        if row is None:
            raise ValueError(
                f"Chunk result row not found in S3 file for job_id={job_id} chunk_id={chunk_id} chunk_index={chunk_index}"
            )

        total += row["sum"]
        chunk_details.append({
            "chunk_id": chunk_id,
            "chunk_index": chunk_index,
            "n0": row["n0"],
            "n1": row["n1"],
            "chunk_sum": str(row["sum"]),
            "s3_key": key
        })

    zmain = Decimal(2) * total
    ztotal = zmain + gabcke

    final_payload = {
        "job_id": job_id,
        "status": "DONE",
        "assembled_at": int(time.time()),
        "number_of_chunks_requested": expected_chunks,
        "chunks_done": len(done_chunks),
        "sum_chunks": str(total),
        "zmain": str(zmain),
        "gabcke": str(gabcke),
        "ztotal": str(ztotal),
        "precalc_s3_key": precalc_key,
        "final_result_s3_key": build_final_result_key(job_id),
        "chunks": chunk_details
    }

    # Optional carry-through fields from jobs/precalc
    if "t" in job_item:
        final_payload["t"] = str(job_item["t"])
    elif "t" in precalc:
        final_payload["t"] = str(precalc["t"])

    if "N" in job_item:
        final_payload["N"] = int(job_item["N"])
    elif precalc_N is not None:
        final_payload["N"] = int(precalc_N)

    return final_payload


def process_job(job_item: dict) -> dict:
    job_id = str(job_item["job_id"])
    expected_chunks = extract_expected_chunks(job_item)

    log.info("Checking job_id=%s expected_chunks=%s", job_id, expected_chunks)

    chunk_items = query_chunks_for_job(job_id)
    done_count = sum(1 for c in chunk_items if c.get("status") == "DONE")
    failed_count = sum(1 for c in chunk_items if c.get("status") == "FAILED")

    log.info(
        "Job chunk status summary: job_id=%s total_chunks_seen=%s done=%s failed=%s expected=%s",
        job_id, len(chunk_items), done_count, failed_count, expected_chunks
    )

    if failed_count > 0:
        update_job_status(
            job_id,
            "FAILED",
            {
                "chunks_done": done_count,
                "chunks_failed": failed_count,
                "assembly_error": f"{failed_count} chunks in FAILED status"
            }
        )
        return {
            "job_id": job_id,
            "action": "marked_failed",
            "done_chunks": done_count,
            "failed_chunks": failed_count
        }

    if done_count != expected_chunks:
        return {
            "job_id": job_id,
            "action": "waiting",
            "done_chunks": done_count,
            "expected_chunks": expected_chunks
        }

    claimed = try_mark_job_assembling(job_id)
    if not claimed:
        log.info("Job already claimed by another invocation: job_id=%s", job_id)
        return {
            "job_id": job_id,
            "action": "skipped_already_claimed"
        }

    try:
        final_payload = assemble_job(job_item, chunk_items)
        final_key = build_final_result_key(job_id)

        write_s3_json(RESULTS_BUCKET, final_key, final_payload)

        update_job_status(
            job_id,
            "DONE",
            {
                "chunks_done": final_payload["chunks_done"],
                "sum_chunks": final_payload["sum_chunks"],
                "zmain": final_payload["zmain"],
                "gabcke": final_payload["gabcke"],
                "ztotal": final_payload["ztotal"],
                "final_result_s3_key": final_key,
                "assembled_at": final_payload["assembled_at"]
            }
        )

        log.info("Job assembled successfully: job_id=%s final_key=%s", job_id, final_key)

        return {
            "job_id": job_id,
            "action": "assembled",
            "final_result_s3_key": final_key,
            "ztotal": final_payload["ztotal"]
        }

    except Exception as e:
        log.exception("Assembly failed for job_id=%s: %s", job_id, e)
        update_job_status(job_id, "FAILED", {"assembly_error": str(e)})
        return {
            "job_id": job_id,
            "action": "assembly_failed",
            "error": str(e)
        }


# -----------------------------------------------------------------------------
# LAMBDA ENTRYPOINT
# -----------------------------------------------------------------------------
def lambda_handler(event, context):
    log.info("Monitor assembler started")
    log.info("JOBS_TABLE=%s", JOBS_TABLE)
    log.info("CHUNKS_TABLE=%s", CHUNKS_TABLE)
    log.info("RESULTS_BUCKET=%s", RESULTS_BUCKET)
    log.info("PRECALC_PREFIX=%s", PRECALC_PREFIX)
    log.info("CHUNK_RESULTS_PREFIX=%s", CHUNK_RESULTS_PREFIX)
    log.info("MAX_JOBS_PER_RUN=%s", MAX_JOBS_PER_RUN)

    jobs = scan_jobs_in_chunks_created(MAX_JOBS_PER_RUN)
    log.info("Found %s job(s) in CHUNKS_CREATED", len(jobs))

    results = []
    for job_item in jobs:
        try:
            results.append(process_job(job_item))
        except Exception as e:
            job_id = str(job_item.get("job_id", "UNKNOWN"))
            log.exception("Unexpected error while processing job_id=%s: %s", job_id, e)
            results.append({
                "job_id": job_id,
                "action": "unexpected_error",
                "error": str(e)
            })

    return {
        "statusCode": 200,
        "body": json.dumps({
            "processed_jobs": len(results),
            "results": results
        })
    }