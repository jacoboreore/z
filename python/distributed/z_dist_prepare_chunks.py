import json
import os
import math
import datetime
from urllib.parse import unquote_plus

import boto3
from botocore.exceptions import ClientError
from boto3.dynamodb.conditions import Key


AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
JOBS_TABLE = os.environ["JOBS_TABLE"]
CHUNKS_TABLE = os.environ["CHUNKS_TABLE"]

# Optional:
# if true, Lambda 2 will also enqueue chunk messages to SQS after writing DynamoDB
ENABLE_SQS_DISPATCH = os.environ.get("ENABLE_SQS_DISPATCH", "false").lower() == "true"
QUEUE_URL = os.environ.get("QUEUE_URL", "")

s3 = boto3.client("s3", region_name=AWS_REGION)
dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)
sqs = boto3.client("sqs", region_name=AWS_REGION) if ENABLE_SQS_DISPATCH and QUEUE_URL else None

def parse_precalc_text(text: str) -> dict:
    """
    Parse pre_calc.txt of the form:

    N=11894
    t=888818888
    theta=...
    gabcke=...
    useC1C2=1
    prec=256
    """
    data = {}

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        data[k.strip()] = v.strip()

    required = ["N", "t", "theta", "gabcke", "useC1C2", "prec"]
    missing = [k for k in required if k not in data]
    if missing:
        raise ValueError(f"pre_calc.txt missing required keys: {missing}")

    data["N"] = int(data["N"])
    data["useC1C2"] = int(data["useC1C2"])
    data["prec"] = int(data["prec"])

    if data["N"] <= 0:
        raise ValueError("N must be > 0")

    return data

def extract_job_id_from_key(key: str) -> str:
    """
    Expected key:
      jobs/<job_id>/pre_calc.txt
    """
    parts = key.strip("/").split("/")
    if len(parts) < 3:
        raise ValueError(f"Unexpected S3 key format: {key}")
    if parts[0] != "jobs":
        raise ValueError(f"Expected key to start with jobs/: {key}")
    return parts[1]

def build_ranges(N: int, number_of_chunks: int) -> list[dict]:
    """
    Build inclusive chunk ranges matching the working C++ code:

      for (n = n0; n <= n1; ++n)

    Example:
      N=11894, chunks=3
      [1,3965], [3966,7930], [7931,11894]
    """
    if N <= 0:
        raise ValueError("N must be > 0")
    if number_of_chunks <= 0:
        raise ValueError("number_of_chunks must be > 0")

    chunk_size = math.ceil(N / number_of_chunks)
    ranges = []

    n0 = 1
    chunk_index = 0

    while n0 <= N:
        n1 = min(n0 + chunk_size - 1, N)

        ranges.append({
            "chunk_id": f"chunk-{chunk_index:06d}",
            "chunk_index": chunk_index,
            "n0": n0,
            "n1": n1,
            "size": n1 - n0 + 1
        })

        n0 = n1 + 1
        chunk_index += 1

    validate_ranges(ranges, N)
    return ranges

def validate_ranges(ranges: list[dict], N: int) -> None:
    if not ranges:
        raise ValueError("No ranges generated")

    expected_n0 = 1

    for idx, row in enumerate(ranges):
        n0 = row["n0"]
        n1 = row["n1"]

        if n0 != expected_n0:
            raise ValueError(
                f"Gap or overlap detected at chunk {idx}: expected n0={expected_n0}, got {n0}"
            )

        if n1 < n0:
            raise ValueError(
                f"Invalid chunk {idx}: n1={n1} < n0={n0}"
            )

        expected_n0 = n1 + 1

    if ranges[-1]["n1"] != N:
        raise ValueError(
            f"Coverage error: last chunk ends at {ranges[-1]['n1']} but N={N}"
        )

def get_existing_job(job_id: str):
    resp = jobs_table.get_item(Key={"job_id": job_id})
    return resp.get("Item")

def put_job_if_new(
    job_id: str,
    pre: dict,
    number_of_cpus: int,
    number_of_chunks_requested: int,
    number_of_chunks_created: int,
    bucket: str,
    key: str
) -> bool:
    """
    Returns True if inserted, False if already existed.
    Helps make the Lambda idempotent for duplicate S3 events.
    """
    now = datetime.datetime.utcnow().isoformat() + "Z"

    item = {
        "job_id": job_id,
        "created_at": now,
        "status": "CHUNKS_CREATED",
        "t": pre["t"],                 # keep as string for full precision
        "N": pre["N"],
        "theta": pre["theta"],         # keep as string
        "gabcke": pre["gabcke"],       # keep as string
        "useC1C2": pre["useC1C2"],
        "prec": pre["prec"],
        "number_of_cpus": number_of_cpus,
        "number_of_chunks_requested": number_of_chunks_requested,
        "number_of_chunks_created": number_of_chunks_created,
        "completed_chunks": 0,
        "failed_chunks": 0,
        "precalc_bucket": bucket,
        "precalc_key": key,
    }

    try:
        jobs_table.put_item(
            Item=item,
            ConditionExpression="attribute_not_exists(job_id)"
        )
        return True
    except ClientError as e:
        if e.response["Error"]["Code"] == "ConditionalCheckFailedException":
            return False
        raise

def put_chunks(
    job_id: str,
    pre: dict,
    number_of_cpus: int,
    ranges: list[dict]
) -> None:
    now = datetime.datetime.utcnow().isoformat() + "Z"

    with chunks_table.batch_writer() as batch:
        for row in ranges:
            batch.put_item(
                Item={
                    "job_id": job_id,
                    "chunk_id": row["chunk_id"],
                    "created_at": now,
                    "status": "PENDING",
                    "t": pre["t"],
                    "N": pre["N"],
                    "theta": pre["theta"],
                    "gabcke": pre["gabcke"],
                    "useC1C2": pre["useC1C2"],
                    "prec": pre["prec"],
                    "number_of_cpus": number_of_cpus,
                    "chunk_index": row["chunk_index"],
                    "n0": row["n0"],
                    "n1": row["n1"],
                    "size": row["size"]
                }
            )

def dispatch_chunks_to_sqs(job_id: str, ranges: list[dict], pre: dict, number_of_cpus: int):
    if not ENABLE_SQS_DISPATCH or not sqs or not QUEUE_URL:
        return {"dispatched": False, "count": 0}

    count = 0

    for row in ranges:
        message = {
            "job_id": job_id,
            "chunk_id": row["chunk_id"],
            "chunk_index": row["chunk_index"],
            "n0": row["n0"],
            "n1": row["n1"],
            "size": row["size"],
            "t": pre["t"],
            "N": pre["N"],
            "theta": pre["theta"],
            "gabcke": pre["gabcke"],
            "useC1C2": pre["useC1C2"],
            "prec": pre["prec"],
            "number_of_cpus": number_of_cpus
        }

        sqs.send_message(
            QueueUrl=QUEUE_URL,
            MessageBody=json.dumps(message),
            MessageAttributes={
                "job_id": {"DataType": "String", "StringValue": job_id},
                "chunk_id": {"DataType": "String", "StringValue": row["chunk_id"]},
            }
        )
        count += 1

    return {"dispatched": True, "count": count}

def process_record(record: dict) -> dict:
    event_name = record.get("eventName", "")
    if not event_name.startswith("ObjectCreated:"):
        return {"skipped": True, "reason": f"Unsupported eventName={event_name}"}

    bucket = record["s3"]["bucket"]["name"]
    key = unquote_plus(record["s3"]["object"]["key"])

    if not key.endswith("pre_calc.txt"):
        return {"skipped": True, "reason": f"Not a pre_calc file: {key}"}

    # Example key: jobs/<job_id>/pre_calc.txt
    job_id = extract_job_id_from_key(key)

    # Pull the file
    obj = s3.get_object(Bucket=bucket, Key=key)
    text = obj["Body"].read().decode("utf-8")

    pre = parse_precalc_text(text)

    # Optional override: if you want CPUs/chunk-count as metadata on the object later,
    # you can read x-amz-meta-* here. For now use simple defaults or env vars.
    # Better still: Lambda 1 can upload a job manifest too.
    number_of_cpus = int(os.environ.get("DEFAULT_NUMBER_OF_CPUS", "8"))
    number_of_chunks = int(os.environ.get("DEFAULT_NUMBER_OF_CHUNKS", "4"))

    ranges = build_ranges(pre["N"], number_of_chunks)

    inserted = put_job_if_new(
        job_id=job_id,
        pre=pre,
        number_of_cpus=number_of_cpus,
        number_of_chunks_requested=number_of_chunks,
        number_of_chunks_created=len(ranges),
        bucket=bucket,
        key=key
    )

    if not inserted:
        return {
            "skipped": True,
            "reason": f"Job already exists, duplicate S3 event ignored for {job_id}",
            "job_id": job_id
        }

    put_chunks(
        job_id=job_id,
        pre=pre,
        number_of_cpus=number_of_cpus,
        ranges=ranges
    )

    dispatch_info = dispatch_chunks_to_sqs(
        job_id=job_id,
        ranges=ranges,
        pre=pre,
        number_of_cpus=number_of_cpus
    )

    return {
        "job_id": job_id,
        "bucket": bucket,
        "key": key,
        "N": pre["N"],
        "t": pre["t"],
        "number_of_cpus": number_of_cpus,
        "number_of_chunks_requested": number_of_chunks,
        "number_of_chunks_created": len(ranges),
        "sample_ranges": ranges[:5],
        "sqs": dispatch_info
    }

def lambda_handler(event, context):
    try:
        if "Records" not in event:
            return {
                "statusCode": 400,
                "body": json.dumps({"error": "Expected S3 event with Records"})
            }

        results = []
        for record in event["Records"]:
            results.append(process_record(record))

        return {
            "statusCode": 200,
            "body": json.dumps({"results": results}, default=str)
        }

    except Exception as e:
        return {
            "statusCode": 500,
            "body": json.dumps({"error": str(e)})
        }