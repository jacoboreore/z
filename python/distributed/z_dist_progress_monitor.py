import os
import json
import time
from decimal import Decimal

import boto3
from boto3.dynamodb.conditions import Key

AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
JOBS_TABLE = os.environ.get("JOBS_TABLE", "valley-scanner-dev-jobs")
CHUNKS_TABLE = os.environ.get("CHUNKS_TABLE", "valley-scanner-dev-chunks")

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)


def json_safe(obj):
    if isinstance(obj, Decimal):
        return int(obj) if obj % 1 == 0 else float(obj)
    raise TypeError(f"Not JSON serializable: {type(obj)}")


def now_epoch():
    return int(time.time())


def get_expected_chunks(job):
    for k in ("number_of_chunks_requested", "num_chunks", "chunk_count", "expected_chunks"):
        if k in job:
            return int(job[k])
    return None


def query_chunks(job_id):
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


def lambda_handler(event, context):
    job_id = event.get("job_id")

    if not job_id and "queryStringParameters" in event:
        qs = event.get("queryStringParameters") or {}
        job_id = qs.get("job_id")

    if not job_id:
        return {
            "statusCode": 400,
            "body": json.dumps({"error": "Missing job_id"})
        }

    job_resp = jobs_table.get_item(Key={"job_id": str(job_id)})
    job = job_resp.get("Item")

    if not job:
        return {
            "statusCode": 404,
            "body": json.dumps({"error": "Job not found", "job_id": job_id})
        }

    chunks = query_chunks(job_id)

    expected = get_expected_chunks(job)
    total_records = len(chunks)

    counts = {}
    for c in chunks:
        status = c.get("status", "UNKNOWN")
        counts[status] = counts.get(status, 0) + 1

    done = counts.get("DONE", 0)
    failed = counts.get("FAILED", 0)
    running = counts.get("RUNNING", 0) + counts.get("PROCESSING", 0)
    pending = counts.get("PENDING", 0) + counts.get("CREATED", 0)

    denominator = expected or total_records or 1
    progress_pct = round((done / denominator) * 100, 2)

    created_at = int(job.get("created_at", job.get("createdAt", 0)) or 0)
    updated_at = int(job.get("updated_at", 0) or 0)
    elapsed_seconds = now_epoch() - created_at if created_at else None

    response = {
        "job_id": str(job_id),
        "job_status": job.get("status"),
        "expected_chunks": expected,
        "chunk_records_found": total_records,
        "chunks_done": done,
        "chunks_failed": failed,
        "chunks_running": running,
        "chunks_pending_or_created": pending,
        "status_counts": counts,
        "progress_pct": progress_pct,
        "created_at": created_at or None,
        "updated_at": updated_at or None,
        "elapsed_seconds": elapsed_seconds,
        "final_result_s3_key": job.get("final_result_s3_key"),
        "ztotal": str(job.get("ztotal")) if "ztotal" in job else None,
        "assembler_error": job.get("assembler_error"),
    }

    return {
        "statusCode": 200,
        "headers": {"Content-Type": "application/json"},
        "body": json.dumps(response, default=json_safe)
    }