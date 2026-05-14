import os
import json
import time
import logging
import boto3

from boto3.dynamodb.conditions import Key
from botocore.exceptions import ClientError

AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

JOBS_TABLE = os.environ.get("JOBS_TABLE", "valley-scanner-dev-jobs")
CHUNKS_TABLE = os.environ.get("CHUNKS_TABLE", "valley-scanner-dev-chunks")
STATUS_INDEX = os.environ.get("STATUS_INDEX", "status-index")

ASSEMBLER_LAMBDA = os.environ.get("ASSEMBLER_LAMBDA", "z_dist_monitor_assembler2_phase2")

MAX_JOBS_PER_RUN = int(os.environ.get("MAX_JOBS_PER_RUN", "25"))
MAX_ASYNC_INVOCATIONS_PER_RUN = int(os.environ.get("MAX_ASYNC_INVOCATIONS_PER_RUN", "10"))

log = logging.getLogger()
log.setLevel(logging.INFO)

for handler in log.handlers:
    handler.setLevel(logging.INFO)

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
lambda_client = boto3.client("lambda", region_name=AWS_REGION)

jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)


def now_epoch() -> int:
    return int(time.time())


def query_jobs_in_chunks_created(limit: int) -> list[dict]:
    """
    Query the GSI instead of scanning the whole jobs table.
    Assumes GSI:
      PK = status
      SK = updated_at
    """
    items = []
    resp = jobs_table.query(
        IndexName=STATUS_INDEX,
        KeyConditionExpression=Key("status").eq("CHUNKS_CREATED"),
        Limit=limit
    )
    items.extend(resp.get("Items", []))

    while "LastEvaluatedKey" in resp and len(items) < limit:
        resp = jobs_table.query(
            IndexName=STATUS_INDEX,
            KeyConditionExpression=Key("status").eq("CHUNKS_CREATED"),
            ExclusiveStartKey=resp["LastEvaluatedKey"],
            Limit=max(1, limit - len(items))
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


def extract_expected_chunks(job_item: dict) -> int:
    for key in ("number_of_chunks_requested", "num_chunks", "chunk_count", "expected_chunks"):
        if key in job_item:
            return int(job_item[key])
    raise KeyError(
        "Job item does not contain number_of_chunks_requested/num_chunks/chunk_count/expected_chunks"
    )


def summarize_chunk_status(chunk_items: list[dict]) -> dict:
    done = sum(1 for c in chunk_items if c.get("status") == "DONE")
    failed = sum(1 for c in chunk_items if c.get("status") == "FAILED")
    total = len(chunk_items)
    return {
        "total": total,
        "done": done,
        "failed": failed,
    }


def update_job_status(job_id: str, status: str, extra: dict | None = None) -> None:
    expr_names = {"#s": "status"}
    expr_values = {
        ":status": status,
        ":updated_at": now_epoch()
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


def invoke_assembler_async(job_id: str) -> dict:
    payload = {
        "job_id": str(job_id)
    }

    resp = lambda_client.invoke(
        FunctionName=ASSEMBLER_LAMBDA,
        InvocationType="Event",
        Payload=json.dumps(payload).encode("utf-8")
    )

    return {
        "status_code": resp.get("StatusCode"),
    }


def process_job(job_item: dict) -> dict:
    job_id = str(job_item["job_id"])
    expected_chunks = extract_expected_chunks(job_item)

    log.info("Checking job_id=%s expected_chunks=%s", job_id, expected_chunks)

    chunk_items = query_chunks_for_job(job_id)
    summary = summarize_chunk_status(chunk_items)

    log.info(
        "Chunk summary job_id=%s total=%s done=%s failed=%s expected=%s",
        job_id,
        summary["total"],
        summary["done"],
        summary["failed"],
        expected_chunks,
    )

    if summary["failed"] > 0:
        update_job_status(
            job_id,
            "FAILED",
            {
                "chunks_done": summary["done"],
                "chunks_failed": summary["failed"],
                "monitor_error": f"{summary['failed']} chunks in FAILED status",
            }
        )
        return {
            "job_id": job_id,
            "action": "marked_failed",
            "done_chunks": summary["done"],
            "failed_chunks": summary["failed"],
        }

    if summary["done"] != expected_chunks:
        return {
            "job_id": job_id,
            "action": "waiting",
            "done_chunks": summary["done"],
            "expected_chunks": expected_chunks,
        }

    # All chunks complete. Let the single-job assembler take over.
    resp = invoke_assembler_async(job_id)

    update_job_status(
        job_id,
        "ASSEMBLY_QUEUED",
        {
            "chunks_done": summary["done"],
            "chunks_failed": summary["failed"],
            "assembler_invoked_at": now_epoch(),
            "assembler_invoke_status_code": resp.get("status_code"),
        }
    )

    return {
        "job_id": job_id,
        "action": "assembler_invoked",
        "status_code": resp.get("status_code"),
        "done_chunks": summary["done"],
    }


def lambda_handler(event, context):
    log.info("Monitor ready jobs started")
    log.info("JOBS_TABLE=%s", JOBS_TABLE)
    log.info("CHUNKS_TABLE=%s", CHUNKS_TABLE)
    log.info("status_index=%s", STATUS_INDEX)
    log.info("ASSEMBLER_LAMBDA=%s", ASSEMBLER_LAMBDA)
    log.info("MAX_JOBS_PER_RUN=%s", MAX_JOBS_PER_RUN)

    jobs = query_jobs_in_chunks_created(MAX_JOBS_PER_RUN)
    log.info("Found %s job(s) in CHUNKS_CREATED", len(jobs))

    results = []
    invoked = 0

    for job_item in jobs:
        if invoked >= MAX_ASYNC_INVOCATIONS_PER_RUN:
            log.info("Reached MAX_ASYNC_INVOCATIONS_PER_RUN=%s", MAX_ASYNC_INVOCATIONS_PER_RUN)
            break

        try:
            result = process_job(job_item)
            results.append(result)

            if result.get("action") == "assembler_invoked":
                invoked += 1

        except Exception as e:
            job_id = str(job_item.get("job_id", "UNKNOWN"))
            log.exception("Unexpected error while processing job_id=%s: %s", job_id, e)
            results.append({
                "job_id": job_id,
                "action": "unexpected_error",
                "error": str(e),
            })

    return {
        "statusCode": 200,
        "body": json.dumps({
            "jobs_seen": len(jobs),
            "assemblers_invoked": invoked,
            "results": results,
        })
    }