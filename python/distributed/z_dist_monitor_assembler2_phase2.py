import json
import os
import time
import logging
from decimal import Decimal, getcontext

import boto3
from boto3.dynamodb.conditions import Key
from botocore.exceptions import ClientError

getcontext().prec = 100

AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

JOBS_TABLE = os.environ.get("JOBS_TABLE", "valley-scanner-dev-jobs")
CHUNKS_TABLE = os.environ.get("CHUNKS_TABLE", "valley-scanner-dev-chunks")

RESULTS_BUCKET = os.environ.get("RESULTS_BUCKET", "jacob-z")

PRECALC_PREFIX = os.environ.get("PRECALC_PREFIX", "jobs")
CHUNK_RESULTS_PREFIX = os.environ.get("CHUNK_RESULTS_PREFIX", "z-dist-results")
FINAL_RESULT_FILENAME = os.environ.get("FINAL_RESULT_FILENAME", "final_result.json")

ASSEMBLING_STALE_SECONDS = int(os.environ.get("ASSEMBLING_STALE_SECONDS", "3600"))

log = logging.getLogger()
log.setLevel(logging.INFO)

for handler in log.handlers:
    handler.setLevel(logging.INFO)

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
s3 = boto3.client("s3", region_name=AWS_REGION)

jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)


def now_epoch() -> int:
    return int(time.time())


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


def get_job(job_id: str) -> dict | None:
    resp = jobs_table.get_item(Key={"job_id": str(job_id)})
    return resp.get("Item")


def try_claim_or_reclaim_job(job_id: str) -> bool:
    """
    Claim a job for assembly.

    Allowed transitions:
      CHUNKS_CREATED  -> ASSEMBLING
      ASSEMBLY_QUEUED -> ASSEMBLING
      ASSEMBLING(stale) -> ASSEMBLING
    """
    now = int(time.time())

    # First try the normal claim from ASSEMBLY_QUEUED.
    for ready_status in ("ASSEMBLY_QUEUED", "CHUNKS_CREATED"):
        try:
            jobs_table.update_item(
                Key={"job_id": str(job_id)},
                UpdateExpression="SET #s = :assembling, updated_at = :updated_at",
                ConditionExpression="#s = :ready_status",
                ExpressionAttributeNames={"#s": "status"},
                ExpressionAttributeValues={
                    ":assembling": "ASSEMBLING",
                    ":ready_status": ready_status,
                    ":updated_at": now
                }
            )
            log.info("Claimed job_id=%s from status=%s", job_id, ready_status)
            return True
        except ClientError as e:
            if e.response["Error"]["Code"] != "ConditionalCheckFailedException":
                raise

    # If claim failed, check for stale ASSEMBLING.
    job = get_job(job_id)
    if not job:
        log.warning("Job not found while trying to claim job_id=%s", job_id)
        return False

    if job.get("status") != "ASSEMBLING":
        log.info("Job not claimable job_id=%s status=%s", job_id, job.get("status"))
        return False

    age = now - int(job.get("updated_at", 0))
    if age < ASSEMBLING_STALE_SECONDS:
        log.info("Job already assembling and not stale job_id=%s age=%s", job_id, age)
        return False

    jobs_table.update_item(
        Key={"job_id": str(job_id)},
        UpdateExpression="SET #s = :assembling, updated_at = :updated_at",
        ConditionExpression="#s = :assembling_current",
        ExpressionAttributeNames={"#s": "status"},
        ExpressionAttributeValues={
            ":assembling": "ASSEMBLING",
            ":assembling_current": "ASSEMBLING",
            ":updated_at": now
        }
    )

    log.info("Reclaimed stale ASSEMBLING job_id=%s age=%s", job_id, age)
    return True


def query_chunks_for_job(job_id: str) -> list[dict]:
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
    raise KeyError(
        "Job item does not contain number_of_chunks_requested/num_chunks/chunk_count/expected_chunks"
    )


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

    log.info("Expected chunks=%s done=%s failed=%s total=%s",
         expected_chunks, len(done_chunks), len(failed_chunks), len(chunk_items))

    log.info("Reading precalc: s3://%s/%s", RESULTS_BUCKET, precalc_key)

    for chunk in ordered_chunks:
        chunk_id = str(chunk["chunk_id"])
        chunk_index = int(chunk["chunk_index"])
        key = build_chunk_result_key(job_id, chunk_id)

        log.info("Reading chunk_index=%s chunk_id=%s key=s3://%s/%s",
             chunk_index, chunk_id, RESULTS_BUCKET, key)
        
        content = read_s3_text(RESULTS_BUCKET, key)
        rows = parse_chunk_results(content)

        row = next((r for r in rows if r["chunk_id"] == chunk_index), None)
        if row is None:
            raise ValueError(
                f"Chunk result row not found in S3 file for job_id={job_id} "
                f"chunk_id={chunk_id} chunk_index={chunk_index}"
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
        "assembled_at": now_epoch(),
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

    if "t" in job_item:
        final_payload["t"] = str(job_item["t"])
    elif "t" in precalc:
        final_payload["t"] = str(precalc["t"])

    if "N" in job_item:
        final_payload["N"] = int(job_item["N"])
    elif precalc_N is not None:
        final_payload["N"] = int(precalc_N)

    return final_payload


def lambda_handler(event, context):
    job_id = str(event["job_id"])
    log.info("===== ASSEMBLER START job_id=%s =====", job_id)
    log.info("event=%s", json.dumps(event))
    log.info("RESULTS_BUCKET=%s PRECALC_PREFIX=%s CHUNK_RESULTS_PREFIX=%s",
             RESULTS_BUCKET, PRECALC_PREFIX, CHUNK_RESULTS_PREFIX)

    try:
        if not try_claim_or_reclaim_job(job_id):
            log.info("Could not claim job_id=%s", job_id)
            return {"statusCode": 200, "body": json.dumps({"job_id": job_id, "action": "not_claimed"})}

        job_item = get_job(job_id)
        log.info("job_item=%s", json.dumps(job_item, default=str))

        if not job_item:
            raise RuntimeError(f"Job not found after claim job_id={job_id}")

        chunk_items = query_chunks_for_job(job_id)
        log.info("chunks_found=%s", len(chunk_items))
        log.info("chunk_statuses=%s", json.dumps(
            {s: sum(1 for c in chunk_items if c.get("status") == s)
             for s in set(c.get("status", "MISSING") for c in chunk_items)}
        ))

        final_payload = assemble_job(job_item, chunk_items)
        final_key = build_final_result_key(job_id)

        log.info("Writing final result to s3://%s/%s", RESULTS_BUCKET, final_key)
        write_s3_json(RESULTS_BUCKET, final_key, final_payload)

        update_job_status(job_id, "DONE", {
            "chunks_done": final_payload["chunks_done"],
            "sum_chunks": final_payload["sum_chunks"],
            "zmain": final_payload["zmain"],
            "gabcke": final_payload["gabcke"],
            "ztotal": final_payload["ztotal"],
            "final_result_s3_key": final_key,
            "assembled_at": final_payload["assembled_at"],
            "assembler_error": "",
        })

        log.info("===== ASSEMBLER DONE job_id=%s ztotal=%s =====", job_id, final_payload["ztotal"])
        return {"statusCode": 200, "body": json.dumps({
            "job_id": job_id,
            "action": "assembled",
            "final_result_s3_key": final_key,
            "ztotal": final_payload["ztotal"]
        })}

    except Exception as e:
        log.exception("===== ASSEMBLER FAILED job_id=%s error=%s =====", job_id, e)
        try:
            update_job_status(job_id, "ASSEMBLY_FAILED", {
                "assembler_error": str(e)[:1000],
                "assembler_failed_at": now_epoch(),
            })
        except Exception as update_err:
            log.exception("Could not mark ASSEMBLY_FAILED job_id=%s error=%s", job_id, update_err)
        raise