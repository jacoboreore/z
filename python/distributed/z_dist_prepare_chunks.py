import json
import os
import math
import datetime
import logging
import traceback
from urllib.parse import unquote_plus

import boto3
from botocore.exceptions import ClientError


AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
JOBS_TABLE = os.environ["JOBS_TABLE"]
CHUNKS_TABLE = os.environ["CHUNKS_TABLE"]

ENABLE_SQS_DISPATCH = os.environ.get("ENABLE_SQS_DISPATCH", "false").lower() == "true"
QUEUE_URL = os.environ.get("QUEUE_URL", "")

DEFAULT_NUMBER_OF_CPUS = int(os.environ.get("DEFAULT_NUMBER_OF_CPUS", "96"))
DEFAULT_NUMBER_OF_CHUNKS = int(os.environ.get("DEFAULT_NUMBER_OF_CHUNKS", "10"))

# Logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

s3 = boto3.client("s3", region_name=AWS_REGION)
dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
jobs_table = dynamodb.Table(JOBS_TABLE)
chunks_table = dynamodb.Table(CHUNKS_TABLE)
sqs = boto3.client("sqs", region_name=AWS_REGION) if ENABLE_SQS_DISPATCH and QUEUE_URL else None


logger.info(
    "Lambda cold start: region=%s jobs_table=%s chunks_table=%s enable_sqs=%s queue_url_set=%s default_cpus=%s default_chunks=%s",
    AWS_REGION,
    JOBS_TABLE,
    CHUNKS_TABLE,
    ENABLE_SQS_DISPATCH,
    bool(QUEUE_URL),
    DEFAULT_NUMBER_OF_CPUS,
    DEFAULT_NUMBER_OF_CHUNKS,
)


def parse_precalc_text(text: str) -> dict:
    logger.info("parse_precalc_text: input length=%s chars", len(text))

    data = {}

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        data[k.strip()] = v.strip()

    logger.info("parse_precalc_text: parsed keys=%s", sorted(list(data.keys())))

    required = ["N", "t", "theta", "gabcke", "useC1C2", "prec"]
    missing = [k for k in required if k not in data]
    if missing:
        logger.error("parse_precalc_text: missing required keys=%s", missing)
        raise ValueError(f"pre_calc.txt missing required keys: {missing}")

    data["N"] = int(data["N"])
    data["useC1C2"] = int(data["useC1C2"])
    data["prec"] = int(data["prec"])

    if data["N"] <= 0:
        logger.error("parse_precalc_text: invalid N=%s", data["N"])
        raise ValueError("N must be > 0")

    logger.info(
        "parse_precalc_text: success N=%s t=%s useC1C2=%s prec=%s",
        data["N"],
        data["t"],
        data["useC1C2"],
        data["prec"],
    )
    return data


def extract_job_id_from_key(key: str) -> str:
    logger.info("extract_job_id_from_key: raw key=%s", key)

    parts = key.strip("/").split("/")
    if len(parts) < 3:
        logger.error("extract_job_id_from_key: unexpected key format=%s", key)
        raise ValueError(f"Unexpected S3 key format: {key}")
    if parts[0] != "jobs":
        logger.error("extract_job_id_from_key: key does not start with jobs/ key=%s", key)
        raise ValueError(f"Expected key to start with jobs/: {key}")

    job_id = parts[1]
    logger.info("extract_job_id_from_key: extracted job_id=%s", job_id)
    return job_id


def build_ranges(N: int, number_of_chunks: int) -> list[dict]:
    logger.info("build_ranges: N=%s number_of_chunks=%s", N, number_of_chunks)

    if N <= 0:
        logger.error("build_ranges: invalid N=%s", N)
        raise ValueError("N must be > 0")
    if number_of_chunks <= 0:
        logger.error("build_ranges: invalid number_of_chunks=%s", number_of_chunks)
        raise ValueError("number_of_chunks must be > 0")

    chunk_size = math.ceil(N / number_of_chunks)
    logger.info("build_ranges: computed chunk_size=%s", chunk_size)

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

    if ranges:
        logger.info(
            "build_ranges: generated %s ranges first=%s last=%s",
            len(ranges),
            ranges[0],
            ranges[-1],
        )
    else:
        logger.warning("build_ranges: generated zero ranges")

    return ranges


def validate_ranges(ranges: list[dict], N: int) -> None:
    logger.info("validate_ranges: validating %s ranges for N=%s", len(ranges), N)

    if not ranges:
        logger.error("validate_ranges: no ranges generated")
        raise ValueError("No ranges generated")

    expected_n0 = 1

    for idx, row in enumerate(ranges):
        n0 = row["n0"]
        n1 = row["n1"]

        if n0 != expected_n0:
            logger.error(
                "validate_ranges: gap/overlap at chunk=%s expected_n0=%s actual_n0=%s row=%s",
                idx, expected_n0, n0, row
            )
            raise ValueError(
                f"Gap or overlap detected at chunk {idx}: expected n0={expected_n0}, got {n0}"
            )

        if n1 < n0:
            logger.error(
                "validate_ranges: invalid chunk=%s n0=%s n1=%s row=%s",
                idx, n0, n1, row
            )
            raise ValueError(f"Invalid chunk {idx}: n1={n1} < n0={n0}")

        expected_n0 = n1 + 1

    if ranges[-1]["n1"] != N:
        logger.error(
            "validate_ranges: coverage error last_n1=%s N=%s",
            ranges[-1]["n1"], N
        )
        raise ValueError(f"Coverage error: last chunk ends at {ranges[-1]['n1']} but N={N}")

    logger.info("validate_ranges: success")


def get_existing_job(job_id: str):
    logger.info("get_existing_job: checking job_id=%s", job_id)
    resp = jobs_table.get_item(Key={"job_id": job_id})
    item = resp.get("Item")
    logger.info("get_existing_job: exists=%s for job_id=%s", bool(item), job_id)
    return item


def put_job_if_new(
    job_id: str,
    pre: dict,
    number_of_cpus: int,
    number_of_chunks_requested: int,
    number_of_chunks_created: int,
    bucket: str,
    key: str
) -> bool:
    logger.info(
        "put_job_if_new: job_id=%s bucket=%s key=%s cpus=%s chunks_requested=%s chunks_created=%s",
        job_id, bucket, key, number_of_cpus, number_of_chunks_requested, number_of_chunks_created
    )

    now = datetime.datetime.utcnow().isoformat() + "Z"

    item = {
        "job_id": job_id,
        "created_at": now,
        "status": "CHUNKS_CREATED",
        "t": pre["t"],
        "N": pre["N"],
        "theta": pre["theta"],
        "gabcke": pre["gabcke"],
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
        logger.info("put_job_if_new: inserted new job_id=%s", job_id)
        return True
    except ClientError as e:
        code = e.response["Error"]["Code"]
        if code == "ConditionalCheckFailedException":
            logger.warning("put_job_if_new: job already exists job_id=%s", job_id)
            return False

        logger.exception("put_job_if_new: dynamodb error for job_id=%s", job_id)
        raise


def put_chunks(
    job_id: str,
    pre: dict,
    number_of_cpus: int,
    ranges: list[dict]
) -> None:
    logger.info("put_chunks: writing %s chunks for job_id=%s", len(ranges), job_id)

    now = datetime.datetime.utcnow().isoformat() + "Z"

    with chunks_table.batch_writer() as batch:
        for idx, row in enumerate(ranges):
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

            if idx < 3 or idx == len(ranges) - 1:
                logger.info("put_chunks: queued chunk row=%s", row)

    logger.info("put_chunks: finished writing chunks for job_id=%s", job_id)


def dispatch_chunks_to_sqs(job_id: str, ranges: list[dict], pre: dict, number_of_cpus: int):
    logger.info(
        "dispatch_chunks_to_sqs: enable_sqs=%s queue_url_set=%s job_id=%s range_count=%s",
        ENABLE_SQS_DISPATCH,
        bool(QUEUE_URL),
        job_id,
        len(ranges)
    )

    if not ENABLE_SQS_DISPATCH or not sqs or not QUEUE_URL:
        logger.info("dispatch_chunks_to_sqs: skipped")
        return {"dispatched": False, "count": 0}

    count = 0

    for idx, row in enumerate(ranges):
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

        if idx < 3 or idx == len(ranges) - 1:
            logger.info(
                "dispatch_chunks_to_sqs: sent chunk_id=%s chunk_index=%s",
                row["chunk_id"], row["chunk_index"]
            )

    logger.info("dispatch_chunks_to_sqs: dispatched count=%s for job_id=%s", count, job_id)
    return {"dispatched": True, "count": count}


def process_record(record: dict) -> dict:
    logger.info("process_record: raw record=%s", json.dumps(record, default=str))

    event_name = record.get("eventName", "")
    logger.info("process_record: event_name=%s", event_name)

    # Support both classic S3 notifications and EventBridge-style S3 events
    if "s3" in record:
        logger.info("process_record: detected classic S3 notification shape")
        bucket = record["s3"]["bucket"]["name"]
        key = unquote_plus(record["s3"]["object"]["key"])
    elif "detail" in record:
        logger.info("process_record: detected EventBridge shape")
        detail = record["detail"]
        bucket = detail["bucket"]["name"]
        key = unquote_plus(detail["object"]["key"])
        # EventBridge detail-type is different, so eventName may be absent
        if not event_name:
            event_name = "ObjectCreated:EventBridge"
    else:
        logger.error("process_record: unsupported event structure keys=%s", list(record.keys()))
        return {"skipped": True, "reason": "Unsupported event structure"}

    logger.info("process_record: bucket=%s key=%s", bucket, key)

    if not event_name.startswith("ObjectCreated:"):
        logger.warning("process_record: unsupported eventName=%s", event_name)
        return {"skipped": True, "reason": f"Unsupported eventName={event_name}"}

    if not key.endswith("pre_calc.txt"):
        logger.warning("process_record: not a pre_calc file key=%s", key)
        return {"skipped": True, "reason": f"Not a pre_calc file: {key}"}

    job_id = extract_job_id_from_key(key)

    logger.info("process_record: fetching S3 object bucket=%s key=%s", bucket, key)
    obj = s3.get_object(Bucket=bucket, Key=key)
    text = obj["Body"].read().decode("utf-8")
    logger.info("process_record: fetched S3 object length=%s chars", len(text))

    pre = parse_precalc_text(text)

    number_of_cpus = DEFAULT_NUMBER_OF_CPUS
    number_of_chunks = DEFAULT_NUMBER_OF_CHUNKS

    logger.info(
        "process_record: using number_of_cpus=%s number_of_chunks=%s for job_id=%s",
        number_of_cpus, number_of_chunks, job_id
    )

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
        logger.warning("process_record: duplicate job detected job_id=%s", job_id)
        return {
            "skipped": True,
            "reason": f"Job already exists, duplicate event ignored for {job_id}",
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

    result = {
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

    logger.info("process_record: success result=%s", json.dumps(result, default=str))
    return result


def lambda_handler(event, context):
    logger.info("lambda_handler: invoked")
    logger.info("lambda_handler: event=%s", json.dumps(event, default=str))

    try:
        results = []

        # Support direct EventBridge event (single event, no Records)
        if "Records" in event:
            logger.info("lambda_handler: processing Records count=%s", len(event["Records"]))
            for idx, record in enumerate(event["Records"]):
                logger.info("lambda_handler: processing record index=%s", idx)
                results.append(process_record(record))
        else:
            logger.info("lambda_handler: no Records found, treating full event as one record")
            results.append(process_record(event))

        final = {
            "statusCode": 200,
            "body": json.dumps({"results": results}, default=str)
        }

        logger.info("lambda_handler: completed successfully final=%s", json.dumps(final, default=str))
        return final

    except Exception as e:
        logger.error("lambda_handler: exception=%s", str(e))
        logger.error("lambda_handler: traceback=%s", traceback.format_exc())

        return {
            "statusCode": 500,
            "body": json.dumps({
                "error": str(e),
                "traceback": traceback.format_exc()
            })
        }