import os
import json
import logging
from datetime import datetime, timezone
import boto3
from botocore.exceptions import ClientError
from urllib.parse import parse_qs

logger = logging.getLogger()
logger.setLevel(logging.INFO)

S3 = boto3.client("s3")
BUCKET = os.environ["ARTIFACTS_BUCKET"]
BASE_PREFIX = os.getenv("ARTIFACTS_BASE_PREFIX", "jobs/") 

def _resp(status: int, body):
    return {
        "statusCode": status,
        "headers": {
            "Content-Type": "application/json",
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Headers": "Content-Type,Authorization",
            "Access-Control-Allow-Methods": "GET,OPTIONS",
        },
        "body": json.dumps(body),
    }

def lambda_handler(event, _ctx):
    if event.get("requestContext", {}).get("http", {}).get("method") == "OPTIONS":
        return _resp(200, {"ok": True})

    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    job_id = body.get("jobId")

    if not job_id or job_id is None:
        return {"statusCode": 400, "body": json.dumps({"error": "job_id are required"})}
    
    prefix = f"{BASE_PREFIX}{job_id}/output/"
    logger.info("Listing artifacts bucket=%s prefix=%s", BUCKET, prefix)

    try:
        keys = []
        paginator = S3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=BUCKET, Prefix=prefix):
            for obj in page.get("Contents", []):
                keys.append({
                    "Key": obj["Key"],
                    "Size": obj["Size"],
                    "LastModified": obj["LastModified"].astimezone(timezone.utc).isoformat()
                })
    except ClientError as e:
        logger.exception("List failed")
        return _resp(500, {"error": "list_failed", "detail": e.response.get("Error", {}).get("Message")})

    # presign each
    out = []
    for k in keys:
        try:
            url = S3.generate_presigned_url(
                ClientMethod="get_object",
                Params={"Bucket": BUCKET, "Key": k["Key"]},
                ExpiresIn=300,
            )
            out.append({
                "key": k["Key"],
                "size": k["Size"],
                "lastModified": k["LastModified"],
                "url": url,
                "expiresInSeconds": 300,
            })
        except ClientError as e:
            logger.warning("Presign failed for %s: %s", k["Key"], e)
            out.append({
                "key": k["Key"],
                "size": k["Size"],
                "lastModified": k["LastModified"],
                "url": None,
                "error": "presign_failed",
            })

    return _resp(200, {"jobId": job_id, "bucket": BUCKET, "prefix": prefix, "files": out})
