# monitor_job.py
import os, json, boto3
from datetime import datetime, timezone
from botocore.exceptions import ClientError

EC2 = boto3.client("ec2")
S3 = boto3.client("s3")

S3_BUCKET = os.environ.get("S3_BUCKET", "jacob-z")   # made overridable
TIMEOUT_MIN = int(os.environ.get("TIMEOUT_MIN", "180"))  # 3h default

def _resp(status: int, body: dict):
    return {
        "statusCode": status,
        "headers": {
            "Content-Type": "application/json",
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Headers": "Content-Type,Authorization",
            "Access-Control-Allow-Methods": "GET,POST,OPTIONS",
        },
        "body": json.dumps(body),
    }

def _find_instance(job_id: str):
    res = EC2.describe_instances(
        Filters=[
            {"Name": "tag:z:jobId", "Values": [job_id]},
            {"Name": "instance-state-name", "Values": [
                "pending","running","stopping","stopped","shutting-down","terminated"
            ]},
        ]
    )
    for r in res.get("Reservations", []):
        for i in r.get("Instances", []):
            return i
    return None

def _s3_read_last_line_tail(bucket, key, chunk_size=8192, max_chunks=8):
    try:
        head = S3.head_object(Bucket=bucket, Key=key)
        size = head["ContentLength"]
        if size == 0:
            return ""

        buf = b""
        start = max(0, size - chunk_size)
        end = size - 1
        chunks = 0

        while True:
            resp = S3.get_object(Bucket=bucket, Key=key, Range=f"bytes={start}-{end}")
            buf = resp["Body"].read() + buf
            text = buf.decode("utf-8", "replace").rstrip("\r\n \t\0")
            # walk backward to find the first non-empty line
            for ln in reversed(text.splitlines()):
                if ln.strip():
                    return ln

            chunks += 1
            if start == 0 or chunks >= max_chunks:
                return ""  # nothing non-empty found
            end = start - 1
            start = max(0, start - chunk_size)

    except S3.exceptions.NoSuchKey as s3e:
        return f"Error with S3 key {s3e}"
    except Exception as e:
        return f"Generic S3 exception {bucket} {key} {e}"

def _read_progress_json(bucket: str, job_id: str):
    """Read jobs/{jobId}/progress.json; return dict or None."""
    key = f"jobs/{job_id}/progress.json"
    try:
        obj = S3.get_object(Bucket=bucket, Key=key)
        raw = obj["Body"].read()
        data = json.loads(raw.decode("utf-8", "replace") or "{}")
        # minimal validation / normalization
        out = {
            "jobId": data.get("jobId") or job_id,
            "processed": data.get("processed"),
            "zeros": data.get("zeros") or [],
            "progressUpdatedAt": data.get("updatedAt")
        }
        # Auto-fill updatedAt if missing (use object LastModified)
        if not out["progressUpdatedAt"]:
            lm = obj.get("LastModified")
            if lm:
                out["progressUpdatedAt"] = lm.astimezone(timezone.utc).isoformat()
        return out
    except S3.exceptions.NoSuchKey:
        return None
    except ClientError as e:
        # permission or other S3 error: log but don't fail status
        print(f"progress.json read error: {e}")
        return None
    except Exception as e:
        print(f"progress.json parse error: {e}")
        return None

def lambda_handler(event, context):
    # CORS preflight
    if event.get("requestContext", {}).get("http", {}).get("method") == "OPTIONS":
        return _resp(200, {"ok": True})

    # Accept both POST {jobId} and GET ?jobId=
    body = event.get("body")
    if isinstance(body, str) and body.strip():
        try:
            body = json.loads(body)
        except json.JSONDecodeError:
            body = {}
    elif body is None:
        body = {}

    job_id = body.get("jobId")
    if not job_id:
        # GET ?jobId=...
        raw_qs = event.get("rawQueryString") or ""
        if raw_qs:
            # crude parse to avoid importing urllib for v2 events
            for part in raw_qs.split("&"):
                if part.startswith("jobId="):
                    job_id = part.split("=", 1)[1]
                    job_id = job_id and boto3.utils.unquote_str(job_id)
                    break

    if not job_id:
        return _resp(400, {"error": "Specify a Job ID"})

    # Always try to read progress snapshot
    progress = _read_progress_json(S3_BUCKET, job_id)

    inst = _find_instance(job_id)
    if not inst:
        last = _s3_read_last_line_tail(S3_BUCKET, f"jobs/{job_id}/output/run.log")
        if "JOB_DONE" in last:
            status = "succeeded"
        elif "JOB_FAILED" in last:
            status = "failed"
        else:
            status = "unknown"
        return _resp(200, {
            "status": status,
            "message": last,
            "jobId": job_id,
            **({"processed": progress["processed"], "zeros": progress["zeros"], "progressUpdatedAt": progress["progressUpdatedAt"]} if progress else {})
        })

    state = inst["State"]["Name"]
    launch_time = inst["LaunchTime"]  # UTC datetime
    age_min = (datetime.now(timezone.utc) - launch_time).total_seconds() / 60.0

    if state in ("terminated", "shutting-down"):
        last = _s3_read_last_line_tail(S3_BUCKET, f"jobs/{job_id}/output/run.log")
        status = "succeeded" if "JOB_DONE" in last else ("failed" if "JOB_FAILED" in last else "stopped")
        return _resp(200, {
            "status": status,
            "message": last,
            "jobId": job_id,
            **({"processed": progress["processed"], "zeros": progress["zeros"], "progressUpdatedAt": progress["progressUpdatedAt"]} if progress else {})
        })

    # Optional: auto-terminate if over timeout
    if age_min > TIMEOUT_MIN:
        EC2.terminate_instances(InstanceIds=[inst["InstanceId"]])
        return _resp(200, {"status": "terminating", "reason": "timeout", "jobId": job_id})

    # Running
    return _resp(200, {
        "status": "running",
        "instanceId": inst["InstanceId"],
        "state": state,
        "ageMinutes": round(age_min, 1),
        "jobId": job_id,
        **({"processed": progress["processed"], "zeros": progress["zeros"], "progressUpdatedAt": progress["progressUpdatedAt"]} if progress else {})
    })
