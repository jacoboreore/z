# monitor_job.py
import os, json, boto3
from datetime import datetime, timezone, timedelta

EC2 = boto3.client("ec2")
S3  = boto3.client("s3")

S3_BUCKET   = "jacob-z"
TIMEOUT_MIN = int(os.environ.get("TIMEOUT_MIN", "180"))  # 3h default

def _find_instance(job_id):
    res = EC2.describe_instances(
        Filters=[
            {"Name":"tag:z:jobId","Values":[job_id]},
            {"Name":"instance-state-name","Values":["pending","running","stopping","stopped","shutting-down","terminated"]}
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
        return f"Generic S3 exception {e}"


def lambda_handler(event, context):

    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    job_id = body.get("jobId")

    if not job_id or job_id is None:
        return {"statusCode": 400, "body": json.dumps({"error": "Specify a Job ID"})}

    inst = _find_instance(job_id)
    if not inst:
        # check S3 run.log
        last = _s3_read_last_line_tail(S3_BUCKET, f"jobs/{job_id}/run.log")
        if "JOB_DONE" in last:
            status = "succeeded"
        elif "JOB_FAILED" in last:
            status = "failed"
        else:
            status = "unknown"
        return {"statusCode": 200, "body": json.dumps({"status": status, "message": last})}

    state = inst["State"]["Name"]
    launch_time = inst["LaunchTime"]  # datetime in UTC
    age_min = (datetime.now(timezone.utc) - launch_time).total_seconds() / 60.0

    if state in ("terminated", "shutting-down"):
        last = _s3_read_last_line_tail(S3_BUCKET, f"jobs/{job_id}/run.log")
        status = "succeeded" if "JOB_DONE" in last else ("failed" if "JOB_FAILED" in last else "stopped")
        return {"statusCode": 200, "body": json.dumps({"status": status, "message": last})}

    # Optional: kill if over timeout
    if age_min > TIMEOUT_MIN:
        EC2.terminate_instances(InstanceIds=[inst["InstanceId"]])
        return {"statusCode": 200, "body": json.dumps({"status": "terminating", "reason": "timeout"})}

    return {"statusCode": 200, "body": json.dumps({
        "status": "running",
        "instanceId": inst["InstanceId"],
        "state": state,
        "ageMinutes": round(age_min, 1)
    })}
