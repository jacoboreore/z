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

def _s3_read_last_line(bucket, key):
    try:
        obj = S3.get_object(Bucket=bucket, Key=key)
        data = obj["Body"].read().decode("utf-8", errors="replace")
        lines = [ln for ln in data.splitlines() if ln.strip()]
        return lines[-1] if lines else ""
    except S3.exceptions.NoSuchKey:
        return ""
    except Exception:
        return ""

def handler(event, context):

    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    job_id = body.get("jobId")

    # job_id = (event.get("queryStringParameters") or {}).get("jobId")
    # if not job_id:
    #     return {"statusCode": 400, "body": json.dumps({"error": "jobId is required"})}

    inst = _find_instance(job_id)
    if not inst:
        # check S3 run.log
        last = _s3_read_last_line(S3_BUCKET, f"jobs/{job_id}/run.log")
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
        last = _s3_read_last_line(S3_BUCKET, f"jobs/{job_id}/run.log")
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
