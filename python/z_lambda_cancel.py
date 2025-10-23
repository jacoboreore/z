import os, json, boto3

EC2 = boto3.client("ec2")
S3  = boto3.client("s3")

S3_BUCKET = os.environ.get("S3_BUCKET", "")  # optional: to mark cancel in run.log

def _find_instance(job_id: str):
    res = EC2.describe_instances(
        Filters=[
            {"Name": "tag:z:jobId", "Values": [job_id]},
            {"Name": "instance-state-name", "Values": [
                "pending", "running", "stopping", "stopped", "shutting-down"
            ]},
        ]
    )
    for r in res.get("Reservations", []):
        for i in r.get("Instances", []):
            return i
    return None

def _append_cancel_mark(job_id: str):
    # Optional: append a line so your monitor can show "canceled"
    # Safe no-op if bucket is unset or write fails.
    if not S3_BUCKET:
        return
    key = f"jobs/{job_id}/run.log"
    try:
        # Prefix with a newline in case the file exists
        S3.put_object(
            Bucket=S3_BUCKET,
            Key=key,
            Body=b"\nJOB_CANCELED by API\n",
        )
    except Exception:
        pass

def lambda_handler(event, context):
    body = event.get("body")
    if isinstance(body, str):
        try:
            body = json.loads(body or "{}")
        except Exception:
            body = {}
    elif body is None:
        body = {}

    job_id = (body.get("jobId") or "").strip()
    if not job_id:
        return {"statusCode": 400, "body": json.dumps({"error": "Missing jobId"})}

    inst = _find_instance(job_id)
    if not inst:
        # Still mark the log; maybe it already terminated
        _append_cancel_mark(job_id)
        return {"statusCode": 200, "body": json.dumps({
            "status": "not_found",
            "message": "No active instance for this jobId"
        })}

    iid   = inst["InstanceId"]
    state = inst["State"]["Name"]

    # Already on the way down?
    if state in ("shutting-down", "terminated", "stopping", "stopped"):
        _append_cancel_mark(job_id)
        return {"statusCode": 200, "body": json.dumps({
            "status": "terminating",
            "instanceId": iid,
            "state": state
        })}

    # Terminate
    EC2.terminate_instances(InstanceIds=[iid])
    _append_cancel_mark(job_id)

    return {"statusCode": 200, "body": json.dumps({
        "status": "terminating",
        "instanceId": iid
    })}
