import os, json, boto3
EC2 = boto3.client("ec2")

def handler(event, context):
    body = json.loads(event.get("body") or "{}")
    job_id = body.get("jobId")
    if not job_id:
        return {"statusCode": 400, "body": json.dumps({"error": "jobId is required"})}

    res = EC2.describe_instances(Filters=[{"Name":"tag:z:jobId","Values":[job_id]}])
    ids = []
    for r in res.get("Reservations", []):
        for i in r.get("Instances", []):
            ids.append(i["InstanceId"])
    if not ids:
        return {"statusCode": 200, "body": json.dumps({"ok": True, "message":"No instances found (already gone)."})}

    EC2.terminate_instances(InstanceIds=ids)
    return {"statusCode": 200, "body": json.dumps({"ok": True, "terminated": ids})}
