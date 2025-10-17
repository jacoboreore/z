import os, json, re, boto3

S3_BUCKET = os.environ["S3_BUCKET"]
URL_TTL   = int(os.environ.get("URL_TTL", "900"))  # 15 min
s3 = boto3.client("s3")

def _payload(event):
    # Accept both direct-invoke and API GW shapes
    body = {}
    if isinstance(event.get("body"), str):
        try: body = json.loads(event["body"])
        except: body = {}
    elif isinstance(event.get("body"), dict):
        body = event["body"]
    
    qs = event.get("queryStringParameters") or {}
    body.update({k:v for k,v in qs.items() if v is not None})
    
    return body

def lambda_handler(event, context):
    p = _payload(event)
    job_id = (p.get("jobId") or "").strip()
    if not re.fullmatch(r"[0-9a-fA-F-]{36}", job_id):
        return {"statusCode": 400, "body": json.dumps({"error":"invalid jobId"})}

    prefix = f"jobs/{job_id}/"
    resp = s3.list_objects_v2(Bucket=S3_BUCKET, Prefix=prefix, MaxKeys=500)
    keys = {c["Key"] for c in resp.get("Contents", [])}

    want = [
        f"{prefix}run.log",
        f"{prefix}output/valley_walk.csv",
        f"{prefix}output/refined_zeros.csv",
        f"{prefix}output/valley_walk_for_plotting.png",
    ]

    files = []
    for k in sorted(keys):
        url = s3.generate_presigned_url(
            "get_object",
            Params={"Bucket": S3_BUCKET, "Key": k},
            ExpiresIn=URL_TTL
        )
        files.append({"key": k, "url": url})

    selected = []
    for k in want:
        if k in keys:
            url = s3.generate_presigned_url(
                "get_object", Params={"Bucket": S3_BUCKET, "Key": k}, ExpiresIn=URL_TTL
            )
            selected.append({"key": k, "url": url})

    return {"statusCode": 200, "body": json.dumps({"files": files, "selected": selected})}


## ENV VARS
# S3_BUCKET=jacob-z
# URL_TTL=900