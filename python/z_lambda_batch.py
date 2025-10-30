import os, json, base64, uuid, boto3
from string import Template
import logging


EC2 = boto3.client("ec2")
S3 = boto3.client('s3')

# --- Required env vars ---
SUBNET_ID             = os.environ.get("SUBNET_ID")
SECURITY_GROUP_ID     = os.environ["SECURITY_GROUP_ID"]
INSTANCE_PROFILE_NAME = os.environ["INSTANCE_PROFILE_NAME"]
S3_BUCKET             = os.environ["S3_BUCKET"]
AWS_REGION            = os.environ.get("AWS_REGION", "us-east-1")

# --- autosize CPU ---
INSTANCE_TYPE_SMALL   = os.environ.get("INSTANCE_TYPE_SMALL", "c7i.8xlarge")
INSTANCE_TYPE_LARGE   = os.environ.get("INSTANCE_TYPE_LARGE", "c7i.16xlarge")
INSTANCE_TYPE_XLARGE   = os.environ.get("INSTANCE_TYPE_XLARGE", "c7i.24xlarge")
INSTANCE_TYPE_HUGE   = os.environ.get("INSTANCE_TYPE_XLARGE", "c7i.48xlarge")

# Optional: name the thresholds explicitly (as floats)
THRESH_SMALL_MAX  = 1e8     # t <= 1e10 → SMALL
THRESH_LARGE_MAX  = 1e12     # (1e10, 1e14] → LARGE
THRESH_XLARGE_MAX = 1e16     # (1e14, 1e18) → XLARGE
THRESH_HUGE_MAX   = 1e21     # [1e18, 1e21] → HUGe (and above → HUGe)

#T_THRESHOLD           = float(os.environ.get("T_THRESHOLD", "1e10"))
USE_SPOT              = os.environ.get("USE_SPOT", "false").lower() == "true"

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def load_user_data_from_s3(bucket: str, key: str, variables: dict) -> str:
    try:
        obj = S3.get_object(Bucket=bucket, Key=key)
        text = obj["Body"].read().decode("utf-8")
    except S3.exceptions.NoSuchKey:
        raise RuntimeError(f"User-data template not found: s3://{bucket}/{key}")
    except Exception as e:
        raise RuntimeError(f"Failed to fetch user-data template: {e}")

    # Use safe_substitute so $JOB_ID, $USR, ${ZEROS_CSV}, etc. are left as-is if not provided
    tmpl = Template(text)
    return tmpl.safe_substitute(variables)

def _pick_instance_type(t_val: float) -> str:
    """
    Ranges:
      - t <= 1e10                     → SMALL
      - 1e10 < t <= 1e14              → LARGE
      - 1e14 < t <  1e18              → XLARGE
      - 1e18 <= t <= 1e21 (and above) → HUGE
    """
    if t_val <= THRESH_SMALL_MAX:
        choice = INSTANCE_TYPE_SMALL
    elif t_val <= THRESH_LARGE_MAX:
        choice = INSTANCE_TYPE_LARGE
    elif t_val < THRESH_XLARGE_MAX:  # strictly less than 1e18 to avoid gap
        choice = INSTANCE_TYPE_XLARGE
    else:
        # includes [1e18, 1e21] and anything > 1e21
        choice = INSTANCE_TYPE_HUGE

    logger.info(f"_pick_instance_type: t={t_val} → {choice}")
    return choice

def lambda_handler(event, context):
    # Expect {"image":"ami-...", "t": <number|string>, "zeros": <number|string>}
    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    image = body.get("image")
    t     = body.get("t")
    zeros = body.get("zeros")  # total zeros target

    if not image or t is None or zeros is None:
        return {"statusCode": 400, "body": json.dumps({"error": "image, zeros and t are required"})}

    try:
        t_val = float(t)
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "t must be numeric"})}

    try:
        zeros_val = int(float(zeros))  # accept "1000" or "1000.0"
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "zeros must be numeric"})}

    job_id = str(uuid.uuid4())
    instance_type = _pick_instance_type(t_val)

    vars_map = {
        "S3_BUCKET": S3_BUCKET,
        "job_id": job_id,
        "t_val": t_val,
        "zeros_val": zeros_val,
    }

    template_key = "z-commands/all_commands"
    user_data = load_user_data_from_s3(S3_BUCKET, template_key, vars_map)

    params = {
        "ImageId": image,
        "InstanceType": instance_type,
        "MinCount": 1, "MaxCount": 1,
        "UserData": base64.b64encode(user_data.encode()).decode(),
        "IamInstanceProfile": {"Name": INSTANCE_PROFILE_NAME},
        "InstanceInitiatedShutdownBehavior": "terminate",
        "TagSpecifications": [{
            "ResourceType": "instance",
            "Tags": [
                {"Key":"Name",   "Value": f"z-job-{job_id}"},
                {"Key":"z:jobId","Value": job_id},
                {"Key":"z:t",    "Value": str(t_val)},
            ],
        }],
        "NetworkInterfaces": [{
            "DeviceIndex": 0,
            "SubnetId": SUBNET_ID,
            "Groups": [SECURITY_GROUP_ID],
            "AssociatePublicIpAddress": True
        }],
    }

    if USE_SPOT:
        params["InstanceMarketOptions"] = {
            "MarketType": "spot",
            "SpotOptions": {"SpotInstanceType": "one-time"}
        }

    resp = EC2.run_instances(**params)
    instance_id = resp["Instances"][0]["InstanceId"]

    return {
        "statusCode": 200,
        "body": json.dumps({
            "jobId": job_id,
            "instanceId": instance_id,
            "instanceType": instance_type,
            "message": "Launched. Schedulers (10s) enabled for JSON+S3; job will upload and halt."
        })
    }
