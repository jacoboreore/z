import os, json, base64, uuid, boto3
from string import Template
import logging

EC2 = boto3.client("ec2")
S3  = boto3.client("s3")

# --- Required env vars ---
SUBNET_ID             = os.environ.get("SUBNET_ID")
SECURITY_GROUP_ID     = os.environ["SECURITY_GROUP_ID"]
INSTANCE_PROFILE_NAME = os.environ["INSTANCE_PROFILE_NAME"]
S3_BUCKET             = os.environ["S3_BUCKET"]
AWS_REGION            = os.environ.get("AWS_REGION", "us-east-1")

# --- autosize CPU ---
INSTANCE_TYPE_SMALL   = os.environ.get("INSTANCE_TYPE_SMALL",  "c7i.8xlarge")
INSTANCE_TYPE_LARGE   = os.environ.get("INSTANCE_TYPE_LARGE",  "c7i.16xlarge")
INSTANCE_TYPE_XLARGE  = os.environ.get("INSTANCE_TYPE_XLARGE", "c7i.24xlarge")
INSTANCE_TYPE_HUGE    = os.environ.get("INSTANCE_TYPE_HUGE",   "c7i.48xlarge")

# Thresholds
THRESH_SMALL_MAX   = 1e8
THRESH_LARGE_MAX   = 1e12
THRESH_XLARGE_MAX  = 1e16
THRESH_HUGE_MAX    = 1e21

USE_SPOT = os.environ.get("USE_SPOT", "false").lower() == "true"

# --- Hard limits ---
T_MAX_HARD       = 9_999_999_999_999_999
ZEROS_MAX_HARD   = 25
ZEROS_MIN_HARD   = 1

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

    tmpl = Template(text)
    return tmpl.safe_substitute(variables)

def _pick_instance_type(t_val: float) -> str:
    if t_val <= THRESH_SMALL_MAX:
        return INSTANCE_TYPE_SMALL
    elif t_val <= THRESH_LARGE_MAX:
        return INSTANCE_TYPE_LARGE
    elif t_val < THRESH_XLARGE_MAX:
        return INSTANCE_TYPE_XLARGE
    else:
        return INSTANCE_TYPE_HUGE

def _validate_limits(t_val: float, zeros_val: int):
    """Return error message if limits exceeded, else None."""
    if t_val < 0 or t_val > T_MAX_HARD:
        return f"t={t_val} is outside valid range (0 ≤ t ≤ {T_MAX_HARD})"
    if zeros_val < ZEROS_MIN_HARD or zeros_val > ZEROS_MAX_HARD:
        return f"zeros={zeros_val} is outside valid range ({ZEROS_MIN_HARD}–{ZEROS_MAX_HARD})"
    return None

def lambda_handler(event, context):
    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    image = body.get("image")
    t     = body.get("t")
    zeros = body.get("zeros")

    if not image or t is None or zeros is None:
        return {"statusCode": 400, "body": json.dumps({"error": "image, zeros and t are required"})}

    try:
        t_val = float(t)
        zeros_val = int(float(zeros))
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "t and zeros must be numeric"})}

    # --- Validate limits before doing anything else ---
    err = _validate_limits(t_val, zeros_val)
    if err:
        logger.warning(f"Rejected job: {err}")
        return {
            "statusCode": 400,
            "body": json.dumps({
                "error": err,
                "limits": {
                    "t_max": T_MAX_HARD,
                    "zeros_min": ZEROS_MIN_HARD,
                    "zeros_max": ZEROS_MAX_HARD
                }
            })
        }

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
                {"Key": "Name",    "Value": f"z-job-{job_id}"},
                {"Key": "z:jobId", "Value": job_id},
                {"Key": "z:t",     "Value": str(t_val)},
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
            "message": "Launched successfully. Scheduler active for JSON+S3 output."
        })
    }
