import os
import json
import uuid
import boto3
from string import Template
from decimal import Decimal
import logging

EC2 = boto3.client("ec2")
S3 = boto3.client("s3")

SUBNET_ID = os.environ["SUBNET_ID"]
SECURITY_GROUP_ID = os.environ["SECURITY_GROUP_ID"]
INSTANCE_PROFILE_NAME = os.environ["INSTANCE_PROFILE_NAME"]
S3_BUCKET = os.environ["S3_BUCKET"]
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

INSTANCE_TYPE_SMALL = os.environ.get("INSTANCE_TYPE_SMALL", "c7i.2xlarge")
INSTANCE_TYPE_LARGE = os.environ.get("INSTANCE_TYPE_LARGE", "c7i.2xlarge")
INSTANCE_TYPE_XLARGE = os.environ.get("INSTANCE_TYPE_XLARGE", "c7i.2xlarge")
INSTANCE_TYPE_HUGE = os.environ.get("INSTANCE_TYPE_HUGE", "c7i.2xlarge")

THRESH_SMALL_MAX = Decimal("1e8")
THRESH_LARGE_MAX = Decimal("1e12")
THRESH_XLARGE_MAX = Decimal("1e16")
THRESH_HUGE_MAX = Decimal("1e21")

USE_SPOT = os.environ.get("USE_SPOT", "false").lower() == "true"

T_MAX_HARD = Decimal("9999999999999999999999999")

logger = logging.getLogger()
logger.setLevel(logging.INFO)

for handler in logger.handlers:
    handler.setLevel(logging.INFO)


def load_user_data_from_s3(bucket: str, key: str, variables: dict) -> str:
    obj = S3.get_object(Bucket=bucket, Key=key)
    text = obj["Body"].read().decode("utf-8")
    tmpl = Template(text)
    return tmpl.safe_substitute(variables)

def pick_instance_type(t_val: Decimal) -> str:
    if t_val <= THRESH_SMALL_MAX:
        return INSTANCE_TYPE_SMALL
    elif t_val <= THRESH_LARGE_MAX:
        return INSTANCE_TYPE_LARGE
    elif t_val <= THRESH_XLARGE_MAX:
        return INSTANCE_TYPE_XLARGE
    else:
        return INSTANCE_TYPE_HUGE

def validate_t(t_val: Decimal):
    if t_val <= Decimal("14") or t_val > T_MAX_HARD:
        raise ValueError(f"t={t_val} is outside valid range (14 < t <= {T_MAX_HARD})")

def lambda_handler(event, context):
    body = event.get("body", event)
    if isinstance(body, str):
        body = json.loads(body)

    image = body.get("image")
    t = body.get("t")

    if not image or t is None:
        return {
            "statusCode": 400,
            "body": json.dumps({"error": "image and t are required"})
        }

    try:
        t_val = Decimal(str(t))
        validate_t(t_val)
    except Exception as e:
        return {
            "statusCode": 400,
            "body": json.dumps({"error": f"invalid t: {e}"})
        }

    job_id = str(uuid.uuid4())
    instance_type = pick_instance_type(t_val)

    vars_map = {
        "S3_BUCKET": S3_BUCKET,
        "job_id": job_id,
        "t_val": str(t_val),
    }

    template_key = "z-commands/precalc_commands"
    user_data = load_user_data_from_s3(S3_BUCKET, template_key, vars_map)

    params = {
        "ImageId": image,
        "InstanceType": instance_type,
        "MinCount": 1,
        "MaxCount": 1,
        "UserData": user_data,  # boto3 auto-encodes
        "IamInstanceProfile": {"Name": INSTANCE_PROFILE_NAME},
        "InstanceInitiatedShutdownBehavior": "terminate",
        "TagSpecifications": [{
            "ResourceType": "instance",
            "Tags": [
                {"Key": "Name", "Value": f"z-precalc-{job_id}"},
                {"Key": "z:jobId", "Value": job_id},
                {"Key": "z:t", "Value": str(t_val)},
                {"Key": "z:stage", "Value": "precalc"},
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
            "job_id": job_id,
            "instance_id": instance_id,
            "instance_type": instance_type,
            "s3_prefix": f"jobs/{job_id}/",
            "message": "Pre-calc instance launched"
        })
    }