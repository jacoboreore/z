# lambda_function.py
import os, json, base64, uuid, boto3

EC2 = boto3.client("ec2")

# ----- set env vars the Lambda expects -----
# os.environ["SUBNET_ID"]             = "subnet-084d809e9bb9891bc"
# os.environ["SECURITY_GROUP_ID"]     = "sg-0a4c92ac491102653"
# os.environ["INSTANCE_PROFILE_NAME"] = "z-batch-runner"
# os.environ["S3_BUCKET"]             = "jacob-z"
# os.environ["AWS_REGION"]            = "us-east-1"      # or your region

# # Optional autosize + spot + dry-run
# os.environ["INSTANCE_TYPE_SMALL"]   = "c7i.2xlarge"    # ~8 vCPU
# os.environ["INSTANCE_TYPE_LARGE"]   = "c7i.16xlarge"   # ~64 vCPU
# os.environ["T_THRESHOLD"]           = "1e14"
# os.environ["USE_SPOT"]              = "false"
# os.environ["DRY_RUN"]               = "true"           # first run: permission check only

# --- Required env vars ---
SUBNET_ID             = os.environ.get("SUBNET_ID")
SECURITY_GROUP_ID     = os.environ["SECURITY_GROUP_ID"]
INSTANCE_PROFILE_NAME = os.environ["INSTANCE_PROFILE_NAME"]     # role attached to EC2
S3_BUCKET             = os.environ["S3_BUCKET"]
AWS_REGION            = os.environ.get("AWS_REGION", "us-east-1")

# --- Optional autosize and options ---
INSTANCE_TYPE_SMALL   = os.environ.get("INSTANCE_TYPE_SMALL", "c7i.8xlarge")    # ~32 vCPU
INSTANCE_TYPE_LARGE   = os.environ.get("INSTANCE_TYPE_LARGE", "c7i.16xlarge")   # ~64 vCPU
T_THRESHOLD           = float(os.environ.get("T_THRESHOLD", "1e10"))
USE_SPOT              = os.environ.get("USE_SPOT", "false").lower() == "true"

def _pick_instance_type(t_val: float) -> str:
    return INSTANCE_TYPE_LARGE if t_val >= T_THRESHOLD else INSTANCE_TYPE_SMALL

def handler(event, context):
    # Expect {"image": "ami-...", "t": <number|string>}
    body = event.get("body")
    if isinstance(body, str):
        body = json.loads(body or "{}")
    elif body is None:
        body = {}

    image = body.get("image")
    t     = body.get("t")

    if not image or t is None:
        return {"statusCode": 400, "body": json.dumps({"error": "image and t are required"})}

    try:
        t_val = float(t)
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "t must be numeric"})}

    job_id = str(uuid.uuid4())
    instance_type = _pick_instance_type(t_val)

    # user-data: detect OS user (ubuntu/ec2-user), run binary, upload to S3, terminate
    user_data = f"""#!/bin/bash -xe
exec > /var/log/user-data.log 2>&1

JOB_ID="{job_id}"
T_VAL="{t_val}"
BUCKET="{S3_BUCKET}"
PREFIX="jobs/$JOB_ID"

# Detect default user
USR="ubuntu"
id ubuntu >/dev/null 2>&1 || USR="ec2-user"

# Ensure awscli + jq
# if ! command -v aws >/dev/null 2>&1; then
#   if command -v apt-get >/dev/null 2>&1; then
#     apt-get update -y || true
#     DEBIAN_FRONTEND=noninteractive apt-get install -y awscli jq || true
#   elif command -v yum >/dev/null 2>&1; then
#     yum install -y awscli jq || true
#   fi
# fi

# Workspace
mkdir -p /home/$USR/jacob/output || true
chown -R $USR:$USR /home/$USR || true

#######################################################################
##################### Run the app #####################################
#######################################################################

cd /home/$USR/jacob || exit 1
set +e
sudo -u $USR /home/$USR/jacob/bin/z valley_walk "$T_VAL" 10 $JOB_ID >&1 | tee -a /home/$USR/jacob/output/run.log
RC=${{PIPESTATUS[0]}}
set -e


#######################################################################
##################### S3 reports and work #############################
#######################################################################

# Ship artifacts
aws s3 cp /var/log/user-data.log s3://$BUCKET/$PREFIX/user-data.log || true
[ -d /home/$USR/jacob/output ] && aws s3 cp --recursive /home/$USR/jacob/output s3://$BUCKET/$PREFIX/output/ || true
[ -d /home/$USR/jacob/logs   ] && aws s3 cp --recursive /home/$USR/jacob/logs   s3://$BUCKET/$PREFIX/logs/   || true

# Marker + copy
if [ "$RC" -eq 0 ]; then
  echo "JOB_DONE"  | tee -a /home/$USR/jacob/output/run.log
else
  echo "JOB_FAILED (rc=$RC)" | tee -a /home/$USR/jacob/output/run.log
fi
aws s3 cp /home/$USR/jacob/output/run.log s3://$BUCKET/$PREFIX/run.log || true

# Terminate the instance on shutdown
shutdown -h now
"""

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
        # Attach network interface with public IP for easy troubleshooting
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
            "message": "Launched. Instance will run, upload results to S3, then terminate."
        })
    }
