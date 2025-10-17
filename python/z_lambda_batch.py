# lambda_function.py
import os, json, base64, uuid, boto3

EC2 = boto3.client("ec2")

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

def lambda_handler(event, context):
    # Expect {"image": "ami-...", "t": <number|string>}
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
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "t must be numeric"})}

    try:
        zeros_val = float(zeros)
    except Exception:
        return {"statusCode": 400, "body": json.dumps({"error": "zeros must be numeric"})}

    job_id = str(uuid.uuid4())
    instance_type = _pick_instance_type(t_val)

    # user-data: detect OS user (ubuntu/ec2-user), run binary, upload to S3, terminate
    user_data = f"""#!/bin/bash -xe
exec > /var/log/user-data.log 2>&1

JOB_ID="{job_id}"
T_VAL="{t_val}"
ZEROS_VAL="{zeros_val}"
BUCKET="{S3_BUCKET}"
PREFIX="jobs/$JOB_ID"

# Detect default user
USR="ubuntu"
id ubuntu >/dev/null 2>&1 || USR="ec2-user"

# AWSCLI IS ALREADY INSTALLED IN THE INSTANCE, THEN AN IMAGE WAS GENERATED

# Workspace
mkdir -p /home/$USR/jacob/output || true
chown -R $USR:$USR /home/$USR || true

#######################################################################
##################### Run the app #####################################
#######################################################################

cd /home/$USR/jacob || exit 1
set +e
sudo -u $USR /home/$USR/jacob/bin/z valley_walk "$T_VAL" "$ZEROS_VAL" $JOB_ID >&1 | tee -a /home/$USR/jacob/output/run.log
RC=${{PIPESTATUS[0]}}
set -e

#######################################################################
######################## Plotting #####################################
#######################################################################

latest=$(find output -maxdepth 1 -type f -name '*walk_plot*.csv' -printf '%T@ %p\n' | sort -nr | head -n1 | cut -d' ' -f2-)
[ -n "$latest" ] && python3 plotting_mountains.py "$latest" || echo "No CSVs in output/"

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
