import os

# ----- set env vars the Lambda expects -----
os.environ["SUBNET_ID"]             = "subnet-084d809e9bb9891bc"
os.environ["SECURITY_GROUP_ID"]     = "sg-0a4c92ac491102653"
os.environ["INSTANCE_PROFILE_NAME"] = "z-batch-runner"
os.environ["S3_BUCKET"]             = "jacob-z"
os.environ["AWS_REGION"]            = "us-east-1"      

# Optional autosize + spot + dry-run
os.environ["INSTANCE_TYPE_SMALL"]   = "c7i.2xlarge"    # ~8 vCPU
os.environ["INSTANCE_TYPE_LARGE"]   = "c7i.16xlarge"   # ~64 vCPU
os.environ["T_THRESHOLD"]           = "1e14"
os.environ["USE_SPOT"]              = "false"
os.environ["DRY_RUN"]               = "true"           # first run: permission check only