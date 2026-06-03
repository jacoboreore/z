import json
import os
import re
import logging
from typing import Dict, Any

import boto3
from botocore.exceptions import ClientError

logger = logging.getLogger()
logger.setLevel(logging.INFO)

REGION = os.getenv("AWS_REGION", os.getenv("AWS_DEFAULT_REGION", "us-east-1"))

# Env vars
USER_POOL_ID = os.environ["USER_POOL_ID"]
APP_CLIENT_ID = os.environ["APP_CLIENT_ID"]  # App client WITHOUT secret; ADMIN_NO_SRP_AUTH enabled
ARTIFACTS_BUCKET = os.environ["ARTIFACTS_BUCKET"]
DEFAULT_ARTIFACT_KEY = os.getenv("DEFAULT_ARTIFACT_KEY", "sample/report.txt")

# Optional JSON allowlist in env; else fall back to demo values
try:
    ALLOWLIST: Dict[str, str] = json.loads(os.getenv("ALLOWLIST_JSON", "{}"))
except json.JSONDecodeError:
    ALLOWLIST = {}

if not ALLOWLIST:
    logger.error(f"Allowed users not set")

# uncomment for testing
#     ALLOWLIST = {
#         "alice@example.com": "Alice-Temp#123",
#         "bob@example.com": "Bob-Temp#456",
#         "carol@example.com": "Carol-Temp#789",
#     }

cognito = boto3.client("cognito-idp", region_name=REGION)
S3 = boto3.client("s3", region_name=REGION)

def _json(status: int, body: Dict[str, Any], cors: bool = True):
    headers = {"Content-Type": "application/json"}
    if cors:
        headers.update({
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Headers": "Content-Type,Authorization",
            "Access-Control-Allow-Methods": "GET,POST,OPTIONS",
        })
    return {"statusCode": status, "headers": headers, "body": json.dumps(body)}

def _ensure_user(email: str):
    try:
        cognito.admin_get_user(UserPoolId=USER_POOL_ID, Username=email)
        return
    except ClientError as e:
        if e.response["Error"]["Code"] != "UserNotFoundException":
            raise
    
    # Create user suppressed (no email)
    cognito.admin_create_user(
        UserPoolId=USER_POOL_ID,
        Username=email,
        MessageAction="SUPPRESS",
        UserAttributes=[
            {"Name": "email", "Value": email},
            {"Name": "email_verified", "Value": "true"},
        ],
    )

def _admin_auth(username: str, password: str) -> Dict[str, Any]:
    resp = cognito.admin_initiate_auth(
        UserPoolId=USER_POOL_ID,
        ClientId=APP_CLIENT_ID,
        AuthFlow="ADMIN_NO_SRP_AUTH",
        AuthParameters={"USERNAME": username, "PASSWORD": password},
    )
    
    return resp.get("AuthenticationResult", {})

def _bearer(event) -> str | None:
    headers = event.get("headers") or {}
    token = headers.get("Authorization") or headers.get("authorization") or ""
    m = re.match(r"^Bearer\s+(.+)$", token)
    return m.group(1) if m else None

def lambda_handler(event, context):
    method = event.get("httpMethod", "GET")
    path = event.get("path", "/")

    logger.info(f"PATH EVENT {path}")

    if method == "OPTIONS":
        return _json(200, {"ok": True})

    try:
        try:
            payload = json.loads(event.get("body") or "{}")
        except json.JSONDecodeError:
            return _json(400, {"error": "invalid_json"})

        email = (payload.get("email") or "").strip()
        password = payload.get("password") or ""
        
        if not email or not password:
            return _json(400, {"error": "email and password required"})

        expected = ALLOWLIST.get(email)
        if expected is None or expected != password:
            return _json(401, {"error": "invalid credentials"})

        _ensure_user(email)
        # Set/Reset password permanent
        cognito.admin_set_user_password(
            UserPoolId=USER_POOL_ID, Username=email, Password=password, Permanent=True
        )

        tokens = _admin_auth(email, password)
        return _json(200, {"email": email, "tokens": tokens})
        
    except ClientError as e:
        logger.exception("AWS error")
        return _json(500, {"error": e.response.get("Error", {}).get("Code", "aws_error"), "detail": e.response.get("Error", {}).get("Message")})
    except Exception as e:
        logger.exception("Unhandled error")
        return _json(500, {"error": "server_error", "detail": str(e)})

