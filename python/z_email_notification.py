import os
import json
import boto3
from datetime import datetime, timezone

sns = boto3.client("sns")
TOPIC_ARN = "arn:aws:sns:us-east-1:407875882825:z-mail-alert-topic-jacob"

def lambda_handler(event, context):
    message = {
        "alert": "Valley Scanner calculation triggered",
        "time_utc": datetime.now(timezone.utc).isoformat(),
        "event": event
    }

    sns.publish(
        TopicArn=TOPIC_ARN,
        Subject="VS calculation triggered",
        Message=json.dumps(message, indent=2, default=str)
    )

    return {
        "ok": True,
        "message": "Email notification sent"
    }
