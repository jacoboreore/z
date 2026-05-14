import json
import os
import sys
import uuid
import datetime
import time
from decimal import Decimal
from typing import Any, Dict, List

import boto3
from botocore.exceptions import ClientError

AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
REFINE_JOBS_TABLE = os.environ.get("REFINE_JOBS_TABLE", "vs-refine-jobs")
REFINE_SAMPLES_TABLE = os.environ.get("REFINE_SAMPLES_TABLE", "vs-refine-evals")
DISPATCH_LAMBDA_NAME = os.environ.get("DISPATCH_LAMBDA_NAME", "z_dist_precalc")
AMI_IMAGE = os.environ.get("AMI_IMAGE", "ami-0fabf13310cc5147f")

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
lambda_client = boto3.client("lambda", region_name=AWS_REGION)

jobs_table = dynamodb.Table(REFINE_JOBS_TABLE)
samples_table = dynamodb.Table(REFINE_SAMPLES_TABLE)


def utc_now() -> str:
    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def make_eval_id(stage_name: str, label: str) -> str:
    return f"{stage_name}#{label}"


def load_json_file(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def validate_stage1_payload(data: Dict[str, Any]) -> None:
    required_top = ["action", "stage_name", "refine_id", "input", "metadata", "requested_points"]
    for key in required_top:
        if key not in data:
            raise ValueError(f"Missing required top-level field: {key}")

    if data["action"] != "queue_z_evaluations":
        raise ValueError(f"Unexpected action: {data['action']}")

    if not isinstance(data["requested_points"], list) or len(data["requested_points"]) == 0:
        raise ValueError("requested_points must be a non-empty list")

    declared_sample_count = int(data.get("metadata", {}).get("sample_count", 0))
    actual_sample_count = len(data["requested_points"])

    if declared_sample_count and declared_sample_count != actual_sample_count:
        raise ValueError(
            f"metadata.sample_count={declared_sample_count} does not match "
            f"requested_points length={actual_sample_count}"
        )

    labels_seen = set()
    for point in data["requested_points"]:
        if "label" not in point or "t" not in point:
            raise ValueError(f"Invalid point entry: {point}")
        if point["label"] in labels_seen:
            raise ValueError(f"Duplicate sample label: {point['label']}")
        labels_seen.add(point["label"])


def to_decimal_or_none(value: Any):
    if value is None:
        return None
    try:
        return Decimal(str(value))
    except Exception:
        return str(value)


def put_refine_job(data: Dict[str, Any]) -> None:
    now = int(time.time())

    item = {
        "refine_id": data["refine_id"],
        "stage_name": data["stage_name"],
        "status": "REQUESTED",
        "action": data.get("action", ""),
        "reason": data.get("reason", ""),
        "input_t0": str(data["input"].get("t0", "")),
        "mpfr_prec": int(data.get("config", {}).get("mpfr_prec", 0)),
        "rs_log_digits": int(data.get("config", {}).get("rs_log_digits", 0)),
        "h": str(data.get("metadata", {}).get("h", "")),
        "spacing_est": str(data.get("metadata", {}).get("spacing_est", "")),
        "sample_count": int(data.get("metadata", {}).get("sample_count", 0)),
        "requested_points": data["requested_points"],
        "created_at": now,
        "updated_at": now,
    }

    jobs_table.put_item(Item=item)


def put_sample_items(data: Dict[str, Any]) -> None:
    now = int(time.time())

    with samples_table.batch_writer() as batch:
        for point in data["requested_points"]:
            label = str(point["label"])
            t_value = str(point["t"])
            eval_id = make_eval_id(data["stage_name"], label)

            batch.put_item(
                Item={
                    "refine_id": data["refine_id"],
                    "eval_id": eval_id,
                    "label": label,
                    "stage_name": data["stage_name"],
                    "t": t_value,
                    "status": "PENDING",
                    "source": "refine_stage1_request",
                    "requested_at": now,
                    "updated_at": now,
                }
            )


def invoke_dispatch_lambda(refine_id: str, stage_name: str, point: Dict[str, Any]) -> Dict[str, Any]:
    payload = {
        "t": str(point["t"]),
        "image": AMI_IMAGE,
    }

    resp = lambda_client.invoke(
        FunctionName=DISPATCH_LAMBDA_NAME,
        InvocationType="RequestResponse",  # easier for testing/audit
        Payload=json.dumps(payload).encode("utf-8"),
    )

    response_payload = resp["Payload"].read().decode("utf-8") if "Payload" in resp else ""
    return {
        "status_code": resp.get("StatusCode"),
        "function_error": resp.get("FunctionError"),
        "response_payload": response_payload,
    }


def update_sample_dispatch_status(refine_id: str, sample_id: str, dispatch_result: Dict[str, Any]) -> None:
    now = int(time.time())

    samples_table.update_item(
        Key={"refine_id": refine_id, "eval_id": sample_id},
        UpdateExpression=(
            "SET #st = :st, dispatched_at = :dt, updated_at = :ut, "
            "lambda_status_code = :sc, lambda_function_error = :fe, lambda_response = :lr"
        ),
        ExpressionAttributeNames={"#st": "status"},
        ExpressionAttributeValues={
            ":st": "DISPATCHED" if not dispatch_result.get("function_error") else "DISPATCH_ERROR",
            ":dt": now,
            ":ut": now,
            ":sc": dispatch_result.get("status_code"),
            ":fe": dispatch_result.get("function_error"),
            ":lr": dispatch_result.get("response_payload", ""),
        },
    )


def update_job_dispatched(refine_id: str) -> None:
    now = int(time.time())

    jobs_table.update_item(
        Key={"refine_id": refine_id},
        UpdateExpression="SET #st = :st, updated_at = :ut",
        ExpressionAttributeNames={"#st": "status"},
        ExpressionAttributeValues={
            ":st": "DISPATCHED",
            ":ut": now,
        },
    )


def update_job_error(refine_id: str, error_message: str) -> None:
    now = int(time.time())
    
    jobs_table.update_item(
        Key={"refine_id": refine_id},
        UpdateExpression="SET #st = :st, updated_at = :ut, error_message = :em",
        ExpressionAttributeNames={"#st": "status"},
        ExpressionAttributeValues={
            ":st": "ERROR",
            ":ut": now,
            ":em": error_message[:4000],
        },
    )


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <stage1_json>")
        return 1

    json_path = sys.argv[1]

    try:
        data = load_json_file(json_path)
        validate_stage1_payload(data)

        print(f"[INFO] Loaded stage1 JSON for refine_id={data['refine_id']}")
        put_refine_job(data)
        print("[INFO] Refinement job written to DynamoDB")

        put_sample_items(data)
        print("[INFO] Sample items written to DynamoDB")

        for point in data["requested_points"]:
            label = str(point["label"])
            eval_id = make_eval_id(data["stage_name"], label)

            print(f"[INFO] Dispatching sample {label} t={point['t']}")
            dispatch_result = invoke_dispatch_lambda(
                refine_id=data["refine_id"],
                stage_name=data["stage_name"],
                point=point,
            )
            update_sample_dispatch_status(data["refine_id"], eval_id, dispatch_result)
            print(
                f"[INFO] Sample {label} dispatch status="
                f"{dispatch_result.get('status_code')} "
                f"function_error={dispatch_result.get('function_error')}"
            )

        update_job_dispatched(data["refine_id"])
        print(f"[INFO] Job {data['refine_id']} marked DISPATCHED")
        return 0

    except Exception as e:
        refine_id = None
        try:
            maybe_data = load_json_file(json_path)
            refine_id = maybe_data.get("refine_id")
        except Exception:
            pass

        print(f"[ERROR] {e}", file=sys.stderr)
        if refine_id:
            try:
                update_job_error(refine_id, str(e))
            except Exception as inner:
                print(f"[WARN] Failed to update Dynamo error state: {inner}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())