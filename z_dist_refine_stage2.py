import json
import os
import sys
import time
import subprocess
import datetime
from decimal import Decimal
from typing import Any, Dict, List, Optional, Tuple

import boto3
from boto3.dynamodb.conditions import Key


AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
REFINE_JOBS_TABLE = os.environ.get("REFINE_JOBS_TABLE", "vs-refine-jobs")
REFINE_SAMPLES_TABLE = os.environ.get("REFINE_SAMPLES_TABLE", "vs-refine-evals")
Z_JOBS_TABLE = os.environ.get("Z_JOBS_TABLE", "valley-scanner-dev-jobs")

POLL_SECONDS = int(os.environ.get("POLL_SECONDS", "10"))
MAX_WAIT_SECONDS = int(os.environ.get("MAX_WAIT_SECONDS", "3600"))

# Next-stage handoff
STAGE2_JSON_DIR = os.environ.get("STAGE2_JSON_DIR", "output")
STAGE2_CPP_BIN = os.environ.get("STAGE2_CPP_BIN", "./bin/z")

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
refine_jobs_table = dynamodb.Table(REFINE_JOBS_TABLE)
samples_table = dynamodb.Table(REFINE_SAMPLES_TABLE)
z_jobs_table = dynamodb.Table(Z_JOBS_TABLE)

def make_json_safe(obj):
    if isinstance(obj, Decimal):
        # Preserve integers as int when possible
        if obj == obj.to_integral_value():
            return int(obj)
        return float(obj)

    if isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}

    if isinstance(obj, list):
        return [make_json_safe(v) for v in obj]

    return obj

def utc_now() -> str:
    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def log(msg: str) -> None:
    print(f"[{utc_now()}] [INFO] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[{utc_now()}] [WARN] {msg}", flush=True)


def err(msg: str) -> None:
    print(f"[{utc_now()}] [ERROR] {msg}", file=sys.stderr, flush=True)


def get_refine_job(refine_id: str) -> Optional[Dict[str, Any]]:
    resp = refine_jobs_table.get_item(Key={"refine_id": refine_id})
    return resp.get("Item")


def get_eval_rows(refine_id: str) -> List[Dict[str, Any]]:
    resp = samples_table.query(
        KeyConditionExpression=Key("refine_id").eq(refine_id)
    )
    return resp.get("Items", [])


def get_z_job(job_id: str) -> Optional[Dict[str, Any]]:
    resp = z_jobs_table.get_item(Key={"job_id": job_id})
    return resp.get("Item")


def update_job_status(refine_id: str, status: str, extra: Optional[Dict[str, Any]] = None) -> None:
    extra = extra or {}
    now = int(time.time())

    update_expr = "SET #st = :st, updated_at = :ut"
    expr_names = {"#st": "status"}
    expr_vals = {
        ":st": status,
        ":ut": now,
    }

    for i, (k, v) in enumerate(extra.items(), start=1):
        name_key = f"#k{i}"
        val_key = f":v{i}"
        update_expr += f", {name_key} = {val_key}"
        expr_names[name_key] = k
        expr_vals[val_key] = v

    refine_jobs_table.update_item(
        Key={"refine_id": refine_id},
        UpdateExpression=update_expr,
        ExpressionAttributeNames=expr_names,
        ExpressionAttributeValues=expr_vals,
    )


def update_eval_row(refine_id: str, eval_id: str, extra: Dict[str, Any]) -> None:
    now = now = int(time.time())

    update_expr = "SET updated_at = :ut"
    expr_names = {}
    expr_vals = {":ut": now}

    for i, (k, v) in enumerate(extra.items(), start=1):
        name_key = f"#k{i}"
        val_key = f":v{i}"
        update_expr += f", {name_key} = {val_key}"
        expr_names[name_key] = k
        expr_vals[val_key] = v

    samples_table.update_item(
        Key={"refine_id": refine_id, "eval_id": eval_id},
        UpdateExpression=update_expr,
        ExpressionAttributeNames=expr_names if expr_names else None,
        ExpressionAttributeValues=expr_vals,
    )


def launch_successful(item: Dict[str, Any]) -> bool:
    status = item.get("status", "")
    lambda_function_error = item.get("lambda_function_error")
    return status != "DISPATCH_ERROR" and not lambda_function_error


def eval_failed(item: Dict[str, Any]) -> bool:
    return item.get("status") in {"ERROR", "FAILED", "DISPATCH_ERROR"}


def parse_lambda_response_for_job_id(raw: Any) -> str:
    """
    Expected shape from your sample:
    lambda_response = "\"{""statusCode"": 200, ""body"": ""{\\""job_id\\"": \\""uuid\\""...}""}\""

    We decode:
    1) outer JSON string/object
    2) outer.body if present
    3) inner JSON body
    4) inner['job_id']
    """
    if raw is None:
        raise ValueError("lambda_response is missing")

    text = raw if isinstance(raw, str) else json.dumps(raw)

    # First decode
    first = json.loads(text)

    # Sometimes first decode returns stringified JSON again
    if isinstance(first, str):
        first = json.loads(first)

    if not isinstance(first, dict):
        raise ValueError(f"Unexpected lambda_response shape after first decode: {type(first)}")

    body = first.get("body")
    if body is None:
        raise ValueError("lambda_response does not contain body")

    if isinstance(body, str):
        inner = json.loads(body)
    elif isinstance(body, dict):
        inner = body
    else:
        raise ValueError(f"Unexpected body type in lambda_response: {type(body)}")

    job_id = inner.get("job_id")
    if not job_id:
        raise ValueError("Could not find job_id in lambda_response body")

    return str(job_id)


def ensure_z_job_ids(refine_id: str, eval_rows: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Returns mapping eval_id -> z_job_id.
    If parse fails for any eval, abort.
    Also persists z_job_id back into vs-refine-evals for convenience.
    """
    mapping: Dict[str, str] = {}

    for row in eval_rows:
        if "label" not in row or not row["label"]:
            raise RuntimeError(f"Eval row missing label: {row}")
        if "eval_id" not in row or not row["eval_id"]:
            raise RuntimeError(f"Eval row missing eval_id: {row}")

        eval_id = row["eval_id"]
        label = row.get("label", "")

        if eval_failed(row):
            raise RuntimeError(f"Eval row already failed: {eval_id}")

        if not launch_successful(row):
            raise RuntimeError(f"Eval row not launch-successful yet: eval_id={eval_id} label={label}")

        z_job_id = row.get("z_job_id")
        if not z_job_id:
            z_job_id = parse_lambda_response_for_job_id(row.get("lambda_response"))
            update_eval_row(refine_id, eval_id, {"z_job_id": z_job_id})

        mapping[eval_id] = z_job_id

    return mapping


def z_job_done(item: Dict[str, Any]) -> bool:
    return item.get("status") == "DONE" and item.get("ztotal") not in (None, "")


def z_job_failed(item: Dict[str, Any]) -> bool:
    return item.get("status") in {"ERROR", "FAILED", "ABORTED"}


def normalize_z_string(value: Any) -> str:
    """
    Your sample ztotal values appear like:
      "'-1.18361511...."
    So strip whitespace and a leading apostrophe if present.
    """
    s = str(value).strip()
    if s.startswith("'"):
        s = s[1:]
    return s


def build_stage2_payload(refine_job: Dict[str, Any],
                         eval_rows: List[Dict[str, Any]],
                         z_jobs_by_eval: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:

    for row in eval_rows:
        if "label" not in row or not row["label"]:
            raise RuntimeError(f"Eval row missing label: {row}")
        if "eval_id" not in row or not row["eval_id"]:
            raise RuntimeError(f"Eval row missing eval_id: {row}")

    rows_by_label = {row["label"]: row for row in eval_rows}

    sample_count = int(refine_job.get("sample_count", 0))
    if sample_count <= 0:
        sample_count = len(eval_rows)

    # center is always mandatory
    if "center" not in rows_by_label:
        raise RuntimeError("Missing required label: center")

    # support both 3-point and 5-point layouts
    if {"left2", "left1", "center", "right1", "right2"}.issubset(rows_by_label.keys()):
        ordered_labels = ["left2", "left1", "center", "right1", "right2"]
    elif {"left1", "center", "right1"}.issubset(rows_by_label.keys()):
        ordered_labels = ["left1", "center", "right1"]
    else:
        raise RuntimeError(
            f"Unsupported label set for refinement payload: {sorted(rows_by_label.keys())}"
        )

    if len(ordered_labels) != sample_count:
        raise RuntimeError(
            f"Refine job sample_count={sample_count} but resolved labels={ordered_labels}"
        )

    def eval_entry(label: str) -> Dict[str, Any]:
        row = rows_by_label[label]
        eval_id = row["eval_id"]

        if eval_id not in z_jobs_by_eval:
            raise RuntimeError(f"Missing downstream Z job for eval_id={eval_id} label={label}")

        zjob = z_jobs_by_eval[eval_id]

        return {
            "label": label,
            "t": str(row.get("t", "")),
            "z_value": normalize_z_string(zjob.get("ztotal", "")),
            "status": "COMPLETED",
            "z_job_id": zjob.get("job_id", ""),
            "z_job_status": zjob.get("status", "")
        }

    return {
        "action": "decide_next_refinement_step",
        "stage_name": "refine_stage_2",
        "reason": "stage1_samples_completed",
        "refine_id": refine_job["refine_id"],
        "input": {
            "t0": refine_job.get("input_t0", "")
        },
        "config": {
            "mpfr_prec": refine_job.get("mpfr_prec", 0),
            "rs_log_digits": refine_job.get("rs_log_digits", 0)
        },
        "metadata": {
            "h": refine_job.get("h", ""),
            "spacing_est": refine_job.get("spacing_est", ""),
            "sample_count": sample_count,
            "generated_at": utc_now()
        },
        "evaluations": [eval_entry(label) for label in ordered_labels]
    }


def write_stage2_input_json(refine_id: str, payload: Dict[str, Any]) -> str:
    os.makedirs(STAGE2_JSON_DIR, exist_ok=True)
    path = os.path.join(STAGE2_JSON_DIR, f"{refine_id}_stage2_input.json")

    safe_payload = make_json_safe(payload)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(safe_payload, f, indent=2)

    return path


def stage2_output_path(refine_id: str) -> str:
    os.makedirs(STAGE2_JSON_DIR, exist_ok=True)
    return os.path.join(STAGE2_JSON_DIR, f"{refine_id}_stage2_output.json")


def run_stage2_cpp(input_json_path: str, output_json_path: str) -> int:
    cmd = [STAGE2_CPP_BIN, "z_distributed_refine_stage2", input_json_path, output_json_path, "refinement_job"]
    log(f"Launching C++ stage2: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout, end="", flush=True)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr, flush=True)

    return result.returncode


def main() -> int:
    if len(sys.argv) != 2:
        err(f"Usage: {sys.argv[0]} <refine_id>")
        return 1

    refine_id = sys.argv[1]
    log(f"Starting monitor for refine_id={refine_id}")

    refine_job = get_refine_job(refine_id)
    if not refine_job:
        err(f"Refine job not found: refine_id={refine_id}")
        return 2

    log(f"Loaded refine job with status={refine_job.get('status')}")
    update_job_status(refine_id, "WAITING_EVALS")

    start = time.time()

    while True:
        elapsed = int(time.time() - start)
        if elapsed > MAX_WAIT_SECONDS:
            err(f"Timeout waiting for eval completion refine_id={refine_id}")
            update_job_status(refine_id, "TIMEOUT", {"monitor_error": "Timed out waiting for completed Z jobs"})
            return 3

        eval_rows = get_eval_rows(refine_id)
        log(f"Fetched {len(eval_rows)} eval rows for refine_id={refine_id}")

        expected_samples = int(refine_job.get("sample_count", 0) or 0)
        if expected_samples <= 0:
            expected_samples = 3

        if len(eval_rows) < expected_samples:
            warn(f"Expected {expected_samples} eval rows, got {len(eval_rows)}. Waiting...")
            time.sleep(POLL_SECONDS)
            continue

        failed = [row for row in eval_rows if eval_failed(row)]
        if failed:
            err(f"Detected failed eval rows for refine_id={refine_id}: {[r.get('eval_id') for r in failed]}")
            update_job_status(
                refine_id,
                "ERROR",
                {"monitor_error": f"Failed evals: {[r.get('eval_id') for r in failed]}"}
            )
            return 4

        try:
            z_job_ids = ensure_z_job_ids(refine_id, eval_rows)
        except Exception as e:
            err(f"Failed extracting downstream z job ids for refine_id={refine_id}: {e}")
            update_job_status(refine_id, "ERROR", {"monitor_error": f"z_job_id parse failure: {str(e)}"})
            return 5

        log(f"Resolved z_job_ids for refine_id={refine_id}: {z_job_ids}")

        z_jobs_by_eval: Dict[str, Dict[str, Any]] = {}
        pending = []
        failed_z = []

        for eval_id, job_id in z_job_ids.items():
            zjob = get_z_job(job_id)
            if not zjob:
                warn(f"Z job not found yet for eval_id={eval_id} job_id={job_id}")
                pending.append(eval_id)
                continue

            z_jobs_by_eval[eval_id] = zjob
            status = zjob.get("status")

            log(f"Z job eval_id={eval_id} job_id={job_id} status={status}")

            if z_job_failed(zjob):
                failed_z.append((eval_id, job_id, status))
            elif not z_job_done(zjob):
                pending.append(eval_id)

        if failed_z:
            err(f"Detected failed Z jobs for refine_id={refine_id}: {failed_z}")
            update_job_status(
                refine_id,
                "ERROR",
                {"monitor_error": f"Failed Z jobs: {failed_z}"}
            )
            return 6

        if pending:
            log(f"Still waiting for Z jobs refine_id={refine_id}: {pending}")
            time.sleep(POLL_SECONDS)
            continue

        # At this point all 3 downstream Z jobs are DONE
        log(f"All downstream Z jobs completed for refine_id={refine_id}")
        update_job_status(refine_id, "STAGE1_COMPLETED")

        stage2_payload = build_stage2_payload(refine_job, eval_rows, z_jobs_by_eval)
        input_json_path = write_stage2_input_json(refine_id, stage2_payload)
        output_json_path = stage2_output_path(refine_id)

        log(f"Stage2 input JSON written to {input_json_path}")
        log(f"Stage2 output JSON path will be {output_json_path}")

        update_job_status(
            refine_id,
            "STAGE2_READY",
            {
                "stage2_input_json_path": input_json_path,
                "stage2_output_json_path": output_json_path
            }
        )

        rc = run_stage2_cpp(input_json_path, output_json_path)
        if rc != 0:
            err(f"C++ stage2 returned rc={rc} for refine_id={refine_id}")
            update_job_status(refine_id, "ERROR", {"monitor_error": f"Stage2 C++ failed rc={rc}"})
            return 7

        log(f"C++ stage2 completed successfully for refine_id={refine_id}")
        update_job_status(refine_id, "STAGE2_LAUNCHED")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
