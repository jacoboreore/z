#!/usr/bin/env python3
"""
z_dist_refine_loop.py

Distributed refinement loop controller for Valley Scanner.

Design rules
------------
- C++ owns all math and all refinement decisions.
- Python only orchestrates:
  * read C++ output JSON
  * persist job/sample state
  * dispatch requested Z evaluations
  * monitor DynamoDB until downstream Z jobs finish
  * build next C++ input JSON
  * invoke C++ again
  * repeat until terminal action

Typical flow
------------
1) Initial C++ stage writes an output JSON with:
      action=queue_z_evaluations
      refine_id
      stage_name
      requested_points

2) This script consumes that JSON and dispatches the requested Z(t) evaluations.

3) It waits until all downstream Z jobs are DONE.

4) It builds the input JSON expected by the C++ decision engine and runs:
      ./bin/z z_distributed_refine_stage2 <input_json> <output_json> refinement_job

5) It reads the new C++ output JSON:
      - if action == queue_z_evaluations: loop again
      - otherwise: stop and mark terminal state

Notes
-----
- This script does NOT decide tolerance, Newton acceptance, bracketing, etc.
- It supports generic requested_points emitted by C++:
    * 1 point (e.g. mid_bracket)
    * 3 points (e.g. left1/center/right1)
    * wider scans, or other future labels
"""

from __future__ import annotations

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
from botocore.exceptions import ClientError


# -----------------------------------------------------------------------------
# Environment / configuration
# -----------------------------------------------------------------------------

AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")

REFINE_JOBS_TABLE = os.environ.get("REFINE_JOBS_TABLE", "vs-refine-jobs")
REFINE_SAMPLES_TABLE = os.environ.get("REFINE_SAMPLES_TABLE", "vs-refine-evals")
Z_JOBS_TABLE = os.environ.get("Z_JOBS_TABLE", "valley-scanner-dev-jobs")

DISPATCH_LAMBDA_NAME = os.environ.get("DISPATCH_LAMBDA_NAME", "z_dist_precalc")
AMI_IMAGE = os.environ.get("AMI_IMAGE", "ami-0fabf13310cc5147f")

POLL_SECONDS = int(os.environ.get("POLL_SECONDS", "10"))
MAX_WAIT_SECONDS = int(os.environ.get("MAX_WAIT_SECONDS", "3600"))
MAX_ITERATIONS = int(os.environ.get("MAX_ITERATIONS", "8"))

REFINE_JSON_DIR = os.environ.get("REFINE_JSON_DIR", "output")
STAGE2_CPP_BIN = os.environ.get("STAGE2_CPP_BIN", "./bin/z")

# Terminal actions returned by C++.
# Extend this set later if the refiner grows more terminal outcomes.
TERMINAL_ACTIONS = {
    "accept_root",
    "wait",   # optional / non-final in theory, but here treated as no-more-dispatch
    "error",
    "abort",
    "failed",
}

dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)
lambda_client = boto3.client("lambda", region_name=AWS_REGION)

refine_jobs_table = dynamodb.Table(REFINE_JOBS_TABLE)
samples_table = dynamodb.Table(REFINE_SAMPLES_TABLE)
z_jobs_table = dynamodb.Table(Z_JOBS_TABLE)


# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def log(msg: str) -> None:
    print(f"[{utc_now()}] [INFO] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[{utc_now()}] [WARN] {msg}", flush=True)


def err(msg: str) -> None:
    print(f"[{utc_now()}] [ERROR] {msg}", file=sys.stderr, flush=True)


def banner(title: str) -> None:
    line = "=" * 78
    print(f"\n{line}\n[{utc_now()}] {title}\n{line}", flush=True)

def print_best_z_panel(current_payload: Dict[str, Any], z_jobs_by_eval: Dict[str, Dict[str, Any]]) -> None:
    best = None

    for point in current_payload.get("requested_points", []):
        label = str(point.get("label", ""))
        eval_id = make_eval_id(current_payload["stage_name"], label)
        zjob = z_jobs_by_eval.get(eval_id)

        if not zjob:
            continue

        z_str = normalize_z_string(zjob.get("ztotal", ""))
        t_str = str(point.get("t", ""))

        try:
            abs_z = abs(float(z_str))
        except Exception:
            continue

        if best is None or abs_z < best[0]:
            best = (abs_z, z_str, t_str, label)

    if best:
        print("\n++++++++++++++++ JACOB BEST ++++++++++++++++++++", flush=True)
        print("\n" + "-" * 70, flush=True)
        print(f"BEST Z SO FAR | abs_z={best[0]:.17g} | z={best[1]}", flush=True)
        print(f"t={best[2]}", flush=True)
        print(f"label={best[3]}", flush=True)
        print("-" * 70 + "\n", flush=True)
        
# -----------------------------------------------------------------------------
# JSON helpers
# -----------------------------------------------------------------------------

def make_json_safe(obj: Any) -> Any:
    """
    Convert Dynamo Decimal values and nested structures into JSON-safe objects.

    We preserve integers as int when exact, otherwise use float for readability.
    This is orchestration metadata only; the actual high-precision numeric strings
    for t and Z values are kept as strings from the C++ side.
    """
    if isinstance(obj, Decimal):
        if obj == obj.to_integral_value():
            return int(obj)
        return float(obj)

    if isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}

    if isinstance(obj, list):
        return [make_json_safe(v) for v in obj]

    return obj


def load_json_file(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json_file(path: str, payload: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(make_json_safe(payload), f, indent=2)
        f.write("\n")


# -----------------------------------------------------------------------------
# Dynamo helpers
# -----------------------------------------------------------------------------

def get_refine_job(refine_id: str) -> Optional[Dict[str, Any]]:
    resp = refine_jobs_table.get_item(Key={"refine_id": refine_id})
    return resp.get("Item")


def get_eval_rows_for_stage(refine_id: str, stage_name: str) -> List[Dict[str, Any]]:
    """
    Return all eval rows for a refine_id that belong to one logical stage iteration.

    The table partition key is refine_id, so we query all rows then filter locally.
    This is acceptable here because per-refinement row count is modest.
    """
    resp = samples_table.query(KeyConditionExpression=Key("refine_id").eq(refine_id))
    items = resp.get("Items", [])
    return [row for row in items if row.get("stage_name") == stage_name]


def get_z_job(job_id: str) -> Optional[Dict[str, Any]]:
    resp = z_jobs_table.get_item(Key={"job_id": job_id})
    return resp.get("Item")


def update_job_status(refine_id: str, status: str, extra: Optional[Dict[str, Any]] = None) -> None:
    """
    Update top-level status for the refinement job.
    """
    extra = extra or {}
    now = int(time.time())

    update_expr = "SET #st = :st, updated_at = :ut"
    expr_names = {"#st": "status"}
    expr_vals = {
        ":st": status,
        ":ut": now,
    }

    for i, (k, v) in enumerate(extra.items(), start=1):
        nk = f"#k{i}"
        vk = f":v{i}"
        update_expr += f", {nk} = {vk}"
        expr_names[nk] = k
        expr_vals[vk] = v

    refine_jobs_table.update_item(
        Key={"refine_id": refine_id},
        UpdateExpression=update_expr,
        ExpressionAttributeNames=expr_names,
        ExpressionAttributeValues=expr_vals,
    )


def put_or_update_refine_job(initial_payload: Dict[str, Any]) -> None:
    """
    Create or refresh the top-level refinement job record.

    This is called at loop start so a fresh run can resume or initialize the job.
    """
    now = int(time.time())

    refine_id = initial_payload["refine_id"]

    item = {
        "refine_id": refine_id,
        "stage_name": initial_payload.get("stage_name", ""),
        "status": "LOOP_STARTED",
        "action": initial_payload.get("action", ""),
        "reason": initial_payload.get("reason", ""),
        "input_t0": str(initial_payload.get("input", {}).get("t0", "")),
        "requested_points": initial_payload.get("requested_points", []),
        "created_at": now,
        "updated_at": now,
    }

    existing = get_refine_job(refine_id)
    if existing:
        update_job_status(
            refine_id,
            "LOOP_STARTED",
            {
                "stage_name": initial_payload.get("stage_name", ""),
                "action": initial_payload.get("action", ""),
                "reason": initial_payload.get("reason", ""),
                "requested_points": initial_payload.get("requested_points", []),
            },
        )
    else:
        refine_jobs_table.put_item(Item=item)


def update_eval_row(refine_id: str, eval_id: str, extra: Dict[str, Any]) -> None:
    """
    Patch one eval row in Dynamo.
    """
    now = int(time.time())

    update_expr = "SET updated_at = :ut"
    expr_names = {}
    expr_vals = {":ut": now}

    for i, (k, v) in enumerate(extra.items(), start=1):
        nk = f"#k{i}"
        vk = f":v{i}"
        update_expr += f", {nk} = {vk}"
        expr_names[nk] = k
        expr_vals[vk] = v

    kwargs = {
        "Key": {"refine_id": refine_id, "eval_id": eval_id},
        "UpdateExpression": update_expr,
        "ExpressionAttributeValues": expr_vals,
    }
    if expr_names:
        kwargs["ExpressionAttributeNames"] = expr_names

    samples_table.update_item(**kwargs)


def put_sample_items(refine_id: str, stage_name: str, requested_points: List[Dict[str, Any]]) -> None:
    """
    Persist the requested points for one stage iteration.

    eval_id is stage_name#label so labels remain readable and stage-local.
    """
    now = int(time.time())

    with samples_table.batch_writer() as batch:
        for point in requested_points:
            label = str(point["label"])
            t_value = str(point["t"])
            eval_id = make_eval_id(stage_name, label)

            batch.put_item(
                Item={
                    "refine_id": refine_id,
                    "eval_id": eval_id,
                    "label": label,
                    "stage_name": stage_name,
                    "t": t_value,
                    "status": "PENDING",
                    "source": "refine_loop_request",
                    "requested_at": now,
                    "updated_at": now,
                }
            )


# -----------------------------------------------------------------------------
# Small orchestration helpers
# -----------------------------------------------------------------------------

def make_eval_id(stage_name: str, label: str) -> str:
    return f"{stage_name}#{label}"


def normalize_z_string(value: Any) -> str:
    """
    Downstream ztotal sometimes arrives like:
        "'-1.18361511..."
    Strip whitespace and a leading apostrophe if present.
    """
    s = str(value).strip()
    if s.startswith("'"):
        s = s[1:]
    return s


def parse_lambda_response_for_job_id(raw: Any) -> str:
    """
    Decode the lambda response payload and extract job_id.

    Expected outer shapes seen so far:
      - JSON string
      - JSON object with body field
      - nested JSON string in body
    """
    if raw is None:
        raise ValueError("lambda_response is missing")

    text = raw if isinstance(raw, str) else json.dumps(raw)
    first = json.loads(text)

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


def launch_successful(item: Dict[str, Any]) -> bool:
    status = item.get("status", "")
    lambda_function_error = item.get("lambda_function_error")
    return status != "DISPATCH_ERROR" and not lambda_function_error


def eval_failed(item: Dict[str, Any]) -> bool:
    return item.get("status") in {"ERROR", "FAILED", "DISPATCH_ERROR"}


def z_job_done(item: Dict[str, Any]) -> bool:
    return item.get("status") == "DONE" and item.get("ztotal") not in (None, "")


def z_job_failed(item: Dict[str, Any]) -> bool:
    return item.get("status") in {"ERROR", "FAILED", "ABORTED"}


# -----------------------------------------------------------------------------
# Validation of C++ output
# -----------------------------------------------------------------------------

def validate_queue_payload(data: Dict[str, Any]) -> None:
    """
    Validate the minimum structure emitted by C++ when it requests more Z values.
    """
    required_top = ["action", "stage_name", "refine_id", "requested_points"]
    for key in required_top:
        if key not in data:
            raise ValueError(f"Missing required top-level field: {key}")

    if data["action"] != "queue_z_evaluations":
        raise ValueError(f"Unexpected action for queue payload: {data['action']}")

    requested_points = data["requested_points"]
    if not isinstance(requested_points, list) or len(requested_points) == 0:
        raise ValueError("requested_points must be a non-empty list")

    labels_seen = set()
    for point in requested_points:
        if "label" not in point or "t" not in point:
            raise ValueError(f"Invalid point entry: {point}")
        label = str(point["label"])
        if label in labels_seen:
            raise ValueError(f"Duplicate sample label: {label}")
        labels_seen.add(label)


# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------

def invoke_dispatch_lambda(point: Dict[str, Any]) -> Dict[str, Any]:
    """
    Call the existing Lambda that starts the distributed Z calculation for one t.

    Only orchestration here:
    - Python does not calculate Z
    - Python does not decide what t should be
    """
    payload = {
        "t": str(point["t"]),
        "image": AMI_IMAGE,
    }

    resp = lambda_client.invoke(
        FunctionName=DISPATCH_LAMBDA_NAME,
        InvocationType="RequestResponse",
        Payload=json.dumps(payload).encode("utf-8"),
    )

    response_payload = resp["Payload"].read().decode("utf-8") if "Payload" in resp else ""
    return {
        "status_code": resp.get("StatusCode"),
        "function_error": resp.get("FunctionError"),
        "response_payload": response_payload,
    }


def dispatch_requested_points(refine_id: str, stage_name: str, requested_points: List[Dict[str, Any]]) -> None:
    """
    Persist sample rows and dispatch all requested points for one stage iteration.
    """
    banner(f"Dispatch stage: {stage_name} | refine_id={refine_id}")
    put_sample_items(refine_id, stage_name, requested_points)
    log(f"Stored {len(requested_points)} requested point(s) in Dynamo")

    for idx, point in enumerate(requested_points, start=1):
        label = str(point["label"])
        t_value = str(point["t"])
        eval_id = make_eval_id(stage_name, label)

        log(f"[{idx}/{len(requested_points)}] Dispatching label={label} t={t_value}")
        dispatch_result = invoke_dispatch_lambda(point)

        update_eval_row(
            refine_id,
            eval_id,
            {
                "status": "DISPATCHED" if not dispatch_result.get("function_error") else "DISPATCH_ERROR",
                "dispatched_at": utc_now(),
                "lambda_status_code": dispatch_result.get("status_code"),
                "lambda_function_error": dispatch_result.get("function_error"),
                "lambda_response": dispatch_result.get("response_payload", ""),
            }
        )

        log(
            f"Dispatch result label={label} "
            f"status_code={dispatch_result.get('status_code')} "
            f"function_error={dispatch_result.get('function_error')}"
        )


# -----------------------------------------------------------------------------
# Waiting for Z completion
# -----------------------------------------------------------------------------

def ensure_z_job_ids(refine_id: str, eval_rows: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Resolve downstream Z job_id for each eval row and persist it if missing.
    """
    mapping: Dict[str, str] = {}

    for row in eval_rows:
        label = row.get("label", "")
        eval_id = row.get("eval_id", "")

        if not label or not eval_id:
            raise RuntimeError(f"Malformed eval row: {row}")

        if eval_failed(row):
            raise RuntimeError(f"Eval row already failed: eval_id={eval_id}")

        if not launch_successful(row):
            raise RuntimeError(f"Eval row not launch-successful yet: eval_id={eval_id} label={label}")

        z_job_id = row.get("z_job_id")
        if not z_job_id:
            z_job_id = parse_lambda_response_for_job_id(row.get("lambda_response"))
            update_eval_row(refine_id, eval_id, {"z_job_id": z_job_id})

        mapping[eval_id] = z_job_id

    return mapping


def wait_for_stage_results(refine_id: str, stage_name: str, expected_count: int) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """
    Wait until all downstream Z jobs for one stage are complete.

    Returns
    -------
    eval_rows:
        The stage-local eval rows from the refine samples table.
    z_jobs_by_eval:
        Mapping eval_id -> Z job item from the main Z jobs table.
    """
    banner(f"Waiting for stage results: {stage_name} | refine_id={refine_id}")
    start = time.time()

    while True:
        elapsed = int(time.time() - start)
        if elapsed > MAX_WAIT_SECONDS:
            raise TimeoutError(f"Timed out waiting for completed Z jobs for stage={stage_name}")

        eval_rows = get_eval_rows_for_stage(refine_id, stage_name)
        log(f"Fetched {len(eval_rows)} eval row(s) for stage={stage_name}")

        if len(eval_rows) < expected_count:
            warn(f"Expected {expected_count} eval rows, got {len(eval_rows)}. Waiting...")
            time.sleep(POLL_SECONDS)
            continue

        failed = [row for row in eval_rows if eval_failed(row)]
        if failed:
            raise RuntimeError(f"Detected failed eval rows: {[r.get('eval_id') for r in failed]}")

        z_job_ids = ensure_z_job_ids(refine_id, eval_rows)
        log(f"Resolved z_job_ids: {z_job_ids}")

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
            raise RuntimeError(f"Detected failed Z jobs: {failed_z}")

        if pending:
            log(f"Still waiting for stage={stage_name}: pending eval_ids={pending}")
            time.sleep(POLL_SECONDS)
            continue

        log(f"All downstream Z jobs completed for stage={stage_name}")
        return eval_rows, z_jobs_by_eval


# -----------------------------------------------------------------------------
# Build C++ input
# -----------------------------------------------------------------------------

def build_cpp_decision_input(
    refine_job: Dict[str, Any],
    current_request_payload: Dict[str, Any],
    eval_rows: List[Dict[str, Any]],
    z_jobs_by_eval: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Build the JSON consumed by the C++ decision engine.

    Important:
    - We preserve the order of current_request_payload['requested_points']
      so C++ receives the evaluations in the same logical order it asked for.
    - Python still does not decide what labels "mean".
      It only carries them forward.
    """
    rows_by_eval_id = {row["eval_id"]: row for row in eval_rows}
    requested_points = current_request_payload["requested_points"]

    evaluations = []
    for point in requested_points:
        label = str(point["label"])
        eval_id = make_eval_id(current_request_payload["stage_name"], label)

        if eval_id not in rows_by_eval_id:
            raise RuntimeError(f"Missing eval row for label={label} eval_id={eval_id}")

        row = rows_by_eval_id[eval_id]
        zjob = z_jobs_by_eval.get(eval_id)
        if not zjob:
            raise RuntimeError(f"Missing downstream Z job for eval_id={eval_id}")

        evaluations.append(
            {
                "label": label,
                "t": str(row.get("t", "")),
                "z_value": normalize_z_string(zjob.get("ztotal", "")),
                "status": "COMPLETED",
                "z_job_id": zjob.get("job_id", ""),
                "z_job_status": zjob.get("status", ""),
            }
        )

    return {
        "action": "decide_next_refinement_step",
        "stage_name": current_request_payload.get("stage_name", ""),
        "reason": "requested_samples_completed",
        "refine_id": refine_job["refine_id"],
        "input": {
            "t0": refine_job.get("input_t0", "") or current_request_payload.get("input", {}).get("t0", "")
        },
        "metadata": {
            **make_json_safe(current_request_payload.get("metadata", {})),
            "generated_at": utc_now(),
            "sample_count": len(evaluations),
        },
        "evaluations": evaluations,
    }


def iteration_input_json_path(refine_id: str, iteration_no: int) -> str:
    os.makedirs(REFINE_JSON_DIR, exist_ok=True)
    return os.path.join(REFINE_JSON_DIR, f"{refine_id}_iter{iteration_no:03d}_decision_input.json")


def iteration_output_json_path(refine_id: str, iteration_no: int) -> str:
    os.makedirs(REFINE_JSON_DIR, exist_ok=True)
    return os.path.join(REFINE_JSON_DIR, f"{refine_id}_iter{iteration_no:03d}_decision_output.json")


# -----------------------------------------------------------------------------
# Run C++
# -----------------------------------------------------------------------------

def run_cpp_decision(input_json_path: str, output_json_path: str) -> int:
    """
    Invoke the C++ decision engine.

    The executable contract is assumed to be the same one your current stage2
    Python already uses.
    """
    cmd = [STAGE2_CPP_BIN, "z_distributed_refine_stage2", input_json_path, output_json_path, "refinement_job"]
    banner("Launching C++ decision engine")
    log(f"Command: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout, end="", flush=True)
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr, flush=True)

    return result.returncode


# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

def refinement_loop(initial_json_path: str) -> int:
    """
    Main orchestration loop.

    Loop invariant:
    - current_payload is always the latest C++ output payload.
    - if current_payload.action == queue_z_evaluations, Python dispatches and waits.
    - then Python builds the next input for C++ and obtains the next output payload.
    - loop ends when C++ returns a terminal action.
    """
    banner("Distributed refinement loop started")

    current_payload = load_json_file(initial_json_path)
    validate_queue_payload(current_payload)

    refine_id = current_payload["refine_id"]
    put_or_update_refine_job(current_payload)

    log(f"Loaded initial payload for refine_id={refine_id}")
    update_job_status(
        refine_id,
        "RUNNING",
        {
            "current_stage_name": current_payload.get("stage_name", ""),
            "current_action": current_payload.get("action", ""),
            "initial_json_path": initial_json_path,
        },
    )

    iteration_no = 1
    while True:
        action = current_payload.get("action", "")
        base_stage_name = current_payload.get("stage_name", f"stage_iter_{iteration_no}")
        stage_name = f"{base_stage_name}_iter_{iteration_no:03d}"
        current_payload["stage_name"] = stage_name
        requested_points = current_payload.get("requested_points", [])

        banner(f"Iteration {iteration_no} | stage_name={stage_name} | action={action}")
        if iteration_no > MAX_ITERATIONS:
            update_job_status(
                refine_id,
                "FAILED_TO_CONVERGE",
                {
                    "monitor_error": f"Reached MAX_ITERATIONS={MAX_ITERATIONS}",
                    "last_stage_name": stage_name,
                    "last_action": action,
                },
            )
            err(f"FAILED_TO_CONVERGE: reached MAX_ITERATIONS={MAX_ITERATIONS}")
            return 9

        if action != "queue_z_evaluations":
            # Terminal or unexpected action returned directly from C++.
            log(f"C++ returned non-queue action={action}; stopping loop")
            update_job_status(
                refine_id,
                "COMPLETED" if action in TERMINAL_ACTIONS else "STOPPED",
                {
                    "final_action": action,
                    "final_stage_name": stage_name,
                    "final_reason": current_payload.get("reason", ""),
                },
            )
            return 0

        validate_queue_payload(current_payload)

        # Persist current queue request on job row for easy visual inspection.
        update_job_status(
            refine_id,
            "DISPATCHING",
            {
                "current_stage_name": stage_name,
                "current_action": action,
                "requested_points": requested_points,
                "requested_count": len(requested_points),
            },
        )

        # 1) Dispatch the points requested by C++.
        dispatch_requested_points(refine_id, stage_name, requested_points)

        # 2) Wait for all requested downstream Z jobs to finish.
        update_job_status(refine_id, "WAITING_EVALS", {"current_stage_name": stage_name})
        eval_rows, z_jobs_by_eval = wait_for_stage_results(refine_id, stage_name, len(requested_points))
        print_best_z_panel(current_payload, z_jobs_by_eval)

        # 3) Build next input JSON for C++.
        refine_job = get_refine_job(refine_id)
        if not refine_job:
            raise RuntimeError(f"Refine job disappeared from Dynamo: refine_id={refine_id}")

        decision_input = build_cpp_decision_input(refine_job, current_payload, eval_rows, z_jobs_by_eval)
        input_json_path = iteration_input_json_path(refine_id, iteration_no)
        output_json_path = iteration_output_json_path(refine_id, iteration_no)

        write_json_file(input_json_path, decision_input)
        log(f"Wrote C++ decision input JSON: {input_json_path}")
        log(f"Planned C++ decision output JSON: {output_json_path}")

        update_job_status(
            refine_id,
            "DECIDING",
            {
                "decision_input_json_path": input_json_path,
                "decision_output_json_path": output_json_path,
                "current_stage_name": stage_name,
            },
        )

        # 4) Run C++ decision engine.
        if os.path.exists(output_json_path):
            os.remove(output_json_path)

        rc = run_cpp_decision(input_json_path, output_json_path)

        if not os.path.exists(output_json_path):
            update_job_status(
                refine_id,
                "ERROR",
                {
                    "monitor_error": f"C++ returned rc={rc} but output file was not created",
                    "missing_output_json": output_json_path,
                    "failed_iteration": iteration_no,
                    "failed_stage_name": stage_name,
                },
            )
            err(f"C++ output JSON missing: {output_json_path}")
            return 8

        if rc != 0:
            update_job_status(
                refine_id,
                "ERROR",
                {
                    "monitor_error": f"C++ decision failed rc={rc}",
                    "failed_iteration": iteration_no,
                    "failed_stage_name": stage_name,
                },
            )
            err(f"C++ decision returned rc={rc}")
            return 7

        # 5) Load next C++ output and either loop or stop.
        current_payload = load_json_file(output_json_path)
        next_action = current_payload.get("action", "")
        next_stage = current_payload.get("stage_name", "")

        log(f"C++ produced next action={next_action} next_stage={next_stage}")
        update_job_status(
            refine_id,
            "DECISION_READY",
            {
                "latest_cpp_action": next_action,
                "latest_cpp_stage_name": next_stage,
                "latest_cpp_output_json": output_json_path,
                "latest_cpp_reason": current_payload.get("reason", ""),
            },
        )

        if next_action in TERMINAL_ACTIONS and next_action != "queue_z_evaluations":
            banner(f"Terminal action reached: {next_action}")
            update_job_status(
                refine_id,
                "COMPLETED",
                {
                    "final_action": next_action,
                    "final_stage_name": next_stage,
                    "final_reason": current_payload.get("reason", ""),
                    "final_output_json": output_json_path,
                },
            )
            return 0

        iteration_no += 1


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) != 2:
        err(f"Usage: {sys.argv[0]} <initial_cpp_output_json>")
        return 1

    initial_json_path = sys.argv[1]

    try:
        return refinement_loop(initial_json_path)

    except TimeoutError as e:
        err(str(e))
        try:
            data = load_json_file(initial_json_path)
            refine_id = data.get("refine_id")
            if refine_id:
                update_job_status(refine_id, "TIMEOUT", {"monitor_error": str(e)})
        except Exception:
            pass
        return 3

    except Exception as e:
        err(str(e))
        try:
            data = load_json_file(initial_json_path)
            refine_id = data.get("refine_id")
            if refine_id:
                update_job_status(refine_id, "ERROR", {"monitor_error": str(e)})
        except Exception:
            pass
        return 2


if __name__ == "__main__":
    raise SystemExit(main())