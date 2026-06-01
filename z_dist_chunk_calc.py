import json
import os
import time
import shutil
import logging
import tempfile
import subprocess
from decimal import Decimal, getcontext
import boto3
from botocore.exceptions import ClientError

getcontext().prec = 100

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
AWS_REGION = os.environ.get("AWS_REGION", "us-east-1")
QUEUE_URL = os.environ.get("QUEUE_URL", "https://sqs.us-east-1.amazonaws.com/407875882825/valley-scanner-dev-work")

CPP_EXE = os.environ.get("CPP_EXE", "./bin/z")
POLL_WAIT_SECONDS = int(os.environ.get("POLL_WAIT_SECONDS", "20"))
VISIBILITY_TIMEOUT = int(os.environ.get("VISIBILITY_TIMEOUT", "900"))
SLEEP_EMPTY_SECONDS = int(os.environ.get("SLEEP_EMPTY_SECONDS", "5"))

KEEP_WORKDIR = os.environ.get("KEEP_WORKDIR", "false").lower() == "true"
BASE_WORKDIR = os.environ.get("BASE_WORKDIR", "./worker_runs")
LOGS_DIR = os.environ.get("LOGS_DIR", "./logs")

ENABLE_DYNAMO_STATUS = os.environ.get("ENABLE_DYNAMO_STATUS", "true").lower() == "true"
CHUNKS_TABLE = os.environ.get("CHUNKS_TABLE", "valley-scanner-dev-chunks")

ENABLE_S3_UPLOAD = os.environ.get("ENABLE_S3_UPLOAD", "true").lower() == "true"
RESULTS_BUCKET = os.environ.get("RESULTS_BUCKET", "jacob-z")
RESULTS_PREFIX = os.environ.get("RESULTS_PREFIX", "z-dist-results")

# -----------------------------------------------------------------------------
# ROOT LOGGING
# -----------------------------------------------------------------------------
import logging
import sys

# Root logger
root_log = logging.getLogger()
root_log.setLevel(logging.INFO)

# Remove existing handlers completely
for h in root_log.handlers[:]:
    root_log.removeHandler(h)

# Create explicit console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)

formatter = logging.Formatter(
    '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)

console_handler.setFormatter(formatter)

root_log.addHandler(console_handler)

# Optional named logger
logger = logging.getLogger("valley-scanner")
logger.setLevel(logging.INFO)

logger.info("Logging initialized correctly")

# -----------------------------------------------------------------------------
# AWS CLIENTS
# -----------------------------------------------------------------------------
sqs = boto3.client("sqs", region_name=AWS_REGION)
dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION) if ENABLE_DYNAMO_STATUS else None
s3 = boto3.client("s3", region_name=AWS_REGION) if ENABLE_S3_UPLOAD else None

chunks_table = dynamodb.Table(CHUNKS_TABLE) if ENABLE_DYNAMO_STATUS and CHUNKS_TABLE else None


# -----------------------------------------------------------------------------
# LOGGER HELPERS
# -----------------------------------------------------------------------------
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def build_chunk_logger(chunk_index: int) -> tuple[logging.Logger, str]:
    """
    Creates a dedicated logger for a chunk.
    File path: ./logs/log_chunk_xxxxxx.log
    """
    ensure_dir(LOGS_DIR)
    logfile = os.path.join(LOGS_DIR, f"log_chunk_{chunk_index:06d}.log")

    logger_name = f"valley-worker-chunk-{chunk_index:06d}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    # Avoid duplicate handlers if logger is reused
    if logger.handlers:
        return logger, logfile

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")

    file_handler = logging.FileHandler(logfile, encoding="utf-8")
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger, logfile


def close_logger(logger: logging.Logger) -> None:
    handlers = list(logger.handlers)
    for handler in handlers:
        try:
            handler.flush()
            handler.close()
        except Exception:
            pass
        logger.removeHandler(handler)


# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
def parse_chunk_results(path: str) -> list[dict]:
    results = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            row = {}
            parts = line.split(",")
            for p in parts:
                k, v = p.split("=", 1)
                row[k] = v

            results.append({
                "chunk_id": int(row["chunk_id"]),
                "n0": int(row["n0"]),
                "n1": int(row["n1"]),
                "sum": Decimal(row["sum"].replace("E", "e"))
            })

    results.sort(key=lambda x: x["chunk_id"])
    return results


def update_chunk_status(
    job_id: str,
    chunk_id: str,
    status: str,
    extra: dict | None = None,
    logger: logging.Logger | None = None
) -> None:
    if not chunks_table:
        return

    item = {
        "job_id": str(job_id),
        "chunk_id": str(chunk_id),
        "status": status,
        "updated_at": int(time.time())
    }
    if extra:
        item.update(extra)

    try:
        chunks_table.put_item(Item=item)
        if logger:
            logger.info(
                "Dynamo chunk status updated: job_id=%s chunk_id=%s status=%s",
                job_id, chunk_id, status
            )
    except ClientError as e:
        if logger:
            logger.exception("Failed to update chunk status: %s", e)


def upload_file_to_s3(
    local_path: str,
    bucket: str,
    key: str,
    logger: logging.Logger | None = None
) -> None:
    if not s3:
        return
    s3.upload_file(local_path, bucket, key)
    if logger:
        logger.info("Uploaded to S3: s3://%s/%s", bucket, key)


def run_cmd(cmd: list[str], cwd: str, env: dict | None = None, logger: logging.Logger | None = None) -> None:
    if logger:
        logger.info("Running command: %s", " ".join(cmd))

    result = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True
    )

    if logger:
        logger.info("Command return code: %s", result.returncode)
        if result.stdout:
            logger.info("Command stdout:\n%s", result.stdout.strip())
        if result.stderr:
            logger.info("Command stderr:\n%s", result.stderr.strip())

    result.check_returncode()


def validate_message(body: dict) -> dict:
    required = [
        "job_id",
        "chunk_id",
        "chunk_index",
        "n0",
        "n1",
        "size",
        "t",
        "N",
        "theta",
        "gabcke",
        "useC1C2",
        "prec",
        "number_of_cpus",
    ]
    missing = [k for k in required if k not in body]
    if missing:
        raise ValueError(f"Missing required fields in SQS message: {missing}")

    msg = {
        "job_id": str(body["job_id"]),
        "chunk_id": str(body["chunk_id"]),          # tracking id, e.g. chunk-000001
        "chunk_index": int(body["chunk_index"]),    # numeric id for z.exe
        "n0": int(body["n0"]),
        "n1": int(body["n1"]),
        "size": int(body["size"]),
        "t": str(body["t"]).strip(),
        "N": int(body["N"]),
        "theta": str(body["theta"]).strip(),
        "gabcke": str(body["gabcke"]).strip(),
        "useC1C2": int(body["useC1C2"]),
        "prec": int(body["prec"]),
        "number_of_cpus": int(body["number_of_cpus"]),
    }

    if msg["n0"] <= 0 or msg["n1"] < msg["n0"]:
        raise ValueError(f"Invalid range: n0={msg['n0']}, n1={msg['n1']}")

    expected_size = msg["n1"] - msg["n0"] + 1
    if msg["size"] != expected_size:
        raise ValueError(
            f"Invalid size in message: size={msg['size']} expected={expected_size}"
        )

    return msg


def build_workdir(job_id: str, chunk_id: str) -> str:
    ensure_dir(BASE_WORKDIR)
    safe_chunk_id = chunk_id.replace("/", "_")
    return tempfile.mkdtemp(prefix=f"job_{job_id}_{safe_chunk_id}_", dir=BASE_WORKDIR)


def write_precalc_from_message(pre_file: str, payload: dict, logger: logging.Logger | None = None) -> None:
    """
    Reconstruct the pre_calc.txt expected by z_distributed_chunk
    using the values already computed upstream.
    """
    content = [
        f"N={payload['N']}",
        f"theta={payload['theta']}",
        f"gabcke={payload['gabcke']}",
        f"useC1C2={payload['useC1C2']}",
        f"prec={payload['prec']}",
        f"t={payload['t']}",
    ]

    with open(pre_file, "w", encoding="utf-8") as f:
        f.write("\n".join(content))
        f.write("\n")

    if logger:
        logger.info("Created local pre_calc file: %s", pre_file)
        logger.info("pre_calc.txt content:\n%s", "\n".join(content))


def build_child_env(payload: dict, logger: logging.Logger | None = None) -> dict:
    env = os.environ.copy()
    cpus = str(payload["number_of_cpus"])
    env["OMP_NUM_THREADS"] = cpus
    env["NUMBER_OF_CPUS"] = cpus

    if logger:
        logger.info(
            "Child process env prepared: OMP_NUM_THREADS=%s NUMBER_OF_CPUS=%s",
            env["OMP_NUM_THREADS"],
            env["NUMBER_OF_CPUS"]
        )

    return env


# -----------------------------------------------------------------------------
# CORE PROCESSING
# -----------------------------------------------------------------------------
def process_message(msg: dict) -> None:
    receipt_handle = msg["ReceiptHandle"]
    raw_body = msg["Body"]

    # Temporary logger until chunk_index is known
    temp_logger = root_log
    chunk_logger = None

    try:
        temp_logger.info("Received raw SQS message body: %s", raw_body)
        body = json.loads(raw_body)
        payload = validate_message(body)

        job_id = payload["job_id"]
        chunk_id = payload["chunk_id"]
        chunk_index = payload["chunk_index"]
        n0 = payload["n0"]
        n1 = payload["n1"]

        chunk_logger, logfile = build_chunk_logger(chunk_index)
        chunk_logger.info("============================================================")
        chunk_logger.info("Starting processing for chunk_index=%s", chunk_index)
        chunk_logger.info("Chunk log file: %s", logfile)
        chunk_logger.info("Raw SQS message body: %s", raw_body)
        chunk_logger.info("Validated payload: %s", json.dumps(payload))

        workdir = build_workdir(job_id, chunk_id)
        output_dir = os.path.join(workdir, "output")
        ensure_dir(output_dir)

        pre_file = os.path.join(output_dir, "pre_calc.txt")
        chunk_file = os.path.join(workdir, "chunk_results.txt")
        summary_path = os.path.join(workdir, f"summary_{chunk_id}.json")

        exe_name = os.path.basename(CPP_EXE)
        local_exe = os.path.join(workdir, exe_name)

        chunk_logger.info("Working directory: %s", workdir)
        chunk_logger.info("Output directory: %s", output_dir)
        chunk_logger.info("Local executable path: %s", local_exe)
        chunk_logger.info("Chunk output file path: %s", chunk_file)
        chunk_logger.info("Summary file path: %s", summary_path)

        update_chunk_status(
            job_id,
            chunk_id,
            "RUNNING",
            {
                "chunk_index": chunk_index,
                "t": payload["t"],
                "n0": n0,
                "n1": n1,
                "size": payload["size"],
                "N": payload["N"],
                "theta": payload["theta"],
                "gabcke": payload["gabcke"],
                "useC1C2": payload["useC1C2"],
                "prec": payload["prec"],
                "number_of_cpus": payload["number_of_cpus"],
            },
            logger=chunk_logger
        )

        chunk_logger.info("Copying executable from %s to %s", CPP_EXE, local_exe)
        shutil.copy2(CPP_EXE, local_exe)
        os.chmod(local_exe, 0o755)
        shutil.copytree("./conf", os.path.join(workdir, "conf"), dirs_exist_ok=True)
        chunk_logger.info("Executable copied and chmod applied")

        write_precalc_from_message(pre_file, payload, logger=chunk_logger)
        child_env = build_child_env(payload, logger=chunk_logger)

        cmd = [
            local_exe,
            "z_distributed_chunk",
            pre_file,
            str(chunk_index),
            str(n0),
            str(n1),
            chunk_file,
            job_id,
        ]

        chunk_logger.info(
            "Invoking chunk executable with args: mode=%s pre_file=%s chunk_index=%s n0=%s n1=%s chunk_file=%s job_id=%s",
            "z_distributed_chunk",
            pre_file,
            chunk_index,
            n0,
            n1,
            chunk_file,
            job_id
        )

        run_cmd(cmd, cwd=workdir, env=child_env, logger=chunk_logger)

        if not os.path.exists(chunk_file):
            raise FileNotFoundError(f"chunk_results.txt not found after chunk execution: {chunk_file}")

        chunk_logger.info("Chunk results file found: %s", chunk_file)

        with open(chunk_file, "r", encoding="utf-8") as f:
            raw_chunk_results = f.read()
        chunk_logger.info("chunk_results.txt content:\n%s", raw_chunk_results.strip())

        results = parse_chunk_results(chunk_file)
        chunk_logger.info("Parsed %s result rows from chunk_results.txt", len(results))

        this_row = next((r for r in results if r["chunk_id"] == chunk_index), None)
        if this_row is None:
            raise ValueError(f"No result found for chunk_index={chunk_index} in chunk_results.txt")

        chunk_sum = this_row["sum"]
        chunk_logger.info(
            "Matched result row: chunk_id=%s n0=%s n1=%s sum=%s",
            this_row["chunk_id"],
            this_row["n0"],
            this_row["n1"],
            str(chunk_sum)
        )

        summary = {
            "job_id": job_id,
            "chunk_id": chunk_id,
            "chunk_index": chunk_index,
            "t": payload["t"],
            "N": payload["N"],
            "theta": payload["theta"],
            "gabcke": payload["gabcke"],
            "useC1C2": payload["useC1C2"],
            "prec": payload["prec"],
            "number_of_cpus": payload["number_of_cpus"],
            "n0": n0,
            "n1": n1,
            "size": payload["size"],
            "chunk_sum": str(chunk_sum),
            "status": "DONE",
        }

        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

        chunk_logger.info("Summary JSON written: %s", summary_path)
        chunk_logger.info("Summary content: %s", json.dumps(summary))

        if ENABLE_S3_UPLOAD:
            prefix = f"{RESULTS_PREFIX}/{job_id}/{chunk_id}"
            chunk_logger.info("S3 upload enabled, prefix=%s", prefix)
            upload_file_to_s3(pre_file, RESULTS_BUCKET, f"{prefix}/pre_calc.txt", logger=chunk_logger)
            upload_file_to_s3(chunk_file, RESULTS_BUCKET, f"{prefix}/chunk_results.txt", logger=chunk_logger)
            upload_file_to_s3(summary_path, RESULTS_BUCKET, f"{prefix}/summary.json", logger=chunk_logger)

        update_chunk_status(
            job_id,
            chunk_id,
            "DONE",
            {
                "chunk_index": chunk_index,
                "t": payload["t"],
                "n0": n0,
                "n1": n1,
                "size": payload["size"],
                "N": payload["N"],
                "theta": payload["theta"],
                "gabcke": payload["gabcke"],
                "chunk_sum": str(chunk_sum),
                "useC1C2": payload["useC1C2"],
                "prec": payload["prec"],
                "number_of_cpus": payload["number_of_cpus"],
            },
            logger=chunk_logger
        )

        sqs.delete_message(QueueUrl=QUEUE_URL, ReceiptHandle=receipt_handle)
        chunk_logger.info("Deleted SQS message successfully: job_id=%s chunk_id=%s", job_id, chunk_id)
        chunk_logger.info("Chunk processing finished successfully")
        chunk_logger.info("============================================================")

    except Exception as e:
        logger = chunk_logger if chunk_logger else temp_logger
        logger.exception("Chunk processing failed: %s", e)

        try:
            body = json.loads(raw_body)
            payload = validate_message(body)
            update_chunk_status(
                payload["job_id"],
                payload["chunk_id"],
                "FAILED",
                {
                    "chunk_index": payload["chunk_index"],
                    "error": str(e),
                },
                logger=logger if isinstance(logger, logging.Logger) else None
            )
        except Exception as inner:
            logger.exception("Failed to update FAILED status after processing error: %s", inner)

        raise

    finally:
        try:
            if 'workdir' in locals():
                if KEEP_WORKDIR:
                    if chunk_logger:
                        chunk_logger.info("Keeping workdir for inspection: %s", workdir)
                else:
                    shutil.rmtree(workdir, ignore_errors=True)
                    if chunk_logger:
                        chunk_logger.info("Workdir removed: %s", workdir)
        finally:
            if chunk_logger:
                close_logger(chunk_logger)


# -----------------------------------------------------------------------------
# POLL LOOP
# -----------------------------------------------------------------------------
def poll_once() -> bool:
    resp = sqs.receive_message(
        QueueUrl=QUEUE_URL,
        MaxNumberOfMessages=1,
        WaitTimeSeconds=POLL_WAIT_SECONDS,
        VisibilityTimeout=VISIBILITY_TIMEOUT,
        AttributeNames=["All"],
        MessageAttributeNames=["All"],
    )

    messages = resp.get("Messages", [])
    if not messages:
        return False

    process_message(messages[0])
    return True


def main() -> None:
    ensure_dir(LOGS_DIR)
    ensure_dir(BASE_WORKDIR)

    root_log.info("Valley Scanner worker started")
    root_log.info("QUEUE_URL=%s", QUEUE_URL)
    root_log.info("CPP_EXE=%s", CPP_EXE)
    root_log.info("AWS_REGION=%s", AWS_REGION)
    root_log.info("LOGS_DIR=%s", LOGS_DIR)
    root_log.info("BASE_WORKDIR=%s", BASE_WORKDIR)
    root_log.info("ENABLE_DYNAMO_STATUS=%s", ENABLE_DYNAMO_STATUS)
    root_log.info("ENABLE_S3_UPLOAD=%s", ENABLE_S3_UPLOAD)

    while True:
        try:
            got_one = poll_once()
            if not got_one:
                root_log.info("No messages available, sleeping %s sec", SLEEP_EMPTY_SECONDS)
                time.sleep(SLEEP_EMPTY_SECONDS)
        except KeyboardInterrupt:
            root_log.info("Worker interrupted, exiting")
            break
        except Exception as e:
            root_log.exception("Worker loop error: %s", e)
            time.sleep(3)


if __name__ == "__main__":
    main()