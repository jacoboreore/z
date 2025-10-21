#!/usr/bin/env python3
import argparse, csv, json, os, sys, tempfile
from typing import List, Tuple

def read_csv_safely(path: str) -> List[Tuple[float, float, float]]:
    """
    Reads a CSV that may be concurrently appended to.
    Drops an incomplete trailing line if present.
    Expects header: t,absZ,N
    Returns a list of (t, absZ, N) as floats.
    """
    if not os.path.exists(path):
        return []

    # Read raw to detect incomplete last line robustly
    with open(path, "rb") as f:
        data = f.read()

    # Decode forgivingly (in case writer wrote partial UTF-8 chunk)
    text = data.decode("utf-8", errors="replace")

    # Keep only complete lines (those that end with newline)
    # splitlines(keepends=True) preserves end markers
    lines = text.splitlines(keepends=True)
    complete_lines = [ln for ln in lines if ln.endswith("\n") or ln.endswith("\r\n")]

    if not complete_lines:
        return []

    # Join back and parse as CSV
    joined = "".join(complete_lines)
    rdr = csv.reader(joined.splitlines())
    rows = list(rdr)
    if not rows:
        return []

    # Header check (best-effort)
    header = [h.strip() for h in rows[0]]
    start_idx = 1
    if len(header) >= 3 and header[0].lower() == "t" and header[1].lower() in ("absz", "|z|") and header[2].lower() == "n":
        pass  # header is present
    else:
        # no header — treat all as data
        start_idx = 0

    out = []
    for r in rows[start_idx:]:
        if len(r) < 3:
            continue
        try:
            t = float(r[0])
            absZ = float(r[1])
            N = float(r[2])
            out.append((t, absZ, N))
        except Exception:
            # skip malformed/partial records
            continue
    return out

def write_json_atomically(obj, out_path: str):
    """
    Writes JSON atomically: write to temp in same dir, then replace.
    Ensures readers never see a partial file.
    """
    d = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(d, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=d, delete=False) as tmp:
        tmp_path = tmp.name
        json.dump(obj, tmp, ensure_ascii=False, separators=(",", ":"), indent=2)
        tmp.flush()
        os.fsync(tmp.fileno())
    os.replace(tmp_path, out_path)

def main():
    ap = argparse.ArgumentParser(description="Convert zeros.csv to compact JSON for FE polling.")
    ap.add_argument("--zeros-csv", required=True, help="Path to zeros.csv")
    ap.add_argument("--out", required=True, help="Path to write JSON (e.g., /home/USER/jacob/progress.json)")
    ap.add_argument("--total", type=int, required=True, help="Total zeros expected for this job (denominator)")
    ap.add_argument("--job-id", required=True, help="Job ID string")
    args = ap.parse_args()

    zeros = read_csv_safely(args.zeros_csv)
    processed = len(zeros)

    payload = {
        "jobId": args.job_id,
        "processed": f"{processed}/{args.total}",
        "zeros": [[t, absZ, N] for (t, absZ, N) in zeros]
    }

    write_json_atomically(payload, args.out)

if __name__ == "__main__":
    # never crash the timer; exit 0 on handled errors
    try:
        main()
    except Exception as e:
        # print to stderr for journald
        print(f"zeros_to_json error: {e}", file=sys.stderr)
        sys.exit(0)
