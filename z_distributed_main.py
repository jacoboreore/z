import subprocess
import math
import os
from decimal import Decimal, getcontext

getcontext().prec = 100

JOB_ID = "101"
CPP_EXE = "./z.exe"
PRE_FILE = "output/pre_calc.txt"
CHUNK_FILE = "chunk_results.txt"

def parse_precalc(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            data[k] = v
    data["N"] = int(data["N"])
    return data


def build_ranges(N, num_chunks):
    chunk_size = math.ceil(N / num_chunks)
    ranges = []

    n0 = 1
    chunk_id = 0
    while n0 <= N:
        n1 = min(n0 + chunk_size - 1, N)
        ranges.append((chunk_id, n0, n1))
        n0 = n1 + 1
        chunk_id += 1

    return ranges


def run_precalc(t_str):
    cmd = [CPP_EXE, "z_distributed_precalc", t_str, JOB_ID]
    print("Running:", " ".join(str(cmd)))
    subprocess.run(cmd, check=True)


def run_chunk(chunk_id, n0, n1):
    cmd = [CPP_EXE, "z_distributed_chunk", PRE_FILE, str(chunk_id), str(n0), str(n1), CHUNK_FILE, JOB_ID]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_chunk_results(path):
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


def main():
    t_str = input("Enter t: ").strip()
    num_chunks = int(input("Enter number of chunks: ").strip())

    if os.path.exists(CHUNK_FILE):
        os.remove(CHUNK_FILE)

    run_precalc(t_str)
    pre = parse_precalc(PRE_FILE)
    
    gabcke = Decimal(pre["gabcke"].replace("E", "e"))

    print(f"\nComputed N = {pre['N']}")
    ranges = build_ranges(pre["N"], num_chunks)

    print("\nChunk ranges:")
    for chunk_id, n0, n1 in ranges:
        print(f"  chunk {chunk_id}: [{n0}, {n1}]")

    for chunk_id, n0, n1 in ranges:
        run_chunk(chunk_id, n0, n1)

    results = parse_chunk_results(CHUNK_FILE)

    total = Decimal(0)
    for row in results:
        total += row["sum"]

    zmain = Decimal(2) * total
    ztotal = zmain + gabcke

    print("\nFinal assembly:")
    print("sum(chunks) =", total)
    print("Zmain = 2 * sum(chunks) =", zmain)
    print("Gabcke =", gabcke)
    print("Ztotal =", ztotal)

if __name__ == "__main__":
    main()