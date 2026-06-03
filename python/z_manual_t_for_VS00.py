#!/usr/bin/env python3
from decimal import Decimal, getcontext
import sys

getcontext().prec = 150

CFG_CONST = Decimal("0")8888833333335555555812.089999999999999996669330926
VALLEY_CONSTANT_STEP_DIV = Decimal("100")
M = Decimal("0.5")  # critical line real part

def compute_N_from_y(t: Decimal) -> Decimal:
    return t * t + M * M

def compute_y_from_N(N: Decimal) -> Decimal:
    return (N - M * M).sqrt()

def constant_step(N: Decimal) -> Decimal:
    if CFG_CONST == 0:
        return N.sqrt() / VALLEY_CONSTANT_STEP_DIV
    return CFG_CONST

def generate_vs_points(t0: str, count: int = 10):
    t = Decimal(t0)
    N = compute_N_from_y(t)

    step = constant_step(N)

    print(f"t0   = {t}")
    print(f"N0   = {N}")
    print(f"step = {step}")
    print()

    for i in range(1, count + 1):
        N += step
        t = compute_y_from_N(N)
        #print(f"{i:02d}  t = {t}")
        #print(f"    N = {N}")
        print(f"{t}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python vs_points.py <t0> [count]")
        sys.exit(1)

    t0 = sys.argv[1]
    count = int(sys.argv[2]) if len(sys.argv) >= 3 else 10
    generate_vs_points(t0, count)