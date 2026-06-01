#!/usr/bin/env python3
"""
Independent zero-check script using python-flint (ARB/ACB ball arithmetic).

- Validates ζ(σ + i t) with rigorous upper bounds on |ζ|.
- Optionally performs a few Newton/Halley steps constrained to the vertical line.
- Prints |ζ|_UB and a conservative lower bound for |ζ'| to hint at simple zeros.

Usage (module):
    from validate_zero import validate_zeta_zero
    validate_zeta_zero(0.5, "1e12", prec_bits=200)

CLI (optional):
    python validate_zero.py 0.5 1e12 --prec 200 --iters 3

Requirements:
    python -m pip install python-flint
"""

from flint import acb, arb, ctx
import math
import sys
from typing import Tuple, Union, Optional

NumberLike = Union[str, float, int, arb]

def set_precision_bits(bits: int) -> int:
    """
    Set python-flint decimal precision based on binary bits.
    Returns the dps selected (with +20 guard).
    """
    dps = int(math.ceil(bits * math.log10(2))) + 20  # +20 guard digits
    ctx.dps = dps
    return dps

def abs_ball(x):
    return x.abs() if hasattr(x, "abs") else abs(x)

def ball_is_small(a: arb, thresh: float) -> bool:
    """
    Return True if the ball 'a' is consistent with |a| < thresh,
    or if its radius already straddles zero (mid≈0 within radius).
    """
    mid = a.mid()
    rad = a.rad()
    return (rad >= abs(mid)) or (abs(mid) + rad < thresh)

def safe_upper(a: arb) -> float:
    """Upper bound of a nonnegative ball magnitude."""
    return float(a.mid()) + float(a.rad())

def safe_lower_nonneg(a: arb) -> float:
    """Lower bound of a magnitude ball, clamped to 0."""
    return max(0.0, float(a.mid()) - float(a.rad()))

def zeta_and_deriv(s: acb):
    """
    Return (ζ(s), ζ'(s)). Uses s.zeta(k) when available, else module function.
    """
    if hasattr(s, "zeta"):
        z0 = s.zeta()
        try:
            z1 = s.zeta(1)
        except TypeError:
            from flint import zeta as zeta_fn
            z1 = zeta_fn(s, 1)
        return z0, z1
    from flint import zeta as zeta_fn
    return zeta_fn(s), zeta_fn(s, 1)

def zeta_second(s: acb):
    """Return ζ''(s) if available."""
    if hasattr(s, "zeta"):
        try:
            return s.zeta(2)
        except TypeError:
            from flint import zeta as zeta_fn
            return zeta_fn(s, 2)
    from flint import zeta as zeta_fn
    return zeta_fn(s, 2)

def spacing_est(t_f: float) -> float:
    """
    Local spacing estimate s ≈ 2π / log(t / 2π).
    Guarded so log argument ≥ exp(1) to avoid negative/zero logs.
    """
    t_guard = max(t_f, 10.0)
    x = t_guard / (2.0 * math.pi)
    # ensure denominator is positive and not tiny
    denom = max(math.log(max(x, math.e)), 1.0)
    return (2.0 * math.pi) / denom

def validate_zeta_zero(
    sigma: float,
    t: NumberLike,
    prec_bits: int = 256,
    zero_threshold: float = 1e-10,
    newton_iters: int = 6,
    try_halley: bool = True
) -> Tuple[acb, bool, arb]:
    """
    Validate ζ(σ + i t) using python-flint and optionally refine t along the
    vertical line (Newton/Halley with step caps and backtracking).

    Returns:
        (zeta_val, zero_likely, refined_t)

    zero_likely is based on a ball test:
        rad ≥ |mid|  OR  |mid| + rad < zero_threshold
    """
    set_precision_bits(prec_bits)

    # keep ball arithmetic throughout
    t_arb = arb(str(t)) if not isinstance(t, arb) else t
    s = acb(arb(sigma), t_arb)
    zeta_val, zeta_der = zeta_and_deriv(s)
    a = abs_ball(zeta_val)

    zero_likely = ball_is_small(a, zero_threshold)
    print(f"s = {sigma} + {t} i")
    print(f"zeta(s)     = {zeta_val}")
    print(f"|zeta(s)|_UB ≤ {safe_upper(a)}")

    refined_t = arb(t_arb)

    if not zero_likely:
        for it in range(newton_iters):
            s = acb(arb(sigma), refined_t)
            z0, z1 = zeta_and_deriv(s)
            a0 = abs_ball(z0)

            if ball_is_small(a0, zero_threshold):
                print(f"[newton] converged in {it} iters")
                zero_likely = True
                break

            # Guard: derivative must exclude 0
            a1 = abs_ball(z1)
            if safe_lower_nonneg(a1) == 0.0:
                # bump precision and retry once
                ctx.dps += 10
                z0, z1 = zeta_and_deriv(acb(arb(sigma), refined_t))
                a1 = abs_ball(z1)
                if safe_lower_nonneg(a1) == 0.0:
                    # tiny nudge if still too flat
                    refined_t -= arb("1e-12")
                    continue

            # Newton step on vertical line: dt = Im(z0/z1)
            q = z0 / z1
            dt = q.imag  # arb

            # Optional Halley correction
            if try_halley:
                z2 = zeta_second(s)
                # F(t)=ζ(σ+i t); F' = i ζ'(s); F'' = -ζ''(s)
                F   = z0
                Fp  = acb(arb(0), arb(1)) * z1
                Fpp = -z2
                denom = 2*(Fp*Fp) - F*Fpp
                if not ball_is_small(abs_ball(denom), 1e-40):
                    halley_dt = (2*F*Fp) / denom
                    dt = halley_dt.real  # stay on vertical line

            # Cap step by spacing
            t_float = float(refined_t.mid())
            cap = 0.3 * max(1e-12, spacing_est(t_float))
            dt_f = float(dt)
            if abs(dt_f) > cap:
                dt *= (cap / abs(dt_f))

            # Backtrack if |ζ| doesn't drop (using UB)
            old_ub = safe_upper(a0)
            ok = False
            trial_dt = arb(dt)
            for _ in range(5):
                trial_t = refined_t - trial_dt
                z_try, _ = zeta_and_deriv(acb(arb(sigma), trial_t))
                ub = safe_upper(abs_ball(z_try))
                if ub < old_ub:
                    refined_t = trial_t
                    ok = True
                    break
                trial_dt *= arb("0.5")
            if not ok:
                refined_t -= arb("1e-12")

        # Final report
        s_final = acb(arb(sigma), refined_t)
        z_final, z1_final = zeta_and_deriv(s_final)
        a_final = abs_ball(z_final)
        print(f"refined t   = {refined_t}")
        print(f"zeta(final) = {z_final}")
        print(f"|zeta|_UB   = {safe_upper(a_final)}")
        zero_likely = ball_is_small(a_final, zero_threshold)

        # (optional) simple-zero hint: derivative LB > 0
        a1_final = abs_ball(z1_final)
        print(f"|zeta'|_LB  = {safe_lower_nonneg(a1_final)}")

    return zeta_val, zero_likely, refined_t


# ---------------------- Optional CLI ----------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python validate_zero.py <sigma> <t> [--prec BITS] [--iters N] [--thr THRESH]")
        sys.exit(2)
    sigma = float(sys.argv[1])
    t_in = sys.argv[2]
    prec = 256
    iters = 3
    thr = 1e-10
    for i, arg in enumerate(sys.argv[3:], start=3):
        if arg == "--prec" and i+1 < len(sys.argv):
            prec = int(sys.argv[i+1])
        if arg == "--iters" and i+1 < len(sys.argv):
            iters = int(sys.argv[i+1])
        if arg == "--thr" and i+1 < len(sys.argv):
            thr = float(sys.argv[i+1])
    validate_zeta_zero(sigma, t_in, prec_bits=prec, zero_threshold=thr, newton_iters=iters)
