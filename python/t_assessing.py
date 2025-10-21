from flint import acb, arb, ctx
import math

def set_precision_bits(bits: int):
    dps = int(math.ceil(bits * math.log10(2))) + 20  # +20 guard
    ctx.dps = dps

def abs_ball(x):
    return x.abs() if hasattr(x, "abs") else abs(x)

def ball_is_small(a: arb, thresh: float) -> bool:
    mid = a.mid()
    rad = a.rad()
    return (rad >= abs(mid)) or (abs(mid) + rad < thresh)

def safe_upper(a: arb) -> float:
    return float(a.mid()) + float(a.rad())

def safe_lower_nonneg(a: arb) -> float:
    # lower bound for nonnegative balls (clamped at 0)
    return max(0.0, float(a.mid()) - float(a.rad()))

def zeta_and_deriv(s: acb):
    if hasattr(s, "zeta"):
        z0 = s.zeta()
        try:
            z1 = s.zeta(1)     # first derivative wrt s
        except TypeError:
            from flint import zeta as zeta_fn
            z1 = zeta_fn(s, 1)
        return z0, z1
    from flint import zeta as zeta_fn
    return zeta_fn(s), zeta_fn(s, 1)

def zeta_second(s: acb):
    # second derivative if available
    if hasattr(s, "zeta"):
        try:
            return s.zeta(2)
        except TypeError:
            from flint import zeta as zeta_fn
            return zeta_fn(s, 2)
    from flint import zeta as zeta_fn
    return zeta_fn(s, 2)

def spacing_est(t_f: float) -> float:
    return (2.0*math.pi) / math.log(t_f/(2.0*math.pi))

def validate_zeta_zero(sigma, t, prec_bits=256, zero_threshold=1e-10, newton_iters=6, try_halley=True):
    """
    Check zeta(sigma + i t) using python-flint and optionally refine t (Newton on the vertical line).
    Returns: (zeta_val, zero_likely, refined_t)
    """
    set_precision_bits(prec_bits)

    # make sure t is an arb to keep ball arithmetic throughout
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

            # guard: derivative must not include 0
            a1 = abs_ball(z1)
            if safe_lower_nonneg(a1) == 0.0:
                # derivative too flat; increase precision a bit and retry once
                ctx.dps += 10
                z0, z1 = zeta_and_deriv(acb(arb(sigma), refined_t))
                a1 = abs_ball(z1)
                if safe_lower_nonneg(a1) == 0.0:
                    # bail to very small nudge instead of dividing by near-zero
                    refined_t -= arb("1e-12")
                    continue

            # Newton step on vertical line: dt = Im(z0/z1)
            q = z0 / z1
            dt = q.imag  # arb

            # optional Halley step if derivative is healthy
            if try_halley:
                z2 = zeta_second(s)
                # F(t)=zeta; F' = i z1; F'' = - z2
                F  = z0
                Fp = acb(arb(0), arb(1)) * z1     # i*z1
                Fpp= -z2                          # -z2
                denom = 2*(Fp*Fp) - F*Fpp
                if not ball_is_small(abs_ball(denom), 1e-40):
                    halley_dt = (2*F*Fp) / denom
                    # stick to real step; take real part
                    dt = halley_dt.real

            # cap step by spacing
            t_float = float(refined_t.mid())
            cap = 0.3 * max(1e-12, spacing_est(max(t_float, 10.0)))
            dt_f = float(dt)
            if abs(dt_f) > cap:
                dt *= (cap / abs(dt_f))

            # backtrack if |zeta| doesn't drop (by UB)
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
                # tiny nudge and continue
                refined_t -= arb("1e-12")

        # final report
        s_final = acb(arb(sigma), refined_t)
        z_final, z1_final = zeta_and_deriv(s_final)
        a_final = abs_ball(z_final)
        print(f"refined t   = {refined_t}")
        print(f"zeta(final) = {z_final}")
        print(f"|zeta|_UB   = {safe_upper(a_final)}")
        zero_likely = ball_is_small(a_final, zero_threshold)

        # (optional) simple-zero check: derivative LB > 0
        a1_final = abs_ball(z1_final)
        print(f"|zeta'|_LB  = {safe_lower_nonneg(a1_final)}")

    return zeta_val, zero_likely, refined_t


sigma = 0.5
#t_large = "4400961943973.052602062600039768164996"
#t_large = "4400961943973.05260206260003976816499647134924840351522867308276"
t_large = "88888888898.2973446644705945344"

# start with comfortable precision and a stricter threshold
validate_zeta_zero(sigma, t_large, prec_bits=150, zero_threshold=1e-8, newton_iters=1)
