#include <cmath>
#include <algorithm>
#include <sstream>
#include "refine_newton.h"
#include "confirm_zero.h" // est_spacing, use_threaded_for_t
// #include "z_calculator_small.h"   // compute_Z_small(...)
#include "z_calculator.h"
// #include "z_calculator_rs_gabcke.h"
#include "settings.h"
#include "logs.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ==================== helpers ==================== */

static inline int mpfr_sign(const mpfr_t x) { return mpfr_sgn(x); }
static inline bool opposite_sign(const mpfr_t a, const mpfr_t b)
{
    int sa = mpfr_sign(a), sb = mpfr_sign(b);
    return sa && sb && (sa != sb);
}

static inline double stage_step_tol(const mpfr_t t, double tol_abs)
{
    const double s = est_spacing(t);
    const double k = std::clamp(std::sqrt(std::max(tol_abs, 1e-20)), 1e-4, 0.2);
    const double raw = k * s;
    return std::clamp(raw, CFG().h_min, 0.2 * s);
}

// choose h for derivative (scale with spacing + precision)
static inline double choose_h(const mpfr_t t)
{
    const double s = est_spacing(t);
    const double eps = std::ldexp(1.0, -(int)CFG().mpfr_prec); // 2^(-P)
    const double mag = std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
    const double h0 = std::pow(eps, 1.0 / 3.0) * mag;
    return std::clamp(h0, 1e-9, 0.1 * s);
}

/* ============== lightweight Pegasus burst ============== */

// refine a pre-existing bracket [a,b] with fa*fb<0]; a few steps only
static bool pegasus_refine_bracket(
    mpfr_t root,
    mpfr_t a, mpfr_t b, mpfr_t fa, mpfr_t fb,
    double tol_abs, double tol_step,
    int max_it, int use_threaded_z)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    mpfr_t z, fz, l, s, num, den, tmp;
    mpfr_inits2(P, z, fz, l, s, num, den, tmp, (mpfr_ptr)0);

    if (!opposite_sign(fa, fb))
    {
        mpfr_set(root, a, MPFR_RNDN);
        mpfr_clears(z, fz, l, s, num, den, tmp, (mpfr_ptr)0);
        return false;
    }

    bool ok = false;
    for (int it = 0; it < max_it; ++it)
    {
        // z = a - fa*(b-a)/(fb-fa)
        mpfr_sub(l, b, a, MPFR_RNDN);
        mpfr_sub(den, fb, fa, MPFR_RNDN);
        if (mpfr_zero_p(den))
            break;
        mpfr_mul(num, fa, l, MPFR_RNDN);
        mpfr_div(s, num, den, MPFR_RNDN);
        mpfr_sub(z, a, s, MPFR_RNDN);

        compute_Z(fz, z);

        // termination checks
        mpfr_t absfz;
        mpfr_init2(absfz, P);
        mpfr_abs(absfz, fz, MPFR_RNDN);
        const double absfz_d = mpfr_get_d(absfz, MPFR_RNDN);
        mpfr_clear(absfz);
        if (absfz_d <= tol_abs)
        {
            mpfr_set(root, z, MPFR_RNDN);
            ok = true;
            break;
        }

        mpfr_t absl;
        mpfr_init2(absl, P);
        mpfr_abs(absl, l, MPFR_RNDN);
        const double width = mpfr_get_d(absl, MPFR_RNDN);
        mpfr_clear(absl);
        if (width <= tol_step)
        {
            mpfr_set(root, z, MPFR_RNDN);
            ok = true;
            break;
        }

        // update bracket – scale the UNCHANGED endpoint (Pegasus); Illinois fallback on 0 denom
        if (opposite_sign(fz, fb))
        {
            // root in [z,b] -> shift a=z, scale fa
            mpfr_add(tmp, fb, fz, MPFR_RNDN);
            if (!mpfr_zero_p(tmp))
            {
                mpfr_div(num, fb, tmp, MPFR_RNDN); // m = fb/(fb+fz)
                mpfr_mul(fa, fa, num, MPFR_RNDN);  // fa = m*fa
            }
            else
            {
                mpfr_div_ui(fa, fa, 2, MPFR_RNDN); // Illinois
            }
            mpfr_set(a, z, MPFR_RNDN);
            mpfr_set(fa, fz, MPFR_RNDN);
        }
        else
        {
            // root in [a,z] -> shift b=z, scale fb
            mpfr_add(tmp, fa, fz, MPFR_RNDN);
            if (!mpfr_zero_p(tmp))
            {
                mpfr_div(num, fa, tmp, MPFR_RNDN); // m = fa/(fa+fz)
                mpfr_mul(fb, fb, num, MPFR_RNDN);  // fb = m*fb
            }
            else
            {
                mpfr_div_ui(fb, fb, 2, MPFR_RNDN); // Illinois
            }
            mpfr_set(b, z, MPFR_RNDN);
            mpfr_set(fb, fz, MPFR_RNDN);
        }

        mpfr_set(root, z, MPFR_RNDN); // best-so-far
    }

    mpfr_clears(z, fz, l, s, num, den, tmp, (mpfr_ptr)0);
    return ok;
}

/* ==================== derivatives ==================== */

// 3-point central: Z'(t) ≈ (Z(t+h) - Z(t-h)) / (2h)
static inline void diff_Z_3pt(mpfr_t dZ, const mpfr_t t, double h, int use_threaded_z)
{
    mpfr_t tp, tm, Zp, Zm;
    mpfr_inits2(CFG().mpfr_prec, tp, tm, Zp, Zm, (mpfr_ptr)0);
    mpfr_set_d(tp, h, MPFR_RNDN);
    mpfr_add(tp, t, tp, MPFR_RNDN);
    mpfr_set_d(tm, -h, MPFR_RNDN);
    mpfr_add(tm, t, tm, MPFR_RNDN);

    compute_Z(Zp, tp);
    compute_Z(Zm, tm);

    mpfr_t num;
    mpfr_init2(num, CFG().mpfr_prec);
    mpfr_sub(num, Zp, Zm, MPFR_RNDN);
    mpfr_div_d(dZ, num, 2.0 * h, MPFR_RNDN);
    mpfr_clears(tp, tm, Zp, Zm, num, (mpfr_ptr)0);
}

// 5-point: Z'(t) ≈ [ -Z(t+2h)+8Z(t+h)-8Z(t-h)+Z(t-2h) ]/(12h)
static inline void diff_Z_5pt(mpfr_t dZ, const mpfr_t t, double h, int use_threaded_z)
{
    mpfr_t tp2, tp1, tm1, tm2, Zp2, Zp1, Zm1, Zm2;
    mpfr_inits2(CFG().mpfr_prec, tp2, tp1, tm1, tm2, Zp2, Zp1, Zm1, Zm2, (mpfr_ptr)0);

    mpfr_set_d(tp2, 2 * h, MPFR_RNDN);
    mpfr_add(tp2, t, tp2, MPFR_RNDN);
    mpfr_set_d(tp1, h, MPFR_RNDN);
    mpfr_add(tp1, t, tp1, MPFR_RNDN);
    mpfr_set_d(tm1, -h, MPFR_RNDN);
    mpfr_add(tm1, t, tm1, MPFR_RNDN);
    mpfr_set_d(tm2, -2 * h, MPFR_RNDN);
    mpfr_add(tm2, t, tm2, MPFR_RNDN);

    compute_Z(Zp2, tp2);
    compute_Z(Zp1, tp1);
    compute_Z(Zm1, tm1);
    compute_Z(Zm2, tm2);

    mpfr_t num, term;
    mpfr_inits2(CFG().mpfr_prec, num, term, (mpfr_ptr)0);
    mpfr_neg(num, Zp2, MPFR_RNDN);
    mpfr_mul_ui(term, Zp1, 8, MPFR_RNDN);
    mpfr_add(num, num, term, MPFR_RNDN);
    mpfr_mul_ui(term, Zm1, 8, MPFR_RNDN);
    mpfr_sub(num, num, term, MPFR_RNDN);
    mpfr_add(num, num, Zm2, MPFR_RNDN);
    mpfr_div_d(dZ, num, 12.0 * h, MPFR_RNDN);

    mpfr_clears(tp2, tp1, tm1, tm2, Zp2, Zp1, Zm1, Zm2, num, term, (mpfr_ptr)0);
}

/* ==================== step control ==================== */

// Backtrack on Δ; return false if step shrinks too much without improvement
static inline bool backtrack_step(
    mpfr_t t,
    const mpfr_t Zcur,
    const mpfr_t delta,
    int use_threaded_z,
    double min_fraction = 1.0 / 64.0,
    double min_abs_rel = 1e-18)
{
    mpfr_t trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled;
    mpfr_inits2(CFG().mpfr_prec, trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled, (mpfr_ptr)0);

    mpfr_set_ui(frac, 1, MPFR_RNDN);
    double frac_d = 1.0;

    mpfr_abs(abs_cur, Zcur, MPFR_RNDN);

    bool ok = false;
    for (int k = 0; k < 12; ++k)
    {
        mpfr_mul(scaled, delta, frac, MPFR_RNDN);
        mpfr_abs(abs_scaled, scaled, MPFR_RNDN);

        const double tiny = mpfr_get_d(abs_scaled, MPFR_RNDN);
        const double tmag = std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
        if (tiny < min_abs_rel * tmag)
            break;

        mpfr_sub(trial_t, t, scaled, MPFR_RNDN);
        compute_Z(trial_Z, trial_t);

        mpfr_t abs_trial;
        mpfr_init2(abs_trial, CFG().mpfr_prec);
        mpfr_abs(abs_trial, trial_Z, MPFR_RNDN);

        if (mpfr_cmp(abs_trial, abs_cur) < 0)
        {
            mpfr_set(t, trial_t, MPFR_RNDN);
            ok = true;
            mpfr_clear(abs_trial);
            break;
        }
        mpfr_clear(abs_trial);

        frac_d *= 0.5;
        if (frac_d < min_fraction)
            break;
        mpfr_div_ui(frac, frac, 2, MPFR_RNDN);
    }

    mpfr_clears(trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled, (mpfr_ptr)0);
    return ok;
}

// If stuck (dZ ~ 0 and backtracking failed), try an adaptive nudge:
// sample t ± h0 and move to the side with smaller |Z|.
static inline bool descent_nudge(mpfr_t t, int use_threaded_z)
{
    const double s = est_spacing(t);
    const double h0 = std::max(1e-12, std::min(0.2 * s, choose_h(t)));

    mpfr_t tp, tm, Zp, Zm;
    mpfr_inits2(CFG().mpfr_prec, tp, tm, Zp, Zm, (mpfr_ptr)0);

    mpfr_set_d(tp, h0, MPFR_RNDN);
    mpfr_add(tp, t, tp, MPFR_RNDN);
    mpfr_set_d(tm, -h0, MPFR_RNDN);
    mpfr_add(tm, t, tm, MPFR_RNDN);

    compute_Z(Zp, tp);
    compute_Z(Zm, tm);

    mpfr_t aP, aM;
    mpfr_inits2(CFG().mpfr_prec, aP, aM, (mpfr_ptr)0);
    mpfr_abs(aP, Zp, MPFR_RNDN);
    mpfr_abs(aM, Zm, MPFR_RNDN);

    bool moved = false;
    if (mpfr_cmp(aP, aM) < 0)
    {
        mpfr_set(t, tp, MPFR_RNDN);
        moved = true;
    }
    else if (mpfr_cmp(aM, aP) < 0)
    {
        mpfr_set(t, tm, MPFR_RNDN);
        moved = true;
    }

    mpfr_clears(tp, tm, Zp, Zm, aP, aM, (mpfr_ptr)0);
    return moved;
}

/* ==================== Newton driver ==================== */

bool newton_refine_root(
    mpfr_t t_out,
    const mpfr_t t0,
    double tol_abs,
    double tol_step,
    int max_iter,
    int use_threaded_z)
{
    mpfr_t t, Z, dZ, delta;
    mpfr_inits2(CFG().mpfr_prec, t, Z, dZ, delta, (mpfr_ptr)0);
    mpfr_set(t, t0, MPFR_RNDN);

    // secant safeguard state
    mpfr_t t_prev, Z_prev;
    mpfr_inits2(CFG().mpfr_prec, t_prev, Z_prev, (mpfr_ptr)0);
    bool have_prev = false;

    // optional bracket state (for Pegasus burst)
    mpfr_t br_a, br_b, br_fa, br_fb;
    bool have_bracket = false;
    mpfr_inits2(CFG().mpfr_prec, br_a, br_b, br_fa, br_fb, (mpfr_ptr)0);

    for (int it = 0; it < max_iter; ++it)
    {
        LOGMP_DEBUG("-> Newton single iter = %d", it);
        compute_Z(Z, t);

        // keep/update a bracket when we see a sign change
        if (have_prev && opposite_sign(Z_prev, Z))
        {
            if (!have_bracket)
            {
                mpfr_set(br_a, t_prev, MPFR_RNDN);
                mpfr_set(br_b, t, MPFR_RNDN);
                mpfr_set(br_fa, Z_prev, MPFR_RNDN);
                mpfr_set(br_fb, Z, MPFR_RNDN);
                have_bracket = true;
            }
            else
            {
                mpfr_t new_w, old_w;
                mpfr_inits2(CFG().mpfr_prec, new_w, old_w, (mpfr_ptr)0);
                mpfr_sub(new_w, t, t_prev, MPFR_RNDN);
                mpfr_abs(new_w, new_w, MPFR_RNDN);
                mpfr_sub(old_w, br_b, br_a, MPFR_RNDN);
                mpfr_abs(old_w, old_w, MPFR_RNDN);
                if (mpfr_cmp(new_w, old_w) < 0)
                {
                    mpfr_set(br_a, t_prev, MPFR_RNDN);
                    mpfr_set(br_b, t, MPFR_RNDN);
                    mpfr_set(br_fa, Z_prev, MPFR_RNDN);
                    mpfr_set(br_fb, Z, MPFR_RNDN);
                }
                mpfr_clears(new_w, old_w, (mpfr_ptr)0);
            }
        }

        // stop if |Z| small enough
        mpfr_t absZ;
        mpfr_init2(absZ, CFG().mpfr_prec);
        mpfr_abs(absZ, Z, MPFR_RNDN);
        const double absZ_d = mpfr_get_d(absZ, MPFR_RNDN);
        mpfr_clear(absZ);
        if (absZ_d <= tol_abs)
        {
            mpfr_set(t_out, t, MPFR_RNDN);
            mpfr_clears(br_a, br_b, br_fa, br_fb, (mpfr_ptr)0);
            mpfr_clears(t, Z, dZ, delta, t_prev, Z_prev, (mpfr_ptr)0);
            return true;
        }

        // derivative: 3pt then 5pt if needed
        double h = choose_h(t);
        diff_Z_3pt(dZ, t, h, use_threaded_z);

        bool dZ_ok = !mpfr_zero_p(dZ);
        if (dZ_ok)
        {
            mpfr_t absd;
            mpfr_init2(absd, CFG().mpfr_prec);
            mpfr_abs(absd, dZ, MPFR_RNDN);
            if (mpfr_cmp_d(absd, 1e-30) < 0)
                dZ_ok = false;
            mpfr_clear(absd);
        }
        if (!dZ_ok)
        {
            const double h_big = std::min(2.0 * h, 0.2 * est_spacing(t));
            diff_Z_5pt(dZ, t, h_big, use_threaded_z);
            dZ_ok = !mpfr_zero_p(dZ);
        }

        // delta = Z/dZ or secant fallback
        if (!dZ_ok)
        {
            if (have_prev)
            {
                mpfr_t num, den;
                mpfr_inits2(CFG().mpfr_prec, num, den, (mpfr_ptr)0);
                mpfr_sub(den, Z, Z_prev, MPFR_RNDN);
                if (mpfr_zero_p(den))
                {
                    mpfr_clears(num, den, (mpfr_ptr)0);
                    break;
                }
                mpfr_sub(num, t, t_prev, MPFR_RNDN);
                mpfr_mul(num, num, Z, MPFR_RNDN);
                mpfr_div(delta, num, den, MPFR_RNDN);
                mpfr_clears(num, den, (mpfr_ptr)0);
            }
            else
            {
                if (!descent_nudge(t, use_threaded_z))
                {
                    mpfr_add_d(t, t, 1e-9, MPFR_RNDN);
                }
                continue;
            }
        }
        else
        {
            mpfr_div(delta, Z, dZ, MPFR_RNDN);
        }

        // clamp step within one valley
        const double s = est_spacing(t);
        const double step_max = 0.7 * s;
        double d_abs = std::fabs(mpfr_get_d(delta, MPFR_RNDN));
        if (d_abs > step_max)
        {
            const double scale = step_max / d_abs;
            mpfr_mul_d(delta, delta, scale, MPFR_RNDN);
            d_abs = step_max;
        }

        // accept if already tiny step
        if (d_abs <= tol_step)
        {
            mpfr_set(t_out, t, MPFR_RNDN);
            mpfr_clears(br_a, br_b, br_fa, br_fb, (mpfr_ptr)0);
            mpfr_clears(t, Z, dZ, delta, t_prev, Z_prev, (mpfr_ptr)0);
            return true;
        }

        // save state for secant
        mpfr_set(t_prev, t, MPFR_RNDN);
        mpfr_set(Z_prev, Z, MPFR_RNDN);
        have_prev = true;

        // backtracking; if it fails and we have a bracket, try a tiny Pegasus burst
        if (!backtrack_step(t, Z, delta, use_threaded_z))
        {
            if (have_bracket)
            {
                mpfr_t z_br;
                mpfr_init2(z_br, CFG().mpfr_prec);
                const bool ok = pegasus_refine_bracket(
                    z_br, br_a, br_b, br_fa, br_fb,
                    /*tol_abs*/ tol_abs, /*tol_step*/ stage_step_tol(t, tol_abs),
                    /*max_it*/ 8, use_threaded_z);
                if (ok)
                {
                    mpfr_set(t, z_br, MPFR_RNDN);
                    mpfr_clear(z_br);
                    continue;
                }
                mpfr_clear(z_br);
            }
            if (!descent_nudge(t, use_threaded_z))
            {
                mpfr_sub(t, t, delta, MPFR_RNDN);
            }
        }
    }

    // no convergence: return last t
    mpfr_set(t_out, t, MPFR_RNDN);
    mpfr_clears(br_a, br_b, br_fa, br_fb, (mpfr_ptr)0);
    mpfr_clears(t, Z, dZ, delta, t_prev, Z_prev, (mpfr_ptr)0);
    return false;
}
