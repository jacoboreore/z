#include <cmath>
#include <algorithm>
#include "refine_newton.h"
#include "config.h"
#include "confirm_zero.h"

extern "C++" {
    void compute_Z_small(mpfr_t Z, const mpfr_t t);
    void compute_Z_threaded(mpfr_t Z, mpfr_t t);
}

/* ==================== Utilities ==================== */

// Pick engine based on N = floor(sqrt(t/2π))
// static inline bool use_threaded_for_t(const mpfr_t t) {
//     mpfr_t pi, two_pi, tmp; mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
//     mpfr_const_pi(pi, MPFR_RNDN);
//     mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
//     mpfr_div(tmp, t, two_pi, MPFR_RNDN);
//     mpfr_sqrt(tmp, tmp, MPFR_RNDN);
//     long N = mpfr_get_si(tmp, MPFR_RNDN);
//     mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
//     return (N >= 8);
// }

static inline void eval_Z(mpfr_t Z, mpfr_t t, int use_threaded_z) {
    if (use_threaded_z && use_threaded_for_t(t)) {
        mpfr_t tc; mpfr_init2(tc, MPFR_PRECISION);
        mpfr_set(tc, t, MPFR_RNDN);
        compute_Z_threaded(Z, tc);
        mpfr_clear(tc);
    } else {
        compute_Z_small(Z, t);
    }
}

// π/log(t/2π) with protection
// static inline double est_spacing(const mpfr_t t) {
//     mpfr_t pi, two_pi, x, logx, s;
//     mpfr_inits2(MPFR_PRECISION, pi, two_pi, x, logx, s, (mpfr_ptr)0);
//     mpfr_const_pi(pi, MPFR_RNDN);
//     mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
//     mpfr_div(x, t, two_pi, MPFR_RNDN);
//     mpfr_t minx; mpfr_init2(minx, MPFR_PRECISION); mpfr_set_d(minx, 1e-12, MPFR_RNDN);
//     if (mpfr_cmp(x, minx) < 0) mpfr_set(x, minx, MPFR_RNDN);
//     mpfr_log(logx, x, MPFR_RNDN);
//     if (mpfr_cmp_d(logx, 1e-12) <= 0) { mpfr_clears(pi,two_pi,x,logx,s,minx,(mpfr_ptr)0); return 0.25; }
//     mpfr_div(s, pi, logx, MPFR_RNDN);
//     double out = mpfr_get_d(s, MPFR_RNDN);
//     mpfr_clears(pi,two_pi,x,logx,s,minx,(mpfr_ptr)0);
//     return out;
// }

// Choose h for derivative (scale with spacing + precision)
static inline double choose_h(const mpfr_t t) {
    const double s   = est_spacing(t);
    const double eps = std::ldexp(1.0, -(int)MPFR_PRECISION);   // 2^(-P)
    const double h0  = std::pow(eps, 1.0/3.0) * std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
    return std::clamp(h0, 1e-9, 0.1 * s);
}

/* ==================== Derivative ==================== */

// 3-point central: Z'(t) ≈ (Z(t+h) - Z(t-h)) / (2h)
static inline void diff_Z_3pt(mpfr_t dZ, const mpfr_t t, double h, int use_threaded_z) {
    mpfr_t tp, tm, Zp, Zm; mpfr_inits2(MPFR_PRECISION, tp, tm, Zp, Zm, (mpfr_ptr)0);
    mpfr_set_d(tp,  h, MPFR_RNDN); mpfr_add(tp, t, tp, MPFR_RNDN);
    mpfr_set_d(tm, -h, MPFR_RNDN); mpfr_add(tm, t, tm, MPFR_RNDN);
    eval_Z(Zp, tp, use_threaded_z);
    eval_Z(Zm, tm, use_threaded_z);
    mpfr_t num; mpfr_init2(num, MPFR_PRECISION);
    mpfr_sub(num, Zp, Zm, MPFR_RNDN);
    mpfr_div_d(dZ, num, 2.0*h, MPFR_RNDN);
    mpfr_clears(tp,tm,Zp,Zm,num,(mpfr_ptr)0);
}

// 5-point: Z'(t) ≈ [ -Z(t+2h)+8Z(t+h)-8Z(t-h)+Z(t-2h) ]/(12h)
static inline void diff_Z_5pt(mpfr_t dZ, const mpfr_t t, double h, int use_threaded_z) {
    mpfr_t tp2, tp1, tm1, tm2, Zp2, Zp1, Zm1, Zm2;
    mpfr_inits2(MPFR_PRECISION, tp2, tp1, tm1, tm2, Zp2, Zp1, Zm1, Zm2, (mpfr_ptr)0);
    mpfr_set_d(tp2, 2*h, MPFR_RNDN); mpfr_add(tp2, t, tp2, MPFR_RNDN);
    mpfr_set_d(tp1,   h, MPFR_RNDN); mpfr_add(tp1, t, tp1, MPFR_RNDN);
    mpfr_set_d(tm1,  -h, MPFR_RNDN); mpfr_add(tm1, t, tm1, MPFR_RNDN);
    mpfr_set_d(tm2, -2*h, MPFR_RNDN); mpfr_add(tm2, t, tm2, MPFR_RNDN);
    eval_Z(Zp2, tp2, use_threaded_z);
    eval_Z(Zp1, tp1, use_threaded_z);
    eval_Z(Zm1, tm1, use_threaded_z);
    eval_Z(Zm2, tm2, use_threaded_z);
    mpfr_t num, term; mpfr_inits2(MPFR_PRECISION, num, term, (mpfr_ptr)0);
    mpfr_neg(num, Zp2, MPFR_RNDN);
    mpfr_mul_ui(term, Zp1, 8, MPFR_RNDN); mpfr_add(num, num, term, MPFR_RNDN);
    mpfr_mul_ui(term, Zm1, 8, MPFR_RNDN); mpfr_sub(num, num, term, MPFR_RNDN);
    mpfr_add(num, num, Zm2, MPFR_RNDN);
    mpfr_div_d(dZ, num, 12.0*h, MPFR_RNDN);
    mpfr_clears(tp2,tp1,tm1,tm2,Zp2,Zp1,Zm1,Zm2,num,term,(mpfr_ptr)0);
}

/* ==================== Step control ==================== */

// Backtrack on Δ; return false if step shrinks too much without improvement
static inline bool backtrack_step(
    mpfr_t t,
    const mpfr_t Zcur,
    const mpfr_t delta,
    int use_threaded_z,
    double min_fraction = 1.0/64.0,
    double min_abs_rel  = 1e-18)
{
    mpfr_t trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled;
    mpfr_inits2(MPFR_PRECISION, trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled, (mpfr_ptr)0);

    mpfr_set_ui(frac, 1, MPFR_RNDN);  // λ
    double frac_d = 1.0;

    mpfr_abs(abs_cur, Zcur, MPFR_RNDN);

    bool ok = false;
    for (int k = 0; k < 12; ++k) {
        // scaled = λ * delta ; trial_t = t - scaled
        mpfr_mul(scaled, delta, frac, MPFR_RNDN);
        mpfr_abs(abs_scaled, scaled, MPFR_RNDN);

        // too tiny absolute movement? bail
        double tiny = mpfr_get_d(abs_scaled, MPFR_RNDN);
        const double tmag = std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
        if (tiny < min_abs_rel * tmag) break;

        mpfr_sub(trial_t, t, scaled, MPFR_RNDN);
        eval_Z(trial_Z, trial_t, use_threaded_z);

        mpfr_t abs_trial; mpfr_init2(abs_trial, MPFR_PRECISION);
        mpfr_abs(abs_trial, trial_Z, MPFR_RNDN);

        if (mpfr_cmp(abs_trial, abs_cur) < 0) {
            mpfr_set(t, trial_t, MPFR_RNDN);
            ok = true;
            mpfr_clear(abs_trial);
            break;
        }
        mpfr_clear(abs_trial);

        // halve λ
        frac_d *= 0.5;
        if (frac_d < min_fraction) break;
        mpfr_div_ui(frac, frac, 2, MPFR_RNDN);
    }

    mpfr_clears(trial_t, trial_Z, scaled, frac, abs_cur, abs_scaled, (mpfr_ptr)0);
    return ok;
}

// If we’re stuck (dZ ~ 0 and backtracking failed), try an adaptive nudge:
// sample t ± h0 and move to the side with smaller |Z|.
static inline bool descent_nudge(mpfr_t t, int use_threaded_z) {
    const double s  = est_spacing(t);
    const double h0 = std::max(1e-12, std::min(0.2*s, choose_h(t)));

    mpfr_t tp, tm, Zp, Zm; mpfr_inits2(MPFR_PRECISION, tp, tm, Zp, Zm, (mpfr_ptr)0);
    mpfr_set_d(tp,  h0, MPFR_RNDN); mpfr_add(tp, t, tp, MPFR_RNDN);
    mpfr_set_d(tm, -h0, MPFR_RNDN); mpfr_add(tm, t, tm, MPFR_RNDN);
    eval_Z(Zp, tp, use_threaded_z);
    eval_Z(Zm, tm, use_threaded_z);

    mpfr_t aP, aM; mpfr_inits2(MPFR_PRECISION, aP, aM, (mpfr_ptr)0);
    mpfr_abs(aP, Zp, MPFR_RNDN);
    mpfr_abs(aM, Zm, MPFR_RNDN);
    bool moved = false;
    if (mpfr_cmp(aP, aM) < 0) { mpfr_set(t, tp, MPFR_RNDN); moved = true; }
    else if (mpfr_cmp(aM, aP) < 0) { mpfr_set(t, tm, MPFR_RNDN); moved = true; }

    mpfr_clears(tp,tm,Zp,Zm,aP,aM,(mpfr_ptr)0);
    return moved;
}

/* ==================== Newton driver ==================== */

bool newton_refine_root(
    mpfr_t t_out,
    const mpfr_t t0,
    double tol_abs,
    double tol_step,
    int    max_iter,
    int    use_threaded_z)
{
    mpfr_t t, Z, dZ, delta;
    mpfr_inits2(MPFR_PRECISION, t, Z, dZ, delta, (mpfr_ptr)0);
    mpfr_set(t, t0, MPFR_RNDN);

    // Secant safeguard state
    mpfr_t t_prev, Z_prev; mpfr_inits2(MPFR_PRECISION, t_prev, Z_prev, (mpfr_ptr)0);
    bool have_prev = false;

    for (int it = 0; it < max_iter; ++it) {
        // Evaluate Z(t)
        eval_Z(Z, t, use_threaded_z);

        // Stop if |Z| small enough
        mpfr_t absZ; mpfr_init2(absZ, MPFR_PRECISION);
        mpfr_abs(absZ, Z, MPFR_RNDN);
        const double absZ_d = mpfr_get_d(absZ, MPFR_RNDN);
        mpfr_clear(absZ);
        if (absZ_d <= tol_abs) {
            mpfr_set(t_out, t, MPFR_RNDN);
            mpfr_clears(t,Z,dZ,delta,t_prev,Z_prev,(mpfr_ptr)0);
            return true;
        }

        // Derivative: start with 3pt; if too small/unstable, retry with 5pt and/or larger h
        double h = choose_h(t);
        diff_Z_3pt(dZ, t, h, use_threaded_z);

        bool dZ_ok = !mpfr_zero_p(dZ);
        if (dZ_ok) {
            // If |dZ| is *extremely* tiny, try to improve estimate
            mpfr_t absd; mpfr_init2(absd, MPFR_PRECISION); mpfr_abs(absd, dZ, MPFR_RNDN);
            if (mpfr_cmp_d(absd, 1e-30) < 0) dZ_ok = false;
            mpfr_clear(absd);
        }
        if (!dZ_ok) {
            // upgrade scheme and/or enlarge h a bit
            diff_Z_5pt(dZ, t, std::min(2.0*h, 0.2*est_spacing(t)), use_threaded_z);
            dZ_ok = !mpfr_zero_p(dZ);
        }

        // Compute Newton step delta = Z/dZ or fallbacks
        if (!dZ_ok) {
            if (have_prev) {
                // Secant delta = Z * (t - t_prev) / (Z - Z_prev)
                mpfr_t num, den; mpfr_inits2(MPFR_PRECISION, num, den, (mpfr_ptr)0);
                mpfr_sub(den, Z, Z_prev, MPFR_RNDN);
                if (mpfr_zero_p(den)) { mpfr_clears(num,den,(mpfr_ptr)0); break; }
                mpfr_sub(num, t, t_prev, MPFR_RNDN);
                mpfr_mul(num, num, Z, MPFR_RNDN);
                mpfr_div(delta, num, den, MPFR_RNDN);
                mpfr_clears(num, den, (mpfr_ptr)0);
            } else {
                // No derivative and no secant: try a descent nudge
                if (!descent_nudge(t, use_threaded_z)) {
                    // give a microscopic push and continue
                    mpfr_add_d(t, t, 1e-9, MPFR_RNDN);
                }
                continue; // iterate again
            }
        } else {
            mpfr_div(delta, Z, dZ, MPFR_RNDN);
        }

        // Clamp step size to stay within one valley
        const double s = est_spacing(t);
        const double step_max = 0.7 * s;
        const double d_d = std::fabs(mpfr_get_d(delta, MPFR_RNDN));
        if (d_d > step_max) {
            const double scale = step_max / d_d;
            mpfr_mul_d(delta, delta, scale, MPFR_RNDN);
        }

        // If step already tiny, accept
        if (std::fabs(mpfr_get_d(delta, MPFR_RNDN)) <= tol_step) {
            mpfr_set(t_out, t, MPFR_RNDN);
            mpfr_clears(t,Z,dZ,delta,t_prev,Z_prev,(mpfr_ptr)0);
            return true;
        }

        // Save state for secant fallback
        mpfr_set(t_prev, t, MPFR_RNDN);
        mpfr_set(Z_prev, Z, MPFR_RNDN);
        have_prev = true;

        // Backtracking on Δ to guarantee |Z| descent
        if (!backtrack_step(t, Z, delta, use_threaded_z)) {
            // Try a descent nudge (adaptive), then continue loop
            if (!descent_nudge(t, use_threaded_z)) {
                // As a final fallback, take the damped step anyway
                mpfr_sub(t, t, delta, MPFR_RNDN);
            }
        }
    }

    // No convergence: return last t
    mpfr_set(t_out, t, MPFR_RNDN);
    mpfr_clears(t, Z, dZ, delta, t_prev, Z_prev, (mpfr_ptr)0);
    return false;
}
