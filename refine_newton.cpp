#include "refine_newton.h"
#include "config.h"
//#include "z_calculator_thread.h"
//#include "z_calculator_small.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

extern "C++" {
        void compute_Z_small(mpfr_t Z, const mpfr_t t);
        void compute_Z_threaded(mpfr_t Z, mpfr_t t);
}

// --- Helpers ---------------------------------------------------------

// Decide whether to use threaded RS or small-Z at this t
static inline bool use_threaded_for_t(const mpfr_t t) {
    // N = floor(sqrt(t/2π)) — if N>=8 we prefer the threaded RS path
    mpfr_t pi, two_pi, tmp; 
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
    
    return (N >= 8);
}

// Evaluate Z(t) using the best engine for this t and setting
static inline void eval_Z(mpfr_t Z, mpfr_t t, int use_threaded_z) {
    
        if (use_threaded_z && use_threaded_for_t(t)) {
                // compute_Z_threaded takes non-const mpfr_t for t → copy
                mpfr_t tc; mpfr_init2(tc, MPFR_PRECISION);
                mpfr_set(tc, t, MPFR_RNDN);
                compute_Z_threaded(Z, tc);
                mpfr_clear(tc);
        } else {
                compute_Z_small(Z, t);
        }
}

// Estimate local zero spacing ~ π / log(t/2π) (used to clamp Newton steps)
static inline double est_spacing(const mpfr_t t) {
    mpfr_t pi, two_pi, x, logx, s;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, x, logx, s, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);     // 2π
    mpfr_div(x, t, two_pi, MPFR_RNDN);         // x = t / (2π)

    // Clamp x from below to avoid log(<=0) / tiny logs
    mpfr_t minx; mpfr_init2(minx, MPFR_PRECISION);
    mpfr_set_d(minx, 1e-12, MPFR_RNDN);
    if (mpfr_cmp(x, minx) < 0) mpfr_set(x, minx, MPFR_RNDN);

    mpfr_log(logx, x, MPFR_RNDN);              // log(t/2π)
    // If log too small/nonpositive, return safe cap
    if (mpfr_cmp_d(logx, 1e-12) <= 0) {
        mpfr_clears(pi, two_pi, x, logx, s, minx, (mpfr_ptr)0);
        return 0.25;
    }

    mpfr_div(s, pi, logx, MPFR_RNDN);          // s = π / log(t/2π)
    double out = mpfr_get_d(s, MPFR_RNDN);

    mpfr_clears(pi, two_pi, x, logx, s, minx, (mpfr_ptr)0);
    return out;
}

// Choose a finite-difference step for derivative of Z
static inline double choose_h(const mpfr_t t) {
    double s = est_spacing(t);
    double eps = std::ldexp(1.0, -(int)MPFR_PRECISION);
    double h_candidate = std::pow(eps, 1.0/3.0) * std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
    double h = std::clamp(h_candidate, 1e-9, 0.1 * s);
    return h;
}

// 5-point symmetric derivative of Z at t:
// Z'(t) ≈ [ -Z(t+2h) + 8 Z(t+h) - 8 Z(t-h) + Z(t-2h) ] / (12 h)
static inline void diff_Z(
        mpfr_t dZ, 
        const mpfr_t t, 
        double h, 
        int use_threaded_z) {

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
    mpfr_neg(num, Zp2, MPFR_RNDN);           // -Z(t+2h)
    mpfr_mul_ui(term, Zp1, 8, MPFR_RNDN);    // +8 Z(t+h)
    mpfr_add(num, num, term, MPFR_RNDN);
    mpfr_mul_ui(term, Zm1, 8, MPFR_RNDN);    // -8 Z(t-h)
    mpfr_sub(num, num, term, MPFR_RNDN);
    mpfr_add(num, num, Zm2, MPFR_RNDN);      // + Z(t-2h)

    mpfr_div_d(dZ, num, 12.0 * h, MPFR_RNDN);

    mpfr_clears(tp2,tp1,tm1,tm2,Zp2,Zp1,Zm1,Zm2,num,term,(mpfr_ptr)0);
}

// Backtracking line search on Δt (halve until |Z| decreases or min step)
static inline bool backtrack_step(
        mpfr_t t, 
        const mpfr_t Zcur,
        const mpfr_t delta, 
        int use_threaded_z,
        double min_fraction = 1.0/64.0){

    mpfr_t trial_t, trial_Z, scaled, frac;
    mpfr_inits2(MPFR_PRECISION, trial_t, trial_Z, scaled, frac, (mpfr_ptr)0);

    mpfr_set_ui(frac, 1, MPFR_RNDN);             // λ = 1
    double frac_d = 1.0;

    // compute baseline |Zcur|
    mpfr_t abs_cur; mpfr_init2(abs_cur, MPFR_PRECISION);
    mpfr_abs(abs_cur, Zcur, MPFR_RNDN);

    bool ok = false;
    for (int k = 0; k < 8; ++k) { // up to 8 halvings
        // scaled = λ * delta
        mpfr_mul(scaled, delta, frac, MPFR_RNDN);
        // trial_t = t - scaled
        mpfr_sub(trial_t, t, scaled, MPFR_RNDN);
        // Z(trial_t)
        eval_Z(trial_Z, trial_t, use_threaded_z);

        mpfr_t abs_trial; mpfr_init2(abs_trial, MPFR_PRECISION);
        mpfr_abs(abs_trial, trial_Z, MPFR_RNDN);

        if (mpfr_cmp(abs_trial, abs_cur) < 0) {
            // accept
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

    mpfr_clears(trial_t, trial_Z, scaled, frac, abs_cur, (mpfr_ptr)0);
    return ok;
}

// --- Newton driver ---------------------------------------------------

bool newton_refine_root(
        mpfr_t t_out, 
        const mpfr_t t0,
        double tol_abs, 
        double tol_step,
        int max_iter, 
        int use_threaded_z){
    
    mpfr_t t, Z, dZ, delta;
    mpfr_inits2(MPFR_PRECISION, t, Z, dZ, delta, (mpfr_ptr)0);
    mpfr_set(t, t0, MPFR_RNDN);

    // For secant fallback
    mpfr_t t_prev, Z_prev; 
    mpfr_inits2(MPFR_PRECISION, t_prev, Z_prev, (mpfr_ptr)0);
    
    bool have_prev = false;

    for (int it = 0; it < max_iter; ++it) {
        // Z(t)
        eval_Z(Z, t, use_threaded_z);

        // |Z| stop
        mpfr_t absZ; 
        mpfr_init2(absZ, MPFR_PRECISION);

        mpfr_abs(absZ, Z, MPFR_RNDN);
        double absZ_d = mpfr_get_d(absZ, MPFR_RNDN);
        mpfr_clear(absZ);
        
        if (absZ_d <= tol_abs) { 
                mpfr_set(t_out, t, MPFR_RNDN); 
                mpfr_clears(t,Z,dZ,delta,t_prev,Z_prev,(mpfr_ptr)0); 
                return true; 
        }

        // derivative via 5-pt stencil
        double h = choose_h(t);
        diff_Z(dZ, t, h, use_threaded_z);

        // delta = Z / dZ
        if (mpfr_zero_p(dZ)) {
            // fallback to secant if we have previous
            if (have_prev) {
                // delta = Z * (t - t_prev) / (Z - Z_prev)
                mpfr_t num, den; 
                mpfr_inits2(MPFR_PRECISION, num, den, (mpfr_ptr)0);
                
                mpfr_sub(den, Z, Z_prev, MPFR_RNDN);
                
                if (mpfr_zero_p(den)) { 
                        mpfr_clears(num,den,(mpfr_ptr)0); 
                        break; 
                }

                mpfr_sub(num, t, t_prev, MPFR_RNDN);
                mpfr_mul(num, num, Z, MPFR_RNDN);
                mpfr_div(delta, num, den, MPFR_RNDN);
                mpfr_clears(num, den, (mpfr_ptr)0);
            } else {
                // nudge t by a tiny step and continue
                mpfr_add_d(t, t, 1e-6, MPFR_RNDN);
                continue;
            }
        } else {
            mpfr_div(delta, Z, dZ, MPFR_RNDN);
        }

        // step clamp using local spacing
        double s = est_spacing(t);
        double step_max = 0.7 * s; // keep it within one valley
        double d_d = std::fabs(mpfr_get_d(delta, MPFR_RNDN));
        if (d_d > step_max) {
            double scale = step_max / d_d;
            mpfr_mul_d(delta, delta, scale, MPFR_RNDN);
        }

        // |Δt| stop
        if (std::fabs(mpfr_get_d(delta, MPFR_RNDN)) <= tol_step) {
            mpfr_set(t_out, t, MPFR_RNDN);
            mpfr_clears(t,Z,dZ,delta,t_prev,Z_prev,(mpfr_ptr)0);
            return true;
        }

        // Keep previous for possible secant fallback
        mpfr_set(t_prev, t, MPFR_RNDN);
        mpfr_set(Z_prev, Z, MPFR_RNDN);
        have_prev = true;

        // backtracking line search on Δ to guarantee descent in |Z|
        bool ok = backtrack_step(t, Z, delta, use_threaded_z);
        if (!ok) {
            // last resort: take the damped step anyway
            mpfr_sub(t, t, delta, MPFR_RNDN);
        }
    }

    // No convergence within max_iter; return last iterate
    mpfr_set(t_out, t, MPFR_RNDN);
    mpfr_clears(t, Z, dZ, delta, t_prev, Z_prev, (mpfr_ptr)0);

    return false;
}
