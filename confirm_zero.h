// confirm_zero.h
#pragma once
#include <mpfr.h>

// Estimate local spacing s ≈ π / log(t/2π)
static inline double est_spacing(const mpfr_t t){
    mpfr_t pi, two_pi, x, logx, s;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, x, logx, s, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(x, t, two_pi, MPFR_RNDN);
    // clamp x≥1e-12
    if (mpfr_cmp_d(x, 1e-12) < 0) mpfr_set_d(x, 1e-12, MPFR_RNDN);
    mpfr_log(logx, x, MPFR_RNDN);
    if (mpfr_cmp_d(logx, 1e-12) <= 0){ mpfr_clears(pi,two_pi,x,logx,s,(mpfr_ptr)0); return 0.25; }
    mpfr_div(s, pi, logx, MPFR_RNDN);
    double out = mpfr_get_d(s, MPFR_RNDN);
    mpfr_clears(pi,two_pi,x,logx,s,(mpfr_ptr)0);
    return out;
}

// extern Z evaluators you already have:
extern "C++" void compute_Z_small(mpfr_t Z, const mpfr_t t);
extern "C++" void compute_Z_threaded(mpfr_t Z, mpfr_t t);

static inline bool use_threaded_for_t(const mpfr_t t) {
    mpfr_t pi, two_pi, tmp; mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
    return (N >= 8);
}


// Return true only if you see a **sign change** at ±ε
static inline bool confirm_zero_like(const mpfr_t t, int use_threaded_z){
    double s = est_spacing(t);
    double eps = 0.2 * s;  // check 20% of expected spacing

    mpfr_t tp, tm, Zp, Zm;
    mpfr_inits2(MPFR_PRECISION, tp, tm, Zp, Zm, (mpfr_ptr)0);

    mpfr_add_d(tp, t,  eps, MPFR_RNDN);
    mpfr_add_d(tm, t, -eps, MPFR_RNDN);

    if (use_threaded_z && use_threaded_for_t(t)) {
        mpfr_t tpc, tmc; mpfr_inits2(MPFR_PRECISION, tpc, tmc, (mpfr_ptr)0);
        mpfr_set(tpc, tp, MPFR_RNDN);
        mpfr_set(tmc, tm, MPFR_RNDN);
        compute_Z_threaded(Zp, tpc);
        compute_Z_threaded(Zm, tmc);
        mpfr_clears(tpc, tmc, (mpfr_ptr)0);
    } else {
        compute_Z_small(Zp, tp);
        compute_Z_small(Zm, tm);
    }

    int sp = mpfr_sgn(Zp);
    int sm = mpfr_sgn(Zm);

    mpfr_clears(tp, tm, Zp, Zm, (mpfr_ptr)0);
    return (sp * sm) < 0;  // opposite signs → very likely a real zero nearby
}
