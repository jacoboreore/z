// z_ball_verbose.cpp
#include "z_ball.h"
#include "config.h"
#include <cmath>

extern "C++" {
    void compute_Z_small(mpfr_t Z, const mpfr_t t);
    void compute_Z_threaded(mpfr_t Z, mpfr_t t);
    void remainder_bound(mpfr_t R, const mpfr_t t);
}

// Helpers (same as before)
static inline bool use_threaded_for_t(const mpfr_t t){
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN); if (N<1) N=1;
    mpfr_clears(pi,two_pi,tmp,(mpfr_ptr)0);
    return (N >= 8);
}
static inline long rs_N(const mpfr_t t){
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN); if (N<1) N=1;
    mpfr_clears(pi,two_pi,tmp,(mpfr_ptr)0);
    return N;
}
static inline void rounding_radius(mpfr_t out, long N, double k=4.0){
    mpfr_prec_t P = MPFR_PRECISION;
    mpfr_t sN, ulp, kk;
    mpfr_inits2(P, sN, ulp, kk, (mpfr_ptr)0);
    mpfr_set_ui(sN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sN, sN, MPFR_RNDN);
    mpfr_mul_ui(sN, sN, 2, MPFR_RNDN); // 2*sqrt(N)
    mpfr_set_ui(ulp, 1, MPFR_RNDN);
    mpfr_div_2si(ulp, ulp, P, MPFR_RNDN); // 2^{-P}
    mpfr_set_d(kk, k, MPFR_RNDN);
    mpfr_mul(out, sN, ulp, MPFR_RNDN);
    mpfr_mul(out, out, kk, MPFR_RNDN);
    mpfr_clears(sN, ulp, kk, (mpfr_ptr)0);
}

// OPTIONAL: bound from theta truncation (difference big vs extended asymptotic)
extern "C++" void compute_theta_big(mpfr_t theta, const mpfr_t t);
extern "C++" void compute_theta_extended(mpfr_t theta, const mpfr_t t); // add 1/t^3,...,1/t^9

static inline void theta_radius(mpfr_t out, const mpfr_t t){
    // rad_theta ≈ |Δθ| * 2*sqrt(N)
    long N = rs_N(t);
    mpfr_t th1, th2, dth, factor;
    mpfr_inits2(MPFR_PRECISION, th1, th2, dth, factor, (mpfr_ptr)0);
    compute_theta_big(th1, t);
    compute_theta_extended(th2, t);  // your “exact-asymptotic” up to 1/t^9
    mpfr_sub(dth, th1, th2, MPFR_RNDN);
    mpfr_abs(dth, dth, MPFR_RNDN);
    mpfr_set_ui(factor, (unsigned long)(2.0*std::sqrt((double)N)), MPFR_RNDN);
    mpfr_mul(out, dth, factor, MPFR_RNDN);
    mpfr_clears(th1, th2, dth, factor, (mpfr_ptr)0);
}

void compute_Z_ball_verbose(ZBallV& B, const mpfr_t t){
    // midpoint
    if (use_threaded_for_t(t)) {
        mpfr_t tc; mpfr_init2(tc, MPFR_PRECISION);
        mpfr_set(tc, t, MPFR_RNDN);
        compute_Z_threaded(B.mid, tc);
        mpfr_clear(tc);
    } else {
        compute_Z_small(B.mid, t);
    }

    // components
    remainder_bound(B.rad_R, t);
    rounding_radius(B.rad_round, rs_N(t), /*k=*/4.0);

    // If you have compute_theta_extended implemented, uncomment:
    // theta_radius(B.rad_theta, t);
    // else keep it at 0.

    // total radius
    mpfr_set(B.rad_tot, B.rad_R, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_round, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_theta, MPFR_RNDN);
}
