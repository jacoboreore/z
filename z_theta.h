#pragma once
#include <mpfr.h>

// forward decls of the two implementations
void compute_theta_via_lggamma(mpfr_t theta, const mpfr_t t);
void compute_theta(mpfr_t theta, const mpfr_t t);

// Provide the FULL definition here (inline, not static)
inline void compute_theta_big(mpfr_t theta, const mpfr_t t) {
    
    // (you can still use MPFR_PRECISION from config.h if included before)
    mpfr_t T; mpfr_init2(T, MPFR_PRECISION);
    mpfr_set_d(T, 100.0, MPFR_RNDN);
    if (mpfr_cmp(t, T) < 0) {
        compute_theta_via_lggamma(theta, t);
    }
    else{
        compute_theta(theta, t);
    }
    
    mpfr_clear(T);
}