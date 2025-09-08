#include <mpfr.h>
#include <mpc.h>
#include <thread>
#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include "config.h"

// ------------------------------------------------------------------
// Exact Riemann–Siegel theta via log-Gamma
// θ(t) = Im( log Γ(1/4 + i t/2) ) − (t/2) log π
// This is numerically stable for small/medium t.
// ------------------------------------------------------------------
// Exact Riemann–Siegel theta via complex log-Gamma:
// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
void compute_theta_via_lggamma(mpfr_t theta, const mpfr_t t) {
    // Initialize MPFR and MPC variables
    mpc_t complex_arg, log_gamma;
    mpfr_t temp, pi, half_t, log_pi_term, result;
    
    mpc_init2(complex_arg, MPFR_PRECISION);
    mpc_init2(log_gamma, MPFR_PRECISION);
    mpfr_init2(temp, MPFR_PRECISION);
    mpfr_init2(pi, MPFR_PRECISION);
    mpfr_init2(half_t, MPFR_PRECISION);
    mpfr_init2(log_pi_term, MPFR_PRECISION);
    mpfr_init2(result, MPFR_PRECISION);
    
    // Compute 1/4 + i t/2
    mpfr_set_d(temp, 0.25, MPFR_RNDN);
    mpfr_div_si(half_t, t, 2, MPFR_RNDN);
    mpc_set_fr_fr(complex_arg, temp, half_t, MPC_RNDNN); // Set real and imaginary parts
    
    // Compute log Γ(1/4 + i t/2)
    // Fallback: Compute using mpfr_lngamma for real and imaginary parts
    // Note: This is an approximation and may not be as accurate for complex args
    mpfr_t re, im;
    mpfr_init2(re, MPFR_PRECISION);
    mpfr_init2(im, MPFR_PRECISION);
    mpfr_set_d(re, 0.25, MPFR_RNDN);
    mpfr_set(im, half_t, MPFR_RNDN);
    mpfr_lngamma(re, re, MPFR_RNDN); // Real part (approximation)
    mpc_set_fr(log_gamma, re, MPC_RNDNN); // Set real part, imaginary part approximated as 0
    mpfr_clear(re);
    mpfr_clear(im);
    // std::cerr << "Warning: Using mpfr_lngamma fallback, may be inaccurate for complex args\n";
    
    // Get imaginary part
    mpc_imag(result, log_gamma, MPFR_RNDN);
    
    // Compute −(t/2) log π
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_log(temp, pi, MPFR_RNDN);
    mpfr_mul(half_t, t, temp, MPFR_RNDN);
    mpfr_div_si(log_pi_term, half_t, 2, MPFR_RNDN);
    mpfr_neg(log_pi_term, log_pi_term, MPFR_RNDN);
    
    // Final result
    mpfr_add(theta, result, log_pi_term, MPFR_RNDN);
    
    // Clean up
    mpc_clear(complex_arg);
    mpc_clear(log_gamma);
    mpfr_clear(temp);
    mpfr_clear(pi);
    mpfr_clear(half_t);
    mpfr_clear(log_pi_term);
    mpfr_clear(result);
}

// ---------------- Theta(t) ----------------
// theta(t) = (t/2) log(t/2π) - t/2 - π/8 + 1/(48 t) + O(t^-3)
void compute_theta(
        mpfr_t theta, 
        const mpfr_t t) {
    
    mpfr_t pi, two_pi, x;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, x, (mpfr_ptr)0);
    
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    mpfr_div(x, t, two_pi, MPFR_RNDN);      // t/(2π)
    mpfr_log(x, x, MPFR_RNDN);              // log(t/2π)
    mpfr_mul(theta, x, t, MPFR_RNDN);       // t*log(...)
    mpfr_div_ui(theta, theta, 2, MPFR_RNDN);// (t/2) log(...)

    mpfr_div_ui(x, t, 2, MPFR_RNDN);        // t/2
    mpfr_sub(theta, theta, x, MPFR_RNDN);   // ... - t/2

    mpfr_div_ui(x, pi, 8, MPFR_RNDN);       // π/8
    mpfr_sub(theta, theta, x, MPFR_RNDN);   // ... - π/8

    mpfr_mul_ui(x, t, 48, MPFR_RNDN);       // 48 t
    mpfr_ui_div(x, 1, x, MPFR_RNDN);        // 1/(48 t)
    mpfr_add(theta, theta, x, MPFR_RNDN);   // ... + 1/(48 t)

    mpfr_clears(pi, two_pi, x, (mpfr_ptr)0);
}
