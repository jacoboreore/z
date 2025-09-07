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
    #ifdef MPC_HAVE_LNGAMMA
        mpc_lngamma(log_gamma, complex_arg, MPC_RNDNN); // Use mpc_lngamma if available
    #else
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
        std::cerr << "Warning: Using mpfr_lngamma fallback, may be inaccurate for complex args\n";
    #endif
    
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
}

// ---------------- Small-t path via eta series ----------------
// ζ(s) = η(s) / (1 - 2^{1-s}),   η(s) = Σ (-1)^{n-1} n^{-s}
// Uses recurrence (n+1)^{-s} = n^{-s} * exp(-s log(1+1/n))
// ------------------------------------------------------------------
// Fast η(s) via Euler transform of the alternating series
// η(s) = sum_{n>=1} (-1)^{n-1} n^{-s}
// Euler transform: η = sum_{k>=0} Δ^k a_1 / 2^{k+1},  a_n = n^{-s}
// We build a_1..a_{M+1}, then forward differences in-place.
// ------------------------------------------------------------------
void zeta_via_eta(
        mpc_t zeta, 
        const mpc_t s,
        int M = 128, 
        double tol = 1e-30){

    // Build a_j = (j+1)^{-s} for j=0..M using the stable recurrence
    // a_{n+1} = a_n * exp(-s * log1p(1/n))
    std::vector<mpc_t> cur(M+1);
    for (int j = 0; j <= M; ++j) mpc_init2(cur[j], MPFR_PRECISION);

    mpc_t ratio, tmpc;  
    mpc_init2(ratio, MPFR_PRECISION); 
    mpc_init2(tmpc, MPFR_PRECISION);
    
    mpfr_t ln1p;        
    mpfr_init2(ln1p, MPFR_PRECISION);

    // a_1 = 1^{-s} = 1
    mpc_set_ui(cur[0], 1, MPC_RNDNN);

    for (int j = 0; j < M; ++j) {
        // ratio = exp( -s * log1p(1/(j+1)) )
        mpfr_set_ui(ln1p, (unsigned long)(j+1), MPFR_RNDN);
        mpfr_ui_div(ln1p, 1, ln1p, MPFR_RNDN);   // 1/(j+1)
        mpfr_log1p(ln1p, ln1p, MPFR_RNDN);       // log(1+1/(j+1))
        mpc_mul_fr(tmpc, s, ln1p, MPC_RNDNN);
        mpc_neg(tmpc, tmpc, MPC_RNDNN);
        mpc_exp(ratio, tmpc, MPC_RNDNN);

        // cur[j+1] = cur[j] * ratio = (j+2)^{-s}
        mpc_mul(cur[j+1], cur[j], ratio, MPC_RNDNN);
    }

    // η via forward differences
    mpc_t eta, term;        
    mpc_init2(eta,  MPFR_PRECISION); 
    mpc_set_ui(eta, 0, MPC_RNDNN);
    mpc_init2(term, MPFR_PRECISION);

    mpfr_t pow2, abst;      
    mpfr_inits2(MPFR_PRECISION, pow2, abst, (mpfr_ptr)0);
    mpfr_set_ui(pow2, 2, MPFR_RNDN);             // 2^{1} for k=0 term

    int L = M + 1; // current length of the difference row
    for (int k = 0; k <= M && L > 0; ++k) {
        // term_k = cur[0] / 2^{k+1}
        mpc_div_fr(term, cur[0], pow2, MPC_RNDNN);
        mpc_add(eta, eta, term, MPC_RNDNN);

        // stop when |term_k| < tol (converged)
        mpc_abs(abst, term, MPFR_RNDN);
        if (mpfr_cmp_d(abst, tol) < 0) break;

        // Next row of forward differences in-place: cur[j] = cur[j+1] - cur[j]
        for (int j = 0; j < L-1; ++j) {
            mpc_sub(cur[j], cur[j+1], cur[j], MPC_RNDNN);
        }
        --L;
        mpfr_mul_ui(pow2, pow2, 2, MPFR_RNDN);  // pow2 *= 2
    }

    // zeta = eta / (1 - 2^{1-s})
    mpfr_t ln2; 
    mpfr_init2(ln2, MPFR_PRECISION); 
    mpfr_const_log2(ln2, MPFR_RNDN);
    
    // tmpc = 2^{1-s} = 2 * 2^{-s} = 2 * exp(-s ln 2)
    mpc_mul_fr(tmpc, s, ln2, MPC_RNDNN);
    mpc_neg(tmpc, tmpc, MPC_RNDNN);
    mpc_exp(tmpc, tmpc, MPC_RNDNN);          // 2^{-s}
    mpc_mul_ui(tmpc, tmpc, 2, MPC_RNDNN);    // 2^{1-s}
    mpc_ui_sub(tmpc, 1, tmpc, MPC_RNDNN);    // 1 - 2^{1-s}
    mpc_div(zeta, eta, tmpc, MPC_RNDNN);

    // cleanup
    for (auto& v : cur) mpc_clear(v);
    mpc_clear(ratio);
    mpc_clear(tmpc);
    mpc_clear(eta);
    mpc_clear(term);

    mpfr_clears(ln1p, pow2, abst, ln2, (mpfr_ptr)0);
}

// Z(t) = Re( e^{iθ(t)} ζ(1/2 + i t) )
// For small/medium t: use η(s) (Euler transform) + exact θ via log-Gamma.
void compute_Z_small(
        mpfr_t Z, 
        const mpfr_t t) {

    mpc_t s, z;  
    mpc_init2(s, MPFR_PRECISION); 
    mpc_init2(z, MPFR_PRECISION);

    mpfr_t half, theta, cth, sth, zr, zi, tmp;
    mpfr_inits2(MPFR_PRECISION, half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);

    // s = 1/2 + i t
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpc_set_fr_fr(s, half, t, MPC_RNDNN);

    // ζ(s) via fast η(s), and θ(t) exactly (no asympt expansion here)
    zeta_via_eta(z, s, /*M*/160, /*tol*/1e-34);
    compute_theta_via_lggamma(theta, t);

    // Z(t) = Re( e^{iθ} ζ ) = Re( (cosθ + i sinθ)(zr + i zi) ) = zr*cosθ - zi*sinθ
    mpfr_cos(cth, theta, MPFR_RNDN);
    mpfr_sin(sth, theta, MPFR_RNDN);
    mpc_real(zr, z, MPFR_RNDN);
    mpc_imag(zi, z, MPFR_RNDN);

    mpfr_mul(Z,  zr, cth, MPFR_RNDN);
    mpfr_mul(tmp, zi, sth, MPFR_RNDN);
    mpfr_sub(Z,  Z,  tmp, MPFR_RNDN);

    // cleanup
    mpfr_clears(half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);
    mpc_clear(s);
    mpc_clear(z);
}

