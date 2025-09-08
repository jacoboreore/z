#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include <iostream>
#include "z_calculator_thread.h"
#include "z_calculator_small.h"
#include "config.h"
#include "z_theta.h"
#include "logs.h"

// ---------------- Unified compute_Z ----------------
void compute_Z(mpfr_t Z, mpfr_t t) {
    
    // choose N
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    
    if (N < 1) N = 1;
    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);

    if (N < 8) {
        compute_Z_small(Z, t);
    } else {
        printf("USING THREADING...\n");
        compute_Z_threaded(Z, t);
    }
}

// ---------------- Prefilter ----------------
// Compute partial RS sum up to K=⌊αN⌋, bound the tail + remainder.
// Discard t if |2 S_K| - 2*T_K - R > τ (i.e., cannot be near zero).
bool z_prefilter(const mpfr_t t, double alpha, double tau) {
    
    if (alpha <= 0.0) alpha = 0.25;
    if (alpha >= 1.0) alpha = 0.8;

    mpfr_t pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(Nmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(Nmp, Nmp, MPFR_RNDN);
    long N = mpfr_get_si(Nmp, MPFR_RNDN);
    
    if (N < 1) N = 1;

    if (N < 8) {
        // Option A: always keep near tiny N
        // return true;

        // Option B (better): cheap small-Z gate
        mpfr_t Zs; mpfr_init2(Zs, MPFR_PRECISION);
        compute_Z_small(Zs, t);
        double az = std::fabs(mpfr_get_d(Zs, MPFR_RNDN));
        mpfr_clear(Zs);
        if (az <= 1e-2) return true;   // near zero → keep
        if (az >= 0.3)  return false;  // clearly far → drop
        return true;                   // ambiguous? favor recall
    }

    long K = (long)std::floor(alpha * N);
    if (K < 4) K = (N >= 4) ? 4 : N;
    if (K > N) K = N;

    compute_theta_big(theta, t);

    // Start at n=1
    mpfr_set_zero(S, 0);
    mpfr_set_ui(n, 1, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, t, ln_n, MPFR_RNDN);
    mpfr_sub(phi, theta, phi, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif
    mpfr_add(S, S, cphi, MPFR_RNDN); // term at n=1

    mpfr_t new_c, new_s, sqrt_i1, i1, tmp;
    mpfr_inits2(MPFR_PRECISION, new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

    for (long i = 1; i < K; ++i) {
        // Δφ = t * log1p(1/i)
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(q, 1, n, MPFR_RNDN);          // 1/i
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate (φ -> φ - Δφ)
        
        // new_c = c*cd + s*sd
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        // new_s = s*cd - c*sd
        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        mpfr_set_ui(i1, (unsigned long)(i+1), MPFR_RNDN);
        mpfr_sqrt(sqrt_i1, i1, MPFR_RNDN);
        mpfr_div(new_c, cphi, sqrt_i1, MPFR_RNDN);
        mpfr_add(S, S, new_c, MPFR_RNDN);
        
    }

    mpfr_clears(new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

    mpfr_t twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau;
    mpfr_inits2(MPFR_PRECISION, twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau, (mpfr_ptr)0);

    mpfr_mul_ui(twoS, S, 2, MPFR_RNDN);
    mpfr_abs(abs_twoS, twoS, MPFR_RNDN);

    mpfr_set_ui(Nmp, (unsigned long)N, MPFR_RNDN);
    mpfr_set_ui(Kmp, (unsigned long)K, MPFR_RNDN);
    mpfr_sqrt(sqrtN, Nmp, MPFR_RNDN);
    mpfr_sqrt(sqrtK, Kmp, MPFR_RNDN);
    mpfr_sub(tail, sqrtN, sqrtK, MPFR_RNDN);
    mpfr_mul_ui(tail, tail, 2, MPFR_RNDN); // T_K
    mpfr_mul_ui(tail, tail, 2, MPFR_RNDN); // 2*T_K for Z

    remainder_bound(R, t);

    // L = |2S_K| - 2*T_K - R
    // mpfr_sub(L, abs_twoS, tail, MPFR_RNDN);
    // mpfr_sub(L, L, R, MPFR_RNDN);

    // mpfr_set_d(Tau, tau, MPFR_RNDN);
    // bool keep = (mpfr_cmp(L, Tau) <= 0);
    double tau_boost = 0.50; // tune: 0.4–0.7 works well
    double tau_eff = tau + tau_boost / std::sqrt((double)N);

    mpfr_t TauEff; mpfr_init2(TauEff, MPFR_PRECISION);
    mpfr_set_d(TauEff, tau_eff, MPFR_RNDN);
    bool keep = (mpfr_cmp(L, TauEff) <= 0);
    mpfr_clear(TauEff);

    std::cerr << "[prefilter] N=" << N << " K=" << K << " alpha=" << alpha << " tau=" << tau << "\n";
    log_mpfr("[prefilter] |2S_K|", abs_twoS, 18);
    log_mpfr("[prefilter] 2*T_K", tail, 18);
    log_mpfr("[prefilter] R", R, 18);
    log_mpfr("[prefilter] L", L, 18);
    std::cerr << "[prefilter] keep=" << (keep ? "YES" : "NO") << "\n";

    mpfr_clears(pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd,
                twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau, (mpfr_ptr)0);

    return keep;
}