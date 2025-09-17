// z_calculator_thread.cpp  (single-thread RS with Euler–Maclaurin tail)
#include <mpfr.h>
#include <cmath>
#include <iostream>
#include "z_theta.h"   // compute_theta_stirling_opt(...)
#include "logs.h"      // //log_mpfr(...)
#include "settings.h"  // CFG()

// // Simple diagnostic tail bound you already had
// void remainder_bound(mpfr_t R, const mpfr_t t)
// {
//     const double C = CFG().rs_tail_C;
//     mpfr_t r, tmp;
//     mpfr_inits2(CFG().mpfr_prec, r, tmp, (mpfr_ptr)0);
//     mpfr_sqrt(tmp, t, MPFR_RNDN); // √t
//     mpfr_sqrt(r, tmp, MPFR_RNDN); // t^(1/4)
//     mpfr_ui_div(R, 1, r, MPFR_RNDN);
//     mpfr_mul_d(R, R, C, MPFR_RNDN);
//     mpfr_clears(r, tmp, (mpfr_ptr)0);
// }

// Slow but simple direct sum — useful to sanity check S for small N.
static void rs_sum_direct(mpfr_t S, mpfr_srcptr t, mpfr_srcptr theta, long N)
{
    mpfr_set_zero(S, 0);
    mpfr_t n, lnn, phi, cphi, sqrt_n, term;
    mpfr_inits2(CFG().mpfr_prec, n, lnn, phi, cphi, sqrt_n, term, (mpfr_ptr)0);

    for (long i = 1; i <= N; ++i) {
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_log(lnn, n, MPFR_RNDN);
        mpfr_mul(phi, t, lnn, MPFR_RNDN);
        mpfr_sub(phi, theta, phi, MPFR_RNDN);
        mpfr_cos(cphi, phi, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
        mpfr_div(term, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(S, S, term, MPFR_RNDN);
    }

    mpfr_clears(n, lnn, phi, cphi, sqrt_n, term, (mpfr_ptr)0);
}

// Z(t) with RS truncation + Euler–Maclaurin tail (integral + -1/12 f′)
void compute_Z_threaded(mpfr_t Z, mpfr_t t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;

    // N = floor( sqrt(t/(2π)) ), at least 1
    mpfr_t pi, two_pi, tmp, theta;
    mpfr_inits2(P, pi, two_pi, tmp, theta, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    if (N < 1) N = 1;

    // θ(t): optimized Stirling variant
    compute_theta_stirling_opt(theta, t);

    // Logging header
    // if (CFG().rs_log_debug) {
    //     //log_mpfr("[RS] t", t, CFG().rs_log_digits);
    //     //log_mpfr("[RS] theta", theta, CFG().rs_log_digits);
    //     std::cerr << "[RS] N=" << N << "\n";
    // }

    // ----- Main RS sum via stable rotation recurrence -----
    mpfr_t S;                      mpfr_init2(S, P);               mpfr_set_zero(S, 0);
    mpfr_t cphi, sphi;             mpfr_inits2(P, cphi, sphi, (mpfr_ptr)0);
    mpfr_t dphi, cd, sd;           mpfr_inits2(P, dphi, cd, sd, (mpfr_ptr)0);
    mpfr_t new_c, new_s, tmp1;     mpfr_inits2(P, new_c, new_s, tmp1, (mpfr_ptr)0);
    mpfr_t i_mp, sqrt_i;           mpfr_inits2(P, i_mp, sqrt_i, (mpfr_ptr)0);

    // seed at n = 1 → φ = θ
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sphi, cphi, theta, MPFR_RNDN);
#else
    mpfr_cos(cphi, theta, MPFR_RNDN);
    mpfr_sin(sphi, theta, MPFR_RNDN);
#endif
    mpfr_add(S, S, cphi, MPFR_RNDN);

    // advance i = 1..N-1 using Δφ = t*log1p(1/i)
    for (long i = 1; i <= N - 1; ++i) {
        mpfr_set_ui(i_mp, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(dphi, 1, i_mp, MPFR_RNDN);
        mpfr_log1p(dphi, dphi, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate by -Δφ
        // new_c = c*cd + s*sd
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp1,  sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp1, MPFR_RNDN);
        // new_s = s*cd − c*sd
        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp1,  cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp1, MPFR_RNDN);
        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        // add term at n = i+1
        mpfr_set_ui(i_mp, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_i, i_mp, MPFR_RNDN);
        mpfr_div(tmp1, cphi, sqrt_i, MPFR_RNDN);
        mpfr_add(S, S, tmp1, MPFR_RNDN);
    }

    // Optionally verify rotation sum matches a direct sum (cheap at small N)
    if (CFG().rs_log_debug) {
        mpfr_t S_dir; mpfr_init2(S_dir, P);
        rs_sum_direct(S_dir, t, theta, N);
        // //log_mpfr("[RS] S_rec", S,     CFG().rs_log_digits);
        // //log_mpfr("[RS] S_dir", S_dir, CFG().rs_log_digits);
        mpfr_clear(S_dir);
    }

    // ----- Assemble Z with Euler–Maclaurin tail -----
    // Start with Z_main = 2*S
    mpfr_mul_ui(Z, S, 2, MPFR_RNDN);

    // Build φ_N, cos φ_N, sin φ_N, √N
    mpfr_t nN, lnN, phiN, cN, sN, sqrtN;
    mpfr_inits2(P, nN, lnN, phiN, cN, sN, sqrtN, (mpfr_ptr)0);
    mpfr_set_ui(nN, (unsigned long)N, MPFR_RNDN);
    mpfr_log(lnN, nN, MPFR_RNDN);
    mpfr_mul(phiN, t, lnN, MPFR_RNDN);
    mpfr_sub(phiN, theta, phiN, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sN, cN, phiN, MPFR_RNDN); // sN=sin φN, cN=cos φN
#else
    mpfr_cos(cN, phiN, MPFR_RNDN);
    mpfr_sin(sN, phiN, MPFR_RNDN);
#endif
    mpfr_sqrt(sqrtN, nN, MPFR_RNDN);

    // Half-weight last term: subtract f(N) once (because Z_main was 2*S)
    // f(N) = cos φ_N / √N
    mpfr_t termN; mpfr_init2(termN, P);
    mpfr_div(termN, cN, sqrtN, MPFR_RNDN);
    mpfr_sub(Z, Z, termN, MPFR_RNDN);

    // Tail integral I = ( (1/2)C + t S ) / ( √N * (t^2 + 1/4) )
    mpfr_t t2, denom, num_I, I;
    mpfr_inits2(P, t2, denom, num_I, I, (mpfr_ptr)0);
    mpfr_mul(t2, t, t, MPFR_RNDN);
    mpfr_set_d(denom, 0.25, MPFR_RNDN);
    mpfr_add(denom, denom, t2, MPFR_RNDN);                 // denom = t^2 + 1/4
    mpfr_mul_d(num_I, cN, 0.5, MPFR_RNDN);                 // (1/2)C
    mpfr_mul(tmp1, sN, t, MPFR_RNDN);                      // t S
    mpfr_add(num_I, num_I, tmp1, MPFR_RNDN);               // (1/2)C + tS
    mpfr_div(I, num_I, sqrtN, MPFR_RNDN);
    mpfr_div(I, I, denom, MPFR_RNDN);

    // f'(N) = N^{-3/2} ( -1/2 C + t S )
    mpfr_t fprime, N32;
    mpfr_inits2(P, fprime, N32, (mpfr_ptr)0);
    mpfr_mul(N32, nN, sqrtN, MPFR_RNDN);                   // N^{3/2} = N * √N
    mpfr_mul_d(fprime, cN, -0.5, MPFR_RNDN);               // -1/2 C
    mpfr_add(fprime, fprime, tmp1, MPFR_RNDN);             // -1/2 C + t S
    mpfr_div(fprime, fprime, N32, MPFR_RNDN);              // divide by N^{3/2}

    // ΔZ_tail = 2*I  - (1/6)*f'(N)
    mpfr_t dZ1, dZ2;
    mpfr_inits2(P, dZ1, dZ2, (mpfr_ptr)0);
    mpfr_mul_ui(dZ1, I, 2, MPFR_RNDN);                     // 2*I
    mpfr_mul_d(dZ2, fprime, -1.0/6.0, MPFR_RNDN);          // -(1/6) f'(N)
    mpfr_add(Z, Z, dZ1, MPFR_RNDN);
    mpfr_add(Z, Z, dZ2, MPFR_RNDN);

    // Diagnostics
    if (CFG().rs_log_debug) {
        mpfr_t R; mpfr_init2(R, P);
        remainder_bound(R, t);
        // //log_mpfr("[RS] |R_bound|", R, CFG().rs_log_digits);
        // //log_mpfr("[RS] Z", Z, CFG().rs_log_digits);
        mpfr_clear(R);
    }

    // Cleanup
    mpfr_clears(pi, two_pi, tmp, theta, S, cphi, sphi, dphi, cd, sd,
                new_c, new_s, tmp1, i_mp, sqrt_i,
                nN, lnN, phiN, cN, sN, sqrtN, termN,
                t2, denom, num_I, I, fprime, N32, dZ1, dZ2, (mpfr_ptr)0);
}
