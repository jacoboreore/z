// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: z_theta.cpp
// Description:
//   θ(t) evaluation utilities for the Riemann–Siegel theta function.
//   Provides Stirling-based approximations with fixed K terms or with
//   an optimally truncated Bernoulli series. The public entry points are:
//     - compute_theta_stirling_K(theta, t, K)
//     - compute_theta_stirling(theta, t)           // uses CFG().theta_stirling_terms
//     - compute_theta_stirling_opt(theta, t)       // optimal truncation strategy
// -----------------------------------------------------------------------------

#include <mpfr.h>
#include "settings.h"

// Rational coefficients c_m for:
//   θ(t) = (t/2) log(t/(2π)) − t/2 − π/8 + Σ_{m=1..K} c_m / t^{2m−1}
// Values derive from Bernoulli numbers; table capped at m=9.
struct Rat {
    long long num;
    long long den;
};

static const Rat C[] = {
    {0, 1},                             // placeholder for 0-index
    {1, 48},                            // m = 1
    {7, 5760},                          // m = 2
    {31, 80640},                        // m = 3
    {127, 430080},                      // m = 4
    {511, 1216512},                     // m = 5
    {1414477, 319334400},               // m = 6
    {8191, 3670016},                    // m = 7
    {118518239, 135895449600LL},        // m = 8
    {5749691557LL, 1634992128000LL}     // m = 9
};

/**
 * @brief Stirling-series approximation of θ(t) with K correction terms.
 *
 * Implements
 *   θ(t) = (t/2) log(t/(2π)) − t/2 − π/8 + Σ_{m=1..K} c_m / t^{2m−1},
 * with coefficients c_m taken from the table above (valid up to m=9).
 *
 * @param[out] theta  Result θ(t).
 * @param[in]  t      Height parameter (t > 0).
 * @param[in]  K      Number of correction terms (clamped to [0, 9]).
 */
void compute_theta_stirling_K(mpfr_t theta, const mpfr_t t, int K)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    if (K < 0) K = 0;
    if (K > 9) K = 9; // extend C[] to raise this cap

    mpfr_t pi, two_pi, log_term, base, inv_t, pwr, term, tmp;
    mpfr_inits2(P, pi, two_pi, log_term, base, inv_t, pwr, term, tmp, (mpfr_ptr)0);

    // π and 2π
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    // base = t/(2π), log_term = (t/2) * log(t/(2π))
    mpfr_div(base, t, two_pi, MPFR_RNDN);
    mpfr_log(log_term, base, MPFR_RNDN);
    mpfr_mul(log_term, log_term, t, MPFR_RNDN);
    mpfr_div_ui(log_term, log_term, 2, MPFR_RNDN);

    // theta = log_term − t/2 − π/8
    mpfr_set(theta, log_term, MPFR_RNDN);
    mpfr_div_ui(tmp, t, 2, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);
    mpfr_div_ui(tmp, pi, 8, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);

    // series Σ c_m / t^{2m−1}
    mpfr_ui_div(inv_t, 1, t, MPFR_RNDN); // 1/t
    mpfr_set(pwr, inv_t, MPFR_RNDN);     // pwr = 1/t  for m=1

    for (int m = 1; m <= K; ++m)
    {
        // Ensure pwr = 1 / t^{2m−1}
        if (m > 1) {
            mpfr_mul(pwr, pwr, inv_t, MPFR_RNDN);
            mpfr_mul(pwr, pwr, inv_t, MPFR_RNDN);
        }
        mpfr_set_si(term, C[m].num, MPFR_RNDN);
        mpfr_div_si(term, term, C[m].den, MPFR_RNDN); // term = c_m
        mpfr_mul(term, term, pwr, MPFR_RNDN);         // c_m / t^{2m−1}
        mpfr_add(theta, theta, term, MPFR_RNDN);
    }

    mpfr_clears(pi, two_pi, log_term, base, inv_t, pwr, term, tmp, (mpfr_ptr)0);
}

/**
 * @brief Convenience wrapper for Stirling-series θ(t).
 *
 * Uses K = CFG().theta_stirling_terms, clamped to the table range.
 *
 * @param[out] theta  Result θ(t).
 * @param[in]  t      Height parameter (t > 0).
 */
void compute_theta_stirling(mpfr_t theta, const mpfr_t t)
{
    compute_theta_stirling_K(theta, t, CFG().theta_stirling_terms);
}

// -----------------------------------------------------------------------------
// Optional Lanczos-based θ(t) for small/medium t.
// Implementation retained for future use; requires <mpc.h>.
// -----------------------------------------------------------------------------
/*
// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π  (Lanczos, g=7, N=9)
// void compute_theta_lanczos(mpfr_t theta, const mpfr_t t) { ... }
*/

// Bernoulli numbers B_2, B_4, ..., B_20 (numerators/denominators).
// Extend if additional headroom is required.
static const long BNUM[] = { 1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611 };
static const long BDEN[] = { 6, 30, 42, 30, 66, 2730, 6, 510, 798, 330 };

/**
 * @brief θ(t) via Stirling with optimal truncation of the Bernoulli series.
 *
 * Computes
 *    θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
 * using a Stirling expansion evaluated at precision P+64, truncating when the
 * term magnitude falls below a cutoff ≈ 2^{-(P+10)}, then rounds back to P.
 *
 * Notes:
 *  - Works directly with z = a + i b, where a = 1/4 and b = t/2.
 *  - Uses the identity Im(z^{-(2k−1)}) = − r^{-(2k−1)} sin((2k−1) arg z).
 *
 * @param[out] theta Result θ(t) at CFG().mpfr_prec precision.
 * @param[in]  t     Height parameter (t > 0).
 */
void compute_theta_stirling_opt(mpfr_t theta, mpfr_srcptr t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    const mpfr_prec_t Q = P + 64;

    // z = a + i b with a = 1/4, b = t/2
    mpfr_t a, b, r, argz, logr;
    mpfr_inits2(Q, a, b, r, argz, logr, (mpfr_ptr)0);
    mpfr_set_d(a, 0.25, MPFR_RNDN);
    mpfr_div_ui(b, t, 2, MPFR_RNDN); // b = t/2

    // r = |z|, argz = atan2(b, a), log r
    mpfr_hypot(r, a, b, MPFR_RNDN);
    mpfr_atan2(argz, b, a, MPFR_RNDN);
    mpfr_log(logr, r, MPFR_RNDN);

    // s_im = (a−1/2) * argz + b*log r − b
    mpfr_t half, amhalf, sim, tmp;
    mpfr_inits2(Q, half, amhalf, sim, tmp, (mpfr_ptr)0);
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpfr_sub(amhalf, a, half, MPFR_RNDN);      // a−1/2 = −1/4
    mpfr_mul(sim, amhalf, argz, MPFR_RNDN);    // (a−1/2)*argz
    mpfr_mul(tmp, b, logr, MPFR_RNDN);         // b*log r
    mpfr_add(sim, sim, tmp, MPFR_RNDN);
    mpfr_sub(sim, sim, b, MPFR_RNDN);          // ... − b

    // Bernoulli series: Σ B_{2k}/[2k(2k−1)] * Im(z^{−(2k−1)}),
    // with Im(z^{−(2k−1)}) = − r^{−(2k−1)} sin((2k−1) arg z).
    mpfr_t inv, inv_r2, rpow, ang, s, coeff, mag, cutoff;
    mpfr_inits2(Q, inv, inv_r2, rpow, ang, s, coeff, mag, cutoff, (mpfr_ptr)0);

    mpfr_ui_div(inv, 1, r, MPFR_RNDN);     // 1/r
    mpfr_mul(inv_r2, inv, inv, MPFR_RNDN); // 1/r^2
    mpfr_set(rpow, inv, MPFR_RNDN);        // r^{−1}

    // cutoff ≈ 2^{−(P+10)} as a magnitude stop
    mpfr_set_ui(cutoff, 1, MPFR_RNDN);
    mpfr_div_2si(cutoff, cutoff, P + 10, MPFR_RNDN);

    const int KMAX = (int)(sizeof(BNUM) / sizeof(BNUM[0]));
    for (int k = 1; k <= KMAX; ++k)
    {
        // magnitude check: |B_{2k}|/(2k(2k−1)) * r^{−(2k−1)} < cutoff ?
        mpfr_set_si(mag, BNUM[k-1], MPFR_RNDN);
        mpfr_div_si(mag, mag, BDEN[k-1], MPFR_RNDN);
        mpfr_abs(mag, mag, MPFR_RNDN);
        mpfr_div_ui(mag, mag, (unsigned)(2 * k * (2 * k - 1)), MPFR_RNDN);
        mpfr_mul(mag, mag, rpow, MPFR_RNDN);
        if (mpfr_cmp(mag, cutoff) < 0)
            break;

        // ang = (2k−1) * arg z ; s = sin(ang)
        mpfr_mul_ui(ang, argz, 2 * k - 1, MPFR_RNDN);
        mpfr_sin(s, ang, MPFR_RNDN);

        // coeff = B_{2k}/[2k(2k−1)] * r^{−(2k−1)} * (−sin(ang))
        mpfr_set_si(coeff, BNUM[k-1], MPFR_RNDN);
        mpfr_div_si(coeff, coeff, BDEN[k-1], MPFR_RNDN);
        mpfr_div_ui(coeff, coeff, (unsigned)(2 * k * (2 * k - 1)), MPFR_RNDN);
        mpfr_mul(coeff, coeff, rpow, MPFR_RNDN);
        mpfr_neg(s, s, MPFR_RNDN);
        mpfr_mul(coeff, coeff, s, MPFR_RNDN);
        mpfr_add(sim, sim, coeff, MPFR_RNDN);

        // advance rpow *= r^{−2}
        mpfr_mul(rpow, rpow, inv_r2, MPFR_RNDN);
    }

    // subtract (t/2) log π
    mpfr_t pi_val, logpi, tp2;
    mpfr_inits2(Q, pi_val, logpi, tp2, (mpfr_ptr)0);
    mpfr_const_pi(pi_val, MPFR_RNDN);
    mpfr_log(logpi, pi_val, MPFR_RNDN);
    mpfr_mul(tp2, t, logpi, MPFR_RNDN);
    mpfr_div_ui(tp2, tp2, 2, MPFR_RNDN);
    mpfr_sub(sim, sim, tp2, MPFR_RNDN);

    // return at precision P
    mpfr_prec_round(theta, P, MPFR_RNDN);
    mpfr_set(theta, sim, MPFR_RNDN);

    mpfr_clears(a, b, r, argz, logr, half, amhalf, sim, tmp,
                inv, inv_r2, rpow, ang, s, coeff, mag, cutoff,
                pi_val, logpi, tp2, (mpfr_ptr)0);
}
