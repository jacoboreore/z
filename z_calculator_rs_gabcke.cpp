// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: z_calculator_rs_gabcke.cpp
// Description:
//   - High-precision Hardy Z(t) via the Riemann–Siegel main sum plus
//     Gabcke-style remainder (first correction term R₁ with optional C₁/C₂).
//   - Fast θ(t) using log Γ via Spouge with reflection (stable for Im z large).
//   - Parallel rotation kernel for RS main sum (OpenMP), with dynamic blocking.
//
// Public entry points:
//   - void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t)
//
// File-local helpers (selected):
//   - static long rs_N_from_t(...)
//   - static void Z_rs_main(...), Z_rs_main_parallel(...)
//   - static void compute_theta_spouge(...)
//   - static void R_gabcke(...)
//   - static void mpc_loggamma_spouge_shift(...), _reflect(...)
// -----------------------------------------------------------------------------

#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <omp.h>
#include <sstream>

#include "settings.h"
#include "logs.h"
#include "utils.h"
#include "z_calculator_rs_gabcke.h"

/**
 * @brief Convert MPFR precision in bits to approximate decimal digits.
 * @param P MPFR precision in bits.
 * @return Approximate decimal digits.
 */
static inline long dec_digits_from_bits(mpfr_prec_t P)
{
    return (long)std::ceil(P * 0.3010299956639811952137389L);
}

// -----------------------------------------------------------------------------
// RS parameterization
// -----------------------------------------------------------------------------

/**
 * @brief Compute RS parameters from t:
 *        a = sqrt(t/(2π)), N = floor(a), p = frac(a).
 *
 * @param t      Height parameter.
 * @param P      MPFR precision in bits.
 * @param a_out  Output a.
 * @param p_out  Output p.
 * @param C      Cached constants.
 * @return N = floor(a), clamped to at least 1.
 */
static long rs_N_from_t(const mpfr_t t, mpfr_prec_t P, mpfr_t a_out, mpfr_t p_out, const Constants &C)
{
    mpfr_t tmp;
    mpfr_init2(tmp, P);
    mpfr_div(tmp, t, C.two_pi, MPFR_RNDN);
    mpfr_sqrt(a_out, tmp, MPFR_RNDN); // a = sqrt(t/(2π))
    long N = mpfr_get_si(a_out, MPFR_RNDN);
    if (N < 1)
        N = 1;
    mpfr_set_si(tmp, N, MPFR_RNDN);
    mpfr_sub(p_out, a_out, tmp, MPFR_RNDN); // p = a - floor(a)
    mpfr_clear(tmp);
    return N;
}

// -----------------------------------------------------------------------------
// RS main sum (threaded)
// -----------------------------------------------------------------------------

/**
 * @brief Parallel Riemann–Siegel main sum using rotation updates.
 *
 * Rotates cos/sin by Δφ = t*log1p(1/i) to avoid log(n) recomputation.
 * Uses per-thread partial accumulators and static scheduling for determinism.
 *
 * @param Z           Output main-sum contribution.
 * @param t           Height parameter.
 * @param theta       θ(t).
 * @param P           MPFR precision in bits.
 * @param C           Cached constants.
 * @param max_threads Optional thread cap (≤0 → omp_get_max_threads()).
 * @param block_size  Optional block size (≤0 → dynamic ~16*sqrt(N), clamped).
 */
static void Z_rs_main_parallel(
    mpfr_t Z,
    const mpfr_t t,
    const mpfr_t theta,
    mpfr_prec_t P,
    const Constants &C,
    int max_threads = 0,
    long block_size = 0)
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    if (max_threads <= 0)
        max_threads = omp_get_max_threads();

    // Dynamic block: ~16*sqrt(N), clamped.
    if (block_size <= 0)
    {
        long dyn = static_cast<long>(std::sqrt(static_cast<long double>(N)) * 16.0L);
        if (dyn < 4096)
            dyn = 4096;
        if (dyn > 32768)
            dyn = 32768;
        block_size = dyn;
    }
    const long num_blocks = (N + block_size - 1) / block_size;

    mpfr_set_zero(Z, 0);

    std::vector<mpfr_t> partial(max_threads);
    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_init2(partial[i], P);
        mpfr_set_zero(partial[i], 0);
    }

#pragma omp parallel for schedule(static)
    for (long b = 0; b < num_blocks; ++b)
    {
        const int tid = omp_get_thread_num();
        mpfr_t Sloc;
        mpfr_init2(Sloc, P);
        mpfr_set_zero(Sloc, 0);

        long n0 = b * block_size + 1;
        long n1 = std::min(n0 + block_size - 1, N);

        // Local rotation state for this block
        mpfr_t n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp;
        mpfr_inits2(P, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);

        mpfr_set_ui(n_mp, (unsigned long)n0, MPFR_RNDN);
        mpfr_log(phi0, n_mp, MPFR_RNDN);
        mpfr_mul(phi0, phi0, t, MPFR_RNDN);
        mpfr_neg(phi0, phi0, MPFR_RNDN);
        mpfr_add(phi0, phi0, theta, MPFR_RNDN);

        mpfr_sin_cos(sphi, cphi, phi0, MPFR_RNDN);

        // n = n0
        mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);

        for (long i = n0; i < n1; ++i)
        {
            // rotate by Δφ = t*log1p(1/i)
            mpfr_set_ui(n_mp, (unsigned long)i, MPFR_RNDN);
            mpfr_ui_div(dphi, 1, n_mp, MPFR_RNDN);
            mpfr_log1p(dphi, dphi, MPFR_RNDN);
            mpfr_mul(dphi, dphi, t, MPFR_RNDN);
            mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

            mpfr_t new_c, new_s;
            mpfr_inits2(P, new_c, new_s, (mpfr_ptr)0);
            mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
            mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

            mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
            mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

            mpfr_set(cphi, new_c, MPFR_RNDN);
            mpfr_set(sphi, new_s, MPFR_RNDN);
            mpfr_clears(new_c, new_s, (mpfr_ptr)0);

            // contribute n = i + 1
            mpfr_set_ui(n_mp, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
            mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
            mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);
        }

        mpfr_add(partial[tid], partial[tid], Sloc, MPFR_RNDN);
        mpfr_clears(Sloc, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);
    }

    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_add(Z, Z, partial[i], MPFR_RNDN);
        mpfr_clear(partial[i]);
    }

    // factor 2 · Re part already baked as cos(...) / sqrt(n) sum
    mpfr_mul_ui(Z, Z, 2, MPFR_RNDN);
    mpfr_clears(a, p, (mpfr_ptr)0);
}

// -----------------------------------------------------------------------------
// Spouge parameter heuristic
// -----------------------------------------------------------------------------

/**
 * @brief Choose Spouge parameter a from precision, conservatively.
 * @param P MPFR precision in bits.
 * @return Spouge parameter a in [30,80] for typical high-precision ranges.
 */
int spouge_a_for_prec(mpfr_prec_t P)
{
    int digs = (int)dec_digits_from_bits(P);
    int a = (int)std::lround(0.45 * digs);
    if (a < 30)
        a = 30;
    if (a > 80)
        a = 80;
    return a;
}

// -----------------------------------------------------------------------------
// log Γ via Spouge (shifted) and reflection
// -----------------------------------------------------------------------------

/**
 * @brief log Γ via Spouge on the right half-plane after shifting by (a−1).
 *
 * @param out Output complex log Γ(z).
 * @param z   Input point.
 * @param S   Spouge coefficient table.
 * @param P   MPFR precision in bits.
 */
static void mpc_loggamma_spouge_shift(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P)
{
    mpc_t w, s, term;
    mpc_init2(w, P);
    mpc_init2(s, P);
    mpc_init2(term, P);
    mpc_set(w, z, MPC_RNDNN);
    mpfr_add_si(mpc_realref(w), mpc_realref(w), S.a - 1, MPFR_RNDN);

    mpc_set_fr(s, &S.c[0], MPC_RNDNN);
    for (int k = 1; k < S.a; ++k)
    {
        mpc_set(term, z, MPC_RNDNN);
        mpfr_add_si(mpc_realref(term), mpc_realref(term), k - 1, MPFR_RNDN);
        mpc_ui_div(term, 1, term, MPC_RNDNN);
        mpc_mul_fr(term, term, &S.c[k], MPC_RNDNN);
        mpc_add(s, s, term, MPC_RNDNN);
    }

    mpc_log(out, w, MPC_RNDNN);

    mpc_t zmh;
    mpc_init2(zmh, P);
    mpc_set(zmh, z, MPC_RNDNN);
    mpfr_sub_d(mpc_realref(zmh), mpc_realref(zmh), 0.5, MPFR_RNDN);

    mpc_mul(out, zmh, out, MPC_RNDNN);
    mpc_sub(out, out, w, MPC_RNDNN);

    mpc_log(term, s, MPC_RNDNN);
    mpc_add(out, out, term, MPC_RNDNN);

    mpc_clear(zmh);
    mpc_clear(w);
    mpc_clear(s);
    mpc_clear(term);
}

/**
 * @brief Reflection identity for log Γ:
 *        log Γ(z) = log π − log(sin(π z)) − log Γ(1−z).
 *
 * @param out Output complex log Γ(z).
 * @param z   Input point.
 * @param S   Spouge coefficient table.
 * @param P   MPFR precision in bits.
 * @param C   Cached constants (π, log π).
 */
static void mpc_loggamma_spouge_reflect(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P, const Constants &C)
{
    mpc_t one_minus_z, lg1mz, sin_pi_z, t;
    mpc_init2(one_minus_z, P);
    mpc_init2(lg1mz, P);
    mpc_init2(sin_pi_z, P);
    mpc_init2(t, P);

    mpc_set_si(one_minus_z, 1, MPC_RNDNN);
    mpc_sub(one_minus_z, one_minus_z, z, MPC_RNDNN);

    mpc_loggamma_spouge_shift(lg1mz, one_minus_z, S, P);

    // sin(π z)
    mpc_set_fr(t, C.pi, MPC_RNDNN);
    mpc_mul(t, t, z, MPC_RNDNN);
    mpc_sin(sin_pi_z, t, MPC_RNDNN);

    // out = log π − log(sin π z) − log Γ(1−z)
    mpc_log(out, sin_pi_z, MPC_RNDNN);
    mpc_neg(out, out, MPC_RNDNN);
    mpc_add_fr(out, out, C.log_pi, MPC_RNDNN);
    mpc_sub(out, out, lg1mz, MPC_RNDNN);

    mpc_clear(one_minus_z);
    mpc_clear(lg1mz);
    mpc_clear(sin_pi_z);
    mpc_clear(t);
}

// -----------------------------------------------------------------------------
// θ(t)
// -----------------------------------------------------------------------------

/**
 * @brief Compute θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π.
 *
 * Uses Spouge + reflection for stability at large t.
 *
 * @param theta Output θ(t).
 * @param t     Height parameter.
 * @param P     MPFR precision in bits.
 * @param S     Spouge coefficient table.
 * @param C     Cached constants.
 */
void compute_theta_spouge(mpfr_t theta, const mpfr_t t, mpfr_prec_t P, const Spouge &S, const Constants &C)
{
    mpc_t z, lg;
    mpc_init2(z, P);
    mpc_init2(lg, P);

    mpfr_t tr;
    mpfr_init2(tr, P);
    mpfr_div_ui(tr, t, 2, MPFR_RNDN);
    mpc_set_d_d(z, 0.25, 0.0, MPC_RNDNN);
    mpfr_set(mpc_imagref(z), tr, MPFR_RNDN);

    mpc_loggamma_spouge_reflect(lg, z, S, P, C);

    mpfr_set(theta, mpc_imagref(lg), MPFR_RNDN);
    mpfr_t tmp;
    mpfr_init2(tmp, P);
    mpfr_mul(tmp, C.log_pi, tr, MPFR_RNDN); // (t/2) log π
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);

    mpfr_clear(tr);
    mpfr_clear(tmp);
    mpc_clear(z);
    mpc_clear(lg);
}

// -----------------------------------------------------------------------------
// ψ(p) and derivatives
// -----------------------------------------------------------------------------

/**
 * @brief ψ(p) = cos(2π(p² − p − 1/16)) / cos(2π p), with gentle handling near poles.
 *
 * @param out Output ψ(p).
 * @param p   Fractional part of a.
 * @param P   MPFR precision in bits.
 * @param C   Cached constants.
 */
static void psi_p(mpfr_t out, const mpfr_t p, mpfr_prec_t P, const Constants &C)
{
    mpfr_t num_arg, den_arg, den, num;
    mpfr_inits2(P, num_arg, den_arg, den, num, (mpfr_ptr)0);

    // den = cos(2π p)
    mpfr_mul(den_arg, C.two_pi, p, MPFR_RNDN);
    mpfr_cos(den, den_arg, MPFR_RNDN);

    mpfr_t eps;
    mpfr_init2(eps, P);
    mpfr_set_d(eps, 1e-30, MPFR_RNDN);

    if (mpfr_cmpabs(den, eps) > 0)
    {
        // Regular branch: direct ratio
        mpfr_t t1;
        mpfr_init2(t1, P);
        mpfr_sqr(t1, p, MPFR_RNDN);            // p^2
        mpfr_sub(t1, t1, p, MPFR_RNDN);        // p^2 − p
        mpfr_sub_d(t1, t1, 0.0625, MPFR_RNDN); // −1/16
        mpfr_mul(num_arg, C.two_pi, t1, MPFR_RNDN);
        mpfr_cos(num, num_arg, MPFR_RNDN);
        mpfr_div(out, num, den, MPFR_RNDN);
        mpfr_clear(t1);
    }
    else
    {
        // Near poles of cos(2π p): linearized proxy to avoid blowups
        mpfr_t q, half, d1, d2;
        mpfr_inits2(P, q, half, d1, d2, (mpfr_ptr)0);
        mpfr_set_d(half, 0.5, MPFR_RNDN);
        mpfr_sub_d(d1, p, 0.25, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_sub_d(d2, p, 0.75, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        if (mpfr_cmp(d1, d2) < 0)
        {
            mpfr_sub_d(q, p, 0.25, MPFR_RNDN);
            mpfr_sub(out, half, q, MPFR_RNDN);
        }
        else
        {
            mpfr_sub_d(q, p, 0.75, MPFR_RNDN);
            mpfr_add(out, half, q, MPFR_RNDN);
        }
        mpfr_clears(q, half, d1, d2, (mpfr_ptr)0);
    }

    mpfr_clear(eps);
    mpfr_clears(num_arg, den_arg, den, num, (mpfr_ptr)0);
}

/**
 * @brief Central differences for ψ''(p), ψ'''(p), ψ^{(6)}(p) using a 9-point stencil.
 *
 * @param d2  Output ψ''(p).
 * @param d3  Output ψ'''(p).
 * @param d6  Output ψ^{(6)}(p).
 * @param p   Input p.
 * @param P   MPFR precision in bits.
 * @param W   Precomputed weights (r=4).
 * @param h   Step size (small; chosen upstream).
 * @param C   Cached constants.
 */
static void psi_derivs_central(mpfr_t d2, mpfr_t d3, mpfr_t d6,
                               const mpfr_t p, mpfr_prec_t P,
                               const DiffWeights &W, const mpfr_t h,
                               const Constants &C)
{
    const int r = W.r;
    std::vector<mpfr_t> f(2 * r + 1);
    for (auto &v : f)
        mpfr_init2(v, P);

    // sample ψ at p + (j−r)h
    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t arg, step;
        mpfr_inits2(P, arg, step, (mpfr_ptr)0);
        mpfr_set(arg, p, MPFR_RNDN);
        mpfr_mul_si(step, h, j - r, MPFR_RNDN);
        mpfr_add(arg, arg, step, MPFR_RNDN);
        psi_p(f[j], arg, P, C);
        mpfr_clears(arg, step, (mpfr_ptr)0);
    }

    mpfr_set_zero(d2, 0);
    mpfr_set_zero(d3, 0);
    mpfr_set_zero(d6, 0);
    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t t;
        mpfr_init2(t, P);
        mpfr_mul(t, &W.w2[j], f[j], MPFR_RNDN);
        mpfr_add(d2, d2, t, MPFR_RNDN);
        mpfr_mul(t, &W.w3[j], f[j], MPFR_RNDN);
        mpfr_add(d3, d3, t, MPFR_RNDN);
        mpfr_mul(t, &W.w6[j], f[j], MPFR_RNDN);
        mpfr_add(d6, d6, t, MPFR_RNDN);
        mpfr_clear(t);
    }

    // normalize by h^m
    mpfr_t h2, h3, h6p;
    mpfr_inits2(P, h2, h3, h6p, (mpfr_ptr)0);
    mpfr_mul(h2, h, h, MPFR_RNDN);
    mpfr_mul(h3, h2, h, MPFR_RNDN);
    mpfr_mul(h6p, h3, h3, MPFR_RNDN);
    mpfr_div(d2, d2, h2, MPFR_RNDN);
    mpfr_div(d3, d3, h3, MPFR_RNDN);
    mpfr_div(d6, d6, h6p, MPFR_RNDN);

    mpfr_clears(h2, h3, h6p, (mpfr_ptr)0);
    for (auto &v : f)
        mpfr_clear(v);
}

// -----------------------------------------------------------------------------
// Gabcke correction R₁(t)
// -----------------------------------------------------------------------------

/**
 * @brief Gabcke-style correction term R₁(t) with optional higher-order components.
 *
 * Form:
 *   R₁ ≈ A0 * [ C0 + C1*(2π/t)^(1/2) + C2*(2π/t) ], sign (-1)^(N-1)
 * where A0 = (2π/t)^(1/4), C0 = ψ(p), and C1, C2 use ψ-derivatives.
 *
 * @param R        Output correction.
 * @param t        Height parameter.
 * @param P        MPFR precision in bits.
 * @param useC1C2  Enable C1/C2 terms if true.
 * @param W_opt    Optional prebuilt DiffWeights; if null, built locally.
 * @param C        Cached constants.
 */
void R_gabcke(mpfr_t R, const mpfr_t t, mpfr_prec_t P,
                     bool useC1C2,
                     DiffWeights *W_opt,
                     const Constants &C)
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    // A0 = (2π/t)^(1/4)
    mpfr_t A0;
    mpfr_init2(A0, P);
    mpfr_div(A0, C.two_pi, t, MPFR_RNDN);
    mpfr_rootn_ui(A0, A0, 4, MPFR_RNDN);

    // C0 = ψ(p)
    mpfr_t C0;
    mpfr_init2(C0, P);
    psi_p(C0, p, P, C);
    mpfr_set(R, C0, MPFR_RNDN);

    if (useC1C2)
    {
        // pick h relative to distance from 1/4 or 3/4, clamp to [1e-6, 1e-3]
        mpfr_t d1, quarter, threeq, d2, h;
        mpfr_inits2(P, d1, quarter, threeq, d2, h, (mpfr_ptr)0);
        mpfr_set_d(quarter, 0.25, MPFR_RNDN);
        mpfr_set_d(threeq, 0.75, MPFR_RNDN);
        mpfr_sub(d1, p, quarter, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_sub(d2, p, threeq, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        if (mpfr_cmp(d1, d2) < 0)
            mpfr_mul_d(h, d1, 0.1, MPFR_RNDN);
        else
            mpfr_mul_d(h, d2, 0.1, MPFR_RNDN);
        if (mpfr_cmp_d(h, 1e-6) < 0)
            mpfr_set_d(h, 1e-6, MPFR_RNDN);
        if (mpfr_cmp_d(h, 1e-3) > 0)
            mpfr_set_d(h, 1e-3, MPFR_RNDN);

        DiffWeights localW;
        DiffWeights *W = W_opt ? W_opt : &localW;
        if (W->w2.empty())
            W->build(P);

        mpfr_t psi2, psi3, psi6;
        mpfr_inits2(P, psi2, psi3, psi6, (mpfr_ptr)0);
        psi_derivs_central(psi2, psi3, psi6, p, P, *W, h, C);

        // C1 = -(1/(96 π^2)) ψ'''(p)
        mpfr_t C1;
        mpfr_init2(C1, P);
        mpfr_set(C1, C.pi2, MPFR_RNDN);
        mpfr_mul_ui(C1, C1, 96, MPFR_RNDN);
        mpfr_div(C1, psi3, C1, MPFR_RNDN);
        mpfr_neg(C1, C1, MPFR_RNDN);

        // C2 = ψ^{(6)}/(18432 π^4) + ψ''/(64 π^2)
        mpfr_t C2, term1, term2, denom;
        mpfr_inits2(P, C2, term1, term2, denom, (mpfr_ptr)0);
        mpfr_mul_ui(denom, C.pi4, 18432, MPFR_RNDN); // π^4 * 18432
        mpfr_div(term1, psi6, denom, MPFR_RNDN);
        mpfr_mul_ui(denom, C.pi2, 64, MPFR_RNDN); // π^2 * 64
        mpfr_div(term2, psi2, denom, MPFR_RNDN);
        mpfr_add(C2, term1, term2, MPFR_RNDN);

        // Apply scale factors (2π/t)^(1/2) and (2π/t)
        mpfr_t s12, s1, add;
        mpfr_inits2(P, s12, s1, add, (mpfr_ptr)0);
        mpfr_div(s1, C.two_pi, t, MPFR_RNDN);
        mpfr_sqrt(s12, s1, MPFR_RNDN);

        mpfr_mul(add, C1, s12, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);
        mpfr_mul(add, C2, s1, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);

        mpfr_clears(psi2, psi3, psi6, C1, C2, term1, term2, denom, s12, s1, add,
                    d1, d2, quarter, threeq, h, (mpfr_ptr)0);
    }

    // final amplitude and sign
    mpfr_mul(R, R, A0, MPFR_RNDN);
    if ((N % 2) == 0)
        mpfr_neg(R, R, MPFR_RNDN);

    mpfr_clears(a, p, A0, C0, (mpfr_ptr)0);
}

// -----------------------------------------------------------------------------
// RS main sum (not threaded)
// -----------------------------------------------------------------------------

/**
 * @brief Serial Riemann–Siegel main sum using rotation updates (compact form).
 *
 * @param Z     Output main-sum contribution.
 * @param t     Height parameter.
 * @param theta θ(t).
 * @param P     MPFR precision in bits.
 * @param C     Cached constants.
 */
static void Z_rs_main(mpfr_t Z, const mpfr_t t, const mpfr_t theta, mpfr_prec_t P, const Constants &C)
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    mpfr_t S;
    mpfr_init2(S, P);
    mpfr_set_zero(S, 0);

    mpfr_t cphi, sphi;
    mpfr_inits2(P, cphi, sphi, (mpfr_ptr)0);
    mpfr_sin_cos(sphi, cphi, theta, MPFR_RNDN);
    mpfr_add(S, S, cphi, MPFR_RNDN); // n=1

    mpfr_t dphi, sd, cd, new_c, new_s, tmp;
    mpfr_inits2(P, dphi, sd, cd, new_c, new_s, tmp, (mpfr_ptr)0);
    mpfr_t i_mp, sqrt_n;
    mpfr_inits2(P, i_mp, sqrt_n, (mpfr_ptr)0);

    for (long i = 1; i < N; ++i)
    {
        mpfr_set_ui(i_mp, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(dphi, 1, i_mp, MPFR_RNDN);
        mpfr_log1p(dphi, dphi, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        mpfr_set_ui(i_mp, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_n, i_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(S, S, tmp, MPFR_RNDN);
    }

    mpfr_mul_ui(Z, S, 2, MPFR_RNDN);
    mpfr_clears(a, p, S, cphi, sphi, dphi, sd, cd, new_c, new_s, tmp, i_mp, sqrt_n, (mpfr_ptr)0);
}

/**
 * @brief Combine RS main sum and Gabcke correction:
 *        Z_tot = Z_main + R_gabcke.
 *
 * Threading policy follows CFG().use_threaded_z unless Newton refinement
 * requests threads (CFG().use_threaded_newton), in which case RS runs serially
 * to avoid oversubscription.
 *
 * @param Ztot       Output total Z(t).
 * @param Zmain_out  Optional output for main-sum component (nullable).
 * @param R1_out     Optional output for correction component (nullable).
 * @param t          Height parameter.
 * @param P          MPFR precision in bits.
 * @param useC1C2    Enable C1/C2 terms in Gabcke correction.
 * @param W_opt      Optional prebuilt DiffWeights for reuse.
 * @param S          Spouge coefficients.
 * @param C          Cached constants.
 */
static void Z_rs_with_gabcke(
    mpfr_t Ztot,
    mpfr_t Zmain_out,
    mpfr_t R1_out,
    const mpfr_t t,
    mpfr_prec_t P,
    bool useC1C2,
    DiffWeights *W_opt,
    const Spouge &S,
    const Constants &C)
{
    // Guard: t must be positive
    if (mpfr_cmp_ui(t, 0) <= 0)
    {
        mpfr_set_nan(Ztot);
        if (Zmain_out)
            mpfr_set_nan(Zmain_out);
        if (R1_out)
            mpfr_set_nan(R1_out);
        LOGMP("Z_rs_with_gabcke: t must be > 0");
        return;
    }

    mpfr_t theta;
    mpfr_init2(theta, P);
    compute_theta_spouge(theta, t, P, S, C);

    mpfr_t Zmain;
    mpfr_init2(Zmain, P);

    if (CFG().use_threaded_newton)
    {
        // Let Newton have the threads; keep RS single-thread to reduce oversubscription.
        LOGMP("Disabling parallel RS–Gabcke to parallelize Newton...");
        Z_rs_main(Zmain, t, theta, P, C);
    }
    else
    {
        if (CFG().use_threaded_z)
        {
            LOGMP_DEBUG("Processing parallelized RS–Gabcke...");
            Z_rs_main_parallel(Zmain, t, theta, P, C, /*threads*/ 0, /*block_size*/ 0);
        }
        else
        {
            LOGMP_DEBUG("Processing single-thread RS–Gabcke...");
            Z_rs_main(Zmain, t, theta, P, C);
        }
    }

    mpfr_t R1;
    mpfr_init2(R1, P);
    R_gabcke(R1, t, P, useC1C2, W_opt, C);

    mpfr_add(Ztot, Zmain, R1, MPFR_RNDN);

    if (Zmain_out)
        mpfr_set(Zmain_out, Zmain, MPFR_RNDN);
    if (R1_out)
        mpfr_set(R1_out, R1, MPFR_RNDN);

    mpfr_clears(theta, Zmain, R1, (mpfr_ptr)0);
}

/**
 * @brief Compute Hardy’s Z(t) using:
 *   - θ(t) from Spouge-based log Γ (with reflection),
 *   - RS main sum via rotations (parallel if enabled),
 *   - Gabcke R₁ correction (with optional C₁, C₂).
 *
 * Honors CFG().mpfr_prec, CFG().use_threaded_z, CFG().use_threaded_newton.
 * Keeps the experimental loop to facilitate future batching, tests or vectorization.
 *
 * @param Z Output Z(t).
 * @param t Input height parameter.
 */
void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t)
{
    mpfr_prec_t P = CFG().mpfr_prec;

    // Future-friendly: keep vector/loop if batching more t’s later.
    std::vector<mpfr_srcptr> tvals;
    tvals.push_back(t);

    mpfr_set_default_prec(P);

    // Spouge & constants
    Spouge S;
    int a = spouge_a_for_prec(P);
    S.build(a, P);
    Constants C(P);

    DiffWeights W; // Builds lazily if C1/C2 enabled

    for (auto &ts : tvals)
    {
        mpfr_t tval;
        mpfr_init2(tval, P);
        mpfr_set(tval, ts, MPFR_RNDN);

        mpfr_t Ztot, Zmain, R1;
        mpfr_inits2(P, Ztot, Zmain, R1, (mpfr_ptr)0);

        bool useC1C2 = true; // full correction by default
        Z_rs_with_gabcke(Ztot, Zmain, R1, tval, P, useC1C2, &W, S, C);

        mpfr_set(Z, Ztot, MPFR_RNDN);

        mpfr_clears(tval, Ztot, Zmain, R1, (mpfr_ptr)0);
    }

    tvals.clear();
}
