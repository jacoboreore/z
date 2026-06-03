// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: z_small.cpp
// Description:
//   - Stable η(s) (Dirichlet eta) evaluation via binomial/Euler transform.
//   - Small-t Z(t) path using ζ(s) from η(s) and θ(t) from Stirling wrapper.
//
// Public entry points:
//   - void compute_Z_small(mpfr_t Z, const mpfr_t t)
//
// File-local helpers:
//   - static void zeta_via_eta(mpc_t zeta, const mpc_t s, int M_max, double tol)
// -----------------------------------------------------------------------------

#include <mpfr.h>
#include <mpc.h>
#include <vector>
#include <cmath>

#include "settings.h"
#include "z_theta.h"

// -----------------------------------------------------------------------------
// η(s) via binomial/Euler transform (stable)
//
//   η(s) = Σ_{m ≥ 0} [ 2^{-(m+1)} * Σ_{k=0}^m (-1)^k * C(m,k) * (k+1)^{-s} ]
//
// Precompute a[k] = (k+1)^{-s} using the stable recurrence:
//   a_{k+1} = a_k * exp( -s * log1p(1/(k+1)) )
//
// Accumulate outer m until the current m-term magnitude < tol after a warmup.
// M_max is a hard cap on the number of outer terms.
//
// Notes:
//   - This converges for Re(s) > 0 (including the critical line Re(s)=1/2).
//   - Denominator 1 - 2^{1-s} is handled at the end to recover ζ(s).
// -----------------------------------------------------------------------------

/**
 * @brief Compute ζ(s) via η(s) and the relation ζ(s) = η(s) / (1 - 2^{1-s}).
 *
 * @param[out] zeta   Complex result ζ(s).
 * @param[in]  s      Complex input s (mpc_t, precision from CFG().mpfr_prec).
 * @param[in]  M_max  Hard cap for outer m-sum (default 2000).
 * @param[in]  tol    Absolute tolerance for |mterm| after warmup (default 1e-40).
 */
static void zeta_via_eta(
    mpc_t zeta,
    const mpc_t s,
    int M_max = 2000,
    double tol = 1e-40)
{
    const int M_min = 12; // avoid early termination
    const mpfr_prec_t P = CFG().mpfr_prec;

    // --- Precompute a[k] = (k+1)^{-s} for k=0..M_max
    std::vector<mpc_t> a(static_cast<size_t>(M_max) + 1);
    for (int k = 0; k <= M_max; ++k) mpc_init2(a[static_cast<size_t>(k)], P);

    // a[0] = 1^{-s} = 1
    mpc_set_ui(a[0], 1, MPC_RNDNN);

    mpc_t ratio, tmpc;
    mpc_init2(ratio, P);
    mpc_init2(tmpc, P);

    mpfr_t ln1p;
    mpfr_init2(ln1p, P);

    for (int k = 0; k < M_max; ++k)
    {
        // ratio = exp( -s * log1p(1/(k+1)) )
        mpfr_set_ui(ln1p, static_cast<unsigned long>(k + 1), MPFR_RNDN);
        mpfr_ui_div(ln1p, 1, ln1p, MPFR_RNDN); // 1/(k+1)
        mpfr_log1p(ln1p, ln1p, MPFR_RNDN);
        mpc_mul_fr(tmpc, s, ln1p, MPC_RNDNN);
        mpc_neg(tmpc, tmpc, MPC_RNDNN);
        mpc_exp(ratio, tmpc, MPC_RNDNN);

        // a[k+1] = a[k] * ratio = (k+2)^{-s}
        mpc_mul(a[static_cast<size_t>(k + 1)], a[static_cast<size_t>(k)], ratio, MPC_RNDNN);
    }

    // --- Outer sum over m
    mpc_t eta, inner, mterm;
    mpc_init2(eta, P);
    mpc_set_ui(eta, 0, MPC_RNDNN);
    mpc_init2(inner, P);
    mpc_init2(mterm, P);

    mpfr_t Cm, pow2, abst; // C(m,0)=1 and pow2=2^(m+1)
    mpfr_inits2(P, Cm, pow2, abst, (mpfr_ptr)0);
    mpfr_set_ui(pow2, 2, MPFR_RNDN); // 2^(0+1) initially

    for (int m = 0; m <= M_max; ++m)
    {
        // inner = Σ_{k=0}^m (-1)^k C(m,k) a[k]
        mpfr_set_ui(Cm, 1, MPFR_RNDN);
        mpc_set_ui(inner, 0, MPC_RNDNN);

        for (int k = 0; k <= m; ++k)
        {
            // mterm = C(m,k) * a[k] * ((-1)^k)
            mpc_mul_fr(mterm, a[static_cast<size_t>(k)], Cm, MPC_RNDNN);
            if (k & 1) mpc_neg(mterm, mterm, MPC_RNDNN);
            mpc_add(inner, inner, mterm, MPC_RNDNN);

            if (k < m)
            {
                // Cm *= (m-k)/(k+1)
                mpfr_mul_ui(Cm, Cm, static_cast<unsigned long>(m - k), MPFR_RNDN);
                mpfr_div_ui(Cm, Cm, static_cast<unsigned long>(k + 1), MPFR_RNDN);
            }
        }

        // mterm = inner / 2^(m+1)
        mpc_div_fr(mterm, inner, pow2, MPC_RNDNN);
        mpc_add(eta, eta, mterm, MPC_RNDNN);

        // stop when |mterm| < tol (after warmup)
        mpc_abs(abst, mterm, MPFR_RNDN);
        if (m >= M_min && mpfr_cmp_d(abst, tol) < 0)
            break;

        // pow2 *= 2  (advance outer factor)
        mpfr_mul_ui(pow2, pow2, 2, MPFR_RNDN);
    }

    // --- ζ(s) = η(s) / (1 - 2^{1-s})
    mpfr_t ln2;
    mpfr_init2(ln2, P);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_const_log2(ln2, MPFR_RNDN);
#else
    mpfr_set_ui(ln2, 2, MPFR_RNDN);
    mpfr_log(ln2, ln2, MPFR_RNDN);
#endif

    // tmpc = 2^{1-s} = 2 * exp(-s * ln 2)
    mpc_mul_fr(tmpc, s, ln2, MPC_RNDNN);
    mpc_neg(tmpc, tmpc, MPC_RNDNN);
    mpc_exp(tmpc, tmpc, MPC_RNDNN);
    mpc_mul_ui(tmpc, tmpc, 2, MPC_RNDNN);

    // zeta = eta / (1 - 2^{1-s})
    mpc_ui_sub(mterm, 1, tmpc, MPC_RNDNN); // mterm = 1 - 2^{1-s}
    // (If |mterm| is extremely small, callers may want to raise precision or switch paths.)
    mpc_div(zeta, eta, mterm, MPC_RNDNN);

    // cleanup
    for (auto &v : a) mpc_clear(v);
    mpc_clear(ratio);
    mpc_clear(tmpc);
    mpc_clear(eta);
    mpc_clear(inner);
    mpc_clear(mterm);
    mpfr_clears(ln1p, Cm, pow2, abst, ln2, (mpfr_ptr)0);
}

// -----------------------------------------------------------------------------
// Z(t) = Re( e^{i θ(t)} ζ(1/2 + i t) )
//
// Uses η(s) (above) to evaluate ζ(s) on the critical line, and θ(t) from the
// Stirling-based θ wrapper. Intended for small/medium t, or as a fallback path.
// -----------------------------------------------------------------------------

/**
 * @brief Compute Hardy’s Z(t) for moderate t via η(s) and θ(t).
 *
 * Z(t) = Re( e^{i θ(t)} ζ(1/2 + i t) ) = Re( (cosθ + i sinθ)(zr + i zi) )
 *      = zr·cosθ − zi·sinθ
 *
 * @param[out] Z  Real result Z(t).
 * @param[in]  t  Height parameter (MPFR).
 */
void compute_Z_small(mpfr_t Z, const mpfr_t t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;

    mpc_t s, z;
    mpc_init2(s, P);
    mpc_init2(z, P);

    mpfr_t half, theta, cth, sth, zr, zi, tmp;
    mpfr_inits2(P, half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);

    // s = 1/2 + i t
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpc_set_fr_fr(s, half, t, MPC_RNDNN);

    // ζ(s) via Dirichlet η(s)
    zeta_via_eta(z, s, /*M_max*/ 2000, /*tol*/ 1e-40);

    // θ(t): Stirling-based wrapper (terms controlled by CFG().theta_stirling_terms)
    compute_theta_stirling(theta, t);

    // Z(t) = zr*cosθ − zi*sinθ
    mpfr_cos(cth, theta, MPFR_RNDN);
    mpfr_sin(sth, theta, MPFR_RNDN);
    mpfr_set(zr, mpc_realref(z), MPFR_RNDN);
    mpfr_set(zi, mpc_imagref(z), MPFR_RNDN);

    mpfr_mul(Z, zr, cth, MPFR_RNDN);
    mpfr_mul(tmp, zi, sth, MPFR_RNDN);
    mpfr_sub(Z, Z, tmp, MPFR_RNDN);

    mpfr_clears(half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);
    mpc_clear(s);
    mpc_clear(z);
}
