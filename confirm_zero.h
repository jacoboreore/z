// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: confirm_zero.h
// Purpose:
//   Local spacing estimate and quick zero-likeness checks around a candidate t.
//   • est_spacing(t):   s ≈ π / log(t/2π) with safe clamps.
//   • use_threaded_for_t(t): decide whether RS main sum is large enough to
//                            benefit from threaded evaluation.
//   • confirm_zero_like(t, use_threaded_z): probe Z(t ± ε) and require a
//                            sign change to flag a likely zero nearby.
//
// -----------------------------------------------------------------------------

#pragma once

#include <mpfr.h>

#include "settings.h"
#include "z_calculator.h"       // compute_Z
#include "z_calculator_small.h" // compute_Z_small

/**
 * @brief Estimate local zero spacing s ≈ π / log(t / 2π).
 *
 * @param t MPFR value: imaginary part of s = 1/2 + i t.
 * @return Spacing estimate as double.
 *
 * @details
 *   - Guards against log(≤0) by clamping x = t/(2π) to at least 1e-12.
 *   - If log(x) ≤ 1e-12 after clamping, returns a safe fallback 0.25.
 */
static inline double est_spacing(const mpfr_t t)
{
    mpfr_t pi, two_pi, x, logx, s;
    mpfr_inits2(CFG().mpfr_prec, pi, two_pi, x, logx, s, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(x, t, two_pi, MPFR_RNDN);

    // clamp x ≥ 1e-12 to avoid log(≤0)
    if (mpfr_cmp_d(x, 1e-12) < 0)
        mpfr_set_d(x, 1e-12, MPFR_RNDN);

    mpfr_log(logx, x, MPFR_RNDN);
    // if logx too small/nonpositive, return a safe fallback spacing
    if (mpfr_cmp_d(logx, 1e-12) <= 0)
    {
        mpfr_clears(pi, two_pi, x, logx, s, (mpfr_ptr)0);
        return 0.25;
    }

    mpfr_div(s, pi, logx, MPFR_RNDN); // s = π / logx
    double out = mpfr_get_d(s, MPFR_RNDN);
    mpfr_clears(pi, two_pi, x, logx, s, (mpfr_ptr)0);
    return out;
}

/**
 * @brief Heuristic to decide if threaded RS evaluation is worthwhile at t.
 *
 * @param t MPFR value of the abscissa.
 * @return true if N = ⌊√(t / 2π)⌋ ≥ 8, else false.
 *
 * @details
 *   - Computes N as in the main-sum length of Riemann–Siegel (order of √(t/2π)).
 *   - Threshold 8 is empirical: above it, the threaded overhead tends to pay off.
 */
static inline bool use_threaded_for_t(const mpfr_t t)
{
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(CFG().mpfr_prec, pi, two_pi, tmp, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);

    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
    return (N >= 8);
}

/**
 * @brief Quick “zero-likeness” test: require a sign change across t ± ε.
 *
 * @param t              MPFR abscissa to test.
 * @param use_threaded_z Hint (0/1). If non-zero and use_threaded_for_t(t) is true,
 *                       use the full compute_Z path; otherwise use compute_Z_small.
 * @return true iff sign(Z(t + ε)) · sign(Z(t − ε)) < 0.
 *
 * @details
 *   - ε := max(1e-9, 0.2·est_spacing(t)) provides scale-aware probing while
 *     avoiding underflow at small t.
 *   - For the threaded path, temporaries are copied because some backends take
 *     non-const mpfr_t inputs.
 *
 * @note This is a heuristic screen, not a certified proof of a zero.
 */
static inline bool confirm_zero_like(const mpfr_t t, int use_threaded_z)
{
    const double s = est_spacing(t);
    const double eps = std::max(1e-9, 0.2 * s); // probe at 20% spacing, with a tiny floor

    mpfr_t tp, tm, Zp, Zm;
    mpfr_inits2(CFG().mpfr_prec, tp, tm, Zp, Zm, (mpfr_ptr)0);

    mpfr_add_d(tp, t, eps, MPFR_RNDN);
    mpfr_add_d(tm, t, -eps, MPFR_RNDN);

    if (use_threaded_z && use_threaded_for_t(t))
    {
        // compute_Z may modify inputs; pass temporaries
        mpfr_t tpc, tmc;
        mpfr_inits2(CFG().mpfr_prec, tpc, tmc, (mpfr_ptr)0);
        mpfr_set(tpc, tp, MPFR_RNDN);
        mpfr_set(tmc, tm, MPFR_RNDN);
        compute_Z(Zp, tpc);
        compute_Z(Zm, tmc);
        mpfr_clears(tpc, tmc, (mpfr_ptr)0);
    }
    else
    {
        compute_Z_small(Zp, tp);
        compute_Z_small(Zm, tm);
    }

    const int sp = mpfr_sgn(Zp);
    const int sm = mpfr_sgn(Zm);

    mpfr_clears(tp, tm, Zp, Zm, (mpfr_ptr)0);
    return (sp * sm) < 0; // opposite signs → likely a zero nearby
}
