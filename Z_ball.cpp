// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: z_ball.cpp
// Purpose:
//   Compute a validated “ball” for Hardy’s Z(t): a midpoint value plus
//   conservative radius contributions for the RS tail, rounding, and θ(t)
//   truncation. Chooses the small-ζ/η path or the RS–Gabcke path based on t.
//
// -----------------------------------------------------------------------------

#include <cmath>
#include <mpfr.h>

#include "z_ball.h"
#include "settings.h"
#include "z_calculator.h"
#include "z_theta.h"
#include "z_calculator_small.h"
#include "z_calculator_rs_gabcke.h"
#include "utils.h"

// -----------------------------------------------------------------------------
// Heuristics and component bounds
// -----------------------------------------------------------------------------

/**
 * @brief Decide whether to use the threaded RS–Gabcke path for a given t.
 *
 * Uses the RS length threshold N(t) ≥ 8 as a simple heuristic. For small t,
 * the small-ζ/η path is preferred.
 *
 * @param t  Height parameter.
 * @return true if RS–Gabcke (possibly threaded) should be used.
 */
static inline bool use_threaded_for_t(const mpfr_t t)
{
    return rs_N(t) >= 8;
}

/**
 * @brief Rounding-error radius for the RS main sum accumulation.
 *
 * Model:  k · 2√N · ulp(1)  with ulp(1) = 2^{1-P}, where P is MPFR precision.
 * This is a conservative bound that grows with √N and a tunable factor k.
 *
 * @param[out] out  Radius contribution due to rounding.
 * @param      N    RS length.
 * @param      k    Safety multiplier (default 4.0).
 */
static inline void rounding_radius(mpfr_t out, long N, double k = 4.0)
{
    mpfr_prec_t P = CFG().mpfr_prec;
    mpfr_t sN, ulp, kk;
    mpfr_inits2(P, sN, ulp, kk, (mpfr_ptr)0);

    // 2√N
    mpfr_set_ui(sN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sN, sN, MPFR_RNDN);
    mpfr_mul_ui(sN, sN, 2, MPFR_RNDN);

    // ulp(1) = 2^{1-P}
    mpfr_set_ui(ulp, 1, MPFR_RNDN);
    mpfr_div_2si(ulp, ulp, (long)P - 1, MPFR_RNDN);

    // scale by k
    mpfr_set_d(kk, k, MPFR_RNDN);

    // out = 2√N · ulp(1) · k
    mpfr_mul(out, sN, ulp, MPFR_RNDN);
    mpfr_mul(out, out, kk, MPFR_RNDN);

    mpfr_clears(sN, ulp, kk, (mpfr_ptr)0);
}

/**
 * @brief θ(t) truncation-radius estimate from Stirling-series depth.
 *
 * Compute θ_K and θ_{K+2} (capped at the available table), then use
 * |Δθ| · 2√N as an induced uncertainty on Z via the phase factor.
 *
 * @param[out] out  Radius contribution due to θ(t) truncation.
 * @param      t    Height parameter.
 */
static inline void theta_radius(mpfr_t out, const mpfr_t t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    const long N = rs_N(t);

    // Choose two truncations: K and K+2 (slightly deeper)
    int K  = CFG().theta_stirling_terms; // e.g., 8
    int K2 = K + 2;
    if (K2 > 9)
        K2 = 9; // cap to table

    mpfr_t thK, thK2, dth, sN, factor;
    mpfr_inits2(P, thK, thK2, dth, sN, factor, (mpfr_ptr)0);

    compute_theta_stirling_K(thK,  t, K);
    compute_theta_stirling_K(thK2, t, K2);

    // Δθ
    mpfr_sub(dth, thK, thK2, MPFR_RNDN);
    mpfr_abs(dth, dth, MPFR_RNDN);

    // factor = 2√N
    mpfr_set_ui(sN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sN, sN, MPFR_RNDN);
    mpfr_mul_ui(factor, sN, 2, MPFR_RNDN);

    // out = |Δθ| * 2√N
    mpfr_mul(out, dth, factor, MPFR_RNDN);

    mpfr_clears(thK, thK2, dth, sN, factor, (mpfr_ptr)0);
}

// -----------------------------------------------------------------------------
// Public: compute the validated ball around Z(t)
// -----------------------------------------------------------------------------

/**
 * @brief Compute a verbose Z-ball at height t.
 *
 * Fills:
 *   B.mid       ← midpoint evaluation of Z(t) via small-ζ/η or RS–Gabcke.
 *   B.rad_R     ← remainder (RS tail) bound (via remainder_bound()).
 *   B.rad_round ← rounding-error radius.
 *   B.rad_theta ← θ(t) truncation radius.
 *   B.rad_tot   ← total radius (R + round + theta).
 *
 * Selection:
 *   If rs_N(t) ≥ 8 → compute_Z_rs_gabcke(B.mid, t) (threading per CFG()).
 *   Otherwise     → compute_Z_small(B.mid, t).
 *
 * @param[out] B  Output structure holding midpoint and radii.
 * @param      t  Height where the ball is computed.
 */
void compute_Z_ball_verbose(ZBallV &B, const mpfr_t t)
{
    // Midpoint evaluation
    if (use_threaded_for_t(t))
    {
        // Local copy allows internal routines to adjust rounding modes safely.
        mpfr_t tc;
        mpfr_init2(tc, CFG().mpfr_prec);
        mpfr_set(tc, t, MPFR_RNDN);
        compute_Z_rs_gabcke(B.mid, tc);
        mpfr_clear(tc);
    }
    else
    {
        compute_Z_small(B.mid, t);
    }

    // Component radii
    remainder_bound(B.rad_R, t);
    rounding_radius(B.rad_round, rs_N(t), /*k=*/4.0);
    theta_radius(B.rad_theta, t);

    // Total radius
    mpfr_set(B.rad_tot, B.rad_R, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_round, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_theta, MPFR_RNDN);
}
