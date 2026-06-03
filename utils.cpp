// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: best_state.cpp
// Purpose:
//   Thread-safe tracking of the best (smallest) observed |Z(t)| during scans,
//   plus simple remainder/radius utilities and a safe MPFR→long double caster.
//
// -----------------------------------------------------------------------------

#include <cassert>
#include <limits>
#include <cmath>
#include <mutex>
#include <mpfr.h>

#include "utils.h"
#include "logs.h"
#include "settings.h"

// -----------------------------------------------------------------------------
// Thread-safe singletons (state + mutex)
// -----------------------------------------------------------------------------

/**
 * @brief Obtain the process-wide mutex guarding BestZeroState.
 *
 * @return Reference to a function-local static mutex.
 */
static std::mutex &best_mutex()
{
    static std::mutex m;
    return m;
}

/**
 * @brief Obtain the process-wide BestZeroState singleton.
 *
 * The returned object is initialized by best_init(P) with the desired precision.
 *
 * @return Reference to a function-local static BestZeroState.
 */
BestZeroState &best_state()
{
    static BestZeroState s;
    return s;
}

/**
 * @brief Initialize (or reinitialize) the best-state storage with precision P.
 *
 * If already initialized at the same precision, this is a no-op.
 * If initialized at a different precision, previous MPFR objects are cleared
 * and re-created at the new precision.
 *
 * Sentinel values:
 *   best_absZ ← 9,999,999.0 (large) so any real |Z| improves it.
 *   best_t    ← +0
 *
 * @param P  MPFR working precision (in bits).
 */
void best_init(mpfr_prec_t P)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());

    if (s.initialized && s.P == P)
        return;

    if (s.initialized)
    {
        mpfr_clear(s.best_absZ);
        mpfr_clear(s.best_t);
        s.initialized = false;
        s.P = 0;
    }

    mpfr_init2(s.best_absZ, P);
    mpfr_init2(s.best_t, P);
    s.P = P;
    s.initialized = true;

    // Sentinels
    mpfr_set_d(s.best_absZ, 9'999'999.0, MPFR_RNDN); // start huge so any real |Z| wins
    mpfr_set_zero(s.best_t, 0);                      // +0
}

/**
 * @brief Reset the current “best” holders to their sentinel values.
 *
 * Requires a prior call to best_init(P).
 */
void best_reset()
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(P) must be called before best_reset()");
    mpfr_set_d(s.best_absZ, 9'999'999.0, MPFR_RNDN);
    mpfr_set_zero(s.best_t, 0);
}

// -----------------------------------------------------------------------------
// Update / Read
// -----------------------------------------------------------------------------

/**
 * @brief Compare |Z| against the current best and update if strictly smaller.
 *
 * Logs a message (left intact) when a new best is found.
 *
 * @param Z  Current Z(t) value (MPFR real).
 * @param t  Current t (height) associated with Z.
 */
void best_update(mpfr_srcptr Z, mpfr_srcptr t)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(P) must be called before best_update()");

    mpfr_t absZ;
    mpfr_init2(absZ, s.P);

    mpfr_abs(absZ, Z, MPFR_RNDN);
    if (mpfr_cmp(absZ, s.best_absZ) < 0)
    {
        LOGMP("(*) New better has been calculated t = %.*Rf  |Z| = %.*Rg",
              CFG().rs_log_digits, t,
              CFG().rs_log_digits, Z);

        mpfr_set(s.best_absZ, absZ, MPFR_RNDN);
        mpfr_set(s.best_t, t, MPFR_RNDN);
    }

    mpfr_clear(absZ);
}

/**
 * @brief Retrieve the current best pair (t_best, |Z|_best).
 *
 * The caller must ensure out_t and out_absZ have appropriate precision.
 *
 * @param[out] out_t     Best t so far.
 * @param[out] out_absZ  Smallest |Z| observed so far.
 */
void best_get(mpfr_t out_t, mpfr_t out_absZ)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(P) must be called before best_get()");

    mpfr_set(out_t, s.best_t, MPFR_RNDN);
    mpfr_set(out_absZ, s.best_absZ, MPFR_RNDN);
}

/**
 * @brief Convenience logger for the current best values (diagnostics).
 */
void get_best_t_so_far()
{
    mpfr_t best_t, best_absZ;
    mpfr_inits2(CFG().mpfr_prec, best_t, best_absZ, (mpfr_ptr)0);
    best_get(best_t, best_absZ);

    LOGMP("Best calculations t = %.*Rf |Z| = %.*Rg",
          CFG().rs_log_digits,
          best_t,
          CFG().rs_log_digits,
          best_absZ);

    mpfr_clears(best_t, best_absZ, (mpfr_ptr)0);
}

// -----------------------------------------------------------------------------
// Bounds & casting helpers
// -----------------------------------------------------------------------------

/**
 * @brief Simple diagnostic RS-tail bound R(t) ≈ C · t^{-1/4}.
 *
 * Uses the configurable constant CFG().rs_tail_C as a multiplier.
 *
 * @param[out] R  Output bound on |Remainder(t)|.
 * @param      t  Height.
 */
void remainder_bound(mpfr_t R, const mpfr_t t)
{
    const double C = CFG().rs_tail_C;
    mpfr_t r, tmp;
    mpfr_inits2(CFG().mpfr_prec, r, tmp, (mpfr_ptr)0);

    mpfr_sqrt(tmp, t, MPFR_RNDN); // √t
    mpfr_sqrt(r, tmp, MPFR_RNDN); // t^(1/4)
    mpfr_ui_div(R, 1, r, MPFR_RNDN);
    mpfr_mul_d(R, R, C, MPFR_RNDN);

    mpfr_clears(r, tmp, (mpfr_ptr)0);
}

/**
 * @brief Safe conversion from mpfr_t to long double with special-value handling.
 *
 * Behavior:
 *   - NaN → quiet_NaN()
 *   - +∞  → +infinity()
 *   - −∞  → −infinity()
 *   - Finite → mpfr_get_ld(x, rnd), with final sanity check and sign-aware clamp.
 *
 * @param x    MPFR real value.
 * @param rnd  Rounding mode (default MPFR_RNDN).
 * @return     Converted long double (with special cases handled).
 */
long double mpfr_to_long_double_safe(const mpfr_t x, mpfr_rnd_t rnd /*= MPFR_RNDN*/)
{
    if (!mpfr_number_p(x))
    {
        if (mpfr_nan_p(x))
            return std::numeric_limits<long double>::quiet_NaN();
        if (mpfr_inf_p(x) > 0)
            return std::numeric_limits<long double>::infinity();
        if (mpfr_inf_p(x) < 0)
            return -std::numeric_limits<long double>::infinity();
    }

    long double v = mpfr_get_ld(x, rnd);

    // Clamp if it overflowed the long double range during conversion.
    if (!std::isfinite(v))
    {
        return (mpfr_sgn(x) >= 0)
                   ? std::numeric_limits<long double>::infinity()
                   : -std::numeric_limits<long double>::infinity();
    }
    return v;
}
