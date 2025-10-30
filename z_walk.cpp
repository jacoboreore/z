// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// General Description:
// Core valley-walk routines to locate neighborhoods of zeros of the Riemann
// zeta function on the critical line. This module adapts the sampling step
// based on |Z| descent, optionally refines candidates via Newton, and emits
// CSV records for downstream plotting/analysis.
// -----------------------------------------------------------------------------

#include <iostream>
#include <mpfr.h>
#include <thread>
#include <mutex>
#include <cstdint>
#include <cmath>

#include "z_calculator.h"
#include "custom_y.h"
#include "settings.h"
#include "logs.h"
#include "confirm_zero.h"
#include "utils.h"
#include "refine_newton.h"
#include "refine_newton_thread.h"
#include "file.h"
#include "cli.h"

std::mutex print_mutex;
bool soft_landing = false;

/**
 * @brief Adaptive step controller for the valley walk.
 *
 * If the new |Z| is still large (>|tolerance_to_step|), the step is increased
 * aggressively (up to @cfg valley_max_step). Once inside the tolerance band,
 * the step resets for a soft landing and subsequently scales by the ratio
 * |Z_prev| / |Z_cur| if progress is too slow (ratio < @cfg valley_max_Zs_relation).
 *
 * @param z_prev     Previous Z(t) value.
 * @param z_cur      Current Z(t) value.
 * @param smart_step [in/out] Step in N-space; sign is preserved.
 */
static inline void smart_step(const mpfr_t z_prev,
                              const mpfr_t z_cur,
                              mpfr_t smart_step)
{
    const double tolerance_to_step = CFG().valley_tolerance_to_step;
    const double max_step = CFG().valley_max_step;
    const double max_Zs_relation = CFG().valley_max_Zs_relation;

    mpfr_t abs_prev, abs_cur, tmp, aggressive_factor, ts, ms, grow, mzr;
    mpfr_inits2(CFG().mpfr_prec, abs_prev, abs_cur, tmp, aggressive_factor, ts, ms, grow, mzr, (mpfr_ptr)0);

    mpfr_abs(abs_prev, z_prev, MPFR_RNDN);
    mpfr_abs(abs_cur, z_cur, MPFR_RNDN);

    mpfr_set_d(aggressive_factor, CFG().valley_aggressive_factor, MPFR_RNDN);
    mpfr_set_d(ms, max_step, MPFR_RNDN);
    mpfr_set_d(ts, tolerance_to_step, MPFR_RNDN);
    mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);
    mpfr_set_d(mzr, max_Zs_relation, MPFR_RNDN);

    // If |Z| is far from zero, increase the step aggressively (bounded by max_step).
    if (mpfr_cmp(ts, abs_cur) < 0)
    {
        LOGMP("Applying aggressive step: |Z| far from zero.");
        if (mpfr_cmp(ms, smart_step) > 0)
            mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
    }
    else
    {
        // Entering soft-landing mode: reset the step magnitude once.
        if (!soft_landing)
        {
            LOGMP("Resetting to soft landing step...");
            if (mpfr_sgn(smart_step) < 0)
            {
                mpfr_set(smart_step, grow, MPFR_RNDN);
                mpfr_neg(smart_step, smart_step, MPFR_RNDN);
            }
            else
            {
                mpfr_set(smart_step, grow, MPFR_RNDN);
            }
            soft_landing = true;
        }
        else
        {
            // If progress is slow (|Z_prev|/|Z_cur| < threshold), scale the step up.
            mpfr_div(tmp, abs_prev, abs_cur, MPFR_RNDN);
            LOGMP("Z_prev / Z_cur ratio = %.20Rf", tmp);

            if (mpfr_cmp(mzr, tmp) > 0)
            {
                LOGMP("Increasing step aggressively: from STEP %.20Rf", smart_step);
                mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
                LOGMP("to STEP %.20Rf", smart_step);
            }
        }
    }

    mpfr_clears(abs_prev, abs_cur, tmp, aggressive_factor, ts, ms, grow, mzr, (mpfr_ptr)0);
}

/**
 * @brief Compute a constant step for the mountain walk.
 *
 * If @cfg valley_constant_step == 0, the step is computed as sqrt(N) / divisor,
 * with divisor = @cfg valley_constant_step_div. Otherwise the configured constant
 * is used directly.
 *
 * @param constant_step_out [out] Resulting constant step.
 * @param N                 Current N value (used when computing by formula).
 */
void constant_step(mpfr_t constant_step_out, const mpfr_t N)
{
    const double cfg_const = CFG().valley_constant_step;
    const double cfg_div = CFG().valley_constant_step_div;

    if (cfg_const == 0.0)
    {
        LOGMP("Computing constant step from N = %.*Rg", CFG().rs_log_digits, N);

        mpfr_t calculated_step, divisor;
        mpfr_inits2(CFG().mpfr_prec, calculated_step, divisor, (mpfr_ptr)0);

        mpfr_set_d(divisor, cfg_div, MPFR_RNDN);
        mpfr_sqrt(calculated_step, N, MPFR_RNDN);
        mpfr_div(calculated_step, calculated_step, divisor, MPFR_RNDN);
        mpfr_set(constant_step_out, calculated_step, MPFR_RNDN);

        LOGMP("Constant step = %.12Rf\n", constant_step_out);
        mpfr_clears(calculated_step, divisor, (mpfr_ptr)0);
    }
    else
    {
        mpfr_set_d(constant_step_out, cfg_const, MPFR_RNDN);
        LOGMP("Constant step from configuration = %.12Rf", constant_step_out);
    }
}

/**
 * @brief Single descent from a starting height using adaptive stepping.
 *
 * Uses N-space stepping to decrease |Z| until a bounce (valley bottom)
 * is detected. This function does not perform Newton refinement.
 *
 * @param start_y [in/out] Starting t; updated internally during the walk.
 */
void from_y(mpfr_t start_y)
{
    mpfr_set_default_prec(CFG().mpfr_prec);

    const int max_iters = 1000000;
    int iter = 0;

    mpfr_t N, N_test, z_prev, z_cur, step, grow, t;
    mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, (mpfr_ptr)0);

    LOGMP("Parameters\n---> start_y %.12Rf\n", start_y);

    compute_Z(z_prev, start_y);
    mpfr_abs(z_prev, z_prev, MPFR_RNDN);
    LOGMP("|Z| at start: %.20Rf", z_prev);

    mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);

    compute_N_from_y(N, start_y);
    LOGMP("N from starting y: %.20Rf", N);

    mpfr_add(N_test, N, grow, MPFR_RNDN);
    LOGMP("N + grow: %.20Rf (grow: %.20Rf)", N_test, grow);

    compute_y_from_N(t, N_test);
    compute_Z(z_cur, t);
    mpfr_abs(z_cur, z_cur, MPFR_RNDN);
    LOGMP("|Z| at N + grow: %.20Rf", z_cur);

    // Choose direction: move where |Z| decreases.
    const bool go_down = (mpfr_cmp(z_cur, z_prev) < 0);

    // Initial step: ±grow in the chosen direction.
    mpfr_set_d(step, go_down ? CFG().valley_grow : -1.0 * CFG().valley_grow, MPFR_RNDN);

    LOGMP("Initial direction: %s", go_down ? "DOWN" : "UP");

    while (iter++ < max_iters)
    {
        // Move N by 'step'
        mpfr_add(N, N, step, MPFR_RNDN);

        // Evaluate |Z| at new N
        compute_y_from_N(t, N);
        compute_Z(z_cur, t);
        mpfr_abs(z_cur, z_cur, MPFR_RNDN);

        LOGMP("Scan: |Z| = %.50Rf  t = %.30Rf  STEP = %.10Rf", z_cur, t, step);
        LOGMP_DEBUG("N = %.30Rf", N);

        if (mpfr_cmp(z_cur, z_prev) < 0)
        {
            smart_step(z_prev, z_cur, step);
            mpfr_set(z_prev, z_cur, MPFR_RNDN);
            continue;
        }
        else
        {
            LOGMP("Valley bottom detected (bounce). |Z| = %.20Rf  t = %.20Rf", z_cur, t);
            break;
        }
    }

    mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, (mpfr_ptr)0);
}

/**
 * @brief Multi-pass valley walk with optional Newton refinement.
 * * Deprecating since spacing skips zeros *
 *
 * Iterates over spaced starting points (using a spacing heuristic) and, for each,
 * performs a valley walk. Upon detecting a bounce, triggers Newton refinement
 * (threaded or single-threaded based on configuration). If the refined |Z| is
 * below @cfg newton_tol_abs, the loop can skip further passes for that valley.
 *
 * @param start_y [in/out] Initial t; updated to the last refined t.
 */
void from_y_walk(mpfr_t start_y)
{
    mpfr_set_default_prec(CFG().mpfr_prec);

    const int max_iters = 100; // safety cap per inner walk
    int iter = 0;

    LOGMP("Valley walk initial parameters\n---> start_y %.12Rf\n", start_y);

    for (double valley_cycle = 0; valley_cycle < 5; valley_cycle++)
    {
        mpfr_t N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev;
        mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev, (mpfr_ptr)0);

        // Next y via spacing (spacing can occasionally miss zeros; proceed regardless).
        calculate_next_y_from_spacing(y_next, start_y);
        mpfr_set(start_y, y_next, MPFR_RNDN);

        LOGMP("--------------------------------");
        LOGMP("Starting cycle: next y = %.12Rf", start_y);

        for (double insist = 0; insist < CFG().valley_insistence; insist++)
        {
            iter = 0;

            compute_Z(z_prev, start_y);
            mpfr_abs(z_prev, z_prev, MPFR_RNDN);
            LOGMP("|Z| at start: %.20Rf", z_prev);

            mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);

            compute_N_from_y(N, start_y);
            LOGMP("N from starting y: %.20Rf", N);

            mpfr_add(N_test, N, grow, MPFR_RNDN);
            LOGMP("N + grow: %.20Rf (grow: %.20Rf)", N_test, grow);

            compute_y_from_N(t, N_test);
            compute_Z(z_cur, t);
            mpfr_abs(z_cur, z_cur, MPFR_RNDN);
            LOGMP("|Z| at N + grow: %.20Rf", z_cur);

            const bool go_down = (mpfr_cmp(z_cur, z_prev) < 0);

            // Initial step: ±grow in the chosen direction.
            mpfr_set_d(step, go_down ? CFG().valley_grow : -1.0 * CFG().valley_grow, MPFR_RNDN);
            LOGMP("Initial direction: %s", go_down ? "DOWN" : "UP");

            while (iter++ < max_iters)
            {
                mpfr_add(N, N, step, MPFR_RNDN);

                // Store last t
                mpfr_set(t_prev, t, MPFR_RNDN);

                compute_y_from_N(t, N);
                compute_Z(z_cur, t);
                mpfr_abs(z_cur, z_cur, MPFR_RNDN);

                LOGMP("Scan: |Z| = %.50Rf  t = %.30Rf  STEP = %.10Rf", z_cur, t, step);
                LOGMP_DEBUG("N = %.30Rf", N);

                if (mpfr_cmp(z_cur, z_prev) < 0)
                {
                    smart_step(z_prev, z_cur, step);
                    mpfr_set(z_prev, z_cur, MPFR_RNDN);
                    continue;
                }
                else
                {
                    LOGMP("Valley bottom detected. |Z| = %.20Rf  t = %.20Rf", z_cur, t);
                    LOGMP("Refining from t = %.20Rf ...", t_prev);

                    if (CFG().use_threaded_newton)
                    {
                        LOGMP("Using parallelized Newton...");
                        newton_refine_root_thread(
                            t_star,
                            t_prev,
                            CFG().newton_tol_abs,
                            CFG().newton_tol_step,
                            CFG().newton_max_iter,
                            CFG().use_threaded_z ? 1 : 0);
                    }
                    else
                    {
                        LOGMP("Using single-threaded Newton...");
                        newton_refine_root(
                            t_star,
                            t_prev,
                            CFG().newton_tol_abs,
                            CFG().newton_tol_step,
                            CFG().newton_max_iter,
                            CFG().use_threaded_z ? 1 : 0);
                    }

                    // Update start_y to refined t for subsequent passes.
                    mpfr_set(start_y, t_star, MPFR_RNDN);
                    mpfr_set(t, t_star, MPFR_RNDN);
                    break;
                }
            }

            // Check refined |Z|; if good enough, skip further passes for this valley.
            mpfr_t newton_tol_abs_mp, z_from_bounce;
            mpfr_inits2(CFG().mpfr_prec, newton_tol_abs_mp, z_from_bounce, (mpfr_ptr)0);
            mpfr_set_d(newton_tol_abs_mp, CFG().newton_tol_abs, MPFR_RNDN);

            compute_Z(z_from_bounce, start_y);
            mpfr_abs(z_from_bounce, z_from_bounce, MPFR_RNDN);

            LOGMPF("Refined t = %.*Rg  |Z| = %.*Rg",
                   CFG().rs_log_digits, start_y,
                   CFG().rs_log_digits, z_from_bounce);

            if (mpfr_cmp(newton_tol_abs_mp, z_from_bounce) < 0)
            {
                LOGMP("Valley scanner pass deemed sufficient.");
                insist = CFG().valley_insistence + 1; // force exit
            }
            else
            {
                LOGMP("Passing the valley scanner again...");
            }
            mpfr_clears(newton_tol_abs_mp, z_from_bounce, (mpfr_ptr)0);
        }

        mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev, (mpfr_ptr)0);
    }
}

/**
 * @brief Pure valley-walk loop for plotting and CSV output (no smart-step growth).
 *
 * Executes multiple valley-walk cycles, emitting per-scan samples and minima
 * markers for plotting. On a detected bounce, triggers Newton refinement and records the
 * refined (t, |Z|) to the process CSV. The step is derived from constant_step().
 *
 * @param start_y  [in/out] Starting t; updated between cycles.
 * @param steps    Number of valley cycles to perform.
 * @param opt      CLI/options (unused here; reserved for future behavior).
 */
void from_y_walk_pure_valley_walk(mpfr_t start_y, const double steps, Options opt)
{
    (void)opt; // reserved
    mpfr_set_default_prec(CFG().mpfr_prec);

    const int max_iters = 10000;
    int iter = 0;
    int main_iter = 0;
    bool ok = false;

    // Per-process output file names
    std::string process_file = zcsv::make_process_filename("output", "valley_walk");
    std::string walk_file = zcsv::make_process_filename("output", "valley_walk_plot");

    LOGMPF("Valley walk parameters\n---> start_y %.12Rf\n---> Valley count %.0f\n", start_y, steps);

    while (main_iter++ < steps)
    {
        mpfr_t N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce;
        mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce, (mpfr_ptr)0);

        LOGMP("--------------------------------");
        LOGMP("Starting cycle");

        compute_Z(z_prev, start_y);
        mpfr_abs(z_prev, z_prev, MPFR_RNDN);
        LOGMP("|Z| at start: %.20Rf", z_prev);

        mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);

        compute_N_from_y(N, start_y);
        LOGMP("N from starting y: %.20Rf", N);

        mpfr_add(N_test, N, grow, MPFR_RNDN);
        LOGMP("N + grow: %.20Rf (grow: %.20Rf)", N_test, grow);

        compute_y_from_N(t, N_test);
        compute_Z(z_cur, t);
        mpfr_abs(z_cur, z_cur, MPFR_RNDN);
        LOGMP("|Z| at N + grow: %.20Rf", z_cur);

        const bool go_down = (mpfr_cmp(z_cur, z_prev) < 0);

        // Initialize step from constant-step heuristic
        constant_step(step, N);

        iter = 0;
        if (!go_down)
        {
            LOGMP("Currently on an ascent; progressing until descent begins...");
            do
            {
                mpfr_add(N, N, step, MPFR_RNDN);
                compute_y_from_N(t, N);
                compute_Z(z_cur, t);
                mpfr_abs(z_cur, z_cur, MPFR_RNDN);

                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t, z_cur, N, nullptr, 0, 2, CFG().rs_log_digits);
                    if (!ok)
                        LOGMPF("ERROR writing CSV (walk).");
                }
                if (mpfr_cmp(z_cur, z_prev) > 0)
                {
                    mpfr_set(z_prev, z_cur, MPFR_RNDN);
                }

                // Start descent when |Z| decreases.
                if (mpfr_cmp(z_cur, z_prev) < 0)
                    break;

            } while (iter++ < max_iters);
        }

        LOGMP("Descending...");

        iter = 0;
        while (iter++ < max_iters)
        {
            mpfr_add(N, N, step, MPFR_RNDN);

            // Store last t
            mpfr_set(t_prev, t, MPFR_RNDN);

            compute_y_from_N(t, N);
            compute_Z(z_cur, t);
            mpfr_abs(z_cur, z_cur, MPFR_RNDN);

            LOGMP("Scan: |Z| = %.50Rf  t = %.30Rf  STEP = %.10Rf", z_cur, t, step);
            LOGMP("N = %.30Rf", N);

            if (mpfr_cmp(z_cur, z_prev) < 0)
            {
                mpfr_set(z_prev, z_cur, MPFR_RNDN);
                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t, z_cur, N, nullptr, 0, 2, CFG().rs_log_digits);
                    if (!ok)
                        LOGMPF("ERROR writing CSV (walk).");
                }
                continue;
            }
            else
            {
                LOGMP("Valley bottom detected. |Z| = %.20Rf  t = %.20Rf", z_cur, t);
                LOGMP("Refining from t = %.20Rf ...", t_prev);

                // Mark a minima at t_prev for plotting
                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t_prev, z_prev, N, nullptr, 1, 2, CFG().rs_log_digits);
                    if (!ok)
                        LOGMPF("ERROR writing CSV (minima).");
                }

                if (CFG().use_threaded_newton)
                {
                    LOGMP("Using parallelized Newton...");
                    newton_refine_root_thread(
                        t_refined,
                        t_prev,
                        CFG().newton_tol_abs,
                        CFG().newton_tol_step,
                        CFG().newton_max_iter,
                        CFG().use_threaded_z ? 1 : 0);
                }
                else
                {
                    LOGMP("Using single-threaded Newton...");
                    newton_refine_root(
                        t_refined,
                        t_prev,
                        CFG().newton_tol_abs,
                        CFG().newton_tol_step,
                        CFG().newton_max_iter,
                        CFG().use_threaded_z ? 1 : 0);
                }

                compute_Z(z_from_bounce, t_refined);
                mpfr_abs(z_from_bounce, z_from_bounce, MPFR_RNDN);

                LOGMPF("Refined t = %.*Rg  |Z| = %.*Rg",
                       CFG().rs_log_digits, t_refined,
                       CFG().rs_log_digits, z_from_bounce);

                ok = zcsv::append_record_csv(process_file, t_refined, z_from_bounce, N, nullptr, 2, 1, CFG().rs_log_digits);
                if (!ok)
                    LOGMPF("ERROR writing CSV (process).");

                // Prepare next cycle starting point
                mpfr_set(start_y, t, MPFR_RNDN);
                break;
            }
        }

        mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce, (mpfr_ptr)0);
    }
}
