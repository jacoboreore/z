// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: main.cpp
// Purpose:
//   CLI entry-point coordinating configuration, logging, and execution modes
//   (Z, Valley Scanner/Walk, Z-ball, Prefilter, Newton refinement, progressive
//   refinement from CSV, Odlyzko–Schönhage helpers).
//
// -----------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <gmp.h>
#include <mpfr.h>
#include <algorithm>

#include "cli.h"
#include "settings.h"
#include "z_calculator.h"
#include "z_calculator_rs_gabcke.h"
#include "z_ball.h"
#include "refine_newton.h"
#include "refine_newton_thread.h"
#include "z_walk.h"
#include "logs.h"
#include "utils.h"
#include "z_calculator_os.hpp"
#include "confirm_zero.h"
#include "crash_handlers.hpp"
#include "file.h"
#include "custom_y.h"

/**
 * @brief Program entry point for Valley Scanner CLI.
 *
 * Parses command-line options, initializes configuration and logging, and
 * dispatches to the selected computation mode. This front controller wires up
 * precision, crash handlers, and runtime diagnostics (MPFR/GMP/MPC versions).
 *
 * ### Modes
 * - `Z`:
 *    Compute Hardy's Z(t) at a single @p t and exit.
 * - `Valley_Scanner`:
 *    Convert the input "y-like" value to an internal @p t (via @c from_y ),
 *    run scanner heuristics, and report the best-so-far @p t (global tracker).
 * - `Valley_Walk`:
 *    Perform a pure valley walk from the input using a given step size, then
 *    report the best-so-far @p t.
 * - `Z_ball`:
 *    Compute an interval/bound for Z via @c compute_Z_ball and print it.
 * - `Pre_filter`:
 *    Quickly test if @p t is worth refining (coarse filter).
 * - `Refine`:
 *    Newton-type refinement (single-thread or threaded) around @p t using
 *    tolerances from configuration; then report best-so-far.
 * - `Refine_Progressive`:
 *    Read triplets (t, |Z|, N) from CSV, iteratively tighten tolerances and
 *    refine each row, writing refined rows + spacing to an output CSV.
 * - `Z_OS`:
 *    Use Odlyzko–Schönhage helpers to analyze a local valley block and suggest
 *    a refined step, printing diagnostics.
 *
 * @param argc Argument count from the OS.
 * @param argv Argument vector (UTF-8 expected on POSIX; platform dependent on Windows).
 * @return int Exit code:
 *   - 0 on success
 *   - 80 on CLI parse failure
 *   - 81 on invalid numeric input for @p t
 *
 * @note Global precision (MPFR) is established from @c CFG().mpfr_prec .
 * @note Global “best t” tracker is initialized to support threaded refinement.
 * @note Logging verbosity is controlled by @c CFG().rs_log_debug / @c rs_main_log.
 */
int main(int argc, char **argv)
{
    // Parse command line (cli.h)
    Options opt;
    if (!parse_cli(argc, argv, opt))
        return 80;

    // Load configuration; defaults come from settings.h. Optional file: conf/z.conf
    config_init("conf/z.conf");

    // Make logging flags globally available
    g_main_log = CFG().rs_main_log;
    g_enable_log_debug = CFG().rs_log_debug;

    // Print config if debug logs are enabled
    print_config(CFG());

    // Print API/runtime versions (MPFR/GMP/MPC)
    api_version();

#ifdef OS_USE_MPFR
    mpfr::mpreal::set_default_prec(CFG().mpfr_prec);
    LOGMPF("Precision OS backend = MPFR, prec=%u bits", (unsigned)mpfr::mpreal::get_default_prec());
#else
    LOGMPF("Precision OS backend = long double, digits10=%d", std::numeric_limits<long double>::digits10);
#endif

    // Initialize global "best t" since threaded refiners may discover better points
    best_init(CFG().mpfr_prec);

    // Odlyzko–Schönhage ad-hoc calculations precision (platform aware)
    os_set_default_prec_bits(CFG().mpfr_prec);

    // Crash handler prints a stack trace on fatal errors
    install_crash_handlers();

    // Initialize logger; always tee to console
    init_logger("logs", "z_log_" + opt.job_id, /*also_stdout=*/true);
    LOGMPF("LOG FOR THIS PROCESS:\n%s", log_filename().c_str());

    // Prepare t
    mpfr_t t;
    mpfr_init2(t, CFG().mpfr_prec);

    if (opt.mode != Options::Mode::Refine_Progressive)
    {
        if (mpfr_set_str(t, opt.t_str.c_str(), 10, MPFR_RNDN) != 0)
        {
            std::cerr << "[fatal] invalid number, exiting with 81 t\n";
            mpfr_clear(t);
            return 81;
        }
    }

    // Dispatch by mode
    if (opt.mode == Options::Mode::Z)
    {
        /**
         * @brief Compute Z(t) once and exit.
         *
         * Uses @c compute_Z with current precision. No output printing is done
         * here besides internal logs; adapt as desired for UI integration.
         */
        mpfr_t Z;
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z(Z, t);
        mpfr_clear(Z);
    }
    else if (opt.mode == Options::Mode::Valley_Scanner)
    {
        /**
         * @brief Convert input "y" to @p t and run the scanner.
         *
         * After scanning, prints the best-so-far t via @c get_best_t_so_far .
         */
        from_y(t);
        get_best_t_so_far();
    }
    else if (opt.mode == Options::Mode::Valley_Walk)
    {
        /**
         * @brief Perform a pure valley walk given a step size.
         *
         * @details
         *   - Step size comes from @c opt.step_str .
         *   - @c from_y_walk_pure_valley_walk may emit walking logs/CSV when
         *     enabled in configuration.
         *   - After walking, best-so-far is logged.
         */
        const double steps = std::stod(opt.step_str.c_str());
        from_y_walk_pure_valley_walk(t, steps, opt);
        get_best_t_so_far();
    }
    else if (opt.mode == Options::Mode::Z_ball)
    {
        /**
         * @brief Compute a validated Z-ball around @p t and print midpoint Z.
         *
         * Uses @c compute_Z_ball to build a certified interval; prints the
         * midpoint value for quick inspection.
         */
        mpfr_t Z;
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z_ball(Z, t);
        mpfr_printf("Z_ball Z(%.*Rf) = %.20Re\n", 15, t, Z);
        mpfr_clear(Z);
    }
    else if (opt.mode == Options::Mode::Pre_filter)
    {
        /**
         * @brief Coarse pre-filter to decide whether @p t is worth refining.
         *
         * Prints 1 (true) or 0 (false).
         */
        const bool keep = z_prefilter(t);
        mpfr_printf("Keep t as close to a zero (int 1 = true): %d\n", keep ? 1 : 0);
    }
    else if (opt.mode == Options::Mode::Refine)
    {
        /**
         * @brief Newton-type refinement around @p t.
         *
         * Chooses between single-thread and threaded driver using
         * @c CFG().use_threaded_newton . Tolerances and iteration caps are
         * taken from configuration. The global best tracker may be updated
         * during the process; it is printed at the end.
         */
        mpfr_t t_star;
        mpfr_init2(t_star, CFG().mpfr_prec);

        bool ok = false;

        if (CFG().use_threaded_newton)
        {
            LOGMP("Using parallelized Newton...");
            ok = newton_refine_root_thread(
                t_star,
                t,
                CFG().newton_tol_abs,
                CFG().newton_tol_step,
                CFG().newton_max_iter,
                CFG().use_threaded_z ? 1 : 0);
        }
        else
        {
            LOGMP("Using single-thread Newton...");
            ok = newton_refine_root(
                t_star,
                t,
                CFG().newton_tol_abs,
                CFG().newton_tol_step,
                CFG().newton_max_iter,
                CFG().use_threaded_z ? 1 : 0);
        }

        (void)ok; // result is reflected by best tracker/logs
        get_best_t_so_far();

        mpfr_clear(t_star);
    }
    else if (opt.mode == Options::Mode::Refine_Progressive)
    {
        /**
         * @brief Progressive refinement of many rows from a CSV.
         *
         * Reads a CSV of triplets (t, |Z|, N) and for each row:
         *   1) Re-runs a short Newton burst with progressively tighter tolerances.
         *   2) Recomputes |Z|(t*) and N(t*), and the spacing vs. previous t*.
         *   3) Appends (t*, |Z|(t*), N(t*), spacing) to an output CSV specific
         *      to this process.
         *
         * Input file path is taken from @c opt.file_str .
         * Output file name is generated via @c zcsv::make_process_filename .
         */
        bool ok = false;

        // Output file for this process
        const std::string refine_file = zcsv::make_process_filename("output", "refined");
        LOGMPF("Refined values will be in this file: %s", refine_file.c_str());

        mpfr_t spacing, prev_t;
        mpfr_inits2(CFG().mpfr_prec, spacing, prev_t, (mpfr_ptr)0);

        mpfr_set_d(spacing, 0, MPFR_RNDN);
        mpfr_set_d(prev_t, 0, MPFR_RNDN);

        ok = zcsv::read_triplets_csv(
            opt.file_str,
            refine_file,
            spacing,
            prev_t,
            CFG().mpfr_prec,
            [](const mpfr_t t_from_file,
               const mpfr_t absZ,
               const mpfr_t N,
               size_t line_no,
               const std::string to_file,
               mpfr_t spacing,
               mpfr_t prev_t)
            {
                (void)N; // not needed during refinement burst

                long double newton_tol_abs = CFG().newton_tol_abs;
                long double newton_tol_step = CFG().newton_tol_step;

                LOGMPF("PROCESSING row=%zu t=%.*Rg |Z|=%.*Rg  N=%.*Rg",
                       line_no,
                       CFG().rs_log_digits, t_from_file,
                       CFG().rs_log_digits, absZ,
                       CFG().rs_log_digits, N);

                mpfr_t t_star, best_absZ, NN, working_t;
                mpfr_inits2(CFG().mpfr_prec, t_star, best_absZ, NN, working_t, (mpfr_ptr)0);

                // local working copy of t
                mpfr_set(working_t, t_from_file, MPFR_RNDN);

                // A short progressive-tightening sequence
                for (int i = 0; i < CFG().newton_max_iter; ++i)
                {
                    LOGMP_DEBUG("newton_tol_abs  = %.100Lf", newton_tol_abs);
                    LOGMP_DEBUG("newton_tol_step = %.100Lf", newton_tol_step);

                    if (CFG().use_threaded_newton)
                    {
                        LOGMP_DEBUG("Using parallelized Newton...");
                        newton_refine_root_thread(
                            t_star,
                            working_t,
                            newton_tol_abs,
                            newton_tol_step,
                            /*max_iter*/ 2,
                            CFG().use_threaded_z ? 1 : 0);
                    }
                    else
                    {
                        LOGMP_DEBUG("Using single-thread Newton...");
                        newton_refine_root(
                            t_star,
                            working_t,
                            newton_tol_abs,
                            newton_tol_step,
                            /*max_iter*/ 2,
                            CFG().use_threaded_z ? 1 : 0);
                    }

                    mpfr_set(working_t, t_star, MPFR_RNDN);

                    // tighten tolerances geometrically
                    newton_tol_abs *= 1e-2L;
                    newton_tol_step *= 1e-2L;
                }

                // Compute refined |Z| and N, then spacing vs previous
                compute_Z(best_absZ, working_t);
                mpfr_abs(best_absZ, best_absZ, MPFR_RNDN);
                compute_N_from_y(NN, working_t);

                LOGMPF("Writing to file refined t=%.*Rg N=%.*Rg |Z|=%.*Rg",
                       CFG().rs_log_digits, working_t,
                       CFG().rs_log_digits, NN,
                       CFG().rs_log_digits, best_absZ);

                // Δt (spacing)
                mpfr_sub(spacing, working_t, prev_t, MPFR_RNDN);
                mpfr_set(prev_t, working_t, MPFR_RNDN);

                const bool write_ok =
                    zcsv::append_record_csv(to_file, working_t, best_absZ, NN, spacing,
                                            /*minExp*/ 2, /*maxExp*/ 3, CFG().rs_log_digits);
                if (!write_ok)
                {
                    LOGMPF("Error writing CSV");
                }

                mpfr_clears(t_star, best_absZ, NN, working_t, (mpfr_ptr)0);
            });

        if (!ok)
        {
            LOGMPF("Can't read CSV, or completed with errors.\n"
                   "Ensure consistency if edition was manual, and confirm correct path");
        }

        mpfr_clears(spacing, prev_t, (mpfr_ptr)0);
    }
    else if (opt.mode == Options::Mode::Z_OS)
    {
        /**
         * @brief Odlyzko–Schönhage valley analysis (platform helper).
         *
         * Converts the MPFR @p t to a native type, scans a valley block, prints
         * a summary, and proposes a refined step. Also prints a high-precision
         * MPFR rendering of the new t for reproducibility.
         */
        osvalley::OSConfig oc;

        mpfr_t new_t_from_os, step_os;
        mpfr_init2(new_t_from_os, 256);
        mpfr_init2(step_os, 256);

        const osvalley::Real t0 = mpfr_to_long_double_safe(t, MPFR_RNDN);

        auto blk = osvalley::scan_valley_OS(t0, oc);
        auto ev = osquick::analyze_block(blk,
                                         /*absZ_tol=*/(osvalley::Real)1e-7,
                                         /*ratio_tol=*/(osvalley::Real)0.1);

        osquick::print_summary(ev);

        mpfr_set_ld(step_os, ev.suggested_step, MPFR_RNDN);

        const osvalley::Real t_new = t0 + ev.suggested_step;
        std::cout << "t_refined OS=" << t_new << "\n";

        mpfr_add(new_t_from_os, t, step_os, MPFR_RNDN);
        mpfr_printf("t os = %.*Rf\n", 30, new_t_from_os);

        mpfr_clears(new_t_from_os, step_os, (mpfr_ptr)0);
    }
    else
    {
        /**
         * @brief Fallback mode when no valid operation is requested.
         */
        std::printf("NOTHING TO DO");
    }

    mpfr_clear(t);

    LOGMPF("JOB_DONE");

    // Flush/close the logger gracefully
    shutdown_logger();

    return 0;
}
