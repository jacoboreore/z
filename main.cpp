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

int main(int argc, char **argv)
{
    // parse command line params cli.h
    Options opt;
    if (!parse_cli(argc, argv, opt))
        return 80;

    // load configuration, optional: relative path conf/z.conf.
    // Defaults are in settings.h
    config_init("conf/z.conf");

    // global availability for logs
    g_main_log = CFG().rs_main_log;
    g_enable_log_debug = CFG().rs_log_debug;

    // print config if debug logs are enabled
    print_config(CFG());

    // print API version, MPFR, MPC, GMP
    // Current version as for 07/01/2025

    // MPFR runtime : 4.2.2
    // MPFR headers : 4.2.2
    // GMP  runtime : 6.3.0
    // GMP  headers : 6.3.0
    // MPC  runtime : (no mpc_get_version)
    // MPC  headers : 1.3.1
    api_version();

#ifdef OS_USE_MPFR
    mpfr::mpreal::set_default_prec(CFG().mpfr_prec);
    LOGMPF("Precision OS backend = MPFR, prec=%u bits", (unsigned)mpfr::mpreal::get_default_prec());
#else
    LOGMPF("Precision OS backend = long double, digits 10=%d", std::numeric_limits<long double>::digits10);
#endif

    // Initializing GLOBAL t since we will use threads, sometimes a refined t pops up during the calculations
    best_init(CFG().mpfr_prec);

    // Odlyzko–Schönhage ad-hoc calculations precision.
    // The platform (Opertating System and compiler) would reduce precision due to allowed decimals.
    os_set_default_prec_bits(CFG().mpfr_prec);

    // Crash handler, printing a trace if calculations go wrong
    install_crash_handlers();

    // Init the logger, always print to the console
    init_logger("logs", "z_log_" + opt.job_id, /*also_stdout=*/true);
    LOGMPF("LOG FOR THIS PROCESS:\n%s", log_filename().c_str());

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

    // Calculations based on inputs
    if (opt.mode == Options::Mode::Z)
    {
        mpfr_t Z;
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z(Z, t);
        mpfr_clear(Z);
    }
    else if (opt.mode == Options::Mode::Valley_Scanner)
    {
        from_y(t);
        get_best_t_so_far();
    }
    else if (opt.mode == Options::Mode::Valley_Walk)
    {
        double steps = std::stod(opt.step_str.c_str());
        from_y_walk_pure_valley_walk(t, steps, opt);
        get_best_t_so_far();
    }
    else if (opt.mode == Options::Mode::Z_ball)
    {
        mpfr_t Z;
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z_ball(Z, t);
        mpfr_printf("Z_ball Z(%.*Rf) = %.20Re\n", 15, t, Z);
        mpfr_clear(Z);
    }
    else if (opt.mode == Options::Mode::Pre_filter)
    {
        bool result = z_prefilter(t);
        mpfr_printf("Keep t as close to a zero (int 1 = true): %d\n", result);
    }
    else if (opt.mode == Options::Mode::Refine)
    {
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
            LOGMP("Using Single thread Newton...");
            ok = newton_refine_root(
                t_star,
                t,
                CFG().newton_tol_abs,
                CFG().newton_tol_step,
                CFG().newton_max_iter,
                CFG().use_threaded_z ? 1 : 0);
        }

        get_best_t_so_far();

        mpfr_clear(t_star);
    }
    else if (opt.mode == Options::Mode::Refine_Progressive)
    {
        bool ok = false;

        // get the file name for this process
        std::string refine_file = zcsv::make_process_filename("output", "refined_" + opt.job_id);
        LOGMPF("Refined values will be in this file: %s", refine_file.c_str());

        mpfr_t spacing, prev_t;
        mpfr_inits2(CFG().mpfr_prec, spacing, prev_t, (mpfr_ptr)0);

        mpfr_set_d(spacing, 0, MPFR_RNDN);
        mpfr_set_d(prev_t, 0, MPFR_RNDN);
        
        // open calculations file
        ok = zcsv::read_triplets_csv(
            opt.file_str,
            refine_file,
            spacing,
            prev_t,
            CFG().mpfr_prec,
            [](const mpfr_t t_from_file, const mpfr_t absZ, const mpfr_t N, size_t line_no, const std::string to_file, mpfr_t spacing, mpfr_t prev_t)
            {
                long double newton_tol_abs = CFG().newton_tol_abs;
                long double newton_tol_step = CFG().newton_tol_step;

                LOGMPF("PROCESSING row=%zu t=%.*Rg |Z|=%.*Rg  N=%.*Rg",
                       line_no,
                       CFG().rs_log_digits, t_from_file,
                       CFG().rs_log_digits, absZ,
                       CFG().rs_log_digits, N);

                mpfr_t t_star, best_t, best_absZ, NN, working_t;
                mpfr_inits2(CFG().mpfr_prec, t_star, best_t, best_absZ, NN, working_t, (mpfr_ptr)0);

                // a copy from t
                mpfr_set(working_t, t_from_file, MPFR_RNDN);

                for (int i = 0; i < CFG().newton_max_iter; i++)
                {
                    LOGMP_DEBUG("newton_tol_abs = %.50Lf", newton_tol_abs);
                    LOGMP_DEBUG("newton_tol_step = %.50Lf", newton_tol_step);
                    if (CFG().use_threaded_newton)
                    {
                        LOGMP_DEBUG("Using parallelized Newton...");
                        newton_refine_root_thread(
                            t_star,
                            working_t,
                            newton_tol_abs,
                            newton_tol_step,
                            2,
                            CFG().use_threaded_z ? 1 : 0);
                    }
                    else
                    {
                        LOGMP_DEBUG("Using Single thread Newton...");
                        newton_refine_root(
                            t_star,
                            working_t,
                            newton_tol_abs,
                            newton_tol_step,
                            2,
                            CFG().use_threaded_z ? 1 : 0);
                    }

                    mpfr_set(working_t, t_star, MPFR_RNDN);

                    // increase the precision to refine
                    newton_tol_abs = newton_tol_abs * (1e-2);
                    newton_tol_step = newton_tol_step * (1e-2);
                }

                // write the redined value
                compute_Z(best_absZ, working_t);
                mpfr_abs(best_absZ, best_absZ, MPFR_RNDN);
                compute_N_from_y(NN, working_t);

                LOGMPF("Writting to file refined t=%.*Rg N=%.*Rg |Z|=%.*Rg",
                       CFG().rs_log_digits, working_t,
                       CFG().rs_log_digits, NN,
                       CFG().rs_log_digits, best_absZ);

                // calculate spacing 
                mpfr_sub(spacing, working_t, prev_t, MPFR_RNDN);
                mpfr_set(prev_t, working_t, MPFR_RNDN);
                
                bool write_ok = zcsv::append_record_csv(to_file, working_t, best_absZ, NN, spacing, 2, 3, CFG().rs_log_digits);
                if (!write_ok)
                {
                    LOGMPF("Error writing csv");
                }

                mpfr_clears(t_star, best_t, best_absZ, NN, working_t, (mpfr_ptr)0);
            });

        if (!ok)
        {
            LOGMPF("Can't read CSV, or completed with errors.\nEnsure consistency if edition was manual, and confirm correct path");
        }

        mpfr_clears(spacing, (mpfr_ptr)0);
    }
    else if (opt.mode == Options::Mode::Z_OS)
    {

        osvalley::OSConfig oc;
        mpfr_t new_t_from_os, step_os;

        mpfr_init2(new_t_from_os, 256);
        mpfr_init2(step_os, 256);

        const osvalley::Real t0 = mpfr_to_long_double_safe(t, MPFR_RNDN);
        auto blk = osvalley::scan_valley_OS(t0, oc);
        auto ev = osquick::analyze_block(blk,
                                         /*absZ_tol=*/(osvalley::Real)0.0000001,
                                         /*ratio_tol=*/(osvalley::Real)0.1);

        osquick::print_summary(ev);

        mpfr_set_ld(step_os, ev.suggested_step, MPFR_RNDN);

        osvalley::Real t_new = t0 + ev.suggested_step;
        std::cout << "t_refined OS=" << t_new << "\n";

        mpfr_add(new_t_from_os, t, step_os, MPFR_RNDN);

        mpfr_printf("t os = %.*Rf\n", 30, new_t_from_os);
    }
    else
    {
        printf("NOTHING TO DO");
    }

    mpfr_clear(t);

    mpfr_set_default_prec(CFG().mpfr_prec);

    // sanity close for the log
    shutdown_logger();

    return 0;
}
