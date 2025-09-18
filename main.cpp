#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <gmp.h>
#include <mpfr.h>

#include <algorithm>

#include "cli.h"
#include "settings.h"
#include "z_calculator.h" // compute_Z(...), rs_N_for_log(...), est_spacing(...)
#include "z_calculator_rs_gabcke.h"
#include "z_ball.h" // ZBallV, zballv_init/clear, compute_Z_ball_verbose, zballv_contains_zero
#include "refine_newton.h"
#include "refine_newton thread.h" // newton_refine_root(...)
// #include "confirm_zero.h"  // confirm_zero_like(...)
#include "z_walk.h"
#include "logs.h"
#include "utils.h"
#include "z_calculator_os.hpp"
#include "confirm_zero.h"
#include "crash_handlers.hpp"

int main(int argc, char **argv)
{
    // parse command line params cli.h
    Options opt;
    if (!parse_cli(argc, argv, opt))
        return 80;

    // load configuration, optional: file path relative to the executable conf/z.conf.
    // Defaults are in settings.h
    config_init("conf/z.conf");

    // global availability for logs
    g_main_log = CFG().rs_main_log;
    g_enable_log_debug = CFG().rs_log_debug;

    // print config if debug logs are enabled
    print_config(CFG());

    // print API version, MPFR, MPC, GMP
    // Current version as for 07/01/2025
    // Use these versions to avoid miscalculations or function crashing

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

    // Initializing GLOBAL best_y since we will use threads, sometimes a refined t pops up during the calculations
    best_init(CFG().mpfr_prec);

    // Odlyzko–Schönhage ad-hoc calculations precision.
    // The platform (Opertating System and compiler) would reduce precision due to allowed decimals.
    os_set_default_prec_bits(CFG().mpfr_prec);

    // Crash handler, printing a trace if calculations go wrong
    install_crash_handlers();

    mpfr_t t;
    mpfr_init2(t, CFG().mpfr_prec);

    if (mpfr_set_str(t, opt.t_str.c_str(), 10, MPFR_RNDN) != 0)
    {
        std::cerr << "[fatal] invalid number, exiting with 81 t\n";
        mpfr_clear(t);
        return 81;
    }

    // Calculations based on inputs
    if (opt.mode == Options::Mode::Z)
    {
        mpfr_t Z;
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z(Z, t);
        mpfr_clear(Z);
    }
    else if (opt.mode == Options::Mode::ValleyScanner)
    {
        y_confirmed(t);
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
    else if (opt.mode == Options::Mode::z_OS)
    {

        osvalley::OSConfig oc;

        const osvalley::Real t0 = mpfr_to_long_double_safe(t, MPFR_RNDN);
        auto blk = osvalley::scan_valley_OS(t0, oc);
        auto ev = osquick::analyze_block(blk,
                                         /*absZ_tol=*/(osvalley::Real)0.0000001,
                                         /*ratio_tol=*/(osvalley::Real)0.1);

        osquick::print_summary(ev);

        // // Simple policy for your valley scanner:
        // if (!ev.near)
        // {
        //     // you're “just far”: scale your coarse jump by ev.suggested_step
        //     // (clamp if you like: e.g. to 0.8*spacing)
        //     printf("far");

        //     const osvalley::Real max_j = 0.8L * ev.spacing;
        //     const osvalley::Real step =
        //         std::clamp(ev.suggested_step, -max_j, +max_j);
        //     // advance t by `step` and resample, or adjust your own stepper
        // }
        // else
        // {
        //     printf("close");
        //     // either |Z| small, sign flip in block, or dt small relative to spacing:
        //     // hand off to your sign-change/valley logic + high-prec Newton polish.
        // }

        // bool yeah = osquick::os_keep_for_valley(ev, 1e-2L,
        //                                         0.5L,
        //                                         0.5L,
        //                                         10.0L);

        // printf("\n----> OSSSSSSSSSS %d\n", yeah);

        // // build an OS block around your current t_cur
        // auto block = osvalley::scan_valley_OS(t0, oc);

        // // get spacing however you already do
        // osvalley::Real spacing = est_spacing(t); // π/log(t/2π)

        // // ask OS for a ratio-driven step
        // auto step = osvalley::os_ratio_step_from_block(block, t0, spacing,
        //                                                /*R_target=*/osvalley::Real(5.0));

        // // your valley-scanner policy:
        // if (step.flip_ahead)
        // {
        //     LOGMP_DEBUG("OS approximatting the valley");
        //     // hit the valley edge; refine around t_cur + step.dt_suggest
        // }
        // else if (step.ratio < 1.5)
        // {
        //     LOGMP_DEBUG("OS improve");
        //     // improvement too small -> increase step next time (or re-run with larger R_target)
        // }
        // else
        // {
        //     // take the suggested step
        //     // t0 += step.dt_suggest;
        // }

        // // 1) Make OS block
        // block = osvalley::scan_valley_OS(t0, oc);

        // // 2) Classify
        // osvalley::OSDecision d = osvalley::os_classify_block(block, oc);

        // // 3) Act
        // switch (d.verdict)
        // {
        // case osvalley::OSVerdict::ACCEPT:
        //     printf("\n Run Newton\n");
        //     // Run Newton (full or small RS) from d.t_seed
        //     // newton_refine_root(t_out, d.t_seed, …);
        //     break;
        // case osvalley::OSVerdict::REJECT:
        //     // Skip this t (junk); increase your stride aggressively
        //     printf("\n DROP\n");
        //     break;
        // case osvalley::OSVerdict::RESCAN:
        //     printf("\n More samples\n");
        //     // Keep but maybe: increase samples_per_spacing OR slide t0 ± spacing/2
        //     break;
        // }
    }
    else
    {
        printf("NOTHING TO DO");
    }

    mpfr_clear(t);

    mpfr_set_default_prec(CFG().mpfr_prec);

    return 0;
}
