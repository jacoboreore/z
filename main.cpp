#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <gmp.h>
#include <mpfr.h>

#include "cli.h"
#include "settings.h"
#include "z_calculator.h"  // compute_Z(...), rs_N_for_log(...), est_spacing(...)
#include "z_calculator_rs_gabcke.h"
#include "z_ball.h"        // ZBallV, zballv_init/clear, compute_Z_ball_verbose, zballv_contains_zero
#include "refine_newton.h" // newton_refine_root(...)
//#include "confirm_zero.h"  // confirm_zero_like(...)
#include "z_walk.h"
#include "logs.h"
#include "utils.h"

void get_best_t_so_far(){
    
    // Print the Global best t across all points
    mpfr_t best_t, best_absZ;
    mpfr_inits2(CFG().mpfr_prec, best_t, best_absZ, (mpfr_ptr)0);
    best_get(best_t, best_absZ);

    std::ostringstream oss;
    oss << "BEST |Z| so far, |Z|=" << mp(best_absZ, CFG().rs_log_digits)
        << " at t=" << mp(best_t, CFG().rs_log_digits);

    LOGPF(oss.str(), CFG().rs_log_digits);

    mpfr_clears(best_t, best_absZ, (mpfr_ptr)0);
}

int main(int argc, char **argv)
{

    // print API version, MPFR, MPC, GMP
    // Current version as for 07/01/2005
    
    // MPFR runtime : 4.2.2
    // MPFR headers : 4.2.2
    // GMP  runtime : 6.3.0
    // GMP  headers : 6.3.0
    // MPC  runtime : (no mpc_get_version)
    // MPC  headers : 1.3.1

    api_version();

    // parse command line params cli.h
    Options opt;
    if (!parse_cli(argc, argv, opt))
        return 80;

    // load configuration, optional: file path relative to the executable conf/z.conf.  
    // Defaults are in settings.h
    config_init("conf/z.conf");
    std::cout << "\nCONFIG\n";
    print_config(CFG());

    // global availability for logs
    g_enable_log = CFG().rs_log_debug;

    // Initializing GLOBAL  best_y
    best_init(CFG().mpfr_prec);

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
        mpfr_t step;
        mpfr_init2(step, CFG().mpfr_prec);
        if (mpfr_set_str(step, opt.step_str.c_str(), 10, MPFR_RNDN) != 0)
        {
            std::cerr << "[fatal] mpfr_set_str failed for step\n";
            mpfr_clears(t, step, (mpfr_ptr)0);
            return 65;
        }

        y_confirmed(t);
        mpfr_clears(step, (mpfr_ptr)0);
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
        bool result = z_prefilter(t, CFG().prefilter_alpha, CFG().prefilter_tau);
        mpfr_printf("Keep t as close to a zero (int 1 = true): %d\n", result); 
    }
    else  if (opt.mode == Options::Mode::Refine)
    {
        mpfr_t t_star;
        mpfr_init2(t_star, CFG().mpfr_prec);

        bool ok = newton_refine_root(
            t_star, 
            t,
            CFG().newton_tol_abs,
            CFG().newton_tol_step,
            CFG().newton_max_iter,
            CFG().use_threaded_z ? 1 : 0);

        mpfr_printf("refine[%s]: t_in=%.15Rf  t_out=%.30Rf\n", ok ? "ok" : "stopped", t, t_star);

        ok = refine_in_stages(t_star, t, CFG().use_threaded_z);
        mpfr_printf("refine stages [%s]: t_in=%.15Rf  t_out=%.30Rf\n", ok ? "ok" : "stopped", t, t_star);

        get_best_t_so_far();

        mpfr_clear(t_star);
    } else
    {
        printf("NOTHING TO DO");
    }

    mpfr_clear(t);

    mpfr_set_default_prec(CFG().mpfr_prec);

    return 0;
}


