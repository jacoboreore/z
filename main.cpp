#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <gmp.h>
#include <mpfr.h>

#include "cli.h"
#include "settings.h"
#include "z_calculator.h"   // compute_Z(...), rs_N_for_log(...), est_spacing(...)
#include "z_ball.h"         // ZBallV, zballv_init/clear, compute_Z_ball_verbose, zballv_contains_zero
#include "refine_newton.h"  // newton_refine_root(...)
#include "confirm_zero.h"   // confirm_zero_like(...)

int main(int argc, char** argv) {

    // print API version, MPFR, MPC, GMP
    api_version();
    
    // parse command line params cli.h
    Options opt;
    if (!parse_cli(argc, argv, opt)) return 80;

    // load configuration, optional: pass a file path.  Defaults are in settings.h
    config_init("conf/z.conf");
    std::cout << "\nCONFIG\n";
    print_config(CFG());

    mpfr_t t; 
    mpfr_init2(t, CFG().mpfr_prec);

    if (mpfr_set_str(t, opt.t_str.c_str(), 10, MPFR_RNDN) != 0) {
        std::cerr << "[fatal] invalid number, exiting with 81 t\n"; 
        mpfr_clear(t); 
        return 81;
    }

    if (opt.mode == Options::Mode::Z) {
        
        mpfr_t Z; 
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z(Z, t);
        mpfr_printf("Z(%.*Rf) = %.30Re\n", 15, t, Z);
        mpfr_clear(Z);

    } else if (opt.mode == Options::Mode::ValleyScanner) {
        
        mpfr_t step; 
        mpfr_init2(step, CFG().mpfr_prec);
        if (mpfr_set_str(step, opt.step_str.c_str(), 10, MPFR_RNDN) != 0) {
            std::cerr << "[fatal] mpfr_set_str failed for step\n";
            mpfr_clears(t, step, (mpfr_ptr)0); return 65;
        }
        
        // TODO: not implemented
        mpfr_clears(step, (mpfr_ptr)0);

    } else if (opt.mode == Options::Mode::Z_ball) {
        
        mpfr_t Z; 
        mpfr_init2(Z, CFG().mpfr_prec);
        compute_Z_ball(Z, t);
        mpfr_printf("Z_ball Z(%.*Rf) = %.20Re\n", 15, t, Z);
        mpfr_clear(Z);
        
    } else { // Refine
        
        mpfr_t t_star; 
        mpfr_init2(t_star, CFG().mpfr_prec);
        
        bool ok = newton_refine_root(
            t_star, t,
            CFG().newton_tol_abs,
            CFG().newton_tol_step,
            CFG().newton_max_iter,
            CFG().use_threaded_z ? 1 : 0
        );

        mpfr_printf("refine[%s]: t_in=%.15Rf  t_out=%.30Rf\n", ok ? "ok" : "stopped", t, t_star);
        mpfr_clear(t_star);
    }

    mpfr_clear(t);

    mpfr_set_default_prec(CFG().mpfr_prec);

    return 0;
}

// #include "z_calculator_small.h"
// #include "z_calculator_thread.h"
// #include "z_theta.h"
// #include "settings.h"
// #include "cli.h"
// #include <cstdio>
// #include <cmath>
// #include <iostream>


// static inline void print_val(const char* tag, const mpfr_t x, int d=30){
//     std::printf("%-14s = ", tag);
//     mpfr_out_str(stdout, 10, d, x, MPFR_RNDN);
//     std::printf("\n");
// }

// void rs_sum_direct(mpfr_t out, const mpfr_t t) {
//     mpfr_t theta, term, n, ln_n, sqrt_n, phi, cphi;
//     mpfr_inits2(CFG().mpfr_prec, theta, term, n, ln_n, sqrt_n, phi, cphi, (mpfr_ptr)0);

//     compute_theta_stirling(theta, t, CFG().mpfr_prec, /*K=*/CFG().theta_stirling_terms);

//     // N = floor(sqrt(t / 2π))
//     mpfr_t pi, two_pi, tmp;
//     mpfr_inits2(CFG().mpfr_prec, pi, two_pi, tmp, (mpfr_ptr)0);
//     mpfr_const_pi(pi, MPFR_RNDN);
//     mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
//     mpfr_div(tmp, t, two_pi, MPFR_RNDN);
//     mpfr_sqrt(tmp, tmp, MPFR_RNDN);
//     long N = mpfr_get_si(tmp, MPFR_RNDN);
//     if (N < 1) N = 1;
//     mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);

//     mpfr_set_zero(out, 0);
//     for (long i = 1; i <= N; ++i) {
//         mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
//         mpfr_log(ln_n, n, MPFR_RNDN);
//         mpfr_mul(phi, t, ln_n, MPFR_RNDN);
//         mpfr_sub(phi, theta, phi, MPFR_RNDN);
//         mpfr_cos(cphi, phi, MPFR_RNDN);
//         mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
//         mpfr_div(term, cphi, sqrt_n, MPFR_RNDN);
//         mpfr_add(out, out, term, MPFR_RNDN);
//     }

//     // Z ≈ 2*sum  (optionally half-weight the last term if you’ve chosen that convention)
//     mpfr_mul_ui(out, out, 2, MPFR_RNDN);
// }

// void verify_Z_at(char* t_str) {
//     mpfr_t t, Zthr, Zs, diff;
//     mpfr_inits2(CFG().mpfr_prec, t, Zthr, Zs, diff, (mpfr_ptr)0);
//     mpfr_set_str(t, t_str, 10, MPFR_RNDN);

//     compute_Z_threaded(Zthr, t);
//     compute_Z_small(Zs, t);
//     mpfr_sub(diff, Zthr, Zs, MPFR_RNDN);

//     std::cout << "\n[Z verification]\n";
//     print_val("t            ", t, 20);
//     print_val("Z_threaded   ", Zthr, CFG().rs_log_digits);
//     print_val("Z_small      ", Zs,   CFG().rs_log_digits);
//     print_val("Δ (thr-small)", diff, CFG().rs_log_digits);

//     // print rough t^{-1/4} bound with your C
//     mpfr_t tmp; mpfr_init2(tmp, CFG().mpfr_prec);
//     mpfr_sqrt(tmp, t, MPFR_RNDN);
//     mpfr_sqrt(tmp, tmp, MPFR_RNDN);
//     mpfr_ui_div(tmp, 1, tmp, MPFR_RNDN);
//     double b = CFG().rs_tail_C * mpfr_get_d(tmp, MPFR_RNDN);
//     std::cout << "~bound C t^-1/4 = "  << b << "\n";

//     mpfr_clears(t, Zthr, Zs, diff, tmp, (mpfr_ptr)0);
// }

// int main(int argc, char** argv) {

    
//     // print API version, MPFR, MPC, GMP
//     api_version();
    
//     // parse command line params cli.h
//     Options opt;
//     if (!parse_cli(argc, argv, opt)) return 80;

//     // load configuration, optional: pass a file path.  Defaults are in settings.h
//     config_init("conf/z.conf");
//     std::cout << "\nCONFIG\n";
//     print_config(CFG());

//     mpfr_t t; mpfr_init2(t, CFG().mpfr_prec);
//     mpfr_set_str(t, "7900044.63199511046623912381512550604351774", 10, MPFR_RNDN);
//     verify_Z_at("7900044.63199511046623912381512550604351774");
//     mpfr_clear(t);

//     // in your verification harness
// mpfr_t Zthr, Zs, Zdir, diff1, diff2;
// mpfr_inits2(CFG().mpfr_prec, Zthr, Zs, Zdir, diff1, diff2, (mpfr_ptr)0);

// compute_Z_threaded(Zthr, t);     // now using compute_theta_exact internally
// compute_Z_small(Zs, t);          // also using compute_theta_exact internally
// rs_sum_direct(Zdir, t);          // direct reference

// mpfr_sub(diff1, Zthr, Zdir, MPFR_RNDN);
// mpfr_sub(diff2, Zs,   Zdir, MPFR_RNDN);
// print_val("Z_threaded",  Zthr, CFG().rs_log_digits);
// print_val("Z_small   ",  Zs,   CFG().rs_log_digits);
// print_val("Z_direct  ",  Zdir, CFG().rs_log_digits);
// print_val("thr - direct", diff1, CFG().rs_log_digits);
// print_val("small- direct", diff2, CFG().rs_log_digits);

// }