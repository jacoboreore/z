// #include <mpfr.h>
// #include <cmath>
// #include <cstdio>
// #include <vector>
// #include <string>
// #include <utility>
// //#include "z_walk.h"
// #include "config.h"
// #include "z_calculator.h"
// #include "refine_newton.h"

// // int main(int argc, char* argv[]) {
// //     mpfr_set_default_prec(MPFR_PRECISION);

// //     // Iterate through and print each argument
// //     std::cout << "Arguments:" << std::endl;
// //     for (int i = 0; i < argc; ++i) {
// //         std::cout << "  argv[" << i << "]: " << argv[i] << std::endl;
// //     }

// //     mpfr_t t, Z;
// //     mpfr_inits2(MPFR_PRECISION, t, Z, (mpfr_ptr)0);

// //     // First zero (sanity check)
// //     mpfr_set_str(t, argv[1], 10, MPFR_RNDN);
// //     mpfr_printf("%.10RF\n", t); // %.10RF prints the number with 10 digits after the decimal point
// //     compute_Z(Z, t);
// //     mpfr_printf("Z(%.15Rf) = %.20Re\n", t, Z);

// //     // Your big y sample
// //     //mpfr_set_str(t, "4400961943970.5696075625388437339818624451523392290921263943249583984239027506119881", 50, MPFR_RNDN);
// //     //compute_Z(Z, t);
// //     //mpfr_printf("Z(%.15Rf) = %.10Re\n", t, Z);

// //     mpfr_clears(t, Z, (mpfr_ptr)0);
// //     return 0;
// // }

// // helper: mpfr -> double
// static inline double dval(const mpfr_t x) { return mpfr_get_d(x, MPFR_RNDN); }

// struct Case { const char* label; const char* t_str; };

// static void run_suite(double alpha, double tau,
//                       double near_eps = 1e-3, double far_min = 0.5)
// {
//     printf("=== z_prefilter tests  alpha=%.3f  tau=%.3f  (near_eps=%.1e, far_min=%.2f)\n",
//            alpha, tau, near_eps, far_min);

//     // Mix of near-zeros and clearly-off points
//     // First few nontrivial zeros (approx):
//     // 14.1347251417, 21.0220396388, 25.0108575801, 30.4248761259, 32.9350615877, 37.5861781588
//     std::vector<Case> cases = {
//         {"zero-ish", "14.134725141734693"},
//         {"off",      "14.20"},
//         {"zero-ish", "21.022039638771555"},
//         {"off",      "21.10"},
//         {"zero-ish", "25.010857580145689"},
//         {"off",      "25.20"},
//         {"zero-ish", "30.424876125859513"},
//         {"off",      "30.60"},
//         {"zero-ish", "32.935061587739190"},
//         {"off",      "33.10"},
//         {"zero-ish", "37.586178158825671"},
//         {"off",      "37.80"},
//         // some mid t where RS sum is nontrivial
//         {"mid",      "400"},
//         {"mid",      "500"},
//         {"mid",      "800"},
//         {"mid",      "1000"}
//     };

//     int tp=0, fn=0, tn=0, fp=0, amb=0;

//     mpfr_t t, Z; mpfr_inits2(256, t, Z, (mpfr_ptr)0);

//     for (auto &c : cases) {
//         mpfr_set_str(t, c.t_str, 10, MPFR_RNDN);

//         // Prefilter decision
//         bool keep = z_prefilter(t, alpha, tau);

//         // Ground truth via small-Z
//         compute_Z(Z, t);
//         double z  = dval(Z);
//         double az = std::fabs(z);

//         // classify: near / far / ambiguous (for evaluation only)
//         enum { NEAR, FAR, AMB } clazz = AMB;
//         if (az <= near_eps) clazz = NEAR;
//         else if (az >= far_min) clazz = FAR;
//         else clazz = AMB;

//         // update counters
//         if (clazz == NEAR) { if (keep) ++tp; else ++fn; }
//         else if (clazz == FAR) { if (keep) ++fp; else ++tn; }
//         else { ++amb; }

//         printf("t=%-18s  |Z|=%-10.6g  prefilter=%-5s  class=%s\n",
//                c.t_str, az, keep ? "keep" : "drop",
//                clazz==NEAR ? "near" : (clazz==FAR ? "far" : "amb"));
//     }

//     mpfr_clears(t, Z, (mpfr_ptr)0);

//     printf("\nConfusion (excluding 'amb'):\n");
//     printf("  TP(near & keep)   = %d\n", tp);
//     printf("  FN(near & drop)   = %d\n", fn);
//     printf("  TN(far  & drop)   = %d\n", tn);
//     printf("  FP(far  & keep)   = %d\n", fp);

//     int tot = tp+tn+fp+fn;
//     if (tot > 0) {
//         double recall = (tp+fn) ? double(tp)/double(tp+fn) : 1.0;
//         double precision = (tp+fp) ? double(tp)/double(tp+fp) : 1.0;
//         printf("  recall (zeros caught)   = %.3f\n", recall);
//         printf("  precision (keeps legit) = %.3f\n", precision);
//     }
//     printf("  amb (ignored in stats)  = %d\n", amb);
//     printf("===============================================================\n\n");
// }

// int main() {
//     // // deterministic precision
//     // mpfr_set_default_prec(256);

//     // // Try a couple of (alpha, tau) settings
//     // run_suite(/*alpha*/0.50, /*tau*/0.30);
//     // run_suite(/*alpha*/0.65, /*tau*/0.30);
//     // run_suite(/*alpha*/0.50, /*tau*/0.20);


//       // MPFR
//   printf("MPFR runtime : %s\n", mpfr_get_version());
// #ifdef MPFR_VERSION_STRING
//   printf("MPFR headers : %s\n", MPFR_VERSION_STRING);
// #else
//   printf("MPFR headers : (no MPFR_VERSION_STRING macro)\n");
// #endif

//   // GMP
//   printf("GMP  runtime : %s\n", gmp_version);
//   printf("GMP  headers : %d.%d.%d\n",
//          __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);

//   // MPC
// #ifdef mpc_get_version
//   printf("MPC  runtime : %s\n", mpc_get_version());
// #else
//   printf("MPC  runtime : (no mpc_get_version symbol)\n");
// #endif
// #ifdef MPC_VERSION_STRING
//   printf("MPC  headers : %s\n", MPC_VERSION_STRING);
// #else
//   printf("MPC  headers : (no MPC_VERSION_STRING macro)\n");
// #endif

  
//     mpfr_t t, Z; 
//     mpfr_inits2(256, t, Z, (mpfr_ptr)0);
//     mpfr_t t0, t_star; 
//     mpfr_inits2(MPFR_PRECISION, t0, t_star, (mpfr_ptr)0);

//     // Start near the 1st nontrivial zero
//     mpfr_set_str(t0, "300000035.0600766", 10, MPFR_RNDN);

//     bool ok = newton_refine_root(
//         t_star, 
//         t0,
//         /*tol_abs=*/1e-20,
//         /*tol_step=*/1e-14,
//         /*max_iter=*/20,
//         /*use_threaded_z=*/1
//     );

//     printf("Newton %s\n", ok ? "converged" : "stopped");
//     printf("t = "); 
//     mpfr_printf("%.30Rf\n", t0); 
//     printf("\n");
//     mpfr_printf("%.30Rf\n", t_star); 
//     printf("\n");

//     compute_Z(Z, t_star);

//     mpfr_printf("%.30Rf\n", Z); 
//     printf("\n");

//     compute_Z(Z, t0);

//     mpfr_printf("Z original %.30Rf\n", Z); 
//     printf("\n");

//     mpfr_clears(t0, t_star, (mpfr_ptr)0);

//     return 0;
// }

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <gmp.h>
#include <mpfr.h>

#include "config.h"
#include "z_calculator.h"   // compute_Z(...), rs_N_for_log(...), est_spacing(...)
#include "z_ball.h"         // ZBallV, zballv_init/clear, compute_Z_ball_verbose, zballv_contains_zero
#include "refine_newton.h"  // newton_refine_root(...)
#include "confirm_zero.h"   // confirm_zero_like(...)

// small helper for pretty MPFR prints (decimal, base=10)
static inline void print_val(const char* tag, const mpfr_t x, int digits = 24) {
    std::cout << std::left << std::setw(12) << tag << " = ";
    mpfr_out_str(stdout, 10, digits, x, MPFR_RNDN);
    std::cout << "\n";
}

static long rs_N_for_log(const mpfr_t t) {
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
    return (N < 1 ? 1 : N);
}

int main(int argc, char** argv) {
    mpfr_set_default_prec(MPFR_PRECISION);

    // ---- version banner -----------------------------------------------------
    std::cout << "MPFR runtime : " << mpfr_get_version() << "\n";
#ifdef MPFR_VERSION_STRING
    std::cout << "MPFR headers : " << MPFR_VERSION_STRING << "\n";
#else
    std::cout << "MPFR headers : (no MPFR_VERSION_STRING)\n";
#endif
    std::cout << "GMP  runtime : " << gmp_version << "\n";
    std::cout << "GMP  headers : "
              << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << "\n";
#ifdef mpc_get_version
    std::cout << "MPC  runtime : " << mpc_get_version() << "\n";
#else
    std::cout << "MPC  runtime : (no mpc_get_version)\n";
#endif
#ifdef MPC_VERSION_STRING
    std::cout << "MPC  headers : " << MPC_VERSION_STRING << "\n";
#else
    std::cout << "MPC  headers : (no MPC_VERSION_STRING)\n";
#endif
    std::cout << "------------------------------------------------------------------\n";

    // ---- pick t -------------------------------------------------------------
    mpfr_t t; mpfr_init2(t, MPFR_PRECISION);
    if (argc > 1) {
        if (mpfr_set_str(t, argv[1], 10, MPFR_RNDN) != 0) {
            std::cerr << "Could not parse t from argv[1]; using default.\n";
            mpfr_set_str(t, "300000033.1200961947076831218132392973205613", 10, MPFR_RNDN);
        }
    } else {
        mpfr_set_str(t, "300000033.1200961947076831218132392973205613", 10, MPFR_RNDN);
    }

    // basic info
    mpfr_printf("t = %.30Rf\n", t);
    long N = rs_N_for_log(t);
    std::cout << "Riemann-Siegel N = " << N << "\n";
    double spacing = est_spacing(t);
    std::cout << std::scientific << "est_spacing(t) ~ " << spacing << "  (≈ pi / log(t/2pi))\n";

    // baseline Z(t)
    mpfr_t Z0; mpfr_init2(Z0, MPFR_PRECISION);
    compute_Z(Z0, t);
    print_val("Z(t)", Z0, 24);

    // ---- build Z-ball at t --------------------------------------------------
    ZBallV b; zballv_init(b, MPFR_PRECISION);
    compute_Z_ball_verbose(b, t);

    mpfr_t abs_mid; mpfr_init2(abs_mid, MPFR_PRECISION);
    mpfr_abs(abs_mid, b.mid, MPFR_RNDN);

    // compute total radius = rad_R + rad_round + rad_theta
    mpfr_t rad_total; mpfr_init2(rad_total, MPFR_PRECISION);
    mpfr_set(rad_total, b.rad_R, MPFR_RNDN);
    mpfr_add(rad_total, rad_total, b.rad_round, MPFR_RNDN);
    mpfr_add(rad_total, rad_total, b.rad_theta, MPFR_RNDN);

    std::cout << "\n[ball @ t]\n";
    print_val("mid",        b.mid,      24);
    print_val("|mid|",      abs_mid,    24);
    print_val("rad_R",      b.rad_R,    24);
    print_val("rad_round",  b.rad_round,24);
    print_val("rad_theta",  b.rad_theta,24);
    print_val("rad_total",  rad_total,  24);

    bool contains = (mpfr_cmp(abs_mid, rad_total) <= 0);
    std::cout << "\ncontains_zero (ball test): " << (contains ? "YES" : "NO") << "\n";

    // ---- sign-flip sanity (±0.2 * spacing window) --------------------------
    bool signflip = confirm_zero_like(t, /*use_threaded_z=*/1);
    std::cout << "confirm_zero_like(±0.2·spacing): " << (signflip ? "YES" : "NO") << "\n";

    // ---- optionally refine with Newton if not yet certified -----------------
    mpfr_t t_star; mpfr_init2(t_star, MPFR_PRECISION);
    mpfr_set(t_star, t, MPFR_RNDN);

    if (!contains) {
        bool ok = newton_refine_root(
            t_star, t,
            /*tol_abs=*/1e-5,     // sensible with ~50–80 bit precision
            /*tol_step=*/1e-12,
            /*max_iter=*/20,
            /*use_threaded_z=*/1
        );
        std::cout << "\nNewton " << (ok ? "converged" : "stopped") << "\n";
        mpfr_printf("t* = %.30Rf\n", t_star);

        mpfr_t Z1; mpfr_init2(Z1, MPFR_PRECISION);
        compute_Z(Z1, t_star);
        print_val("Z(t*)", Z1, 24);

        // rebuild ball at t*
        zballv_clear(b); zballv_init(b, MPFR_PRECISION);
        compute_Z_ball_verbose(b, t_star);

        mpfr_abs(abs_mid, b.mid, MPFR_RNDN);
        mpfr_set(rad_total, b.rad_R, MPFR_RNDN);
        mpfr_add(rad_total, rad_total, b.rad_round, MPFR_RNDN);
        mpfr_add(rad_total, rad_total, b.rad_theta, MPFR_RNDN);

        std::cout << "\n[ball @ t*]\n";
        print_val("|mid|",     abs_mid,   24);
        print_val("rad_total", rad_total, 24);
        bool contains2 = (mpfr_cmp(abs_mid, rad_total) <= 0);
        std::cout << "contains_zero (ball test): " << (contains2 ? "YES" : "NO") << "\n";

        mpfr_clear(Z1);
        contains = contains2;   // carry forward for final decision
        mpfr_set(t, t_star, MPFR_RNDN);
    } else {
        std::cout << "\nAlready certified at input t.\n";
    }

    // ---- final decision: ball + signflip -----------------------------------
    bool accept = contains && signflip;
    std::cout << "\nFINAL DECISION: " << (accept ? "ACCEPT (likely zero nearby)" : "REJECT") << "\n";

    // tidy
    mpfr_clears(t, t_star, Z0, abs_mid, rad_total, (mpfr_ptr)0);
    zballv_clear(b);
    return 0;
}
