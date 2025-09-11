#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "z_calculator_thread.h"
#include "z_calculator_small.h"
#include "confirm_zero.h"
#include "refine_newton.h"
#include "settings.h"
#include "z_theta.h"
#include "z_ball.h"
#include "logs.h"


// small helper for pretty MPFR prints (decimal, base=10)
static inline void print_val(const char *tag, const mpfr_t x, int digits = 24)
{
    std::cout << std::left << std::setw(12) << tag << " = ";
    mpfr_out_str(stdout, 10, digits, x, MPFR_RNDN);
    std::cout << "\n";
}

static long rs_N(const mpfr_t t)
{
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(CFG().mpfr_prec, pi, two_pi, tmp, (mpfr_ptr)0);
    
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);

    if (N < 1)
        N = 1;

    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);

    return N;
}

// ---------------- Unified compute_Z ----------------
void compute_Z(mpfr_t Z, mpfr_t t)
{
    long N = rs_N(t);

    if (N < 8)
    {
        compute_Z_small(Z, t);
    }
    else
    {
        // printf("USING THREADING...\n");  // keep or gate by a verbose flag
        compute_Z_threaded(Z, t);
    }
}


void compute_Z_ball(mpfr_t Z, mpfr_t t)
{
    long N = rs_N(t);

    mpfr_t Z1, t_star, Z0, abs_mid, rad_total, t_cur;
    mpfr_inits2(CFG().mpfr_prec, Z1, t_star, Z0, abs_mid, rad_total, t_cur, (mpfr_ptr)0);

    ZBallV b;
    zballv_init(b, CFG().mpfr_prec);

    mpfr_printf("t = %.30Rf\n", t);
    std::cout << "Riemann-Siegel N = " << N << "\n";
    std::cout << std::scientific << "est_spacing(t) = " << est_spacing(t)
              << "  (≈ pi / log(t/2pi))\n";

    compute_Z(Z0, t);
    print_val("Z(t)", Z0, 50);

    compute_Z_ball_verbose(b, t);

    mpfr_abs(abs_mid, b.mid, MPFR_RNDN);
    mpfr_set(rad_total, b.rad_R, MPFR_RNDN);
    mpfr_add(rad_total, rad_total, b.rad_round, MPFR_RNDN);
    mpfr_add(rad_total, rad_total, b.rad_theta, MPFR_RNDN);

    std::cout << "\n[ball @ t]\n";
    print_val("mid", b.mid, 24);
    print_val("|mid|", abs_mid, 24);
    print_val("rad_R", b.rad_R, 24);
    print_val("rad_round", b.rad_round, 24);
    print_val("rad_theta", b.rad_theta, 24);
    print_val("rad_total", rad_total, 24);

    bool contains = (mpfr_cmp(abs_mid, rad_total) <= 0);
    std::cout << "\ncontains_zero (ball test): " << (contains ? "YES" : "NO") << "\n";

    mpfr_set(t_cur, t, MPFR_RNDN);

    if (!contains)
    {
        bool ok = newton_refine_root(
            t_star, t,
            CFG().newton_tol_abs,
            CFG().newton_tol_step,
            CFG().newton_max_iter,
            CFG().use_threaded_z);

        std::cout << "\nNewton " << (ok ? "converged" : "stopped (fail)") << "\n";
        mpfr_printf("t* = %.30Rf\n", t_star);

        compute_Z(Z1, t_star);
        print_val("Z(t*)", Z1, 24);

        // reuse the same ball object safely by clear+init
        zballv_clear(b);
        zballv_init(b, CFG().mpfr_prec);
        compute_Z_ball_verbose(b, t_star);

        mpfr_abs(abs_mid, b.mid, MPFR_RNDN);
        mpfr_set(rad_total, b.rad_R, MPFR_RNDN);
        mpfr_add(rad_total, rad_total, b.rad_round, MPFR_RNDN);
        mpfr_add(rad_total, rad_total, b.rad_theta, MPFR_RNDN);

        std::cout << "\n[ball @ t*]\n";
        print_val("|mid|", abs_mid, 24);
        print_val("rad_total", rad_total, 24);
        bool contains2 = (mpfr_cmp(abs_mid, rad_total) <= 0);
        std::cout << "contains_zero (ball test): " << (contains2 ? "YES" : "NO") << "\n";

        contains = contains2;
        mpfr_set(t_cur, t_star, MPFR_RNDN);
    }
    else
    {
        std::cout << "\nAlready certified at input t.\n";
    }

    bool signflip = confirm_zero_like(t_cur, /*use_threaded_z=*/1);
    std::cout << "confirm_zero_like(±0.2·spacing): " << (signflip ? "YES" : "NO") << "\n";

    bool accept = contains && signflip;
    std::cout << "\nFINAL DECISION: " << (accept ? "ACCEPT (likely zero nearby)" : "REJECT") << "\n";

    mpfr_clears(Z1, t_star, Z0, abs_mid, rad_total, t_cur, (mpfr_ptr)0);
    zballv_clear(b);
}


bool z_prefilter(const mpfr_t t, double alpha, double tau)
{
    // Clamp alpha to sensible bounds
    if (alpha <= 0.0) alpha = 0.25;
    if (alpha >= 1.0) alpha = 0.8;

    // --- Compute N first with a tiny, self-contained MPFR scope (no leaks on early return)
    long N;
    {
        mpfr_t pi, two_pi, tmp;
        mpfr_inits2(CFG().mpfr_prec, pi, two_pi, tmp, (mpfr_ptr)0);
        mpfr_const_pi(pi, MPFR_RNDN);
        mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
        mpfr_div(tmp, t, two_pi, MPFR_RNDN);
        mpfr_sqrt(tmp, tmp, MPFR_RNDN);
        N = mpfr_get_si(tmp, MPFR_RNDN);
        if (N < 1) N = 1;
        mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);
    }

    // --- Small-N fast gate (no large allocations before we might return)
    if (N < 8)
    {
        mpfr_t Zs;
        mpfr_init2(Zs, CFG().mpfr_prec);
        compute_Z_small(Zs, t);
        double az = std::fabs(mpfr_get_d(Zs, MPFR_RNDN));
        mpfr_clear(Zs);

        if (az <= 1e-2)  return true;   // near zero → keep
        if (az >= 0.3)   return false;  // clearly far → drop
        return true;                     // ambiguous? favor recall
    }

    // --- Choose K = floor(alpha*N), bounded to [min(4,N), N]
    long K = (long)std::floor(alpha * (double)N);
    if (K < 4) K = (N >= 4 ? 4 : N);
    if (K > N) K = N;

    // --- Now allocate all MPFR vars we need for the main path
    mpfr_t theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd;
    mpfr_inits2(CFG().mpfr_prec, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, (mpfr_ptr)0);

    // theta via Stirling wrapper (uses CFG().theta_stirling_terms)
    compute_theta_stirling(theta, t);

    // Start at n = 1
    mpfr_set_zero(S, 0);
    mpfr_set_ui(n, 1, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, t, ln_n, MPFR_RNDN);
    mpfr_sub(phi, theta, phi, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif
    mpfr_add(S, S, cphi, MPFR_RNDN); // term at n=1

    // Rotation accumulators
    mpfr_t new_c, new_s, sqrt_i1, i1, tmp;
    mpfr_inits2(CFG().mpfr_prec, new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

    for (long i = 1; i < K; ++i)
    {
        // Δφ = t * log1p(1/i)
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(q, 1, n, MPFR_RNDN);        // q = 1/i
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // Rotate: (cphi, sphi) ← (cphi, sphi) * R(−Δφ)
        // new_c = c*cd + s*sd
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        // new_s = s*cd − c*sd
        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        // Add term (n = i+1): cos(phi) / sqrt(i+1)
        mpfr_set_ui(i1, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_i1, i1, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_i1, MPFR_RNDN);
        mpfr_add(S, S, tmp, MPFR_RNDN);
    }

    mpfr_clears(new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

    // --- Bounds and decision
    mpfr_t twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, TauEff;
    mpfr_inits2(CFG().mpfr_prec, twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, TauEff, (mpfr_ptr)0);

    // |2 S_K|
    mpfr_mul_ui(twoS, S, 2, MPFR_RNDN);
    mpfr_abs(abs_twoS, twoS, MPFR_RNDN);

    // 2*T_K = 4 (√N − √K)
    mpfr_set_ui(sqrtN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sqrtN, sqrtN, MPFR_RNDN);
    mpfr_set_ui(sqrtK, (unsigned long)K, MPFR_RNDN);
    mpfr_sqrt(sqrtK, sqrtK, MPFR_RNDN);
    mpfr_sub(tail, sqrtN, sqrtK, MPFR_RNDN);
    mpfr_mul_ui(tail, tail, 4, MPFR_RNDN); // tail now holds 2*T_K for Z

    // RS remainder bound
    remainder_bound(R, t);

    // L = |2 S_K| - 2*T_K - R
    mpfr_sub(L, abs_twoS, tail, MPFR_RNDN);
    mpfr_sub(L, L, R, MPFR_RNDN);

    // τ_eff = τ + c / √N  (slightly relaxed for small N)
    double tau_boost = 0.50; // tune: 0.4–0.7 works well
    double tau_eff = tau + tau_boost / std::sqrt((double)N);
    mpfr_set_d(TauEff, tau_eff, MPFR_RNDN);

    bool keep = (mpfr_cmp(L, TauEff) <= 0);

    // --- Logs
    std::cerr << "[prefilter] N=" << N << " K=" << K
              << " alpha=" << alpha << " tau=" << tau << "\n";
    log_mpfr("[prefilter] |2S_K|", abs_twoS, 18);
    log_mpfr("[prefilter] 2*T_K", tail, 18);
    log_mpfr("[prefilter] R", R, 18);
    log_mpfr("[prefilter] L", L, 18);
    std::cerr << "[prefilter] keep=" << (keep ? "YES" : "NO") << "\n";

    // --- Cleanup
    mpfr_clears(theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, (mpfr_ptr)0);
    mpfr_clears(twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, TauEff, (mpfr_ptr)0);

    return keep;
}
