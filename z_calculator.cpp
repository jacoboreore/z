#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include "z_calculator_rs_gabcke.h"
#include "z_calculator_small.h"
#include "confirm_zero.h"
#include "refine_newton.h"
#include "settings.h"
#include "z_theta.h"
#include "z_ball.h"
#include "logs.h"
#include "utils.h"

// small helper for pretty MPFR prints (decimal, base=10)
static inline void print_val(const char *tag, const mpfr_t x, int digits = 24)
{
    std::cout << std::left << std::setw(12) << tag << " = ";
    mpfr_out_str(stdout, 10, digits, x, MPFR_RNDN);
    std::cout << "\n";
}

long rs_N(const mpfr_t t)
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

    if (N < 20)
    {
        compute_Z_small(Z, t);
    }
    else
    {
        compute_Z_rs_gabcke(Z, t);
    }

    // record best-so-far across the entire app
    best_update(Z, t);
}

// Z ball calculation for rigor certification
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

// ---------- quick theta for prefilter (3-term Stirling) ----------
static inline void theta_quick(mpfr_t theta, const mpfr_t t)
{
    const mpfr_prec_t P = mpfr_get_prec(t);
    mpfr_t pi, two_pi, u, half_t, term1, pi_over_8, invt, invt3, tmp;
    mpfr_inits2(P, pi, two_pi, u, half_t, term1, pi_over_8, invt, invt3, tmp, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    // u = t / (2π)
    mpfr_div(u, t, two_pi, MPFR_RNDN);

    // half_t = t/2
    mpfr_div_ui(half_t, t, 2, MPFR_RNDN);

    // term1 = (t/2)*log(t/(2π)) - (t/2)
    mpfr_log(term1, u, MPFR_RNDN);
    mpfr_mul(term1, term1, half_t, MPFR_RNDN);
    mpfr_sub(term1, term1, half_t, MPFR_RNDN);

    // theta = term1 - π/8
    mpfr_div_ui(pi_over_8, pi, 8, MPFR_RNDN);
    mpfr_sub(theta, term1, pi_over_8, MPFR_RNDN);

    // + 1/(48 t) + 7/(5760 t^3)
    mpfr_ui_div(invt, 1, t, MPFR_RNDN);
    mpfr_mul(invt3, invt, invt, MPFR_RNDN);
    mpfr_mul(invt3, invt3, invt, MPFR_RNDN);

    // tmp = 1/48 * 1/t
    mpfr_set_ui(tmp, 1, MPFR_RNDN);
    mpfr_div_ui(tmp, tmp, 48, MPFR_RNDN);
    mpfr_mul(tmp, tmp, invt, MPFR_RNDN);
    mpfr_add(theta, theta, tmp, MPFR_RNDN);

    // tmp = 7/5760 * 1/t^3
    mpfr_set_ui(tmp, 7, MPFR_RNDN);
    mpfr_div_ui(tmp, tmp, 5760, MPFR_RNDN);
    mpfr_mul(tmp, tmp, invt3, MPFR_RNDN);
    mpfr_add(theta, theta, tmp, MPFR_RNDN);

    mpfr_clears(pi, two_pi, u, half_t, term1, pi_over_8, invt, invt3, tmp, (mpfr_ptr)0);
}

// ---------- 3-term fast log1p for tiny x (|x| <= 1e-3) ----------
static inline void fast_log1p_3(mpfr_t out, const mpfr_t x)
{
    const double xd = mpfr_get_d(x, MPFR_RNDN);
    if (std::fabs(xd) <= 1e-3)
    {
        const mpfr_prec_t P = mpfr_get_prec(x);
        mpfr_t x2, x3, tmp;
        mpfr_inits2(P, x2, x3, tmp, (mpfr_ptr)0);
        mpfr_mul(x2, x, x, MPFR_RNDN);  // x^2
        mpfr_mul(x3, x2, x, MPFR_RNDN); // x^3
        // out = x - x^2/2 + x^3/3
        mpfr_div_ui(tmp, x2, 2, MPFR_RNDN); // x^2/2
        mpfr_neg(tmp, tmp, MPFR_RNDN);      // -x^2/2
        mpfr_add(out, x, tmp, MPFR_RNDN);   // x - x^2/2
        mpfr_div_ui(tmp, x3, 3, MPFR_RNDN); // x^3/3
        mpfr_add(out, out, tmp, MPFR_RNDN); // ... + x^3/3
        mpfr_clears(x2, x3, tmp, (mpfr_ptr)0);
    }
    else
    {
        mpfr_log1p(out, x, MPFR_RNDN);
    }
}

bool z_prefilter(const mpfr_t t)
{
    // pull knobs only from CFG()
    double alpha = CFG().prefilter_alpha;
    double tau = CFG().prefilter_tau;
    double tau_boost = CFG().prefilter_tau_boost;

    if (!(alpha > 0.0 && alpha < 1.0))
        alpha = 0.98;

    const long N = rs_N(t);

    // cheap small-N gate
    if (N < 8)
    {
        mpfr_t Zs;
        mpfr_init2(Zs, CFG().mpfr_prec);
        compute_Z_small(Zs, t);
        const double az = std::fabs(mpfr_get_d(Zs, MPFR_RNDN));

        LOGMP("N = %.*Rf |Z_small| = %.*Rg Keep = %s",
              N,
              az,
              ((az < 0.3) ? "YES" : "NO-ish"));

        mpfr_clear(Zs);
        if (az >= 0.3)
            return false;
        if (az <= 1e-2)
            return true;
        return true; // ambiguous -> keep
    }

    // K = floor(alpha * N)
    long K = (long)std::floor(alpha * (double)N);
    if (K < 4)
        K = std::min<long>(N, 4);
    if (K > N)
        K = N;

    // fast theta
    mpfr_t theta;
    mpfr_init2(theta, CFG().mpfr_prec);
    theta_quick(theta, t);

    // sum S_K
    mpfr_t S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd;
    mpfr_inits2(CFG().mpfr_prec, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, (mpfr_ptr)0);

    mpfr_set_zero(S, 0);
    mpfr_set_ui(n, 1, MPFR_RNDN);
    mpfr_set_zero(ln_n, 0);
    mpfr_sub(phi, theta, ln_n, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif
    mpfr_add(S, S, cphi, MPFR_RNDN); // n=1 term

    mpfr_t new_c, new_s, sqrt_i1, i1, tmp;
    mpfr_inits2(CFG().mpfr_prec, new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

    // internal tunables (can be moved to CFG() later if you want)
    const int BLOCK = 512;
    const int FASTLOG1P_MIN_I = 1000;
    const double R_SAFETY_MULTIPLIER = 1.05;

    // std::cerr << "[prefilter] N=" << N << "  K=" << K
    //           << "  alpha=" << alpha
    //           << "  tau=" << tau
    //           << "  tau_boost=" << tau_boost << "\n";

    auto emit_snapshot = [&](long i_or_K, const mpfr_t Sk)
    {
        // compute snapshot metrics
        mpfr_t twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, UB, TauEff;
        mpfr_inits2(CFG().mpfr_prec, twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, UB, TauEff, (mpfr_ptr)0);

        // |2S|
        mpfr_mul_ui(twoS, Sk, 2, MPFR_RNDN);
        mpfr_abs(abs_twoS, twoS, MPFR_RNDN);

        // 2*T_K with Ki = i_or_K
        mpfr_set_ui(sqrtN, (unsigned long)N, MPFR_RNDN);
        mpfr_sqrt(sqrtN, sqrtN, MPFR_RNDN);
        mpfr_set_ui(sqrtK, (unsigned long)i_or_K, MPFR_RNDN);
        mpfr_sqrt(sqrtK, sqrtK, MPFR_RNDN);
        mpfr_sub(tail, sqrtN, sqrtK, MPFR_RNDN);
        mpfr_mul_ui(tail, tail, 4, MPFR_RNDN); // 2*T

        remainder_bound(R, t);
        mpfr_mul_d(R, R, R_SAFETY_MULTIPLIER, MPFR_RNDN);

        // L = |2S| - (2T + R)
        mpfr_sub(L, abs_twoS, tail, MPFR_RNDN);
        mpfr_sub(L, L, R, MPFR_RNDN);

        // UB = |2S| + (2T + R)
        mpfr_add(UB, abs_twoS, tail, MPFR_RNDN);
        mpfr_add(UB, UB, R, MPFR_RNDN);

        const double tau_eff = tau + tau_boost / std::sqrt((double)N);
        mpfr_set_d(TauEff, tau_eff, MPFR_RNDN);

        // keep if L <= tau_eff
        const bool keep_here = (mpfr_cmp_d(L, tau_eff) <= 0);

        mpfr_clears(twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, UB, TauEff, (mpfr_ptr)0);
        return keep_here;
    };

    bool decided = false;
    bool keep = true; // default favoring recall

    for (long i = 1; i < K; ++i)
    {
        // Δφ = t * log1p(1/i)
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(q, 1, n, MPFR_RNDN);
        if (i >= FASTLOG1P_MIN_I)
            fast_log1p_3(dphi, q);
        else
            mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate (c,s) by -Δφ
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        // term for n = i+1
        mpfr_set_ui(i1, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_i1, i1, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_i1, MPFR_RNDN);
        mpfr_add(S, S, tmp, MPFR_RNDN);

        if ((i % BLOCK) == 0)
        {
            const bool keep_here = emit_snapshot(i + 1, S);
            if (!keep_here)
            {
                keep = false;
                decided = true;
                break;
            }
        }
    }

    if (!decided)
    {
        // final decision
        keep = emit_snapshot(K, S);
    }

    LOGMP("[prefilter] DECISION: %s", (keep ? "KEEP" : "DROP"));

    mpfr_clears(theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);
    return keep;
}
