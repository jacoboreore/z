#include <cmath>
#include "z_ball.h"
#include "settings.h"
#include "z_calculator.h"
#include "z_theta.h"
#include "z_calculator_small.h"
#include "z_calculator_rs_gabcke.h"
#include "utils.h"

static inline bool use_threaded_for_t(const mpfr_t t)
{
    return rs_N(t) >= 8;
}

static inline void rounding_radius(mpfr_t out, long N, double k = 4.0)
{
    mpfr_prec_t P = CFG().mpfr_prec;
    mpfr_t sN, ulp, kk;
    mpfr_inits2(P, sN, ulp, kk, (mpfr_ptr)0);

    mpfr_set_ui(sN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sN, sN, MPFR_RNDN);
    mpfr_mul_ui(sN, sN, 2, MPFR_RNDN);

    mpfr_set_ui(ulp, 1, MPFR_RNDN);
    mpfr_div_2si(ulp, ulp, (long)P - 1, MPFR_RNDN); // 2^{1-P}

    mpfr_set_d(kk, k, MPFR_RNDN);

    mpfr_mul(out, sN, ulp, MPFR_RNDN);
    mpfr_mul(out, out, kk, MPFR_RNDN);

    mpfr_clears(sN, ulp, kk, (mpfr_ptr)0);
}

// z_ball.cpp (replace your theta_radius)
static inline void theta_radius(mpfr_t out, const mpfr_t t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    const long N = rs_N(t);

    // Choose two truncations: K (from config) and K2 (a bit deeper)
    int K = CFG().theta_stirling_terms; // e.g., 8
    int K2 = K + 2;                     // one or two orders deeper
    if (K2 > 9)
        K2 = 9; // cap to table

    mpfr_t thK, thK2, dth, sN, factor;
    mpfr_inits2(P, thK, thK2, dth, sN, factor, (mpfr_ptr)0);

    compute_theta_stirling_K(thK, t, K);
    compute_theta_stirling_K(thK2, t, K2);

    mpfr_sub(dth, thK, thK2, MPFR_RNDN);
    mpfr_abs(dth, dth, MPFR_RNDN);

    // factor = 2*sqrt(N)
    mpfr_set_ui(sN, (unsigned long)N, MPFR_RNDN);
    mpfr_sqrt(sN, sN, MPFR_RNDN);
    mpfr_mul_ui(factor, sN, 2, MPFR_RNDN);

    // out = |Δθ| * 2√N
    mpfr_mul(out, dth, factor, MPFR_RNDN);

    mpfr_clears(thK, thK2, dth, sN, factor, (mpfr_ptr)0);
}

void compute_Z_ball_verbose(ZBallV &B, const mpfr_t t)
{
    // midpoint
    if (use_threaded_for_t(t))
    {
        mpfr_t tc;
        mpfr_init2(tc, CFG().mpfr_prec);
        mpfr_set(tc, t, MPFR_RNDN);
        compute_Z_rs_gabcke(B.mid, tc);
        mpfr_clear(tc);
    }
    else
    {
        compute_Z_small(B.mid, t);
    }

    // components
    remainder_bound(B.rad_R, t);
    rounding_radius(B.rad_round, rs_N(t), /*k=*/4.0);

    // If you have compute_theta_extended implemented, uncomment:
    // theta_radius(B.rad_theta, t);
    // else keep it at 0.

    // total radius
    mpfr_set(B.rad_tot, B.rad_R, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_round, MPFR_RNDN);
    mpfr_add(B.rad_tot, B.rad_tot, B.rad_theta, MPFR_RNDN);
}
