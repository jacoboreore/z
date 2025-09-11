#pragma once
#include <mpfr.h>

struct ZBallV
{
    mpfr_t mid;       // midpoint
    mpfr_t rad_tot;   // total radius (used for contains-zero test)
    mpfr_t rad_R;     // analytic RS remainder bound
    mpfr_t rad_round; // rounding accumulation bound
    mpfr_t rad_theta; // phase (theta) truncation bound (optional; see below)
};

static inline void zballv_init(ZBallV &b, mpfr_prec_t P)
{
    mpfr_init2(b.mid, P);
    mpfr_init2(b.rad_tot, P);
    mpfr_init2(b.rad_R, P);
    mpfr_init2(b.rad_round, P);
    mpfr_init2(b.rad_theta, P);
    mpfr_set_zero(b.mid, 0);
    mpfr_set_zero(b.rad_tot, 0);
    mpfr_set_zero(b.rad_R, 0);
    mpfr_set_zero(b.rad_round, 0);
    mpfr_set_zero(b.rad_theta, 0);
}

static inline void zballv_clear(ZBallV &b)
{
    mpfr_clear(b.mid);
    mpfr_clear(b.rad_tot);
    mpfr_clear(b.rad_R);
    mpfr_clear(b.rad_round);
    mpfr_clear(b.rad_theta);
}

static inline bool zballv_contains_zero(const ZBallV &b)
{
    mpfr_t a;
    mpfr_init2(a, mpfr_get_prec(b.mid));
    mpfr_abs(a, b.mid, MPFR_RNDN);

    bool ok = (mpfr_cmp(a, b.rad_tot) <= 0);

    mpfr_clear(a);

    return ok;
}

void compute_Z_ball_verbose(ZBallV &B, const mpfr_t t);