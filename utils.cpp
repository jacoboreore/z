#include "utils.h"
#include "logs.h"
#include "settings.h"
#include <cassert>
#include <limits>
#include <cmath>

static std::mutex &best_mutex()
{
    static std::mutex m;
    return m;
}

BestZeroState &best_state()
{
    static BestZeroState s;
    return s;
}

void best_init(mpfr_prec_t P)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());

    if (s.initialized && s.P == P)
        return;

    if (s.initialized)
    {
        mpfr_clear(s.best_absZ);
        mpfr_clear(s.best_t);
        s.initialized = false;
        s.P = 0;
    }

    mpfr_init2(s.best_absZ, P);
    mpfr_init2(s.best_t, P);
    s.P = P;
    s.initialized = true;

    // sentinel values
    mpfr_set_d(s.best_absZ, 9'999'999.0, MPFR_RNDN); // start huge so any real |Z| wins
    mpfr_set_zero(s.best_t, 0);                      // 0 with sign=+0
}

void best_reset()
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(P) must be called before best_reset()");
    mpfr_set_d(s.best_absZ, 9'999'999.0, MPFR_RNDN);
    mpfr_set_zero(s.best_t, 0);
}

void best_update(mpfr_srcptr Z, mpfr_srcptr t)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(P) must be called before best_update()");

    mpfr_t absZ;
    mpfr_init2(absZ, s.P);

    mpfr_abs(absZ, Z, MPFR_RNDN);
    if (mpfr_cmp(absZ, s.best_absZ) < 0)
    {
        LOGMP("(*) New better has been calculated t = %.*Rf  |Z| = %.*Rg",
              CFG().rs_log_digits, t,
              CFG().rs_log_digits, Z);

        mpfr_set(s.best_absZ, absZ, MPFR_RNDN);
        mpfr_set(s.best_t, t, MPFR_RNDN);
    }

    mpfr_clear(absZ);
}

void best_get(mpfr_t out_t, mpfr_t out_absZ)
{
    auto &s = best_state();
    std::lock_guard<std::mutex> lock(best_mutex());
    assert(s.initialized && "best_init(Precision) must be called before best_get()");

    // Caller must have set appropriate precision on out_t/out_absZ
    mpfr_set(out_t, s.best_t, MPFR_RNDN);
    mpfr_set(out_absZ, s.best_absZ, MPFR_RNDN);
}

void get_best_t_so_far()
{
    mpfr_t best_t, best_absZ;
    mpfr_inits2(CFG().mpfr_prec, best_t, best_absZ, (mpfr_ptr)0);
    best_get(best_t, best_absZ);

    LOGMP("Best calculations t = %.*Rf |Z| = %.*Rg",
          CFG().rs_log_digits,
          best_t,
          CFG().rs_log_digits,
          best_absZ);

    mpfr_clears(best_t, best_absZ, (mpfr_ptr)0);
}

// Simple diagnostic tail bound you already had
void remainder_bound(mpfr_t R, const mpfr_t t)
{
    const double C = CFG().rs_tail_C;
    mpfr_t r, tmp;
    mpfr_inits2(CFG().mpfr_prec, r, tmp, (mpfr_ptr)0);
    mpfr_sqrt(tmp, t, MPFR_RNDN); // √t
    mpfr_sqrt(r, tmp, MPFR_RNDN); // t^(1/4)
    mpfr_ui_div(R, 1, r, MPFR_RNDN);
    mpfr_mul_d(R, R, C, MPFR_RNDN);
    mpfr_clears(r, tmp, (mpfr_ptr)0);
}

long double mpfr_to_long_double_safe(const mpfr_t x,
                                            mpfr_rnd_t rnd = MPFR_RNDN) {
    if (!mpfr_number_p(x)) {
        if (mpfr_nan_p(x))   return std::numeric_limits<long double>::quiet_NaN();
        if (mpfr_inf_p(x) > 0) return std::numeric_limits<long double>::infinity();
        if (mpfr_inf_p(x) < 0) return -std::numeric_limits<long double>::infinity();
    }

    long double v = mpfr_get_ld(x, rnd);

    // Optional: clamp if it overflowed the long double range
    if (!std::isfinite(v)) {
        return (mpfr_sgn(x) >= 0)
            ? std::numeric_limits<long double>::infinity()
            : -std::numeric_limits<long double>::infinity();
    }
    return v;
}