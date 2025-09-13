// utils.cpp
#include "utils.h"
#include "logs.h"
#include <cassert>

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
        LOGMP("New better has been calculated t=%.*Rf  |Z|=%.*Rg",
               24, t,
               20, Z);
               
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