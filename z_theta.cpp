#include "settings.h"
#include <algorithm>
#include <mpfr.h>
#include <mpc.h>

// Coefficients c_m up to m=9 for:
// θ(t) = t/2 log(t/(2π)) - t/2 - π/8 + Σ c_m / t^{2m-1}
//
// Values are the standard ones derived from Bernoulli numbers.
// c1 = 1/48
// c2 = 7/5760
// c3 = 31/80640
// c4 = 127/430080
// c5 = 511/1216512
// c6 = 1414477/319334400
// c7 = 8191/3670016
// c8 = 118518239/135895449600
// c9 = 5749691557/1634992128000
//
// (You can trim/extend as you like.)
struct Rat
{
    long long num;
    long long den;
};

static const Rat C[] = {
    {0, 1}, // dummy for 0-index
    {1, 48},
    {7, 5760},
    {31, 80640},
    {127, 430080},
    {511, 1216512},
    {1414477, 319334400},
    {8191, 3670016},
    {118518239, 135895449600LL},
    {5749691557LL, 1634992128000LL},
};

void compute_theta_stirling_K(mpfr_t theta, const mpfr_t t, int K)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    if (K < 0)
        K = 0;
    if (K > 9)
        K = 9; // table cap; extend table if needed

    mpfr_t pi, two_pi, log_term, base, inv_t, pwr, term, tmp;
    mpfr_inits2(P, pi, two_pi, log_term, base, inv_t, pwr, term, tmp, (mpfr_ptr)0);

    // π and 2π
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    // base = t/(2π)
    mpfr_div(base, t, two_pi, MPFR_RNDN);

    // log_term = (t/2) * log( t/(2π) )
    mpfr_log(log_term, base, MPFR_RNDN);
    mpfr_mul(log_term, log_term, t, MPFR_RNDN);
    mpfr_div_ui(log_term, log_term, 2, MPFR_RNDN);

    // theta = log_term - t/2 - π/8
    mpfr_set(theta, log_term, MPFR_RNDN);
    mpfr_div_ui(tmp, t, 2, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);
    mpfr_div_ui(tmp, pi, 8, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);

    // series Σ c_m / t^{2m-1}
    mpfr_ui_div(inv_t, 1, t, MPFR_RNDN); // 1/t
    mpfr_set(pwr, inv_t, MPFR_RNDN);     // pwr = 1/t

    for (int m = 1; m <= K; ++m)
    {
        // Ensure pwr = 1 / t^{2m-1}
        if (m > 1)
        {
            // multiply by (1/t)^2 each step
            mpfr_mul(pwr, pwr, inv_t, MPFR_RNDN);
            mpfr_mul(pwr, pwr, inv_t, MPFR_RNDN);
        }
        mpfr_set_si(term, C[m].num, MPFR_RNDN);
        mpfr_div_si(term, term, C[m].den, MPFR_RNDN); // term = c_m (as mpfr)
        mpfr_mul(term, term, pwr, MPFR_RNDN);         // c_m / t^{2m-1}
        mpfr_add(theta, theta, term, MPFR_RNDN);
    }

    mpfr_clears(pi, two_pi, log_term, base, inv_t, pwr, term, tmp, (mpfr_ptr)0);
}

void compute_theta_stirling(mpfr_t theta, const mpfr_t t)
{
    compute_theta_stirling_K(theta, t, CFG().theta_stirling_terms);
}
