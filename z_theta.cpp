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

// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
// Lanczos with g=7, N=9 (classic coefficients). Accurate for small/medium t.
// void compute_theta_lanczos(mpfr_t theta, const mpfr_t t)
// {
//     const mpfr_prec_t P = CFG().mpfr_prec;

//     // ---- constants ----
//     mpfr_t half, quarter, logpi, ln2pi_half;
//     mpfr_inits2(P, half, quarter, logpi, ln2pi_half, (mpfr_ptr)0);
//     mpfr_set_d(half, 0.5, MPFR_RNDN);
//     mpfr_set_d(quarter, 0.25, MPFR_RNDN);
//     mpfr_log(logpi, mpfr_const_pi(), MPFR_RNDN); // needs MPFR >= 4.1; fallback below if not

//     // ln(sqrt(2π)) = 0.5*ln(2π)
//     mpfr_t two_pi, ln_two_pi;
//     mpfr_inits2(P, two_pi, ln_two_pi, (mpfr_ptr)0);
//     mpfr_const_pi(two_pi, MPFR_RNDN);
//     mpfr_mul_ui(two_pi, two_pi, 2, MPFR_RNDN);
//     mpfr_log(ln_two_pi, two_pi, MPFR_RNDN);
//     mpfr_div_ui(ln2pi_half, ln_two_pi, 2, MPFR_RNDN);

//     // ---- Lanczos coefficients for g = 7, N = 9 ----
//     // p[0..8] (real)
//     static const char *PSTR[9] = {
//         "0.99999999999980993227684700473478",
//         "676.520368121885098567009190444019",
//         "-1259.13921672240287047156078755283",
//         "771.3234287776530788486528258894",
//         "-176.61502916214059906584551354",
//         "12.507343278686904814458936853",
//         "-0.13857109526572011689554707",
//         "9.984369578019570859563e-6",
//         "1.50563273514931155834e-7"};
//     mpfr_t p[9];
//     for (int i = 0; i < 9; ++i)
//     {
//         mpfr_init2(p[i], P);
//         mpfr_set_str(p[i], PSTR[i], 10, MPFR_RNDN);
//     }

//     // z = 1/4 + i*(t/2)
//     mpc_t z, sum, term, denom, x, z_minus_half, t1, t2, t3, lg;
//     mpc_init2(z, P);
//     mpc_init2(sum, P);
//     mpc_init2(term, P);
//     mpc_init2(denom, P);
//     mpc_init2(x, P);
//     mpc_init2(z_minus_half, P);
//     mpc_init2(t1, P);
//     mpc_init2(t2, P);
//     mpc_init2(t3, P);
//     mpc_init2(lg, P);

//     mpfr_t t_half;
//     mpfr_init2(t_half, P);
//     mpfr_div_ui(t_half, t, 2, MPFR_RNDN);
//     mpc_set_fr_fr(z, quarter, t_half, MPC_RNDNN);

//     // A_g(z) = p0 + sum_{k=1..8} p_k / (z + k)
//     mpc_set_fr_fr(sum, p[0], mpfr_const_zero(), MPC_RNDNN); // real p0
//     for (int k = 1; k <= 8; ++k)
//     {
//         // denom = z + k
//         mpc_set(denom, z, MPC_RNDNN);
//         mpfr_t kfr;
//         mpfr_init2(kfr, P);
//         mpfr_set_ui(kfr, (unsigned long)k, MPFR_RNDN);
//         mpc_add_fr(denom, denom, kfr, MPC_RNDNN);
//         // term = p_k / denom
//         mpc_set_fr_fr(term, p[k], mpfr_const_zero(), MPC_RNDNN);
//         mpc_div(term, term, denom, MPC_RNDNN);
//         mpc_add(sum, sum, term, MPC_RNDNN);
//         mpfr_clear(kfr);
//     }

//     // x = z + (g - 1/2), with g=7 -> (g-1/2)=6.5
//     mpfr_t gmh;
//     mpfr_init2(gmh, P);
//     mpfr_set_d(gmh, 6.5, MPFR_RNDN);
//     mpc_set(x, z, MPC_RNDNN);
//     mpc_add_fr(x, x, gmh, MPC_RNDNN);

//     // z_minus_half = z - 1/2
//     mpc_set(z_minus_half, z, MPC_RNDNN);
//     mpc_sub_fr(z_minus_half, z_minus_half, half, MPC_RNDNN);

//     // lg = (z-1/2)*log(x) - x + log(sqrt(2π)) + log(sum)
//     mpc_log(t1, x, MPC_RNDNN); // t1 = log(x)
//     mpc_mul(t1, z_minus_half, t1, MPC_RNDNN);
//     mpc_neg(t2, x, MPC_RNDNN);   // t2 = -x
//     mpc_log(t3, sum, MPC_RNDNN); // t3 = log(sum)
//     mpc_add(lg, t1, t2, MPC_RNDNN);
//     mpc_add_fr(lg, lg, ln2pi_half, MPC_RNDNN);
//     mpc_add(lg, lg, t3, MPC_RNDNN);

//     // theta = Im(lg) - (t/2) * log(pi)
//     mpfr_t im;
//     mpfr_init2(im, P);
//     mpc_imag(im, lg, MPFR_RNDN);
//     mpfr_t corr;
//     mpfr_init2(corr, P);
//     mpfr_mul(corr, t_half, logpi, MPFR_RNDN);
//     mpfr_sub(theta, im, corr, MPFR_RNDN);

//     // cleanup
//     mpfr_clears(half, quarter, logpi, ln2pi_half, two_pi, ln_two_pi, t_half, gmh, im, corr, (mpfr_ptr)0);
//     for (int i = 0; i < 9; ++i)
//         mpfr_clear(p[i]);
//     mpc_clear(z);
//     mpc_clear(sum);
//     mpc_clear(term);
//     mpc_clear(denom);
//     mpc_clear(x);
//     mpc_clear(z_minus_half);
//     mpc_clear(t1);
//     mpc_clear(t2);
//     mpc_clear(t3);
//     mpc_clear(lg);
// }

// Bernoulli numbers B_2, B_4, ..., B_20  (extend if you want more headroom)
static const long BNUM[] = {
    1,   -1,     1,     -1,      5,
   -691,  7,  -3617,  43867, -174611
};
static const long BDEN[] = {
     6,   30,    42,     30,     66,
  2730,    6,    510,    798,     330
};
// These give B_2..B_20. You can add B_22 etc. if you like.

// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
// Use Stirling with optimal truncation; evaluate at P+64 then return at P.
void compute_theta_stirling_opt(mpfr_t theta, mpfr_srcptr t)
{
    const mpfr_prec_t P = CFG().mpfr_prec;
    const mpfr_prec_t Q = P + 64;

    // z = a + i b with a = 1/4, b = t/2
    mpfr_t a, b, r, argz, logr;
    mpfr_inits2(Q, a, b, r, argz, logr, (mpfr_ptr)0);
    mpfr_set_d(a, 0.25, MPFR_RNDN);
    mpfr_div_ui(b, t, 2, MPFR_RNDN);                 // b = t/2

    // r = |z|, argz = atan2(b,a), log r
    mpfr_hypot(r, a, b, MPFR_RNDN);
    mpfr_atan2(argz, b, a, MPFR_RNDN);
    mpfr_log(logr, r, MPFR_RNDN);

    // s_im = (a-1/2)*argz + b*log r - b
    mpfr_t half, amhalf, sim, tmp;
    mpfr_inits2(Q, half, amhalf, sim, tmp, (mpfr_ptr)0);
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpfr_sub(amhalf, a, half, MPFR_RNDN);            // a-1/2 = -1/4
    mpfr_mul(sim, amhalf, argz, MPFR_RNDN);          // (-1/4)*argz
    mpfr_mul(tmp, b, logr, MPFR_RNDN);               // b*log r
    mpfr_add(sim, sim, tmp, MPFR_RNDN);
    mpfr_sub(sim, sim, b, MPFR_RNDN);                // ... - b

    // Bernoulli series: ∑ B_{2k}/[2k(2k-1)] * Im(z^{-(2k-1)})
    // Im(z^{-(2k-1)}) = - r^{-(2k-1)} sin((2k-1) argz)
    mpfr_t inv, inv_r2, rpow, ang, s, coeff, mag, cutoff;
    mpfr_inits2(Q, inv, inv_r2, rpow, ang, s, coeff, mag, cutoff, (mpfr_ptr)0);

    mpfr_ui_div(inv, 1, r, MPFR_RNDN);               // inv = 1/r
    mpfr_mul(inv_r2, inv, inv, MPFR_RNDN);           // 1/r^2
    mpfr_set(rpow, inv, MPFR_RNDN);                  // r^{-1}

    // cutoff ≈ 2^{-(P+10)} as a simple magnitude stop
    mpfr_set_ui(cutoff, 1, MPFR_RNDN);
    mpfr_div_2si(cutoff, cutoff, P + 10, MPFR_RNDN);

    const int KMAX = (int)(sizeof(BNUM)/sizeof(BNUM[0]));
    for (int k = 1; k <= KMAX; ++k)
    {
        // magnitude check: |B_{2k}|/(2k(2k-1)) * r^{-(2k-1)} < cutoff ?
        mpfr_set_si(mag, BNUM[k-1], MPFR_RNDN);
        mpfr_div_si(mag, mag, BDEN[k-1], MPFR_RNDN);
        mpfr_abs(mag, mag, MPFR_RNDN);
        mpfr_div_ui(mag, mag, (unsigned)(2*k*(2*k-1)), MPFR_RNDN);
        mpfr_mul(mag, mag, rpow, MPFR_RNDN);
        if (mpfr_cmp(mag, cutoff) < 0)
            break;

        // ang = (2k-1)*argz ; s = sin(ang)
        mpfr_mul_ui(ang, argz, 2*k - 1, MPFR_RNDN);
        mpfr_sin(s, ang, MPFR_RNDN);

        // coeff = B_{2k}/[2k(2k-1)] * r^{-(2k-1)} * (-sin(ang))
        mpfr_set_si(coeff, BNUM[k-1], MPFR_RNDN);
        mpfr_div_si(coeff, coeff, BDEN[k-1], MPFR_RNDN);
        mpfr_div_ui(coeff, coeff, (unsigned)(2*k*(2*k-1)), MPFR_RNDN);
        mpfr_mul(coeff, coeff, rpow, MPFR_RNDN);
        mpfr_neg(s, s, MPFR_RNDN);
        mpfr_mul(coeff, coeff, s, MPFR_RNDN);
        mpfr_add(sim, sim, coeff, MPFR_RNDN);

        // advance rpow *= r^{-2}
        mpfr_mul(rpow, rpow, inv_r2, MPFR_RNDN);
    }

    // subtract (t/2) log π
    mpfr_t pi2, logpi, tp2;
    mpfr_inits2(Q, pi2, logpi, tp2, (mpfr_ptr)0);
    mpfr_const_pi(pi2, MPFR_RNDN);
    mpfr_log(logpi, pi2, MPFR_RNDN);
    mpfr_mul(tp2, t, logpi, MPFR_RNDN);
    mpfr_div_ui(tp2, tp2, 2, MPFR_RNDN);
    mpfr_sub(sim, sim, tp2, MPFR_RNDN);

    // return at precision P
    mpfr_prec_round(theta, P, MPFR_RNDN);
    mpfr_set(theta, sim, MPFR_RNDN);

    mpfr_clears(a,b,r,argz,logr,half,amhalf,sim,tmp,
                inv,inv_r2,rpow,ang,s,coeff,mag,cutoff,
                pi2,logpi,tp2, (mpfr_ptr)0);
}