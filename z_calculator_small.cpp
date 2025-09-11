#include <mpfr.h>
#include <mpc.h>
#include <thread>
#include <vector>
#include <cmath>
#include <iostream>
#include <atomic>
#include "settings.h"
#include "z_theta.h"

// Fast η(s) via Euler transform of the alternating series
// η(s) = sum_{n>=1} (-1)^{n-1} n^{-s}
// Euler transform: η = sum_{k>=0} Δ^k a_1 / 2^{k+1},  a_n = n^{-s}
// We build a_1..a_{M+1}, then forward differences in-place.
static void zeta_via_eta(
    mpc_t zeta,
    const mpc_t s,
    int M = 128,
    double tol = 1e-30)
{
    // Build a_j = (j+1)^{-s} for j=0..M using the stable recurrence
    // a_{n+1} = a_n * exp(-s * log1p(1/n))
    std::vector<mpc_t> cur(M + 1);
    for (int j = 0; j <= M; ++j)
        mpc_init2(cur[j], CFG().mpfr_prec);

    mpc_t ratio, tmpc;
    mpc_init2(ratio, CFG().mpfr_prec);
    mpc_init2(tmpc, CFG().mpfr_prec);

    mpfr_t ln1p;
    mpfr_init2(ln1p, CFG().mpfr_prec);

    // a_1 = 1^{-s} = 1
    mpc_set_ui(cur[0], 1, MPC_RNDNN);

    for (int j = 0; j < M; ++j)
    {
        // ratio = exp( -s * log1p(1/(j+1)) )
        mpfr_set_ui(ln1p, (unsigned long)(j + 1), MPFR_RNDN);
        mpfr_ui_div(ln1p, 1, ln1p, MPFR_RNDN); // 1/(j+1)
        mpfr_log1p(ln1p, ln1p, MPFR_RNDN);     // log(1 + 1/(j+1))
        mpc_mul_fr(tmpc, s, ln1p, MPC_RNDNN);  // s * log(...)
        mpc_neg(tmpc, tmpc, MPC_RNDNN);        // -s * log(...)
        mpc_exp(ratio, tmpc, MPC_RNDNN);       // exp(-s * log(...))

        // cur[j+1] = cur[j] * ratio = (j+2)^{-s}
        mpc_mul(cur[j + 1], cur[j], ratio, MPC_RNDNN);
    }

    // η via forward differences
    mpc_t eta, term;
    mpc_init2(eta, CFG().mpfr_prec);
    mpc_set_ui(eta, 0, MPC_RNDNN);
    mpc_init2(term, CFG().mpfr_prec);

    mpfr_t pow2, abst;
    mpfr_inits2(CFG().mpfr_prec, pow2, abst, (mpfr_ptr)0);
    mpfr_set_ui(pow2, 2, MPFR_RNDN); // 2^{1} for k=0 term

    int L = M + 1; // current length of the difference row
    for (int k = 0; k <= M && L > 0; ++k)
    {
        // term_k = cur[0] / 2^{k+1}
        mpc_div_fr(term, cur[0], pow2, MPC_RNDNN);
        mpc_add(eta, eta, term, MPC_RNDNN);

        // stop when |term_k| < tol (converged)
        mpc_abs(abst, term, MPFR_RNDN);
        if (mpfr_cmp_d(abst, tol) < 0)
            break;

        // Next row of forward differences in-place: cur[j] = cur[j+1] - cur[j]
        for (int j = 0; j < L - 1; ++j)
        {
            mpc_sub(cur[j], cur[j + 1], cur[j], MPC_RNDNN);
        }
        --L;
        mpfr_mul_ui(pow2, pow2, 2, MPFR_RNDN); // pow2 *= 2
    }

    // zeta = eta / (1 - 2^{1-s})
    mpfr_t ln2;
    mpfr_init2(ln2, CFG().mpfr_prec);
    // Portable ln(2): avoid relying on mpfr_const_log2 availability
    mpfr_set_ui(ln2, 2, MPFR_RNDN);
    mpfr_log(ln2, ln2, MPFR_RNDN);

    // tmpc = 2^{1-s} = 2 * 2^{-s} = 2 * exp(-s ln 2)
    mpc_mul_fr(tmpc, s, ln2, MPC_RNDNN);
    mpc_neg(tmpc, tmpc, MPC_RNDNN);
    mpc_exp(tmpc, tmpc, MPC_RNDNN);       // 2^{-s}
    mpc_mul_ui(tmpc, tmpc, 2, MPC_RNDNN); // 2^{1-s}
    mpc_ui_sub(tmpc, 1, tmpc, MPC_RNDNN); // 1 - 2^{1-s}
    mpc_div(zeta, eta, tmpc, MPC_RNDNN);

    // cleanup
    for (auto &v : cur)
        mpc_clear(v);
    mpc_clear(ratio);
    mpc_clear(tmpc);
    mpc_clear(eta);
    mpc_clear(term);

    mpfr_clears(ln1p, pow2, abst, ln2, (mpfr_ptr)0);
}

// Z(t) = Re( e^{iθ(t)} ζ(1/2 + i t) )
// For small/medium t: use η(s) (Euler transform) + θ via Stirling wrapper.
void compute_Z_small(mpfr_t Z, const mpfr_t t)
{
    mpc_t s, z;
    mpc_init2(s, CFG().mpfr_prec);
    mpc_init2(z, CFG().mpfr_prec);

    mpfr_t half, theta, cth, sth, zr, zi, tmp;
    mpfr_inits2(CFG().mpfr_prec, half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);

    // s = 1/2 + i t
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpc_set_fr_fr(s, half, t, MPC_RNDNN);

    // ζ(s) via fast η(s), and θ(t) via Stirling (configurable K)
    zeta_via_eta(z, s, /*M*/ 160, /*tol*/ 1e-34);
    compute_theta_stirling(theta, t);

    // Z(t) = Re( e^{iθ} ζ ) = Re( (cosθ + i sinθ)(zr + i zi) ) = zr*cosθ - zi*sinθ
    mpfr_cos(cth, theta, MPFR_RNDN);
    mpfr_sin(sth, theta, MPFR_RNDN);
    mpfr_set(zr, mpc_realref(z), MPFR_RNDN);
    mpfr_set(zi, mpc_imagref(z), MPFR_RNDN);

    mpfr_mul(Z, zr, cth, MPFR_RNDN);
    mpfr_mul(tmp, zi, sth, MPFR_RNDN);
    mpfr_sub(Z, Z, tmp, MPFR_RNDN);

    // cleanup
    mpfr_clears(half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);
    mpc_clear(s);
    mpc_clear(z);
}