#include <mpfr.h>
#include <mpc.h>
#include <vector>
#include <cmath>
#include "settings.h"
#include "z_theta.h"

// ------------------------------------------------------------
// η(s) via binomial/Euler transform (stable)
//   η(s) = sum_{m>=0} [ 1/2^{m+1} * sum_{k=0}^m (-1)^k C(m,k) (k+1)^{-s} ]
// precompute a[k] = (k+1)^{-s} using the stable recurrence
//   a_{k+1} = a_k * exp( -s * log1p(1/(k+1)) )
// Then accumulate outer m until the m-term is < tol (after a few warmups).
// M is interpreted as a *hard cap* (M_max).
// ------------------------------------------------------------
static void zeta_via_eta(
    mpc_t zeta,
    const mpc_t s,
    int M_max = 2000,
    double tol = 1e-40)
{
    const int M_min = 12; // don't stop too early

    mpfr_prec_t P = CFG().mpfr_prec;

    // --- Precompute a[k] = (k+1)^{-s} for k=0..M_max
    std::vector<mpc_t> a(M_max + 1);

    for (int k = 0; k <= M_max; ++k)
        mpc_init2(a[k], P);

    // a[0] = 1^{-s} = 1
    mpc_set_ui(a[0], 1, MPC_RNDNN);

    mpc_t ratio, tmpc;
    mpc_init2(ratio, P);
    mpc_init2(tmpc, P);

    mpfr_t ln1p;
    mpfr_init2(ln1p, P);

    for (int k = 0; k < M_max; ++k)
    {
        // ratio = exp( -s * log1p(1/(k+1)) )
        mpfr_set_ui(ln1p, (unsigned long)(k + 1), MPFR_RNDN);
        mpfr_ui_div(ln1p, 1, ln1p, MPFR_RNDN); // 1/(k+1)
        mpfr_log1p(ln1p, ln1p, MPFR_RNDN);
        mpc_mul_fr(tmpc, s, ln1p, MPC_RNDNN);
        mpc_neg(tmpc, tmpc, MPC_RNDNN);
        mpc_exp(ratio, tmpc, MPC_RNDNN);

        // a[k+1] = a[k] * ratio = (k+2)^{-s}
        mpc_mul(a[k + 1], a[k], ratio, MPC_RNDNN);
    }

    // --- Outer sum over m
    mpc_t eta, inner, mterm;
    mpc_init2(eta, P);
    mpc_set_ui(eta, 0, MPC_RNDNN);
    mpc_init2(inner, P);
    mpc_init2(mterm, P);

    mpfr_t Cm, pow2, abst; // C(m,0)=1 and pow2=2^(m+1)
    mpfr_inits2(P, Cm, pow2, abst, (mpfr_ptr)0);

    mpfr_set_ui(pow2, 2, MPFR_RNDN);

    for (int m = 0; m <= M_max; ++m)
    {
        // inner = sum_{k=0}^m (-1)^k C(m,k) a[k]
        // compute C(m,k) iteratively in mpfr: C(m,0)=1; C(m,k+1)=C(m,k)*(m-k)/(k+1)
        mpfr_set_ui(Cm, 1, MPFR_RNDN);
        mpc_set_ui(inner, 0, MPC_RNDNN);

        for (int k = 0; k <= m; ++k)
        {
            // term = C(m,k) * a[k] * ((-1)^k)
            mpc_mul_fr(mterm, a[k], Cm, MPC_RNDNN);
            if (k & 1)
                mpc_neg(mterm, mterm, MPC_RNDNN);
            mpc_add(inner, inner, mterm, MPC_RNDNN);

            if (k < m)
            {
                // Cm *= (m-k)/(k+1)
                mpfr_mul_ui(Cm, Cm, (unsigned long)(m - k), MPFR_RNDN);
                mpfr_div_ui(Cm, Cm, (unsigned long)(k + 1), MPFR_RNDN);
            }
        }

        // mterm = inner / 2^(m+1)
        mpc_div_fr(mterm, inner, pow2, MPC_RNDNN);
        mpc_add(eta, eta, mterm, MPC_RNDNN);

        // stop criterion: |mterm| < tol after a warmup
        mpc_abs(abst, mterm, MPFR_RNDN);
        if (m >= M_min && mpfr_cmp_d(abst, tol) < 0)
            break;

        // pow2 *= 2  (next outer term)
        mpfr_mul_ui(pow2, pow2, 2, MPFR_RNDN);
    }

    // --- ζ(s) = η(s) / (1 - 2^{1-s})
    mpfr_t ln2;
    mpfr_init2(ln2, P);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_const_log2(ln2, MPFR_RNDN);
#else
    mpfr_set_ui(ln2, 2, MPFR_RNDN);
    mpfr_log(ln2, ln2, MPFR_RNDN);
#endif

    // tmpc = 2^{1-s} = 2 * exp(-s * ln 2)
    mpc_mul_fr(tmpc, s, ln2, MPC_RNDNN);
    mpc_neg(tmpc, tmpc, MPC_RNDNN);
    mpc_exp(tmpc, tmpc, MPC_RNDNN);
    mpc_mul_ui(tmpc, tmpc, 2, MPC_RNDNN);

    // zeta = eta / (1 - 2^{1-s})
    mpc_ui_sub(mterm, 1, tmpc, MPC_RNDNN); // mterm = 1 - 2^{1-s}
    mpc_div(zeta, eta, mterm, MPC_RNDNN);

    // cleanup
    for (auto &v : a)
        mpc_clear(v);

    mpc_clear(ratio);
    mpc_clear(tmpc);
    mpc_clear(eta);
    mpc_clear(inner);
    mpc_clear(mterm);
    mpfr_clears(ln1p, Cm, pow2, abst, ln2, (mpfr_ptr)0);
}

// ------------------------------------------------------------
// Z(t) = Re( e^{iθ(t)} ζ(1/2 + i t) )
// Use η(s) from above and θ from the Stirling wrapper.
// ------------------------------------------------------------
void compute_Z_small(mpfr_t Z, const mpfr_t t)
{
    mpfr_prec_t P = CFG().mpfr_prec;

    mpc_t s, z;
    mpc_init2(s, P);
    mpc_init2(z, P);
    
    mpfr_t half, theta, cth, sth, zr, zi, tmp;
    mpfr_inits2(P, half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);

    // s = 1/2 + i t
    mpfr_set_d(half, 0.5, MPFR_RNDN);
    mpc_set_fr_fr(s, half, t, MPC_RNDNN);

    // ζ(s) via binomial/Euler η(s)
    // (You can tune max M and tol in settings if you add them.)
    zeta_via_eta(z, s, /*M_max*/ 2000, /*tol*/ 1e-40);

    // θ(t): currently using your Stirling-based function
    compute_theta_stirling(theta, t);

    // Z(t) = Re( (cosθ + i sinθ)(zr + i zi) ) = zr*cosθ - zi*sinθ
    mpfr_cos(cth, theta, MPFR_RNDN);
    mpfr_sin(sth, theta, MPFR_RNDN);
    mpfr_set(zr, mpc_realref(z), MPFR_RNDN);
    mpfr_set(zi, mpc_imagref(z), MPFR_RNDN);

    mpfr_mul(Z, zr, cth, MPFR_RNDN);
    mpfr_mul(tmp, zi, sth, MPFR_RNDN);
    mpfr_sub(Z, Z, tmp, MPFR_RNDN);

    mpfr_clears(half, theta, cth, sth, zr, zi, tmp, (mpfr_ptr)0);
    mpc_clear(s);
    mpc_clear(z);
}
