#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <omp.h>
#include <sstream>
#include "settings.h"
#include "logs.h"
#include "utils.h"

// -------------------------------
// Constants: precompute π, 2π, π², π⁴, log π once per precision
// -------------------------------
struct Constants
{
    mpfr_t pi, two_pi, pi2, pi4, log_pi;
    mpfr_prec_t P;

    explicit Constants(mpfr_prec_t prec = 0) : P(0)
    {
        mpfr_init2(pi, 2); // temporary small; set in init
        mpfr_init2(two_pi, 2);
        mpfr_init2(pi2, 2);
        mpfr_init2(pi4, 2);
        mpfr_init2(log_pi, 2);
        if (prec)
            init(prec);
    }

    void init(mpfr_prec_t prec)
    {
        if (P == prec)
            return;
        P = prec;
        mpfr_set_prec(pi, P);
        mpfr_set_prec(two_pi, P);
        mpfr_set_prec(pi2, P);
        mpfr_set_prec(pi4, P);
        mpfr_set_prec(log_pi, P);

        mpfr_const_pi(pi, MPFR_RNDN);          // π
        mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN); // 2π
        mpfr_mul(pi2, pi, pi, MPFR_RNDN);      // π²
        mpfr_mul(pi4, pi2, pi2, MPFR_RNDN);    // π⁴
        mpfr_log(log_pi, pi, MPFR_RNDN);       // log π
    }

    ~Constants()
    {
        mpfr_clear(pi);
        mpfr_clear(two_pi);
        mpfr_clear(pi2);
        mpfr_clear(pi4);
        mpfr_clear(log_pi);
    }
};

// -------------------------------
// Spouge coeffs c_k  (real)
// c0 = sqrt(2π)
// c_k = (-1)^(k-1)/(k-1)! * (a-k)^(k-1/2) * e^(a-k),  k=1..a-1
// -------------------------------
struct Spouge
{
    int a;
    std::vector<__mpfr_struct> c; // store POD, not mpfr_t
    mpfr_t sqrt2pi;

    Spouge() : a(0)
    {
        mpfr_init2(sqrt2pi, 128);
        mpfr_const_pi(sqrt2pi, MPFR_RNDN);
        mpfr_mul_ui(sqrt2pi, sqrt2pi, 2, MPFR_RNDN);
        mpfr_sqrt(sqrt2pi, sqrt2pi, MPFR_RNDN);
    }
    ~Spouge()
    {
        for (auto &ci : c)
            mpfr_clear(&ci);
        mpfr_clear(sqrt2pi);
    }

    void build(int A, mpfr_prec_t P)
    {
        a = A;
        for (auto &ci : c)
            mpfr_clear(&ci); // in case build() called twice
        c.clear();
        c.resize(a);
        for (int k = 0; k < a; ++k)
            mpfr_init2(&c[k], P);
        mpfr_set(&c[0], sqrt2pi, MPFR_RNDN);

        mpfr_t tmp1, tmp2, tmp3;
        mpfr_inits2(P, tmp1, tmp2, tmp3, (mpfr_ptr)0);

        for (int k = 1; k < a; ++k)
        {
            int sgn = (k % 2 == 1) ? +1 : -1;

            mpfr_set_si(tmp1, a - k, MPFR_RNDN); // (a-k)
            mpfr_set_si(tmp2, k, MPFR_RNDN);
            mpfr_sub_d(tmp2, tmp2, 0.5, MPFR_RNDN); // k - 1/2
            mpfr_pow(tmp1, tmp1, tmp2, MPFR_RNDN);  // (a-k)^(k-1/2)

            mpfr_set_si(tmp2, a - k, MPFR_RNDN);
            mpfr_exp(tmp2, tmp2, MPFR_RNDN);       // e^(a-k)
            mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN); // *

            if (sgn < 0)
                mpfr_neg(tmp1, tmp1, MPFR_RNDN);

            mpfr_fac_ui(tmp2, k - 1, MPFR_RNDN);    // (k-1)!
            mpfr_div(&c[k], tmp1, tmp2, MPFR_RNDN); // <- address!
        }
        mpfr_clears(tmp1, tmp2, tmp3, (mpfr_ptr)0);
    }
};

// -------------------------------
// Optional: 9-point 8th-order central difference weights
// Solve once per derivative order m ∈ {2,3,6}.
// -------------------------------
struct DiffWeights
{
    int r;
    std::vector<__mpfr_struct> w2, w3, w6;
    mpfr_prec_t P;

    DiffWeights() : r(4), P(0) {}
    ~DiffWeights() { clear(); }

    void clear()
    {
        for (auto *arr : {&w2, &w3, &w6})
        {
            for (auto &v : *arr)
                mpfr_clear(&v);
            arr->clear();
        }
    }

    void build(mpfr_prec_t prec)
    {
        P = prec;
        clear();
        w2.resize(2 * r + 1);
        w3.resize(2 * r + 1);
        w6.resize(2 * r + 1);
        for (int i = 0; i < 2 * r + 1; ++i)
        {
            mpfr_init2(&w2[i], P);
            mpfr_init2(&w3[i], P);
            mpfr_init2(&w6[i], P);
        }
        compute_for_m(2, w2);
        compute_for_m(3, w3);
        compute_for_m(6, w6);
    }

    void compute_for_m(int m, std::vector<__mpfr_struct> &out)
    {
        const int n = 2 * r;
        auto idx = [&](int i, int j)
        { return i * (2 * r + 1) + j; };

        std::vector<__mpfr_struct> A((n + 1) * (2 * r + 1)), b(n + 1), x(2 * r + 1);
        for (auto &v : A)
            mpfr_init2(&v, P);
        for (auto &v : b)
            mpfr_init2(&v, P);
        for (auto &v : x)
            mpfr_init2(&v, P);

        // A[i,j] = (j-r)^i
        for (int i = 0; i <= n; ++i)
        {
            for (int j = 0; j <= 2 * r; ++j)
            {
                int jr = j - r;
                if (i == 0)
                {
                    mpfr_set_ui(&A[idx(i, j)], 1, MPFR_RNDN);
                }
                else
                {
                    mpfr_set_si(&A[idx(i, j)], jr, MPFR_RNDN);
                    for (int k = 1; k < i; ++k)
                        mpfr_mul_si(&A[idx(i, j)], &A[idx(i, j)], jr, MPFR_RNDN);
                }
            }
        }
        for (int i = 0; i <= n; ++i)
            mpfr_set_zero(&b[i], 0);
        mpfr_fac_ui(&b[m], m, MPFR_RNDN);

        // Augmented matrix [A|b]
        const int ROWS = n + 1, COLS = 2 * r + 2;
        std::vector<__mpfr_struct> M(ROWS * COLS);
        for (auto &v : M)
            mpfr_init2(&v, P);

        for (int i = 0; i < ROWS; ++i)
        {
            for (int j = 0; j < 2 * r + 1; ++j)
                mpfr_set(&M[i * COLS + j], &A[idx(i, j)], MPFR_RNDN);
            mpfr_set(&M[i * COLS + (COLS - 1)], &b[i], MPFR_RNDN);
        }

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < ROWS; ++i)
        {
            int p = i;
            for (int r2 = i + 1; r2 < ROWS; ++r2)
                if (mpfr_cmpabs(&M[r2 * COLS + i], &M[p * COLS + i]) > 0)
                    p = r2;
            if (p != i)
                for (int j = 0; j < COLS; ++j)
                    mpfr_swap(&M[p * COLS + j], &M[i * COLS + j]);

            mpfr_t piv;
            mpfr_init2(piv, P);
            mpfr_set(piv, &M[i * COLS + i], MPFR_RNDN);
            for (int j = i; j < COLS; ++j)
                mpfr_div(&M[i * COLS + j], &M[i * COLS + j], piv, MPFR_RNDN);
            mpfr_clear(piv);

            for (int r2 = 0; r2 < ROWS; ++r2)
            {
                if (r2 == i)
                    continue;
                mpfr_t factor;
                mpfr_init2(factor, P);
                mpfr_set(factor, &M[r2 * COLS + i], MPFR_RNDN);
                if (mpfr_zero_p(factor) == 0)
                {
                    for (int j = i; j < COLS; ++j)
                    {
                        mpfr_t tmp;
                        mpfr_init2(tmp, P);
                        mpfr_mul(tmp, factor, &M[i * COLS + j], MPFR_RNDN);
                        mpfr_sub(&M[r2 * COLS + j], &M[r2 * COLS + j], tmp, MPFR_RNDN);
                        mpfr_clear(tmp);
                    }
                }
                mpfr_clear(factor);
            }
        }

        for (int j = 0; j < 2 * r + 1; ++j)
            mpfr_set(&x[j], &M[j * COLS + (COLS - 1)], MPFR_RNDN);
        for (int j = 0; j < 2 * r + 1; ++j)
            mpfr_set(&out[j], &x[j], MPFR_RNDN);

        for (auto &v : A)
            mpfr_clear(&v);
        for (auto &v : b)
            mpfr_clear(&v);
        for (auto &v : x)
            mpfr_clear(&v);
        for (auto &v : M)
            mpfr_clear(&v);
    }
};

// -------------------------------
// precision helpers
// -------------------------------
static inline long dec_digits_from_bits(mpfr_prec_t P)
{
    // ~ P * log10(2)
    return (long)std::ceil(P * 0.3010299956639811952137389L);
}

// -------------------------------
// RS parameters and helpers
// a = sqrt(t/(2π)),  N = floor(a),  p = frac(a)
// -------------------------------
static long rs_N_from_t(const mpfr_t t, mpfr_prec_t P, mpfr_t a_out, mpfr_t p_out, const Constants &C)
{
    mpfr_t tmp;
    mpfr_init2(tmp, P);
    mpfr_div(tmp, t, C.two_pi, MPFR_RNDN);
    mpfr_sqrt(a_out, tmp, MPFR_RNDN); // a = sqrt(t/(2π))
    long N = mpfr_get_si(a_out, MPFR_RNDN);
    if (N < 1)
        N = 1;
    // p = a - floor(a)
    mpfr_set_si(tmp, N, MPFR_RNDN);
    mpfr_sub(p_out, a_out, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
    return N;
}

// Parallel version of Z_rs_main (unchanged logic, but uses Constants via rs_N_from_t)
static void Z_rs_main_parallel(
    mpfr_t Z,
    const mpfr_t t,
    const mpfr_t theta,
    mpfr_prec_t P,
    const Constants &C,
    int max_threads = 0,
    long block_size = 0) // <-- default 0 to enable dynamic sizing
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    if (max_threads <= 0)
        max_threads = omp_get_max_threads();

    // --- Dynamic block size: sqrt(N)*16, clamped to [4096, 32768]
    if (block_size <= 0)
    {
        long dyn = static_cast<long>(std::sqrt(static_cast<long double>(N)) * 16.0L);
        if (dyn < 4096)
            dyn = 4096;
        if (dyn > 32768)
            dyn = 32768;
        block_size = dyn;
    }

    const long num_blocks = (N + block_size - 1) / block_size;

    mpfr_set_zero(Z, 0);

    std::vector<mpfr_t> partial(max_threads);
    for (int i = 0; i < max_threads; ++i)
        mpfr_init2(partial[i], P), mpfr_set_zero(partial[i], 0);

#pragma omp parallel for schedule(static)
    for (long b = 0; b < num_blocks; ++b)
    {
        const int tid = omp_get_thread_num();
        mpfr_t Sloc;
        mpfr_init2(Sloc, P);
        mpfr_set_zero(Sloc, 0);

        long n0 = b * block_size + 1;
        long n1 = std::min(n0 + block_size - 1, N);

        mpfr_t n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp;
        mpfr_inits2(P, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);

        mpfr_set_ui(n_mp, (unsigned long)n0, MPFR_RNDN);
        mpfr_log(phi0, n_mp, MPFR_RNDN);
        mpfr_mul(phi0, phi0, t, MPFR_RNDN);
        mpfr_neg(phi0, phi0, MPFR_RNDN);
        mpfr_add(phi0, phi0, theta, MPFR_RNDN);

        mpfr_sin_cos(sphi, cphi, phi0, MPFR_RNDN);

        mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);

        for (long i = n0; i < n1; ++i)
        {
            mpfr_set_ui(n_mp, (unsigned long)i, MPFR_RNDN);
            mpfr_ui_div(dphi, 1, n_mp, MPFR_RNDN);
            mpfr_log1p(dphi, dphi, MPFR_RNDN);
            mpfr_mul(dphi, dphi, t, MPFR_RNDN);

            mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

            mpfr_t new_c, new_s;
            mpfr_inits2(P, new_c, new_s, (mpfr_ptr)0);
            mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
            mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

            mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
            mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

            mpfr_set(cphi, new_c, MPFR_RNDN);
            mpfr_set(sphi, new_s, MPFR_RNDN);
            mpfr_clears(new_c, new_s, (mpfr_ptr)0);

            mpfr_set_ui(n_mp, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
            mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
            mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);
        }

        mpfr_add(partial[tid], partial[tid], Sloc, MPFR_RNDN);
        mpfr_clears(Sloc, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);
    }

    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_add(Z, Z, partial[i], MPFR_RNDN);
        mpfr_clear(partial[i]);
    }

    mpfr_mul_ui(Z, Z, 2, MPFR_RNDN);
    mpfr_clears(a, p, (mpfr_ptr)0);
}

static inline int spouge_a_for_prec(mpfr_prec_t P)
{
    int digs = (int)dec_digits_from_bits(P);
    int a = (int)std::lround(0.45 * digs);
    if (a < 30)
        a = 30;
    if (a > 80)
        a = 80;
    return a;
}

// -------------------------------
// log Γ via Spouge (right half) and reflection
// -------------------------------
static void mpc_loggamma_spouge_shift(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P)
{
    mpc_t w, s, term;
    mpc_init2(w, P);
    mpc_init2(s, P);
    mpc_init2(term, P);
    mpc_set(w, z, MPC_RNDNN);
    mpfr_add_si(mpc_realref(w), mpc_realref(w), S.a - 1, MPFR_RNDN);

    mpc_set_fr(s, &S.c[0], MPC_RNDNN);
    for (int k = 1; k < S.a; ++k)
    {
        mpc_set(term, z, MPC_RNDNN);
        mpfr_add_si(mpc_realref(term), mpc_realref(term), k - 1, MPFR_RNDN);
        mpc_ui_div(term, 1, term, MPC_RNDNN);
        mpc_mul_fr(term, term, &S.c[k], MPC_RNDNN);
        mpc_add(s, s, term, MPC_RNDNN);
    }

    mpc_log(out, w, MPC_RNDNN);
    mpc_t zmh;
    mpc_init2(zmh, P);
    mpc_set(zmh, z, MPC_RNDNN);
    mpfr_sub_d(mpc_realref(zmh), mpc_realref(zmh), 0.5, MPFR_RNDN);

    mpc_mul(out, zmh, out, MPC_RNDNN);
    mpc_sub(out, out, w, MPC_RNDNN);

    mpc_log(term, s, MPC_RNDNN);
    mpc_add(out, out, term, MPC_RNDNN);

    mpc_clear(zmh);
    mpc_clear(w);
    mpc_clear(s);
    mpc_clear(term);
}

static void mpc_loggamma_spouge_reflect(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P, const Constants &C)
{
    mpc_t one_minus_z, lg1mz, sin_pi_z, t;
    mpc_init2(one_minus_z, P);
    mpc_init2(lg1mz, P);
    mpc_init2(sin_pi_z, P);
    mpc_init2(t, P);

    mpc_set_si(one_minus_z, 1, MPC_RNDNN);
    mpc_sub(one_minus_z, one_minus_z, z, MPC_RNDNN);

    mpc_loggamma_spouge_shift(lg1mz, one_minus_z, S, P);

    // t = π (real)
    mpc_set_fr(t, C.pi, MPC_RNDNN);
    mpc_mul(t, t, z, MPC_RNDNN);
    mpc_sin(sin_pi_z, t, MPC_RNDNN);

    // out = log π − log(sin π z) − lg(1−z)
    mpc_log(out, sin_pi_z, MPC_RNDNN);
    mpc_neg(out, out, MPC_RNDNN);
    mpc_add_fr(out, out, C.log_pi, MPC_RNDNN);
    mpc_sub(out, out, lg1mz, MPC_RNDNN);

    mpc_clear(one_minus_z);
    mpc_clear(lg1mz);
    mpc_clear(sin_pi_z);
    mpc_clear(t);
}

// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
static void compute_theta_spouge(mpfr_t theta, const mpfr_t t, mpfr_prec_t P, const Spouge &S, const Constants &C)
{
    mpc_t z, lg;
    mpc_init2(z, P);
    mpc_init2(lg, P);

    mpfr_t tr;
    mpfr_init2(tr, P);
    mpfr_div_ui(tr, t, 2, MPFR_RNDN);
    mpc_set_d_d(z, 0.25, 0.0, MPC_RNDNN);
    mpfr_set(mpc_imagref(z), tr, MPFR_RNDN);

    mpc_loggamma_spouge_reflect(lg, z, S, P, C);

    mpfr_set(theta, mpc_imagref(lg), MPFR_RNDN);
    mpfr_t tmp;
    mpfr_init2(tmp, P);
    mpfr_mul(tmp, C.log_pi, tr, MPFR_RNDN); // (t/2) log π
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);

    mpfr_clear(tr);
    mpfr_clear(tmp);
    mpc_clear(z);
    mpc_clear(lg);
}

// -------------------------------
// ψ(p) = cos(2π(p^2 - p - 1/16)) / cos(2π p), with safety near 1/4, 3/4
// -------------------------------
static void psi_p(mpfr_t out, const mpfr_t p, mpfr_prec_t P, const Constants &C)
{
    mpfr_t num_arg, den_arg, den, num;
    mpfr_inits2(P, num_arg, den_arg, den, num, (mpfr_ptr)0);

    // den = cos(2π p)
    mpfr_mul(den_arg, C.two_pi, p, MPFR_RNDN);
    mpfr_cos(den, den_arg, MPFR_RNDN);

    mpfr_t eps;
    mpfr_init2(eps, P);
    mpfr_set_d(eps, 1e-30, MPFR_RNDN);

    if (mpfr_cmpabs(den, eps) > 0)
    {
        // num = cos(2π(p^2 - p - 1/16))
        mpfr_t t1;
        mpfr_init2(t1, P);
        mpfr_sqr(t1, p, MPFR_RNDN);            // p^2
        mpfr_sub(t1, t1, p, MPFR_RNDN);        // p^2 - p
        mpfr_sub_d(t1, t1, 0.0625, MPFR_RNDN); // -1/16
        mpfr_mul(num_arg, C.two_pi, t1, MPFR_RNDN);
        mpfr_cos(num, num_arg, MPFR_RNDN);
        mpfr_div(out, num, den, MPFR_RNDN);
        mpfr_clear(t1);
    }
    else
    {
        mpfr_t q, half, d1, d2;
        mpfr_inits2(P, q, half, d1, d2, (mpfr_ptr)0);
        mpfr_set_d(half, 0.5, MPFR_RNDN);
        mpfr_sub_d(d1, p, 0.25, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_sub_d(d2, p, 0.75, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        if (mpfr_cmp(d1, d2) < 0)
        {
            mpfr_sub_d(q, p, 0.25, MPFR_RNDN);
            mpfr_sub(out, half, q, MPFR_RNDN);
        }
        else
        {
            mpfr_sub_d(q, p, 0.75, MPFR_RNDN);
            mpfr_add(out, half, q, MPFR_RNDN);
        }
        mpfr_clears(q, half, d1, d2, (mpfr_ptr)0);
    }

    mpfr_clear(eps);
    mpfr_clears(num_arg, den_arg, den, num, (mpfr_ptr)0);
}

// central derivatives of ψ at p, 9-point, step h; returns ψ'', ψ''', ψ''''''.
static void psi_derivs_central(mpfr_t d2, mpfr_t d3, mpfr_t d6,
                               const mpfr_t p, mpfr_prec_t P,
                               const DiffWeights &W, const mpfr_t h,
                               const Constants &C)
{
    const int r = W.r;
    std::vector<mpfr_t> f(2 * r + 1);
    for (auto &v : f)
        mpfr_init2(v, P);

    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t arg, step;
        mpfr_inits2(P, arg, step, (mpfr_ptr)0);
        mpfr_set(arg, p, MPFR_RNDN);
        mpfr_mul_si(step, h, j - r, MPFR_RNDN);
        mpfr_add(arg, arg, step, MPFR_RNDN);
        psi_p(f[j], arg, P, C);
        mpfr_clears(arg, step, (mpfr_ptr)0);
    }

    mpfr_set_zero(d2, 0);
    mpfr_set_zero(d3, 0);
    mpfr_set_zero(d6, 0);
    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t t;
        mpfr_init2(t, P);
        mpfr_mul(t, &W.w2[j], f[j], MPFR_RNDN);
        mpfr_add(d2, d2, t, MPFR_RNDN);
        mpfr_mul(t, &W.w3[j], f[j], MPFR_RNDN);
        mpfr_add(d3, d3, t, MPFR_RNDN);
        mpfr_mul(t, &W.w6[j], f[j], MPFR_RNDN);
        mpfr_add(d6, d6, t, MPFR_RNDN);
        mpfr_clear(t);
    }

    mpfr_t h2, h3, h6;
    mpfr_inits2(P, h2, h3, h6, (mpfr_ptr)0);
    mpfr_mul(h2, h, h, MPFR_RNDN);
    mpfr_mul(h3, h2, h, MPFR_RNDN);
    mpfr_mul(h6, h3, h3, MPFR_RNDN);

    mpfr_div(d2, d2, h2, MPFR_RNDN);
    mpfr_div(d3, d3, h3, MPFR_RNDN);
    mpfr_div(d6, d6, h6, MPFR_RNDN);

    mpfr_clears(h2, h3, h6, (mpfr_ptr)0);
    for (auto &v : f)
        mpfr_clear(v);
}

// -------------------------------
// Gabcke correction R1: A0 * [ C0 + C1*(2π/t)^(1/2) + C2*(2π/t) ]
// with sign (-1)^(N-1)
// -------------------------------
static void R_gabcke(mpfr_t R, const mpfr_t t, mpfr_prec_t P,
                     bool useC1C2,
                     DiffWeights *W_opt,
                     const Constants &C)
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    // A0 = (2π/t)^(1/4)
    mpfr_t A0;
    mpfr_init2(A0, P);
    mpfr_div(A0, C.two_pi, t, MPFR_RNDN);
    mpfr_rootn_ui(A0, A0, 4, MPFR_RNDN);

    // C0 = ψ(p)
    mpfr_t C0;
    mpfr_init2(C0, P);
    psi_p(C0, p, P, C);
    mpfr_set(R, C0, MPFR_RNDN);

    if (useC1C2)
    {
        mpfr_t d1, quarter, threeq, d2, h;
        mpfr_inits2(P, d1, quarter, threeq, d2, h, (mpfr_ptr)0);
        mpfr_set_d(quarter, 0.25, MPFR_RNDN);
        mpfr_set_d(threeq, 0.75, MPFR_RNDN);
        mpfr_sub(d1, p, quarter, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_sub(d2, p, threeq, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        if (mpfr_cmp(d1, d2) < 0)
            mpfr_mul_d(h, d1, 0.1, MPFR_RNDN);
        else
            mpfr_mul_d(h, d2, 0.1, MPFR_RNDN);
        if (mpfr_cmp_d(h, 1e-6) < 0)
            mpfr_set_d(h, 1e-6, MPFR_RNDN);
        if (mpfr_cmp_d(h, 1e-3) > 0)
            mpfr_set_d(h, 1e-3, MPFR_RNDN);

        DiffWeights localW;
        DiffWeights *W = W_opt ? W_opt : &localW;
        if (W->w2.empty())
            W->build(P);

        mpfr_t psi2, psi3, psi6;
        mpfr_inits2(P, psi2, psi3, psi6, (mpfr_ptr)0);
        psi_derivs_central(psi2, psi3, psi6, p, P, *W, h, C);

        // C1 = -(1/(96π^2)) ψ'''(p)
        mpfr_t C1;
        mpfr_init2(C1, P);
        mpfr_mul(C1, C.pi2, C.pi2, MPFR_RNDN); // (π²)*(π²) -> temporary; we only need π²*96, do simpler:
        // simpler: C1 = psi3 / (96 * π²)
        mpfr_set(C1, C.pi2, MPFR_RNDN);
        mpfr_mul_ui(C1, C1, 96, MPFR_RNDN);
        mpfr_div(C1, psi3, C1, MPFR_RNDN);
        mpfr_neg(C1, C1, MPFR_RNDN);

        // C2 = ψ^{(6)}/(18432π^4) + ψ''/(64π^2)   (kept exactly as your stable version)
        mpfr_t C2, term1, term2, denom;
        mpfr_inits2(P, C2, term1, term2, denom, (mpfr_ptr)0);

        // term1 = psi6 / (18432 * π^4)
        mpfr_mul_ui(denom, C.pi4, 18432, MPFR_RNDN);
        mpfr_div(term1, psi6, denom, MPFR_RNDN);

        // term2 = psi2 / (64 * π^2)
        mpfr_mul_ui(denom, C.pi2, 64, MPFR_RNDN);
        mpfr_div(term2, psi2, denom, MPFR_RNDN);

        mpfr_add(C2, term1, term2, MPFR_RNDN);

        // scale factors (2π/t)^(1/2) and (2π/t)
        mpfr_t s12, s1, add;
        mpfr_inits2(P, s12, s1, add, (mpfr_ptr)0);
        mpfr_div(s1, C.two_pi, t, MPFR_RNDN);
        mpfr_sqrt(s12, s1, MPFR_RNDN);

        mpfr_mul(add, C1, s12, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);
        mpfr_mul(add, C2, s1, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);

        mpfr_clears(psi2, psi3, psi6, C1, C2, term1, term2, denom, s12, s1, add, d1, d2, quarter, threeq, h, (mpfr_ptr)0);
    }

    mpfr_mul(R, R, A0, MPFR_RNDN);
    if ((N % 2) == 0)
        mpfr_neg(R, R, MPFR_RNDN);

    mpfr_clears(a, p, A0, C0, (mpfr_ptr)0);
}

// -------------------------------
// Main RS sum via rotations, NO half-weight
// -------------------------------
static void Z_rs_main(mpfr_t Z, const mpfr_t t, const mpfr_t theta, mpfr_prec_t P, const Constants &C)
{
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p, C);

    mpfr_t S;
    mpfr_init2(S, P);
    mpfr_set_zero(S, 0);

    mpfr_t cphi, sphi;
    mpfr_inits2(P, cphi, sphi, (mpfr_ptr)0);
    mpfr_sin_cos(sphi, cphi, theta, MPFR_RNDN);
    mpfr_add(S, S, cphi, MPFR_RNDN);

    mpfr_t dphi, sd, cd, new_c, new_s, tmp;
    mpfr_inits2(P, dphi, sd, cd, new_c, new_s, tmp, (mpfr_ptr)0);
    mpfr_t i_mp, sqrt_n;
    mpfr_inits2(P, i_mp, sqrt_n, (mpfr_ptr)0);

    for (long i = 1; i < N; ++i)
    {
        mpfr_set_ui(i_mp, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(dphi, 1, i_mp, MPFR_RNDN);
        mpfr_log1p(dphi, dphi, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        mpfr_set_ui(i_mp, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_n, i_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(S, S, tmp, MPFR_RNDN);
    }

    mpfr_mul_ui(Z, S, 2, MPFR_RNDN);

    mpfr_clears(a, p, S, cphi, sphi, dphi, sd, cd, new_c, new_s, tmp, i_mp, sqrt_n, (mpfr_ptr)0);
}

// -------------------------------
// Top-level: Z_tot = Z_main + R_gabcke
// -------------------------------
static void Z_rs_with_gabcke(
    mpfr_t Ztot,
    mpfr_t Zmain_out,
    mpfr_t R1_out,
    const mpfr_t t,
    mpfr_prec_t P,
    bool useC1C2,
    DiffWeights *W_opt,
    const Spouge &S,
    const Constants &C)
{
    mpfr_t theta;
    mpfr_init2(theta, P);
    compute_theta_spouge(theta, t, P, S, C);

    mpfr_t Zmain;
    mpfr_init2(Zmain, P);

    if (CFG().use_threaded_newton)
    {
        LOGMP("Disablig parallel RS-Gabcke to paralelize Newton...");
        LOGMP("Processing Single Thread RS-Gabcke...");
        Z_rs_main(Zmain, t, theta, P, C);
    }
    else
    {
        if (CFG().use_threaded_z)
        {
            LOGMP_DEBUG("Processing parallelized RS-Gabcke...");
            Z_rs_main_parallel(Zmain, t, theta, P, C, /*threads*/ 0, /*block_size*/ 0);
        }
        else
        {
            LOGMP_DEBUG("Processing Single Thread RS-Gabcke...");
            Z_rs_main(Zmain, t, theta, P, C);
        }
    }

    mpfr_t R1;
    mpfr_init2(R1, P);
    R_gabcke(R1, t, P, useC1C2, W_opt, C);

    mpfr_add(Ztot, Zmain, R1, MPFR_RNDN);

    if (Zmain_out)
        mpfr_set(Zmain_out, Zmain, MPFR_RNDN);
    if (R1_out)
        mpfr_set(R1_out, R1, MPFR_RNDN);

    mpfr_clears(theta, Zmain, R1, (mpfr_ptr)0);
}

void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t)
{
    mpfr_prec_t P = CFG().mpfr_prec;
    std::vector<mpfr_srcptr> tvals;
    tvals.push_back(t);

    mpfr_set_default_prec(P);

    // prepare Spouge coeffs
    Spouge S;
    int a = spouge_a_for_prec(P);
    S.build(a, P);

    // precompute π-related constants once
    Constants C(P);

    DiffWeights W; // build on first use if useC1C2=true

    for (auto &ts : tvals)
    {
        mpfr_t tval;
        mpfr_init2(tval, P);
        mpfr_set(tval, ts, MPFR_RNDN);

        mpfr_t Ztot, Zmain, R1;
        mpfr_inits2(P, Ztot, Zmain, R1, (mpfr_ptr)0);

        bool useC1C2 = true;
        Z_rs_with_gabcke(Ztot, Zmain, R1, tval, P, useC1C2, &W, S, C);

        mpfr_set(Z, Ztot, MPFR_RNDN);

        mpfr_clears(tval, Ztot, Zmain, R1, (mpfr_ptr)0);
    }

    tvals.clear();
}
