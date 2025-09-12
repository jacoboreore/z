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

// local vars
std::ostringstream oss;

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
            mpfr_exp(tmp2, tmp2, MPFR_RNDN); // e^(a-k)
            mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

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
static void mpfr_two_pi(mpfr_t x, mpfr_prec_t P)
{
    mpfr_const_pi(x, MPFR_RNDN);
    mpfr_mul_ui(x, x, 2, MPFR_RNDN);
}

static long rs_N_from_t(const mpfr_t t, mpfr_prec_t P, mpfr_t a_out, mpfr_t p_out)
{
    mpfr_t two_pi, tmp;
    mpfr_inits2(P, two_pi, tmp, (mpfr_ptr)0);
    mpfr_two_pi(two_pi, P);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(a_out, tmp, MPFR_RNDN); // a = sqrt(t/(2π))
    long N = mpfr_get_si(a_out, MPFR_RNDN);
    if (N < 1)
        N = 1;
    // p = a - floor(a)
    mpfr_set_si(tmp, N, MPFR_RNDN);
    mpfr_sub(p_out, a_out, tmp, MPFR_RNDN);
    mpfr_clears(two_pi, tmp, (mpfr_ptr)0);
    return N;
}

// Parallel version of Z_rs_main
static void Z_rs_main_parallel(
    mpfr_t Z,
    const mpfr_t t,
    const mpfr_t theta,
    mpfr_prec_t P,
    int max_threads = 0,
    long block_size = 8192)
{

    //LOG("\nParallel execution, large t Riemann Siegel, Gabcke correction...");
    // N and p only to get N; a/p not needed here beyond that
    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p);

    if (max_threads <= 0)
        max_threads = omp_get_max_threads();

    if (block_size <= 0)
        block_size = 8192; // tune: 4k–32k

    const long num_blocks = (N + block_size - 1) / block_size;

    // Global accumulator
    mpfr_set_zero(Z, 0);

    // Per-thread partial sums (to avoid many critical sections)
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

        // block range [n0, n1]
        long n0 = b * block_size + 1;
        long n1 = std::min(n0 + block_size - 1, N);

        // thread-local MPFR temps
        mpfr_t n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp;
        mpfr_inits2(P, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);

        // φ(n0) = θ - t*log(n0)
        mpfr_set_ui(n_mp, (unsigned long)n0, MPFR_RNDN);
        mpfr_log(phi0, n_mp, MPFR_RNDN);
        mpfr_mul(phi0, phi0, t, MPFR_RNDN);
        mpfr_neg(phi0, phi0, MPFR_RNDN);
        mpfr_add(phi0, phi0, theta, MPFR_RNDN);

        mpfr_sin_cos(sphi, cphi, phi0, MPFR_RNDN);

        // term at n0
        mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);

        // walk inside the block
        for (long i = n0; i < n1; ++i)
        {
            // Δφ = t * log1p(1/i)
            mpfr_set_ui(n_mp, (unsigned long)i, MPFR_RNDN);
            mpfr_ui_div(dphi, 1, n_mp, MPFR_RNDN);
            mpfr_log1p(dphi, dphi, MPFR_RNDN);
            mpfr_mul(dphi, dphi, t, MPFR_RNDN);

            // rotate (cphi,sphi) by -Δφ
            mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

            mpfr_t new_c, new_s;
            mpfr_inits2(P, new_c, new_s, (mpfr_ptr)0);
            // new_c =  cphi*cd + sphi*sd
            mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
            mpfr_add(new_c, new_c, tmp, MPFR_RNDN);
            // new_s =  sphi*cd - cphi*sd
            mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
            mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
            mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

            mpfr_set(cphi, new_c, MPFR_RNDN);
            mpfr_set(sphi, new_s, MPFR_RNDN);
            mpfr_clears(new_c, new_s, (mpfr_ptr)0);

            // add term at n=i+1
            mpfr_set_ui(n_mp, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
            mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
            mpfr_add(Sloc, Sloc, tmp, MPFR_RNDN);
        }

        // accumulate to per-thread bin
        mpfr_add(partial[tid], partial[tid], Sloc, MPFR_RNDN);

        mpfr_clears(Sloc, n_mp, sqrt_n, phi0, cphi, sphi, dphi, sd, cd, tmp, (mpfr_ptr)0);
    }

    // reduce per-thread partials
    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_add(Z, Z, partial[i], MPFR_RNDN);
        mpfr_clear(partial[i]);
    }

    // Z = 2S (no half-weight)
    mpfr_mul_ui(Z, Z, 2, MPFR_RNDN);

    mpfr_clears(a, p, (mpfr_ptr)0);
}

static inline int spouge_a_for_prec(mpfr_prec_t P)
{
    // heuristic: a ≈ 0.45 * decimal_digits, cap [30, 80]
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
// log Γ(z) ≈ (z-1/2)log(z+a-1) − (z+a-1) + log( c0 + Σ c_k/(z+k−1) )
// For reflection: log Γ(z) = log π − log sin(π z) − log Γ(1−z)
// -------------------------------
static void mpc_loggamma_spouge_shift(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P)
{
    mpc_t w, s, term;
    mpc_init2(w, P);
    mpc_init2(s, P);
    mpc_init2(term, P);
    mpc_set(w, z, MPC_RNDNN);
    // w = z + (a-1)
    mpfr_add_si(mpc_realref(w), mpc_realref(w), S.a - 1, MPFR_RNDN);

    // s = c0 + Σ c_k / (z + k - 1)
    mpc_set_fr(s, &S.c[0], MPC_RNDNN);

    for (int k = 1; k < S.a; ++k)
    {
        mpc_set(term, z, MPC_RNDNN);
        mpfr_add_si(mpc_realref(term), mpc_realref(term), k - 1, MPFR_RNDN); // z + (k-1)
        mpc_ui_div(term, 1, term, MPC_RNDNN);                                // 1/(z+k-1)
        mpc_mul_fr(term, term, &S.c[k], MPC_RNDNN);                          // c_k/(...)
        mpc_add(s, s, term, MPC_RNDNN);
    }

    // (z - 1/2) * log(w) - w + log(s)
    mpc_log(out, w, MPC_RNDNN); // out = log(w)
    mpc_t zmh;
    mpc_init2(zmh, P);
    mpc_set(zmh, z, MPC_RNDNN);
    mpfr_sub_d(mpc_realref(zmh), mpc_realref(zmh), 0.5, MPFR_RNDN); // z - 1/2

    mpc_mul(out, zmh, out, MPC_RNDNN); // (z-1/2)log(w)
    mpc_sub(out, out, w, MPC_RNDNN);   // ... - w

    mpc_log(term, s, MPC_RNDNN); // log(s)
    mpc_add(out, out, term, MPC_RNDNN);

    mpc_clear(zmh);
    mpc_clear(w);
    mpc_clear(s);
    mpc_clear(term);
}

static void mpc_loggamma_spouge_reflect(mpc_t out, const mpc_t z, const Spouge &S, mpfr_prec_t P)
{
    mpc_t one_minus_z, lg1mz, sin_pi_z, t;
    mpc_init2(one_minus_z, P);
    mpc_init2(lg1mz, P);
    mpc_init2(sin_pi_z, P);
    mpc_init2(t, P);

    // 1 - z
    mpc_set_si(one_minus_z, 1, MPC_RNDNN);
    mpc_sub(one_minus_z, one_minus_z, z, MPC_RNDNN);

    // lg(1 - z)
    mpc_loggamma_spouge_shift(lg1mz, one_minus_z, S, P);

    // t = π  (as a complex with imag=0)
    mpfr_t pi;
    mpfr_init2(pi, P);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpc_set_fr(t, pi, MPC_RNDNN);
    // sin(π z)
    mpc_mul(t, t, z, MPC_RNDNN);
    mpc_sin(sin_pi_z, t, MPC_RNDNN);

    // out = log(pi) - log(sin π z) - lg(1-z)
    mpc_log(out, sin_pi_z, MPC_RNDNN);
    mpc_neg(out, out, MPC_RNDNN); // -log(sin π z)

    // + log π
    mpfr_log(pi, pi, MPFR_RNDN);
    mpc_add_fr(out, out, pi, MPC_RNDNN);

    // - lg(1-z)
    mpc_sub(out, out, lg1mz, MPC_RNDNN);

    mpfr_clear(pi);

    mpc_clear(one_minus_z);
    mpc_clear(lg1mz);
    mpc_clear(sin_pi_z);
    mpc_clear(t);
}

// θ(t) = Im(log Γ(1/4 + i t/2)) − (t/2) log π
static void compute_theta_spouge(mpfr_t theta, const mpfr_t t, mpfr_prec_t P, const Spouge &S)
{
    mpc_t z, lg;
    mpc_init2(z, P);
    mpc_init2(lg, P);

    // z = 1/4 + i t/2
    mpfr_t tr;
    mpfr_init2(tr, P);
    mpfr_div_ui(tr, t, 2, MPFR_RNDN);
    mpc_set_d_d(z, 0.25, 0.0, MPC_RNDNN);
    mpfr_set(mpc_imagref(z), tr, MPFR_RNDN);

    mpc_loggamma_spouge_reflect(lg, z, S, P);

    // theta = Im(lg) - (t/2) log π
    mpfr_set(theta, mpc_imagref(lg), MPFR_RNDN);
    mpfr_t tmp;
    mpfr_init2(tmp, P);
    mpfr_const_pi(tmp, MPFR_RNDN);
    mpfr_log(tmp, tmp, MPFR_RNDN);
    mpfr_mul(tmp, tmp, tr, MPFR_RNDN);
    mpfr_sub(theta, theta, tmp, MPFR_RNDN);

    mpfr_clear(tr);
    mpfr_clear(tmp);

    mpc_clear(z);
    mpc_clear(lg);
}

// -------------------------------
// ψ(p) = cos(2π(p^2 - p - 1/16)) / cos(2π p), with safety near 1/4, 3/4
// -------------------------------
static void psi_p(mpfr_t out, const mpfr_t p, mpfr_prec_t P)
{
    mpfr_t two_pi, num_arg, den_arg, den, num;
    mpfr_inits2(P, two_pi, num_arg, den_arg, den, num, (mpfr_ptr)0);
    mpfr_two_pi(two_pi, P);

    // den = cos(2π p)
    mpfr_mul(den_arg, two_pi, p, MPFR_RNDN);
    mpfr_cos(den, den_arg, MPFR_RNDN);

    // if |den| not tiny, do direct
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
        mpfr_mul(num_arg, two_pi, t1, MPFR_RNDN);
        mpfr_cos(num, num_arg, MPFR_RNDN);
        mpfr_div(out, num, den, MPFR_RNDN);
        mpfr_clear(t1);
    }
    else
    {
        // near removable singularity: linearized fallbacks
        mpfr_t q, half;
        mpfr_inits2(P, q, half, (mpfr_ptr)0);
        mpfr_set_d(half, 0.5, MPFR_RNDN);
        // choose 1/4 or 3/4
        mpfr_t d1, d2;
        mpfr_inits2(P, d1, d2, (mpfr_ptr)0);
        mpfr_sub_d(d1, p, 0.25, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_sub_d(d2, p, 0.75, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        if (mpfr_cmp(d1, d2) < 0)
        {
            // ~ 1/2 - (p-1/4)
            mpfr_sub_d(q, p, 0.25, MPFR_RNDN);
            mpfr_sub(out, half, q, MPFR_RNDN);
        }
        else
        {
            // ~ 1/2 + (p-3/4)
            mpfr_sub_d(q, p, 0.75, MPFR_RNDN);
            mpfr_add(out, half, q, MPFR_RNDN);
        }
        mpfr_clears(q, half, d1, d2, (mpfr_ptr)0);
    }

    mpfr_clear(eps);

    mpfr_clears(two_pi, num_arg, den_arg, den, num, (mpfr_ptr)0);
}



// central derivatives of ψ at p, 9-point, step h; returns ψ'', ψ''', ψ''''''.
static void psi_derivs_central(mpfr_t d2, mpfr_t d3, mpfr_t d6,
                               const mpfr_t p, mpfr_prec_t P,
                               const DiffWeights &W, const mpfr_t h)
{
    const int r = W.r;
    std::vector<mpfr_t> f(2 * r + 1);
    for (auto &v : f)
        mpfr_init2(v, P);

    // f_j = ψ(p + (j-r) h)
    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t arg;
        mpfr_init2(arg, P);
        mpfr_set(arg, p, MPFR_RNDN);
        mpfr_t step;
        mpfr_init2(step, P);
        mpfr_mul_si(step, h, j - r, MPFR_RNDN);
        mpfr_add(arg, arg, step, MPFR_RNDN);
        psi_p(f[j], arg, P);
        mpfr_clears(arg, step, (mpfr_ptr)0);
    }
    // dot products
    mpfr_set_zero(d2, 0);
    mpfr_set_zero(d3, 0);
    mpfr_set_zero(d6, 0);
    for (int j = 0; j < 2 * r + 1; j++)
    {
        mpfr_t t;
        mpfr_init2(t, P);
        // ψ''
        mpfr_mul(t, &W.w2[j], f[j], MPFR_RNDN);
        mpfr_add(d2, d2, t, MPFR_RNDN);
        // ψ'''
        mpfr_mul(t, &W.w3[j], f[j], MPFR_RNDN);
        mpfr_add(d3, d3, t, MPFR_RNDN);
        // ψ''''''
        mpfr_mul(t, &W.w6[j], f[j], MPFR_RNDN);
        mpfr_add(d6, d6, t, MPFR_RNDN);
        mpfr_clear(t);
    }
    // divide by h^m
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
                     DiffWeights *W_opt = nullptr)
{
    mpfr_t two_pi, a, p;
    mpfr_inits2(P, two_pi, a, p, (mpfr_ptr)0);
    mpfr_two_pi(two_pi, P);
    // a, p, N
    long N = rs_N_from_t(t, P, a, p);

    // A0 = (2π/t)^(1/4)
    mpfr_t A0;
    mpfr_init2(A0, P);
    mpfr_div(A0, two_pi, t, MPFR_RNDN);
    mpfr_rootn_ui(A0, A0, 4, MPFR_RNDN);

    // C0
    mpfr_t C0;
    mpfr_init2(C0, P);
    psi_p(C0, p, P);

    mpfr_set(R, C0, MPFR_RNDN);

    if (useC1C2)
    {
        // step size h
        mpfr_t d1, quarter, threeq;
        mpfr_inits2(P, d1, quarter, threeq, (mpfr_ptr)0);
        mpfr_set_d(quarter, 0.25, MPFR_RNDN);
        mpfr_set_d(threeq, 0.75, MPFR_RNDN);
        mpfr_sub(d1, p, quarter, MPFR_RNDN);
        mpfr_abs(d1, d1, MPFR_RNDN);
        mpfr_t d2;
        mpfr_init2(d2, P);
        mpfr_sub(d2, p, threeq, MPFR_RNDN);
        mpfr_abs(d2, d2, MPFR_RNDN);
        mpfr_t h;
        mpfr_init2(h, P);
        if (mpfr_cmp(d1, d2) < 0)
        {
            // near 1/4
            mpfr_mul_d(h, d1, 0.1, MPFR_RNDN);
        }
        else
        {
            mpfr_mul_d(h, d2, 0.1, MPFR_RNDN);
        }
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
        psi_derivs_central(psi2, psi3, psi6, p, P, *W, h);

        // C1 = -(1/(96π^2)) ψ'''(p)
        mpfr_t C1;
        mpfr_init2(C1, P);
        mpfr_const_pi(C1, MPFR_RNDN);
        mpfr_mul(C1, C1, C1, MPFR_RNDN);    // π^2
        mpfr_mul_ui(C1, C1, 96, MPFR_RNDN); // 96π^2
        mpfr_div(C1, psi3, C1, MPFR_RNDN);
        mpfr_neg(C1, C1, MPFR_RNDN);

        // C2 = ψ^{(6)}/(18432π^4) + ψ''/(64π^2)
        mpfr_t C2, tmp;
        mpfr_inits2(P, C2, tmp, (mpfr_ptr)0);
        mpfr_const_pi(tmp, MPFR_RNDN);      // tmp = π
        mpfr_mul(tmp, tmp, tmp, MPFR_RNDN); // π^2
        mpfr_t pi2;
        mpfr_init2(pi2, P);
        mpfr_set(pi2, tmp, MPFR_RNDN);
        mpfr_mul(tmp, tmp, tmp, MPFR_RNDN); // π^4
        // ψ''/(64π^2)
        mpfr_t term2;
        mpfr_init2(term2, P);
        mpfr_mul_ui(term2, pi2, 64, MPFR_RNDN);
        mpfr_div(term2, psi2, term2, MPFR_RNDN);
        // ψ^{(6)}/(18432 π^4)
        mpfr_t term1;
        mpfr_init2(term1, P);
        mpfr_mul_ui(term1, tmp, 18432, MPFR_RNDN);
        mpfr_div(term1, psi6, term1, MPFR_RNDN);
        mpfr_add(C2, term1, term2, MPFR_RNDN);

        // scale factors (2π/t)^(1/2) and (2π/t)
        mpfr_t s12, s1;
        mpfr_inits2(P, s12, s1, (mpfr_ptr)0);
        mpfr_div(s1, two_pi, t, MPFR_RNDN); // (2π)/t
        mpfr_sqrt(s12, s1, MPFR_RNDN);      // ^1/2
        // R += C1*s12 + C2*s1
        mpfr_t add;
        mpfr_init2(add, P);
        mpfr_mul(add, C1, s12, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);
        mpfr_mul(add, C2, s1, MPFR_RNDN);
        mpfr_add(R, R, add, MPFR_RNDN);

        mpfr_clears(psi2, psi3, psi6, C1, C2, tmp, pi2, term1, term2, s12, s1, add, (mpfr_ptr)0);
        mpfr_clears(d1, d2, quarter, threeq, h, (mpfr_ptr)0);
    }

    // R *= A0
    mpfr_mul(R, R, A0, MPFR_RNDN);

    // sign (-1)^(N-1)
    if ((N % 2) == 0)
        mpfr_neg(R, R, MPFR_RNDN);

    mpfr_clears(two_pi, a, p, A0, C0, (mpfr_ptr)0);
}

// -------------------------------
// Main RS sum via rotations, NO half-weight
// Z = 2 Σ_{n=1}^N cos(θ - t log n)/√n
// with stable rotation using Δφ = t log(1 + 1/i)
// -------------------------------
static void Z_rs_main(mpfr_t Z, const mpfr_t t, const mpfr_t theta, mpfr_prec_t P)
{

    mpfr_t a, p;
    mpfr_inits2(P, a, p, (mpfr_ptr)0);
    long N = rs_N_from_t(t, P, a, p);

    mpfr_t S;
    mpfr_init2(S, P);
    mpfr_set_zero(S, 0);

    // seed n=1: φ=θ
    mpfr_t cphi, sphi;
    mpfr_inits2(P, cphi, sphi, (mpfr_ptr)0);
    mpfr_sin_cos(sphi, cphi, theta, MPFR_RNDN);
    // n=1 term
    mpfr_add(S, S, cphi, MPFR_RNDN);

    mpfr_t dphi, sd, cd, new_c, new_s, tmp;
    mpfr_inits2(P, dphi, sd, cd, new_c, new_s, tmp, (mpfr_ptr)0);
    mpfr_t i_mp, sqrt_n;
    mpfr_inits2(P, i_mp, sqrt_n, (mpfr_ptr)0);

    for (long i = 1; i < N; ++i)
    {
        // Δφ = t*log1p(1/i)
        mpfr_set_ui(i_mp, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(dphi, 1, i_mp, MPFR_RNDN);
        mpfr_log1p(dphi, dphi, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

        // rotate by -Δφ
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        // add term at n=i+1
        mpfr_set_ui(i_mp, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_sqrt(sqrt_n, i_mp, MPFR_RNDN);
        mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(S, S, tmp, MPFR_RNDN);
    }

    // Z = 2S (no half-weight)
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
    bool parallel)
{

    mpfr_t theta;
    mpfr_init2(theta, P);

    compute_theta_spouge(theta, t, P, S);

    mpfr_t Zmain;
    mpfr_init2(Zmain, P);

    // parallel optimization
    if (parallel)
        // threads=0 → uses omp_get_max_threads().
        // Tune block_size (4k–32k). Larger blocks = fewer log(n0) calls, smaller overhead; smaller blocks = more parallelism.
        Z_rs_main_parallel(Zmain, t, theta, P, /*threads*/ 0, /*block_size*/ 8192);
    else
        // single thread
        Z_rs_main(Zmain, t, theta, P);

    mpfr_t R1;
    mpfr_init2(R1, P);
    R_gabcke(R1, t, P, useC1C2, W_opt);

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

    // Examples for calculations
    // tvals = {
    //     "3047.70881764088986878223663882",
    //     "5045.94167550201313249597140915",
    //     "7040.82831205880106241115667351",
    //     "88888891.29471578310681783197410193108012386",
    //     "888888883.41171178275288939532420360893015575",
    //     "10000"};

    mpfr_set_default_prec(P);

    //oss << "mpfr_prec = " << P << " bits (" << dec_digits_from_bits(P) << " d.p.)\n";
    //LOGP(oss.str(), CFG().rs_log_digits);    

    // prepare Spouge coeffs
    Spouge S;
    int a = spouge_a_for_prec(P);
    S.build(a, P);

    // optional derivative weights once
    DiffWeights W; // build on first use if useC1C2=true

    for (auto &ts : tvals)
    {
        mpfr_t tval;
        mpfr_init2(tval, P);
        //mpfr_set_str(tval, ts.c_str(), 10, MPFR_RNDN);
        mpfr_set(tval, ts, MPFR_RNDN);

        mpfr_t Ztot, Zmain, R1;
        mpfr_inits2(P, Ztot, Zmain, R1, (mpfr_ptr)0);

        bool useC1C2 = true;
        Z_rs_with_gabcke(Ztot, Zmain, R1, tval, P, useC1C2, &W, S, true);

        oss << "t = " << mp(tval, CFG().rs_log_digits, 'g') << " " << "Ztot = " << mp(Ztot, CFG().rs_log_digits, 'g') << "\n";   // or 'e' / 'f'
        LOGP(oss.str(), CFG().rs_log_digits);
        //std::cout << "a(spouge) = " << a << "\n";
        //print_mpfr("Z_main", Zmain, 30);
        //print_mpfr("R1", R1, 30);
        //print_mpfr("Z_tot", Ztot, 30);

        // last Z calculated from the array is returned, this loop while helps to test many t
        mpfr_set(Z, Ztot, MPFR_RNDN);

        mpfr_clears(tval, Ztot, Zmain, R1, (mpfr_ptr)0);
    }

    tvals.clear();
}