#include "z_theta.h"
#include "logs.h"
#include "settings.h"
#include "z_calculator.h"
#include <mpfr.h>
#include <mpc.h>
#include <thread>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <atomic>

struct WorkerArgs
{
    mpfr_srcptr t;
    mpfr_srcptr theta;
    long i0, i1;
    int id;
    mpfr_t partial;
};

static void pairwise_add_inplace(std::vector<mpfr_ptr> &v, mpfr_t out)
{
    if (v.empty())
    {
        mpfr_set_zero(out, 0);
        return;
    }
    size_t n = v.size();
    while (n > 1)
    {
        size_t m = 0;
        for (size_t i = 0; i + 1 < n; i += 2)
        {
            mpfr_add(v[i], v[i], v[i + 1], MPFR_RNDN);
            v[m++] = v[i];
        }
        if (n & 1)
            v[m++] = v[n - 1];
        n = m;
    }
    mpfr_add(out, out, v[0], MPFR_RNDN);
}

static inline void rs_block_sum(
    mpfr_t out, mpfr_srcptr t, mpfr_srcptr theta, long i0, long i1)
{
    // out += sum_{i=i0..i1} cos(theta - t*log i)/sqrt(i)
    mpfr_t n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, new_c, new_s;
    mpfr_inits2(CFG().mpfr_prec,
                n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, new_c, new_s, (mpfr_ptr)0);

    // Seed at i0
    mpfr_set_ui(n, (unsigned long)i0, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, t, ln_n, MPFR_RNDN);
    mpfr_sub(phi, theta, phi, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif

    for (long i = i0; i <= i1; ++i)
    {
        // add cos(phi)/sqrt(i)
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
        mpfr_div(term, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(out, out, term, MPFR_RNDN);

        if (i == i1)
            break;

        // advance phase: Δφ = t*log1p(1/i)
        mpfr_ui_div(q, 1, n, MPFR_RNDN);
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate by -Δφ
        // new_c = c*cd + s*sd
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);
        // new_s = s*cd − c*sd
        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);
    }

    mpfr_clears(n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, new_c, new_s, (mpfr_ptr)0);
}

// configurable remainder bound
void remainder_bound(mpfr_t R, const mpfr_t t)
{
    const double C = CFG().rs_tail_C; // <--- configurable
    mpfr_t r, tmp;
    mpfr_inits2(CFG().mpfr_prec, r, tmp, (mpfr_ptr)0);

    mpfr_sqrt(tmp, t, MPFR_RNDN); // √t
    mpfr_sqrt(r, tmp, MPFR_RNDN); // t^(1/4)
    mpfr_ui_div(R, 1, r, MPFR_RNDN);
    mpfr_mul_d(R, R, C, MPFR_RNDN);

    mpfr_clears(r, tmp, (mpfr_ptr)0);
}

static void rs_worker(WorkerArgs *A)
{
    mpfr_t tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp;
    mpfr_inits2(CFG().mpfr_prec, tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);

    mpfr_set(tloc, A->t, MPFR_RNDN);
    mpfr_set(thetaloc, A->theta, MPFR_RNDN);

    mpfr_set_ui(n, (unsigned long)A->i0, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, tloc, ln_n, MPFR_RNDN);
    mpfr_sub(phi, thetaloc, phi, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif
    mpfr_set_zero(A->partial, 0);

    for (long i = A->i0; i <= A->i1; ++i)
    {
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
        mpfr_div(term, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(A->partial, A->partial, term, MPFR_RNDN);

        if (i == A->i1)
            break;

        mpfr_ui_div(q, 1, n, MPFR_RNDN);
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, tloc, MPFR_RNDN);
#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR == 4 && MPFR_VERSION_MINOR >= 1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        mpfr_mul(tmp, cphi, cd, MPFR_RNDN);
        mpfr_mul(term, sphi, sd, MPFR_RNDN);
        mpfr_add(phi, tmp, term, MPFR_RNDN); // new_c
        mpfr_mul(tmp, sphi, cd, MPFR_RNDN);
        mpfr_mul(term, cphi, sd, MPFR_RNDN);
        mpfr_sub(sphi, tmp, term, MPFR_RNDN); // new_s
        mpfr_set(cphi, phi, MPFR_RNDN);
    }

    mpfr_clears(tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);
}

// heuristic: pick a block size that gives ~8 blocks per thread (cap at 256 min)
static inline long choose_block(long N, unsigned nth)
{
    long tgt = (long)std::max(8u, 8u); // stylistic
    long blk = N / std::max<long>((long)nth * tgt, 1);

    if (blk < 256)
        blk = std::min<long>(N, 256);

    return blk;
}

void compute_Z_threaded(mpfr_t Z, mpfr_t t)
{
    // --- N and θ
    long N = rs_N(t);

    mpfr_t theta;
    mpfr_init2(theta, CFG().mpfr_prec);

    compute_theta_stirling(theta, t);

    if (CFG().rs_log_debug)
    {
        log_mpfr("[RS] t", t, CFG().rs_log_digits);
        log_mpfr("[RS] theta", theta, CFG().rs_log_digits);
        std::cerr << "[RS] N=" << N << "\n";
    }

    // --- thread count
    unsigned hw = std::thread::hardware_concurrency();
    unsigned maxth = (CFG().max_threads > 0) ? (unsigned)CFG().max_threads : (hw ? hw : 1u);
    unsigned nth = std::clamp<unsigned>(hw ? hw : 1u, CFG().min_threads, maxth);
    if (N < CFG().rs_min_N_for_threads)
        nth = 1u;
    if ((int)nth > N)
        nth = (unsigned)N;

    // --- dynamic block scheduler
    std::atomic<long> next(1);
    const long block = choose_block(N, nth);

    if (CFG().rs_log_debug)
    {
        std::cerr << "[RS] threads=" << nth << " block=" << block << "\n";
    }

    std::vector<std::thread> pool;
    pool.reserve(nth);

    // each thread keeps its own partial
    std::vector<mpfr_t> partials(nth);
    for (unsigned k = 0; k < nth; ++k)
    {
        mpfr_init2(partials[k], CFG().mpfr_prec);
        mpfr_set_zero(partials[k], 0);
    }

    auto worker = [&](unsigned tid)
    {
        mpfr_t local;
        mpfr_init2(local, CFG().mpfr_prec);
        mpfr_set_zero(local, 0);

        for (;;)
        {
            long i0 = next.fetch_add(block, std::memory_order_relaxed);
            if (i0 > N)
                break;
            long i1 = i0 + block - 1;
            if (i1 > N)
                i1 = N;

            rs_block_sum(local, t, theta, i0, i1);
        }

        mpfr_add(partials[tid], partials[tid], local, MPFR_RNDN);
        mpfr_clear(local);
    };

    for (unsigned k = 0; k < nth; ++k)
        pool.emplace_back(worker, k);

    for (auto &th : pool)
        th.join();

    // --- reduction
    mpfr_t S;
    mpfr_init2(S, CFG().mpfr_prec);
    mpfr_set_zero(S, 0);

    std::vector<mpfr_ptr> parts;
    parts.reserve(nth);

    for (unsigned k = 0; k < nth; ++k)
        parts.push_back(partials[k]);

    pairwise_add_inplace(parts, S);

    for (unsigned k = 0; k < nth; ++k)
        mpfr_clear(partials[k]);

    // Z = 2*S (+ optional ½-weight last term)
    mpfr_mul_ui(Z, S, 2, MPFR_RNDN);

    if (CFG().rs_half_weight_last && N >= 1)
    {
        mpfr_t nN, lnN, phiN, cphiN, sqrtN, termN;
        mpfr_inits2(CFG().mpfr_prec, nN, lnN, phiN, cphiN, sqrtN, termN, (mpfr_ptr)0);
        mpfr_set_ui(nN, (unsigned long)N, MPFR_RNDN);
        mpfr_log(lnN, nN, MPFR_RNDN);
        mpfr_mul(phiN, t, lnN, MPFR_RNDN);
        mpfr_sub(phiN, theta, phiN, MPFR_RNDN);
        mpfr_cos(cphiN, phiN, MPFR_RNDN);
        mpfr_sqrt(sqrtN, nN, MPFR_RNDN);
        mpfr_div(termN, cphiN, sqrtN, MPFR_RNDN);
        mpfr_sub(Z, Z, termN, MPFR_RNDN);
        mpfr_clears(nN, lnN, phiN, cphiN, sqrtN, termN, (mpfr_ptr)0);
    }

    // remainder (diagnostic)
    if (CFG().rs_log_debug)
    {
        mpfr_t R;
        mpfr_init2(R, CFG().mpfr_prec);
        remainder_bound(R, t);
        log_mpfr("[RS] |R_bound|", R, CFG().rs_log_digits);
        log_mpfr("[RS] Z", Z, CFG().rs_log_digits);
        mpfr_clear(R);
    }

    mpfr_clears(theta, S, (mpfr_ptr)0);
}
