#include <mpfr.h>
#include <mpc.h>
#include <thread>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <atomic>
#include "config.h"


struct WorkerArgs {
    mpfr_srcptr t;      // shared read-only
    mpfr_srcptr theta;  // shared read-only
    long i0, i1;
    mpfr_prec_t prec;
    int id;
    mpfr_t partial;     // output
};


static inline void log_mpfr(
    const char* tag, 
    const mpfr_t x, 
    int digits=30) {

    std::cerr << tag << " = ";
    mpfr_out_str(stderr, 10, digits, x, MPFR_RNDN);
    std::cerr << "\n";

}

// ---------------- Theta(t) ----------------
// theta(t) = (t/2) log(t/2π) - t/2 - π/8 + 1/(48 t) + O(t^-3)
void compute_theta(
        mpfr_t theta, 
        const mpfr_t t) {
    
    mpfr_t pi, two_pi, x;

    mpfr_inits2(MPFR_PRECISION, pi, two_pi, x, (mpfr_ptr)0);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    mpfr_div(x, t, two_pi, MPFR_RNDN);      // t/(2π)
    mpfr_log(x, x, MPFR_RNDN);              // log(t/2π)
    mpfr_mul(theta, x, t, MPFR_RNDN);       // t*log(...)
    mpfr_div_ui(theta, theta, 2, MPFR_RNDN);// (t/2) log(...)

    mpfr_div_ui(x, t, 2, MPFR_RNDN);        // t/2
    mpfr_sub(theta, theta, x, MPFR_RNDN);   // ... - t/2

    mpfr_div_ui(x, pi, 8, MPFR_RNDN);       // π/8
    mpfr_sub(theta, theta, x, MPFR_RNDN);   // ... - π/8

    mpfr_mul_ui(x, t, 48, MPFR_RNDN);       // 48 t
    mpfr_ui_div(x, 1, x, MPFR_RNDN);        // 1/(48 t)
    mpfr_add(theta, theta, x, MPFR_RNDN);   // ... + 1/(48 t)

    mpfr_clears(pi, two_pi, x, (mpfr_ptr)0);
}

// ---------------- Remainder bound ----------------
// Safe bound R(t) ≲ C * t^(-1/4). Tweak C if you later add sharper corrections.
static inline void remainder_bound(mpfr_t R, const mpfr_t t, double C = 0.05) {
    mpfr_t r, tmp;
    
    mpfr_inits2(MPFR_PRECISION, r, tmp, (mpfr_ptr)0);
    mpfr_sqrt(tmp, t, MPFR_RNDN);   // √t
    mpfr_sqrt(r, tmp, MPFR_RNDN);   // t^(1/4)
    mpfr_ui_div(R, 1, r, MPFR_RNDN);
    mpfr_mul_d(R, R, C, MPFR_RNDN);
    mpfr_clears(r, tmp, (mpfr_ptr)0);
}

// ---------------- Prefilter ----------------
// Compute partial RS sum up to K=⌊αN⌋, bound the tail + remainder.
// Discard t if |2 S_K| - 2*T_K - R > τ (i.e., cannot be near zero).
bool z_prefilter(const mpfr_t t, double alpha, double tau) {
    
    if (alpha <= 0.0) alpha = 0.25;
    if (alpha >= 1.0) alpha = 0.8;

    mpfr_t pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, (mpfr_ptr)0);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(Nmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(Nmp, Nmp, MPFR_RNDN);
    long N = mpfr_get_si(Nmp, MPFR_RNDN);
    if (N < 1) N = 1;

    long K = (long)std::floor(alpha * N);
    if (K < 4) K = (N >= 4) ? 4 : N;
    if (K > N) K = N;

    compute_theta(theta, t);

    // Start at n=1
    mpfr_set_zero(S, 0);
    mpfr_set_ui(n, 1, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, t, ln_n, MPFR_RNDN);
    mpfr_sub(phi, theta, phi, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif
    mpfr_add(S, S, cphi, MPFR_RNDN); // term at n=1

    for (long i = 1; i < K; ++i) {
        // Δφ = t * log1p(1/i)
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_ui_div(q, 1, n, MPFR_RNDN);          // 1/i
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, t, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate (φ -> φ - Δφ)
        mpfr_t new_c, new_s, sqrt_i1, i1, tmp;
        mpfr_inits2(MPFR_PRECISION, new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);

        // new_c = c*cd + s*sd
        mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  sphi, sd, MPFR_RNDN);
        mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

        // new_s = s*cd - c*sd
        mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
        mpfr_mul(tmp,  cphi, sd, MPFR_RNDN);
        mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

        mpfr_set(cphi, new_c, MPFR_RNDN);
        mpfr_set(sphi, new_s, MPFR_RNDN);

        mpfr_set_ui(i1, (unsigned long)(i+1), MPFR_RNDN);
        mpfr_sqrt(sqrt_i1, i1, MPFR_RNDN);
        mpfr_div(new_c, cphi, sqrt_i1, MPFR_RNDN);
        mpfr_add(S, S, new_c, MPFR_RNDN);

        mpfr_clears(new_c, new_s, sqrt_i1, i1, tmp, (mpfr_ptr)0);
    }

    mpfr_t twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau;
    mpfr_inits2(MPFR_PRECISION, twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau, (mpfr_ptr)0);

    mpfr_mul_ui(twoS, S, 2, MPFR_RNDN);
    mpfr_abs(abs_twoS, twoS, MPFR_RNDN);

    mpfr_set_ui(Nmp, (unsigned long)N, MPFR_RNDN);
    mpfr_set_ui(Kmp, (unsigned long)K, MPFR_RNDN);
    mpfr_sqrt(sqrtN, Nmp, MPFR_RNDN);
    mpfr_sqrt(sqrtK, Kmp, MPFR_RNDN);
    mpfr_sub(tail, sqrtN, sqrtK, MPFR_RNDN);
    mpfr_mul_ui(tail, tail, 2, MPFR_RNDN); // T_K
    mpfr_mul_ui(tail, tail, 2, MPFR_RNDN); // 2*T_K for Z

    remainder_bound(R, t);

    // L = |2S_K| - 2*T_K - R
    mpfr_sub(L, abs_twoS, tail, MPFR_RNDN);
    mpfr_sub(L, L, R, MPFR_RNDN);

    mpfr_set_d(Tau, tau, MPFR_RNDN);
    bool keep = (mpfr_cmp(L, Tau) <= 0);

    std::cerr << "[prefilter] N=" << N << " K=" << K << " alpha=" << alpha << " tau=" << tau << "\n";
    log_mpfr("[prefilter] |2S_K|", abs_twoS, 18);
    log_mpfr("[prefilter] 2*T_K", tail, 18);
    log_mpfr("[prefilter] R", R, 18);
    log_mpfr("[prefilter] L", L, 18);
    std::cerr << "[prefilter] keep=" << (keep ? "YES" : "NO") << "\n";

    mpfr_clears(pi, two_pi, Nmp, Kmp, theta, S, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd,
                twoS, abs_twoS, tail, sqrtN, sqrtK, R, L, Tau, (mpfr_ptr)0);
    return keep;
}

static void rs_worker(WorkerArgs* A) {
    const long int prec = A->prec;

    mpfr_t tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp;
    mpfr_inits2(MPFR_PRECISION, tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);

    mpfr_set(tloc,    A->t,     MPFR_RNDN);
    mpfr_set(thetaloc,A->theta, MPFR_RNDN);

    std::cerr << "[worker " << A->id << "] range [" << A->i0 << "," << A->i1 << "]\n";

    mpfr_set_ui(n, (unsigned long)A->i0, MPFR_RNDN);
    mpfr_log(ln_n, n, MPFR_RNDN);
    mpfr_mul(phi, tloc, ln_n, MPFR_RNDN);
    mpfr_sub(phi, thetaloc, phi, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
    mpfr_sin_cos(sphi, cphi, phi, MPFR_RNDN);
#else
    mpfr_cos(cphi, phi, MPFR_RNDN);
    mpfr_sin(sphi, phi, MPFR_RNDN);
#endif

    mpfr_set_zero(A->partial, 0);

    for (long i = A->i0; i <= A->i1; ++i) {
        mpfr_set_ui(n, (unsigned long)i, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n, MPFR_RNDN);
        mpfr_div(term, cphi, sqrt_n, MPFR_RNDN);
        mpfr_add(A->partial, A->partial, term, MPFR_RNDN);

        if (i == A->i0 || i == A->i0+1 || i == A->i1-1 || i == A->i1) {
            std::cerr << "[worker " << A->id << "] i=" << i << " term=";
            mpfr_out_str(stderr, 10, 18, term, MPFR_RNDN);
            std::cerr << "  cphi=";
            mpfr_out_str(stderr, 10, 12, cphi, MPFR_RNDN);
            std::cerr << "\n";
        }

        if (i == A->i1) break;

        // Δφ = t * log1p(1/i)
        mpfr_ui_div(q, 1, n, MPFR_RNDN);
        mpfr_log1p(dphi, q, MPFR_RNDN);
        mpfr_mul(dphi, dphi, tloc, MPFR_RNDN);

#if MPFR_VERSION_MAJOR > 4 || (MPFR_VERSION_MAJOR==4 && MPFR_VERSION_MINOR>=1)
        mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);
#else
        mpfr_cos(cd, dphi, MPFR_RNDN);
        mpfr_sin(sd, dphi, MPFR_RNDN);
#endif
        // rotate: new_c = c*cd + s*sd ; new_s = s*cd - c*sd
        mpfr_mul(tmp, cphi, cd, MPFR_RNDN);
        mpfr_mul(term, sphi, sd, MPFR_RNDN);
        mpfr_add(phi, tmp, term, MPFR_RNDN);      // phi holds new_c

        mpfr_mul(tmp, sphi, cd, MPFR_RNDN);
        mpfr_mul(term, cphi, sd, MPFR_RNDN);
        mpfr_sub(sphi, tmp, term, MPFR_RNDN);     // sphi becomes new_s

        mpfr_set(cphi, phi, MPFR_RNDN);
    }

    std::cerr << "[worker " << A->id << "] partial=";
    mpfr_out_str(stderr, 10, 18, A->partial, MPFR_RNDN);
    std::cerr << "\n";

    mpfr_clears(tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);
}

void compute_Z_threaded(mpfr_t Z, mpfr_t t) {
    
    unsigned int nthreads = std::thread::hardware_concurrency();

    mpfr_t theta, pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, theta, pi, two_pi, tmp, (mpfr_ptr)0);
    compute_theta(theta, t);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    if (N < 1) N = 1;

    std::cerr << "[main] threads=" << nthreads << "  N=" << N << "\n";
    log_mpfr("[main] t", t, 20);
    log_mpfr("[main] theta", theta, 20);

    // If N is tiny, RS truncation is not balanced; caller should use small-t path.
    if ((int)nthreads > N) nthreads = (unsigned)N;

    std::vector<WorkerArgs> args(nthreads);
    std::vector<std::thread> pool; pool.reserve(nthreads);

    long base = N / (long)nthreads, rem = N % (long)nthreads;
    long cur = 1;
    for (unsigned k = 0; k < nthreads; ++k) {
        long len = base + ((long)k < rem ? 1 : 0);
        args[k].t     = t;
        args[k].theta = theta;
        args[k].i0    = cur;
        args[k].i1    = cur + len - 1;
        args[k].prec  = MPFR_PRECISION;
        args[k].id    = (int)k;
        mpfr_init2(args[k].partial, MPFR_PRECISION);
        mpfr_set_zero(args[k].partial, 0);

        std::cerr << "[main] chunk " << k << ": [" << args[k].i0 << "," << args[k].i1 << "]\n";

        cur += len;
    }

    for (unsigned k = 0; k < nthreads; ++k) pool.emplace_back(rs_worker, &args[k]);
    for (auto& th : pool) th.join();

    mpfr_t sum; mpfr_init2(sum, MPFR_PRECISION); mpfr_set_zero(sum, 0);
    for (unsigned k = 0; k < nthreads; ++k) {
        mpfr_add(sum, sum, args[k].partial, MPFR_RNDN);
        mpfr_clear(args[k].partial);
    }
    mpfr_mul_ui(sum, sum, 2, MPFR_RNDN); // 2 * Σ

    mpfr_t R; mpfr_init2(R, MPFR_PRECISION);
    remainder_bound(R, t);
    if ((N - 1) & 1) mpfr_neg(R, R, MPFR_RNDN);

    mpfr_add(Z, sum, R, MPFR_RNDN);

    log_mpfr("[main] 2*sum", sum, 20);
    log_mpfr("[main] R", R, 20);
    log_mpfr("[main] Z", Z, 20);
    if (mpfr_nan_p(Z))  std::cerr << "[main] Z is NaN!\n";
    if (mpfr_inf_p(Z))  std::cerr << "[main] Z is Inf!\n";

    mpfr_clears(theta, pi, two_pi, tmp, sum, R, (mpfr_ptr)0);
}



