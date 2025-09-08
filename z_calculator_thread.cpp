#include <mpfr.h>
#include <mpc.h>
#include <thread>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <atomic>
#include "config.h"
#include "z_theta.h"
#include "logs.h"

struct WorkerArgs {
    mpfr_srcptr t;      // shared read-only
    mpfr_srcptr theta;  // shared read-only
    long i0, i1;
    int id;
    mpfr_t partial;     // output sum
};

// In-place pairwise reducer: writes partial sums back into v[0..m).
// In-place pairwise reducer operating on pointers.
// v holds mpfr_ptr (pointers to mpfr_t), not mpfr_t by value.
static void pairwise_add_inplace(std::vector<mpfr_ptr>& v, mpfr_t out) {
    if (v.empty()) { mpfr_set_zero(out, 0); return; }

    size_t n = v.size();
    while (n > 1) {
        size_t m = 0;
        for (size_t i = 0; i + 1 < n; i += 2) {
            // accumulate v[i+1] into v[i] (in-place)
            mpfr_add(v[i], v[i], v[i+1], MPFR_RNDN);
            v[m++] = v[i];            // keep pointer to the updated accumulator
        }
        if (n & 1) v[m++] = v[n - 1]; // carry the last (odd) pointer
        n = m;
    }
    // out += total
    mpfr_add(out, out, v[0], MPFR_RNDN);
}

// ---------------- Remainder bound ----------------
// Safe bound R(t) ≲ C * t^(-1/4). Tweak C if you later add sharper corrections.
 void remainder_bound(
        mpfr_t R, 
        const mpfr_t t) {
    
    double C = 0.05;
    mpfr_t r, tmp;
    
    mpfr_inits2(MPFR_PRECISION, r, tmp, (mpfr_ptr)0);
    mpfr_sqrt(tmp, t, MPFR_RNDN);   // √t
    mpfr_sqrt(r, tmp, MPFR_RNDN);   // t^(1/4)
    mpfr_ui_div(R, 1, r, MPFR_RNDN);
    mpfr_mul_d(R, R, C, MPFR_RNDN);
    mpfr_clears(r, tmp, (mpfr_ptr)0);
}

static void rs_worker(WorkerArgs* A) {

    mpfr_t tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp;
    mpfr_inits2(MPFR_PRECISION, tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);

    mpfr_set(tloc,    A->t,     MPFR_RNDN);
    mpfr_set(thetaloc,A->theta, MPFR_RNDN);

    //std::cerr << "[worker " << A->id << "] range [" << A->i0 << "," << A->i1 << "]\n";

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
            //std::cerr << "[worker " << A->id << "] i=" << i << " term=";
            //mpfr_out_str(stderr, 10, 18, term, MPFR_RNDN);
            //std::cerr << "  cphi=";
            //mpfr_out_str(stderr, 10, 12, cphi, MPFR_RNDN);
            //std::cerr << "\n";
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

    //std::cerr << "[worker " << A->id << "] partial=";
    //mpfr_out_str(stderr, 10, 18, A->partial, MPFR_RNDN);
    //std::cerr << "\n";

    mpfr_clears(tloc, thetaloc, n, ln_n, phi, cphi, sphi, q, dphi, cd, sd, sqrt_n, term, tmp, (mpfr_ptr)0);
}

void compute_Z_threaded(
        mpfr_t Z, 
        mpfr_t t) {
    
    unsigned int nthreads = std::thread::hardware_concurrency();

    mpfr_t theta, pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, theta, pi, two_pi, tmp, (mpfr_ptr)0);
    
    //compute_theta(theta, t);
    compute_theta_big(theta, t);

    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    
    if (N < 1) N = 1;
    
    //std::cerr << "[main] threads=" << nthreads << "  N=" << N << "\n";
    log_mpfr("[main] t", t, 20);
    log_mpfr("[main] theta", theta, 20);

    // If N is tiny, RS truncation is not balanced; caller should use small-t path.
    if ((int)nthreads > N) {
        nthreads = (unsigned)N;
    }

    std::vector<WorkerArgs> args(nthreads);
    std::vector<std::thread> pool; 
    pool.reserve(nthreads);

    long base = N / (long)nthreads, rem = N % (long)nthreads;
    long cur = 1;
    
    for (unsigned k = 0; k < nthreads; ++k) {
        long len = base + ((long)k < rem ? 1 : 0);
        args[k].t     = t;
        args[k].theta = theta;
        args[k].i0    = cur;
        args[k].i1    = cur + len - 1;
        args[k].id    = (int)k;

        mpfr_init2(args[k].partial, MPFR_PRECISION);
        mpfr_set_zero(args[k].partial, 0);

        //std::cerr << "[main] chunk " << k << ": [" << args[k].i0 << "," << args[k].i1 << "]\n";

        cur += len;
    }

    for (unsigned k = 0; k < nthreads; ++k) {
        pool.emplace_back(rs_worker, &args[k]);
    }
    
    for (auto& th : pool) {
        th.join();
    }

    // mpfr_t sum; 
    // mpfr_init2(sum, MPFR_PRECISION); 
    // mpfr_set_zero(sum, 0);

    // std::vector<mpfr_t> parts; parts.reserve(nthreads);
    // for (unsigned k=0; k<nthreads; ++k) { 
    //     parts.push_back(args[k].partial);
    // }
    
    // pairwise_add(parts, sum);
    
    // for (unsigned k=0; k<nthreads; ++k) { 
    //     mpfr_clear(args[k].partial);
    // }
    
    // for (unsigned k = 0; k < nthreads; ++k) {
    //     mpfr_add(sum, sum, args[k].partial, MPFR_RNDN);
    //     mpfr_clear(args[k].partial);
    // }

    // mpfr_mul_ui(sum, sum, 2, MPFR_RNDN); // 2 * Σ
    // mpfr_set(Z, sum, MPFR_RNDN);

    mpfr_t sum; 
    mpfr_init2(sum, MPFR_PRECISION); 
    mpfr_set_zero(sum, 0);

    // collect pointers to thread partials
    std::vector<mpfr_ptr> parts; 
    parts.reserve(nthreads);
    for (unsigned k = 0; k < nthreads; ++k) {
        parts.push_back(args[k].partial);  // mpfr_t decays to mpfr_ptr here
    }

    // reduce in-place and then clear the thread partials
    pairwise_add_inplace(parts, sum);
    for (unsigned k = 0; k < nthreads; ++k) {
        mpfr_clear(args[k].partial);
    }

    // final Z = 2 * Σ
    mpfr_mul_ui(sum, sum, 2, MPFR_RNDN);
    mpfr_set(Z, sum, MPFR_RNDN);

    // (optional) report an error bar from the remainder bound
    mpfr_t R; mpfr_init2(R, MPFR_PRECISION);
    remainder_bound(R, t);
    log_mpfr("[main] |R_bound|", R, 20);
    mpfr_clear(R);

    // good: log AFTER setting Z; don't reference R again
    log_mpfr("[main] 2*sum", sum, 20);
    log_mpfr("[main] Z", Z, 20);
    
    if (mpfr_nan_p(Z))  std::cerr << "[main] Z is NaN!\n";
    if (mpfr_inf_p(Z))  std::cerr << "[main] Z is Inf!\n";

    mpfr_clears(theta, pi, two_pi, tmp, sum, (mpfr_ptr)0);
}



