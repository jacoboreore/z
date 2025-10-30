// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: os_valley.hpp
// Purpose:
//   Precision-shimmed OS (Odlyzko–Schönhage) block planner and evaluator
//   for approximate Hardy Z(t) sampling, plus quick analysis/classification
//   helpers for valley scanning and Newton seeding.
// -----------------------------------------------------------------------------

#pragma once

#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#include "logs.h"
#include "settings.h"

// ===== constants & math wrappers usable by both native and MPFR backends =====

/**
 * @brief π for the chosen Real backend.
 * @tparam T Real type (long double or mpfr::mpreal).
 * @return π in type T.
 */
template <class T>
inline T os_pi();

template <>
inline long double os_pi<long double>()
{
    return 3.141592653589793238462643383279502884L;
}

#if OS_HAS_MPFR
template <>
inline mpfr::mpreal os_pi<mpfr::mpreal>()
{
    return mpfr::const_pi();
}
#endif

/**
 * @brief 2π for the chosen Real backend.
 * @tparam T Real type.
 * @return 2π in type T.
 */
template <class T>
inline T os_two_pi() { return (T)2 * os_pi<T>(); }

/**
 * @brief Absolute value with ADL for std and MPFR.
 */
template <class T>
inline T OsAbs(const T &x)
{
    using std::abs;
#if OS_HAS_MPFR
    using mpfr::abs;
#endif
    return abs(x);
}

/**
 * @brief Cosine with ADL for std and MPFR.
 */
template <class T>
inline T OsCos(const T &x)
{
    using std::cos;
#if OS_HAS_MPFR
    using mpfr::cos;
#endif
    return cos(x);
}

/**
 * @brief Sine with ADL for std and MPFR.
 */
template <class T>
inline T OsSin(const T &x)
{
    using std::sin;
#if OS_HAS_MPFR
    using mpfr::sin;
#endif
    return sin(x);
}

/**
 * @brief Natural logarithm with ADL for std and MPFR.
 */
template <class T>
inline T OsLog(const T &x)
{
    using std::log;
#if OS_HAS_MPFR
    using mpfr::log;
#endif
    return log(x);
}

/**
 * @brief Square root with ADL for std and MPFR.
 */
template <class T>
inline T OsSqrt(const T &x)
{
    using std::sqrt;
#if OS_HAS_MPFR
    using mpfr::sqrt;
#endif
    return sqrt(x);
}

/**
 * @brief Floor with ADL for std and MPFR.
 */
template <class T>
inline T OsFloor(const T &x)
{
    using std::floor;
#if OS_HAS_MPFR
    using mpfr::floor;
#endif
    return floor(x);
}

/**
 * @brief Ceil with ADL for std and MPFR.
 */
template <class T>
inline T OsCeil(const T &x)
{
    using std::ceil;
#if OS_HAS_MPFR
    using mpfr::ceil;
#endif
    return ceil(x);
}

/**
 * @brief Round-to-nearest (half-away-from-zero) to int64.
 * @tparam T Real backend type.
 * @param x Input value.
 * @return Rounded 64-bit integer.
 */
template <class T>
inline long long OsLLRound(const T &x)
{
    if (x >= T(0))
        return (long long)OsFloor(x + T(0.5));
    else
        return (long long)OsCeil(x - T(0.5));
}

// =============================== end precision shim ===============================

namespace osvalley
{
    using Real = OS_REAL_T;
    using Cpx = std::complex<Real>;

    // ---------------------------------- config ----------------------------------

    /**
     * @brief OS block configuration (bridged from global CFG()).
     *
     * Values mirror the valley-scanner configuration to maintain consistent
     * thresholds and sampling strategies.
     */
    struct OSConfig
    {
        Real enable_min_t = CFG().enable_min_t;              ///< Only use OS when t0 ≥ this.
        int samples_per_spacing = CFG().samples_per_spacing; ///< Δt = spacing / S.
        std::size_t K_points = CFG().K_points;               ///< Points computed per block (≤ M).
        Real delta_log_max = CFG().delta_log_max;            ///< Bin width in log n (smaller → more accurate).
        int poly_deg_R = CFG().poly_deg_R;                   ///< Taylor degree (2..6 reasonable).
        bool return_main_sum = CFG().return_main_sum;        ///< If true: skip θ factor (debug output).

        // Classification thresholds (mirrors valley cfg)
        Real tol_refine = CFG().tol_refine;             ///< Newton "near enough" tolerance.
        Real z_far = CFG().z_far;                       ///< |Zc| above this is considered far (order-1).
        Real ratio_far = CFG().ratio_far;               ///< |dt_pred|/spacing threshold for off-center.
        Real min_mult = CFG().min_mult;                 ///< Accept if min|Z| ≤ min_mult · tol_refine.
        Real eps_abs = CFG().eps_abs;                   ///< Absolute guard band for trusted sign.
        Real eps_rel = CFG().eps_rel;                   ///< Relative guard band (fraction of median |Z|).
        Real flip_center_frac = CFG().flip_center_frac; ///< Flip must be within ±fraction·spacing of center to auto-accept.
    };

    // ------------------------------- theta series --------------------------------

    /**
     * @brief θ and its derivatives at t0.
     */
    struct ThetaSeries
    {
        Real th0; ///< θ(t0)
        Real th1; ///< θ'(t0)
        Real th2; ///< θ''(t0)
        Real th3; ///< θ'''(t0)
    };

    // ------------------------------- OS plan state --------------------------------

    /**
     * @brief Plan for one OS block centered at t0.
     *
     * Contains sampling geometry, per-degree coefficient bins, FFT buffers,
     * and precomputed phase factors.
     */
    struct OSPlan
    {
        // Chosen parameters
        Real t0{};        ///< Center t0.
        Real dt{};        ///< Step Δt between consecutive outputs.
        Real dlog{};      ///< Bin width in log n.
        std::size_t M{};  ///< FFT length (power-of-two).
        std::size_t K{};  ///< Number of outputs requested.
        int R{};          ///< Taylor degree.
        ThetaSeries TS{}; ///< θ series at t0.

        // Per-degree input to FFT: C_r[m]
        std::vector<std::vector<Cpx>> C; ///< Size R×M.

        // FFT outputs per degree: G_r[k]
        std::vector<std::vector<Cpx>> G; ///< Size R×M.

        // Precomputed exp(+i t0 m Δlog)
        std::vector<Cpx> E_m; ///< Length M.
    };

    // --------------------------- OS ratio / step suggestion -----------------------

    /**
     * @brief Ratio/step suggestion computed from an OS block for valley scanning.
     */
    struct OSRatioStep
    {
        osvalley::Real dt_suggest{}; ///< Suggested signed Δt to move.
        osvalley::Real Z0{};         ///< Z at start index.
        osvalley::Real Z1_pred{};    ///< Z at suggested point (or flip edge).
        osvalley::Real ratio{};      ///< |Z0| / |Z1_pred|.
        bool flip_ahead{};           ///< Whether a sign flip is encountered before reaching the target step.
        int i0{-1};                  ///< Start index used.
        int i1{-1};                  ///< End index chosen.
    };

    // ------------------------------ assembled block -------------------------------

    /**
     * @brief OS block samples and optional pre-θ main sum.
     */
    struct BlockResult
    {
        Real t0{}, dt{};
        std::vector<Real> t; ///< Length K: t[k].
        std::vector<Real> Z; ///< Length K: approximate Z(t[k]).
        std::vector<Cpx> F;  ///< Length K: optional main sum before θ (debug).
    };

    // ---------------------------------- utilities ---------------------------------

    /**
     * @brief Next power-of-two ≥ n (clamped at 1).
     */
    static inline std::size_t next_pow2(std::size_t n)
    {
        if (n <= 1)
            return 1;
        --n;

        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
#if SIZE_MAX > 0xffffffffu
        n |= n >> 32;
#endif
        return n + 1;
    }

    /**
     * @brief In-place iterative Cooley–Tukey FFT (power-of-two length).
     *
     * Forward transform is unnormalized; inverse divides by n.
     * @param a    Input/output vector (complex).
     * @param sign -1 = forward, +1 = inverse.
     */
    static inline void fft_inplace(std::vector<Cpx> &a, int sign /*-1: forward, +1: inverse*/)
    {
        // Iterative Cooley–Tukey, power-of-two length, unnormalized forward
        const std::size_t n = a.size();
        assert(n && ((n & (n - 1)) == 0));

        // Bit-reversal permutation
        for (std::size_t i = 1, j = 0; i < n; i++)
        {
            std::size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }

        // Stages
        for (std::size_t len = 2; len <= n; len <<= 1)
        {
            const Real ang = (sign < 0 ? (Real)-1 : (Real) + 1) * (os_two_pi<Real>() / (Real)len);
            const Cpx wlen = {OsCos(ang), OsSin(ang)};
            for (std::size_t i = 0; i < n; i += len)
            {
                Cpx w{1, 0};
                for (std::size_t j = 0; j < len / 2; ++j)
                {
                    const Cpx u = a[i + j];
                    const Cpx v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }

        if (sign > 0)
        {
            const Real inv = (Real)1 / (Real)n;
            for (auto &x : a)
                x *= inv;
        }
    }

    /**
     * @brief θ(t) ≈ (t/2)log(t/2π) − t/2 − π/8 + 1/(48t) + 7/(5760 t^3) + …
     *        θ'(t) = 0.5·log(t/2π) + O(t^-2), θ''(t)=1/(2t), θ'''(t)=−1/(2t^2).
     *
     * Three-term Stirling is sufficient at large t (e.g., 1e14).
     * @param t0 Expansion center.
     * @return ThetaSeries {θ, θ', θ'', θ'''} at t0.
     */
    static inline ThetaSeries theta_series_at(Real t0)
    {
        const Real two_pi = os_two_pi<Real>();
        const Real u = t0 / two_pi;

        const Real th0 = (Real)0.5 * t0 * OsLog(u) - (Real)0.5 * t0 - os_pi<Real>() / (Real)8 + (Real)1 / ((Real)48 * t0) + (Real)7 / ((Real)5760 * t0 * t0 * t0);

        const Real th1 = (Real)0.5 * OsLog(u);
        const Real th2 = (Real)1 / ((Real)2 * t0);
        const Real th3 = (Real)-1 / ((Real)2 * t0 * t0);
        return {th0, th1, th2, th3};
    }

    /**
     * @brief Evaluate θ(t0 + s) from a local ThetaSeries via a cubic expansion.
     * @param S ThetaSeries at t0.
     * @param s Offset s = t − t0.
     * @return θ(t0 + s).
     */
    static inline Real theta_eval_from_series(const ThetaSeries &S, Real s)
    {
        return S.th0 + S.th1 * s + (S.th2 * (s * s)) / (Real)2 + (S.th3 * (s * s * s)) / (Real)6;
    }

    /**
     * @brief RS truncation length N ≈ floor(sqrt(t/(2π))).
     * @param t Height.
     * @return N as uint64_t, clamped to at least 1.
     */
    static inline std::uint64_t rs_length_N(Real t)
    {
        Real Nf = OsSqrt(t / (os_two_pi<Real>()));
        if (Nf < (Real)1)
            Nf = (Real)1;
        return (std::uint64_t)OsFloor(Nf);
    }

    /**
     * @brief Local zero spacing estimate ≈ π / log(t/2π).
     * @param t Height.
     * @return Spacing estimate.
     */
    static inline Real zero_spacing_est(Real t)
    {
        return os_pi<Real>() / OsLog(t / (os_two_pi<Real>()));
    }

    /**
     * @brief Build an OS plan centered at t0 with K_points outputs.
     *
     * Chooses Δt from spacing/S, selects M as the next power-of-two satisfying
     * Δlog = 2π/(M Δt) ≤ delta_log_max and M ≥ K. Precomputes exp(+i t0 m Δlog).
     *
     * @param t0       Block center.
     * @param K_points Number of outputs to produce.
     * @param cfg      OS configuration.
     * @return OSPlan with geometry and buffers ready for accumulation.
     */
    static inline OSPlan make_plan(Real t0, std::size_t K_points, const OSConfig &cfg)
    {
        OSPlan P;
        P.t0 = t0;
        P.K = K_points;

        // Δt from spacing / S
        const Real spacing = zero_spacing_est(t0);
        P.dt = spacing / (Real)std::max(1, cfg.samples_per_spacing);

        // Pick M so that Δlog = 2π / (M Δt) ≤ delta_log_max and M ≥ K
        const Real need = os_two_pi<Real>() / (P.dt * cfg.delta_log_max); // minimal M (pre pow2)
        std::size_t Mmin = (std::size_t)OsCeil(need);
        P.M = next_pow2(std::max<std::size_t>(P.K, Mmin));
        P.dlog = os_two_pi<Real>() / ((Real)P.M * P.dt);

        P.R = std::max(2, cfg.poly_deg_R);
        P.TS = theta_series_at(P.t0);

        P.C.assign(P.R, std::vector<Cpx>(P.M, Cpx{0, 0}));
        P.G.assign(P.R, std::vector<Cpx>(P.M, Cpx{0, 0}));
        P.E_m.resize(P.M);

        // Precompute E_m = exp(i t0 m Δlog)
        for (std::size_t m = 0; m < P.M; ++m)
        {
            const Real ang = P.t0 * ((Real)m * P.dlog);
            P.E_m[m] = {OsCos(ang), OsSin(ang)};
        }
        return P;
    }

    /**
     * @brief Accumulate per-bin Taylor coefficients C_r[m] for n = 1..N.
     *
     * For each n, computes its log-bin m, the base coefficient, and accumulates
     * the first R terms of the Taylor expansion in Δt around t0.
     *
     * @param P Plan with geometry and buffers.
     * @param N RS truncation length.
     */
    static inline void accumulate_bins(OSPlan &P, std::uint64_t N)
    {
        for (std::uint64_t n = 1; n <= N; ++n)
        {
            const Real x = OsLog((Real)n); // log n
            const Real mreal = x / P.dlog;
            std::int64_t m = (std::int64_t)OsLLRound(mreal);

            // Wrap m into [0, M)
            m %= (std::int64_t)P.M;
            if (m < 0)
                m += (std::int64_t)P.M;

            const Real delta = x - (Real)m * P.dlog; // δ ∈ ~[-Δlog/2, +Δlog/2]

            // Base coefficient: n^{-1/2} * e^{-i t0 log n} * e^{+i t0 m Δlog}
            const Real invsqrt = (Real)1 / OsSqrt((Real)n);
            const Real a = -(P.t0 * x);
            const Cpx phase = {OsCos(a), OsSin(a)}; // e^{-i t0 log n}
            const Cpx c_full = invsqrt * phase * P.E_m[(std::size_t)m];

            // b = (-i Δt δ)
            const Cpx b(0, -(P.dt * delta)); // purely imaginary

            // Accumulate R terms: c_full * b^r / r!
            Cpx pow(1, 0);
            Real invfact = (Real)1; // 1/0!, then *= 1/(r+1)
            for (int r = 0; r < P.R; ++r)
            {
                P.C[r][(std::size_t)m] += c_full * (pow * invfact);
                pow *= b;
                invfact /= (Real)(r + 1);
            }
        }
    }

    /**
     * @brief Run forward FFT (unnormalized) on each degree r to obtain G_r.
     * @param P Plan whose C_r buffers are transformed in-place to G_r.
     */
    static inline void transform(OSPlan &P)
    {
        for (int r = 0; r < P.R; ++r)
        {
            P.G[r] = P.C[r];         // copy to preserve C (debug/inspection)
            fft_inplace(P.G[r], -1); // forward unnormalized
        }
    }

    /**
     * @brief Evaluate a completed OS block to produce {t, Z, F}.
     *
     * For each k, computes F_k ≈ ∑_{r=0}^{R-1} k^r G_r[k]. If return_main_sum
     * is false, multiplies by e^{i θ(t)} and outputs Z_k = 2·Re(F_k e^{iθ}).
     *
     * @param P   Completed plan (G_r available).
     * @param cfg OS configuration (return_main_sum flag).
     * @return BlockResult with t, Z, and F vectors.
     */
    static inline BlockResult evaluate_block(const OSPlan &P, const OSConfig &cfg)
    {
        BlockResult R;
        R.t0 = P.t0;
        R.dt = P.dt;
        R.t.resize(P.K);
        R.Z.resize(P.K);
        R.F.resize(P.K);

        for (std::size_t k = 0; k < P.K; ++k)
        {
            // F_k ≈ sum_{r=0}^{R-1} k^r · G_r[k]
            Cpx Fk(0, 0);
            Real kr = (Real)1;
            for (int r = 0; r < P.R; ++r)
            {
                Fk += kr * P.G[r][k];
                kr *= (Real)k;
            }
            R.F[k] = Fk;

            const Real t = P.t0 + (Real)k * P.dt;
            R.t[k] = t;

            if (cfg.return_main_sum)
            {
                R.Z[k] = (Real)Fk.real(); // raw main sum (debug)
            }
            else
            {
                // Multiply by e^{i θ(t)} and take 2·Re(...)
                const Real s = (Real)k * P.dt; // s = t − t0
                const Real th = theta_eval_from_series(P.TS, s);
                const Cpx eith(OsCos(th), OsSin(th));
                R.Z[k] = (Real)2 * (Fk * eith).real();
            }
        }
        return R;
    }

    /**
     * @brief Full OS pipeline for one block: plan → accumulate → FFT → evaluate.
     *
     * Throws if t0 < enable_min_t to route callers to a small-Z or RS–Gabcke path.
     *
     * @param t0  Block center.
     * @param cfg OS configuration.
     * @return BlockResult with approximate Z(t) samples and optional F(t).
     */
    inline BlockResult scan_valley_OS(Real t0, const OSConfig &cfg)
    {
        if (t0 < cfg.enable_min_t)
            throw std::runtime_error(
                "OS z disabled for small t0; use RS–Gabcke or small-z path, or raise enable_min_t.");

        // Choose RS length once for the block
        const std::uint64_t N = rs_length_N(t0);

        OSPlan P = make_plan(t0, cfg.K_points, cfg);
        accumulate_bins(P, N);
        transform(P);

        // Debug prints 
        // const Real spacing = zero_spacing_est(t0);
        // const Real cov_log = (Real)P.M * P.dlog; // 2π/Δt
        // fprintf(stderr,
        //   "[OS] t0=%.0Lf  N=%llu  M=%zu  R=%d  Δt=%.3Le  spacing≈%.3Le  Δlog=%.3Le  cover_log=%.3Le\n",
        //   (long double)t0, (unsigned long long)N, P.M, P.R, (long double)P.dt,
        //   (long double)spacing, (long double)P.dlog, (long double)cov_log);

        return evaluate_block(P, cfg);
    }

    /**
     * @brief Return sign with hysteresis thresholding.
     * @param x   Value.
     * @return +1 if x>0, −1 if x<0, 0 if x==0.
     */
    static inline int os_sgn(const osvalley::Real &x) { return (x > (Real)0) - (x < (Real)0); }

    /**
     * @brief Compute an OSRatioStep from a block, suggesting the next Δt.
     *
     * Uses central slope and a target ratio to estimate a bounded step, then
     * walks sample-by-sample to stop at the first trusted sign flip if present.
     *
     * @param block      Completed block (t, Z, dt).
     * @param t_cur      Current t.
     * @param spacing    Local zero spacing estimate.
     * @param R_target   Desired improvement factor |Z0|/|Z1| (default 5).
     * @param s_min_frac Minimum step as fraction of spacing (default 0.02).
     * @param s_max_frac Maximum step as fraction of spacing (default 0.7).
     * @param g_floor    Gradient floor to avoid division by ~0 (default 1e-12).
     * @return OSRatioStep with indices, ratio, flip flag, and Δt suggestion.
     */
    OSRatioStep os_ratio_step_from_block(
        const osvalley::BlockResult &block,
        osvalley::Real t_cur,
        osvalley::Real spacing,
        osvalley::Real R_target = osvalley::Real(5),
        osvalley::Real s_min_frac = osvalley::Real(0.02),
        osvalley::Real s_max_frac = osvalley::Real(0.7),
        osvalley::Real g_floor = osvalley::Real(1e-12))
    {
        OSRatioStep out;

        const int K = (int)block.Z.size();
        if (K < 3)
            return out; // Not enough samples

        const osvalley::Real dt = block.dt;

        // Nearest sample index to t_cur
        int i0 = (int)OsLLRound((t_cur - block.t0) / dt);
        i0 = std::clamp(i0, 1, K - 2); // Keep room for central difference
        out.i0 = i0;

        const osvalley::Real Z0 = block.Z[i0];
        out.Z0 = Z0;

        // Central difference for dZ/dt
        const osvalley::Real g = (block.Z[i0 + 1] - block.Z[i0 - 1]) / ((Real)2 * dt);
        const osvalley::Real g_abs = std::max(OsAbs(g), g_floor);

        // Target step length for ratio R: s = (1 − 1/R)·|Z0|/|g|
        osvalley::Real s = ((Real)1 - (Real)1 / R_target) * OsAbs(Z0) / g_abs;

        // Clamp by spacing
        const osvalley::Real s_min = s_min_frac * spacing;
        const osvalley::Real s_max = s_max_frac * spacing;
        s = std::clamp(s, s_min, s_max);

        // Descent direction: d|Z|/dt ≈ sign(Z0) * g → choose step sign = −sign(Z0 * g)
        int dir = -os_sgn(Z0 * g);
        if (dir == 0)
            dir = -os_sgn(Z0); // Fallback

        // Walk along OS samples to locate flip before reaching s
        int max_hops = (int)OsFloor(OsAbs(s / dt));
        max_hops = std::max(max_hops, 1);
        int i = i0;
        bool flip = false;
        for (int hop = 0; hop < max_hops; ++hop)
        {
            const int inext = i + dir;
            if (inext <= 0 || inext >= K)
                break;
            if (block.Z[i] * block.Z[inext] <= (Real)0)
            {
                // Found sign flip; land on edge
                i = inext;
                flip = true;
                break;
            }
            i = inext;
        }

        // If no flip, step to the bounded target index
        if (!flip)
        {
            int i_target = i0 + dir * max_hops;
            i_target = std::clamp(i_target, 0, K - 1);
            i = i_target;
        }

        out.i1 = i;
        out.Z1_pred = block.Z[i];
        out.flip_ahead = flip;

        // Observed ratio from OS samples
        const osvalley::Real Z1abs = std::max(OsAbs(out.Z1_pred), (Real)1e-30);
        out.ratio = OsAbs(Z0) / Z1abs;

        // Suggested Δt (signed)
        out.dt_suggest = (i - i0) * dt;

        return out;
    }

    // ------------------------------- classification -------------------------------

    /**
     * @brief Verdicts from OS-based classification.
     */
    enum class OSVerdict
    {
        ACCEPT,
        REJECT,
        RESCAN
    };

    /**
     * @brief Decision bundle used by valley scanner to gate Newton refinement.
     */
    struct OSDecision
    {
        OSVerdict verdict{OSVerdict::RESCAN}; ///< Decision: ACCEPT / REJECT / RESCAN.
        long double t_seed{};                 ///< Seed for Newton (flip secant or min-based).
        bool has_flip{};                      ///< Whether a flip was detected.
        long double min_absZ{};               ///< Minimum |Z| within the block (FYI).
        long double absZc{};                  ///< |Z| at the center (FYI).
        long double ratio{};                  ///< |dt_pred|/spacing (FYI).
    };

    /**
     * @brief Trusted sign with hysteresis.
     * @param z   Value.
     * @param eta Band threshold.
     * @return +1, 0, or −1.
     */
    inline int trusted_sign(long double z, long double eta)
    {
        if (z > eta)
            return +1;
        if (z < -eta)
            return -1;
        return 0;
    }

    /**
     * @brief Median of absolute values as a robust amplitude scale.
     * @param v Input vector.
     * @return Median(|v|).
     */
    inline long double median_abs(const std::vector<long double> &v)
    {
        if (v.empty())
            return 0;
        std::vector<long double> u(v.begin(), v.end());
        for (auto &x : u)
            x = OsAbs(x);
        auto mid = u.begin() + u.size() / 2;
        std::nth_element(u.begin(), mid, u.end());
        return *mid;
    }

    /**
     * @brief Classify a completed block and emit a decision for refinement.
     *
     * Detects a trusted sign flip near the center, or accepts near-min blocks.
     * Rejects far/off-center cases; otherwise requests a rescan with adjusted t0.
     *
     * @param B  Completed block.
     * @param oc OS configuration thresholds.
     * @return OSDecision containing verdict and suggested Newton seed.
     */
    inline OSDecision os_classify_block(const osvalley::BlockResult &B, const osvalley::OSConfig &oc)
    {
        using osvalley::Real;
        OSDecision out{};

        const int K = (int)B.Z.size();
        if (K < 2)
            return out;
        const int c = K / 2;

        const long double spacing =
            (long double)os_pi<Real>() / OsLog((long double)B.t[c] / ((long double)2 * (long double)os_pi<Real>()));

        // Center stats
        const long double Zc = (long double)B.Z[c];
        const long double absZc = OsAbs(Zc);

        long double dZdt = 0;
        if (K >= 3)
            dZdt = ((long double)B.Z[c + 1] - (long double)B.Z[c - 1]) / ((long double)2 * (long double)B.dt);

        const long double dt_pred = (dZdt != 0 ? Zc / dZdt : 0);
        const long double ratio = (spacing > 0 ? OsAbs(dt_pred / spacing) : 1e9L);

        // Robust guard band for sign trust
        std::vector<long double> zvals(K);
        for (int i = 0; i < K; ++i)
            zvals[i] = (long double)B.Z[i];
        const long double amp = median_abs(zvals);
        const long double eta = std::max((long double)oc.eps_abs, (long double)(oc.eps_rel * amp));

        // Scan for trusted flip and locate nearest to center
        bool flip = false;
        int flip_iL = -1;
        long double t_seed = 0;

        for (int i = 0; i + 1 < K; ++i)
        {
            const int s0 = trusted_sign(zvals[i], eta);
            const int s1 = trusted_sign(zvals[i + 1], eta);
            if (s0 && s1 && s0 != s1)
            {
                flip = true;
                flip_iL = i;

                // Secant seed in [i, i+1]
                const long double z0 = zvals[i], z1 = zvals[i + 1];
                const long double t0 = (long double)B.t0 + (long double)i * (long double)B.dt;
                t_seed = t0 - z0 * ((long double)B.dt / (z1 - z0));
                break; // Take first; near center for well-planned blocks
            }
        }

        // Find min |Z| in block (for near-but-no-flip)
        int j = 0;
        for (int i = 1; i < K; ++i)
            if (OsAbs(zvals[i]) < OsAbs(zvals[j]))
                j = i;

        const long double min_absZ = OsAbs(zvals[j]);
        if (!flip)
        {
            // Fallback Newton seed at min sample if slope is usable
            long double slope = 0;
            if (j > 0 && j + 1 < K)
                slope = (zvals[j + 1] - zvals[j - 1]) / ((long double)2 * (long double)B.dt);
            else if (j + 1 < K)
                slope = (zvals[j + 1] - zvals[j]) / (long double)B.dt;
            else if (j > 0)
                slope = (zvals[j] - zvals[j - 1]) / (long double)B.dt;

            t_seed = ((long double)B.t0 + (long double)j * (long double)B.dt) - (OsAbs(slope) > 1e-12L ? zvals[j] / slope : 0);
        }

        // Final classification
        const bool flip_near_center =
            flip && (OsAbs(((long double)flip_iL + 0.5L - (long double)c) * (long double)B.dt) <= oc.flip_center_frac * spacing);
        const bool near_by_min = (min_absZ <= oc.min_mult * oc.tol_refine);
        const bool far_center = (absZc >= oc.z_far) && (ratio >= oc.ratio_far);

        if (flip_near_center || near_by_min)
        {
            out.verdict = OSVerdict::ACCEPT;
        }
        else if (!flip && !near_by_min && far_center)
        {
            out.verdict = OSVerdict::REJECT;
        }
        else
        {
            out.verdict = OSVerdict::RESCAN;
        }

        out.t_seed = t_seed;
        out.has_flip = flip;
        out.min_absZ = min_absZ;
        out.absZc = absZc;
        out.ratio = ratio;
        return out;
    }

} // namespace osvalley

// =============================== quick analysis API ===============================

namespace osquick
{
    using osvalley::Real;

    /**
     * @brief One-shot block evaluation summary for dashboards / logs.
     */
    struct OSEval
    {
        Real t_center{};
        Real Zc{}, absZc{};
        Real dZdt{};
        Real dt_pred{};  ///< Newton offset from center.
        Real spacing{};  ///< ≈ π/log(t/2π).
        Real ratio{};    ///< |dt_pred| / spacing.
        Real min_absZ{}; ///< min_i |Z[i]|.
        int idx_min{-1}; ///< Index of min_absZ.
        bool has_flip{false};
        bool is_near{false};
        Real suggested_step{}; ///< Suggested step = −dt_pred.
    };

    /**
     * @brief Local zero spacing estimate (duplicate helper for convenience).
     * @param t Height.
     * @return ≈ π / log(t/2π).
     */
    inline Real local_spacing(Real t)
    {
        const Real two_pi = (Real)2 * os_pi<Real>();
        return os_pi<Real>() / OsLog(t / two_pi);
    }

    /**
     * @brief Analyze a BlockResult and compute a compact OSEval summary.
     *
     * Conditions for is_near are inclusive:
     *   - |Zc| ≤ absZ_tol, or
     *   - any sign flip inside the block, or
     *   - |dt_pred|/spacing ≤ ratio_tol, or
     *   - min_i |Z[i]| ≤ absZ_tol.
     *
     * @param b         Completed block.
     * @param absZ_tol  Absolute |Z| threshold (default 0.5).
     * @param ratio_tol Ratio threshold for |dt_pred|/spacing (default 0.30).
     * @return OSEval with derived quantities and suggested_step.
     */
    inline OSEval analyze_block(const osvalley::BlockResult &b,
                                Real absZ_tol = (Real)0.5,
                                Real ratio_tol = (Real)0.30)
    {
        OSEval ev{};
        const int K = (int)b.Z.size();
        const int c = K / 2;

        ev.t_center = b.t[c];
        ev.Zc = b.Z[c];
        ev.absZc = OsAbs((long double)ev.Zc); // Safe cast if Real = long double
        ev.spacing = local_spacing(ev.t_center);

        // Slope at center (guard edges)
        if (K >= 3)
        {
            ev.dZdt = (b.Z[c + 1] - b.Z[c - 1]) / ((Real)2 * b.dt);
        }
        else
        {
            ev.dZdt = (Real)0;
        }

        // One-shot Newton offset (if slope is usable)
        ev.dt_pred = (ev.dZdt != (Real)0) ? (ev.Zc / ev.dZdt) : (Real)0;
        ev.ratio = (ev.spacing > (Real)0 ? OsAbs((long double)(ev.dt_pred / ev.spacing)) : (Real)1e9);

        // Scan for sign flip & min |Z|
        ev.min_absZ = OsAbs((long double)b.Z[0]);
        ev.idx_min = 0;
        for (int i = 1; i < K; ++i)
        {
            const Real zi = b.Z[i];
            if ((b.Z[i - 1] <= (Real)0 && zi > (Real)0) ||
                (b.Z[i - 1] >= (Real)0 && zi < (Real)0))
            {
                ev.has_flip = true;
            }
            const Real a = OsAbs((long double)zi);
            if (a < ev.min_absZ)
            {
                ev.min_absZ = a;
                ev.idx_min = i;
            }
        }

        // Verdict: any condition qualifies as "near enough"
        ev.is_near = (ev.absZc <= absZ_tol) || ev.has_flip || (ev.ratio <= ratio_tol) || (ev.min_absZ <= absZ_tol);
        ev.suggested_step = -ev.dt_pred;

        return ev;
    }

    /**
     * @brief Emit a compact textual summary for OSEval (stdout).
     * @param ev Evaluation bundle.
     */
    inline void print_summary(const OSEval &ev)
    {
        std::cout << std::scientific << std::setprecision(18);
        std::cout << "t_c=" << ev.t_center
                  << "\n|Zc|=" << ev.absZc
                  << "\ndZdt=" << ev.dZdt
                  << "\ndt_pred=" << ev.dt_pred
                  << "\nspacing=" << ev.spacing
                  << "\nratio=" << ev.ratio
                  << "\nmin|Z|=" << ev.min_absZ
                  << "\nflip=" << (ev.has_flip ? 1 : 0)
                  << "\nverdict=" << (ev.is_near ? "NEAR" : "FAR")
                  << "\nsuggested_step=" << ev.suggested_step
                  << "\n";
    }

} // namespace osquick
