#pragma once
#include "logs.h"
#include "settings.h"
#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <cassert>

// ---- constants & math wrappers that work for both backends ----
template <class T>
inline T os_pi();
template <>
inline long double os_pi<long double>() { return 3.141592653589793238462643383279502884L; }
#if OS_HAS_MPFR
template <>
inline mpfr::mpreal os_pi<mpfr::mpreal>() { return mpfr::const_pi(); }
#endif

template <class T>
inline T os_two_pi() { return (T)2 * os_pi<T>(); }

template <class T>
inline T OsAbs(const T &x)
{
    using std::abs;
#if OS_HAS_MPFR
    using mpfr::abs;
#endif
    return abs(x);
}

template <class T>
inline T OsCos(const T &x)
{
    using std::cos;
#if OS_HAS_MPFR
    using mpfr::cos;
#endif
    return cos(x);
}
template <class T>
inline T OsSin(const T &x)
{
    using std::sin;
#if OS_HAS_MPFR
    using mpfr::sin;
#endif
    return sin(x);
}
template <class T>
inline T OsLog(const T &x)
{
    using std::log;
#if OS_HAS_MPFR
    using mpfr::log;
#endif
    return log(x);
}
template <class T>
inline T OsSqrt(const T &x)
{
    using std::sqrt;
#if OS_HAS_MPFR
    using mpfr::sqrt;
#endif
    return sqrt(x);
}
template <class T>
inline T OsFloor(const T &x)
{
    using std::floor;
#if OS_HAS_MPFR
    using mpfr::floor;
#endif
    return floor(x);
}
template <class T>
inline T OsCeil(const T &x)
{
    using std::ceil;
#if OS_HAS_MPFR
    using mpfr::ceil;
#endif
    return ceil(x);
}
// round-to-nearest (half-away-from-zero), returns int64
template <class T>
inline long long OsLLRound(const T &x)
{
    if (x >= T(0))
        return (long long)OsFloor(x + T(0.5));
    else
        return (long long)OsCeil(x - T(0.5));
}
// ======= end precision shim =======

// // ------- numerics backend (long double for stability at t>=1e14) -------
// #ifndef OS_REAL_T
// #define OS_REAL_T long double
// #endif

// #ifndef os_pi
// #define os_pi 3.141592653589793238462643383279502884L
// #endif

namespace osvalley
{
    using Real = OS_REAL_T;
    using Cpx = std::complex<Real>;

    // ---------------------- config ----------------------
    struct OSConfig
    {
        Real enable_min_t = CFG().enable_min_t;        // only use OS when t0 >= this
        int samples_per_spacing = CFG().samples_per_spacing; // Δt = spacing / S
        std::size_t K_points = CFG().K_points;               // points computed per block (≤ M)
        Real delta_log_max = CFG().delta_log_max;            // bin width in log n (smaller → more accurate)
        int poly_deg_R = CFG().poly_deg_R;                   // Taylor degree (2..6 reasonable)
        bool return_main_sum = CFG().return_main_sum;        // if true: skip θ factor (debug)

        // classification thresholds (mirrors valley cfg)
        Real tol_refine = CFG().tol_refine;             // what you pass to Newton as “near enough”
        Real z_far = CFG().z_far;                       // |Zc| above this looks far (order-1)
        Real ratio_far = CFG().ratio_far;               // |dt_pred|/spacing above this looks off-center
        Real min_mult = CFG().min_mult;                 // accept if min|Z| <= min_mult * tol_refine
        Real eps_abs = CFG().eps_abs;                   // absolute guard band for ‘trusted sign’
        Real eps_rel = CFG().eps_rel;                   // relative guard band (fraction of median |Z|)
        Real flip_center_frac = CFG().flip_center_frac; // flip must be within ±0.6*spacing of center to auto-accept
    };

    // θ and its derivatives @ t0 (Stirling side good for t≥1e9)
    struct ThetaSeries
    {
        Real th0, th1, th2, th3; // θ, θ', θ'', θ'''
    };

    // ---------------------- main planner ----------------------
    struct OSPlan
    {
        // chosen parameters
        Real t0{}, dt{}, dlog{}; // center t0, step Δt, bin width Δlog
        std::size_t M{};         // FFT length (power of two)
        std::size_t K{};         // number of outputs requested
        int R{};                 // Taylor degree
        ThetaSeries TS{};        // θ series @ t0

        // per-degree input to FFT: C_r[m]
        std::vector<std::vector<Cpx>> C; // R x M
        // FFT outputs per degree
        std::vector<std::vector<Cpx>> G; // R x M

        // precomputed exp(+i t0 m Δlog)
        std::vector<Cpx> E_m; // length M
    };

    // ---- OS ratio + step suggestion for valley scanning ----
    // Inputs:
    //   block : OS block (t0, dt, vectors t, Z, F)
    //   t_cur : your current t
    //   spacing : local zero spacing ~ pi / log(t/2pi) (use your est_spacing())
    //   R_target : desired improvement factor |Z0|/|Z1| (e.g., 5.0)
    //   cfg clamps: s_min_frac, s_max_frac (fractions of spacing), g_floor
    // Output:
    //   step suggestion, observed ratio from OS samples, and flip detection
    struct OSRatioStep
    {
        osvalley::Real dt_suggest{}; // signed Δt (direction to move)
        osvalley::Real Z0{};         // OS Z at start
        osvalley::Real Z1_pred{};    // OS Z at suggested point (or at flip edge)
        osvalley::Real ratio{};      // |Z0| / |Z1_pred|
        bool flip_ahead{};           // sign flip encountered en route
        int i0{-1}, i1{-1};          // sample indices (start/end)
    };

    
    // ---------------------- assemble F_k and Z_k ----------------------
    struct BlockResult
    {
        Real t0{}, dt{};
        std::vector<Real> t; // K values
        std::vector<Real> Z; // approximate Z(t)
        std::vector<Cpx> F;  // optional: main sum (before θ), length K
    };

    // ---------------------- tiny utilities ----------------------
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

    static inline void fft_inplace(std::vector<Cpx> &a, int sign /*-1: forward, +1: inverse*/)
    {
        // Iterative Cooley–Tukey, power-of-two length, unnormalized forward
        const std::size_t n = a.size();
        assert(n && ((n & (n - 1)) == 0));
        // bit-reverse
        for (std::size_t i = 1, j = 0; i < n; i++)
        {
            std::size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }
        // stages
        for (std::size_t len = 2; len <= n; len <<= 1)
        {
            const Real ang = (sign < 0 ? -1 : +1) * (os_two_pi<Real>() / (Real)len);
            const Cpx wlen = {OsCos(ang), OsSin(ang)};
            for (std::size_t i = 0; i < n; i += len)
            {
                Cpx w{1, 0};
                for (std::size_t j = 0; j < len / 2; ++j)
                {
                    Cpx u = a[i + j];
                    Cpx v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (sign > 0)
        {
            const Real inv = 1 / (Real)n;
            for (auto &x : a)
                x *= inv;
        }
    }

    // θ(t) ≈ (t/2)log(t/2π) - t/2 - π/8 + 1/(48t) + 7/(5760 t^3) + ...
    // θ'(t)= 0.5*log(t/2π) + 1/(48 t^2) + O(t^-4)  (but we prefer exact low orders below)
    static inline ThetaSeries theta_series_at(Real t0)
    {
        const Real two_pi = os_two_pi<Real>();
        // θ0 via 3-term Stirling (more than enough at 1e14)
        Real u = t0 / two_pi;
        Real th0 = 0.5L * t0 * OsLog(u) - 0.5L * t0 - os_pi<Real>() / 8.0L + 1.0L / (48.0L * t0) + 7.0L / (5760.0L * t0 * t0 * t0);
        // exact low derivatives from θ'(t)=1/2 log(t/2π) + O(t^-2)
        Real th1 = 0.5L * OsLog(u);
        Real th2 = 1.0L / (2.0L * t0);
        Real th3 = -1.0L / (2.0L * t0 * t0);
        return {th0, th1, th2, th3};
    }

    static inline Real theta_eval_from_series(const ThetaSeries &S, Real s)
    {
        // s = t - t0
        return S.th0 + S.th1 * s + (S.th2 * (s * s)) / 2.0L + (S.th3 * (s * s * s)) / 6.0L;
    }

    // RS length N ~ floor( sqrt(t/(2π)) )
    static inline std::uint64_t rs_length_N(Real t)
    {
        Real Nf = OsSqrt(t / (os_two_pi<Real>()));
        if (Nf < 1)
            Nf = 1;
        return (std::uint64_t)std::floor(Nf);
    }

    // local zero spacing ~ π / log(t/2π)
    static inline Real zero_spacing_est(Real t)
    {
        return os_pi<Real>() / OsLog(t / (os_two_pi<Real>()));
    }

    // build plan (no n-loop yet)
    static inline OSPlan make_plan(Real t0, std::size_t K_points, const OSConfig &cfg)
    {
        OSPlan P;
        P.t0 = t0;
        P.K = K_points;
        // Δt from spacing / S
        const Real spacing = zero_spacing_est(t0);
        P.dt = spacing / (Real)std::max(1, cfg.samples_per_spacing);

        // pick M so that Δlog = 2π / (M Δt) ≤ delta_log_max and M ≥ K
        const Real need = os_two_pi<Real>() / (P.dt * cfg.delta_log_max); // minimal M (before pow2)
        std::size_t Mmin = (std::size_t)std::ceil(need);
        P.M = next_pow2(std::max<std::size_t>(P.K, Mmin));
        P.dlog = os_two_pi<Real>() / ((Real)P.M * P.dt);

        P.R = std::max(2, cfg.poly_deg_R);
        P.TS = theta_series_at(P.t0);

        P.C.assign(P.R, std::vector<Cpx>(P.M, Cpx{0, 0}));
        P.G.assign(P.R, std::vector<Cpx>(P.M, Cpx{0, 0}));
        P.E_m.resize(P.M);

        // precompute E_m = exp(i t0 m Δlog)
        for (std::size_t m = 0; m < P.M; ++m)
        {
            Real ang = P.t0 * ((Real)m * P.dlog);
            P.E_m[m] = {OsCos(ang), OsSin(ang)};
        }
        return P;
    }

    // ---------------------- build coefficients (O(N*R)) ----------------------
    static inline void accumulate_bins(OSPlan &P, std::uint64_t N)
    {
        const Real inv2 = 0.5L;
        // we’ll round to nearest bin: m = round(log n / Δlog) modulo M
        for (std::uint64_t n = 1; n <= N; ++n)
        {
            Real x = OsLog((Real)n); // log n
            Real mreal = x / P.dlog;
            std::int64_t m = (std::int64_t)llround(mreal);
            // wrap m into [0, M)
            m %= (std::int64_t)P.M;
            if (m < 0)
                m += (std::int64_t)P.M;
            Real delta = x - (Real)m * P.dlog; // δ in ~[-Δlog/2, +Δlog/2]

            // base coefficient: n^{-1/2} * exp( -i t0 * x ) * exp(+i t0 m Δlog)
            Real invsqrt = 1.0L / OsSqrt((Real)n);
            Real a = -P.t0 * x;
            Cpx phase_full = {OsCos(a), OsSin(a)}; // e^{-i t0 log n}
            Cpx c_full = invsqrt * phase_full * P.E_m[(std::size_t)m];

            // b = (-i Δt δ)
            Cpx b(0, -P.dt * delta); // purely imaginary

            // accumulate R terms: c_full * b^r / r!
            Cpx pow = Cpx(1, 0);
            Real invfact = 1.0L; // 1/0! then *= 1/(r+1)
            for (int r = 0; r < P.R; ++r)
            {
                P.C[r][(std::size_t)m] += c_full * (pow * invfact);
                pow *= b;
                invfact /= (Real)(r + 1);
            }
        }
    }

    // ---------------------- FFT per degree ----------------------
    static inline void transform(OSPlan &P)
    {
        for (int r = 0; r < P.R; ++r)
        {
            P.G[r] = P.C[r];         // copy
            fft_inplace(P.G[r], -1); // forward unnormalized
        }
    }

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
            // F_k ≈ sum_{r=0}^{R-1} k^r * G_r[k]
            Cpx Fk(0, 0);
            Real kr = 1;
            for (int r = 0; r < P.R; ++r)
            {
                Fk += kr * P.G[r][k];
                kr *= (Real)k;
            }
            R.F[k] = Fk;

            Real t = P.t0 + (Real)k * P.dt;
            R.t[k] = t;

            if (cfg.return_main_sum)
            {
                R.Z[k] = (Real)Fk.real(); // raw (debug)
            }
            else
            {
                // multiply by e^{i θ(t)} and take 2 Re(...)
                Real s = (Real)k * P.dt; // s = t - t0
                Real th = theta_eval_from_series(P.TS, s);
                Cpx eith(OsCos(th), OsSin(th));
                Real Zk = 2.0L * (Fk * eith).real();
                R.Z[k] = Zk;
            }
        }
        return R;
    }

    inline BlockResult scan_valley_OS(Real t0, const OSConfig &cfg)
    {
        if (t0 < cfg.enable_min_t)
            throw std::runtime_error("OS z disabled for small t0, use 'z' (RS-Gabcke or small z) as param or raise enable_min_t.\n Exiting");

        // choose RS length once for the block
        const std::uint64_t N = rs_length_N(t0);

        OSPlan P = make_plan(t0, cfg.K_points, cfg);
        accumulate_bins(P, N);
        transform(P);

        // quick debug, not used after code refinement, still useful
        //     const Real spacing = zero_spacing_est(t0);
        //     const Real cov_log = (Real)P.M * P.dlog; // 2π/Δt
        //     fprintf(stderr,
        //             "[OS] t0=%.0Lf  N=%llu  M=%zu  R=%d  Δt=%.3Le  spacing≈%.3Le  Δlog=%.3Le  cover_log=%.3Le\n",
        //             t0, (unsigned long long)N, P.M, P.R, (long double)P.dt,
        //             (long double)spacing, (long double)P.dlog, (long double)cov_log);
        
        return evaluate_block(P, cfg);
    }

    
    static inline int sgn(osvalley::Real x) { return (x > 0) - (x < 0); }

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

        const auto K = (int)block.Z.size();
        if (K < 3)
            return out; // not enough samples

        const osvalley::Real dt = block.dt;
        // nearest sample index to t_cur
        int i0 = (int)OsLLRound((t_cur - block.t0) / dt);
        i0 = std::clamp(i0, 1, K - 2); // keep room for central diff
        out.i0 = i0;

        const osvalley::Real Z0 = block.Z[i0];
        out.Z0 = Z0;

        // central difference for dZ/dt
        const osvalley::Real g = (block.Z[i0 + 1] - block.Z[i0 - 1]) / (2 * dt);
        const osvalley::Real g_abs = std::max(std::abs(g), g_floor);

        // target step length for ratio R:  s = (1 - 1/R) * |Z0| / |g|
        osvalley::Real s = (osvalley::Real(1) - osvalley::Real(1) / R_target) * std::abs(Z0) / g_abs;

        // clamp by spacing
        const osvalley::Real s_min = s_min_frac * spacing;
        const osvalley::Real s_max = s_max_frac * spacing;
        s = std::clamp(s, s_min, s_max);

        // descent direction: we want d|Z|/dt * step < 0
        // First-order d|Z|/dt ≈ sign(Z0) * g  => choose step sign = -sign(Z0 * g)
        int dir = -sgn(Z0 * g);
        if (dir == 0)
            dir = -sgn(Z0); // degenerate fallback

        // walk along OS samples to locate flip before reaching s
        int max_hops = (int)std::floor(std::abs(s / dt));
        max_hops = std::max(max_hops, 1);
        int i = i0;
        bool flip = false;
        for (int hop = 0; hop < max_hops; ++hop)
        {
            int inext = i + dir;
            if (inext <= 0 || inext >= K)
                break;
            if (block.Z[i] * block.Z[inext] <= 0)
            {
                // found sign flip; land at the edge (inext)
                i = inext;
                flip = true;
                break;
            }
            i = inext;
        }
        // if no flip, we may still need to “top off” to match s more closely
        if (!flip)
        {
            int i_target = i0 + dir * max_hops;
            i_target = std::clamp(i_target, 0, K - 1);
            i = i_target;
        }

        out.i1 = i;
        out.Z1_pred = block.Z[i];
        out.flip_ahead = flip;

        // observed ratio from OS samples
        const osvalley::Real Z1abs = std::max(std::abs(out.Z1_pred), osvalley::Real(1e-30));
        out.ratio = std::abs(Z0) / Z1abs;

        // suggested Δt (signed)
        out.dt_suggest = (i - i0) * dt;

        return out;
    }

    enum class OSVerdict
    {
        ACCEPT,
        REJECT,
        RESCAN
    };

    struct OSDecision
    {
        OSVerdict verdict{OSVerdict::RESCAN};
        long double t_seed{};   // seed for Newton (flip secant or min-based)
        bool has_flip{};        // FYI
        long double min_absZ{}; // FYI
        long double absZc{};    // FYI
        long double ratio{};    // FYI (|dt_pred|/spacing)
    };

    // helper: trusted sign with hysteresis
    inline int trusted_sign(long double z, long double eta)
    {
        if (z > eta)
            return +1;
        if (z < -eta)
            return -1;
        return 0;
    }

    // median(|Z|) as robust amplitude scale
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

    inline OSDecision os_classify_block(const osvalley::BlockResult &B, const osvalley::OSConfig &oc)
    {
        using osvalley::Real;
        OSDecision out{};

        const int K = (int)B.Z.size();
        if (K < 2)
            return out;
        const int c = K / 2;

        const long double spacing =
            (long double)os_pi<Real>() / OsLog((long double)B.t[c] / (2.0L * (long double)os_pi<Real>()));

        // center stats
        const long double Zc = (long double)B.Z[c];
        const long double absZc = OsAbs(Zc);

        long double dZdt = 0;
        if (K >= 3)
            dZdt = ((long double)B.Z[c + 1] - (long double)B.Z[c - 1]) / (2.0L * (long double)B.dt);

        const long double dt_pred = (dZdt != 0 ? Zc / dZdt : 0);
        const long double ratio = (spacing > 0 ? OsAbs(dt_pred / spacing) : 1e9L);

        // robust guard band for sign trust
        std::vector<long double> zvals(K);
        for (int i = 0; i < K; ++i)
            zvals[i] = (long double)B.Z[i];
        const long double amp = median_abs(zvals);
        const long double eta = std::max((long double)oc.eps_abs, (long double)(oc.eps_rel * amp));

        // scan for trusted flip & locate nearest to center
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

                // secant seed in [i, i+1]
                const long double z0 = zvals[i], z1 = zvals[i + 1];
                const long double t0 = (long double)B.t0 + (long double)i * (long double)B.dt;
                t_seed = t0 - z0 * ((long double)B.dt / (z1 - z0));
                break; // take first; it will be near center in a well-planned block
            }
        }

        // find min |Z| in block (for near-but-no-flip)
        int j = 0;
        for (int i = 1; i < K; ++i)
            if (OsAbs(zvals[i]) < OsAbs(zvals[j]))
                j = i;
        const long double min_absZ = OsAbs(zvals[j]);
        if (!flip)
        {
            // fallback Newton seed at min sample if slope ok
            long double slope = 0;
            if (j > 0 && j + 1 < K)
                slope = (zvals[j + 1] - zvals[j - 1]) / (2.0L * (long double)B.dt);
            else if (j + 1 < K)
                slope = (zvals[j + 1] - zvals[j]) / (long double)B.dt;
            else if (j > 0)
                slope = (zvals[j] - zvals[j - 1]) / (long double)B.dt;
            t_seed = ((long double)B.t0 + (long double)j * (long double)B.dt) -
                     (OsAbs(slope) > 1e-12L ? zvals[j] / slope : 0);
        }

        // classification
        const bool flip_near_center = flip && (OsAbs(((long double)flip_iL + 0.5L - (long double)c) * (long double)B.dt) <= oc.flip_center_frac * spacing);
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

namespace osquick
{

    using osvalley::Real;

    struct OSEval
    {
        Real t_center{};
        Real Zc{}, absZc{};
        Real dZdt{};
        Real dt_pred{};  // Newton offset from center
        Real spacing{};  // ~ pi/log(t/2pi)
        Real ratio{};    // |dt_pred| / spacing
        Real min_absZ{}; // min_i |Z[i]|
        int idx_min{-1}; // where min_absZ occurs
        bool has_flip{false};
        bool near{false};
        Real suggested_step{}; // = -dt_pred
    };

    inline Real local_spacing(Real t)
    {
        // spacing ≈ π / log(t/2π)
        const Real two_pi = (Real)2.0L * os_pi<Real>();
        return os_pi<Real>() / OsLog(t / two_pi);
    }

    inline OSEval analyze_block(const osvalley::BlockResult &b,
                                Real absZ_tol = (Real)0.5L,
                                Real ratio_tol = (Real)0.30L)
    {
        OSEval ev{};
        const int K = (int)b.Z.size();
        const int c = K / 2;

        ev.t_center = b.t[c];
        ev.Zc = b.Z[c];
        ev.absZc = OsAbs((long double)ev.Zc); // safe cast if Real=long double
        ev.spacing = local_spacing(ev.t_center);

        // slope at center (guard edges)
        if (K >= 3)
        {
            ev.dZdt = (b.Z[c + 1] - b.Z[c - 1]) / ((Real)2.0L * b.dt);
        }
        else
        {
            ev.dZdt = (Real)0;
        }

        // one-shot Newton offset (if slope is usable)
        if (ev.dZdt != (Real)0)
        {
            ev.dt_pred = ev.Zc / ev.dZdt;
        }
        else
        {
            ev.dt_pred = (Real)0;
        }
        ev.ratio = (ev.spacing > (Real)0 ? OsAbs((long double)(ev.dt_pred / ev.spacing)) : (Real)1e9L);

        // scan for sign flip & min |Z|
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
            Real a = OsAbs((long double)zi);
            if (a < ev.min_absZ)
            {
                ev.min_absZ = a;
                ev.idx_min = i;
            }
        }

        // verdict: any of these means "near enough"
        ev.near = (ev.absZc <= absZ_tol) || ev.has_flip || (ev.ratio <= ratio_tol) || (ev.min_absZ <= absZ_tol);

        // step suggestion (for your scanner to “home in” quickly)
        ev.suggested_step = -ev.dt_pred;

        return ev;
    }

    inline void print_summary(const OSEval &ev)
    {
        std::cout << std::scientific << std::setprecision(18);
        std::cout << "t_c=" << ev.t_center
                  << "\n|Zc|=" << ev.absZc
                  << "\ndZdt=" << ev.dZdt
                  << "\ndt=" << ev.dt_pred
                  << "\nspacing=" << ev.spacing
                  << "\nratio=" << ev.ratio
                  << "\nmin|Z|=" << ev.min_absZ
                  << "\nflip=" << (ev.has_flip ? 1 : 0)
                  << "\nverdict=" << (ev.near ? "NEAR" : "FAR")
                  << "\nsuggested_step=" << ev.suggested_step
                  << "\n";
    }

    // inline bool os_keep_for_valley(const osquick::OSEval &ev,
    //                                long double tol_refine = 1e-4L,
    //                                long double z_far = 0.5L,
    //                                long double ratio_far = 0.5L,
    //                                long double min_mult = 10.0L)
    // {
    //     // definite keep if we already bracketed a zero in the block
    //     if (ev.has_flip)
    //         return true;

    //     // keep if the best sample in the block is “close-ish”
    //     if (ev.min_absZ <= min_mult * tol_refine)
    //         return true;

    //     // otherwise, if the center looks obviously far and off-center, skip
    //     if (ev.absZc >= z_far && ev.ratio >= ratio_far)
    //         return false;

    //     // uncertain: keep (you can also choose to resample with a smaller block)
    //     return true;
    // }

    // // --- OS sign-change detector & zero seed (for BlockResult) ------------------
    // // Tune from config if you like:
    // struct OSFlipCfg
    // {
    //     long double eps_abs = 1e-3L;   // absolute guard band for "trusting" a sample's sign
    //     long double eps_rel = 0.02L;   // relative guard band (fraction of local amplitude)
    //     long double min_mult = 10.0L;  // accept min|Z| <= min_mult * tol_refine as "near"
    //     long double slope_min = 1e-6L; // avoid divide-by-tiny when estimating from min sample
    // };

    // // Result of scanning one OS block
    // struct OSFlipHit
    // {
    //     bool has_flip = false;               // adjacent trusted samples cross sign
    //     bool near = false;                   // no flip, but min|Z| is small enough
    //     size_t i_left = SIZE_MAX;            // left index of crossing (if has_flip)
    //     long double t_zero = 0;              // secant/interpolated zero seed
    //     long double z_left = 0, z_right = 0; // values around the crossing
    //     size_t i_min = SIZE_MAX;             // index of min|Z| in block
    //     long double zc_abs = 0;              // |Z| at center (optional diagnostic)
    // };

    // // median(|Z|) for a robust local amplitude scale
    // static long double os_median_abs(const std::vector<long double> &Z)
    // {
    //     std::vector<long double> v;
    //     v.reserve(Z.size());
    //     for (auto z : Z)
    //         v.push_back(OsAbs(z));
    //     if (v.empty())
    //         return 0;
    //     auto mid = v.begin() + v.size() / 2;
    //     std::nth_element(v.begin(), mid, v.end());
    //     return *mid;
    // }

    // static inline int trusted_sign(long double z, long double eta)
    // {
    //     if (z > eta)
    //         return +1;
    //     if (z < -eta)
    //         return -1;
    //     return 0; // too close to 0 to trust sign
    // }

    // // Scan an OS block for a sign flip; otherwise produce a good zero seed
    // OSFlipHit os_find_flip(const osvalley::BlockResult &B,
    //                        long double tol_refine, // e.g. 1e-4L (your Newton entry tol)
    //                        const OSFlipCfg &cfg = {})
    // {
    //     OSFlipHit out;

    //     const size_t K = B.Z.size();
    //     if (K < 2)
    //         return out;

    //     // Guard band (hysteresis) for trusting signs
    //     const long double amp = os_median_abs(B.Z);
    //     const long double eta = std::max(cfg.eps_abs, cfg.eps_rel * amp);

    //     // (optional) |Z| at the center sample, handy for logging
    //     const size_t cI = K / 2;
    //     out.zc_abs = OsAbs(B.Z[cI]);

    //     // 1) Look for adjacent trusted sign change and do a linear secant zero
    //     for (size_t i = 0; i + 1 < K; ++i)
    //     {
    //         const long double z0 = B.Z[i];
    //         const long double z1 = B.Z[i + 1];
    //         const int s0 = trusted_sign(z0, eta);
    //         const int s1 = trusted_sign(z1, eta);
    //         if (s0 && s1 && s0 != s1)
    //         {
    //             const long double t0 = B.t0 + (long double)i * B.dt;
    //             // linear crossing in [t_i, t_{i+1}]
    //             const long double t_cross = t0 - z0 * (B.dt / (z1 - z0));

    //             out.has_flip = true;
    //             out.i_left = i;
    //             out.t_zero = t_cross;
    //             out.z_left = z0;
    //             out.z_right = z1;
    //             return out;
    //         }
    //     }

    //     // 2) No trusted flip: use min|Z| sample as "last good approx" to a zero
    //     size_t j = 0;
    //     for (size_t i = 1; i < K; ++i)
    //         if (OsAbs(B.Z[i]) < OsAbs(B.Z[j]))
    //             j = i;

    //     out.i_min = j;

    //     // "near" if min|Z| is within a loose factor of the refine tol
    //     const long double min_absZ = OsAbs(B.Z[j]);
    //     if (min_absZ <= cfg.min_mult * tol_refine)
    //         out.near = true;

    //     // Estimate a zero from min sample via a finite-difference slope
    //     const long double tj = B.t0 + (long double)j * B.dt;
    //     long double slope = 0;
    //     if (j > 0 && j + 1 < K)
    //     {
    //         slope = (B.Z[j + 1] - B.Z[j - 1]) / (2 * B.dt);
    //     }
    //     else if (j + 1 < K)
    //     {
    //         slope = (B.Z[j + 1] - B.Z[j]) / B.dt;
    //     }
    //     else if (j > 0)
    //     {
    //         slope = (B.Z[j] - B.Z[j - 1]) / B.dt;
    //     }
    //     if (OsAbs(slope) > cfg.slope_min)
    //     {
    //         out.t_zero = tj - B.Z[j] / slope; // Newton one-shot from min sample
    //     }
    //     else
    //     {
    //         out.t_zero = tj; // fallback: just use the min sample’s t
    //     }

    //     return out;
    // }

} // namespace osquick
