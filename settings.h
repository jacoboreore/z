#pragma once
#include <cstdint>

// ======= OS precision backend switch (long double ⇄ MPFR) =======
#if defined(OS_USE_MPFR)
#include <mpreal.h>
#define OS_HAS_MPFR 1
using OS_MPReal = mpfr::mpreal;
// set a sane default; override from your settings before use
inline void os_set_default_prec_bits(int bits) { mpfr::mpreal::set_default_prec(bits); }
#else
#define OS_HAS_MPFR 0
using OS_MPReal = long double;
inline void os_set_default_prec_bits(int) { /* no-op */ }
#endif

// ---- pick Real type (keeps your OS_REAL_T override working) ----
#ifndef OS_REAL_T
#define OS_REAL_T OS_MPReal
#endif

using Real = OS_REAL_T;

struct Config
{
    int mpfr_prec = 250; // precision bits
    
    double prefilter_tau_boost = 0.50;
    double prefilter_alpha = 0.50;
    double prefilter_tau = 0.30;

    double newton_tol_abs = 1e-5;
    double newton_tol_step = 1e-14;
    int newton_max_iter = 5;
    bool use_threaded_z = true;
    bool use_threaded_newton = false;

    double ball_mid_threshold = 1e-8;

    double h_min = 1e-9;             // min FD step
    double h_frac_spacing = 0.10;    // h ≤ h_frac_spacing * spacing
    double step_frac_spacing = 0.70; // |Δt| ≤ step_frac_spacing * spacing

    // RS / threading knobs
    bool rs_half_weight_last = false;  // apply ½‑weight on last term
    double rs_tail_C = 0.5;            // remainder constant (tunable)
    int min_threads = 1;               // 1 = no forced oversubscription
    int max_threads = 0;               // 0 = hw_concurrency
    int rs_min_N_for_threads = 64;     // don’t thread tiny N
    double theta_exact_switch_t = 1e9; // reserved for exact theta gate

    // logging
    bool rs_log_debug = true;
    bool rs_main_log = true;
    int rs_log_digits = 24;

    int theta_stirling_terms = 8; // integer term count

    double valley_grow = 2e-8;            // factor to add to N to walk the to the valley
    double valley_aggressive_factor = 10; // factor to multply if the step to walk fastee
    double valley_tolerance_to_step = 1;
    double valley_max_step = 1e10;
    double valley_max_Zs_relation = 2;

    Real enable_min_t = 1e14L;    // only use OS when t0 >= this
    int samples_per_spacing = 8;  // Δt = spacing / S
    std::size_t K_points = 1024;  // points computed per block (≤ M)
    Real delta_log_max = 1e-3L;   // bin width in log n (smaller → more accurate)
    int poly_deg_R = 4;           // Taylor degree (2..6 reasonable)
    bool return_main_sum = false; // if true: skip θ factor (debug)

    // NEW: classification thresholds (can mirror your valley cfg)
    Real tol_refine = 1e-4L;      // what you pass to Newton as “near enough”
    Real z_far = 0.5L;            // |Zc| above this looks far (order-1)
    Real ratio_far = 0.5L;        // |dt_pred|/spacing above this looks off-center
    Real min_mult = 10.0L;        // accept if min|Z| <= min_mult * tol_refine
    Real eps_abs = 1e-3L;         // absolute guard band for ‘trusted sign’
    Real eps_rel = 0.02L;         // relative guard band (fraction of median |Z|)
    Real flip_center_frac = 0.6L; // flip must be within ±0.6*spacing of center to auto-accept

};

// Initialize once (optional; CFG() will lazy‑init if you don’t call this)
void config_init(const char *path = nullptr);

// Get the global config (thread‑safe, lazy‑initialized once)
const Config &CFG();

// Optional: reload config (e.g., after editing file/env)
void config_reload(const char *path = nullptr);

void api_version();
void print_config(const Config &cfg);
void apply_precision(const Config &cfg);