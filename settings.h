// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: settings.h
// Purpose:
//   Configuration loader and reporter for the Z(t) computation suite.
//   - Parses key–value pairs from a config file and environment variables.
//   - Exposes a single global Config instance via CFG().
//   - Applies MPFR precision globally.
//   - Prints library version information (GMP/MPFR/MPC).
//
// -----------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <cstring>
#include <string>

#if defined(OS_USE_MPFR)
#include "mpreal.h"
#define OS_HAS_MPFR 1
using OS_MPReal = mpfr::mpreal;
inline void os_set_default_prec_bits(int bits) { mpfr::mpreal::set_default_prec(bits); }
#else
using OS_MPReal = long double;
inline void os_set_default_prec_bits(int) { /* no-op */ }
#endif

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

    long double newton_tol_abs = 1e-5;
    long double newton_tol_step = 1e-14;
    int newton_max_iter = 5;
    bool use_threaded_z = true;
    bool use_threaded_newton = false;

    double ball_mid_threshold = 1e-8;

    double h_min = 1e-9;
    double h_frac_spacing = 0.10;
    double step_frac_spacing = 0.70;

    // RS / threading knobs
    bool rs_half_weight_last = false;
    double rs_tail_C = 0.5;
    int min_threads = 1;
    int max_threads = 0;
    int rs_min_N_for_threads = 64;
    double theta_exact_switch_t = 1e9;

    // logs
    bool rs_log_debug = true;
    bool rs_main_log = true;
    int rs_log_digits = 24;

    int theta_stirling_terms = 8;

    double valley_grow = 2e-8;
    double valley_aggressive_factor = 10;
    double valley_constant_step = 0;
    double valley_constant_step_div = 50;
    double valley_tolerance_to_step = 1;
    double valley_max_step = 1e10;
    double valley_max_Zs_relation = 2;
    double valley_insistence = 3;

    Real enable_min_t = 1e14L;
    int samples_per_spacing = 8;
    std::size_t K_points = 1024;
    Real delta_log_max = 1e-3L;
    int poly_deg_R = 4;
    bool return_main_sum = false;

    Real tol_refine = 1e-4L;
    Real z_far = 0.5L;
    Real ratio_far = 0.5L;
    Real min_mult = 10.0L;
    Real eps_abs = 1e-3L;
    Real eps_rel = 0.02L;
    Real flip_center_frac = 0.6L;

    // Outputs
    bool walking_output_enabled = true;
    bool walking_for_plotting_enabled = true;

    std::string header_refined = "t,absZ,spacing";
    std::string header_walking = "t,absZ,N";
    std::string header_walking_for_plotting = "t,absZ,N,minima";
};

// Initialize once
void config_init(const char *path = nullptr);

// Get the global config (thread‑safe, lazy‑initialized once)
const Config &CFG();

void config_reload(const char *path = nullptr);

void api_version();
void print_config(const Config &cfg);
void apply_precision(const Config &cfg);