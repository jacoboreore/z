#pragma once


struct Config {
    int mpfr_prec = 250; // bits

    double prefilter_alpha = 0.50;
    double prefilter_tau = 0.30;

    double newton_tol_abs = 1e-5;
    double newton_tol_step = 1e-14;
    int newton_max_iter = 5;
    bool use_threaded_z = true;

    double ball_mid_threshold = 1e-8;

    double h_min = 1e-9; // min FD step
    double h_frac_spacing = 0.10; // h ≤ h_frac_spacing * spacing
    double step_frac_spacing = 0.70; // |Δt| ≤ step_frac_spacing * spacing

    // RS / threading knobs
    bool rs_half_weight_last = false; // apply ½‑weight on last term
    double rs_tail_C = 0.5; // remainder constant (tunable)
    int min_threads = 1; // 1 = no forced oversubscription
    int max_threads = 0; // 0 = hw_concurrency
    int rs_min_N_for_threads = 64; // don’t thread tiny N
    double theta_exact_switch_t = 1e9; // reserved for exact theta gate


    // logging
    bool rs_log_debug = true;
    int rs_log_digits = 24;


    int theta_stirling_terms = 8; // integer term count
};


// Initialize once (optional; CFG() will lazy‑init if you don’t call this)
void config_init(const char* path = nullptr);


// Get the global config (thread‑safe, lazy‑initialized once)
const Config& CFG();


// Optional: reload config (e.g., after editing file/env)
void config_reload(const char* path = nullptr);


void api_version();
void print_config(const Config& cfg);
void apply_precision(const Config& cfg);