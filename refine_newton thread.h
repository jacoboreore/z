#pragma once
#include <mpfr.h>

// Refine a root of Z(t)=0 near the initial guess t0 using Newton-Raphson.
// On success, writes the refined root to t_out and returns true.
// Tolerances:
//  - tol_abs:  stop if |Z(t)| <= tol_abs
//  - tol_step: stop if |Δt|  <= tol_step
// max_iter: hard cap on iterations (e.g. 20)
// use_threaded_z: if non-zero, use compute_Z_threaded where beneficial.
//
// NOTE: This is line-searched Newton with 3/5-point symmetric derivatives.
// It auto-chooses small-Z vs threaded depending on t.
// bool newton_refine_root(
//     mpfr_t t_out,
//     const mpfr_t t0,
//     double tol_abs,
//     double tol_step,
//     int max_iter,
//     int use_threaded_z);

bool newton_refine_root_thread(
    mpfr_t t_out,
    const mpfr_t t0,
    double tol_abs,
    double tol_step,
    int max_iter,
    int use_threaded_z);

// // Coarse-to-fine multi-pass wrapper (mirrors your Python idea)
// bool refine_in_stages(
//     mpfr_t t_out, const mpfr_t t0,
//     int use_threaded_z);
