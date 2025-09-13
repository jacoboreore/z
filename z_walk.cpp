#include <iostream>
#include <mpfr.h>
#include <thread>
#include <mutex>
#include <cstdint>
#include <cmath>
#include "z_calculator.h"
#include "custom_y.h"
#include "settings.h"
#include "logs.h"

std::mutex print_mutex;

// Step adaptation: if z_cur < z_prev → shrink step (refine), else grow (escape).
// Keeps the sign of 'step' so "down" stays negative.
// Bounds keep step in a sane range [1, 1e12] in magnitude.
static inline void smart_step(const mpfr_t z_prev,
                              const mpfr_t z_cur,
                              mpfr_t step,
                              double grow = 10.0,
                              double shrink = 0.5,
                              double step_min = 1.0,
                              double step_max = 1e12)
{
    mpfr_t abs_prev, abs_cur, sgn, tmp;
    mpfr_inits2(CFG().mpfr_prec, abs_prev, abs_cur, sgn, tmp, (mpfr_ptr)0);

    mpfr_abs(abs_prev, z_prev, MPFR_RNDN);
    mpfr_abs(abs_cur, z_cur, MPFR_RNDN);

    // sgn = sign(step) in { -1, 0, +1 }
    int s = mpfr_sgn(step);
    mpfr_set_si(sgn, (s >= 0 ? 1 : -1), MPFR_RNDN);

    // scale = shrink or grow
    const double scale = (mpfr_cmp(abs_cur, abs_prev) < 0) ? shrink : grow;

    // |step| = clamp(|step| * scale, step_min, step_max)
    mpfr_abs(tmp, step, MPFR_RNDN);
    mpfr_mul_d(tmp, tmp, scale, MPFR_RNDN);

    if (mpfr_cmp_d(tmp, step_min) < 0)
        mpfr_set_d(tmp, step_min, MPFR_RNDN);

    if (mpfr_cmp_d(tmp, step_max) > 0)
        mpfr_set_d(tmp, step_max, MPFR_RNDN);

    // step ← sgn * |step|
    mpfr_mul(step, tmp, sgn, MPFR_RNDN);

    mpfr_clears(abs_prev, abs_cur, sgn, tmp, (mpfr_ptr)0);
}

// Walk N downward in big leaps to see if Z keeps improving.
// All MPFR state is local to this thread. No set_default_prec here.
static void backward_search_uj(
    uint64_t N0, mpfr_prec_t P,
    uint64_t steps = 1000000ULL,
    uint64_t stride = 100000ULL)
{
    mpfr_t N, y, z, best_z;
    mpfr_inits2(P, N, y, z, best_z, (mpfr_ptr)0);

    // N ← N0
    if (N0 > (uint64_t)std::numeric_limits<unsigned long>::max())
    {
        // fall back via double if uj isn’t available; or split into loops
        mpfr_set_d(N, (double)N0, MPFR_RNDN);
    }
    else
    {
        mpfr_set_ui(N, (unsigned long)N0, MPFR_RNDN);
    }

    mpfr_set_ui(best_z, 1000, MPFR_RNDN); // Large initial |Z|

    for (uint64_t i = 0; i < steps; ++i)
    {
        // N = N - stride
        mpfr_t s;
        mpfr_init2(s, P);
        if (stride > (uint64_t)std::numeric_limits<unsigned long>::max())
            mpfr_set_d(s, (double)stride, MPFR_RNDN);
        else
            mpfr_set_ui(s, (unsigned long)stride, MPFR_RNDN);
        mpfr_sub(N, N, s, MPFR_RNDN);
        mpfr_clear(s);

        compute_y_from_N(y, N);
        compute_Z(z, y);
        mpfr_abs(z, z, MPFR_RNDN);

        if (mpfr_cmp(z, best_z) < 0)
        {
            mpfr_set(best_z, z, MPFR_RNDN);
            std::lock_guard<std::mutex> lock(print_mutex);
            mpfr_printf("*** BACKWARD: improved at N=%.0Rf, |Z|=%.12Rf, y=%.12Rf\n", N, z, y);
        }
    }

    mpfr_clears(N, y, z, best_z, (mpfr_ptr)0);
}

// int mountain_walk(mpfr_t start)
// {

//     mpfr_set_default_prec(CFG().mpfr_prec);

//     // Initialize big numbers
//     mpfr_t N,
//         calculated_y,
//         corrected_y,
//         z,
//         prev_z,
//         threshold_z;

//     mpfr_inits(
//         N,
//         calculated_y,
//         corrected_y,
//         z,
//         prev_z,
//         threshold_z,
//         (mpfr_ptr)0);

//     // Initialize starting values
//     mpfr_set(N, start, MPFR_RNDN);
//     mpfr_set_ui(prev_z, 0, MPFR_RNDN);
//     mpfr_set_d(threshold_z, 0.5, MPFR_RNDN);

//     long count_attempts = 0;
//     long step = 1;
//     int switch_dir = 1; // 1=up, 0=down

//     while (true)
//     {
//         // Calculate jacob_refined
//         compute_y_from_N(calculated_y, N);
//         compute_Z(z, calculated_y);
//         mpfr_abs(z, z, MPFR_RNDN);

//         if (count_attempts % 10 == 0 && count_attempts != 0)
//         {
//             count_attempts = 0;
//             step *= 100000;
//         }

//         // N = N + step
//         mpfr_t step_mpfr;
//         mpfr_init_set_si(step_mpfr, step, MPFR_RNDN);
//         mpfr_add(N, N, step_mpfr, MPFR_RNDN);
//         mpfr_clear(step_mpfr);

//         if (switch_dir == 1)
//         { // going up
//             if (mpfr_cmp(z, prev_z) > 0)
//             {
//                 count_attempts++;
//             }
//             else
//             {
//                 // gmp_printf("Switching to down at step: %ld, z: %.20Rf\n", step, z);
//                 // mpfr_printf("DOWN: %ld, z: %.20Rf\n", step, z);
//                 switch_dir = 0;
//                 step = 1;
//                 count_attempts = 0;
//             }
//             mpfr_set(prev_z, z, MPFR_RNDN);
//             continue;
//         }

//         if (switch_dir == 0)
//         { // going dow
//             if (mpfr_cmp(z, prev_z) < 0)
//             {
//                 count_attempts++;
//                 mpfr_set(prev_z, z, MPFR_RNDN);
//                 continue;
//             }
//             else
//             {
//                 mpfr_printf("UP AGAIN: %ld, z: %.20Rf y: %.20Rf\n", step, z, calculated_y);
//                 // gmp_printf("Switching again to up at step: %ld, z: %.20Rf\n", step, z);
//                 switch_dir = 1;
//                 step = 1;
//                 count_attempts = 0;
//             }
//         }

//         // Apply correction
//         // y_correction(corrected_y, calculated_y, "newton");
//         // z_abs(z, corrected_y);

//         if (mpfr_cmp(z, threshold_z) < 0)
//         {
//             // gmp_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);
//             //  mpfr_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);
//             mpfr_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);

//             // Launch backward checker in a separate thread
//             mpfr_t N_copy;
//             mpfr_init2(N_copy, CFG().mpfr_prec);
//             mpfr_set(N_copy, N, MPFR_RNDN);

//             std::thread(backward_search, N_copy, CFG().mpfr_prec, 100000000000).detach();
//         }

//         // Reset for next up-direction search
//         switch_dir = 1;
//         step = 1;
//     }

//     // Clean up
//     mpfr_clears(N, calculated_y, corrected_y, z, prev_z, threshold_z, (mpfr_ptr)0);

//     return 0;
// }

// int n_range_walk(mpfr_t min_t, mpfr_t max_t)
// {

//     mpfr_t N,
//         Z,
//         y,
//         ant,
//         best_y,
//         best_z,
//         one;

//     // Convert command-line arguments to mpfr_t for precision
//     mpfr_set_default_prec(CFG().mpfr_prec);
//     int prec = CFG().mpfr_prec;

//     mpfr_inits2(prec,
//                 N,
//                 Z,
//                 y,
//                 ant,
//                 best_y,
//                 best_z,
//                 one,
//                 (mpfr_ptr)0);

//     mpfr_set_ui(one, 1, MPFR_RNDN);    // constant 1
//     mpfr_set_ui(ant, 1000, MPFR_RNDN); // initial best value (arbitrary large)
//     mpfr_set(N, min_t, MPFR_RNDN);

//     mpfr_t threshold;
//     mpfr_init_set_d(threshold, 0.0001, MPFR_RNDN);

//     while (mpfr_cmp(N, max_t) <= 0)
//     {

//         compute_y_from_N(y, N);
//         compute_Z(Z, y);

//         mpfr_abs(Z, Z, MPFR_RNDN);

//         if (mpfr_cmp(Z, threshold) < 0)
//         {
//             mpfr_printf("Best y %.20Rf/%.20Rf/%.20Rf\n", y, Z, N);
//         }

//         // Check if Z < 0.001
//         /*if (mpfr_cmp(Z, threshold) < 0) {
//             //std::cout << std::fixed << std::setprecision(20);
//             if (mpfr_cmp(Z, ant) < 0) {
//                 mpfr_printf("Best N: %.20Rf\n", N);
//                 mpfr_printf("Best Z: %.20Rf\n", best_z);
//                 mpfr_printf("Best y: %.20Rf\n\n", best_y);

//                 mpfr_set(best_z, Z, MPFR_RNDN);
//                 mpfr_set(best_y, y, MPFR_RNDN);
//                 mpfr_set(ant, Z, MPFR_RNDN);

//             }
//         }*/

//         // Increment N by 1: N = N + 1
//         mpfr_add(N, N, one, MPFR_RNDN);
//     }

//     // std::cout << std::fixed << std::setprecision(20);
//     // mpfr_printf("Best N: %.20Rf\n", N);
//     // mpfr_printf("Best Z: %.20Rf\n", best_z);
//     // mpfr_printf("Best y: %.20Rf\n", best_y);

//     mpfr_clears(min_t, max_t, N, Z, y, ant, best_y, best_z, one, (mpfr_ptr)0);

//     return 0;
// }

void y_confirmed(mpfr_t start_y)
{
    mpfr_set_default_prec(CFG().mpfr_prec);

    const int max_iters = 1000000;
    int iter = 0;

    mpfr_t N, z_prev, z_cur, y_next, step, yN, y_best, Np, Nm, yp, ym, zp, zm;
    mpfr_inits2(CFG().mpfr_prec, N, z_prev, z_cur, y_next, step, yN, y_best, Np, Nm, yp, ym, zp, zm, (mpfr_ptr)0);

    // next y by spacing
    //calculate_next_y_from_spacing(y_next, start_y);
    mpfr_set(y_next, start_y, MPFR_RNDN);

    LOGMP("Process params\n---> start_y %.12Rf\nNext y %.12Rf\n", start_y, y_next);

    // |Z| at next y
    compute_Z(z_prev, y_next);
    mpfr_abs(z_prev, z_prev, MPFR_RNDN);
    LOGMP("|Z| at next y: %.12Rf\n", z_prev);

    // N = floor(N(y_next))
    compute_N_from_y(N, y_next);
    mpfr_floor(N, N);
    LOGMP("N from next y: %.0Rf\n", N);

    // Try both directions (N±1) to select initial descent
    mpfr_add_ui(Np, N, 1, MPFR_RNDN);
    compute_y_from_N(yp, Np);
    compute_Z(zp, yp);
    mpfr_abs(zp, zp, MPFR_RNDN);

    mpfr_sub_ui(Nm, N, 1, MPFR_RNDN);
    compute_y_from_N(ym, Nm);
    compute_Z(zm, ym);
    mpfr_abs(zm, zm, MPFR_RNDN);

    // choose direction: go where |Z| is smaller
    bool go_down = (mpfr_cmp(zm, zp) < 0);

    mpfr_set(N, go_down ? Nm : Np, MPFR_RNDN);
    mpfr_set(z_cur, go_down ? zm : zp, MPFR_RNDN);
    mpfr_set(yN, go_down ? ym : yp, MPFR_RNDN);
    mpfr_set(y_best, yN, MPFR_RNDN);

    // initial step: ±1 in the chosen direction
    mpfr_set_si(step, go_down ? -1 : 1, MPFR_RNDN);

    LOGMP("\n\n****** Initial direction: %s ******\n", go_down ? "DOWN" : "UP");

    while (iter++ < max_iters)
    {
        // move N by 'step'
        mpfr_add(N, N, step, MPFR_RNDN);

        // evaluate |Z| at new N
        compute_y_from_N(yN, N);
        compute_Z(z_prev, yN); // reuse z_prev as temp
        mpfr_abs(z_prev, z_prev, MPFR_RNDN);

        // ---- per-step log (Z, N, step) ----
        if (CFG().rs_log_debug)
        {
            //if ((iter % 50) == 0)
            LOGMP("scan: Z=%.12Rf N=%.0Rf step=%.0Rf t=%.20Rf\n", z_prev, N, step, yN);
        }
        // -----------------------------------

        if (mpfr_cmp(z_prev, z_cur) < 0)
        {
            // improvement → record best, tighten step
            mpfr_set(z_cur, z_prev, MPFR_RNDN);
            mpfr_set(y_best, yN, MPFR_RNDN);

            // shrink step magnitude for honing
            smart_step(
                z_prev,
                z_cur,
                step,
                CFG().valley_grow,
                CFG().valley_shrink,
                CFG().valley_step_min,
                CFG().valley_step_max);

            // (optional) occasional refinement log remains if you like:
            // if ((iter % 50) == 0)
            //     mpfr_printf("refine: N=%.0Rf, |Z|=%.12Rf, step=%.0Rf\n", N, z_cur, step);
            continue;
        }
        else
        {
            // got worse → grow step to jump valleys, and flip direction
            smart_step(z_cur, z_prev, step, /*grow=*/10.0, /*shrink=*/0.5);
            mpfr_neg(step, step, MPFR_RNDN); // flip direction once

            // One “bounce” is your valley bottom
            mpfr_printf("\n\n*** valley bottom ~ y=%.25Rf, |Z|=%.25Rf, N=%.0Rf\n", y_best, z_cur, N);
            // mpfr_printf("refine: N=%.0Rf, |Z|=%.12Rf, step=%.0Rf\n", N, z_cur, step);

            // optional: launch backward probe (non-blocking)
            // uint64_t N0 = (uint64_t)std::llround(mpfr_get_d(N, MPFR_RNDN));
            // std::thread(backward_search_uj, N0, CFG().mpfr_prec, 200000ULL, 20000ULL).detach();

            break;
        }
    }

    mpfr_clears(N, z_prev, z_cur, y_next, step, yN, y_best, Np, Nm, yp, ym, zp, zm, (mpfr_ptr)0);

// Use this code for segmentation faults, to debug memory claiming
// #define FREE_MPFR(var)                                           \
//     do                                                           \
//     {                                                            \
//         fprintf(stderr, "[free] about to clear %s ...\n", #var); \
//         fflush(stderr);                                          \
//         mpfr_clear(var);                                         \
//         fprintf(stderr, "[free] cleared %s\n", #var);            \
//         fflush(stderr);                                          \
//     } while (0)

//     // Clear in a deterministic order (reverse order often helps spot issues)
//     FREE_MPFR(zm);
//     FREE_MPFR(zp);
//     FREE_MPFR(ym);
//     FREE_MPFR(yp);
//     FREE_MPFR(Nm);
//     FREE_MPFR(Np);
//     FREE_MPFR(y_best);
//     FREE_MPFR(yN);
//     FREE_MPFR(step);
//     FREE_MPFR(y_next);
//     FREE_MPFR(z_cur);
//     FREE_MPFR(z_prev);
//     FREE_MPFR(N);

// #undef FREE_MPFR
}
