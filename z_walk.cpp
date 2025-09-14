#include <iostream>
#include <mpfr.h>
#include <thread>
#include <mutex>
#include <cstdint>
#include <cmath>
#include "z_calculator.h"
#include "custom_y.h"
#include "custom_y.h"
#include "settings.h"
#include "logs.h"
#include "confirm_zero.h"
#include "utils.h"

std::mutex print_mutex;

bool soft_landing = false;

// Step adaptation: if z_cur < z_prev → shrink step (refine), else grow (escape).
// Keeps the sign of 'step' so "down" stays negative.
// Bounds keep step in a sane range [1, 1e12] in magnitude.
static inline void smart_step(const mpfr_t z_prev,
                              const mpfr_t z_cur,
                              mpfr_t smart_step)
{
    double tolerance_to_step = CFG().valley_tolerance_to_step;
    double max_step = CFG().valley_max_step;
    double max_Zs_relation = CFG().valley_max_Zs_relation;

    mpfr_t abs_prev, abs_cur, sgn, tmp, aggressive_factor, ts, ms, grow, mzr;
    mpfr_inits2(CFG().mpfr_prec, abs_prev, abs_cur, sgn, tmp, aggressive_factor, ts, ms, grow, mzr, (mpfr_ptr)0);

    mpfr_abs(abs_prev, z_prev, MPFR_RNDN);
    mpfr_abs(abs_cur, z_cur, MPFR_RNDN);

    mpfr_set_d(aggressive_factor, CFG().valley_aggressive_factor, MPFR_RNDN);
    mpfr_set_d(ms, max_step, MPFR_RNDN);
    mpfr_set_d(ts, tolerance_to_step, MPFR_RNDN);
    mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);
    mpfr_set_d(mzr, max_Zs_relation, MPFR_RNDN);

    // Evaluating current z, if z > 0.1 (conf), we increase the step aggressively
    // z is far from zero, tolerance 0.1
    if (mpfr_cmp(ts, abs_cur) < 0)
    {
        LOGMP_DEBUG("Applying aggressive step");
        if (mpfr_cmp(ms, smart_step) > 0)
            mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
    }
    else
    {
        if (!soft_landing)
        {
            LOGMP_DEBUG("Reseting to smart step...");
            if (mpfr_sgn(smart_step) < 0)
            {
                mpfr_set(smart_step, grow, MPFR_RNDN);
                mpfr_neg(smart_step, smart_step, MPFR_RNDN);
            }
            else
            {
                mpfr_set(smart_step, grow, MPFR_RNDN);
            }
            soft_landing = true;
        }
        else
        {
            mpfr_div(tmp, abs_prev, abs_cur, MPFR_RNDN);
            LOGMP("Z1/Z relation %.20Rf", tmp);

            // if the relation between the previous Z and the current Z is less than mzr, we are walkig too slowly
            if (mpfr_cmp(mzr, tmp) > 0)
            {
                LOGMP("Increasing aggressively the step to a better z", tmp);
                mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
            }
        }
    }

    // if a z 0.1 or less is detected, we land to the valley softly
    mpfr_clears(abs_prev, abs_cur, sgn, tmp, aggressive_factor, ts, ms, grow, mzr, (mpfr_ptr)0);
}

void y_confirmed(mpfr_t start_y)
{
    mpfr_set_default_prec(CFG().mpfr_prec);

    const int max_iters = 1000000;
    int iter = 0;

    mpfr_t N, N_test, z_prev, z_cur, step, grow, t;
    mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, (mpfr_ptr)0);

    // next y by spacing
    // calculate_next_y_from_spacing(y_next, start_y);
    // mpfr_set(y_next, start_y, MPFR_RNDN);

    LOGMP("Process params\n---> start_y %.12Rf\n\n", start_y);

    // |Z| at next y
    compute_Z(z_prev, start_y);
    mpfr_abs(z_prev, z_prev, MPFR_RNDN);
    LOGMP("|Z| at start: %.20Rf", z_prev);

    mpfr_set_d(grow, CFG().valley_grow, MPFR_RNDN);

    compute_N_from_y(N, start_y);
    LOGMP("N from strting y: %.20Rf\n", N);

    mpfr_add(N_test, N, grow, MPFR_RNDN);
    LOGMP("N+1: %.20Rf (grow: %.20Rf)", N_test, grow);

    compute_y_from_N(t, N_test);
    compute_Z(z_cur, t);
    mpfr_abs(z_cur, z_cur, MPFR_RNDN);

    LOGMP("|Z| at N+1: %.20Rf", z_cur);

    // choose direction: go where |Z| is smaller
    bool go_down = (mpfr_cmp(z_cur, z_prev) < 0);

    // initial step: ±1 in the chosen direction
    mpfr_set_d(step, go_down ? CFG().valley_grow : -1 * CFG().valley_grow, MPFR_RNDN);

    LOGMP("\n\n****** Initial direction: %s ******\n", go_down ? "DOWN" : "UP");

    while (iter++ < max_iters)
    {
        // move N by 'step'
        mpfr_add(N, N, step, MPFR_RNDN);

        // evaluate |Z| at new N
        compute_y_from_N(t, N);
        compute_Z(z_cur, t); // reuse z_prev as temp
        mpfr_abs(z_cur, z_cur, MPFR_RNDN);

        // ---- per-step log (Z, N, step) ----
        if (CFG().rs_main_log)
        {
            LOGMP("Scan: Z = %.50Rf t = %.30Rf STEP = %.10Rf", z_cur, t, step);
            LOGMP_DEBUG("N = %.30Rf\n", N);
        }
        // -----------------------------------

        if (mpfr_cmp(z_cur, z_prev) < 0)
        {
            smart_step(z_prev, z_cur, step);
            mpfr_set(z_prev, z_cur, MPFR_RNDN);

            continue;
        }
        else
        {
            // One “bounce” is your valley bottom
            LOGMP("\n\n*** valley bottom, detected sign change, Mountain position going UP |Z| = %.20Rf t = %.20Rf\n", z_cur, t);
            // mpfr_printf("refine: N=%.0Rf, |Z|=%.12Rf, step=%.0Rf\n", N, z_cur, step);

            // optional: launch backward probe (non-blocking)
            // uint64_t N0 = (uint64_t)std::llround(mpfr_get_d(N, MPFR_RNDN));
            // std::thread(backward_search_uj, N0, CFG().mpfr_prec, 200000ULL, 20000ULL).detach();
            // precision_landing_from_best();
            break;
        }
    }

    mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, (mpfr_ptr)0);
}
