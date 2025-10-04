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
#include "refine_newton.h"
#include "refine_newton_thread.h"
#include "file.h"

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
    // |Z| is too far from a zero
    if (mpfr_cmp(ts, abs_cur) < 0)
    {
        LOGMP("Applying aggressive step, |Z| is too far from zero");
        if (mpfr_cmp(ms, smart_step) > 0)
            mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
    }
    else
    {
        if (!soft_landing)
        {
            LOGMP("------> Reseting to smart step for a soft landing...");
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
                LOGMP("Increasing aggressively the step, from STEP %.20Rf", smart_step);
                mpfr_mul(smart_step, smart_step, aggressive_factor, MPFR_RNDN);
                LOGMP("to STEP %.20Rf", smart_step);
            }
        }
    }

    // if a z 0.1 or less is detected, we land to the valley softly
    mpfr_clears(abs_prev, abs_cur, sgn, tmp, aggressive_factor, ts, ms, grow, mzr, (mpfr_ptr)0);
}

void contant_step(mpfr_t constat_step, const mpfr_t N)
{
    double valley_constant_step = CFG().valley_constant_step;
    double valley_constant_step_div = CFG().valley_constant_step_div;

    if (valley_constant_step == 0)
    {
        LOGMP("Evaluating constant step from N = %.*Rg", CFG().rs_log_digits, N);

        mpfr_t calculated_step, divisor;
        mpfr_inits2(CFG().mpfr_prec, calculated_step, divisor, (mpfr_ptr)0);

        mpfr_set_ui(divisor, valley_constant_step_div, MPFR_RNDN);

        mpfr_sqrt(calculated_step, N, MPFR_RNDN);
        mpfr_div(calculated_step, calculated_step, divisor, MPFR_RNDN);
        mpfr_set(constat_step, calculated_step, MPFR_RNDN);

        LOGMP("Constant step = %.12Rf\n\n", constat_step);

        mpfr_clears(calculated_step, divisor, (mpfr_ptr)0);
    }
    else
    {
        mpfr_set_d(constat_step, valley_constant_step, MPFR_RNDN);
        LOGMP("Constant step from log, step =  %.12Rf", constat_step);
    }
}

void from_y(mpfr_t start_y)
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
        LOGMP("Scan: Z = %.50Rf t = %.30Rf STEP = %.10Rf", z_cur, t, step);
        LOGMP_DEBUG("N = %.30Rf\n", N);

        // -----------------------------------
        if (mpfr_cmp(z_cur, z_prev) < 0)
        {
            smart_step(z_prev, z_cur, step);
            mpfr_set(z_prev, z_cur, MPFR_RNDN);

            continue;
        }
        else
        {
            // bounce
            LOGMP("\n\n*** valley bottom, sign change detected, Mountain position going UP |Z| = %.20Rf t = %.20Rf\n", z_cur, t);

            break;
        }
    }

    mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, (mpfr_ptr)0);
}

void from_y_walk(mpfr_t start_y)
{
    mpfr_set_default_prec(CFG().mpfr_prec);
    const int max_iters = 100; // killer switch
    int iter = 0;
    double newton_tol_step = CFG().newton_tol_step;

    LOGMP("Process y walk initial params\n---> start_y %.12Rf\n\n", start_y);

    for (double count_valley_from_y_walk = 0; count_valley_from_y_walk < 5; count_valley_from_y_walk++)
    {
        mpfr_t N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev;
        mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev, (mpfr_ptr)0);

        // next y by spacing, the spacing formula is skipping zeros
        // tested Odlyzco DB
        //      (103793332903, mpf('30607946001.0414392944755272329218235304200151'))
        //      (103793332904, mpf('30607946001.17507333178113001164219442100034168'))   <------ skipped
        calculate_next_y_from_spacing(y_next, start_y);

        mpfr_set(start_y, y_next, MPFR_RNDN);

        LOGMP("\n\n--------------------------------\n--------Starting cycle----------\n--------------------------------\n");
        LOGMP("Next y = %.12Rf\n\n", start_y);

        for (double count_valley_insistence = 0; count_valley_insistence < CFG().valley_insistence; count_valley_insistence++)
        {
            iter = 0; // killer switch

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

            // FIX: going UP to the mountain result in an already calculated zero
            // if (!go_down)
            // {
            //     LOGMP("\n\n Climbing the mountain, we need to go down...\n ");
            //     mpfr_set_d(step, CFG().valley_grow, MPFR_RNDN);
            //     do
            //     {
            //         mpfr_add(N, N, step, MPFR_RNDN);
            //         compute_y_from_N(t, N);
            //         compute_Z(z_cur, t);
            //         mpfr_abs(z_cur, z_cur, MPFR_RNDN);

            //         if (mpfr_cmp(z_cur, z_prev) > 0)
            //         {
            //             smart_step(z_prev, z_cur, step);
            //             mpfr_set(z_prev, z_cur, MPFR_RNDN);
            //         }

            //         go_down = (mpfr_cmp(z_cur, z_prev) < 0);
            //     } while (!go_down && iter++ < max_iters);

            //     LOGMP("\n\nNow going down...");
            // }

            // initial step: ±1 in the chosen direction
            mpfr_set_d(step, go_down ? CFG().valley_grow : -1 * CFG().valley_grow, MPFR_RNDN);

            LOGMP("\n\n****** Initial direction: %s ******\n", go_down ? "DOWN" : "UP");

            iter = 0; // killer switch
            while (iter++ < max_iters)
            {
                // move N by 'step'
                mpfr_add(N, N, step, MPFR_RNDN);

                // store last t
                mpfr_set(t_prev, t, MPFR_RNDN);

                // evaluate |Z| at new N
                compute_y_from_N(t, N);
                compute_Z(z_cur, t);
                mpfr_abs(z_cur, z_cur, MPFR_RNDN);

                // per-step log (Z, N, step)
                LOGMP("Scan: Z = %.50Rf t = %.30Rf STEP = %.10Rf", z_cur, t, step);
                LOGMP_DEBUG("N = %.30Rf\n", N);

                // -----------------------------------
                if (mpfr_cmp(z_cur, z_prev) < 0)
                {
                    // keep walking to the valley
                    smart_step(z_prev, z_cur, step);
                    mpfr_set(z_prev, z_cur, MPFR_RNDN);

                    continue;
                }
                else
                {
                    // Bounce
                    LOGMP("\n\n*** valley bottom, sign change detected, Mountain position going UP |Z| = %.20Rf t = %.20Rf\n", z_cur, t);
                    LOGMP("refining t = %.20Rf...", t_prev);

                    if (CFG().use_threaded_newton)
                    {
                        LOGMP("Using parallelized Newton...");
                        newton_refine_root_thread(
                            t_star,
                            t_prev,
                            CFG().newton_tol_abs,
                            CFG().newton_tol_step,
                            CFG().newton_max_iter,
                            CFG().use_threaded_z ? 1 : 0);
                    }
                    else
                    {
                        LOGMP("Using Single thread Newton...");
                        newton_refine_root(
                            t_star,
                            t_prev,
                            CFG().newton_tol_abs,
                            CFG().newton_tol_step,
                            CFG().newton_max_iter,
                            CFG().use_threaded_z ? 1 : 0);
                    }

                    // set the refine t to the start_y for either, another pass or find the next t
                    mpfr_set(start_y, t_star, MPFR_RNDN);
                    mpfr_set(t, t_star, MPFR_RNDN);

                    break;
                }
            }

            // if we attained a good zero, we don't need to pass again here
            mpfr_t newton_tol_step, z_from_bounce;
            mpfr_inits2(CFG().mpfr_prec, newton_tol_step, z_from_bounce, (mpfr_ptr)0);
            mpfr_set_d(newton_tol_step, CFG().newton_tol_abs, MPFR_RNDN);

            compute_Z(z_from_bounce, start_y);
            mpfr_abs(z_from_bounce, z_from_bounce, MPFR_RNDN);

            LOGMPF("====================> LOG Refined t = %.*Rg |Z| = %.*Rg",
                   CFG().rs_log_digits,
                   start_y,
                   CFG().rs_log_digits,
                   z_from_bounce);

            if (mpfr_cmp(newton_tol_step, z_from_bounce) < 0)
            {
                LOGMP("No need to pass the valley scanner again...");
                count_valley_insistence = CFG().valley_insistence + 1; // force exiting
            }
            else
            {
                LOGMP("Passing the valley scanner again...");
            }
            mpfr_clears(newton_tol_step, z_from_bounce, (mpfr_ptr)0);
        }

        mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, y_next, t_star, t_prev, (mpfr_ptr)0);
    }
}

void from_y_walk_pure_valley_walk(mpfr_t start_y, const double steps)
{
    mpfr_set_default_prec(CFG().mpfr_prec);
    const int max_iters = 10000;
    int iter = 0;
    int main_iter = 0;
    bool ok = false;
    // double newton_tol_step = CFG().newton_tol_step;

    // get the file name for this process
    std::string process_file = zcsv::make_process_filename("output", "valley_walk");

    // get the file name for the valley walking calculation for plotting
    std::string walk_file = zcsv::make_process_filename("output", "valley_walk_plot");

    LOGMPF("Process y walk initial params\n---> start_y %.12Rf\n---> Valley count %f\n\n******************\n", start_y, steps);

    while (main_iter++ < steps)
    {
        // no overhead here, just restarting everything
        mpfr_t N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce;
        mpfr_inits2(CFG().mpfr_prec, N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce, (mpfr_ptr)0);

        LOGMP("\n\n--------------------------------\n--------Starting cycle----------\n--------------------------------\n");
        // LOGMP("Next y = %.12Rf\n\n", start_y);

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

        contant_step(step, N);

        iter = 0;
        if (!go_down)
        {
            LOGMP("\n\nClimbing the mountain, we need to go down...\n ");
            do
            {
                mpfr_add(N, N, step, MPFR_RNDN);
                compute_y_from_N(t, N);
                compute_Z(z_cur, t);
                mpfr_abs(z_cur, z_cur, MPFR_RNDN);

                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t, z_cur, N, nullptr, 0, 2, CFG().rs_log_digits);
                    if (!ok)
                    {
                        LOGMPF("\n\n******** ERROR writing csv ***********\n\n");
                    }
                }
                if (mpfr_cmp(z_cur, z_prev) > 0)
                {
                    // smart_step(z_prev, z_cur, step);
                    // contant_step(step);
                    mpfr_set(z_prev, z_cur, MPFR_RNDN);
                }

                go_down = (mpfr_cmp(z_cur, z_prev) < 0);
            } while (!go_down && iter++ < max_iters);
        }

        LOGMP("\nGoing down...\n");
        // initial step: ±1 in the chosen direction
        // mpfr_set_d(step, go_down ? CFG().valley_grow : -1 * CFG().valley_grow, MPFR_RNDN);

        // LOGMP("\n\n****** Initial direction: %s ******\n", go_down ? "DOWN" : "UP");

        iter = 0; // killer switch

        while (iter++ < max_iters)
        {
            // move N by 'step'
            mpfr_add(N, N, step, MPFR_RNDN);

            // store last t
            mpfr_set(t_prev, t, MPFR_RNDN);

            // evaluate |Z| at new N
            compute_y_from_N(t, N);
            compute_Z(z_cur, t);
            mpfr_abs(z_cur, z_cur, MPFR_RNDN);

            // per-step log (Z, N, step)
            LOGMP("Scan: Z = %.50Rf t = %.30Rf STEP = %.10Rf", z_cur, t, step);
            LOGMP("N = %.30Rf\n", N);

            // -----------------------------------
            if (mpfr_cmp(z_cur, z_prev) < 0)
            {
                // keep walking to the valley
                mpfr_set(z_prev, z_cur, MPFR_RNDN);
                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t, z_cur, N, nullptr, 0, 2, CFG().rs_log_digits);
                    if (!ok)
                    {
                        LOGMPF("\n\n******** ERROR writing csv ***********\n\n");
                    }
                }
                continue;
            }
            else
            {
                // Bounce
                LOGMP("\n\n*** valley bottom, sign change detected, Mountain position going UP |Z| = %.20Rf t = %.20Rf\n", z_cur, t);
                LOGMP("refining t = %.20Rf...", t_prev);

                // values for plotting, we register a minima at t_prev
                if (CFG().walking_for_plotting_enabled)
                {
                    ok = zcsv::append_record_csv(walk_file, t_prev, z_prev, N, nullptr, 1, 2, CFG().rs_log_digits);
                    if (!ok)
                    {
                        LOGMPF("\n\n******** ERROR writing csv ***********\n\n");
                    }
                }
                if (CFG().use_threaded_newton)
                {
                    LOGMP("Using parallelized Newton...");
                    newton_refine_root_thread(
                        t_refined,
                        t_prev,
                        CFG().newton_tol_abs,
                        CFG().newton_tol_step,
                        CFG().newton_max_iter,
                        CFG().use_threaded_z ? 1 : 0);
                }
                else
                {
                    LOGMP("Using Single thread Newton...");
                    newton_refine_root(
                        t_refined,
                        t_prev,
                        CFG().newton_tol_abs,
                        CFG().newton_tol_step,
                        CFG().newton_max_iter,
                        CFG().use_threaded_z ? 1 : 0);
                }

                compute_Z(z_from_bounce, t_refined);
                mpfr_abs(z_from_bounce, z_from_bounce, MPFR_RNDN);

                // zcsv::append_record_csv(walk_file, t_refined, z_from_bounce, N, CFG().rs_log_digits);

                LOGMPF("====================> LOG Refined t = %.*Rg |Z| = %.*Rg",
                       CFG().rs_log_digits,
                       t_refined,
                       CFG().rs_log_digits,
                       z_from_bounce);

                ok = zcsv::append_record_csv(process_file, t_refined, z_from_bounce, N, nullptr, 2, 1, CFG().rs_log_digits);
                if (!ok)
                {
                    LOGMPF("\n\n******** ERROR writing csv ***********\n\n");
                }

                mpfr_set(start_y, t, MPFR_RNDN);

                break;
            }
        }

        mpfr_clears(N, N_test, z_prev, z_cur, step, grow, t, t_refined, t_prev, z_from_bounce, (mpfr_ptr)0);
    }
}
