#include <iostream>
#include <mpfr.h>
#include <thread>
#include <mutex>
#include "z_calculator.h"
#include "custom_y.h"
#include "config.h"

std::mutex print_mutex;

void smart_step(mpfr_t z1, mpfr_t z2, mpfr_t N, mpfr_t step){

     mpfr_t 
        one,
        ten,
        correction_percentaje,
        bad_threshold;

    mpfr_inits(
        one,
        ten,
        correction_percentaje,
        bad_threshold,
        (mpfr_ptr)0);

    // init 
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_ui(ten, 10, MPFR_RNDN);

    // far from a zero, we need to increase the step 
    if(mpfr_cmp(z2, one) > 0){
        mpfr_mul(step, step, ten, MPFR_RNDN);
    } else {
        mpfr_set_d(bad_threshold, 0.5, MPFR_RNDN);
        if(mpfr_cmp(z2, bad_threshold) > 0){
            mpfr_mul(step, step, ten, MPFR_RNDN);
        }
    }

    mpfr_clears(
        one,
        ten,
        correction_percentaje,
        bad_threshold,
        (mpfr_ptr)0);

}

void backward_search(mpfr_t original_N, int precision, float steps = 10000000000000) {
    mpfr_set_default_prec(precision);

    mpfr_t N, y, z, best_z;
    mpfr_inits2(precision, N, y, z, best_z, (mpfr_ptr)0);
    
    mpfr_set(N, original_N, MPFR_RNDN);
    mpfr_set_ui(best_z, 1000, MPFR_RNDN); // Large initial Z

    for (int i = 1; i <= steps; i+=100000000000) {
        mpfr_sub_ui(N, N, 1, MPFR_RNDN);  // N = N - 1

        compute_y_from_N(y, N);
        compute_Z(z, y);
        mpfr_abs(z, z, MPFR_RNDN);

        if (mpfr_cmp(z, best_z) < 0) {
            mpfr_set(best_z, z, MPFR_RNDN);
            std::lock_guard<std::mutex> lock(print_mutex);
            //mpfr_printf("THREAD: Improved Z found at N: %.0Rf, z: %.20Rf, y: %.20Rf\n", N, z, y);
        }
    }

    mpfr_printf("*** THREAD: Improved Z found at N: %.0Rf, z: %.20Rf, y: %.20Rf\n", N, z, y);

    mpfr_clears(N, y, z, best_z, (mpfr_ptr)0);
}

int mountain_walk(mpfr_t start) {
    
    mpfr_set_default_prec(MPFR_PRECISION);

    // Initialize big numbers
    mpfr_t  N, 
            calculated_y, 
            corrected_y, 
            z, 
            prev_z, 
            threshold_z;

    mpfr_inits(
        N, 
        calculated_y, 
        corrected_y, 
        z, 
        prev_z, 
        threshold_z, 
        (mpfr_ptr)0);

    // Initialize starting values
    mpfr_set(N, start, MPFR_RNDN);
    mpfr_set_ui(prev_z, 0, MPFR_RNDN);
    mpfr_set_d(threshold_z, 0.5, MPFR_RNDN);

    long count_attempts = 0;
    long step = 1;
    int switch_dir = 1; // 1=up, 0=down

    while (true) {
        // Calculate jacob_refined
        compute_y_from_N(calculated_y, N);
        compute_Z(z, calculated_y);
        mpfr_abs(z, z, MPFR_RNDN);

        if (count_attempts % 10 == 0 && count_attempts != 0) {
            count_attempts = 0;
            step *= 100000;
        }

        // N = N + step
        mpfr_t step_mpfr;
        mpfr_init_set_si(step_mpfr, step, MPFR_RNDN);
        mpfr_add(N, N, step_mpfr, MPFR_RNDN);
        mpfr_clear(step_mpfr);

        if (switch_dir == 1) { // going up
            if (mpfr_cmp(z, prev_z) > 0) {
                count_attempts++;
            } else {
                //gmp_printf("Switching to down at step: %ld, z: %.20Rf\n", step, z);
                //mpfr_printf("DOWN: %ld, z: %.20Rf\n", step, z);
                switch_dir = 0;
                step = 1;
                count_attempts = 0;
            }
            mpfr_set(prev_z, z, MPFR_RNDN);
            continue;
        }

        if (switch_dir == 0) { // going dow
            if (mpfr_cmp(z, prev_z) < 0) {
                count_attempts++;
                mpfr_set(prev_z, z, MPFR_RNDN);
                continue;
            } else {
                mpfr_printf("UP AGAIN: %ld, z: %.20Rf y: %.20Rf\n", step, z, calculated_y);
                //gmp_printf("Switching again to up at step: %ld, z: %.20Rf\n", step, z);
                switch_dir = 1;
                step = 1;
                count_attempts = 0;
            }
        }

        // Apply correction
        //y_correction(corrected_y, calculated_y, "newton");
        //z_abs(z, corrected_y);

        if (mpfr_cmp(z, threshold_z) < 0) {
            //gmp_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);
            // mpfr_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);
            mpfr_printf("Found near-zero! y: %.20Rf, N: %.0Rf, z: %.20Rf\n", calculated_y, N, z);

            // Launch backward checker in a separate thread
            mpfr_t N_copy;
            mpfr_init2(N_copy, MPFR_PRECISION);
            mpfr_set(N_copy, N, MPFR_RNDN);

            std::thread(backward_search, N_copy, MPFR_PRECISION, 100000000000).detach();
        } 
        

        // Reset for next up-direction search
        switch_dir = 1;
        step = 1;
    }

    // Clean up
    mpfr_clears(N, calculated_y, corrected_y, z, prev_z, threshold_z, (mpfr_ptr)0);

    return 0;
}

int n_range_walk(mpfr_t min_t, mpfr_t max_t ){
    
    mpfr_t  N, 
            Z, 
            y, 
            ant, 
            best_y, 
            best_z, 
            one;

    // Convert command-line arguments to mpfr_t for precision
    mpfr_set_default_prec(MPFR_PRECISION);
    int prec = MPFR_PRECISION;

	mpfr_inits2(prec, 
        N, 
        Z, 
        y, 
        ant, 
        best_y, 
        best_z, 
        one, 
        (mpfr_ptr) 0);

    mpfr_set_ui(one, 1, MPFR_RNDN); // constant 1
    mpfr_set_ui(ant, 1000, MPFR_RNDN); // initial best value (arbitrary large)
    mpfr_set(N, min_t, MPFR_RNDN);

    mpfr_t threshold;
    mpfr_init_set_d(threshold, 0.0001, MPFR_RNDN);

    while (mpfr_cmp(N, max_t) <= 0) {

        compute_y_from_N(y, N);
        compute_Z(Z, y);

        mpfr_abs(Z, Z, MPFR_RNDN);

        if (mpfr_cmp(Z, threshold) < 0) {
            mpfr_printf("Best y %.20Rf/%.20Rf/%.20Rf\n", y, Z, N);
        }

           // Check if Z < 0.001
        /*if (mpfr_cmp(Z, threshold) < 0) {
            //std::cout << std::fixed << std::setprecision(20);
            if (mpfr_cmp(Z, ant) < 0) {
                mpfr_printf("Best N: %.20Rf\n", N);
                mpfr_printf("Best Z: %.20Rf\n", best_z);
                mpfr_printf("Best y: %.20Rf\n\n", best_y);

                mpfr_set(best_z, Z, MPFR_RNDN);
                mpfr_set(best_y, y, MPFR_RNDN);
                mpfr_set(ant, Z, MPFR_RNDN);
            
            }
        }*/

        // Increment N by 1: N = N + 1
        mpfr_add(N, N, one, MPFR_RNDN);
    }

    //std::cout << std::fixed << std::setprecision(20);
    //mpfr_printf("Best N: %.20Rf\n", N);
    //mpfr_printf("Best Z: %.20Rf\n", best_z);
    //mpfr_printf("Best y: %.20Rf\n", best_y);

    mpfr_clears(min_t, max_t, N, Z, y, ant, best_y, best_z, one, (mpfr_ptr) 0);

    return 0;
}

int y_confirmed(mpfr_t start_y){

    mpfr_set_default_prec(MPFR_PRECISION);
    int prec = MPFR_PRECISION;
    int steps = 1000000, cont_steps = 0;
    int switch_dir; // 1=up, 0=down

    // Initialize big numbers
    mpfr_t  N, 
            z1,
            z2,
            next,
            step,
            calculated_y,
            stored_y_before_switch_direction,
            minus_step;

    mpfr_inits(
        N, 
        z1,
        z2,
        next,
        step,
        calculated_y,
        stored_y_before_switch_direction,
        minus_step,
        (mpfr_ptr)0);

    // init 
    mpfr_set_ui(step, 1, MPFR_RNDN); // constant 1
    mpfr_set_si(minus_step, -1, MPFR_RNDN); // constant 1

    // log process parameters
    mpfr_printf("Process params\n--->start_y %.10Rf\n", start_y);

    // calculate the next y
    calculate_next_y_from_spacing(start_y, next);
    mpfr_printf("Next y %.10Rf\n", next);
    
    compute_Z(z1, next);
    mpfr_abs(z1, z1, MPFR_RNDN);
    mpfr_printf("Z from Next y %.10Rf\n", z1);
    
    // Initializing Mountain WALK
    compute_N_from_y(next, N); 
    mpfr_floor(N, N);
    mpfr_printf("N from Next y %.10Rf\n", N);
    
    mpfr_add(N, N, step, MPFR_RNDN);
    compute_y_from_N(calculated_y, N);
    compute_Z(z2, calculated_y);
    mpfr_abs(z2, z2, MPFR_RNDN);
    mpfr_printf("z from the Next N+1 %.10Rf %.10Rf\n", N, z2);

    if (mpfr_cmp(z2, z1) > 0){
        switch_dir = 1; // going up
        mpfr_printf("STEPPING BACK\n");
    } else {
        switch_dir = 0; // going down
        mpfr_printf("MOVING FORWARD\n");
    }

    mpfr_printf("Going %d z1 %.20Rf z2 %.10Rf\n", switch_dir, z1, z2);
    
    while (cont_steps < steps) {

        if (switch_dir == 1) {
            mpfr_add(N, N, minus_step, MPFR_RNDN);
        }
        else{
            mpfr_add(N, N, step, MPFR_RNDN);
        }
        
        cont_steps++;

        // Calculate jacob_refined
        compute_y_from_N(calculated_y, N);
        compute_Z(z1, calculated_y);
        //compute_Z(z1, calculated_y, MPFR_PRECISION);
        mpfr_abs(z1, z1, MPFR_RNDN);

        if (mpfr_cmp(z1, z2) < 0){
            if(cont_steps % 10 == 0){
                if (switch_dir == 1) {
                    smart_step(z1, z2, N, minus_step);
                    mpfr_printf("y %.10Rf/%.10Rf/%.10Rf\n", z1, calculated_y, minus_step);
                }
                else{
                    smart_step(z1, z2, N, step);
                    mpfr_printf("y %.10Rf/%.10Rf/%.10Rf\n", z1, calculated_y, step);
                }
            }    
            
            mpfr_set(z2, z1, MPFR_RNDN);
            mpfr_set(stored_y_before_switch_direction, calculated_y, MPFR_RNDN);
            
        }
        else {
            mpfr_printf("*** Best y %.10Rf/%.10Rf/%.10Rf\n", stored_y_before_switch_direction, z1, N);
            break;
        }
    }

    // Clean up
    mpfr_clears(
        N, 
        z1,
        z2, 
        next,
        step,
        minus_step,
        calculated_y,
        stored_y_before_switch_direction, (mpfr_ptr)0);

    return 0;
}
