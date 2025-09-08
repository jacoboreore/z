#pragma once

#include <mpfr.h>

void compute_y_from_N(mpfr_t y, mpfr_t N);

void compute_N_from_y(mpfr_t y, mpfr_t N);

void calculate_next_y_from_spacing(mpfr_t y, mpfr_t next);
