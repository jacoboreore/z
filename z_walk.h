#pragma once

// Precision in bits (~77 decimal places for 256)
constexpr int PRECISION = 256;

int mountain_walk(mpfr_t start);

int y_confirmed(mpfr_t start_y);

int n_range_walk(mpfr_t min_t, mpfr_t max_t);
