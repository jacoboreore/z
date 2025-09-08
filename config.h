#ifndef CONFIG_H
#define CONFIG_H

// Set the precision in bits (~77 decimal places for 256)
constexpr mpfr_prec_t  MPFR_PRECISION = 100;

// Threshold for detecting near-zero
constexpr double Z_THRESHOLD = 0.001;

// Step adjustment factor
constexpr int STEP_INCREASE_FACTOR = 1000;

// Initial step
constexpr long INITIAL_STEP = 1;

#endif // CONFIG_H