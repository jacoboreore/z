#pragma once
#include "cli.h"

int mountain_walk(mpfr_t start);

void from_y(mpfr_t start_y);

void from_y_walk(mpfr_t start_y);

void from_y_walk_pure_valley_walk(mpfr_t start_y, const double steps, Options opt);

int n_range_walk(mpfr_t min_t, mpfr_t max_t);
