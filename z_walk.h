// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// General Description:
// Core valley-walk routines declarations to locate neighborhoods of zeros of the Riemann
// zeta function on the critical line. This module adapts the sampling step
// based on |Z| descent, optionally refines candidates via Newton, and emits
// CSV records for downstream plotting/analysis.
// -----------------------------------------------------------------------------

#pragma once
#include "cli.h"

int mountain_walk(mpfr_t start);

void from_y(mpfr_t start_y);

void from_y_walk(mpfr_t start_y);

void from_y_walk_pure_valley_walk(mpfr_t start_y, const double steps, Options opt);

int n_range_walk(mpfr_t min_t, mpfr_t max_t);
