// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: custom_y.h
// Purpose:
//   Mountain-Valley core calculations, transform any given N to t using the 
//   documented formula Real-Complex conversion at m = 1/2.
//
// -----------------------------------------------------------------------------

#pragma once

#include <mpfr.h>

void compute_y_from_N(mpfr_t y, mpfr_t N);

void compute_N_from_y(mpfr_t y, mpfr_t N);

void calculate_next_y_from_spacing(mpfr_t y, mpfr_t next);
