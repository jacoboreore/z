
// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: z_calculator.h
// Description:
//   Unified Z(t) evaluation, ball-certification helper, and a lightweight
//   prefilter. Keeps multiple historical/experimental pathways as testimony
//   of versions tested during development.
//
// -----------------------------------------------------------------------------


#pragma once
#include <mpfr.h>

void compute_Z(mpfr_t Z, mpfr_t t);

void compute_Z_ball(mpfr_t Z, mpfr_t t);

bool z_prefilter(const mpfr_t t);

long rs_N(const mpfr_t t);
