// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: z_small.h
// Description:
//   - Stable η(s) (Dirichlet eta) evaluation via binomial/Euler transform.
//   - Small-t Z(t) path using ζ(s) from η(s) and θ(t) from Stirling wrapper.
//   - Preserves experimental structure for reproducibility/documentation.
//
// Public entry points:
//   - void compute_Z_small(mpfr_t Z, const mpfr_t t)
//
// -----------------------------------------------------------------------------

#pragma once

#include <mpfr.h>

void compute_Z_small(mpfr_t Z, const mpfr_t t);
