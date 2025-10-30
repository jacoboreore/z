// -----------------------------------------------------------------------------
// Module: z_calculator_rs_gabcke.h
// Description:
//   - High-precision Hardy Z(t) via the Riemann–Siegel main sum plus
//     Gabcke-style remainder (first correction term R₁ with optional C₁/C₂).
//   - Fast θ(t) using log Γ via Spouge with reflection (stable for Im z large).
//   - Parallel rotation kernel for RS main sum (OpenMP), with dynamic blocking.
//
// Public entry points:
//   - void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t)
// -----------------------------------------------------------------------------

#pragma once

#include <mpfr.h>

void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t);
