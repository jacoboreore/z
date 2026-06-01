// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: refine_newton.h
// Purpose:
//   Newton–secant–Pegasus hybrid root refiner for Hardy’s Z(t).
//   - Tracks best |Z| via observer callbacks.
//   - Uses adaptive finite-difference derivatives (3-pt/5-pt).
//   - Backtracking with descent nudges; optional Pegasus burst if bracketed.
// -----------------------------------------------------------------------------

#pragma once
#include <mpfr.h>

bool newton_refine_root(
    mpfr_t t_out, const mpfr_t t0,
    double tol_abs, double tol_step,
    int max_iter, int use_threaded_z);
