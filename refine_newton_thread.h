// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.//
//
// File: refine_newton_thread.h
// Purpose:
//   Thread-aware Newton–secant–Pegasus hybrid refiner for Hardy’s Z(t) with
//   MPFR-backed memoization and careful OpenMP usage.
// ============================================================================

#pragma once
#include <mpfr.h>

bool newton_refine_root_thread(
    mpfr_t t_out,
    const mpfr_t t0,
    double tol_abs,
    double tol_step,
    int max_iter,
    int use_threaded_z);


