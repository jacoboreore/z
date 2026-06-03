// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// Module: theta_stirling.h
// Description:
//   θ(t) evaluation header utilities for the Riemann–Siegel theta function.
//   Provides Stirling-based approximations with fixed K terms or with
//   an optimally truncated Bernoulli series. The public entry points are:
//     - compute_theta_stirling_K(theta, t, K)
//     - compute_theta_stirling(theta, t)           // uses CFG().theta_stirling_terms
//     - compute_theta_stirling_opt(theta, t)       // optimal truncation strategy
// -----------------------------------------------------------------------------

#pragma once
#include "settings.h"
#include <mpfr.h>

void compute_theta_stirling(mpfr_t theta, const mpfr_t t);

void compute_theta_stirling_K(mpfr_t theta, const mpfr_t t, int K);

void compute_theta_stirling_opt(mpfr_t theta, mpfr_srcptr t);

// no used but implemented for t values below 2000
// void compute_theta_lanczos(mpfr_t theta, const mpfr_t t);