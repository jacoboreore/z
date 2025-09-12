#pragma once
#include "settings.h"
#include <mpfr.h>

// Riemann–Siegel theta function θ(t)
// Implementation: real Stirling/Atkinson series
// θ(t) = t/2·log(t/2π) − t/2 − π/8 + Σ_{k≥1} a_{2k-1}/t^{2k-1}
// with a1=1/48, a3=7/5760, a5=31/80640, a7=127/430080, a9=511/1216512 …
// K = CFG().theta_stirling_terms (clamped to available coefficients)
// Valid for t ≥ 0; θ(0) = −π/8.
void compute_theta_stirling(mpfr_t theta, const mpfr_t t);

void compute_theta_stirling_K(mpfr_t theta, const mpfr_t t, int K);
