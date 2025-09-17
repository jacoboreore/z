#pragma once
#include <mpfr.h>

// Threaded RS evaluator: Z(t) = 2 * Σ_{n=1..N} cos(θ(t) - t log n)/√n  (+ optional ½-weight last)
void compute_Z_threaded(mpfr_t Z, mpfr_t t);

// Configurable analytic remainder bound |R(t)|
void remainder_bound(mpfr_t R, const mpfr_t t);
