#pragma once

#include <mpfr.h>

void compute_Z_threaded(mpfr_t Z, mpfr_t t);

void remainder_bound(mpfr_t R, const mpfr_t t); 

