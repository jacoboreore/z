#pragma once

#include <mpfr.h>

void compute_Z(mpfr_t Z, mpfr_t t);

bool z_prefilter(const mpfr_t t, double alpha, double tau);
