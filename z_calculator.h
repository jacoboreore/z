#pragma once

#include "settings.h"
#include <mpfr.h>

void compute_Z(mpfr_t Z, mpfr_t t);

void compute_Z_ball(mpfr_t Z, mpfr_t t);

bool z_prefilter(const mpfr_t t, double alpha, double tau);

static long rs_N(const mpfr_t t);
