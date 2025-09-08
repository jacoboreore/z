#pragma once

#include <mpfr.h>
#include <mpc.h>
#include <iostream>

static inline void log_mpfr(
    const char* tag, 
    const mpfr_t x, 
    int digits=30) {

    //std::cerr << tag << " = ";
    //mpfr_out_str(stderr, 10, digits, x, MPFR_RNDN);
    //std::cerr << "\n";
}