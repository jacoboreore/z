#include <mpfr.h>
#include <mpc.h>
#include <cmath>
#include "z_calculator_thread.h"
#include "z_calculator_small.h"
#include "config.h"

// ---------------- Unified compute_Z ----------------
void compute_Z(mpfr_t Z, mpfr_t t) {
    
    // choose N
    mpfr_t pi, two_pi, tmp;
    mpfr_inits2(MPFR_PRECISION, pi, two_pi, tmp, (mpfr_ptr)0);
    
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);
    mpfr_div(tmp, t, two_pi, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    
    long N = mpfr_get_si(tmp, MPFR_RNDN);
    
    if (N < 1) N = 1;
    mpfr_clears(pi, two_pi, tmp, (mpfr_ptr)0);

    if (N < 8) {
        compute_Z_small(Z, t);
    } else {
        compute_Z_threaded(Z, t);
    }
}
