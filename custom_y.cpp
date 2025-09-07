#include <iostream>
#include <mpfr.h>
#include <cmath>
#include <iomanip>
#include "config.h"

void compute_y_from_N(mpfr_t y, mpfr_t N) {
	
    int prec = MPFR_PRECISION;

	mpfr_t cumulative, minus_one;
	mpfr_init2(cumulative, prec);
    mpfr_init2(minus_one, prec);

	mpfr_set_si(minus_one, -1, MPFR_RNDN);
	
	// sqrt(4N - 1) / 2
	mpfr_mul_ui(cumulative, N, 4, MPFR_RNDN);
	mpfr_add(cumulative, cumulative, minus_one, MPFR_RNDN);
	mpfr_sqrt(cumulative, cumulative, MPFR_RNDN);
	mpfr_div_ui(y, cumulative, 2, MPFR_RNDN);

	//std::cout << "calculated y ";
	//mpfr_printf("%.15Rf\n", y);

	mpfr_clears(cumulative, minus_one, (mpfr_ptr) 0);
	
}

void calculate_next_y_from_spacing(mpfr_t y, mpfr_t next) {
    int prec = MPFR_PRECISION;

    mpfr_t pi, two_pi, y_div_twopi, log_term, spacing;

    mpfr_inits2(prec, pi, two_pi, y_div_twopi, log_term, spacing, (mpfr_ptr) 0);

    // two_pi = 2 * pi
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN);

    // y_div_twopi = y / (2 * pi)
    mpfr_div(y_div_twopi, y, two_pi, MPFR_RNDN);

    // log_term = log(y / (2 * pi))
    mpfr_log(log_term, y_div_twopi, MPFR_RNDN);

    // spacing = (2 * pi) / log_term
    mpfr_div(spacing, two_pi, log_term, MPFR_RNDN);

    // next = y + spacing
    mpfr_add(next, y, spacing, MPFR_RNDN);

    // Optional print
    // mpfr_printf("Next y ≈ %.20Rf\n", next);

    mpfr_clears(pi, two_pi, y_div_twopi, log_term, spacing, (mpfr_ptr) 0);
}


void compute_N_from_y(mpfr_t y, mpfr_t N) {
	
    int prec = MPFR_PRECISION;

	mpfr_t cumulative, one, four;
	mpfr_init2(cumulative, prec);
    mpfr_init2(one, prec);
    mpfr_init2(four, prec);
    

	mpfr_set_si(one, 1, MPFR_RNDN);
    mpfr_set_si(four, 4, MPFR_RNDN);
	
	// (4 * y **2 + 1) / 4
	mpfr_mul(cumulative, y, y, MPFR_RNDN);
	mpfr_mul(cumulative, four, cumulative, MPFR_RNDN);
	mpfr_add(cumulative, cumulative, one, MPFR_RNDN);
	mpfr_div(N, cumulative, four, MPFR_RNDN);

	//mpfr_printf("%.15Rf\n", N);

	mpfr_clears(cumulative, one, y, (mpfr_ptr) 0);
	
}