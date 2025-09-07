
#include <iostream>
#include <mpfr.h>
#include <cstring> // For strcmp
#include "z_walk.h"
#include "config.h"
#include "z_calculator.h"


int main(int argc, char* argv[]) {
    mpfr_set_default_prec(MPFR_PRECISION);

    // Iterate through and print each argument
    std::cout << "Arguments:" << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << "  argv[" << i << "]: " << argv[i] << std::endl;
    }

    mpfr_t t, Z;
    mpfr_inits2(MPFR_PRECISION, t, Z, (mpfr_ptr)0);

    // First zero (sanity check)
    mpfr_set_str(t, argv[1], 10, MPFR_RNDN);
    mpfr_printf("%.10RF\n", t); // %.10RF prints the number with 10 digits after the decimal point
    compute_Z(Z, t);
    mpfr_printf("Z(%.15Rf) = %.20Re\n", t, Z);

    // Your big y sample
    //mpfr_set_str(t, "4400961943970.5696075625388437339818624451523392290921263943249583984239027506119881", 50, MPFR_RNDN);
    //compute_Z(Z, t);
    //mpfr_printf("Z(%.15Rf) = %.10Re\n", t, Z);

    mpfr_clears(t, Z, (mpfr_ptr)0);
    return 0;
}