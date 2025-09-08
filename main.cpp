#include <mpfr.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <utility>
//#include "z_walk.h"
#include "config.h"
#include "z_calculator.h"

// int main(int argc, char* argv[]) {
//     mpfr_set_default_prec(MPFR_PRECISION);

//     // Iterate through and print each argument
//     std::cout << "Arguments:" << std::endl;
//     for (int i = 0; i < argc; ++i) {
//         std::cout << "  argv[" << i << "]: " << argv[i] << std::endl;
//     }

//     mpfr_t t, Z;
//     mpfr_inits2(MPFR_PRECISION, t, Z, (mpfr_ptr)0);

//     // First zero (sanity check)
//     mpfr_set_str(t, argv[1], 10, MPFR_RNDN);
//     mpfr_printf("%.10RF\n", t); // %.10RF prints the number with 10 digits after the decimal point
//     compute_Z(Z, t);
//     mpfr_printf("Z(%.15Rf) = %.20Re\n", t, Z);

//     // Your big y sample
//     //mpfr_set_str(t, "4400961943970.5696075625388437339818624451523392290921263943249583984239027506119881", 50, MPFR_RNDN);
//     //compute_Z(Z, t);
//     //mpfr_printf("Z(%.15Rf) = %.10Re\n", t, Z);

//     mpfr_clears(t, Z, (mpfr_ptr)0);
//     return 0;
// }

// helper: mpfr -> double
static inline double dval(const mpfr_t x) { return mpfr_get_d(x, MPFR_RNDN); }

struct Case { const char* label; const char* t_str; };

static void run_suite(double alpha, double tau,
                      double near_eps = 1e-3, double far_min = 0.5)
{
    printf("=== z_prefilter tests  alpha=%.3f  tau=%.3f  (near_eps=%.1e, far_min=%.2f)\n",
           alpha, tau, near_eps, far_min);

    // Mix of near-zeros and clearly-off points
    // First few nontrivial zeros (approx):
    // 14.1347251417, 21.0220396388, 25.0108575801, 30.4248761259, 32.9350615877, 37.5861781588
    std::vector<Case> cases = {
        {"zero-ish", "14.134725141734693"},
        {"off",      "14.20"},
        {"zero-ish", "21.022039638771555"},
        {"off",      "21.10"},
        {"zero-ish", "25.010857580145689"},
        {"off",      "25.20"},
        {"zero-ish", "30.424876125859513"},
        {"off",      "30.60"},
        {"zero-ish", "32.935061587739190"},
        {"off",      "33.10"},
        {"zero-ish", "37.586178158825671"},
        {"off",      "37.80"},
        // some mid t where RS sum is nontrivial
        {"mid",      "400"},
        {"mid",      "500"},
        {"mid",      "800"},
        {"mid",      "1000"}
    };

    int tp=0, fn=0, tn=0, fp=0, amb=0;

    mpfr_t t, Z; mpfr_inits2(256, t, Z, (mpfr_ptr)0);

    for (auto &c : cases) {
        mpfr_set_str(t, c.t_str, 10, MPFR_RNDN);

        // Prefilter decision
        bool keep = z_prefilter(t, alpha, tau);

        // Ground truth via small-Z
        compute_Z(Z, t);
        double z  = dval(Z);
        double az = std::fabs(z);

        // classify: near / far / ambiguous (for evaluation only)
        enum { NEAR, FAR, AMB } clazz = AMB;
        if (az <= near_eps) clazz = NEAR;
        else if (az >= far_min) clazz = FAR;
        else clazz = AMB;

        // update counters
        if (clazz == NEAR) { if (keep) ++tp; else ++fn; }
        else if (clazz == FAR) { if (keep) ++fp; else ++tn; }
        else { ++amb; }

        printf("t=%-18s  |Z|=%-10.6g  prefilter=%-5s  class=%s\n",
               c.t_str, az, keep ? "keep" : "drop",
               clazz==NEAR ? "near" : (clazz==FAR ? "far" : "amb"));
    }

    mpfr_clears(t, Z, (mpfr_ptr)0);

    printf("\nConfusion (excluding 'amb'):\n");
    printf("  TP(near & keep)   = %d\n", tp);
    printf("  FN(near & drop)   = %d\n", fn);
    printf("  TN(far  & drop)   = %d\n", tn);
    printf("  FP(far  & keep)   = %d\n", fp);

    int tot = tp+tn+fp+fn;
    if (tot > 0) {
        double recall = (tp+fn) ? double(tp)/double(tp+fn) : 1.0;
        double precision = (tp+fp) ? double(tp)/double(tp+fp) : 1.0;
        printf("  recall (zeros caught)   = %.3f\n", recall);
        printf("  precision (keeps legit) = %.3f\n", precision);
    }
    printf("  amb (ignored in stats)  = %d\n", amb);
    printf("===============================================================\n\n");
}

int main() {
    // deterministic precision
    mpfr_set_default_prec(256);

    // Try a couple of (alpha, tau) settings
    run_suite(/*alpha*/0.50, /*tau*/0.30);
    run_suite(/*alpha*/0.65, /*tau*/0.30);
    run_suite(/*alpha*/0.50, /*tau*/0.20);

    return 0;
}