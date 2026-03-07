// -----------------------------------------------------------------------------
// Module: z_calculator_rs_gabcke.h
// Description:
//   - High-precision Hardy Z(t) via the Riemann–Siegel main sum plus
//     Gabcke-style remainder (first correction term R₁ with optional C₁/C₂).
//   - Fast θ(t) using log Γ via Spouge with reflection (stable for Im z large).
//   - Parallel rotation kernel for RS main sum (OpenMP), with dynamic blocking.
//
// Public entry points:
//   - void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t)
// -----------------------------------------------------------------------------

#pragma once

#include <mpfr.h>
#include <vector>

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

/**
 * @brief Container for frequently used constants (π, 2π, π², π⁴, log π)
 *        at a given MPFR precision.
 *
 * Reinitialize via init() when precision changes.
 */
struct Constants
{
    mpfr_t pi, two_pi, pi2, pi4, log_pi;
    mpfr_prec_t P;

    /// Construct with optional precision; real initialization via init().
    explicit Constants(mpfr_prec_t prec = 0) : P(0)
    {
        mpfr_init2(pi, 2);
        mpfr_init2(two_pi, 2);
        mpfr_init2(pi2, 2);
        mpfr_init2(pi4, 2);
        mpfr_init2(log_pi, 2);
        if (prec)
            init(prec);
    }

    /**
     * @brief Initialize cached constants to precision @p prec (idempotent).
     * @param prec MPFR precision in bits.
     */
    void init(mpfr_prec_t prec)
    {
        if (P == prec)
            return;
        P = prec;
        mpfr_set_prec(pi, P);
        mpfr_set_prec(two_pi, P);
        mpfr_set_prec(pi2, P);
        mpfr_set_prec(pi4, P);
        mpfr_set_prec(log_pi, P);

        mpfr_const_pi(pi, MPFR_RNDN);          // π
        mpfr_mul_ui(two_pi, pi, 2, MPFR_RNDN); // 2π
        mpfr_mul(pi2, pi, pi, MPFR_RNDN);      // π²
        mpfr_mul(pi4, pi2, pi2, MPFR_RNDN);    // π⁴
        mpfr_log(log_pi, pi, MPFR_RNDN);       // log π
    }

    /// Destroy and clear internal MPFR variables.
    ~Constants()
    {
        mpfr_clear(pi);
        mpfr_clear(two_pi);
        mpfr_clear(pi2);
        mpfr_clear(pi4);
        mpfr_clear(log_pi);
    }
};

// -----------------------------------------------------------------------------
// Spouge coefficients (real) for log Γ on the right half-plane
// -----------------------------------------------------------------------------

/**
 * @brief Spouge series coefficient table for log Γ, right half-plane.
 *
 * Stores coefficients in raw __mpfr_struct (POD) to avoid extra indirection.
 */
struct Spouge
{
    int a;
    std::vector<__mpfr_struct> c; ///< c[0] = sqrt(2π); c_k for k∈[1,a-1]
    mpfr_t sqrt2pi;

    /// Construct with default a=0 and initialize sqrt(2π) at moderate precision.
    Spouge() : a(0)
    {
        mpfr_init2(sqrt2pi, 128);
        mpfr_const_pi(sqrt2pi, MPFR_RNDN);
        mpfr_mul_ui(sqrt2pi, sqrt2pi, 2, MPFR_RNDN);
        mpfr_sqrt(sqrt2pi, sqrt2pi, MPFR_RNDN);
    }

    /// Clear internal MPFR variables and coefficient storage.
    ~Spouge()
    {
        for (auto &ci : c)
            mpfr_clear(&ci);
        mpfr_clear(sqrt2pi);
    }

    /**
     * @brief Build Spouge coefficient table for parameter @p A at precision @p P.
     *        Coefficients c_k for k≥1 are computed in parallel.
     * @param A Spouge parameter a (> 1).
     * @param P MPFR precision in bits.
     */
    void build(int A, mpfr_prec_t P)
    {
        a = A;

        // Re-init c with correct precision
        for (auto &ci : c)
            mpfr_clear(&ci);
        c.clear();
        c.resize(a);
        for (int k = 0; k < a; ++k)
            mpfr_init2(&c[k], P);

        // c0 = sqrt(2π)
        mpfr_t c0;
        mpfr_init2(c0, P);
        mpfr_const_pi(c0, MPFR_RNDN);
        mpfr_mul_ui(c0, c0, 2, MPFR_RNDN);
        mpfr_sqrt(c0, c0, MPFR_RNDN);
        mpfr_set(&c[0], c0, MPFR_RNDN);
        mpfr_clear(c0);

// Crunch c_k in parallel. thread-local temporaries.
#pragma omp parallel
        {
            mpfr_t t1, t2, t3;
            mpfr_inits2(P, t1, t2, t3, (mpfr_ptr)0);

#pragma omp for
            for (int k = 1; k < a; ++k)
            {
                // sign flips every k
                int sgn = (k & 1) ? +1 : -1;

                // (a-k)^(k-1/2)
                mpfr_set_si(t1, a - k, MPFR_RNDN);
                mpfr_set_si(t2, k, MPFR_RNDN);
                mpfr_sub_d(t2, t2, 0.5, MPFR_RNDN);
                mpfr_pow(t1, t1, t2, MPFR_RNDN);

                // * e^(a-k)
                mpfr_set_si(t2, a - k, MPFR_RNDN);
                mpfr_exp(t2, t2, MPFR_RNDN);
                mpfr_mul(t1, t1, t2, MPFR_RNDN);

                if (sgn < 0)
                    mpfr_neg(t1, t1, MPFR_RNDN);

                // divide by (k-1)!
                mpfr_fac_ui(t2, (unsigned long)(k - 1), MPFR_RNDN);
                mpfr_div(&c[k], t1, t2, MPFR_RNDN);
            }

            mpfr_clears(t1, t2, t3, (mpfr_ptr)0);
        }
    }
};

void compute_Z_rs_gabcke(mpfr_t Z, const mpfr_t t);
void compute_theta_spouge(mpfr_t theta, const mpfr_t t, mpfr_prec_t P, const Spouge &S, const Constants &C);
int spouge_a_for_prec(mpfr_prec_t P);
