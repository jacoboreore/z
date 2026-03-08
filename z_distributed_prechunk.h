#pragma once

#include <mpfr.h>
#include <string>

struct PrecalcData
{
    std::int64_t N;
    std::string theta_str;
    std::string t_str;
    std::string gabcke_str;
    int useC1C2;
    int prec;
};

// -----------------------------------------------------------------------------
// RS parameterization
// -----------------------------------------------------------------------------

/**
 * @brief Compute RS parameters from t:
 *        a = sqrt(t/(2π)), N = floor(a), p = frac(a).
 *
 * @param t      Height parameter.
 * @param P      MPFR precision in bits.
 * @param a_out  Output a.
 * @param p_out  Output p.
 * @param C      Cached constants.
 * @return N = floor(a), clamped to at least 1.
 */
static std::int64_t rs_N_from_t_local(const mpfr_t t, mpfr_prec_t P, mpfr_t a_out, mpfr_t p_out, const Constants &C)
{
    mpfr_t tmp;
    mpfr_init2(tmp, P);

    mpfr_div(tmp, t, C.two_pi, MPFR_RNDN);
    mpfr_sqrt(a_out, tmp, MPFR_RNDN); // a = sqrt(t/(2π))

    std::int64_t N = mpfr_get_si(a_out, MPFR_RNDN);

    if (N < 1)
        N = 1;

    mpfr_set_uj(tmp, static_cast<uintmax_t>(N), MPFR_RNDN);
    mpfr_sub(p_out, a_out, tmp, MPFR_RNDN); // p = a - floor(a)
    mpfr_clear(tmp);

    return N;
}

bool write_precalc_file(const std::string &path, const PrecalcData &data);
bool read_precalc_file(const std::string &path, PrecalcData &data);

// Computes theta(t), N and writes to pre_cal.txt-compatible structure.
bool compute_precalc_data(const mpfr_t t, PrecalcData &out);

// Computes S_chunk = sum_{n=n0..n1} cos(theta - t log n)/sqrt(n)
bool compute_rs_chunk_from_precalc(
    const PrecalcData &pre,
    std::int64_t n0,
    std::int64_t n1,
    std::string &chunk_sum_str);

// Append chunk result to output file.
bool append_chunk_result(
    const std::string &output_path,
    int chunk_id,
    std::int64_t n0,
    std::int64_t n1,
    const std::string &chunk_sum_str);