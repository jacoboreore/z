#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <cstdio>

#include <mpfr.h>
#include <mpc.h>

#include "settings.h"
#include "logs.h"
#include "utils.h"
#include "z_calculator_rs_gabcke.h"
#include "z_distributed_prechunk.h"

// Reuse from your existing codebase
struct Constants;
struct Spouge;

// ---------- small helpers ----------

// static std::string mpfr_to_string(const mpfr_t x, int digits = 0)
// {
//     char *buf = nullptr;
//     if (digits <= 0)
//         digits = 50;
//     mpfr_asprintf(&buf, "%.*Re", digits, x);
//     std::string s(buf ? buf : "");
//     mpfr_free_str(buf);
//     return s;
// }

static bool string_to_mpfr(const std::string &s, mpfr_t out)
{
    return mpfr_set_str(out, s.c_str(), 10, MPFR_RNDN) == 0;
}

bool write_precalc_file(const std::string &path, const PrecalcData &data)
{
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    if (!f)
        return false;

    f << "N=" << data.N << "\n";
    f << "prec=" << data.prec << "\n";
    f << "t=" << data.t_str << "\n";
    f << "theta=" << data.theta_str << "\n";
    return true;
}

bool read_precalc_file(const std::string &path, PrecalcData &data)
{

    std::cout << "Reading pre_calc file: " << path;

    std::ifstream f(path);
    if (!f)
        return false;

    std::string line;
    while (std::getline(f, line))
    {
        auto pos = line.find('=');
        if (pos == std::string::npos)
            continue;

        std::string key = line.substr(0, pos);
        std::string val = line.substr(pos + 1);

        if (key == "N")
            data.N = std::stol(val);
        else if (key == "prec")
            data.prec = static_cast<mpfr_prec_t>(std::stol(val));
        else if (key == "t")
            data.t_str = val;
        else if (key == "theta")
            data.theta_str = val;
    }

    return true;
}

bool compute_precalc_data(const mpfr_t t, PrecalcData &out)
{
    mpfr_prec_t P = CFG().mpfr_prec;
    out.prec = P;
    out.t_str = mpfr_to_string(t, 80);

    Constants C(P);
    Spouge S;
    int a = spouge_a_for_prec(P);
    S.build(a, P);

    mpfr_t theta, a_mp, p_mp;
    mpfr_inits2(P, theta, a_mp, p_mp, (mpfr_ptr)0);

    compute_theta_spouge(theta, t, P, S, C);
    long N = rs_N_from_t_local(t, P, a_mp, p_mp, C);

    out.N = N;
    out.theta_str = mpfr_to_string(theta, 80);

    mpfr_clears(theta, a_mp, p_mp, (mpfr_ptr)0);
    return true;
}

bool compute_rs_chunk_from_precalc(
    const PrecalcData &pre,
    long n0,
    long n1,
    std::string &chunk_sum_str)
{
    if (n0 < 1 || n1 < n0)
        return false;

    mpfr_prec_t P = pre.prec;

    mpfr_t t, theta, sum, n_mp, sqrt_n, logn, phi, cphi;
    mpfr_inits2(P, t, theta, sum, n_mp, sqrt_n, logn, phi, cphi, (mpfr_ptr)0);

    if (!string_to_mpfr(pre.t_str, t))
    {
        mpfr_clears(t, theta, sum, n_mp, sqrt_n, logn, phi, cphi, (mpfr_ptr)0);
        return false;
    }

    if (!string_to_mpfr(pre.theta_str, theta))
    {
        mpfr_clears(t, theta, sum, n_mp, sqrt_n, logn, phi, cphi, (mpfr_ptr)0);
        return false;
    }

    mpfr_set_zero(sum, 0);

    for (long n = n0; n <= n1; ++n)
    {
        mpfr_set_si(n_mp, n, MPFR_RNDN);
        mpfr_log(logn, n_mp, MPFR_RNDN);      // log(n)
        mpfr_mul(phi, t, logn, MPFR_RNDN);    // t log n
        mpfr_neg(phi, phi, MPFR_RNDN);        // -t log n
        mpfr_add(phi, phi, theta, MPFR_RNDN); // theta - t log n

        mpfr_cos(cphi, phi, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
        mpfr_div(cphi, cphi, sqrt_n, MPFR_RNDN);

        mpfr_add(sum, sum, cphi, MPFR_RNDN);
    }

    chunk_sum_str = mpfr_to_string(sum, 80);

    mpfr_clears(t, theta, sum, n_mp, sqrt_n, logn, phi, cphi, (mpfr_ptr)0);
    return true;
}

bool append_chunk_result(
    const std::string &output_path,
    int chunk_id,
    long n0,
    long n1,
    const std::string &chunk_sum_str)
{
    std::ofstream f(output_path, std::ios::out | std::ios::app);
    if (!f)
        return false;

    f << "chunk_id=" << chunk_id
      << ",n0=" << n0
      << ",n1=" << n1
      << ",sum=" << chunk_sum_str
      << "\n";

    return true;
}