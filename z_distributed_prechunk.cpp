#include <fstream>
#include <vector>
#include <mpfr.h>
#include <mpc.h>
#include <omp.h>
#include <cmath>

#include "settings.h"
#include "logs.h"
#include "utils.h"
#include "z_calculator_rs_gabcke.h"
#include "z_distributed_prechunk.h"

static bool string_to_mpfr(const std::string &s, mpfr_t out)
{
    return mpfr_set_str(out, s.c_str(), 10, MPFR_RNDN) == 0;
}

static int calculate_num_digits(mpfr_prec_t P)
{
    return static_cast<int>(dec_digits_from_bits(P)) + 10;
}

bool write_precalc_file(const std::string &path, const PrecalcData &data)
{
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    if (!f)
        return false;

    f << "N=" << data.N << "\n";
    f << "t=" << data.t_str << "\n";
    f << "theta=" << data.theta_str << "\n";
    f << "gabcke=" << data.gabcke_str << "\n";
    f << "useC1C2=" << data.useC1C2 << "\n";
    f << "prec=" << CFG().mpfr_prec << "\n";

    return true;
}

bool read_precalc_file(const std::string &path, PrecalcData &data)
{
    std::cout << "Reading pre_calc file: " << path << "\n";

    std::ifstream f(path);
    if (!f)
        return false;

    data = PrecalcData{};

    bool hasN = false;
    bool hasT = false;
    bool hasTheta = false;
    bool hasGabcke = false;
    bool hasUseC1C2 = false;
    bool hasPrec = false;

    std::string line;
    try
    {
        while (std::getline(f, line))
        {
            if (line.empty())
                continue;

            auto pos = line.find('=');
            if (pos == std::string::npos || pos == 0 || pos == line.size() - 1)
                continue;

            std::string key = line.substr(0, pos);
            std::string val = line.substr(pos + 1);

            if (key == "N")
            {
                data.N = std::stoll(val);
                hasN = true;
            }
            else if (key == "t")
            {
                data.t_str = val;
                hasT = true;
            }
            else if (key == "theta")
            {
                data.theta_str = val;
                hasTheta = true;
            }
            else if (key == "gabcke")
            {
                data.gabcke_str = val;
                hasGabcke = true;
            }
            else if (key == "useC1C2")
            {
                data.useC1C2 = std::stoi(val);
                hasUseC1C2 = true;
            }
            else if (key == "prec")
            {
                data.prec = static_cast<mpfr_prec_t>(std::stoll(val));
                hasPrec = true;
            }
        }
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Invalid pre_calc file format: " << ex.what() << "\n";
        return false;
    }

    if (!hasN || !hasT || !hasTheta || !hasGabcke || !hasUseC1C2 || !hasPrec)
        return false;

    if (data.N <= 0)
        return false;

    if (data.t_str.empty() || data.theta_str.empty() || data.gabcke_str.empty())
        return false;

    if (data.useC1C2 != 0 && data.useC1C2 != 1)
        return false;

    if (data.prec < MPFR_PREC_MIN)
        return false;

    // Optional strict validation against current runtime config:
    // if (data.prec != CFG().mpfr_prec)
    // {
    //     std::cerr << "Precision mismatch: file=" << data.prec
    //               << ", runtime=" << CFG().mpfr_prec << "\n";
    //     return false;
    // }

    return true;
}

bool compute_precalc_data(const mpfr_t t, PrecalcData &out)
{
    mpfr_prec_t P = CFG().mpfr_prec;
    int digits = calculate_num_digits(P);

    out.t_str = mpfr_to_string(t, digits);

    Constants C(P);
    Spouge S;
    int a = spouge_a_for_prec(P);
    S.build(a, P);

    mpfr_t theta, a_mp, p_mp;
    mpfr_inits2(P, theta, a_mp, p_mp, (mpfr_ptr)0);

    compute_theta_spouge(theta, t, P, S, C);
    std::int64_t N = rs_N_from_t_local(t, P, a_mp, p_mp, C);

    out.N = N;
    out.theta_str = mpfr_to_string(theta, digits);

    // Gabcke correction
    DiffWeights W;
    mpfr_t R1;
    mpfr_init2(R1, P);

    bool useC1C2 = true;
    R_gabcke(R1, t, P, useC1C2, &W, C);

    out.gabcke_str = mpfr_to_string(R1, digits);
    out.useC1C2 = useC1C2 ? 1 : 0;

    mpfr_clear(R1);

    mpfr_clears(theta, a_mp, p_mp, (mpfr_ptr)0);
    return true;
}

bool compute_rs_chunk_from_precalc_serial(
    const PrecalcData &pre,
    std::int64_t n0,
    std::int64_t n1,
    std::string &chunk_sum_str)
{
    if (n0 < 1 || n1 < n0)
        return false;

    mpfr_prec_t P = CFG().mpfr_prec;
    int digits = calculate_num_digits(P);

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

    for (std::int64_t n = n0; n <= n1; ++n)
    {
        mpfr_set_uj(n_mp, static_cast<uintmax_t>(n), MPFR_RNDN);
        mpfr_log(logn, n_mp, MPFR_RNDN);      // log(n)
        mpfr_mul(phi, t, logn, MPFR_RNDN);    // t log n
        mpfr_neg(phi, phi, MPFR_RNDN);        // -t log n
        mpfr_add(phi, phi, theta, MPFR_RNDN); // theta - t log n

        mpfr_cos(cphi, phi, MPFR_RNDN);
        mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
        mpfr_div(cphi, cphi, sqrt_n, MPFR_RNDN);

        mpfr_add(sum, sum, cphi, MPFR_RNDN);
    }

    chunk_sum_str = mpfr_to_string(sum, digits);

    mpfr_clears(t, theta, sum, n_mp, sqrt_n, logn, phi, cphi, (mpfr_ptr)0);
    return true;
}

bool compute_rs_chunk_from_precalc_threaded(
    const PrecalcData &pre,
    std::int64_t n0,
    std::int64_t n1,
    std::string &chunk_sum_str)
{
    if (n0 < 1 || n1 < n0)
        return false;

    mpfr_prec_t P = CFG().mpfr_prec;
    int digits = calculate_num_digits(P);

    mpfr_t t, theta, sum;
    mpfr_inits2(P, t, theta, sum, (mpfr_ptr)0);

    if (!string_to_mpfr(pre.t_str, t))
    {
        mpfr_clears(t, theta, sum, (mpfr_ptr)0);
        return false;
    }

    if (!string_to_mpfr(pre.theta_str, theta))
    {
        mpfr_clears(t, theta, sum, (mpfr_ptr)0);
        return false;
    }

    mpfr_set_zero(sum, 0);

    int max_threads = omp_get_max_threads();
    std::vector<mpfr_t> partials(max_threads);

    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_init2(partials[i], P);
        mpfr_set_zero(partials[i], 0);
    }

#pragma omp parallel num_threads(max_threads)
    {
        int tid = omp_get_thread_num();

        mpfr_t local_sum, n_mp, sqrt_n, logn, phi, cphi;
        mpfr_t t_local, theta_local;

        mpfr_inits2(P, local_sum, n_mp, sqrt_n, logn, phi, cphi, t_local, theta_local, (mpfr_ptr)0);
        mpfr_set_zero(local_sum, 0);

        mpfr_set(t_local, t, MPFR_RNDN);
        mpfr_set(theta_local, theta, MPFR_RNDN);

        std::int64_t total_terms = n1 - n0 + 1;
        std::int64_t block_size = (total_terms + max_threads - 1) / max_threads;

        std::int64_t start_n = n0 + tid * block_size;
        std::int64_t end_n = start_n + block_size - 1;
        if (end_n > n1)
            end_n = n1;

        if (start_n <= n1)
        {
            for (std::int64_t n = start_n; n <= end_n; ++n)
            {
                mpfr_set_uj(n_mp, static_cast<uintmax_t>(n), MPFR_RNDN);
                mpfr_log(logn, n_mp, MPFR_RNDN);            // log(n)
                mpfr_mul(phi, t_local, logn, MPFR_RNDN);    // t log n
                mpfr_neg(phi, phi, MPFR_RNDN);              // -t log n
                mpfr_add(phi, phi, theta_local, MPFR_RNDN); // theta - t log n

                mpfr_cos(cphi, phi, MPFR_RNDN);
                mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
                mpfr_div(cphi, cphi, sqrt_n, MPFR_RNDN);

                mpfr_add(local_sum, local_sum, cphi, MPFR_RNDN);
            }
        }

        mpfr_set(partials[tid], local_sum, MPFR_RNDN);

        mpfr_clears(local_sum, n_mp, sqrt_n, logn, phi, cphi, t_local, theta_local, (mpfr_ptr)0);
    }

    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_add(sum, sum, partials[i], MPFR_RNDN);
        mpfr_clear(partials[i]);
    }

    chunk_sum_str = mpfr_to_string(sum, digits);

    mpfr_clears(t, theta, sum, (mpfr_ptr)0);
    return true;
}

bool compute_rs_chunk_from_precalc_threaded_rot(
    const PrecalcData &pre,
    std::int64_t n0,
    std::int64_t n1,
    std::string &chunk_sum_str)
{
    if (n0 < 1 || n1 < n0)
        return false;

    mpfr_prec_t P = CFG().mpfr_prec;
    int digits = calculate_num_digits(P);

    mpfr_t t, theta, sum;
    mpfr_inits2(P, t, theta, sum, (mpfr_ptr)0);

    if (!string_to_mpfr(pre.t_str, t))
    {
        mpfr_clears(t, theta, sum, (mpfr_ptr)0);
        return false;
    }

    if (!string_to_mpfr(pre.theta_str, theta))
    {
        mpfr_clears(t, theta, sum, (mpfr_ptr)0);
        return false;
    }

    mpfr_set_zero(sum, 0);

    int max_threads = omp_get_max_threads();
    std::vector<mpfr_t> partials(max_threads);
    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_init2(partials[i], P);
        mpfr_set_zero(partials[i], 0);
    }

    std::int64_t total_terms = n1 - n0 + 1;
    std::int64_t thread_span = (total_terms + max_threads - 1) / max_threads;

    // Keep this modest, like original working code
    // TODO: create config for this
    std::int64_t inner_block_size = CFG().dist_rot_inner_block_size;
    if (inner_block_size < 1024)
        inner_block_size = 1024;

#pragma omp parallel num_threads(max_threads)
    {
        int tid = omp_get_thread_num();

        std::int64_t start_n = n0 + tid * thread_span;
        std::int64_t end_n = start_n + thread_span - 1;
        if (end_n > n1)
            end_n = n1;

        mpfr_t local_sum;
        mpfr_t t_local, theta_local;
        mpfr_t n_mp, sqrt_n;
        mpfr_t phi0, cphi, sphi;
        mpfr_t dphi, sd, cd;
        mpfr_t tmp, new_c, new_s;

        mpfr_inits2(P,
                    local_sum,
                    t_local, theta_local,
                    n_mp, sqrt_n,
                    phi0, cphi, sphi,
                    dphi, sd, cd,
                    tmp, new_c, new_s,
                    (mpfr_ptr)0);

        mpfr_set_zero(local_sum, 0);
        mpfr_set(t_local, t, MPFR_RNDN);
        mpfr_set(theta_local, theta, MPFR_RNDN);

        if (start_n <= n1)
        {
            for (std::int64_t block_start = start_n; block_start <= end_n; block_start += inner_block_size)
            {
                std::int64_t block_end = block_start + inner_block_size - 1;
                if (block_end > end_n)
                    block_end = end_n;

                // Re-anchor phase exactly at block_start
                mpfr_set_uj(n_mp, static_cast<uintmax_t>(block_start), MPFR_RNDN);
                mpfr_log(phi0, n_mp, MPFR_RNDN);
                mpfr_mul(phi0, phi0, t_local, MPFR_RNDN);
                mpfr_neg(phi0, phi0, MPFR_RNDN);
                mpfr_add(phi0, phi0, theta_local, MPFR_RNDN);

                mpfr_sin_cos(sphi, cphi, phi0, MPFR_RNDN);

                // contribute n = block_start
                mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
                mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
                mpfr_add(local_sum, local_sum, tmp, MPFR_RNDN);

                // rotate only within this small block
                for (std::int64_t i = block_start; i < block_end; ++i)
                {
                    mpfr_set_uj(n_mp, static_cast<uintmax_t>(i), MPFR_RNDN);
                    mpfr_ui_div(dphi, 1, n_mp, MPFR_RNDN);
                    mpfr_log1p(dphi, dphi, MPFR_RNDN);
                    mpfr_mul(dphi, dphi, t_local, MPFR_RNDN);

                    mpfr_sin_cos(sd, cd, dphi, MPFR_RNDN);

                    // phi -> phi - dphi
                    mpfr_mul(new_c, cphi, cd, MPFR_RNDN);
                    mpfr_mul(tmp, sphi, sd, MPFR_RNDN);
                    mpfr_add(new_c, new_c, tmp, MPFR_RNDN);

                    mpfr_mul(new_s, sphi, cd, MPFR_RNDN);
                    mpfr_mul(tmp, cphi, sd, MPFR_RNDN);
                    mpfr_sub(new_s, new_s, tmp, MPFR_RNDN);

                    mpfr_set(cphi, new_c, MPFR_RNDN);
                    mpfr_set(sphi, new_s, MPFR_RNDN);

                    mpfr_set_uj(n_mp, static_cast<uintmax_t>(i + 1), MPFR_RNDN);
                    mpfr_sqrt(sqrt_n, n_mp, MPFR_RNDN);
                    mpfr_div(tmp, cphi, sqrt_n, MPFR_RNDN);
                    mpfr_add(local_sum, local_sum, tmp, MPFR_RNDN);
                }
            }
        }

        mpfr_set(partials[tid], local_sum, MPFR_RNDN);

        mpfr_clears(local_sum,
                    t_local, theta_local,
                    n_mp, sqrt_n,
                    phi0, cphi, sphi,
                    dphi, sd, cd,
                    tmp, new_c, new_s,
                    (mpfr_ptr)0);
    }

    for (int i = 0; i < max_threads; ++i)
    {
        mpfr_add(sum, sum, partials[i], MPFR_RNDN);
        mpfr_clear(partials[i]);
    }

    chunk_sum_str = mpfr_to_string(sum, digits);

    mpfr_clears(t, theta, sum, (mpfr_ptr)0);
    return true;
}

bool compute_rs_chunk_from_precalc(const PrecalcData &pre,
                                   std::int64_t n0,
                                   std::int64_t n1,
                                   std::string &chunk_sum_str)
{

    switch (CFG().dist_chunk_mode)
    {
    case 1:
        LOGMPF("SERIAL chunk calculations");
        return compute_rs_chunk_from_precalc_serial(pre, n0, n1, chunk_sum_str);
    case 2:
        LOGMPF("PARALLEL chunk calculations");
        return compute_rs_chunk_from_precalc_threaded(pre, n0, n1, chunk_sum_str);
    case 3:
        LOGMPF("PARALLEL ROTATED chunk calculations");
        return compute_rs_chunk_from_precalc_threaded_rot(pre, n0, n1, chunk_sum_str);
    default:
        LOGMPF("Invalid processing mode, defaulting to PARALLEL ROTATED chunk calculations");
        return compute_rs_chunk_from_precalc_threaded_rot(pre, n0, n1, chunk_sum_str);
    }
}

bool append_chunk_result(
    const std::string &output_path,
    int chunk_id,
    std::int64_t n0,
    std::int64_t n1,
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