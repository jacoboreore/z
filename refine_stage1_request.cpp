// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// File: refine_stage1_request.cpp
// Purpose:
//   First-stage distributed refiner request generator.
//   Given (job_id, t0), compute initial sampling points for remote Z(t)
//   evaluation and write them to a JSON file for Python/auditing.
//
// Usage:
//   refine_stage1_request <job_id> <t0> <output_json>
//
// Example:
//   refine_stage1_request job-001 123456789.123 /tmp/refine_stage1_job-001.json
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <mpfr.h>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "settings.h"
#include "logs.h"
#include "utils.h"
#include "confirm_zero.h"

// ----------------------------------------------------------------------------
// Small helpers
// ----------------------------------------------------------------------------

static inline double choose_h_stage1(const mpfr_t t)
{
    // Same philosophy as current refiner:
    // h scales with precision and local magnitude, but is clamped by spacing.
    const double s = est_spacing(t);
    const double eps = std::ldexp(1.0, -(int)CFG().mpfr_prec); // 2^(-P)
    const double mag = std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
    const double h0 = std::pow(eps, 1.0 / 3.0) * mag;

    return std::clamp(h0, 1e-9, 0.1 * s);
}

static std::string json_escape(const std::string &s)
{
    std::ostringstream oss;
    for (char c : s)
    {
        switch (c)
        {
        case '\"':
            oss << "\\\"";
            break;
        case '\\':
            oss << "\\\\";
            break;
        case '\b':
            oss << "\\b";
            break;
        case '\f':
            oss << "\\f";
            break;
        case '\n':
            oss << "\\n";
            break;
        case '\r':
            oss << "\\r";
            break;
        case '\t':
            oss << "\\t";
            break;
        default:
            if ((unsigned char)c < 0x20)
            {
                oss << "\\u"
                    << std::hex << std::setw(4) << std::setfill('0')
                    << (int)(unsigned char)c
                    << std::dec << std::setfill(' ');
            }
            else
            {
                oss << c;
            }
        }
    }
    return oss.str();
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------

int main_stage1(
    const mpfr_t t,
    const std::string job_id,
    const std::string output_json)
{
    const mpfr_prec_t P = CFG().mpfr_prec;

    mpfr_t t0, h_mp, spacing_mp, t_left, t_right;
    mpfr_inits2(P, t0, h_mp, spacing_mp, t_left, t_right, (mpfr_ptr)0);

    mpfr_set(t0, t, MPFR_RNDN);

    const double spacing = est_spacing(t0);
    const double h = choose_h_stage1(t0);

    mpfr_set_d(h_mp, h, MPFR_RNDN);
    mpfr_set_d(spacing_mp, spacing, MPFR_RNDN);

    mpfr_sub(t_left, t0, h_mp, MPFR_RNDN);
    mpfr_add(t_right, t0, h_mp, MPFR_RNDN);

    const std::string t0_s = mpfr_to_string(t0, CFG().rs_log_digits);
    const std::string t_left_s = mpfr_to_string(t_left, CFG().rs_log_digits);
    const std::string t_right_s = mpfr_to_string(t_right, CFG().rs_log_digits);
    const std::string h_s = mpfr_to_string(h_mp, CFG().rs_log_digits);
    const std::string spacing_s = mpfr_to_string(spacing_mp, CFG().rs_log_digits);

    try
    {
        std::filesystem::path out_path(output_json);
        std::filesystem::path dir = out_path.parent_path();

        if (!dir.empty() && !std::filesystem::exists(dir))
        {
            LOGMP_DEBUG("Creating directory structure: %s", dir.string().c_str());
            std::filesystem::create_directories(dir);
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR creating directory structure for "
                  << output_json << " : " << e.what() << "\n";

        mpfr_clears(t0, h_mp, spacing_mp, t_left, t_right, (mpfr_ptr)0);
        return 4;
    }

    std::ofstream out(output_json, std::ios::out | std::ios::trunc);

    if (!out)
    {
        std::cerr << "ERROR: cannot open output file: " << output_json << "\n";

        mpfr_clears(t0, h_mp, spacing_mp, t_left, t_right, (mpfr_ptr)0);
        return 3;
    }

    out << "{\n";
    out << "  \"action\": \"queue_z_evaluations\",\n";
    out << "  \"stage_name\": \"refine_stage_1\",\n";
    out << "  \"reason\": \"initial_newton_sampling\",\n";
    out << "  \"refine_id\": \"" << json_escape(job_id) << "\",\n";
    out << "  \"input\": {\n";
    out << "    \"t0\": \"" << t0_s << "\"\n";
    out << "  },\n";
    out << "  \"config\": {\n";
    out << "    \"mpfr_prec\": " << (unsigned long)CFG().mpfr_prec << ",\n";
    out << "    \"rs_log_digits\": " << CFG().rs_log_digits << "\n";
    out << "  },\n";
    out << "  \"metadata\": {\n";
    out << "    \"spacing_est\": \"" << spacing_s << "\",\n";
    out << "    \"h\": \"" << h_s << "\",\n";
    out << "    \"sample_count\": 3\n";
    out << "  },\n";
    out << "  \"requested_points\": [\n";
    out << "    {\n";
    out << "      \"label\": \"left1\",\n";
    out << "      \"t\": \"" << t_left_s << "\"\n";
    out << "    },\n";
    out << "    {\n";
    out << "      \"label\": \"center\",\n";
    out << "      \"t\": \"" << t0_s << "\"\n";
    out << "    },\n";
    out << "    {\n";
    out << "      \"label\": \"right1\",\n";
    out << "      \"t\": \"" << t_right_s << "\"\n";
    out << "    }\n";
    out << "  ]\n";
    out << "}\n";

    out.close();

    LOGMP_DEBUG("refine_stage1_request wrote JSON: %s", output_json.c_str());
    LOGMP_DEBUG("refine_id=%s t0=%s h=%s spacing=%s",
                job_id.c_str(),
                t0_s.c_str(),
                h_s.c_str(),
                spacing_s.c_str());

    mpfr_clears(t0, h_mp, spacing_mp, t_left, t_right, (mpfr_ptr)0);

    return 0;
}