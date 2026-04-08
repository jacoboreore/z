// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// File: refine_stage2_decision.cpp
// Purpose:
//   Stage 2 decision engine for distributed refinement.
//   Reads stage2 JSON (3 completed evaluations) and emits the next decision:
//     - accept_root
//     - bracket_detected -> request midpoint evaluation
//     - newton_step -> request 3 new evaluations around next center
//
// Usage:
//   refine_stage2_decision <input_json> <output_json>
//
// Notes:
//   - Uses nlohmann/json for parsing/writing.
//   - Reuses existing mpfr_to_string(), est_spacing(), CFG() style.
//   - Hardcoded tolerances for now, to be parameterized later.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mpfr.h>
#include <nlohmann/json.hpp>

#include "confirm_zero.h"
#include "settings.h"
#include "logs.h"
#include "utils.h"

using json = nlohmann::json;

// -----------------------------------------------------------------------------
// Hardcoded for now; parameterize later
// -----------------------------------------------------------------------------
static constexpr double STAGE2_TOL_ABS = 1e-12;
static constexpr double DERIV_EPS = 1e-30;
static constexpr double MAX_STEP_FACTOR = 0.7;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

static bool load_json_file(const std::string &path, json &j)
{
    std::ifstream in(path);
    if (!in)
        return false;
    in >> j;
    return true;
}

static bool save_json_file(const std::string &path, const json &j)
{
    try
    {
        std::filesystem::path out_path(path);
        std::filesystem::path dir = out_path.parent_path();

        if (!dir.empty() && !std::filesystem::exists(dir))
            std::filesystem::create_directories(dir);

        std::ofstream out(path, std::ios::out | std::ios::trunc);
        if (!out)
            return false;

        out << j.dump(2) << "\n";
        return true;
    }
    catch (...)
    {
        return false;
    }
}

static bool mpfr_from_string(mpfr_t x, const std::string &s)
{
    return mpfr_set_str(x, s.c_str(), 10, MPFR_RNDN) == 0;
}

static double mpfr_abs_to_double(const mpfr_t x)
{
    mpfr_t a;
    mpfr_init2(a, mpfr_get_prec(x));
    mpfr_abs(a, x, MPFR_RNDN);
    double d = mpfr_get_d(a, MPFR_RNDN);
    mpfr_clear(a);
    return d;
}

static int mpfr_sign_of(const mpfr_t x)
{
    return mpfr_sgn(x);
}

static bool opposite_sign(const mpfr_t a, const mpfr_t b)
{
    const int sa = mpfr_sign_of(a);
    const int sb = mpfr_sign_of(b);
    return sa && sb && (sa != sb);
}

static double choose_h_stage(const mpfr_t t)
{
    const double s = est_spacing(t);
    const double eps = std::ldexp(1.0, -(int)CFG().mpfr_prec);
    const double mag = std::max(1.0, std::fabs(mpfr_get_d(t, MPFR_RNDN)));
    const double h0 = std::pow(eps, 1.0 / 3.0) * mag;
    return std::clamp(h0, 1e-9, 0.1 * s);
}

static json make_point(const std::string &label, const mpfr_t t)
{
    return json{
        {"label", label},
        {"t", mpfr_to_string(t, CFG().rs_log_digits)}};
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main_stage2(const std::string input_json, const std::string output_json)
{
    json in;
    if (!load_json_file(input_json, in))
    {
        LOGMPF("ERROR: could not load input json: %s", input_json);
        return 2;
    }

    try
    {
        const std::string refine_id = in.at("refine_id").get<std::string>();
        const std::string t0_str = in.at("input").at("t0").get<std::string>();

        const auto evals = in.at("evaluations");
        if (!evals.is_array() || evals.size() != 3)
        {
            LOGMPF("ERROR: expected exactly 3 evaluations");
            return 3;
        }

        std::string left_t_str, center_t_str, right_t_str;
        std::string left_z_str, center_z_str, right_z_str;
        std::string left_status, center_status, right_status;

        for (const auto &e : evals)
        {
            const std::string label = e.at("label").get<std::string>();
            const std::string t_str = e.at("t").get<std::string>();
            const std::string z_str = e.at("z_value").get<std::string>();
            const std::string st = e.at("status").get<std::string>();

            if (label == "left1")
            {
                left_t_str = t_str;
                left_z_str = z_str;
                left_status = st;
            }
            else if (label == "center")
            {
                center_t_str = t_str;
                center_z_str = z_str;
                center_status = st;
            }
            else if (label == "right1")
            {
                right_t_str = t_str;
                right_z_str = z_str;
                right_status = st;
            }
        }

        if (left_status != "COMPLETED" || center_status != "COMPLETED" || right_status != "COMPLETED")
        {
            json out = {
                {"action", "wait"},
                {"stage_name", "refine_stage_2"},
                {"reason", "not_all_evaluations_completed"},
                {"refine_id", refine_id}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                return 4;
            }
            return 0;
        }

        const mpfr_prec_t P = CFG().mpfr_prec;

        mpfr_t t_left, t_center, t_right;
        mpfr_t z_left, z_center, z_right;
        mpfr_t t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid;
        mpfr_inits2(P,
                    t_left, t_center, t_right,
                    z_left, z_center, z_right,
                    t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                    (mpfr_ptr)0);

        if (!mpfr_from_string(t_left, left_t_str) ||
            !mpfr_from_string(t_center, center_t_str) ||
            !mpfr_from_string(t_right, right_t_str) ||
            !mpfr_from_string(z_left, left_z_str) ||
            !mpfr_from_string(z_center, center_z_str) ||
            !mpfr_from_string(z_right, right_z_str))
        {
            LOGMPF("ERROR: invalid numeric field in input json");
            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 5;
        }

        const double abs_center = mpfr_abs_to_double(z_center);

        // ---------------------------------------------------------------------
        // Decision 1: already good enough
        // ---------------------------------------------------------------------
        if (abs_center <= STAGE2_TOL_ABS)
        {
            json out = {
                {"action", "accept_root"},
                {"stage_name", "refine_stage_2"},
                {"reason", "center_value_within_tolerance"},
                {"refine_id", refine_id},
                {"result", {{"root_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"root_abs_z", mpfr_to_string(z_center, CFG().rs_log_digits)}, {"tolerance_used", STAGE2_TOL_ABS}}}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                mpfr_clears(t_left, t_center, t_right,
                            z_left, z_center, z_right,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                            (mpfr_ptr)0);
                return 6;
            }

            LOGMPF("Stage2 decision: accept_root refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 0;
        }

        // ---------------------------------------------------------------------
        // Decision 2: bracket detected
        // ---------------------------------------------------------------------
        if (opposite_sign(z_left, z_center))
        {
            mpfr_add(t_mid, t_left, t_center, MPFR_RNDN);
            mpfr_div_ui(t_mid, t_mid, 2, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "bracket_detected_left_center"},
                {"refine_id", refine_id},
                {"metadata", {{"bracket_left_t", mpfr_to_string(t_left, CFG().rs_log_digits)}, {"bracket_right_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"sample_count", 1}}},
                {"requested_points", json::array({make_point("mid_bracket", t_mid)})}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                mpfr_clears(t_left, t_center, t_right,
                            z_left, z_center, z_right,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                            (mpfr_ptr)0);
                return 7;
            }

            LOGMPF("Stage2 decision: bracket left-center refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 0;
        }

        if (opposite_sign(z_center, z_right))
        {
            mpfr_add(t_mid, t_center, t_right, MPFR_RNDN);
            mpfr_div_ui(t_mid, t_mid, 2, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "bracket_detected_center_right"},
                {"refine_id", refine_id},
                {"metadata", {{"bracket_left_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"bracket_right_t", mpfr_to_string(t_right, CFG().rs_log_digits)}, {"sample_count", 1}}},
                {"requested_points", json::array({make_point("mid_bracket", t_mid)})}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                mpfr_clears(t_left, t_center, t_right,
                            z_left, z_center, z_right,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                            (mpfr_ptr)0);
                return 8;
            }

            LOGMPF("Stage2 decision: bracket center-right refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 0;
        }

        // ---------------------------------------------------------------------
        // Decision 3: central-difference Newton step
        // dZ ≈ (Zr - Zl) / (tr - tl)
        // t_next = t_center - Z_center / dZ
        // ---------------------------------------------------------------------
        mpfr_sub(ztmp, z_right, z_left, MPFR_RNDN); // numerator
        mpfr_sub(dZ, t_right, t_left, MPFR_RNDN);   // temp denominator
        if (mpfr_zero_p(dZ))
        {
            LOGMPF("ERROR: t_right - t_left == 0");
            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 9;
        }

        mpfr_div(dZ, ztmp, dZ, MPFR_RNDN);

        const double abs_dz = mpfr_abs_to_double(dZ);
        if (abs_dz < DERIV_EPS)
        {
            // Derivative too small -> ask for wider scan around center
            const double h = choose_h_stage(t_center) * 2.0;
            mpfr_set_d(h_mp, h, MPFR_RNDN);

            mpfr_sub(t_req_left, t_center, h_mp, MPFR_RNDN);
            mpfr_add(t_req_right, t_center, h_mp, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "derivative_too_small_request_wider_scan"},
                {"refine_id", refine_id},
                {"metadata", {{"base_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)}, {"sample_count", 3}}},
                {"requested_points", json::array({make_point("left_wide", t_req_left),
                                                  make_point("center", t_center),
                                                  make_point("right_wide", t_req_right)})}};

            if (!save_json_file(output_json, out))
            {
                LOGMPF("ERROR: failed to write output json: %s", output_json);
                mpfr_clears(t_left, t_center, t_right,
                            z_left, z_center, z_right,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                            (mpfr_ptr)0);
                return 10;
            }

            LOGMPF("Stage2 decision: derivative too small, wider scan refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 0;
        }

        // delta = Zc / dZ
        mpfr_div(delta, z_center, dZ, MPFR_RNDN);

        // clamp delta to a fraction of spacing
        const double spacing = est_spacing(t_center);
        const double step_max = MAX_STEP_FACTOR * spacing;
        double delta_abs = std::fabs(mpfr_get_d(delta, MPFR_RNDN));

        if (delta_abs > step_max)
        {
            const double scale = step_max / delta_abs;
            mpfr_mul_d(delta, delta, scale, MPFR_RNDN);
        }

        mpfr_sub(t_next, t_center, delta, MPFR_RNDN);

        const double h_next = choose_h_stage(t_next);
        mpfr_set_d(h_mp, h_next, MPFR_RNDN);

        mpfr_sub(t_req_left, t_next, h_mp, MPFR_RNDN);
        mpfr_add(t_req_right, t_next, h_mp, MPFR_RNDN);

        json out = {
            {"action", "queue_z_evaluations"},
            {"stage_name", "refine_stage_3"},
            {"reason", "newton_step_from_stage2"},
            {"refine_id", refine_id},
            {"input", {{"t0", t0_str}}},
            {"metadata", {{"previous_center_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"previous_center_z", mpfr_to_string(z_center, CFG().rs_log_digits)}, {"derivative_est", mpfr_to_string(dZ, CFG().rs_log_digits)}, {"delta", mpfr_to_string(delta, CFG().rs_log_digits)}, {"spacing_est", mpfr_to_string(t_center, CFG().rs_log_digits)}, // optional placeholder
                          {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)},
                          {"sample_count", 3}}},
            {"requested_points", json::array({make_point("left1", t_req_left),
                                              make_point("center", t_next),
                                              make_point("right1", t_req_right)})}};

        if (!save_json_file(output_json, out))
        {
            std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
            mpfr_clears(t_left, t_center, t_right,
                        z_left, z_center, z_right,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                        (mpfr_ptr)0);
            return 11;
        }

        LOGMPF("Stage2 decision: newton step refine_id=%s", refine_id.c_str());
        LOGMPF("refine_id=%s t_next=%s delta=%s dZ=%s",
               refine_id.c_str(),
               mpfr_to_string(t_next, CFG().rs_log_digits).c_str(),
               mpfr_to_string(delta, CFG().rs_log_digits).c_str(),
               mpfr_to_string(dZ, CFG().rs_log_digits).c_str());

        mpfr_clears(t_left, t_center, t_right,
                    z_left, z_center, z_right,
                    t_next, ztmp, dZ, delta, h_mp, t_req_left, t_req_right, t_mid,
                    (mpfr_ptr)0);
        return 0;
    }
    catch (const std::exception &e)
    {
        LOGMPF("ERROR: exception in stage2 decision: ", e.what());
        return 12;
    }
}