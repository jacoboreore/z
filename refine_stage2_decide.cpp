// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// File: refine_stage2_decision.cpp
// Purpose:
//   Stage 2 decision engine for distributed refinement.
//   Reads stage2 JSON (3 completed evaluations) and emits the next decision:
//     - accept_root
//     - bracket_detected -> request 3-point
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
static constexpr double STAGE2_TOL_ABS = 1e-8;
static constexpr double DERIV_EPS = 1e-30;
static constexpr double MAX_STEP_FACTOR = 0.7;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

struct EvalPoint
{
    std::string label;
    std::string t_str;
    std::string z_str;
    std::string status;
};

static bool have_label(const std::vector<EvalPoint> &pts, const std::string &label)
{
    for (const auto &p : pts)
        if (p.label == label)
            return true;
    return false;
}

static const EvalPoint *find_label(const std::vector<EvalPoint> &pts, const std::string &label)
{
    for (const auto &p : pts)
        if (p.label == label)
            return &p;
    return nullptr;
}

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
    catch (const std::exception &e)
    {
        LOGMPF("ERROR: exception in save_json_file: %s", e.what());
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
        LOGMPF("ERROR: could not load input json: %s", input_json.c_str());
        return 2;
    }

    try
    {
        const std::string refine_id = in.at("refine_id").get<std::string>();
        const std::string t0_str = in.at("input").at("t0").get<std::string>();

        const auto evals = in.at("evaluations");
        if (!evals.is_array() || evals.empty())
        {
            LOGMPF("ERROR: expected non-empty evaluations array");
            return 3;
        }

        std::vector<EvalPoint> pts;
        for (const auto &e : evals)
        {
            pts.push_back(EvalPoint{
                e.at("label").get<std::string>(),
                e.at("t").get<std::string>(),
                e.at("z_value").get<std::string>(),
                e.at("status").get<std::string>()});
        }

        const bool has3 =
            have_label(pts, "left1") &&
            have_label(pts, "center") &&
            have_label(pts, "right1");

        const bool has5 =
            have_label(pts, "left2") &&
            have_label(pts, "left1") &&
            have_label(pts, "center") &&
            have_label(pts, "right1") &&
            have_label(pts, "right2");

        LOGMPF("Stage2 input refine_id=%s eval_count=%d has3=%d has5=%d",
               refine_id.c_str(),
               (int)pts.size(),
               (int)has3,
               (int)has5);

        for (const auto &p : pts)
        {
            LOGMPF("Stage2 eval refine_id=%s label=%s t=%s status=%s",
                   refine_id.c_str(),
                   p.label.c_str(),
                   p.t_str.c_str(),
                   p.status.c_str());
        }

        if (!has3)
        {
            LOGMPF("ERROR: missing required minimum 3-point labels");
            return 3;
        }

        const EvalPoint *p_left1 = find_label(pts, "left1");
        const EvalPoint *p_center = find_label(pts, "center");
        const EvalPoint *p_right1 = find_label(pts, "right1");
        const EvalPoint *p_left2 = find_label(pts, "left2");
        const EvalPoint *p_right2 = find_label(pts, "right2");

        for (const auto &p : pts)
        {
            if (p.status != "COMPLETED")
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

                LOGMPF("Stage2 decision: wait refine_id=%s because not all evals are completed", refine_id.c_str());
                return 0;
            }
        }

        const mpfr_prec_t P = CFG().mpfr_prec;

        mpfr_t t_left2, t_left1, t_center, t_right1, t_right2;
        mpfr_t z_left2, z_left1, z_center, z_right1, z_right2;
        mpfr_t t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom;
        mpfr_inits2(P,
                    t_left2, t_left1, t_center, t_right1, t_right2,
                    z_left2, z_left1, z_center, z_right1, z_right2,
                    t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                    (mpfr_ptr)0);

        if (!mpfr_from_string(t_left1, p_left1->t_str) ||
            !mpfr_from_string(t_center, p_center->t_str) ||
            !mpfr_from_string(t_right1, p_right1->t_str) ||
            !mpfr_from_string(z_left1, p_left1->z_str) ||
            !mpfr_from_string(z_center, p_center->z_str) ||
            !mpfr_from_string(z_right1, p_right1->z_str))
        {
            LOGMPF("ERROR: invalid numeric field in mandatory 3-point input json");
            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 5;
        }

        if (has5)
        {
            if (!mpfr_from_string(t_left2, p_left2->t_str) ||
                !mpfr_from_string(t_right2, p_right2->t_str) ||
                !mpfr_from_string(z_left2, p_left2->z_str) ||
                !mpfr_from_string(z_right2, p_right2->z_str))
            {
                LOGMPF("ERROR: invalid numeric field in 5-point input json");
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 5;
            }
        }

        const double abs_center = mpfr_abs_to_double(z_center);

        LOGMPF("Stage2 refine_id=%s center_abs_z=%.17g", refine_id.c_str(), abs_center);

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
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 6;
            }

            LOGMPF("Stage2 decision: accept_root refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 0;
        }

        // ---------------------------------------------------------------------
        // Decision 2: bracket detected
        // Use immediate neighbors around center
        // ---------------------------------------------------------------------
        if (opposite_sign(z_left1, z_center))
        {
            mpfr_add(t_mid, t_left1, t_center, MPFR_RNDN);
            mpfr_div_ui(t_mid, t_mid, 2, MPFR_RNDN);

            const double h_next = choose_h_stage(t_mid);
            mpfr_set_d(h_mp, h_next, MPFR_RNDN);

            mpfr_sub(t_req_left1, t_mid, h_mp, MPFR_RNDN);
            mpfr_add(t_req_right1, t_mid, h_mp, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "bracket_detected_left_center"},
                {"refine_id", refine_id},
                {"input", {{"t0", t0_str}}},
                {"metadata", {{"bracket_left_t", mpfr_to_string(t_left1, CFG().rs_log_digits)}, {"bracket_right_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"mid_t", mpfr_to_string(t_mid, CFG().rs_log_digits)}, {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)}, {"sample_count", 3}}},
                {"requested_points", json::array({make_point("left1", t_req_left1),
                                                  make_point("center", t_mid),
                                                  make_point("right1", t_req_right1)})}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 7;
            }

            LOGMPF("Stage2 decision: bracket left-center refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 0;
        }

        if (opposite_sign(z_center, z_right1))
        {
            mpfr_add(t_mid, t_center, t_right1, MPFR_RNDN);
            mpfr_div_ui(t_mid, t_mid, 2, MPFR_RNDN);

            const double h_next = choose_h_stage(t_mid);
            mpfr_set_d(h_mp, h_next, MPFR_RNDN);

            mpfr_sub(t_req_left1, t_mid, h_mp, MPFR_RNDN);
            mpfr_add(t_req_right1, t_mid, h_mp, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "bracket_detected_center_right"},
                {"refine_id", refine_id},
                {"input", {{"t0", t0_str}}},
                {"metadata", {{"bracket_left_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"bracket_right_t", mpfr_to_string(t_right1, CFG().rs_log_digits)}, {"mid_t", mpfr_to_string(t_mid, CFG().rs_log_digits)}, {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)}, {"sample_count", 3}}},
                {"requested_points", json::array({make_point("left1", t_req_left1),
                                                  make_point("center", t_mid),
                                                  make_point("right1", t_req_right1)})}};

            if (!save_json_file(output_json, out))
            {
                std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 8;
            }

            LOGMPF("Stage2 decision: bracket center-right refine_id=%s", refine_id.c_str());

            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 0;
        }

        // ---------------------------------------------------------------------
        // Decision 3: derivative estimate
        // 3-point or 5-point depending on available labels
        // ---------------------------------------------------------------------
        bool used_5pt = false;

        if (has5)
        {
            // dZ = [ -Z(+2h) + 8Z(+h) - 8Z(-h) + Z(-2h) ] / (12 h)
            mpfr_sub(h_mp, t_right1, t_center, MPFR_RNDN);
            if (mpfr_zero_p(h_mp))
            {
                LOGMPF("ERROR: right1 - center == 0 in 5-point derivative");
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 9;
            }

            mpfr_neg(dZ, z_right2, MPFR_RNDN);         // -Z(+2h)
            mpfr_mul_ui(ztmp, z_right1, 8, MPFR_RNDN); // +8Z(+h)
            mpfr_add(dZ, dZ, ztmp, MPFR_RNDN);
            mpfr_mul_ui(ztmp, z_left1, 8, MPFR_RNDN); // -8Z(-h)
            mpfr_sub(dZ, dZ, ztmp, MPFR_RNDN);
            mpfr_add(dZ, dZ, z_left2, MPFR_RNDN); // +Z(-2h)

            mpfr_mul_ui(denom, h_mp, 12, MPFR_RNDN);
            mpfr_div(dZ, dZ, denom, MPFR_RNDN);

            used_5pt = true;
            LOGMPF("Stage2 derivative mode refine_id=%s mode=5-point", refine_id.c_str());
        }
        else
        {
            mpfr_sub(ztmp, z_right1, z_left1, MPFR_RNDN);
            mpfr_sub(denom, t_right1, t_left1, MPFR_RNDN);
            if (mpfr_zero_p(denom))
            {
                LOGMPF("ERROR: t_right1 - t_left1 == 0 in 3-point derivative");
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 9;
            }

            mpfr_div(dZ, ztmp, denom, MPFR_RNDN);
            LOGMPF("Stage2 derivative mode refine_id=%s mode=3-point", refine_id.c_str());
        }

        const double abs_dz = mpfr_abs_to_double(dZ);
        LOGMPF("Stage2 derivative refine_id=%s abs_dz=%.17g", refine_id.c_str(), abs_dz);

        if (abs_dz < DERIV_EPS)
        {
            // Ask for 5-point scan around center
            const double h = choose_h_stage(t_center);
            mpfr_set_d(h_mp, h, MPFR_RNDN);
            mpfr_mul_ui(denom, h_mp, 2, MPFR_RNDN);

            mpfr_sub(t_req_left2, t_center, denom, MPFR_RNDN);
            mpfr_sub(t_req_left1, t_center, h_mp, MPFR_RNDN);
            mpfr_add(t_req_right1, t_center, h_mp, MPFR_RNDN);
            mpfr_add(t_req_right2, t_center, denom, MPFR_RNDN);

            json out = {
                {"action", "queue_z_evaluations"},
                {"stage_name", "refine_stage_3"},
                {"reason", "derivative_too_small_request_5_point_scan"},
                {"refine_id", refine_id},
                {"input", {{"t0", t0_str}}},
                {"metadata", {{"base_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)}, {"derivative_mode", used_5pt ? "5-point" : "3-point"}, {"sample_count", 5}}},
                {"requested_points", json::array({make_point("left2", t_req_left2),
                                                  make_point("left1", t_req_left1),
                                                  make_point("center", t_center),
                                                  make_point("right1", t_req_right1),
                                                  make_point("right2", t_req_right2)})}};

            if (!save_json_file(output_json, out))
            {
                LOGMPF("ERROR: failed to write output json: %s", output_json.c_str());
                mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                            z_left2, z_left1, z_center, z_right1, z_right2,
                            t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                            (mpfr_ptr)0);
                return 10;
            }

            LOGMPF("Stage2 decision: derivative too small refine_id=%s -> requesting 5 points", refine_id.c_str());

            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 0;
        }

        // delta = Zc / dZ
        mpfr_div(delta, z_center, dZ, MPFR_RNDN);

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
        mpfr_mul_ui(denom, h_mp, 2, MPFR_RNDN);

        mpfr_sub(t_req_left1, t_next, h_mp, MPFR_RNDN);
        mpfr_add(t_req_right1, t_next, h_mp, MPFR_RNDN);

        json requested_points;
        int next_sample_count = 3;

        if (used_5pt)
        {
            mpfr_sub(t_req_left2, t_next, denom, MPFR_RNDN);
            mpfr_add(t_req_right2, t_next, denom, MPFR_RNDN);

            requested_points = json::array({make_point("left2", t_req_left2),
                                            make_point("left1", t_req_left1),
                                            make_point("center", t_next),
                                            make_point("right1", t_req_right1),
                                            make_point("right2", t_req_right2)});
            next_sample_count = 5;
        }
        else
        {
            requested_points = json::array({make_point("left1", t_req_left1),
                                            make_point("center", t_next),
                                            make_point("right1", t_req_right1)});
            next_sample_count = 3;
        }

        json out = {
            {"action", "queue_z_evaluations"},
            {"stage_name", "refine_stage_3"},
            {"reason", "newton_step_from_stage2"},
            {"refine_id", refine_id},
            {"input", {{"t0", t0_str}}},
            {"metadata", {{"previous_center_t", mpfr_to_string(t_center, CFG().rs_log_digits)}, {"previous_center_z", mpfr_to_string(z_center, CFG().rs_log_digits)}, {"derivative_est", mpfr_to_string(dZ, CFG().rs_log_digits)}, {"delta", mpfr_to_string(delta, CFG().rs_log_digits)}, {"spacing_est", std::to_string(spacing)}, {"h", mpfr_to_string(h_mp, CFG().rs_log_digits)}, {"derivative_mode", used_5pt ? "5-point" : "3-point"}, {"sample_count", next_sample_count}}},
            {"requested_points", requested_points}};

        if (!save_json_file(output_json, out))
        {
            std::cerr << "ERROR: failed to write output json: " << output_json << "\n";
            mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                        z_left2, z_left1, z_center, z_right1, z_right2,
                        t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                        (mpfr_ptr)0);
            return 11;
        }

        LOGMPF("Stage2 decision: newton step refine_id=%s next_sample_count=%d derivative_mode=%s",
               refine_id.c_str(),
               next_sample_count,
               used_5pt ? "5-point" : "3-point");

        LOGMPF("refine_id=%s t_next=%s delta=%s dZ=%s",
               refine_id.c_str(),
               mpfr_to_string(t_next, CFG().rs_log_digits).c_str(),
               mpfr_to_string(delta, CFG().rs_log_digits).c_str(),
               mpfr_to_string(dZ, CFG().rs_log_digits).c_str());

        mpfr_clears(t_left2, t_left1, t_center, t_right1, t_right2,
                    z_left2, z_left1, z_center, z_right1, z_right2,
                    t_next, ztmp, dZ, delta, h_mp, t_req_left2, t_req_left1, t_req_right1, t_req_right2, t_mid, denom,
                    (mpfr_ptr)0);
        return 0;
    }
    catch (const std::exception &e)
    {
        LOGMPF("ERROR: exception in stage2 decision: %s", e.what());
        return 12;
    }
}