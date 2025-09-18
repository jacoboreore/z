#include <mpfr.h>
#ifdef HAVE_MPC
#include <mpc>
#endif
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <mutex>
#include <type_traits>
#include <string>
#include <sstream>
#include <iomanip>
#include "settings.h"
#include "logs.h"

#ifdef OS_USE_MPFR
#include "mpreal.h"
#define HAVE_MPREAL 1
#else
#define HAVE_MPREAL 0
#endif

template <class T>
static inline std::enable_if_t<std::is_integral_v<T>, std::string>
fmt_num(const T &v, int /*digits*/)
{
    return std::to_string(v);
}

static inline std::string fmt_num(const bool &v, int /*digits*/)
{
    return v ? "true" : "false";
}

template <class T>
static inline std::enable_if_t<std::is_floating_point_v<T>, std::string>
fmt_num(const T &v, int digits)
{
    std::ostringstream oss;
    // oss.setf(std::ios::scientific);
    oss << v;
    return oss.str();
}

#if HAVE_MPREAL
static inline std::string fmt_num(const mpfr::mpreal &v, int digits)
{
    // mpfr::mpreal → string with desired digits (base10)
    return v.toString(10, digits);
}
#endif

// Fallback for things like std::string / const char*
static inline std::string fmt_num(const std::string &s, int) { return s; }
static inline std::string fmt_num(const char *s, int) { return s ? std::string(s) : std::string("(null)"); }

template <class T>
static inline void log_kv(const char *key, const T &value, int digits = 18)
{
    const std::string sval = fmt_num(value, digits);
    LOGMP_DEBUG("%-26s = %s", key, sval.c_str());
}

static Config g_cfg;          // the single global instance
static std::once_flag g_once; // guards one-time init
static std::mutex g_mutex;    // guards reload

// ---------- helpers ----------
static inline std::string trim(std::string s)
{
    auto issp = [](unsigned char c)
    { return std::isspace(c); };

    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c)
                                    { return !issp(c); }));

    s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c)
                         { return !issp(c); })
                .base(),
            s.end());

    return s;
}

static inline bool ieq(const std::string &a, const char *b)
{
    size_t blen = std::char_traits<char>::length(b);
    if (a.size() != blen)
        return false;
    for (size_t i = 0; i < blen; ++i)
        if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    return true;
}

static inline bool to_bool(std::string s)
{
    s = trim(s);
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c)
                   { return std::tolower(c); });
    return (s == "1" || s == "true" || s == "yes" || s == "on");
}

static inline void set_kv(Config &c, const std::string &k, const std::string &v)
{
    // All application precision term
    if (ieq(k, "mpfr_prec"))
        c.mpfr_prec = std::max(2, std::stoi(v));

    else if (ieq(k, "prefilter_tau_boost"))
        c.prefilter_tau_boost = std::stod(v);

    else if (ieq(k, "prefilter_alpha"))
        c.prefilter_alpha = std::stod(v);
    else if (ieq(k, "prefilter_tau"))
        c.prefilter_tau = std::stod(v);

    else if (ieq(k, "newton_tol_abs"))
        c.newton_tol_abs = std::stod(v);
    else if (ieq(k, "newton_tol_step"))
        c.newton_tol_step = std::stod(v);
    else if (ieq(k, "newton_max_iter"))
        c.newton_max_iter = std::max(1, std::stoi(v));

    else if (ieq(k, "use_threaded_z"))
        c.use_threaded_z = to_bool(v);
    else if (ieq(k, "use_threaded_newton"))
        c.use_threaded_newton = to_bool(v);

    else if (ieq(k, "ball_mid_threshold"))
        c.ball_mid_threshold = std::stod(v);
    else if (ieq(k, "h_min"))
        c.h_min = std::stod(v);
    else if (ieq(k, "h_frac_spacing"))
        c.h_frac_spacing = std::stod(v);
    else if (ieq(k, "step_frac_spacing"))
        c.step_frac_spacing = std::stod(v);

    else if (ieq(k, "min_threads"))
        c.min_threads = std::stod(v);
    else if (ieq(k, "max_threads"))
        c.max_threads = std::stod(v);
    else if (ieq(k, "rs_min_N_for_threads"))
        c.rs_min_N_for_threads = std::stod(v);
    else if (ieq(k, "rs_tail_C"))
        c.rs_tail_C = std::stod(v);

    else if (ieq(k, "rs_log_debug"))
        c.rs_log_debug = to_bool(v);
    else if (ieq(k, "rs_main_log"))
        c.rs_main_log = to_bool(v);

    else if (ieq(k, "rs_log_digits"))
        c.rs_log_digits = std::stod(v);

    else if (ieq(k, "rs_half_weight_last"))
        c.rs_half_weight_last = to_bool(v);
    else if (ieq(k, "theta_exact_switch_t"))
        c.theta_exact_switch_t = std::stod(v);
    else if (ieq(k, "theta_stirling_terms"))
        c.theta_stirling_terms = std::stod(v);

    else if (ieq(k, "valley_grow"))
        c.valley_grow = std::stod(v);
    else if (ieq(k, "valley_aggressive_factor"))
        c.valley_aggressive_factor = std::stod(v);
    else if (ieq(k, "valley_tolerance_to_step"))
        c.valley_tolerance_to_step = std::stod(v);
    else if (ieq(k, "valley_max_step"))
        c.valley_max_step = std::stod(v);
    else if (ieq(k, "valley_max_Zs_relation"))
        c.valley_max_Zs_relation = std::stod(v);

    else if (ieq(k, "enable_min_t"))
        c.enable_min_t = std::stold(v);
    else if (ieq(k, "samples_per_spacing"))
        c.samples_per_spacing = std::stoi(v);
    else if (ieq(k, "K_points"))
        c.K_points = std::stoi(v);
    else if (ieq(k, "delta_log_max"))
        c.delta_log_max = std::stold(v);
    else if (ieq(k, "poly_deg_R"))
        c.poly_deg_R = std::stoi(v);
    else if (ieq(k, "return_main_sum"))
        c.return_main_sum = to_bool(v);
    else if (ieq(k, "tol_refine"))
        c.tol_refine = std::stold(v);
    else if (ieq(k, "z_far"))
        c.z_far = std::stold(v);
    else if (ieq(k, "ratio_far"))
        c.ratio_far = std::stold(v);
    else if (ieq(k, "min_mult"))
        c.min_mult = std::stold(v);
    else if (ieq(k, "eps_abs"))
        c.eps_abs = std::stold(v);
    else if (ieq(k, "eps_rel"))
        c.eps_rel = std::stold(v);
    else if (ieq(k, "flip_center_frac"))
        c.flip_center_frac = std::stold(v);
}

static void load_file(Config &cfg, const std::string &path)
{

    std::ifstream f;

    std::cout << "\nConfig file: " << path << " \n";

    f.open(path);
    if (!f.is_open())
    {
        std::cerr << "Error while reading config file: " << path;
        return;
    }

    std::string line;
    while (std::getline(f, line))
    {
        auto hash = line.find_first_of("#");
        if (hash != std::string::npos)
            line.erase(hash);
        auto eq = line.find('=');
        if (eq == std::string::npos)
            continue;
        auto k = trim(line.substr(0, eq));
        auto v = trim(line.substr(eq + 1));

        if (!k.empty() && !v.empty())
            set_kv(cfg, k, v);
    }
}

static void load_env(Config &cfg)
{
    auto get = [](const char *name) -> const char *
    { return std::getenv(name); };
#define OV(name, key)            \
    do                           \
    {                            \
        if (auto *s = get(name)) \
            set_kv(cfg, key, s); \
    } while (0)
    OV("ZETA_MPFR_PREC", "mpfr_prec");
    OV("ZETA_PREFILTER_ALPHA", "prefilter_alpha");
    OV("ZETA_PREFILTER_TAU", "prefilter_tau");
    OV("ZETA_NEWTON_TOL_ABS", "newton_tol_abs");
    OV("ZETA_NEWTON_TOL_STEP", "newton_tol_step");
    OV("ZETA_NEWTON_MAX_ITER", "newton_max_iter");
    OV("ZETA_USE_THREADED_Z", "use_threaded_z");
    OV("ZETA_BALL_MID_THRESHOLD", "ball_mid_threshold");
    OV("ZETA_H_MIN", "h_min");
    OV("ZETA_H_FRAC_SPACING", "h_frac_spacing");
    OV("ZETA_STEP_FRAC_SPACING", "step_frac_spacing");

    OV("ZETA_MIN_THREADS", "min_threads");
    OV("ZETA_MAX_THREADS", "max_threads");
    OV("ZETA_RS_MIN_N_FOR_THREADS", "rs_min_N_for_threads");
    OV("ZETA_RS_TAIL_C", "rs_tail_C");
    OV("ZETA_RS_LOG_DEBUG", "rs_log_debug");
    OV("ZETA_RS_LOG_DIGITS", "rs_log_digits");

#undef OV
}

static void load_from_paths(Config &cfg, const char *path)
{
    load_file(cfg, path);
}

// ---------- public API ----------
void print_config(const Config &c)
{
    const int D = c.rs_log_digits > 0 ? c.rs_log_digits : 18;

    LOGMP_DEBUG("\n");

    log_kv("mpfr_prec", c.mpfr_prec);

    // --- prefilter ---
    log_kv("prefilter_alpha", c.prefilter_alpha, D);
    log_kv("prefilter_tau", c.prefilter_tau, D);
    log_kv("prefilter_tau_boost", c.prefilter_tau_boost, D);
    LOGMP_DEBUG("\n");

    // --- newton ---
    log_kv("newton_tol_abs", c.newton_tol_abs, D);
    log_kv("newton_tol_step", c.newton_tol_step, D);
    log_kv("newton_max_iter", c.newton_max_iter);
    LOGMP_DEBUG("\n");

    // --- threading / zero confirm ---
    log_kv("use_threaded_z", c.use_threaded_z);
    log_kv("use_threaded_newton", c.use_threaded_newton);
    log_kv("ball_mid_threshold", c.ball_mid_threshold, D);
    LOGMP_DEBUG("\n");

    // --- step/spacing ---
    log_kv("h_min", c.h_min, D);
    log_kv("h_frac_spacing", c.h_frac_spacing, D);
    log_kv("step_frac_spacing", c.step_frac_spacing, D);
    LOGMP_DEBUG("\n");

    // --- RS core ---
    log_kv("rs_half_weight_last", c.rs_half_weight_last);
    log_kv("rs_tail_C", c.rs_tail_C, D);
    log_kv("min_threads", c.min_threads);
    log_kv("max_threads", c.max_threads);
    log_kv("rs_min_N_for_threads", c.rs_min_N_for_threads);
    log_kv("theta_exact_switch_t", c.theta_exact_switch_t, D);
    LOGMP_DEBUG("\n");

    // --- logging ---
    log_kv("rs_log_debug", c.rs_log_debug);
    log_kv("rs_main_log", c.rs_main_log);
    log_kv("rs_log_digits", c.rs_log_digits);
    LOGMP_DEBUG("\n");

    // --- theta ---
    log_kv("theta_stirling_terms", c.theta_stirling_terms);
    LOGMP_DEBUG("\n");

    // --- valley ---
    log_kv("valley_grow", c.valley_grow, D);
    log_kv("valley_aggressive_factor", c.valley_aggressive_factor, D);
    log_kv("valley_max_step", c.valley_max_step, D);
    log_kv("valley_max_Zs_relation", c.valley_max_Zs_relation, D);
    LOGMP_DEBUG("\n");

    // --- OS / Odlyzko–Schönhage ---
    log_kv("enable_min_t", c.enable_min_t, D);
    log_kv("samples_per_spacing", c.samples_per_spacing);
    log_kv("K_points", c.K_points);
    log_kv("delta_log_max", c.delta_log_max, D);
    log_kv("poly_deg_R", c.poly_deg_R);
    log_kv("return_main_sum", c.return_main_sum);
    log_kv("tol_refine", c.tol_refine, D);
    log_kv("z_far", c.z_far, D);
    log_kv("ratio_far", c.ratio_far, D);
    log_kv("min_mult", c.min_mult, D);
    log_kv("eps_abs", c.eps_abs, D);
    log_kv("eps_rel", c.eps_rel, D);
    log_kv("flip_center_frac", c.flip_center_frac, D);
    LOGMP_DEBUG("\n");
}

void apply_precision(const Config &cfg)
{
    mpfr_set_default_prec(cfg.mpfr_prec);
}

void config_init(const char *path)
{
    std::call_once(g_once, [path]
                   {
        // start from built-in defaults in g_cfg
        load_from_paths(g_cfg, path);
        // TODO: get vars from environment
        //load_env(g_cfg);
        apply_precision(g_cfg); });
}

const Config &CFG()
{
    // lazy-init exactly once; no logging, no re-init
    config_init(nullptr);
    return g_cfg;
}

void config_reload(const char *path)
{
    std::lock_guard<std::mutex> lk(g_mutex);
    Config fresh; // start from defaults
    load_from_paths(fresh, path);
    load_env(fresh);
    g_cfg = fresh; // atomically swap contents
    apply_precision(g_cfg);
}

void api_version()
{
    LOGMP_DEBUG("------------------------------------------------------------------");

    // MPFR
    LOGMP_DEBUG("MPFR runtime : %s", mpfr_get_version());
#ifdef MPFR_VERSION_STRING
    LOGMP_DEBUG("MPFR headers : %s", MPFR_VERSION_STRING);
#else
    LOGMP_DEBUG("MPFR headers : (unknown)");
#endif

    // GMP
    LOGMP_DEBUG("GMP  runtime : %s", gmp_version);
#if defined(__GNU_MP_VERSION) && defined(__GNU_MP_VERSION_MINOR)
#if defined(__GNU_MP_VERSION_PATCHLEVEL)
    LOGMP_DEBUG("GMP  headers : %d.%d.%d",
                __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL);
#else
    LOGMP_DEBUG("GMP  headers : %d.%d",
                __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR);
#endif
#else
    LOGMP_DEBUG("GMP  headers : (unknown)");
#endif

    // MPC (optional)
#if defined(HAVE_MPC) || defined(MPC_VERSION_STRING) || defined(__MPC_H__)
#if defined(MPC_VERSION_STRING)
    LOGMP_DEBUG("MPC  headers : %s", MPC_VERSION_STRING);
#else
    LOGMP_DEBUG("MPC  headers : (unknown)");
#endif
#if defined(MPC_VERSION_MAJOR) || defined(mpc_get_version)
    LOGMP_DEBUG("MPC  runtime : %s", mpc_get_version());
#else
    LOGMP_DEBUG("MPC  runtime : (mpc_get_version unavailable)");
#endif
#else
    LOGMP_DEBUG("MPC  runtime : (libmpc not detected)");
    LOGMP_DEBUG("MPC  headers : (libmpc not detected)");
#endif

    LOGMP_DEBUG("------------------------------------------------------------------");
}
