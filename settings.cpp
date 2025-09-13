#include <mpfr.h>
#include <mpc.h>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <mutex>
#include "settings.h"

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
    else if (ieq(k, "valley_shrink"))
        c.valley_shrink = std::stod(v);
    else if (ieq(k, "valley_step_min"))
        c.valley_step_min = std::stod(v);
    else if (ieq(k, "valley_step_max"))
        c.valley_step_max = std::stod(v);
}

static void load_file(Config &cfg, const std::string &path)
{

    std::ifstream f;

    std::cout << "Config file: " << path;

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
    std::cout.setf(std::ios::boolalpha);
    std::cout
        << "mpfr_prec               = " << c.mpfr_prec << "\n"
        << "prefilter_alpha         = " << c.prefilter_alpha << "\n"
        << "prefilter_tau           = " << c.prefilter_tau << "\n"
        << "prefilter_tau_boost     = " << c.prefilter_tau_boost << "\n"
        << "newton_tol_abs          = " << c.newton_tol_abs << "\n"
        << "newton_tol_step         = " << c.newton_tol_step << "\n"
        << "newton_max_iter         = " << c.newton_max_iter << "\n"
        << "use_threaded_z          = " << c.use_threaded_z << "\n"
        << "use_threaded_newton     = " << c.use_threaded_newton << "\n"
        << "ball_mid_threshold      = " << c.ball_mid_threshold << "\n"
        << "h_min                   = " << c.h_min << "\n"
        << "h_frac_spacing          = " << c.h_frac_spacing << "\n"
        << "step_frac_spacing       = " << c.step_frac_spacing << "\n"
        << "rs_half_weight_last     = " << c.rs_half_weight_last << "\n"
        << "rs_tail_C               = " << c.rs_tail_C << "\n"
        << "min_threads             = " << c.min_threads << "\n"
        << "max_threads             = " << c.max_threads << "\n"
        << "rs_min_N_for_threads    = " << c.rs_min_N_for_threads << "\n"
        << "theta_exact_switch_t    = " << c.theta_exact_switch_t << "\n"
        << "rs_log_debug            = " << c.rs_log_debug << "\n"
        << "rs_log_digits           = " << c.rs_log_digits << "\n"
        << "theta_stirling_terms    = " << c.theta_stirling_terms << "\n"
        << "valley_grow             = " << c.valley_grow << "\n"
        << "valley_shrink           = " << c.valley_shrink << "\n"
        << "valley_step_min         = " << c.valley_step_min << "\n"
        << "valley_step_max         = " << c.valley_step_max << "\n"

        << "------------------------------------------------------------------\n";
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
        load_env(g_cfg);
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
    std::cout << "MPFR runtime : " << mpfr_get_version() << "\n";
#ifdef MPFR_VERSION_STRING
    std::cout << "MPFR headers : " << MPFR_VERSION_STRING << "\n";
#else
    std::cout << "MPFR headers : (no MPFR_VERSION_STRING)\n";
#endif
    std::cout << "GMP  runtime : " << gmp_version << "\n";
    std::cout << "GMP  headers : " << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << "\n";
#ifdef mpc_get_version
    std::cout << "MPC  runtime : " << mpc_get_version() << "\n";
#else
    std::cout << "MPC  runtime : (no mpc_get_version)\n";
#endif
#ifdef MPC_VERSION_STRING
    std::cout << "MPC  headers : " << MPC_VERSION_STRING << "\n";
#else
    std::cout << "MPC  headers : (no MPC_VERSION_STRING)\n";
#endif
    std::cout << "------------------------------------------------------------------\n";
}
