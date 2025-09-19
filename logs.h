// logs.h
#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <string>
#include <mpfr.h>

// Global control flag
inline bool g_main_log = true; // main logger
inline bool g_enable_log_debug = true; // debug logger

inline std::string current_timestamp()
{
    using namespace std::chrono;
    auto now = system_clock::now();
    std::time_t t = system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

// Internal logger (prints one line; msg already formatted)
inline void log_msg_internal(const char *file,
                             int line,
                             const std::string &msg,
                             int /*precision*/ = 15,
                             bool force_log = false,
                             int type = 1) // 1 normal, 2 debug
{
    if (!force_log)
    {
        if (type == 1 && !g_main_log)
        {
            return;
        }

        if (type == 2 && !g_enable_log_debug)
        {
            return;
        }
    }

    std::string module(file);
    if (auto pos = module.find_last_of("/\\"); pos != std::string::npos)
        module = module.substr(pos + 1);

    std::ostringstream oss;
    oss << "[" << current_timestamp() << "] "
        << "[" << module << ":" << line << "] "
        << msg;

    std::cout << oss.str() << std::endl;
}

/* -----------------------------------------------------
   MPFR-aware printf-like logging macros (no va_list).
   Use %Rf / %Re / %Rg with either %.<digits>R? or %.*R?.
   Example:
     LOGMP("t = %.*Rf |Z| = %.*Rg", 24, tval, 20, Ztot);
   ----------------------------------------------------- */

#define LOGMP(fmt, ...)                                         \
    do                                                          \
    {                                                           \
        char *_buf = nullptr;                                   \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);    \
        std::string _s;                                         \
        if (_n >= 0 && _buf)                                    \
        {                                                       \
            _s.assign(_buf, static_cast<size_t>(_n));           \
            mpfr_free_str(_buf);                                \
        }                                                       \
        log_msg_internal(__FILE__, __LINE__, _s, 15, false, 1); \
    } while (0)

#define LOGMPF(fmt, ...)                                       \
    do                                                         \
    {                                                          \
        char *_buf = nullptr;                                  \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);   \
        std::string _s;                                        \
        if (_n >= 0 && _buf)                                   \
        {                                                      \
            _s.assign(_buf, static_cast<size_t>(_n));          \
            mpfr_free_str(_buf);                               \
        }                                                      \
        log_msg_internal(__FILE__, __LINE__, _s, 15, true, 1); \
    } while (0)

#define LOGMP_DEBUG(fmt, ...)                                   \
    do                                                          \
    {                                                           \
        char *_buf = nullptr;                                   \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);    \
        std::string _s;                                         \
        if (_n >= 0 && _buf)                                    \
        {                                                       \
            _s.assign(_buf, static_cast<size_t>(_n));           \
            mpfr_free_str(_buf);                                \
        }                                                       \
        log_msg_internal(__FILE__, __LINE__, _s, 15, false, 2); \
    } while (0)
