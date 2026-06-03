// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: logs.h
// Purpose:
//   Minimal, MPFR-aware logging utilities with optional file output and
//   thread-safe emission. Provides three convenience macros:
//
//     - LOGMP(...)       : normal logs (honors g_main_log)
//     - LOGMPF(...)      : forced logs (always emit; ignores g_main_log)
//     - LOGMP_DEBUG(...) : debug logs (honors g_enable_log_debug)
//
//   All three accept printf-like format strings compatible with mpfr_asprintf,
//   so you can print mpfr_t via %Rg / %Rf, etc.
//
// -----------------------------------------------------------------------------

#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <string>
#include <fstream>
#include <mutex>
#include <filesystem>
#include <cstdint>
#include <mpfr.h>

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

// ================== Global controls ==================

/**
 * @brief Master switch for normal logs (used by LOGMP).
 * @note Forced logs (LOGMPF) ignore this flag.
 */
inline bool g_main_log = true;

/**
 * @brief Master switch for debug logs (used by LOGMP_DEBUG).
 */
inline bool g_enable_log_debug = true;

/** @brief Output file stream for logs (if enabled). */
inline std::ofstream g_log_file;

/** @brief When true, write logs to @ref g_log_file (see init_logger). */
inline bool g_log_to_file = false;

/** @brief When true, also echo logs to stdout (besides file). */
inline bool g_log_also_stdout = true;

/** @brief Global mutex to serialize log writes across threads. */
inline std::mutex g_log_mutex;

/** @brief Absolute path of the active log file (if any). */
inline std::string g_log_filename;

// ================== Date/Time helpers ==================

/**
 * @brief Current local timestamp as "YYYY-MM-DD HH:MM:SS".
 * @return Formatted timestamp string.
 * @threadsafe Uses only local temporaries.
 */
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

/**
 * @brief Return current process identifier as unsigned long.
 */
inline unsigned long current_pid()
{
#if defined(_WIN32)
    return static_cast<unsigned long>(GetCurrentProcessId());
#else
    return static_cast<unsigned long>(getpid());
#endif
}

/**
 * @brief Build a log filename under @p dir with prefix + timestamp + pid.
 *
 * Pattern: <prefix>_YYYYMMDD_HHMMSS_pid<NNN>.log
 *
 * @param dir    Output directory (created if missing).
 * @param prefix Filename prefix (e.g., "valley_scanner" or job-id variant).
 * @return Absolute/relative path as string.
 */
inline std::string make_log_filename(const std::string &dir = "logs",
                                     const std::string &prefix = "valley_scanner")
{
    namespace fs = std::filesystem;
    fs::create_directories(dir);

    // timestamp: YYYYMMDD_HHMMSS
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char tb[32]{};
    std::strftime(tb, sizeof(tb), "%Y%m%d_%H%M%S", &tm);

    std::ostringstream oss;
    oss << prefix << "_" << tb << "_pid" << current_pid() << ".log";
    return (fs::path(dir) / oss.str()).string();
}

// ================== Init / shutdown ==================

/**
 * @brief Initialize logging sinks.
 *
 * Opens a new file at logs/<prefix>_YYYYMMDD_HHMMSS_pidXXXX.log and writes an
 * opening banner line. Also configures whether to mirror to stdout.
 *
 * @param dir          Directory for log files (created if missing).
 * @param prefix       Filename prefix (you often use a job-id here).
 * @param also_stdout  When true, echo logs to stdout in addition to file.
 * @return true if file logging is available (file could be opened).
 *
 * @note Even if file logging fails, logs can still go to stdout if
 *       @p also_stdout is true.
 */
inline bool init_logger(const std::string &dir = "logs",
                        const std::string &prefix = "valley_scanner",
                        bool also_stdout = true)
{
    std::lock_guard<std::mutex> lock(g_log_mutex);
    g_log_filename = make_log_filename(dir, prefix);
    g_log_file.open(g_log_filename, std::ios::out | std::ios::trunc);
    g_log_to_file = g_log_file.is_open();
    g_log_also_stdout = also_stdout;
    if (g_log_to_file)
    {
        // First line marker
        g_log_file << "[" << current_timestamp() << "] [logger] opened " << g_log_filename << "\n";
        g_log_file.flush();
    }
    return g_log_to_file;
}

/**
 * @brief Flush and close the active log file (if any).
 * @post g_log_to_file becomes false.
 */
inline void shutdown_logger()
{
    std::lock_guard<std::mutex> lock(g_log_mutex);
    if (g_log_file.is_open())
    {
        g_log_file << "[" << current_timestamp() << "] [logger] closing\n";
        g_log_file.flush();
        g_log_file.close();
    }
    g_log_to_file = false;
}

/**
 * @brief Get the absolute/relative path of the active log file.
 * @return Reference to internal string (valid until next init_logger).
 */
inline const std::string &log_filename() { return g_log_filename; }

// ================== Core logging backend ==================

/**
 * @brief Internal logger used by the public macros.
 *
 * @param file        Caller’s __FILE__.
 * @param line        Caller’s __LINE__.
 * @param msg         Preformatted log message.
 * @param precision   (Reserved) decimal precision hint for future formatters.
 * @param force_log   If true, bypass g_main_log/g_enable_log_debug checks.
 * @param type        1 = normal/info; 2 = debug.
 *
 * @note This routine is thread-safe. It holds the global mutex while
 *       writing to file/stdout, and flushes the file per line to support
 *       real-time tailing.
 */
inline void log_msg_internal(const char *file,
                             int line,
                             const std::string &msg,
                             int /*precision*/ = 15,
                             bool force_log = false,
                             int type = 1)
{
    try
    {
        if (!force_log)
        {
            if (type == 1 && !g_main_log)
                return;
            if (type == 2 && !g_enable_log_debug)
                return;
        }

        const char *safe_file = (file != nullptr) ? file : "<null-file>";
        std::string module(safe_file);

        if (auto pos = module.find_last_of("/\\"); pos != std::string::npos)
            module = module.substr(pos + 1);

        std::ostringstream oss;
        oss << "[" << current_timestamp() << "] "
            << "[" << module << ":" << line << "] "
            << msg;

        const std::string line_out = oss.str();

        std::lock_guard<std::mutex> lock(g_log_mutex);
        if (g_log_to_file && g_log_file.is_open())
        {
            g_log_file << line_out << '\n';
            g_log_file.flush();
        }
        if (g_log_also_stdout)
        {
            std::cout << line_out << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "[LOGGER-FAIL] line=" << line
                  << " err=" << e.what() << '\n';
    }
    catch (...)
    {
        std::cerr << "[LOGGER-FAIL] line=" << line
                  << " err=unknown\n";
    }
}

// ================== MPFR-aware convenience macros ==================

/**
 * @def LOGMP(fmt, ...)
 * @brief Emit a normal log line using mpfr_asprintf-style formatting.
 * @details Honors @ref g_main_log (no output if false).
 *
 * Usage:
 * @code
 * mpfr_t t; mpfr_init2(t, 200);
 * // ...
 * LOGMP("t=%.*Rf |Z|=%.*Rg", 24, t, 24, absZ);
 * @endcode
 */
#define LOGMP(fmt, ...)                                         \
    do                                                          \
    {                                                           \
        char *_buf = nullptr;                                   \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);    \
        std::string _s;                                         \
        if (_n >= 0 && _buf != nullptr)                         \
        {                                                       \
            _s.assign(_buf, static_cast<size_t>(_n));           \
        }                                                       \
        else                                                    \
        {                                                       \
            _s = "<formatting_failed>";                         \
        }                                                       \
        if (_buf != nullptr)                                    \
        {                                                       \
            mpfr_free_str(_buf);                                \
        }                                                       \
        log_msg_internal(__FILE__, __LINE__, _s, 15, false, 1); \
    } while (0)

/**
 * @def LOGMPF(fmt, ...)
 * @brief Emit a forced log line (always prints regardless of switches).
 * @details Ignores @ref g_main_log; still goes to file/stdout per settings.
 */
#define LOGMPF(fmt, ...)                                       \
    do                                                         \
    {                                                          \
        char *_buf = nullptr;                                  \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);   \
        std::string _s;                                        \
        if (_n >= 0 && _buf != nullptr)                        \
        {                                                      \
            _s.assign(_buf, static_cast<size_t>(_n));          \
        }                                                      \
        else                                                   \
        {                                                      \
            _s = "<formatting_failed>";                        \
        }                                                      \
        if (_buf != nullptr)                                   \
        {                                                      \
            mpfr_free_str(_buf);                               \
        }                                                      \
        log_msg_internal(__FILE__, __LINE__, _s, 15, true, 1); \
    } while (0)

/**
 * @def LOGMP_DEBUG(fmt, ...)
 * @brief Emit a debug log line using mpfr_asprintf-style formatting.
 * @details Honors @ref g_enable_log_debug (no output if false).
 */
#define LOGMP_DEBUG(fmt, ...)                                   \
    do                                                          \
    {                                                           \
        char *_buf = nullptr;                                   \
        int _n = mpfr_asprintf(&_buf, (fmt), ##__VA_ARGS__);    \
        std::string _s;                                         \
        if (_n >= 0 && _buf != nullptr)                         \
        {                                                       \
            _s.assign(_buf, static_cast<size_t>(_n));           \
        }                                                       \
        else                                                    \
        {                                                       \
            _s = "<formatting_failed>";                         \
        }                                                       \
        if (_buf != nullptr)                                    \
        {                                                       \
            mpfr_free_str(_buf);                                \
        }                                                       \
        log_msg_internal(__FILE__, __LINE__, _s, 15, false, 2); \
    } while (0)
