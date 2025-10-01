// logs.h
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
inline bool g_main_log = true;         // main logger on/off
inline bool g_enable_log_debug = true; // debug logger on/off

// File logging state
inline std::ofstream g_log_file;
inline bool g_log_to_file = false;    // write to file?
inline bool g_log_also_stdout = true; // echo to stdout too?
inline std::mutex g_log_mutex;
inline std::string g_log_filename;

// Date-Time helper
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

// Filename helpers, process ID to the name plus the current date
inline unsigned long current_pid()
{
#if defined(_WIN32)
    return static_cast<unsigned long>(GetCurrentProcessId());
#else
    return static_cast<unsigned long>(getpid());
#endif
}

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

// Init / shutdown
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

// expose the filename
inline const std::string &log_filename() { return g_log_filename; }

// Core logging
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
            return;
        if (type == 2 && !g_enable_log_debug)
            return;
    }

    std::string module(file);
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
        g_log_file.flush(); // keep immediate for agent tailing
    }
    if (g_log_also_stdout)
    {
        std::cout << line_out << std::endl;
    }
}

// ================== MPFR-aware macros ==================
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
