#pragma once
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>
#include <string>

// Global control flag
inline bool g_enable_log = true;

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

// Internal logger
inline void log_msg_internal(
    const char *file,
    int line,
    const std::string &msg,
    int precision = 15,
    bool force_log = false)
{
    if (!g_enable_log && !force_log)
        return;

    // Module name from filename
    std::string module(file);
    auto slash = module.find_last_of("/\\");
    if (slash != std::string::npos)
        module = module.substr(slash + 1);

    std::ostringstream oss;
    oss << "[" << current_timestamp() << "] "
        << "[" << module << ":" << line << "] "
        << std::fixed << std::setprecision(precision)
        << msg;

    std::cout << oss.str() << std::endl;
}

// Macros
#define LOG(msg) log_msg_internal(__FILE__, __LINE__, (msg))
#define LOGF(msg) log_msg_internal(__FILE__, __LINE__, (msg), 15, true)
#define LOGP(msg, prec) log_msg_internal(__FILE__, __LINE__, (msg), (prec))
#define LOGPF(msg, prec) log_msg_internal(__FILE__, __LINE__, (msg), (prec), true)
