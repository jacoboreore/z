// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: crash_handlers.hpp
// Purpose:
//   Lightweight cross-platform crash/exception signal handling.
//   • Installs handlers for common fatal signals (SIGSEGV, SIGABRT, SIGFPE, SIGILL).
//   • On POSIX (Unix/macOS), prints a best-effort backtrace to stderr.
//   • Exits the process with 128 + signal code (conventional shell status).
//
// -----------------------------------------------------------------------------

#pragma once

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(__unix__) || defined(__APPLE__)
#include <execinfo.h>
#endif

/**
 * @brief Emit a best-effort stack backtrace to stderr.
 * @details
 *   POSIX: uses backtrace()/backtrace_symbols_fd() from execinfo.h.
 *   Non-POSIX: prints a placeholder message (no backtrace available).
 * @note Intended for use inside a signal handler; avoids heap allocations.
 */
inline void os_dump_stack()
{
#if defined(__unix__) || defined(__APPLE__)
    void *bt[64];
    int n = ::backtrace(bt, 64);
    std::fprintf(stderr, "\n---- BACKTRACE ----\n");
    ::backtrace_symbols_fd(bt, n, fileno(stderr));
    std::fprintf(stderr, "-------------------\n");
#else
    std::fprintf(stderr, "(no backtrace on this platform)\n");
#endif
}

/**
 * @brief Resolve a textual signal name for logging.
 * @param sig Numeric signal identifier (e.g., SIGSEGV).
 * @return Human-readable signal name (best effort; never null).
 * @details
 *   POSIX: uses strsignal().
 *   Non-POSIX: returns a small static mapping for common signals.
 */
inline const char *os_signal_name(int sig)
{
#if defined(__unix__) || defined(__APPLE__)
    const char *name = ::strsignal(sig); // global namespace
    return name ? name : "unknown";
#else
    // Minimal fallback for non-POSIX platforms
    switch (sig)
    {
    case SIGSEGV:
        return "SIGSEGV";
    case SIGABRT:
        return "SIGABRT";
    case SIGFPE:
        return "SIGFPE";
    case SIGILL:
        return "SIGILL";
    case SIGINT:
        return "SIGINT";
    case SIGTERM:
        return "SIGTERM";
    default:
        return "signal";
    }
#endif
}

/**
 * @brief Common crash signal handler: prints signal & backtrace, then terminates.
 * @param sig Received signal number.
 * @note Uses std::_Exit(128+sig) to avoid running static destructors/re-entry.
 */
inline void os_signal_handler(int sig)
{
    const char *name = os_signal_name(sig);
    std::fprintf(stderr, "\n[FATAL] signal %d (%s) received\n", sig, name);
    os_dump_stack();
    std::_Exit(128 + sig);
}

/**
 * @brief Install crash handlers for common fatal signals.
 * @details
 *   POSIX: installs via sigaction with SA_RESETHAND for SIGSEGV/SIGABRT/SIGFPE/SIGILL.
 *   Others: uses std::signal for the same set.
 * @note Call this early in main() to capture faults during initialization.
 */
inline void install_crash_handlers()
{
#if defined(__unix__) || defined(__APPLE__)
    // Prefer sigaction on POSIX for reliable behavior
    struct sigaction sa{};
    sa.sa_handler = os_signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESETHAND;
    ::sigaction(SIGSEGV, &sa, nullptr);
    ::sigaction(SIGABRT, &sa, nullptr);
    ::sigaction(SIGFPE, &sa, nullptr);
    ::sigaction(SIGILL, &sa, nullptr);
#else
    std::signal(SIGSEGV, os_signal_handler);
    std::signal(SIGABRT, os_signal_handler);
    std::signal(SIGFPE, os_signal_handler);
    std::signal(SIGILL, os_signal_handler);
#endif
}
