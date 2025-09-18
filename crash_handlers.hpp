#pragma once
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>   // <-- for ::strsignal on POSIX

#if defined(__unix__) || defined(__APPLE__)
  #include <execinfo.h>
#endif

inline void os_dump_stack() {
#if defined(__unix__) || defined(__APPLE__)
    void* bt[64];
    int n = ::backtrace(bt, 64);
    std::fprintf(stderr, "\n---- BACKTRACE ----\n");
    ::backtrace_symbols_fd(bt, n, fileno(stderr));
    std::fprintf(stderr, "-------------------\n");
#else
    std::fprintf(stderr, "(no backtrace on this platform)\n");
#endif
}

inline const char* os_signal_name(int sig) {
#if defined(__unix__) || defined(__APPLE__)
    const char* name = ::strsignal(sig);  // global namespace
    return name ? name : "unknown";
#else
    // Minimal fallback for non-POSIX platforms
    switch (sig) {
        case SIGSEGV: return "SIGSEGV";
        case SIGABRT: return "SIGABRT";
        case SIGFPE:  return "SIGFPE";
        case SIGILL:  return "SIGILL";
        case SIGINT:  return "SIGINT";
        case SIGTERM: return "SIGTERM";
        default:      return "signal";
    }
#endif
}

inline void os_signal_handler(int sig) {
    const char* name = os_signal_name(sig);
    std::fprintf(stderr, "\n[FATAL] signal %d (%s) received\n", sig, name);
    os_dump_stack();
    std::_Exit(128 + sig);
}

inline void install_crash_handlers() {
#if defined(__unix__) || defined(__APPLE__)
    // Prefer sigaction on POSIX for reliable behavior
    struct sigaction sa{};
    sa.sa_handler = os_signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESETHAND;
    ::sigaction(SIGSEGV, &sa, nullptr);
    ::sigaction(SIGABRT, &sa, nullptr);
    ::sigaction(SIGFPE,  &sa, nullptr);
    ::sigaction(SIGILL,  &sa, nullptr);
#else
    std::signal(SIGSEGV, os_signal_handler);
    std::signal(SIGABRT, os_signal_handler);
    std::signal(SIGFPE,  os_signal_handler);
    std::signal(SIGILL,  os_signal_handler);
#endif
}
