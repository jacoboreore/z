#pragma once
#include <mpfr.h>
#include <string>
#include <sstream>
#include <mutex>

// --- Best-zero global state (singleton) ---
struct BestZeroState {
    mpfr_t best_absZ;   // |Z| of the best-so-far
    mpfr_t best_t;      // corresponding t
    mpfr_prec_t P;      // precision used for the mpfr_t's
    bool initialized;
    BestZeroState() : P(0), initialized(false) {}
};

BestZeroState& best_state();           // access singleton
void best_init(mpfr_prec_t P);         // call once at startup (or when precision changes)
void best_reset();                     // set |Z|=9'999'999 and t=0
void best_update(mpfr_srcptr Z, mpfr_srcptr t); // compare |Z| and update if better
void best_get(mpfr_t out_t, mpfr_t out_absZ);   // copy out a snapshot (thread-safe)

inline std::string mpfr_to_string(mpfr_srcptr x, int digits, char mode = 'g')
{
    // mode: 'f' fixed, 'e' scientific, 'g' general
    char fmt[8]; // e.g. "%.*Rg"
    fmt[0] = '%';
    fmt[1] = '.';
    fmt[2] = '*';
    fmt[3] = 'R';
    fmt[4] = mode;
    fmt[5] = '\0';

    int n = mpfr_snprintf(nullptr, 0, fmt, digits, x);
    std::string s(n + 1, '\0'); // +1 for null
    mpfr_snprintf(s.data(), s.size(), fmt, digits, x);
    s.pop_back(); // drop the null
    return s;
}

// Stream manipulator: usage -> oss << mp(tval, 30, 'g')
struct MpfrPrintable
{
    mpfr_srcptr x;
    int digits;
    char mode;
};

inline MpfrPrintable mp(mpfr_srcptr x, int digits, char mode = 'g') { return {x, digits, mode}; }
inline std::ostream &operator<<(std::ostream &os, const MpfrPrintable &m)
{
    return os << mpfr_to_string(m.x, m.digits, m.mode);
}
