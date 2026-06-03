// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: file.h  (zcsv helpers)
// Purpose:
//   Portable CSV append/read helpers for Valley Scanner outputs, with
//   fsync/commit guarantees on write and simple row-wise parsing utilities.
//   Thread-safe appends via a process-global mutex.
//
// -----------------------------------------------------------------------------

#pragma once
#include <string>
#include <mutex>
#include <filesystem>
#include <cerrno>
#include <cstring>

#include <mpfr.h>

#include "utils.h"
#include "logs.h"

#if defined(_WIN32)
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#define ZC_OPEN _open
#define ZC_WRITE _write
#define ZC_CLOSE _close
#define ZC_SYNC _commit
#define ZC_OFLAGS (_O_WRONLY | _O_CREAT | _O_APPEND | _O_BINARY)
#define ZC_MODE (_S_IREAD | _S_IWRITE)
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#define ZC_OPEN ::open
#define ZC_WRITE ::write
#define ZC_CLOSE ::close
#define ZC_SYNC ::fsync
#define ZC_OFLAGS (O_WRONLY | O_CREAT | O_APPEND)
#define ZC_MODE (0644)
#endif

/**
 * @brief Extract human-readable text from an exception-like object.
 * @tparam E exception type (must support .what() or fallback is returned)
 * @return e.what() when available, otherwise "<no what() available>"
 */
template <class E>
const char *ex_text(const E &e)
{
    if constexpr (requires { e.what(); })
        return e.what();
    else
        return "<no what() available>";
}

namespace zcsv
{
    // -------------------------------------------------------------------------
    // Process identity & filenames
    // -------------------------------------------------------------------------

    /**
     * @brief Get current process id as unsigned long (platform-specific).
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
     * @brief Build a per-process output CSV filename in @p dir with @p prefix.
     *
     * Pattern is currently simplified to "<prefix>.csv" (stable name), but the
     * timestamp and pid scaffolding remains for future use if desired.
     *
     * @param dir    Output directory (created if missing).
     * @param prefix Filename stem (no extension).
     * @return Full path to CSV file.
     */
    inline std::string make_process_filename(const std::string &dir = "output",
                                             const std::string &prefix = "process")
    {
        namespace fs = std::filesystem;
        fs::create_directories(dir);

        // timestamp scaffolding retained for optional use
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
        // oss << prefix << "_" << tb << "_pid" << current_pid() << ".csv";
        oss << prefix << ".csv"; // simple file name

        return (fs::path(dir) / oss.str()).string();
    }

    /**
     * @brief Global CSV mutex for append operations.
     * @return Reference to a process-global mutex.
     */
    inline std::mutex &csv_mutex()
    {
        static std::mutex m;
        return m;
    }

    // -------------------------------------------------------------------------
    // APPEND (single-record, fsync-on-write)
    // -------------------------------------------------------------------------

    /**
     * @brief Append one CSV line in a thread-safe, durable manner.
     *
     * Modes and headers:
     *   - mode=1 (walking)                 → header_walking            → cols: t,absZ,N
     *   - mode=2 (walking for plotting)    → header_walking_for_plotting → cols: t,absZ,N,minima
     *   - mode=3 (refine/progressive out)  → header_refined            → cols: t,absZ,spacing
     *
     * The header strings come from CFG() fields (e.g., header_walking).
     *
     * @param filepath  Target CSV path (parent dir created if needed).
     * @param t         MPFR value to print.
     * @param absZ      |Z|(t) or Z(t) magnitude .
     * @param N         N(t) or counting proxy (ignored in refine mode).
     * @param spacing   Δt vs previous refined t (used only in mode=3).
     * @param minima    For plotting: 1 = last |Z| before bounce, 2 = refined t.
     * @param mode      1=walking, 2=walking_for_plotting, 3=refined.
     * @param digits    Significant digits for mpfr_to_string (default 24).
     * @param write_header_if_new  Emit header if file does not yet exist.
     * @return true on success (write+sync completed), false on error.
     *
     * @threadsafe Uses @ref csv_mutex to serialize concurrent appends.
     * @warning mpfr_to_string must be available from utils.h as used here.
     */
    inline bool append_record_csv(
        const std::string &filepath,
        const mpfr_t t,
        const mpfr_t absZ,
        const mpfr_t N,
        const mpfr_t spacing,
        const int minima, // 1 -> last |Z| before bouncing, 2 -> refined t
        const int mode,   // 1 walking, 2 walking for plotting, 3 refine
        int digits = 24,
        bool write_header_if_new = true)
    {
        std::lock_guard<std::mutex> lock(csv_mutex());

        // Ensure parent directory exists
        try
        {
            std::filesystem::path p(filepath);
            if (p.has_parent_path())
            {
                std::filesystem::create_directories(p.parent_path());
            }
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            LOGMPF("Error writing %s:%d — %s", __func__, __LINE__, ex_text(e));
        }

        const bool new_file = !std::filesystem::exists(filepath);

        // Build the CSV line in memory first (single write → no interleaving).
        // header_refined = t,absZ,spacing
        // header_walking = t,absZ,N
        // header_walking_for_plotting = t,absZ,N,minima

        std::string line;
        if (write_header_if_new && new_file)
        {
            switch (mode)
            {
            case 1:
                // header_walking
                line += CFG().header_walking + "\n";
                break;
            case 2:
                // header_walking_for_plotting
                line += CFG().header_walking_for_plotting + "\n";
                break;
            case 3:
                // refine
                line += CFG().header_refined + "\n";
                break;
            }
            // line += "t,absZ,N\n";
        }

        switch (mode)
        {
        case 1:
            // walking
            line += mpfr_to_string(t, digits, 'g');
            line += ',';
            line += mpfr_to_string(absZ, digits, 'g');
            line += ',';
            line += mpfr_to_string(N, digits, 'g');
            line += '\n';
            break;
        case 2:
            // walking_for_plotting
            line += mpfr_to_string(t, digits, 'g');
            line += ',';
            line += mpfr_to_string(absZ, digits, 'g');
            line += ',';
            line += mpfr_to_string(N, digits, 'g');
            line += ',';
            line += std::to_string(minima);
            line += '\n';
            break;
        case 3:
            // refine
            line += mpfr_to_string(t, digits, 'g');
            line += ',';
            line += mpfr_to_string(absZ, digits, 'g');
            line += ',';
            line += mpfr_to_string(spacing, digits, 'g');
            line += '\n';
            break;
        }

        // Open, write, sync, close
        int fd = ZC_OPEN(filepath.c_str(), ZC_OFLAGS, ZC_MODE);
        if (fd < 0)
        {
            // optional: print or handle errno
            return false;
        }

        const char *data = line.data();
        size_t left = line.size();
        while (left > 0)
        {
            int wrote = ZC_WRITE(fd, data, static_cast<unsigned int>(left));
            if (wrote < 0)
            {
                ZC_CLOSE(fd);
                return false;
            }
            data += wrote;
            left -= wrote;
        }

        // Flush to disk for durability
        if (ZC_SYNC(fd) != 0)
        {
            // even if sync fails, still close and report failure
            ZC_CLOSE(fd);
            return false;
        }

        ZC_CLOSE(fd);
        return true;
    }

    // -------------------------------------------------------------------------
    // READING
    // -------------------------------------------------------------------------

    /** @brief Trim leading/trailing ASCII whitespace in-place. */
    inline void trim_inplace(std::string &s)
    {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a])))
            ++a;
        while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1])))
            --b;

        s.assign(s.data() + a, b - a);
    }

    /**
     * @brief Split a CSV line with three fields into a,b,c (no quotes/escapes).
     * @return true on success; false if fewer than 3 tokens are present.
     * @note Each field is trimmed.
     */
    inline bool split3(const std::string &line, std::string &a, std::string &b, std::string &c)
    {
        std::istringstream iss(line);

        if (!std::getline(iss, a, ','))
            return false;
        if (!std::getline(iss, b, ','))
            return false;
        if (!std::getline(iss, c))
            return false;

        trim_inplace(a);
        trim_inplace(b);
        trim_inplace(c);

        return true;
    }

    /**
     * @brief Parse base-10 string into mpfr_t, with fallback for 'e'→'@' exponent.
     * @param out Destination mpfr_t (precision already set).
     * @param s   Source string (e.g., "1.23e-4").
     * @param rnd MPFR rounding mode.
     * @return true if parse succeeded, false otherwise.
     */
    inline bool mpfr_from_str10(mpfr_t out, const std::string &s, mpfr_rnd_t rnd = MPFR_RNDN)
    {
        if (mpfr_set_str(out, s.c_str(), 10, rnd) == 0)
            return true;
        // fallback: replace 'e'/'E' with '@' (MPFR also accepts mantissa@exp)
        std::string t = s;
        for (char &ch : t)
            if (ch == 'e' || ch == 'E')
                ch = '@';
        return mpfr_set_str(out, t.c_str(), 10, rnd) == 0;
    }

    /**
     * @brief Read a CSV of triplets (t, |Z|, N) and invoke a callback per row.
     *
     * Format:
     *   - Optional header auto-detected if first token is "t" or "T".
     *   - Each subsequent non-empty, non-comment line is split by ',' into
     *     three fields and parsed as base-10 numbers into (t, absZ, N).
     *
     * Callback signature:
     * @code
     *   on_row(t, absZ, N, line_no, to_file, spacing, prev_t);
     * @endcode
     * The callback may append refined rows to @p to_file using append_record_csv.
     *
     * @tparam F        Callable type (will be forwarded).
     * @param filepath  Input CSV path.
     * @param to_file   Output CSV path (used by callback).
     * @param spacing   MPFR var passed through for Δt tracking (callback owns logic).
     * @param prev_t    MPFR var holding previous t* (callback updates it).
     * @param prec_bits Precision for per-line temporaries (t, absZ, N).
     * @param on_row    User callback invoked per parsed row.
     * @param has_header_auto Auto-detect and skip header row if present.
     * @return true if file opened and all parsed rows were processed; false if
     *         open/parse errors occurred (processing continues best-effort).
     *
     * @note Internally reuses mpfr temporaries per line; caller owns @p spacing
     *       and @p prev_t lifetimes.
     */
    template <class F>
    bool read_triplets_csv(
        const std::string &filepath,
        const std::string to_file,
        mpfr_t spacing,
        mpfr_t prev_t,
        const mpfr_prec_t prec_bits,
        F &&on_row,
        bool has_header_auto = true)
    {
        std::ifstream in(filepath);

        if (!in)
        {
            LOGMPF("Error opening the file");
            return false;
        }

        std::string line;
        size_t line_no = 0;
        bool ok_all = true;

        // temp mpfr vars (reused per line)
        mpfr_t t, absZ, N;
        mpfr_init2(t, prec_bits);
        mpfr_init2(absZ, prec_bits);
        mpfr_init2(N, prec_bits);

        // Optional: autodetect and skip header
        if (has_header_auto)
        {
            std::streampos pos = in.tellg();
            if (std::getline(in, line))
            {
                ++line_no;
                std::string a, b, c;
                if (split3(line, a, b, c) &&
                    (a == "t" || a == "T"))
                {
                    // header detected: do nothing (already consumed)
                }
                else
                {
                    // not a header: rewind
                    in.clear();
                    in.seekg(pos);
                    --line_no;
                }
            }
        }

        while (std::getline(in, line))
        {
            ++line_no;
            // skip empty/comment
            std::string s = line;
            trim_inplace(s);
            if (s.empty() || s[0] == '#')
                continue;

            std::string a, b, c;
            if (!split3(s, a, b, c))
            {
                ok_all = false;
                continue;
            }

            if (!mpfr_from_str10(t, a) ||
                !mpfr_from_str10(absZ, b) ||
                !mpfr_from_str10(N, c))
            {
                ok_all = false; // malformed number, skip row
                continue;
            }

            // user callback
            on_row(t, absZ, N, line_no, to_file, spacing, prev_t);
        }

        mpfr_clear(t);
        mpfr_clear(absZ);
        mpfr_clear(N);
        // mpfr_free_cache(); // optional: free MPFR per-thread cache

        return ok_all;
    }

} // namespace zcsv
