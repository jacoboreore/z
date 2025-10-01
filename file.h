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
    // process ID
    inline unsigned long current_pid()
    {
#if defined(_WIN32)
        return static_cast<unsigned long>(GetCurrentProcessId());
#else
        return static_cast<unsigned long>(getpid());
#endif
    }

    inline std::string make_process_filename(const std::string &dir = "output",
                                             const std::string &prefix = "process")
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
        oss << prefix << "_" << tb << "_pid" << current_pid() << ".csv";

        return (fs::path(dir) / oss.str()).string();
    }

    inline std::mutex &csv_mutex()
    {
        static std::mutex m;
        return m;
    }

    // Append a single CSV record: t,absZ,N (all mpfr_t)
    // Thread-safe; opens file, appends, fsync/commit, closes.
    inline bool append_record_csv(const std::string &filepath,
                                  const mpfr_t t,
                                  const mpfr_t absZ,
                                  const mpfr_t N,
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
        std::string line;
        if (write_header_if_new && new_file)
        {
            line += "t,absZ,N\n";
        }

        line += mpfr_to_string(t, digits, 'g');
        line += ',';
        line += mpfr_to_string(absZ, digits, 'g');
        line += ',';
        line += mpfr_to_string(N, digits, 'g');
        line += '\n';

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

    ///////////////////////////////////////////////////////////////////////////////////
    // READING
    ///////////////////////////////////////////////////////////////////////////////////

    // helpers
    inline void trim_inplace(std::string &s)
    {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a])))
            ++a;
        while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1])))
            --b;

        s.assign(s.data() + a, b - a);
    }

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

    // --- main API ------------------------------------------------
    // Calls `on_row(t, absZ, N, line_no)` for each data row.
    // Returns true on full success; false if file open/parse error (continues best-effort).
    template <class F>
    bool read_triplets_csv(
        const std::string &filepath,
        const std::string to_file,
        mpfr_prec_t prec_bits,
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
            on_row(t, absZ, N, line_no, to_file);
        }

        mpfr_clear(t);
        mpfr_clear(absZ);
        mpfr_clear(N);
        // mpfr_free_cache(); // optional: free MPFR per-thread cache

        return ok_all;
    }

} // namespace zcsv
