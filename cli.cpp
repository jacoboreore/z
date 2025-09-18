#include <mpfr.h>
#include <algorithm>
#include <cctype>
#include <iostream>
#include "cli.h"

static bool mpfr_num_ok(const char *s)
{
    mpfr_t tmp;

    mpfr_init2(tmp, 64);
    int ok = (mpfr_set_str(tmp, s, 10, MPFR_RNDN) == 0) && mpfr_number_p(tmp);
    mpfr_clear(tmp);

    return ok;
}

void print_usage(const char *prog)
{
    std::cerr
        << "Usage:\n"
        << "  " << prog << " z <t>\n"
        << "  " << prog << " valley_scanner <t> <count>\n"
        << "  " << prog << " refine <t>\n\n"
        << "Args:\n"
        << "  <t>     : real number (decimal or scientific), the imaginary part\n"
        << "  <count>  : number of zeros to calculate\n";
}

bool parse_cli(int argc, char **argv, Options &out)
{
    if (argc < 3)
    {
        print_usage(argv[0]);
        return false;
    }

    // normalize mode to lowercase
    std::string mode = argv[1];

    std::transform(
        mode.begin(),
        mode.end(),
        mode.begin(),
        [](unsigned char c)
        {
            return std::tolower(c);
        });

    if (mode == "z" || mode == "z_ball")
    {
        if (argc != 3)
        {
            std::cerr << "[error] z requires exactly 1 number (t)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = (mode == "z") ? Options::Mode::Z : Options::Mode::Z_ball;
        out.t_str = argv[2];
        out.step_str.clear();

        return true;
    }
    else if (mode == "pre_filter")
    {

        if (argc != 3)
        {
            std::cerr << "[error] z requires exactly 1 number (t)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::Pre_filter;
        out.t_str = argv[2];
        out.step_str.clear();

        return true;
    }
    else if (mode == "valley_scanner")
    {
        if (argc != 4)
        {
            std::cerr << "[error] valley_scanner requires 2 numbers: <t> <step>\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        if (!mpfr_num_ok(argv[3]))
        {
            std::cerr << "[error] invalid number for step: " << argv[3] << "\n";
            return false;
        }

        out.mode = Options::Mode::ValleyScanner;
        out.t_str = argv[2];
        out.step_str = argv[3];

        return true;
    }
    else if (mode == "refine")
    {
        if (argc != 3)
        {
            std::cerr << "[error] refine requires exactly 1 number (t)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::Refine;
        out.t_str = argv[2];
        out.step_str.clear();

        return true;
    }
    else if (mode == "z_os")
    {
        if (argc != 3)
        {
            std::cerr << "[error] z_os Odlyzko/Schönhage requires 1 numbers: <t> <step>\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::z_OS;
        out.t_str = argv[2];

        return true;
    }
    else
    {
        std::cerr << "[error] unknown mode: " << argv[1] << "\n";
        print_usage(argv[0]);
        return false;
    }
}
