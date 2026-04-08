// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: cli.cpp
// Purpose:
//   Command-line parser for Valley Scanner. Translates user arguments into an
//   Options structure, performing basic validation and help printing.
//
// Supported modes:
//   • z <t> <job_id>
//   • valley_walk <t> <num_zeros> <job_id>
//   • refine <t> <job_id>
//   • pre_filter <t> <job_id>
//   • valley_scanner <t> <job_id>
//   • refine_progressive <file> <job_id>
//   • z_os <t> <job_id>
//
// -----------------------------------------------------------------------------

#include <mpfr.h>
#include <algorithm>
#include <cctype>
#include <iostream>

#include "cli.h"

int local_argc;

/**
 * @brief Validate that a string represents a finite numeric value.
 *
 * @param s Input C-string to test.
 * @return true if successfully parsed as a number by MPFR, false otherwise.
 *
 * @details
 *   Uses mpfr_set_str with 64-bit precision and base 10.
 *   Checks both parse success and mpfr_number_p() (filters NaN/Inf).
 */
static bool mpfr_num_ok(const char *s)
{
    mpfr_t tmp;
    mpfr_init2(tmp, 64);
    int ok = (mpfr_set_str(tmp, s, 10, MPFR_RNDN) == 0) && mpfr_number_p(tmp);
    mpfr_clear(tmp);
    return ok;
}

/**
 * @brief Print full usage instructions for all supported modes.
 *
 * @param prog Executable name (argv[0]).
 *
 * @details
 *   Prints a concise command reference for all CLI modes:
 *     • z — Compute a single Z(t) value.
 *     • z_ball — Variant of z using ball arithmetic or alternate evaluation.
 *     • pre_filter — Run preliminary filter at a given t.
 *     • valley_scanner — Initiate valley scan around t.
 *     • valley_walk — Walk valleys starting from t for N zeros.
 *     • refine — Refine a previously found zero at t.
 *     • refine_progressive — Refine using a CSV or dataset file.
 *     • z_os — Odlyzko/Schönhage Z(t) evaluation mode.
 *
 *   Each mode requires specific arguments, detailed below.
 */
void print_usage(const char *prog)
{
    std::cerr
        << "Received: " << local_argc << " arguments\n"
        << "Usage:\n"
        << "  " << prog << " z <t> <job_id>\n"
        << "      Compute Z(t) at a single point.\n\n"

        << "  " << prog << " z_ball <t> <job_id>\n"
        << "      Compute Z(t) using alternative (ball-based) evaluation.\n\n"

        << "  " << prog << " pre_filter <t>\n"
        << "      Perform a pre-filter evaluation at t (quick test).\n\n"

        << "  " << prog << " valley_scanner <t> <job_id>\n"
        << "      Scan a small region around t for candidate zeros.\n\n"

        << "  " << prog << " valley_walk <t> <num_zeros> <job_id>\n"
        << "      Walk along successive valleys starting from t for the given number of zeros.\n\n"

        << "  " << prog << " refine <t> <job_id>\n"
        << "      Refine a previously estimated zero around t.\n\n"

        << "  " << prog << " refine_progressive <file_path> <job_id>\n"
        << "      Perform a progressive refinement pass using an existing CSV dataset.\n\n"

        << "  " << prog << " z_os <t> <job_id>\n"
        << "      Compute Z(t) using the Odlyzko–Schönhage algorithm.\n\n"

        << "  " << prog << " z_distributed_precalc <t> <job_id>\n"
        << "      Compute precalc for distributed z.\n\n"

        << "  " << prog << " z_distributed_chunk <file> <chunk_id> <start> <end> <job_id>\n"
        << "      Compute precalc for distributed z.\n\n"
        
        << "  " << prog << " z_distributed_refine_stage1 <t> <job_id> <output>\n"
        << "      Compute stage 1 for refinement z.\n\n";

        
}

/**
 * @brief Parse the command-line arguments into an Options struct.
 *
 * @param argc Standard argument count from main().
 * @param argv Standard argument vector from main().
 * @param[out] out Filled on success with parsed fields.
 * @return true on successful parse, false on syntax or numeric errors.
 *
 * @details
 *   - Mode string is case-insensitive (converted to lowercase).
 *   - Checks argument count per subcommand.
 *   - Validates numeric parameters with mpfr_num_ok().
 *   - Always sets job_id from the last argument.
 *   - Prints errors to stderr on invalid usage.
 */
bool parse_cli(int argc, char **argv, Options &out)
{
    // print the arguments:
    std::cout << "Program arguments:" << std::endl;

    // Loop from index 0 to argc - 1
    for (int i = 0; i < argc; i++)
    {
        // argv[0] is typically the name of the executable itself
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
    }

    std::cout << "parsing... \n ";

    if (argc < 3)
    {
        print_usage(argv[0]);
        return false;
    }

    local_argc = argc;

    // Normalize mode to lowercase for comparison
    std::string mode = argv[1];
    std::transform(mode.begin(), mode.end(), mode.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });

    // ---------------------- Z / Z_BALL ----------------------
    if (mode == "z" || mode == "z_ball")
    {
        if (argc != 4)
        {
            std::cerr << "[error] z requires exactly 1 number (t) and a string (job_id)\n";
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
    }
    // ---------------------- PREFILTER -----------------------
    else if (mode == "pre_filter")
    {
        if (argc != 3)
        {
            std::cerr << "[error] pre_filter requires exactly 1 number (t)\n";
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
    }
    // ---------------------- VALLEY_SCANNER ------------------
    else if (mode == "valley_scanner")
    {
        if (argc != 4)
        {
            std::cerr << "[error] valley_scanner requires 1 number (t) and a string (job_id)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::Valley_Scanner;
        out.t_str = argv[2];
    }
    // ---------------------- VALLEY_WALK ---------------------
    else if (mode == "valley_walk")
    {
        if (argc != 5)
        {
            std::cerr << "[error] valley_walk requires 2 numbers (t and number of zeros) and a string (job_id)\n";
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

        out.mode = Options::Mode::Valley_Walk;
        out.t_str = argv[2];
        out.step_str = argv[3];
    }
    // ---------------------- REFINE --------------------------
    else if (mode == "refine")
    {
        if (argc != 4)
        {
            std::cerr << "[error] refine requires exactly 1 number (t) and a string (job_id)\n";
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
    }
    // ---------------------- REFINE_PROGRESSIVE --------------
    else if (mode == "refine_progressive")
    {
        if (argc != 4)
        {
            std::cerr << "[error] refine_progressive requires exactly 1 string (file path) and a string (job_id)\n";
            print_usage(argv[0]);
            return false;
        }

        out.mode = Options::Mode::Refine_Progressive;
        out.file_str = argv[2];
        out.step_str.clear();
    }
    // ---------------------- Z_OS -----------------------------
    else if (mode == "z_os")
    {
        if (argc != 4)
        {
            std::cerr << "[error] z_os requires 1 number (t) and a string (job_id)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::Z_OS;
        out.t_str = argv[2];
    }
    // ---------------------- Z_Distributed_Precalc -----------------------------
    else if (mode == "z_distributed_precalc")
    {
        if (argc != 4)
        {
            std::cerr << "[error] z_distributed_precalc requires 1 number (t) and a string (job_id)\n";
            print_usage(argv[0]);
            return false;
        }

        if (!mpfr_num_ok(argv[2]))
        {
            std::cerr << "[error] invalid number for t: " << argv[2] << "\n";
            return false;
        }

        out.mode = Options::Mode::Z_Distributed_Precalc;
        out.t_str = argv[2];
    }
    // ---------------------- Z_Distributed_Chunk -----------------------------
    else if (mode == "z_distributed_chunk")
    {
        if (argc != 8)
        {
            std::cerr << "[error] z_distributed_chunk requires 1 string (file), 3 numbers (chunk_id, start, end), a string chunk_file, and a string (job_id)\n";
            print_usage(argv[0]);
            return false;
        }

        out.mode = Options::Mode::Z_Distributed_Chunk;
        out.t_str = argv[2];
    }
    // ---------------------- z_distributed_refine_stage1 -----------------------------
    else if (mode == "z_distributed_refine_stage1")
    {
        if (argc != 5)
        {
            std::cerr << "[error] z_distributed_refine_stage1 requires  1 number (t),a string (job_id), an output file (json to python for this stage) \n";
            print_usage(argv[0]);
            return false;
        }

        out.mode = Options::Mode::Z_Distributed_Refine_Stage1;
        out.t_str = argv[2];
    }
    // ---------------------- z_distributed_refine_stage2 -----------------------------
    else if (mode == "z_distributed_refine_stage2")
    {
        if (argc != 4)
        {
            std::cerr << "[error] z_distributed_refine_stage2 requires a string (input_json), a string (output_json) \n";
            print_usage(argv[0]);
            return false;
        }

        out.mode = Options::Mode::Z_Distributed_Refine_Stage2;
    }
    // ---------------------- UNKNOWN --------------------------
    else
    {
        std::cerr << "[error] unknown mode: '" << argv[1] << "'\n";
        print_usage(argv[0]);
        return false;
    }

    // ---------------------- COMMON ---------------------------
    out.job_id = argv[argc - 1];
    return true;
}
