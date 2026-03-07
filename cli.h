// -----------------------------------------------------------------------------
// Valley Scanner © 2025 Jacob Orellana
// All rights reserved. Redistribution and use in source or binary forms
// are permitted only with explicit written permission of the author.
//
// File: cli.h
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

#pragma once
#include <string>

struct Options
{
    enum class Mode
    {
        Z,
        Z_ball,
        Z_OS,
        Valley_Scanner,
        Valley_Walk,
        Refine,
        Refine_Progressive,
        Pre_filter,
        Z_Distributed_Precalc,
        Z_Distributed_Chunk
    } mode;

    std::string t_str;      // raw decimal string for t
    std::string step_str;   // only set for Valley Walk
    std::string file_str;   // file path or other string to reuse this var
    std::string job_id;      // job id to identify logs, outputs, graphics for this job id
    
};

bool parse_cli(int argc, char **argv, Options &out);

void print_usage(const char *prog);
