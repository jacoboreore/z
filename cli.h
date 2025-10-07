// cli.h
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
    } mode;

    std::string t_str;      // raw decimal string for t
    std::string step_str;   // only set for Valley Walk
    std::string file_str;   // file path or other string to reuse this var
    std::string job_id;      // job id to identify logs, outputs, graphics for this job id
    
};

bool parse_cli(int argc, char **argv, Options &out);

void print_usage(const char *prog);
