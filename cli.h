// cli.h
#pragma once
#include <string>

struct Options
{
    enum class Mode
    {
        Z,
        Z_ball,
        z_OS,
        ValleyScanner,
        Refine,
        Pre_filter,
    } mode;

    std::string t_str;    // raw decimal string for t
    std::string step_str; // only set for ValleyScanner
};

bool parse_cli(int argc, char **argv, Options &out);

void print_usage(const char *prog);
