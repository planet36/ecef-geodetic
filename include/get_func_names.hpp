// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// Get the names of the algorithms to test
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "map_func_name_to_func_info.hpp"

#include <algorithm>
#if defined(DEBUG)
#include <cassert>
#endif
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ranges>
#include <string>
#include <vector>

/// get the algorithms named in \a argv from \a first_arg on, or all of them if none are named
/**
* The program exits with a list of the valid names if any name is not a key of
* \c map_func_name_to_func_info.  It also exits if a name is given more than once.
* \param argc the argument count from \c main
* \param argv the argument vector from \c main
* \param first_arg the index in \a argv of the first algorithm name
* \pre \a first_arg is from 1 to \a argc
* \return the names in the order given, or every key of \c map_func_name_to_func_info in
* sorted order
*/
[[nodiscard]] inline std::vector<std::string>
get_func_names(const int argc, char* const* argv, const int first_arg)
{
#if defined(DEBUG)
    assert((first_arg >= 1) && (first_arg <= argc));
#endif

    std::vector<std::string> func_names;

    if (argc > first_arg)
    {
        // use given functions
        func_names.reserve(static_cast<std::size_t>(argc - first_arg));
        for (int i = first_arg; i < argc; ++i)
        {
            func_names.emplace_back(argv[i]);
        }

        // validate func_names
        for (const auto& func_name : func_names)
        {
            if (!map_func_name_to_func_info.contains(func_name))
            {
                warnx("\"%s\" is not a valid function name", func_name.c_str());

                fmt::println(stderr, "Valid function names are:");
                const auto keys = std::views::keys(map_func_name_to_func_info);
                fmt::println(stderr, "  {}", fmt::join(keys, "\n  "));

                std::exit(EXIT_FAILURE);
            }
        }

        // A repeated name would run its algorithm twice, and the accuracy test would
        // overwrite or pool the results of the two runs under one name.
        for (auto it = func_names.cbegin(); it != func_names.cend(); ++it)
        {
            if (std::find(func_names.cbegin(), it, *it) != it)
            {
                errx(EXIT_FAILURE, "\"%s\" is given more than once", it->c_str());
            }
        }
    }
    else
    {
        // use all functions
        func_names.reserve(map_func_name_to_func_info.size());
        for (const auto& [func_name, ignore] : map_func_name_to_func_info)
        {
            func_names.push_back(func_name);
        }
    }

    return func_names;
}
