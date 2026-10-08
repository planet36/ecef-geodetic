// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// Get the number of threads to run each benchmark on
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "parse_int.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <exception>
#include <thread>

/// Get the number of threads to run each benchmark on
/**
* The count is read from the environment variable \c NUM_THREADS, where 0 means one thread
* per hardware thread.  When \c NUM_THREADS is not set, the count is 1.
*
* The program exits if \c NUM_THREADS is not an integer in the interval from 0 to the number
* of hardware threads.
*/
[[nodiscard]] inline int
get_num_threads()
{
    constexpr int min_threads = 1;
    const auto hw_threads = static_cast<int>(std::thread::hardware_concurrency());
    const auto max_threads = std::max(min_threads, hw_threads);

    int num_threads = min_threads;

    try
    {
        // parse_env_int takes these args: name, min, max, default
        num_threads = parse_env_int("NUM_THREADS", 0, max_threads, min_threads);
    }
    catch (const std::exception& ex)
    {
        (void)std::fflush(stdout);
        errx(EXIT_FAILURE, "%s", ex.what());
    }

    if (num_threads == 0)
    {
        num_threads = max_threads;
    }

    return num_threads;
}
