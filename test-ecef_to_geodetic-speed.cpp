// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

// test speed of ECEF-to-Geodetic functions

#include "get_num_threads.hpp"
#include "map_func_name_to_func_info.hpp"
#include "read_coords.hpp"

#include <benchmark/benchmark.h>
#include <concepts>
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <exception>
#include <fmt/ranges.h>
#include <ranges>
#include <string>
#include <vector>

template <std::floating_point T>
void
BM_do_ecef_to_geodetic_test_speed(benchmark::State& BM_state,
                                  const ecef_to_geodetic_func<T>& func,
                                  const std::vector<ECEF<T>>& ecef_vec)
{
    T lat_rad{};
    T lon_rad{};
    T ht{};

    size_t i = 0;
    const size_t n = ecef_vec.size();

    for (auto _ : BM_state) // NOLINT(clang-analyzer-deadcode.DeadStores)
    {
        func(ecef_vec[i].x, ecef_vec[i].y, ecef_vec[i].z, lat_rad, lon_rad, ht);

        benchmark::DoNotOptimize(lat_rad);
        benchmark::DoNotOptimize(lon_rad);
        benchmark::DoNotOptimize(ht);

        ++i;
        if (i == n)
            i = 0;
    }
}

int
main(int argc, char* argv[])
try
{
    // Copied from benchmark.h
    benchmark::MaybeReenterWithoutASLR(argc, argv);
    benchmark::Initialize(&argc, argv);

    /*
    // function names may be given in argv
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
    */

    const int num_threads = get_num_threads();

    std::vector<std::string> func_names;

    if (argc > 1)
    {
        // use given functions
        for (int i = 1; i < argc; ++i)
        {
            func_names.emplace_back(argv[i]);
        }

        // validate func_names
        for (const auto& func_name : func_names)
        {
            // verify the given function names are valid
            if (!map_func_name_to_func_info.contains(func_name))
            {
                fmt::println(stderr, "Error: \"{}\" is not a valid function name.",
                             func_name);

                fmt::println(stderr, "Valid function names are:");
                const auto keys = std::views::keys(map_func_name_to_func_info);
                fmt::println(stderr, "  {}", fmt::join(keys, "\n  "));

                return EXIT_FAILURE;
            }
        }
    }
    else
    {
        // use all functions
        for (const auto& [func_name, ignore] : map_func_name_to_func_info)
        {
            func_names.push_back(func_name);
        }
    }

    std::vector<ECEF<double>> ecef_vec;
    read_coords_ecef(ecef_vec);

    for (const auto& func_name : func_names)
    {
        const auto& func_info = map_func_name_to_func_info.at(func_name);

        benchmark::RegisterBenchmark(func_info.display_name,
                                     &BM_do_ecef_to_geodetic_test_speed<double>,
                                     func_info.func, ecef_vec)->Threads(num_threads);
    }

    benchmark::RunSpecifiedBenchmarks();

    benchmark::Shutdown();

    return 0;
}
catch (const std::exception& ex)
{
    (void)std::fflush(stdout);
    errx(EXIT_FAILURE, "%s", ex.what());
}
