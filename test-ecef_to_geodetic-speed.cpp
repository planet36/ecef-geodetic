// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

// test speed of ECEF-to-Geodetic functions

#include "ecef-coord.hpp"
#include "ecef_to_geodetic-funcs.hpp"
#include "get_func_names.hpp"
#include "get_num_threads.hpp"
#include "map_func_name_to_func_info.hpp"
#include "read_coords.hpp"

#include <benchmark/benchmark.h>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <exception>
#include <string>
#include <vector>

/// time \a func on the points in \a ecef_vec
/**
* Each benchmark iteration converts the next point and wraps back to the first after the
* last, and the results are kept from being optimized away.
* \param BM_state the benchmark state that drives the iterations
* \param func the algorithm under test
* \param ecef_vec the input points
* \pre \a ecef_vec is not empty
*/
template <std::floating_point T>
void
BM_do_ecef_to_geodetic_test_speed(benchmark::State& BM_state,
                                  const ecef_to_geodetic_func<T>& func,
                                  const std::vector<ECEF<T>>& ecef_vec)
{
    T lat_rad{};
    T lon_rad{};
    T ht{};

    std::size_t i = 0;
    const std::size_t n = ecef_vec.size();

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

    const std::vector<std::string> func_names = get_func_names(argc, argv, 1);

    std::vector<ECEF<double>> ecef_vec;
    read_coords_ecef(ecef_vec);

    // With no points, the benchmark loop would read past the end of ecef_vec.
    if (ecef_vec.empty())
        errx(EXIT_FAILURE, "no input coordinates");

    for (const auto& func_name : func_names)
    {
        const auto& func_info = map_func_name_to_func_info.at(func_name);

        // Passing ecef_vec as an extra argument would copy it into every benchmark.  The
        // reference is safe because ecef_vec outlives RunSpecifiedBenchmarks.
        benchmark::RegisterBenchmark(func_info.display_name,
            [func = func_info.func, &ecef_vec](benchmark::State& BM_state)
            {
                BM_do_ecef_to_geodetic_test_speed(BM_state, func, ecef_vec);
            })->Threads(num_threads);
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
