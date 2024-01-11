// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

// test speed of ECEF-to-Geodetic functions

#include "map_func_name_to_func_info.hpp"
#include "read_coords.hpp"

#include <algorithm>
#include <benchmark/benchmark.h>
#include <concepts>
#include <cstdlib>
#include <fmt/ranges.h>
#include <ranges>
#include <string>
#include <thread>
#include <vector>

template <std::floating_point T>
void BM_do_ecef_to_geodetic_test_speed(
	benchmark::State& state,
	const ecef_to_geodetic_func<T>& func,
	const std::vector<ECEF<T>>& ecef_vec)
{
	for (auto _ : state)
	{
		for (const auto& ecef_given : ecef_vec)
		{
			T lat_rad{};
			T lon_rad{};
			T ht{};
			func(ecef_given.x, ecef_given.y, ecef_given.z, lat_rad, lon_rad, ht);
			benchmark::DoNotOptimize(lat_rad);
			benchmark::DoNotOptimize(lon_rad);
			benchmark::DoNotOptimize(ht);
		}
	}

	state.SetItemsProcessed(state.iterations() * ecef_vec.size());
}

int main(int argc, char** argv)
{
	// copied from /usr/include/benchmark/benchmark.h
	benchmark::Initialize(&argc, argv);

	// function names may be given in argv
	//if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;

	// {{{ determine num_threads

	constexpr unsigned int min_threads = 1;
	const unsigned int max_threads = (std::thread::hardware_concurrency() != 0) ? std::thread::hardware_concurrency() : min_threads;
	// https://en.wikipedia.org/wiki/Elvis_operator
	//const unsigned int max_threads = std::thread::hardware_concurrency() ?: min_threads;

	unsigned int num_threads;

	try
	{
		num_threads = static_cast<unsigned int>(std::stoi(std::getenv("NUM_THREADS")));
	}
	catch (...)
	{
		num_threads = max_threads;
	}

	num_threads = std::clamp(num_threads, min_threads, max_threads);

	/*
	if (num_threads > min_threads)
		// Don't use all the cores
		--num_threads;
	*/

	// }}}

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
				fmt::println(stderr, "Error: \"{}\" is not a valid function name.", func_name);

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

		if (num_threads == 1)
			benchmark::RegisterBenchmark(func_name,
					&BM_do_ecef_to_geodetic_test_speed<double>,
					func_info.func, ecef_vec);
		else
			benchmark::RegisterBenchmark(func_name,
					&BM_do_ecef_to_geodetic_test_speed<double>,
					func_info.func, ecef_vec)->Threads(num_threads);
	}

	benchmark::RunSpecifiedBenchmarks();

	benchmark::Shutdown();

	return 0;
}
