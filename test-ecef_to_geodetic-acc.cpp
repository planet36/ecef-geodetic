// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

// test accuracy (and optionally speed) of ECEF-to-Geodetic functions

#include "angle.hpp"
#include "angle_unit.hpp"
#include "ecef-coord.hpp"
#include "ecef_to_geodetic-funcs.hpp"
#include "geodetic-coord.hpp"
#include "geodetic_to_ecef.hpp"
#include "ilog.hpp"
#include "map_func_name_to_func_info.hpp"
#include "read_coords.hpp"
#include "running_stats.hpp"
#include "stats.hpp"

#include <algorithm>
#include <benchmark/benchmark.h>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <execution>
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <gnu/libc-version.h>
#include <map>
#include <mutex>
#include <nlohmann/json.hpp>
#include <numeric>
#include <random>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unistd.h>
#include <utility>
#include <vector>

inline constexpr std::string_view program_version = "2024-01-09";
// https://man7.org/linux/man-pages/man3/gnu_get_libc_version.3.html
const std::string_view glibc_version = gnu_get_libc_version();
// https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
inline constexpr std::string_view gcc_version = __VERSION__;

/// the distance error (m) recorded in place of one that is not finite
/**
* It keeps the statistics finite, and it is large enough to make the algorithm's
* \c ilog10_mean_dist_err 99.
*/
inline constexpr double non_finite_dist_err = 99E97;

/// the \c ilog10_mean_dist_err reported when the mean distance error is exactly 0
inline constexpr int ilog10_mean_dist_err_exact = -99;

/// the \c ilog10_mean_dist_err reported for an inaccurate algorithm
inline constexpr int ilog10_mean_dist_err_inaccurate = 99;

/// the largest \c ilog10_mean_dist_err of an accurate algorithm
inline constexpr int max_ilog10_mean_dist_err_accurate = 2;

/// statistics of the distance errors (m) of one algorithm
template <std::floating_point T>
struct dist_err_stats
{
    T mean{};
    T stdev{};
    T max{};
    T sum{};

    dist_err_stats() = default;

    /// copy the statistics from \a rs
    explicit dist_err_stats(const running_stats<T>& rs) noexcept :
    mean(rs.mean()),
    stdev(rs.standard_deviation()),
    max(rs.max_abs()),
    sum(rs.sum_abs())
    {}
};

/// get the round-trip distance error (m) of \a func at \a ecef_given
/**
* The geodetic result of \a func is converted back to ECEF exactly.
*/
template <std::floating_point T>
auto
round_trip_dist_err(const ecef_to_geodetic_func<T>& func, const ECEF<T>& ecef_given) noexcept
{
    T lat_rad{};
    T lon_rad{};
    T ht{};
    func(ecef_given.x, ecef_given.y, ecef_given.z, lat_rad, lon_rad, ht);
    const Geodetic<angle_unit::radian, T> geod_result{lat_rad, lon_rad, ht};
    const auto ecef_result = geodetic_to_ecef(geod_result);
    const auto dist_err = euclidean_dist(ecef_given, ecef_result);

    return dist_err;
}

/// get the distance error statistics of \a func over \a ecef_vec from running statistics
template <std::floating_point T>
auto
do_ecef_to_geodetic_test_acc_running(const ecef_to_geodetic_func<T>& func,
                                     const std::vector<ECEF<T>>& ecef_vec) noexcept
{
    running_stats<T> rs;

    for (const auto& ecef_given : ecef_vec)
    {
        auto dist_err = round_trip_dist_err(func, ecef_given);
        if (!std::isfinite(dist_err))
        {
            dist_err = non_finite_dist_err;
        }

        rs.push(dist_err);
    }

    return dist_err_stats(rs);
}

/// get the distance error statistics of \a func over \a ecef_vec from every error, sorted
/**
* Summing the errors in ascending order tries to minimize floating-point precision loss.
* However, it reduces the precision loss of the sum by less than 1E-12, and it takes about
* twice as long as \c do_ecef_to_geodetic_test_acc_running.
*/
template <std::floating_point T>
auto
do_ecef_to_geodetic_test_acc_collect(const ecef_to_geodetic_func<T>& func,
                                     const std::vector<ECEF<T>>& ecef_vec)
{
    std::vector<T> dist_errs;
    dist_errs.reserve(ecef_vec.size());
    dist_err_stats<T> stats;

    for (const auto& ecef_given : ecef_vec)
    {
        auto dist_err = round_trip_dist_err(func, ecef_given);
        if (!std::isfinite(dist_err))
        {
            dist_err = non_finite_dist_err;
        }
        dist_errs.push_back(dist_err);
    }

    std::ranges::sort(dist_errs);

    stats.mean = arithmetic_mean_val(dist_errs);
    // the sample standard deviation, which is what running_stats gives
    stats.stdev = stdev_val(dist_errs, true);
    stats.max = max_val(dist_errs);
    stats.sum = sum_val(dist_errs);

    return stats;
}

/// get the distance error statistics of \a func over \a ecef_vec
/**
* \param func the algorithm under test
* \param ecef_vec the input points
* \param collect_dist_err whether to keep every distance error for more precise statistics
*/
template <std::floating_point T>
inline auto
do_ecef_to_geodetic_test_acc(const ecef_to_geodetic_func<T>& func,
                             const std::vector<ECEF<T>>& ecef_vec,
                             const bool collect_dist_err)
{
    if (collect_dist_err)
        return do_ecef_to_geodetic_test_acc_collect(func, ecef_vec);
    else
        return do_ecef_to_geodetic_test_acc_running(func, ecef_vec);
}

#if 1
/// run \a func on every point in \a ecef_vec, keeping the compiler from discarding the results
template <std::floating_point T>
void
do_ecef_to_geodetic_test_speed(const ecef_to_geodetic_func<T>& func,
                               const std::vector<ECEF<T>>& ecef_vec) noexcept
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
#else
auto do_ecef_to_geodetic_test_speed =
    []<std::floating_point T>(const ecef_to_geodetic_func<T>& func,
                              const std::vector<ECEF<T>>& ecef_vec)
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
};
#endif

#define nl (void)putchar('\n')

/// print the command-line usage to stderr
void
print_usage(const char* program_name)
{
    fmt::print(stderr, R"(Usage: {} [OPTION]... [FUNC_NAME]...
Test the accuracy of ECEF-to-Geodetic functions, and print the results as JSON.
Read the input coordinates from stdin.  Test every function unless FUNC_NAMEs are given.

  -a    run the accuracy test
  -1    run the single-point accuracy test
  -s N  run N rounds of the speed test
  -m N  skip functions whose ilog10_mean_dist_err exceeds N
  -g    read geodetic input instead of ECEF
  -t    use multiple threads
  -c    collect every distance error for more precise statistics
  -v    print progress to stderr
)",
               program_name);
}

int
main(int argc, char* argv[])
{
    /*
    ** Methodology for testing a single point:
    ** Given a Geodetic point (g1), perform an exact conversion to ECEF (e1).
    ** For the given algorithm, convert e1 to Geodetic (g2).
    ** Perform an exact conversion of g2 to ECEF (e2).
    ** Calculate the Euclidean distance from e1 to e2.  This is the distance error.
    **
    ** For many Geodetic points, vary the latitude and height.
    ** The longitude calculation is the same in all algorithms (lon = atan2(y, x)), so it's not
    ** meaningful to vary.
    ** The set of test points must include some at the equator, the poles, and non-zero heights.
    */

    using namespace std::literals;
    using json = nlohmann::json;

    constexpr INPUT_DATA_COORD_SYSTEM default_input_data_coord_system =
        INPUT_DATA_COORD_SYSTEM::ECEF;

    bool verbose = false;
    bool do_acc_test = false;
    bool do_single_point_acc_test = false;
    unsigned int num_speed_test_iterations = 0;
    int max_ilog10_mean_dist_err = ilog10_mean_dist_err_inaccurate;
    INPUT_DATA_COORD_SYSTEM input_data_coord_system = default_input_data_coord_system;
    bool use_multiple_threads = false;
    bool collect_dist_err = false;
    json json_output;

    constexpr const char* short_options = "+va1s:m:gtc";

    int c = 0;
    while ((c = getopt(argc, argv, short_options)) != -1)
    {
        switch (c)
        {
        case 'v':
            verbose = true;
            break;

        case 'a':
            do_acc_test = true;
            break;

        case '1':
            do_single_point_acc_test = true;
            break;

        case 's':
            try
            {
                const auto tmp = std::stoll(optarg);
                num_speed_test_iterations =
                    std::saturating_cast<decltype(num_speed_test_iterations)>(tmp);
            }
            catch (const std::invalid_argument& ex)
            {
                errx(EXIT_FAILURE, "invalid argument: %s: \"%s\"", ex.what(), optarg);
            }
            catch (const std::out_of_range& ex)
            {
                errx(EXIT_FAILURE, "out of range: %s: \"%s\"", ex.what(), optarg);
            }
            break;

        case 'm':
            try
            {
                max_ilog10_mean_dist_err = std::stoi(optarg);
            }
            catch (const std::invalid_argument& ex)
            {
                errx(EXIT_FAILURE, "invalid argument: %s: \"%s\"", ex.what(), optarg);
            }
            catch (const std::out_of_range& ex)
            {
                errx(EXIT_FAILURE, "out of range: %s: \"%s\"", ex.what(), optarg);
            }
            break;

        case 'g':
            input_data_coord_system = INPUT_DATA_COORD_SYSTEM::GEODETIC;
            break;

        case 't':
            use_multiple_threads = true;
            break;

        case 'c':
            collect_dist_err = true;
            break;

        default:
            print_usage(argv[0]);
            std::exit(EXIT_FAILURE);
        }
    }

    std::vector<std::string> func_names;

    if (argc > optind)
    {
        // use given functions
        for (int i = optind; i < argc; ++i)
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

                std::exit(EXIT_FAILURE);
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

    // filter out the functions whose hard-coded ilog10_mean_dist_err exceeds
    // max_ilog10_mean_dist_err
    erase_if(func_names,
             [max_ilog10_mean_dist_err](const auto& func_name)
             {
                 const auto it = map_func_name_to_func_info.find(func_name);
                 return it->second.ilog10_mean_dist_err > max_ilog10_mean_dist_err;
             });

    for (const auto& func_name : func_names)
    {
        const auto& func_info = map_func_name_to_func_info.at(func_name);

        json_output["func_names"][func_name]["info"] = {
            {"num_lines", func_info.num_lines},
            {"needs_code_for_corner_cases", func_info.needs_code_for_corner_cases},
            {"ilog10_mean_dist_err", func_info.ilog10_mean_dist_err},
            {"display_name", func_info.display_name},
            {"algo_author", func_info.algo_author},
            {"code_copyright", func_info.code_copyright},
            {"license", func_info.license},
            {"orig_impl_lang", func_info.orig_impl_lang},
            {"url", func_info.url},
            {"citation", func_info.citation},
        };
    }

    json_output["program_version"] = program_version;
    json_output["glibc_version"] = glibc_version;
    json_output["gcc_version"] = gcc_version;
    json_output["num_speed_test_iterations"] = num_speed_test_iterations;
    json_output["max_ilog10_mean_dist_err"] = max_ilog10_mean_dist_err;
    json_output["input_data_coord_system"] = to_string(input_data_coord_system);
    json_output["units_of_measurement"]["acc"] = "meters";
    json_output["units_of_measurement"]["speed"] = "nanoseconds";
    // https://en.cppreference.com/w/cpp/numeric/math/math_errhandling
    json_output["math_errno_set"] = static_cast<bool>(math_errhandling & MATH_ERRNO);
    json_output["math_errexcept_set"] = static_cast<bool>(math_errhandling & MATH_ERREXCEPT);

    std::vector<ECEF<double>> ecef_vec;

    if (do_acc_test || do_single_point_acc_test || (num_speed_test_iterations > 0))
    {
        if (verbose)
            fmt::print(stderr, "# reading input data ... ");
        read_coords(input_data_coord_system, ecef_vec);
        if (verbose)
            fmt::println(stderr, "done");

        // With no points, every mean distance error would be NaN.
        if (ecef_vec.empty())
            errx(EXIT_FAILURE, "no input coordinates");
    }

    json_output["num_input_coords"] = ecef_vec.size();

    if (do_acc_test)
    {
        if (verbose)
            fmt::print(stderr, "# doing accuracy test ... ");

        if (use_multiple_threads)
        {
            std::mutex mtx;
            std::for_each(std::execution::par, std::begin(func_names), std::end(func_names),
                [&](const std::string& func_name)
                {
                    const auto& func_info = map_func_name_to_func_info.at(func_name);

                    const auto stats = do_ecef_to_geodetic_test_acc(func_info.func, ecef_vec,
                                                                    collect_dist_err);

                    int ilog10_mean_dist_err = ilog10(stats.mean);
                    if (stats.mean == 0)
                        ilog10_mean_dist_err = ilog10_mean_dist_err_exact;
                    else if (ilog10_mean_dist_err > max_ilog10_mean_dist_err_accurate)
                        ilog10_mean_dist_err = ilog10_mean_dist_err_inaccurate;

                    std::lock_guard guard{mtx};

                    json_output["func_names"][func_name]["acc"] = {
                        {"mean_dist_err", stats.mean},
                        {"stdev_dist_err", stats.stdev},
                        {"max_dist_err", stats.max},
                        {"sum_dist_err", stats.sum},
                        {"ilog10_mean_dist_err", ilog10_mean_dist_err},
                    };
                });
        }
        else
        {
            for (const auto& func_name : func_names)
            {
                if (verbose)
                    fmt::println(stderr, "# {}", func_name);

                const auto& func_info = map_func_name_to_func_info.at(func_name);

                const auto stats =
                    do_ecef_to_geodetic_test_acc(func_info.func, ecef_vec, collect_dist_err);

                int ilog10_mean_dist_err = ilog10(stats.mean);
                if (stats.mean == 0)
                    ilog10_mean_dist_err = ilog10_mean_dist_err_exact;
                else if (ilog10_mean_dist_err > max_ilog10_mean_dist_err_accurate)
                    ilog10_mean_dist_err = ilog10_mean_dist_err_inaccurate;

                json_output["func_names"][func_name]["acc"] = {
                    {"mean_dist_err", stats.mean},
                    {"stdev_dist_err", stats.stdev},
                    {"max_dist_err", stats.max},
                    {"sum_dist_err", stats.sum},
                    {"ilog10_mean_dist_err", ilog10_mean_dist_err},
                };
            }
        }

        if (verbose)
            fmt::println(stderr, "done");
    }

    if (do_single_point_acc_test)
    {
        // Every function's record for a point starts with the same string and radius.
        std::vector<json> point_records;
        point_records.reserve(ecef_vec.size());
        for (const auto& ecef_given : ecef_vec)
        {
            point_records.push_back(
                {{"ecef", ecef_given.to_string()}, {"radius", L2_norm(ecef_given)}});
        }

        for (const auto& func_name : func_names)
        {
            const auto& func_info = map_func_name_to_func_info.at(func_name);
            auto& acc1 = json_output["func_names"][func_name]["acc1"];

            // Pair each point with its record in point_records.
            for (const auto& [ecef_given, point_record] :
                 std::views::zip(ecef_vec, point_records))
            {
                // Copy the record, since every function starts from the same one.
                json record = point_record;
                record["dist_err"] = round_trip_dist_err(func_info.func, ecef_given);
                // Move the copy into the array, since it is not used again.
                acc1.push_back(std::move(record));
            }
        }
    }

    if (num_speed_test_iterations > 0)
    {
        // This cannot be a map of func name to running stats because the median
        // value must be obtained later.
        std::map<std::string, std::multiset<double>> map_func_name_to_time_per_call;

        std::random_device rd;
        std::seed_seq seeder{rd(), rd(), rd(), rd()};
        std::mt19937_64 rng(seeder);

        // Used for padding the output string
        const std::size_t max_strlen_num_speed_test_iterations =
            std::to_string(num_speed_test_iterations).size();

        // Used to estimate time remaining (minutes)
        running_stats<double> speed_test_iteration_durations;

        while (num_speed_test_iterations > 0)
        {
            const auto iteration_t0 = std::chrono::steady_clock::now();
            // The system clock can be converted to tm, not the steady_clock.
            const auto system_clock_now = std::chrono::system_clock::now();

            if (verbose)
                fmt::println(stderr, "# speed tests remaining: {:{}};  {:%FT%T%z}",
                             num_speed_test_iterations, max_strlen_num_speed_test_iterations,
                             system_clock_now);

            if (speed_test_iteration_durations.num_data_values() > 0)
            {
                if (verbose)
                    fmt::println(stderr, "# est. time remaining: {:.1f} min",
                                 speed_test_iteration_durations.mean() *
                                     num_speed_test_iterations);
            }

            std::shuffle(func_names.begin(), func_names.end(), rng);

            // With multiple threads, each algorithm is about 20% slower on average, presumably
            // because of cache misses.
            if (use_multiple_threads)
            {
                std::mutex mtx;
                std::for_each(
                    std::execution::par, std::begin(func_names), std::end(func_names),
                    [&](const std::string& func_name)
                    {
                        const auto& func_info = map_func_name_to_func_info.at(func_name);

                        const auto t0 = std::chrono::steady_clock::now();

                        do_ecef_to_geodetic_test_speed(func_info.func, ecef_vec);

                        const auto t1 = std::chrono::steady_clock::now();

                        const auto duration = t1 - t0;

                        // (nanoseconds)
                        const double time_per_call = (duration / 1.0ns) / ecef_vec.size();

                        std::lock_guard guard{mtx};

                        map_func_name_to_time_per_call[func_name].insert(time_per_call);
                    });
            }
            else
            {
                if (verbose)
                    fmt::print(stderr, "# ");

                for (const auto& func_name : func_names)
                {
                    if (verbose)
                        fmt::print(stderr, ".");

                    const auto& func_info = map_func_name_to_func_info.at(func_name);

                    const auto t0 = std::chrono::steady_clock::now();

                    do_ecef_to_geodetic_test_speed(func_info.func, ecef_vec);

                    const auto t1 = std::chrono::steady_clock::now();

                    const auto duration = t1 - t0;

                    // (nanoseconds)
                    const double time_per_call = (duration / 1.0ns) / ecef_vec.size();

                    map_func_name_to_time_per_call[func_name].insert(time_per_call);
                }

                if (verbose)
                    fmt::println(stderr, "");
            }

            const auto iteration_t1 = std::chrono::steady_clock::now();

            const auto iteration_duration = iteration_t1 - iteration_t0;

            speed_test_iteration_durations.push(iteration_duration / 1.0min);

            num_speed_test_iterations--;
        }

        if (!map_func_name_to_time_per_call.empty())
        {
            // func_names might have been shuffled
            for (const auto& [func_name, multiset_time_per_call] :
                 map_func_name_to_time_per_call)
            {
                running_stats<double> rs;

                // convert the multiset of speed data to running stats
                rs.push(multiset_time_per_call.cbegin(), multiset_time_per_call.cend());

                json_output["func_names"][func_name]["speed"] = {
                    {"median_time_per_call",
                     fmt::format("{:.1f}", median_val(multiset_time_per_call))},
                    {"mean_time_per_call", fmt::format("{:.1f}", rs.mean())},
                    {"stdev_time_per_call", fmt::format("{:.1f}", rs.standard_deviation())},
                };
            }
        }
    }

    fmt::println("{}", json_output.dump(4));

    return 0;
}
