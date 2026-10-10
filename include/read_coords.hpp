// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// Read ECEF and Geodetic coordinates from stdin
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "angle.hpp"
#include "ecef-coord.hpp"
#include "geodetic-coord.hpp"
#include "geodetic_to_ecef.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/// parse the numbers on one input line
/**
* \param line the input line
* \param line_num the 1-based line number, for error messages
* \return the numbers on the line
* \exception std::invalid_argument the line holds text that does not parse as a number
*/
template <std::floating_point T>
std::vector<T>
read_line_values(const std::string& line, const std::size_t line_num)
{
    std::istringstream values_stream(line);
    std::vector<T> values{std::istream_iterator<T>{values_stream}, std::istream_iterator<T>{}};

    // istream_iterator stops quietly at the first token it cannot parse, such as "nan", "junk",
    // or an overflow, so count the tokens to catch whatever it left behind.
    std::istringstream tokens_stream(line);
    const auto num_tokens = std::distance(std::istream_iterator<std::string>{tokens_stream},
                                          std::istream_iterator<std::string>{});

    if (std::cmp_not_equal(num_tokens, values.size()))
        throw std::invalid_argument(
            std::format("Invalid number on input line {}: {}", line_num, line));

    return values;
}

/// read ECEF coordinates from stdin
/**
* Each line holds X, Y, and Z, or W and Z, in meters.  W is the distance from the Z axis, and a
* line with W and Z puts the point on the prime meridian.  Blank lines are skipped.
* \param[in,out] ecef_vec the vector the points are appended to
* \exception std::invalid_argument a line holds text that does not parse as a number, or
* other than 2 or 3 values
*/
template <std::floating_point T>
void
read_coords_ecef(std::vector<ECEF<T>>& ecef_vec)
{
    std::string input_line;
    std::size_t line_num = 0;
    while (std::getline(std::cin, input_line))
    {
        ++line_num;
        const auto input_vec = read_line_values<T>(input_line, line_num);

        // The parser rejects any line whose tokens are not all numbers, so an empty result
        // means the line is blank.
        if (input_vec.empty())
            continue;

        T x{};
        T y{};
        T z{};

        switch (input_vec.size())
        {
        case 2: // W, Z
            x = std::abs(input_vec.at(0));
            y = 0;
            z = input_vec.at(1);
            break;

        case 3: // X, Y, Z
            x = input_vec.at(0);
            y = input_vec.at(1);
            z = input_vec.at(2);
            break;

        default:
            throw std::invalid_argument(std::format(
                "Invalid input data dimensions on line {}: {}", line_num, input_vec.size()));
            break;
        }

        ecef_vec.emplace_back(x, y, z);
    }
}

/// read Geodetic coordinates from stdin
/**
* Each line holds the latitude, longitude, and height, or the latitude and height.  Angles are
* in degrees and heights in meters, and a line without a longitude puts the point on the prime
* meridian.  Blank lines are skipped.
* \param[in,out] geod_vec the vector the points are appended to
* \exception std::invalid_argument a line holds text that does not parse as a number, or
* other than 2 or 3 values
*/
template <std::floating_point T>
void
read_coords_geod(std::vector<Geodetic<angle_unit::degree, T>>& geod_vec)
{
    std::string input_line;
    std::size_t line_num = 0;
    while (std::getline(std::cin, input_line))
    {
        ++line_num;
        const auto input_vec = read_line_values<T>(input_line, line_num);

        // The parser rejects any line whose tokens are not all numbers, so an empty result
        // means the line is blank.
        if (input_vec.empty())
            continue;

        // input angle unit is degrees
        ang_deg<T> lat;
        ang_deg<T> lon;
        T ht{};

        switch (input_vec.size())
        {
        case 2: // lat, ht
            lat = input_vec.at(0);
            lon = 0;
            ht = input_vec.at(1);
            break;

        case 3: // lat, lon, ht
            lat = input_vec.at(0);
            lon = input_vec.at(1);
            ht = input_vec.at(2);
            break;

        default:
            throw std::invalid_argument(std::format(
                "Invalid input data dimensions on line {}: {}", line_num, input_vec.size()));
            break;
        }

        geod_vec.emplace_back(lat, lon, ht);
    }
}

/// the coordinate system of the input lines
enum struct INPUT_DATA_COORD_SYSTEM
{
    ECEF,
    GEODETIC,
};

/// get the name of the input coordinate system
[[nodiscard]] inline std::string_view
to_string(const INPUT_DATA_COORD_SYSTEM x) noexcept
{
    switch (x)
    {
    case INPUT_DATA_COORD_SYSTEM::ECEF:
        return "ECEF";
        break;
    case INPUT_DATA_COORD_SYSTEM::GEODETIC:
        return "GEODETIC";
        break;
    default:
        return "???";
        break;
    }
}

/// read ECEF or Geodetic coordinates from stdin
/**
* Geodetic input is converted to ECEF.  \c read_coords_ecef and \c read_coords_geod describe
* the line formats.
* \param input_data_coord_system the coordinate system of the input lines
* \param[in,out] ecef_vec the vector the points are appended to
* \exception std::invalid_argument a line holds text that does not parse as a number, or
* other than 2 or 3 values, or \a input_data_coord_system is not a known value
*/
template <std::floating_point T>
void
read_coords(const INPUT_DATA_COORD_SYSTEM input_data_coord_system,
            std::vector<ECEF<T>>& ecef_vec)
{
    switch (input_data_coord_system)
    {
    case INPUT_DATA_COORD_SYSTEM::ECEF:
        read_coords_ecef(ecef_vec);
        break;

    case INPUT_DATA_COORD_SYSTEM::GEODETIC:
        {
            std::vector<Geodetic<angle_unit::degree, T>> geod_vec;
            read_coords_geod(geod_vec);

            ecef_vec.reserve(ecef_vec.size() + geod_vec.size());

            for (const auto& geod : geod_vec)
            {
                ecef_vec.push_back(geodetic_to_ecef(geod));
            }
        }

        break;

    default:
        throw std::invalid_argument(std::format("Invalid INPUT_DATA_COORD_SYSTEM: {}",
                                                std::to_underlying(input_data_coord_system)));
        break;
    }
}
