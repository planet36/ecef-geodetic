// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Read ECEF and Geodetic coordinates from stdin
/**
\file
\author Steven Ward
*/

#pragma once

#include "angle.hpp"
#include "ecef-coord.hpp"
#include "geodetic-coord.hpp"
#include "geodetic_to_ecef.hpp"

#include <cmath>
#include <concepts>
#include <fmt/format.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/// read ECEF coordinates from stdin
template <std::floating_point T>
void
read_coords_ecef(std::vector<ECEF<T>>& ecef_vec)
{
    std::string input_line;
    while (std::getline(std::cin, input_line))
    {
        std::istringstream iss(input_line);

        const std::vector<T> input_vec{std::istream_iterator<T>{iss},
                                       std::istream_iterator<T>{}};

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
            throw std::invalid_argument(fmt::format(
                "Invalid input data dimensions: {}", input_vec.size()));
            break;
        }

        ecef_vec.emplace_back(x, y, z);
    }
}

/// read Geodetic coordinates from stdin
template <std::floating_point T>
void
read_coords_geod(std::vector<Geodetic<angle_unit::degree, T>>& geod_vec)
{
    std::string input_line;
    while (std::getline(std::cin, input_line))
    {
        std::istringstream iss(input_line);

        const std::vector<T> input_vec{std::istream_iterator<T>{iss},
                                       std::istream_iterator<T>{}};

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
            throw std::invalid_argument(fmt::format(
                "Invalid input data dimensions: {}", input_vec.size()));
            break;
        }

        geod_vec.emplace_back(lat, lon, ht);
    }
}

enum struct INPUT_DATA_COORD_SYSTEM
{
    ECEF,
    GEODETIC,
};

[[nodiscard]] std::string_view
to_string(const INPUT_DATA_COORD_SYSTEM& x)
{
    switch (x)
    {
    case INPUT_DATA_COORD_SYSTEM::ECEF    : return "ECEF"    ; break;
    case INPUT_DATA_COORD_SYSTEM::GEODETIC: return "GEODETIC"; break;
    default: return "???"; break;
    }
}

/// read ECEF or Geodetic coordinates from stdin
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

            ecef_vec.reserve(geod_vec.size());

            for (const auto& geod : geod_vec)
            {
                ecef_vec.push_back(geodetic_to_ecef(geod));
            }
        }

        break;

    default:
        throw std::invalid_argument(
            fmt::format("Invalid INPUT_DATA_COORD_SYSTEM: {}",
                        std::to_underlying(input_data_coord_system)));
        break;
    }
}
