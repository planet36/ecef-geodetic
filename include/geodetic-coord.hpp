// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Geodetic coordinate class
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "angle.hpp"

#include <cmath>
#include <concepts>
#include <fmt/base.h>
#include <sstream>
#include <string>

constexpr int geodetic_default_precision = 6; ///< The default precision
constexpr int geodetic_precision_add = 5; ///< The precision to add to the geodetic latitude and longitude

/// Geodetic to string
/**
* \param lat geodetic latitude
* \param lon geodetic longitude
* \param ht ellipsoid height (meters)
* \param precision the number of digits after the decimal place to generate for \a ht
* \note \a lat and \a lon are converted to decimal degrees.
* \note The precision of \a lat and \a lon is <code>\a precision + geodetic_precision_add</code>.
* \return the string representation of the geodetic coordinate
*/
template <angle_unit U, std::floating_point T>
[[nodiscard]] std::string
geodetic_to_str(const angle<U, T>& lat,
                const angle<U, T>& lon,
                const T ht,
                int precision = geodetic_default_precision)
{
    std::string result;

    result += fmt::format("{:.{}f}", lat.to_deg(), precision + geodetic_precision_add);
    result += fmt::format(" {:.{}f}", lon.to_deg(), precision + geodetic_precision_add);

    if (!std::isnan(ht))
        result += fmt::format(" {:.{}f}", ht, precision);

    return result;
}

/// string to Geodetic
/**
* \param s the string representation of the geodetic coordinate
* \note \a lat and \a lon are converted from decimal degrees.
* \param[out] lat geodetic latitude
* \param[out] lon geodetic longitude
* \param[out] ht ellipsoid height (meters)
* \retval true if an error has occurred on the associated stream
*/
template <angle_unit U, std::floating_point T>
bool
str_to_geodetic(const std::string& s, angle<U, T>& lat, angle<U, T>& lon, T& ht)
{
    std::istringstream iss(s);

    T lat_deg{};
    T lon_deg{};

    iss >> lat_deg >> lon_deg;

    lat = ang_deg<T>{lat_deg};
    lon = ang_deg<T>{lon_deg};

    iss >> std::ws;

    //ht = 0;

    // ht is optional
    if (!(iss.eof() || iss.bad()))
        iss >> ht;

    // https://en.cppreference.com/w/cpp/io/basic_ios/operator!
    return !iss;
}

/// Geodetic coordinate
/**
* \sa https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
*/
template <angle_unit U, std::floating_point T>
struct Geodetic
{
    using this_t = Geodetic<U, T>;

    angle<U, T> lat{}; // geodetic latitude
    angle<U, T> lon{}; // geodetic longitude
    T ht{};            // ellipsoid/geodetic height (meters)

    Geodetic() = default;

    constexpr Geodetic(
        const angle<U, T>& _lat,
        const angle<U, T>& _lon,
        const T _ht = 0) :
        lat(_lat),
        lon(_lon),
        ht(_ht)
    {}

    /// conversion ctor
    template <std::floating_point T2>
    constexpr Geodetic(const Geodetic<U, T2>& that) :
        lat(that.lat),
        lon(that.lon),
        ht(that.ht)
    {}

    auto operator<=>(const this_t&) const = default;

    void normalize() { normalize_geodetic(lat, lon); }

    [[nodiscard]] std::string to_string(const int precision = geodetic_default_precision) const
    {
        return geodetic_to_str(lat, lon, ht, precision);
    }
};
