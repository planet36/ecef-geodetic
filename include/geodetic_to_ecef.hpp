// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Geodetic-to-ECEF coordinate conversion
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "angle.hpp"
#include "ecef-coord.hpp"
#include "ellipsoid-wgs84.hpp"
#include "geodetic-coord.hpp"

#include <cmath>
#include <concepts>

/// convert from geodetic to ECEF
/**
* Source:
* NGA.STND.0036_1.0.0_WGS84 2014-07-08
* Equation (4-14)
* \param lat_rad geodetic latitude (radians)
* \param lon_rad geodetic longitude (radians)
* \param ht ellipsoid height (meters)
* \param[out] x X coordinate (meters)
* \param[out] y Y coordinate (meters)
* \param[out] z Z coordinate (meters)
*/
template <std::floating_point T>
void
geodetic_to_ecef(const T lat_rad, const T lon_rad, const T ht, T& x, T& y, T& z)
{
    static constexpr auto& ell = WGS84<T>;

    const auto sin_lat = std::sin(lat_rad);
    const auto cos_lat = std::cos(lat_rad);

    const auto sin_lon = std::sin(lon_rad);
    const auto cos_lon = std::cos(lon_rad);

    const auto Rn = ell.get_Rn(sin_lat);

    x = (Rn + ht) * cos_lat * cos_lon;
    y = (Rn + ht) * cos_lat * sin_lon;
    z = (Rn * (1 - ell.e2) + ht) * sin_lat;
}

/// convert from geodetic to ECEF
/**
* Source:
* NGA.STND.0036_1.0.0_WGS84 2014-07-08
* Equation (4-14)
* \param lat geodetic latitude
* \param lon geodetic longitude
* \param ht ellipsoid height (meters)
* \param[out] x X coordinate (meters)
* \param[out] y Y coordinate (meters)
* \param[out] z Z coordinate (meters)
*/
template <angle_unit U, std::floating_point T>
void
geodetic_to_ecef(const angle<U, T>& lat, const angle<U, T>& lon, const T ht, T& x, T& y, T& z)
{
    static constexpr auto& ell = WGS84<T>;

    const auto sin_lat = sin(lat);
    const auto cos_lat = cos(lat);

    const auto sin_lon = sin(lon);
    const auto cos_lon = cos(lon);

    const auto Rn = ell.get_Rn(sin_lat);

    x = (Rn + ht) * cos_lat * cos_lon;
    y = (Rn + ht) * cos_lat * sin_lon;
    z = (Rn * (1 - ell.e2) + ht) * sin_lat;
}

template <angle_unit U, std::floating_point T>
auto
geodetic_to_ecef(const Geodetic<U, T>& geod)
{
    T x{};
    T y{};
    T z{};
    geodetic_to_ecef(geod.lat, geod.lon, geod.ht, x, y, z);
    return ECEF{x, y, z};
}
