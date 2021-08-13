// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Function for converting coordinates from geodetic to ECEF
/**
\file
\author Steven Ward
*/

#pragma once

#include "wgs84-utils.hpp"

#include <cmath>

/// Convert Geodetic to ECEF
/**
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08
Equation (4-14)
\param lat geodetic latitude (radians)
\param lon geodetic longitude (radians)
\param ht ellipsoid height (meters)
\param[out] x X coordinate (meters)
\param[out] y Y coordinate (meters)
\param[out] z Z coordinate (meters)
*/
void geodetic_to_ecef(const double lat, const double lon, const double ht,
                      double& x, double& y, double& z)
{
	const auto sin_lat = std::sin(lat);
	const auto cos_lat = std::cos(lat);

	const auto sin_lon = std::sin(lon);
	const auto cos_lon = std::cos(lon);

	const auto Rn = get_Rn(sin_lat);

	x = (Rn + ht) * cos_lat * cos_lon;
	y = (Rn + ht) * cos_lat * sin_lon;
	z = (Rn * (1 - WGS84_Ellipsoid::e2) + ht) * sin_lat;
}
