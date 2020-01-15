// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// WGS 84 defining parameters, derived geometric constants, and utility functions
/**
\file
\author Steven Ward
*/

#pragma once

#include <cmath>

namespace WGS84_Ellipsoid
// {{{
{

// defining parameters

/// semi-major axis (equatorial radius of the earth) (meters)
constexpr double a = 6'378'137.0L;

/// flattening factor of the earth
constexpr double f = 1 / 298.257223563L; // (a-b)/a

// derived geometric constants

/// semi-minor axis (polar radius of the earth) (meters)
constexpr auto b = a*(1-f);

/// a squared
constexpr auto a2 = a*a;

/// b squared
constexpr auto b2 = b*b;

/// second flattening
constexpr auto fp = f/(1-f); // (a-b)/b

/// third flattening
constexpr auto n = f/(2-f); // (a-b)/(a+b)

/// first eccentricity squared
constexpr auto e2 = f*(2-f); // (a2-b2)/a2

/// first eccentricity
constexpr auto e = std::sqrt(e2);

/// second eccentricity squared
constexpr auto ep2 = e2/(1-e2); // (a2-b2)/b2

/// second eccentricity
constexpr auto ep = std::sqrt(ep2);

/// third eccentricity squared
constexpr auto epp2 = e2/(2-e2); // (a2-b2)/(a2+b2)

/// third eccentricity
constexpr auto epp = std::sqrt(epp2);

/// linear eccentricity squared
constexpr auto c2 = a2 - b2;

/// linear eccentricity
constexpr auto c = std::sqrt(c2);

/// angular eccentricity
constexpr auto alpha = std::asin(e); // std::acos(b/a)

}
//}}}


/// get the radius of curvature in the prime vertical (meters)
/**
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08
Equation (4-15)
\param sin_lat sine of the geodetic latitude
\return the radius of curvature in the prime vertical (meters)
*/
double get_Rn(const double sin_lat)
{
	const auto d2 = 1 - WGS84_Ellipsoid::e2 * sin_lat * sin_lat;
	const auto d = std::sqrt(d2);

	return WGS84_Ellipsoid::a / d;
}


/// get the ellipsoid height (meters)
/**
Source: Rapp, page 122 (132)

\verbatim
Original equations:
z = (Rn * (1-e2) + h) * sin
w = (Rn + h) * cos

Equatorial case:
h = w / cos - Rn

Polar case:
h = z / sin - Rn * (1-e2)
\endverbatim
\param w distance from the rotational (i.e. Z) axis (meters)
\param z distance above the equatorial (i.e. X-Y) plane (meters)
\param sin_lat sine of the geodetic latitude
\param cos_lat cosine of the geodetic latitude
\param Rn prime vertical radius of curvature (meters)
\return ellipsoid height (meters)
*/
double get_ht(
	const double w, const double z,
	const double sin_lat, const double cos_lat, const double Rn)
{
	// https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
	// cos(45 deg) == 1/sqrt(2)
	if (cos_lat > M_SQRT1_2) // Equatorial
	{
		return w / cos_lat - Rn;
	}
	else // Polar
	{
		return z / sin_lat - Rn * (1 - WGS84_Ellipsoid::e2);
	}
}


/// get the ellipsoid height (meters)
/**
\param w distance from the rotational (i.e. Z) axis (meters)
\param z distance above the equatorial (i.e. X-Y) plane (meters)
\param sin_lat sine of the geodetic latitude
\param cos_lat cosine of the geodetic latitude
\return ellipsoid height (meters)
*/
double get_ht(
	const double w, const double z,
	const double sin_lat, const double cos_lat)
{
	return get_ht(w, z, sin_lat, cos_lat, get_Rn(sin_lat));
}


/// get the ellipsoid height (meters)
/**
\param w distance from the rotational (i.e. Z) axis (meters)
\param z distance above the equatorial (i.e. X-Y) plane (meters)
\param lat_rad geodetic latitude (radians)
\return ellipsoid height (meters)
*/
double get_ht_r(const double w, const double z, const double lat_rad)
{
	const auto sin_lat = std::sin(lat_rad);
	const auto cos_lat = std::cos(lat_rad);

	return get_ht(w, z, sin_lat, cos_lat);
}

