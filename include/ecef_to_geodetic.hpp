// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// ECEF-to-Geodetic coordinate conversion
/**
\file
\author Steven Ward
*/

#pragma once

#include "angle.hpp"
#include "ecef-coord.hpp"
#include "ellipsoid-wgs84.hpp"
#include "geodetic-coord.hpp"

#include <cmath>
#include <concepts>

/// convert from ECEF to geodetic
/**
D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to geodetic
coordinates," in IEEE Transactions on Aerospace and Electronic Systems, vol.
32, no. 1, pp. 473-476, Jan. 1996, doi: 10.1109/7.481290.

U.S. Government work, U.S. copyright does not apply.

Converted to C++ and modified by Steven Ward.  No rights reserved.

\sa https://ieeexplore.ieee.org/document/481290
\param x X coordinate (meters)
\param y Y coordinate (meters)
\param z Z coordinate (meters)
\param[out] lat_rad geodetic latitude (radians)
\param[out] lon_rad geodetic longitude (radians)
\param[out] ht ellipsoid height (meters)
*/
template <std::floating_point T>
void
ecef_to_geodetic(const T x,
                 const T y,
                 const T z,
                 T& lat_rad,
                 T& lon_rad,
                 T& ht)
{
	static constexpr auto& ell = WGS84<T>;

	const auto w2 = x * x + y * y;
	const auto w = std::sqrt(w2);
	const auto z2 = z * z;

	// atan2 returns unintuitive values when given zeros
	// https://en.cppreference.com/w/cpp/numeric/math/atan2
	if (w != 0)
		lon_rad = std::atan2(y, x);
	else // on the axis of rotation
		lon_rad = 0;

	constexpr auto a1 = ell.a * ell.e2;
	constexpr auto a2 = a1 * a1;
	constexpr auto a3 = a1 * ell.e2 / 2;
	constexpr auto a4 = 2.5 * a2;
	constexpr auto a5 = a1 + a3;
	//constexpr auto a6 = (1 - ell.e2);

	const auto r2 = w2 + z2;
	const auto r = std::sqrt(r2);

	const auto s2 = z2 / r2;
	const auto c2 = w2 / r2;
	auto u = a2 / r;
	auto v = a3 - a4 / r;

	T s{};
	T c{};
	T ss{};

	// cos(45 deg)^2 == 0.5
	if (c2 > 0.5) // Equatorial
	{
		s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r);
		lat_rad = std::asin(s);
		ss = s * s;
		c = std::sqrt(1 - ss);
	}
	else // Polar
	{
		c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r);
		lat_rad = std::acos(c);
		ss = 1 - c * c;
		s = std::sqrt(ss);

		if (z < 0)
		{
			lat_rad = -lat_rad;
			s = -s;
		}
	}

	const auto d2 = 1 - ell.e2 * ss;
	const auto Rn = ell.a / std::sqrt(d2);
	const auto Rm = (1 - ell.e2) * Rn / d2;
	const auto rf = (1 - ell.e2) * Rn;
	u = w - Rn * c;
	v = z - rf * s;
	const auto f = c * u + s * v;
	const auto m = c * v - s * u;
	const auto p = m / (Rm + f);

	lat_rad += p;

#if 1
	// custom ht
	ht = f + m * p / 2;
#else
	// common ht
	ht = ell.get_ht(w, z, std::sin(lat_rad), std::cos(lat_rad));
#endif
}

/// convert from ECEF to geodetic
/**
D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to geodetic
coordinates," in IEEE Transactions on Aerospace and Electronic Systems, vol.
32, no. 1, pp. 473-476, Jan. 1996, doi: 10.1109/7.481290.

U.S. Government work, U.S. copyright does not apply.

Converted to C++ and modified by Steven Ward.  No rights reserved.

\sa https://ieeexplore.ieee.org/document/481290
\param x X coordinate (meters)
\param y Y coordinate (meters)
\param z Z coordinate (meters)
\param[out] lat geodetic latitude
\param[out] lon geodetic longitude
\param[out] ht ellipsoid height (meters)
*/
template <angle_unit U, std::floating_point T>
void
ecef_to_geodetic(const T x,
                 const T y,
                 const T z,
                 angle<U, T>& lat,
                 angle<U, T>& lon,
                 T& ht)
{
	static constexpr auto& ell = WGS84<T>;

	const auto w2 = x * x + y * y;
	const auto w = std::sqrt(w2);
	const auto z2 = z * z;

	// atan2 returns unintuitive values when given zeros
	// https://en.cppreference.com/w/cpp/numeric/math/atan2
	if (w != 0)
		lon = a_atan2(y, x);
	else // on the axis of rotation
		lon = 0;

	constexpr auto a1 = ell.a * ell.e2;
	constexpr auto a2 = a1 * a1;
	constexpr auto a3 = a1 * ell.e2 / 2;
	constexpr auto a4 = 2.5 * a2;
	constexpr auto a5 = a1 + a3;
	//constexpr auto a6 = (1 - ell.e2);

	const auto r2 = w2 + z2;
	const auto r = std::sqrt(r2);

	const auto s2 = z2 / r2;
	const auto c2 = w2 / r2;
	auto u = a2 / r;
	auto v = a3 - a4 / r;

	T s{};
	T c{};
	T ss{};

	// cos(45 deg)^2 == 0.5
	if (c2 > 0.5) // Equatorial
	{
		s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r);
		lat = a_asin(s);
		ss = s * s;
		c = std::sqrt(1 - ss);
	}
	else // Polar
	{
		c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r);
		lat = a_acos(c);
		ss = 1 - c * c;
		s = std::sqrt(ss);

		if (z < 0)
		{
			lat = -lat;
			s = -s;
		}
	}

	const auto d2 = 1 - ell.e2 * ss;
	const auto Rn = ell.a / std::sqrt(d2);
	const auto Rm = (1 - ell.e2) * Rn / d2;
	const auto rf = (1 - ell.e2) * Rn;
	u = w - Rn * c;
	v = z - rf * s;
	const auto f = c * u + s * v;
	const auto m = c * v - s * u;
	const auto p = m / (Rm + f);

	lat += p;

#if 1
	// custom ht
	ht = f + m * p / 2;
#else
	// common ht
	ht = ell.get_ht(w, z, sin(lat), cos(lat));
#endif
}

template <std::floating_point T>
auto ecef_to_geodetic(const ECEF<T>& ecef)
{
	T lat_rad{};
	T lon_rad{};
	T ht{};
	ecef_to_geodetic(ecef.x, ecef.y, ecef.z, lat_rad, lon_rad, ht);
	return Geodetic<angle_unit::radian, T>{lat_rad, lon_rad, ht};
}
