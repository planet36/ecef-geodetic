// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// ECEF coordinate class
/**
\file
\author Steven Ward
*/

#pragma once

#include <cmath>
#include <concepts>
#include <fmt/core.h>
#include <sstream>
#include <string>
#include <type_traits>

constexpr int ecef_default_precision = 6; ///< The default precision

/// ECEF to string
/**
\param x X coordinate (meters)
\param y Y coordinate (meters)
\param z Z coordinate (meters)
\param precision the number of digits after the decimal place to generate for \a x, \a y, and \a z
\return the string representation of the ECEF coordinate
*/
template <std::floating_point T>
std::string ecef_to_str(
	const T x,
	const T y,
	const T z,
	int precision = ecef_default_precision)
{
	return fmt::format("{:.{}f} {:.{}f} {:.{}f}",
			x, precision,
			y, precision,
			z, precision);
}

/// string to ECEF
/**
\param s the string representation of the ECEF coordinate
\param[out] x X coordinate (meters)
\param[out] y Y coordinate (meters)
\param[out] z Z coordinate (meters)
\retval true if an error has occurred on the associated stream
*/
template <std::floating_point T>
bool str_to_ecef(
	const std::string& s,
	T& x,
	T& y,
	T& z)
{
	std::istringstream iss(s);

	iss >> x >> y >> z;

	// https://en.cppreference.com/w/cpp/io/basic_ios/operator!
	return !iss;
}
/*

IEEE Std 1278.1-2012
IEEE Standard for Distributed Interactive Simulation - Application Protocols


1.6.3.1 World coordinate system
Locations in the simulated world are identified using a right-handed, geocentric Cartesian coordinate system
called the world coordinate system. The shape of the world is described in NIMA TR 8350.2. The origin of
the coordinate system is the centroid of the World Geodetic System 1984 (WGS 84) reference frame
(ellipsoid) as defined in NIMA TR 8350.2. The axes of this system are labeled X, Y, and Z, with the positive
X-axis passing through the prime meridian at the equator, with the positive Y-axis passing through 90° east
longitude at the Equator and the positive Z-axis passing through the north pole as shown in Figure 1. A
distance of one unit measured in world coordinates corresponds to a distance of 1 m in the simulated world.
A straight line in the world coordinate system is a straight line in the simulated world. This is a rotating
reference frame that rotates on a daily period as the Earth rotates.

*/

/// Earth-Centered Earth-Fixed coordinate
/**
\sa https://en.wikipedia.org/wiki/ECEF

Source:
NGA.STND.0036_1.0.0_WGS84

\verbatim
In Figure 2.1, the origin and axes are defined as follows:

Origin = Earth's center of mass

Z-Axis = The direction of the IERS Reference Pole (IRP). This direction
 corresponds to the direction of the BIH Conventional Terrestrial
 Pole (CTP) (epoch 1984.0) with an uncertainty of 0.005" [1]

X-Axis = Intersection of the IERS Reference Meridian (IRM) and the plane
 passing through the origin and normal to the Z-axis. The IRM is
 coincident with the BIH Zero Meridian (epoch 1984.0) with an
 uncertainty of 0.005" [1]

Y-Axis = Completes a right-handed, Earth-Centered Earth-Fixed (ECEF)
 orthogonal coordinate system

The WGS 84 Coordinate System origin also serves as the geometric center of the
WGS 84 Ellipsoid and the Z-axis serves as the rotational axis of this ellipsoid of revolution
\endverbatim
*/
template <std::floating_point T>
struct ECEF
{
	using this_t = ECEF<T>;

	T x{}; // X coordinate (meters)
	T y{}; // Y coordinate (meters)
	T z{}; // Z coordinate (meters)

	ECEF() = default;

	constexpr ECEF(
		const T _x,
		const T _y,
		const T _z) :
		x(_x),
		y(_y),
		z(_z)
	{}

	/// conversion ctor
	template <std::floating_point T2>
	constexpr ECEF(const ECEF<T2>& that) :
		x(that.x),
		y(that.y),
		z(that.z)
	{}

	auto operator<=>(const this_t&) const = default;

	void normalize()
	{
	}

	std::string to_string(int precision = ecef_default_precision) const
	{
		return ecef_to_str(x, y, z, precision);
	}
};

/// ECEF<T> - ECEF<T2>
template <std::floating_point T, std::floating_point T2>
constexpr auto operator-(const ECEF<T>& p1, const ECEF<T2>& p2)
{
	using result_type = typename std::common_type_t<T, T2>;

	return ECEF<result_type>{
		p1.x - p2.x,
		p1.y - p2.y,
		p1.z - p2.z
	};
}

/// get the L1-norm
/**
\sa https://mathworld.wolfram.com/L1-Norm.html
\sa https://en.wikipedia.org/wiki/Distance#Distance_in_Euclidean_space
\param p1 the ECEF coordinate
\return the L1-norm
*/
template <std::floating_point T>
auto L1_norm(const ECEF<T>& p1)
{
	return
		std::abs(p1.x) +
		std::abs(p1.y) +
		std::abs(p1.z);
}

/// get the L2-norm
/**
\sa https://mathworld.wolfram.com/L2-Norm.html
\sa https://en.wikipedia.org/wiki/Euclidean_distance
\sa https://en.wikipedia.org/wiki/Distance#Distance_in_Euclidean_space
\param p1 the ECEF coordinate
\return the L2-norm
*/
template <std::floating_point T>
auto L2_norm(const ECEF<T>& p1)
{
	return std::sqrt(
		p1.x * p1.x +
		p1.y * p1.y +
		p1.z * p1.z);
}

/// get the Euclidean distance
/**
\sa https://en.wikipedia.org/wiki/Euclidean_distance
\sa https://en.wikipedia.org/wiki/Distance#Distance_in_Euclidean_space
\param p1 the first ECEF coordinate
\param p2 the second ECEF coordinate
\return the Euclidean distance
*/
template <std::floating_point T>
auto euclidean_dist(const ECEF<T>& p1, const ECEF<T>& p2)
{
#if 0
	return L2_norm(ECEF<T>{
		p1.x - p2.x,
		p1.y - p2.y,
		p1.z - p2.z});
#else
	return L2_norm(p1 - p2);
#endif
}
