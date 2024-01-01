// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// angle conversion utilities
/**
\file
\author Steven Ward
*/

#pragma once

#include <cmath>
#include <numbers>

/// milliradians per radian
inline constexpr unsigned short mrad_per_rad = 1'000;

/// radians per revolution
inline constexpr double rad_per_rev = 2 * std::numbers::pi;

/// degrees per revolution
inline constexpr unsigned short deg_per_rev = 360;

/// arcminutes per degree
inline constexpr unsigned short arcmin_per_deg = 60;

/// arcseconds per arcminute
inline constexpr unsigned short arcsec_per_arcmin = 60;

/// degrees per radian
inline constexpr double deg_per_rad = 180 / std::numbers::pi;

inline constexpr unsigned short quadrants_per_rev = 4;
inline constexpr unsigned short sextants_per_rev = 6;
inline constexpr unsigned short octants_per_rev = 8;
inline constexpr unsigned short hexacontades_per_rev = 60;
inline constexpr unsigned short binary_degrees_per_rev = 256;
inline constexpr unsigned short gradians_per_rev = 400;

/// convert to radians from milliradians
constexpr auto
rad_from_mrad(const double x_mrad)
{
	return x_mrad / mrad_per_rad;
}

/// convert to milliradians from radians
constexpr auto
mrad_from_rad(const double x_rad)
{
	return mrad_per_rad * x_rad;
}

/// convert to revolutions from radians
constexpr auto
rev_from_rad(const double x_rad)
{
	return x_rad / rad_per_rev;
}

/// convert to radians from revolutions
constexpr auto
rad_from_rev(const double x_rev)
{
	return rad_per_rev * x_rev;
}

/// convert to revolutions from degrees
constexpr auto
rev_from_deg(const double x_deg)
{
	return x_deg / deg_per_rev;
}

/// convert to degrees from revolutions
constexpr auto
deg_from_rev(const double x_rev)
{
	return deg_per_rev * x_rev;
}

/// convert to degrees from arcminutes
constexpr auto
deg_from_arcmin(const double x_arcmin)
{
	return x_arcmin / arcmin_per_deg;
}

/// convert to arcminutes from degrees
constexpr auto
arcmin_from_deg(const double x_deg)
{
	return arcmin_per_deg * x_deg;
}

/// convert to arcminutes from arcseconds
constexpr auto
arcmin_from_arcsec(const double x_arcsec)
{
	return x_arcsec / arcsec_per_arcmin;
}

/// convert to arcseconds from arcminutes
constexpr auto
arcsec_from_arcmin(const double x_arcmin)
{
	return arcsec_per_arcmin * x_arcmin;
}

/// convert to radians from degrees
constexpr auto
rad_from_deg(const double x_deg)
{
	return x_deg / deg_per_rad;
}

/// convert to degrees from radians
constexpr auto
deg_from_rad(const double x_rad)
{
	return deg_per_rad * x_rad;
}

/// convert from degrees to degees and arcminutes
/**
\note \a arcmin will have the same sign as \a deg if the angle is not zero
\param[in] x_deg the angle (degrees)
\param[out] deg (whole number) degrees
\param[out] arcmin (decimal) arcminutes
*/
constexpr void
deg_to_dm(const double x_deg, double& deg, double& arcmin)
{
	auto tmp = x_deg;
	deg = std::trunc(tmp);

	tmp = arcmin_from_deg(tmp - deg);
	arcmin = tmp;
}

/// convert from degrees to degees, arcminutes, and arcseconds
/**
\note \a arcmin and \a arcsec will have the same sign as \a deg if the angle is not zero
\param[in] x_deg the angle (degrees)
\param[out] deg (whole number) degrees
\param[out] arcmin (whole number) arcminutes
\param[out] arcsec (decimal) arcseconds
*/
constexpr void
deg_to_dms(const double x_deg, double& deg, double& arcmin, double& arcsec)
{
	auto tmp = x_deg;
	deg = std::trunc(tmp);

	tmp = arcmin_from_deg(tmp - deg);
	arcmin = std::trunc(tmp);

	tmp = arcsec_from_arcmin(tmp - arcmin);
	arcsec = tmp;
}

/// convert to degrees from degrees and arcminutes
/**
\param[in] deg (whole number) degrees
\param[in] arcmin (decimal) arcminutes
\return (decimal) degrees
*/
constexpr auto
deg_from_dm(const double deg, const double arcmin)
{
	return deg + deg_from_arcmin(arcmin);
}

/// convert to degrees from degrees, arcminutes, and arcseconds
/**
\param[in] deg (whole number) degrees
\param[in] arcmin (whole number) arcminutes
\param[in] arcsec (decimal) arcseconds
\return (decimal) degrees
*/
constexpr auto
deg_from_dms(const double deg, const double arcmin, const double arcsec)
{
	return deg_from_dm(deg, arcmin + arcmin_from_arcsec(arcsec));
}
