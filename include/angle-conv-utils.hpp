// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// angle conversion utilities
/**
* \file
* \author Steven Ward
*/

#pragma once

#include <cmath>
#include <concepts>
#include <numbers>

/// milliradians per radian
inline constexpr unsigned short mrad_per_rad = 1000;

/// radians per revolution
template <std::floating_point T>
inline constexpr T rad_per_rev_v = 2 * std::numbers::pi_v<T>;

/// radians per revolution, as a \c double
inline constexpr double rad_per_rev = rad_per_rev_v<double>;

/// degrees per revolution
inline constexpr unsigned short deg_per_rev = 360;

/// arcminutes per degree
inline constexpr unsigned short arcmin_per_deg = 60;

/// arcseconds per arcminute
inline constexpr unsigned short arcsec_per_arcmin = 60;

/// degrees per radian
template <std::floating_point T>
inline constexpr T deg_per_rad_v = 180 / std::numbers::pi_v<T>;

/// degrees per radian, as a \c double
inline constexpr double deg_per_rad = deg_per_rad_v<double>;

inline constexpr unsigned short quadrants_per_rev = 4;
inline constexpr unsigned short sextants_per_rev = 6;
inline constexpr unsigned short octants_per_rev = 8;
inline constexpr unsigned short hexacontades_per_rev = 60;
inline constexpr unsigned short binary_degrees_per_rev = 256;
inline constexpr unsigned short gradians_per_rev = 400;

/// convert to radians from milliradians
template <std::floating_point T>
constexpr auto
rad_from_mrad(const T x_mrad)
{
    return x_mrad / mrad_per_rad;
}

/// convert to milliradians from radians
template <std::floating_point T>
constexpr auto
mrad_from_rad(const T x_rad)
{
    return mrad_per_rad * x_rad;
}

/// convert to revolutions from radians
template <std::floating_point T>
constexpr auto
rev_from_rad(const T x_rad)
{
    return x_rad / rad_per_rev_v<T>;
}

/// convert to radians from revolutions
template <std::floating_point T>
constexpr auto
rad_from_rev(const T x_rev)
{
    return rad_per_rev_v<T> * x_rev;
}

/// convert to revolutions from degrees
template <std::floating_point T>
constexpr auto
rev_from_deg(const T x_deg)
{
    return x_deg / deg_per_rev;
}

/// convert to degrees from revolutions
template <std::floating_point T>
constexpr auto
deg_from_rev(const T x_rev)
{
    return deg_per_rev * x_rev;
}

/// convert to degrees from arcminutes
template <std::floating_point T>
constexpr auto
deg_from_arcmin(const T x_arcmin)
{
    return x_arcmin / arcmin_per_deg;
}

/// convert to arcminutes from degrees
template <std::floating_point T>
constexpr auto
arcmin_from_deg(const T x_deg)
{
    return arcmin_per_deg * x_deg;
}

/// convert to arcminutes from arcseconds
template <std::floating_point T>
constexpr auto
arcmin_from_arcsec(const T x_arcsec)
{
    return x_arcsec / arcsec_per_arcmin;
}

/// convert to arcseconds from arcminutes
template <std::floating_point T>
constexpr auto
arcsec_from_arcmin(const T x_arcmin)
{
    return arcsec_per_arcmin * x_arcmin;
}

/// convert to radians from degrees
template <std::floating_point T>
constexpr auto
rad_from_deg(const T x_deg)
{
    return x_deg / deg_per_rad_v<T>;
}

/// convert to degrees from radians
template <std::floating_point T>
constexpr auto
deg_from_rad(const T x_rad)
{
    return deg_per_rad_v<T> * x_rad;
}

/// convert from degrees to degees and arcminutes
/**
* \note \a arcmin will have the same sign as \a deg if the angle is not zero
* \param[in] x_deg the angle (degrees)
* \param[out] deg (whole number) degrees
* \param[out] arcmin (decimal) arcminutes
*/
template <std::floating_point T>
constexpr void
deg_to_dm(const T x_deg, T& deg, T& arcmin)
{
    auto tmp = x_deg;
    deg = std::trunc(tmp);

    tmp = arcmin_from_deg(tmp - deg);
    arcmin = tmp;
}

/// convert from degrees to degees, arcminutes, and arcseconds
/**
* \note \a arcmin and \a arcsec will have the same sign as \a deg if the angle is not zero
* \param[in] x_deg the angle (degrees)
* \param[out] deg (whole number) degrees
* \param[out] arcmin (whole number) arcminutes
* \param[out] arcsec (decimal) arcseconds
*/
template <std::floating_point T>
constexpr void
deg_to_dms(const T x_deg, T& deg, T& arcmin, T& arcsec)
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
* \param[in] deg (whole number) degrees
* \param[in] arcmin (decimal) arcminutes
* \return (decimal) degrees
*/
template <std::floating_point T>
constexpr auto
deg_from_dm(const T deg, const T arcmin)
{
    return deg + deg_from_arcmin(arcmin);
}

/// convert to degrees from degrees, arcminutes, and arcseconds
/**
* \param[in] deg (whole number) degrees
* \param[in] arcmin (whole number) arcminutes
* \param[in] arcsec (decimal) arcseconds
* \return (decimal) degrees
*/
template <std::floating_point T>
constexpr auto
deg_from_dms(const T deg, const T arcmin, const T arcsec)
{
    return deg_from_dm(deg, arcmin + arcmin_from_arcsec(arcsec));
}
