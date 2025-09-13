// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// ellipsoid class
/**
* \file
* \author Steven Ward
* Source:
* NGA.STND.0036_1.0.0_WGS84 2014-07-08
*
* \sa https://nsgreg.nga.mil/doc/view?i=4085
* \sa https://web.archive.org/web/20181220230431/https://earth-info.nga.mil/GandG/publications/NGA_STND_0036_1_0_0_WGS84/NGA.STND.0036_1.0.0_WGS84.pdf
*
* \verbatim
NATIONAL GEOSPATIAL-INTELLIGENCE AGENCY (NGA)
STANDARDIZATION DOCUMENT
DEPARTMENT OF DEFENSE
WORLD GEODETIC SYSTEM 1984
Its Definition and Relationships with Local Geodetic Systems
2014-07-08
Version 1.0.0
\endverbatim
*/

#pragma once

#include <cmath>
#include <concepts>
#include <stdexcept>

/// An ellipsoid and all its defining parameters and derived geometric constants
/**
* Some of the following values are not included:
*
* Table 3.2 Special WGS 84 Parameters
*
* Table 3.3 Other Fundamental Constants and Best Accepted Values
*
* Table 3.4 Special WGS 84 Parameters
*
* (some of) Table 3.5 WGS 84 Ellipsoid Derived Geometric Constants
*
* (some of) Table 3.6 WGS 84 Derived Physical Constants
*
* Table 3.7 WGS 84 Derived Moments of Inertia
*
* \sa https://arxiv.org/pdf/2212.05818.pdf
* \sa https://en.wikipedia.org/wiki/Flattening
* \sa https://en.wikipedia.org/wiki/Eccentricity_(mathematics)
*/
template <std::floating_point T>
struct Ellipsoid
{
    // defining parameters

    /// semi-major axis (equatorial radius of the earth) (meters)
    const T a;

    /// flattening factor of the earth
    /**
    * f == 0: sphere
    * f >  0: oblate ellipsoid
    * f <  0: prolate ellipsoid
    */
    const T f; // (a - b) / a

    /// geocentric gravitational constant (m³/s²)
    const T GM;

    /// nominal mean angular velocity of the earth (rad/s)
    const T omega;

    // derived geometric constants

    /// semi-minor axis (polar radius of the earth) (meters)
    const T b = a * (1 - f);

    /// a squared
    const T a2 = a * a;

    /// b squared
    const T b2 = b * b;

    /// second flattening
    const T fp = f / (1 - f); // (a - b) / b

    /// third flattening
    const T n = f / (2 - f); // (a - b) / (a + b)

    /// first eccentricity squared
    const T e2 = f * (2 - f); // (a2 - b2) / a2

    /// first eccentricity
    const T e = std::sqrt(std::abs(e2));

    /// second eccentricity squared
    const T ep2 = e2 / (1 - e2); // (a2 - b2) / b2

    /// second eccentricity
    const T ep = std::sqrt(std::abs(ep2));

    /// third eccentricity squared
    const T epp2 = e2 / (2 - e2); // (a2 - b2) / (a2 + b2)

    /// third eccentricity
    const T epp = std::sqrt(std::abs(epp2));

    /// linear eccentricity squared
    const T c2 = a2 - b2;

    /// linear eccentricity
    const T c = std::sqrt(std::abs(c2));

    /// angular eccentricity
    const T alpha = std::asin(e); // std::acos(b / a)

    // derived physical constants

    /// γₑ - normal gravity at the equator (on the ellipsoid) (m/s²)
    const T gamma_e = 9.7803253359L;

    /// γₚ - normal gravity at the poles (on the ellipsoid) (m/s²)
    const T gamma_p = 9.8321849379L;

    /// Somigliana's Formula - normal gravity formula constant
    const T k = (1 - f) * gamma_p / gamma_e - 1;

    /// normal gravity formula constant
    const T m = omega * omega * a2 * b / GM;

    Ellipsoid() = delete;

    constexpr Ellipsoid(const T a_,
                        const T f_recip_, // 1 / f
                        const T GM_ = 3.986004418E14L,
                        const T omega_ = 7.292115E-5L) :
    a(a_),
    f(1 / f_recip_),
    GM(GM_),
    omega(omega_)
    {
        if (a <= T{})
            throw std::invalid_argument("Equatorial radius must be positive");

        if (!std::isfinite(a))
            throw std::invalid_argument("Equatorial radius must be finite");

        if (b <= T{})
            throw std::invalid_argument("Polar radius must be positive");

        if (!std::isfinite(b))
            throw std::invalid_argument("Polar radius must be finite");
    }

    /// get the radius of curvature in the prime vertical (meters)
    /**
    * Source:
    * NGA.STND.0036_1.0.0_WGS84 2014-07-08
    * Equation (4-15)
    * \param sin_lat sine of the geodetic latitude
    * \return the radius of curvature in the prime vertical (meters)
    */
    [[nodiscard]] auto get_Rn(const T sin_lat) const
    {
        const auto d2 = 1 - e2 * sin_lat * sin_lat;
        const auto d = std::sqrt(d2);

        return a / d;
    }

    /**
    * Derivation of ellipsoid radius:
    *
    * https://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_center
    *
    * R(θ) = b / √(1 - e² * cos(θ)²)
    *
    * https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    *
    * θ(φ) = atan(tan(φ) * (1 - e²))
    *
    * https://www.wolframalpha.com/input/?i=simplify+cos%28atan%28x%29%29**2
    *
    * cos(atan(x))² = 1 / (x² + 1)
    *
    * cos(θ)² = cos(φ)² / ((1 - e²)² * sin(φ)² + cos(φ)²)
    *
    * R(φ) = b / √(1 - e² * cos(φ)² / ((1 - e²)² * sin(φ)² + cos(φ)²))
    *      = R_N * √((1 - e²)² * sin(φ)² + cos(φ)²)
    *      = R_N * √(1 - e² * sin(φ)² * (2 - e²))
    */

    /// get the ellipsoid radius (meters)
    /**
    * \sa https://www.oc.nps.edu/oc2902w/geodesy/radiigeo.pdf
    * \param sin_lat sine of the geodetic latitude
    * \return the ellipsoid radius (meters)
    */
    [[nodiscard]] auto get_R(const T sin_lat) const
    {
        return get_Rn(sin_lat) * std::sqrt(1 - e2 * sin_lat * sin_lat * (2 - e2));
    }

    /// get the radius of curvature in the meridian (meters)
    /**
    * Source:
    * NGA.STND.0036_1.0.0_WGS84 2014-07-08
    * Page 7-8
    * \param sin_lat sine of the geodetic latitude
    * \return the radius of curvature in the meridian (meters)
    */
    [[nodiscard]] auto get_Rm(const T sin_lat) const
    {
        const auto d2 = 1 - e2 * sin_lat * sin_lat;
        const auto d = std::sqrt(d2);

        return a * (1 - e2) / (d2 * d);
    }

    /// get the normal gravity on the ellipsoid surface (m/s²)
    /**
    * Source:
    * NGA.STND.0036_1.0.0_WGS84 2014-07-08
    * Page 4-1
    * \param sin_lat sine of the geodetic latitude
    * \return the normal gravity on the ellipsoid surface (m/s²)
    */
    [[nodiscard]] auto get_gamma(const T sin_lat) const
    {
        const auto d2 = 1 - e2 * sin_lat * sin_lat;
        const auto d = std::sqrt(d2);

        return gamma_e * (1 + k * sin_lat * sin_lat) / d;
    }

    /// get the normal gravity above the ellipsoid (m/s²)
    /**
    * Source:
    * NGA.STND.0036_1.0.0_WGS84 2014-07-08
    * Page 4-3
    * \param sin_lat sine of the geodetic latitude
    * \param ht ellipsoid height (meters)
    * \return the normal gravity above the ellipsoid (m/s²)
    */
    [[nodiscard]] auto get_gamma_h(const T sin_lat, const T ht) const
    {
        return get_gamma(sin_lat) *
               (1 - 2 * ht * (1 + f + m - 2 * f * sin_lat * sin_lat) / a + 3 * ht * ht / a2);
    }

    /// get the height above the ellipsoid (meters)
    /**
    * Source: Rapp, page 122 (132)
    *
    * \verbatim
    Original equations:
    z = (Rn * (1-e2) + h) * sin
    w = (Rn + h) * cos

    Equatorial case:
    h = w / cos - Rn

    Polar case:
    h = z / sin - Rn * (1-e2)
    \endverbatim
    * \param w distance from the rotational (i.e. Z) axis (meters)
    * \param z distance above the equatorial (i.e. X-Y) plane (meters)
    * \param sin_lat sine of the geodetic latitude
    * \param cos_lat cosine of the geodetic latitude
    * \param Rn prime vertical radius of curvature (meters)
    * \return the height above the ellipsoid (meters)
    */
    [[nodiscard]] auto get_ht(const T w,
                              const T z,
                              const T sin_lat,
                              const T cos_lat,
                              const T Rn) const
    {
        // https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
        // cos(45 deg) == 1/sqrt(2)
        if (cos_lat > M_SQRT1_2) // Equatorial
            return w / cos_lat - Rn;
        else // Polar
            return z / sin_lat - Rn * (1 - e2);
    }

    /// get the height above the ellipsoid (meters)
    /**
    * \param w distance from the rotational (i.e. Z) axis (meters)
    * \param z distance above the equatorial (i.e. X-Y) plane (meters)
    * \param sin_lat sine of the geodetic latitude
    * \param cos_lat cosine of the geodetic latitude
    * \return the height above the ellipsoid (meters)
    */
    [[nodiscard]] auto get_ht(const T w, const T z, const T sin_lat, const T cos_lat) const
    {
        return get_ht(w, z, sin_lat, cos_lat, get_Rn(sin_lat));
    }

    /// get the height above the ellipsoid (meters)
    /**
    * \param w distance from the rotational (i.e. Z) axis (meters)
    * \param z distance above the equatorial (i.e. X-Y) plane (meters)
    * \param lat_rad geodetic latitude (radians)
    * \return the height above the ellipsoid (meters)
    */
    [[nodiscard]] auto get_ht(const T w, const T z, const T lat_rad) const
    {
        const auto sin_lat = std::sin(lat_rad);
        const auto cos_lat = std::cos(lat_rad);

        return get_ht(w, z, sin_lat, cos_lat);
    }

    bool operator==(const Ellipsoid& that) const
    {
        return this->a == that.a && this->f == that.f && this->GM == that.GM &&
               this->omega == that.omega;
    }

    bool operator!=(const Ellipsoid& that) const { return !operator==(that); }
};
