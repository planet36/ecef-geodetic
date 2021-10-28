// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Functions for converting coordinates from ECEF to geodetic
/**
\file
\author Steven Ward
*/

#pragma once

#include "ellipsoid-wgs84.hpp"

#include <cmath>
#include <concepts>
#include <string>

/// the ellipsoid to use
constexpr auto ell = WGS84<double>;

/// get the 2D hypotenuse
/**
\sa https://mathworld.wolfram.com/Norm.html
\sa https://mathworld.wolfram.com/VectorNorm.html
\sa https://mathworld.wolfram.com/L2-Norm.html
\sa https://en.cppreference.com/w/cpp/numeric/math/hypot
\param[in,out] x,y the vector coordinates
\return the 2D hypotenuse
*/
template <std::floating_point T>
auto hypot(const T x, const T y)
{
#if 0
	// SDW: this is a little more accurate, but much slower
	return std::hypot(x, y);
#else
	return std::sqrt(x * x + y * y);
#endif
}

/// make a unit vector (i.e. normalize the vector components)
/**
\pre The magnitude of &lt;\a x, \a y&gt; is non-zero.
\sa https://mathworld.wolfram.com/NormalizedVector.html
\sa https://mathworld.wolfram.com/UnitVector.html
\param[in,out] x,y the vector coordinates
*/
template <std::floating_point T>
void normalize(T& x, T& y)
{
	const auto h = hypot(x, y);
	x /= h;
	y /= h;
}

/// get the cosine of the angle, given the sine of the angle
/**
\pre original angle is in the interval [-90°, +90°]
\pre abs(\a sin_x) <= 1
\param sin_x the sine of the angle
\return the cosine of the angle
*/
template <std::floating_point T>
auto cos_from_sin(const T sin_x)
{
	return std::sqrt(1 - sin_x * sin_x);
}

// these are the lines in the function params and after the function body
//constexpr int _lines_to_ignore = 7+2;

// POW2
#define SQ(X) ((X) * (X))
// POW3
#define CB(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))
#define POW5(X) ((X) * (X) * (X) * (X) * (X))
#define POW6(X) ((X) * (X) * (X) * (X) * (X) * (X))

// common declarations for the ECEF-to-geodetic functions
#define COMMON_FIRST_DECLS \
	const auto w2 = x * x + y * y; \
	const auto w = std::sqrt(w2); \
	[[gnu::unused]] const auto z2 = z * z; \
	lon_rad = std::atan2(y, x); \

// these are the lines in the common first decls
constexpr int _lines_common_first_decls = 4;

// common declarations for the ECEF-to-geodetic functions
// that need code for corner cases (i.e. equator or the poles)
#define COMMON_FIRST_DECLS_CHECKED \
	const auto w2 = x * x + y * y; \
	if (w2 != 0) \
	{ \
		lon_rad = std::atan2(y, x); \
	} \
	else /* on the axis of rotation */ \
	{ \
		lon_rad = 0; \
		if (z == 0) /* center of earth */ \
		{ \
			lat_rad = 0; \
			ht = -ell.a; \
		} \
		else \
		{ \
			lat_rad = std::copysign(M_PI_2, z); \
			ht = std::abs(z) - ell.b; \
		} \
		return; \
	} \
	const auto w = std::sqrt(w2); \
	if (z == 0) /* on the equatorial plane */ \
	{ \
		lat_rad = 0; \
		ht = w - ell.a; \
		return; \
	} \
	[[gnu::unused]] const auto z2 = z * z; \

// these are the lines in the common first decls (checked)
constexpr int _lines_common_first_decls_checked = 28;

//#undef COMMON_FIRST_DECLS
//#define COMMON_FIRST_DECLS COMMON_FIRST_DECLS_CHECKED

struct func_info_t
{
	void (&func_ref)(const double, const double, const double,
	                 double&, double&, double&);
	int                    num_lines;
	const bool             needs_code_for_corner_cases;
	const int              ilog10_mean_dist_err;
	const std::string_view algo_author;
	const std::string_view code_copyright;
	const std::string_view license;
	const std::string_view orig_impl_lang;
	const std::string_view url;
	const std::string_view citation;

	func_info_t(
			void (&_func_ref)(const double, const double, const double,
			                  double&, double&, double&),
			const int              _num_lines,
			const bool             _needs_code_for_corner_cases,
			const int              _ilog10_mean_dist_err,
			const std::string_view _algo_author,
			const std::string_view _code_copyright,
			const std::string_view _license,
			const std::string_view _orig_impl_lang,
			const std::string_view _url,
			const std::string_view _citation
			):

		func_ref                    ( _func_ref                    ) ,
		num_lines                   ( _num_lines                   ) ,
		needs_code_for_corner_cases ( _needs_code_for_corner_cases ) ,
		ilog10_mean_dist_err        ( _ilog10_mean_dist_err        ) ,
		algo_author                 ( _algo_author                 ) ,
		code_copyright              ( _code_copyright              ) ,
		license                     ( _license                     ) ,
		orig_impl_lang              ( _orig_impl_lang              ) ,
		url                         ( _url                         ) ,
		citation                    ( _citation                    )
	{
		if (needs_code_for_corner_cases)
		{
			num_lines += _lines_common_first_decls_checked;
		}
		else
		{
			num_lines += _lines_common_first_decls;
		}
	}
};

/**
* Derivation of naive II algorithm...
*
* Another source:
*
* Geometric Geodesy
* Part 1
* by
* Richard H. Rapp
* The Ohio State University
* April 1991
*
* Page 121 (131)
*/

/**
* Derivation (starting with naive I algorithm):
*
* tan == z * (Rn + ht) / (w * (Rn * (1-e2) + ht))
*
*
* Equatorial:
*
* w / (Rn + ht) == cos
*
* tan == z / (cos * (Rn * (1-e2) + ht))
*
* cos * (Rn * (1-e2) + ht) == w - cos * Rn * e2
*
* tan == z / (w - cos * Rn * e2)
*
*
* Polar:
*
* z / (Rn * (1-e2) + ht) == sin
*
* tan == sin * (Rn + ht) / w
*
* sin * (Rn + ht) == z + sin * Rn * e2
*
* tan == (z + sin * Rn * e2) / w
*/

/**
* Original equation: (from Hirvonen)
*
* Equatorial:
* tan = z / (w - cos * Rn * e2)
*
* Polar:
* tan = (z + sin * Rn * e2) / w
*
* (z + sin * Rn * e2) / w = z / (w - cos * Rn * e2)
* (z + sin * Rn * e2) * (w - cos * Rn * e2) = z * w
* z * w + w * sin * Rn * e2 - z * cos * Rn * e2 - sin * cos * Rn**2 * e4 = z * w
* w * sin * Rn * e2 - z * cos * Rn * e2 - sin * cos * Rn**2 * e4 = 0
* Rn * e2 * (w * sin - z * cos - sin * cos * Rn * e2) = 0
* w * sin - z * cos - sin * cos * Rn * e2 = 0
* Rn * e2 * sin * cos - w * sin + z * cos = 0
*
* let d = sqrt(1 - e2 * sin**2)
* Rn = a / d
*
* Original equation:
* f(φ) = Rn * e2 * sin * cos - w * sin + z * cos = 0
*
* first derivative: (unsimplified)
* f'(φ) = a * e2**2 * sin**2 * cos**2 / d**3 - a * e2 * sin**2 / d + a * e2 * cos**2 / d - w * cos - z * sin = 0
*
* first derivative: (simplified)
* f'(φ) = Rn * (d**2 - (1-e2) / d**2) - w * cos - z * sin = 0
*
* second derivative: (unsimplified)
* f''(φ) = a * e2 * sin * cos * (d**2 - (1-e2) / d**2) / d**3 - 2 * a * e2 * sin * cos * ((1-e2) / d**4 + 1) / d + w * sin - z * cos = 0
*
* second derivative: (simplified)
* f''(φ) = -Rn * e2 * sin * cos * (3 * (1-e2) / d**4 + 1) + w * sin - z * cos = 0
* OR
* f''(φ) = -Rn * e2 * sin * cos * (3 * (1-e2) / d**2 + d**2) / d**2 + w * sin - z * cos = 0
*/

/// get f, f'
template <std::floating_point T>
void get_f_fp(const T w, const T z, const T sin_lat, const T cos_lat,
              T& f, T& fp)
{
	const auto d2 = 1 - WGS84<T>.e2 * sin_lat * sin_lat;
	const auto d = std::sqrt(d2);
	const auto Rn = WGS84<T>.a / d;
	f = Rn * WGS84<T>.e2 * sin_lat * cos_lat - w * sin_lat + z * cos_lat;
	fp = Rn * (d2 - (1 - WGS84<T>.e2) / d2) - w * cos_lat - z * sin_lat;
}

constexpr int _lines_f_fp = 10;

/// get f, f', f''
template <std::floating_point T>
void get_f_fp_fpp(const T w, const T z, const T sin_lat, const T cos_lat,
                  T& f, T& fp, T& fpp)
{
	const auto d2 = 1 - WGS84<T>.e2 * sin_lat * sin_lat;
	const auto d = std::sqrt(d2);
	const auto Rn = WGS84<T>.a / d;
	f = Rn * WGS84<T>.e2 * sin_lat * cos_lat - w * sin_lat + z * cos_lat;
	fp = Rn * (d2 - (1 - WGS84<T>.e2) / d2) - w * cos_lat - z * sin_lat;
	fpp = -Rn * WGS84<T>.e2 * sin_lat * cos_lat * (d2 + 3 * (1 - WGS84<T>.e2) / d2) / d2 + w * sin_lat - z * cos_lat;
}

constexpr int _lines_f_fp_fpp = 11;

template <std::floating_point T>
auto newton_raphson_delta_lat(const T w, const T z,
                              const T sin_lat, const T cos_lat)
{
	T f, fp;
	get_f_fp(w, z, sin_lat, cos_lat, f, fp);
	return f / fp;
}

constexpr int _lines_newton_raphson_delta_lat = 8 + _lines_f_fp;

template <std::floating_point T>
auto householder_delta_lat(const T w, const T z,
                           const T sin_lat, const T cos_lat)
{
	T f, fp, fpp;
	get_f_fp_fpp(w, z, sin_lat, cos_lat, f, fp, fpp);
	return (f / fp) * (1 + 0.5 * f * fpp / (fp * fp));
}

constexpr int _lines_householder_delta_lat = 8 + _lines_f_fp_fpp;

template <std::floating_point T>
auto schroder_delta_lat(const T w, const T z,
                        const T sin_lat, const T cos_lat)
{
	T f, fp, fpp;
	get_f_fp_fpp(w, z, sin_lat, cos_lat, f, fp, fpp);
	return (f * fp) / (fp * fp - f * fpp);
}

constexpr int _lines_schroder_delta_lat = 8 + _lines_f_fp_fpp;

template <std::floating_point T>
auto halley_delta_lat(const T w, const T z,
                      const T sin_lat, const T cos_lat)
{
	T f, fp, fpp;
	get_f_fp_fpp(w, z, sin_lat, cos_lat, f, fp, fpp);
	return (f * fp) / (fp * fp - 0.5 * f * fpp);
}

constexpr int _lines_halley_delta_lat = 8 + _lines_f_fp_fpp;

template <std::floating_point T>
auto ligas_f1(const T w, const T we,
              const T z, const T ze)
{
	return (1 - WGS84<T>.e2) * we * (ze - z) - ze * (we - w);
}

template <std::floating_point T>
auto ligas_f2(const T we, const T ze)
{
	return (1 - WGS84<T>.e2) * we * we + ze * ze - WGS84<T>.b2;
}

template <std::floating_point T>
auto det(const T A[2][2])
{
	return A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

template <std::floating_point T>
void inv(const T A[2][2], T result[2][2])
{
	const auto d = det(A);

	result[0][0] = +A[1][1] / d;
	result[0][1] = -A[0][1] / d;
	result[1][0] = -A[1][0] / d;
	result[1][1] = +A[0][0] / d;
}

template <std::floating_point T>
void mul(const T A[2][2], const T X[2], T result[2])
{
	result[0] = A[0][0] * X[0] + A[0][1] * X[1];
	result[1] = A[1][0] * X[0] + A[1][1] * X[1];
}

template <std::floating_point T>
void ligas_Jacobian(const T w, const T we, const T z, const T ze,
                    T result[2][2])
{
	result[0][0] = (1 - WGS84<T>.e2) * (ze - z) - ze;
	result[0][1] = (1 - WGS84<T>.e2) * we - (we - w);
	result[1][0] = 2 * (1 - WGS84<T>.e2) * we;
	result[1][1] = 2 * ze;
}

constexpr int _lines_ligas_util = 46;

template <std::floating_point T>
auto lin_wang_1995_delta_m(const T w2, const T z2, const T m)
{
	const auto tmp_a = WGS84<T>.a + 2 * m / WGS84<T>.a;
	const auto tmp_b = WGS84<T>.b + 2 * m / WGS84<T>.b;

	const auto f = w2 / (tmp_a * tmp_a) + z2 / (tmp_b * tmp_b) - 1;
	const auto fp = -4 * (w2 / (WGS84<T>.a * CB(tmp_a)) + z2 / (WGS84<T>.b * CB(tmp_b)));

	return f / fp;
}

constexpr int _lines_lin_wang_1995_delta_m = 11;

template <std::floating_point T>
auto shu_2010_delta_k(const T w2, const T z2, const T k)
{
	const auto p = WGS84<T>.a + WGS84<T>.b * k;
	const auto q = WGS84<T>.b + WGS84<T>.a * k;
	const auto p2 = p * p;
	const auto q2 = q * q;

	const auto f = p2 * q2 - w2 * q2 - z2 * p2;
	const auto fp = 2 * (WGS84<T>.b * p * (q2 - z2) + WGS84<T>.a * q * (p2 - w2));

	return f / fp;
}

constexpr int _lines_shu_2010_delta_k = 13;

template <std::floating_point T>
auto wu_2003_delta_t(const T A, const T B, const T C, const T t)
{
	const auto t2 = t * t;
	const auto t3 = t * t * t;
	const auto t4 = t * t * t * t;

	const auto f = (A * t4 + (B + C) * t3 + (B - C) * t - A);
	const auto fp = (4 * A * t3 + 3 * (B + C) * t2 + (B - C));

	return f / fp;
}

constexpr int _lines_wu_2003_delta_t = 12;

namespace borkowski_1989
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// a*e2/(1-f) == b*ep2
	// SDW: NaNs at the poles
	//const auto E = (ell.b * std::abs(z) - ell.a2 * ell.e2) / (ell.a * w);
	//const auto F = (ell.b * std::abs(z) + ell.a2 * ell.e2) / (ell.a * w);
	const auto E = ((1 - ell.f) * std::abs(z) - ell.a * ell.e2) / w;
	const auto F = ((1 - ell.f) * std::abs(z) + ell.a * ell.e2) / w;

	const auto P = 4 * (E * F + 1) / 3;
	const auto Q = 2 * (E * E - F * F);
	const auto D = P * P * P + Q * Q;
	double v;
	if (D < 0)
	{
		const auto sqrt_n_P = std::sqrt(-P);
		v = 2 * sqrt_n_P * std::cos(std::acos(Q / CB(sqrt_n_P)) / 3);
	}
	else
	{
		const auto sqrt_D = std::sqrt(D);
		v = std::cbrt(sqrt_D - Q) - std::cbrt(sqrt_D + Q);
	}
	const auto G = (std::sqrt(E * E + v) + E) / 2;
	const auto t = std::sqrt(G * G + (F - v * G) / (2 * G - E)) - G;

	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
	auto sin_lat = 1 - t * t;
	auto cos_lat = 2 * t * (1 - ell.f);

	if (z < 0)
		sin_lat = -sin_lat;

	lat_rad = std::atan2(sin_lat, cos_lat);

	// http://www.astro.uni.torun.pl/~kb/Papers/ASS/Geod-ASS.htm
	// SDW: this does not work
	//ht = (r - ell.a * t) * cos_lat + (z - ell.b) * sin_lat;

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "K.M. Borkowski",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm",
	/*.citation                    =*/ R"(Borkowski, K.M. Bull. Geodesique (1989) 63: 50. https://doi.org/10.1007/BF02520228
https://link.springer.com/article/10.1007/BF02520228)"
);

}
// }}}

namespace bowring_1976_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// (i = 0) geocentric to parametric
	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.f);
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "B.R. Bowring",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/233496226_Transformation_from_spatial_to_geographical_coordinates",
	/*.citation                    =*/ R"(R. Bowring, B. (1976). Transformation from spatial to geographical coordinates. Survey Review. 23. 323-327. 10.1179/003962676791280626.
https://www.tandfonline.com/doi/abs/10.1179/sre.1976.23.181.323)"
);

}
// }}}

namespace bowring_1976_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// (i = 0) geocentric to parametric
	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.f);
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	// geodetic to parametric
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);

	// (i = 2) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "B.R. Bowring",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/233496226_Transformation_from_spatial_to_geographical_coordinates",
	/*.citation                    =*/ R"(R. Bowring, B. (1976). Transformation from spatial to geographical coordinates. Survey Review. 23. 323-327. 10.1179/003962676791280626.
https://www.tandfonline.com/doi/abs/10.1179/sre.1976.23.181.323)"
);

}
// }}}

namespace bowring_1985_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	// (i = 0) ??? to parametric
	auto sin_lat = z * (1 - ell.f + ell.a * ell.e2 / r);
	// SDW: this is not more accurate
	//auto sin_lat = z * (1 - ell.f) + ell.a * ell.e2 * (z / r);
	auto cos_lat = w;
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "B.R. Bowring",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://gis.stackexchange.com/questions/28446/computational-most-efficient-way-to-convert-cartesian-to-geodetic-coordinates",
	/*.citation                    =*/ R"(Bowring's Method with improved initial guess [1985])"
);

}
// }}}

namespace bowring_1985_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	// (i = 0) ??? to parametric
	auto sin_lat = z * (1 - ell.f + ell.a * ell.e2 / r);
	// SDW: this is not more accurate
	//auto sin_lat = z * (1 - ell.f) + ell.a * ell.e2 * (z / r);
	auto cos_lat = w;
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	// geodetic to parametric
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);

	// (i = 2) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "B.R. Bowring",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://gis.stackexchange.com/questions/28446/computational-most-efficient-way-to-convert-cartesian-to-geodetic-coordinates",
	/*.citation                    =*/ R"(Bowring's Method with improved initial guess [1985])"
);

}
// }}}

namespace bowring_toms_1995_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);
	double aD_b;

	if (r <= 2E6 + (ell.a + ell.b) / 2)
	{
		// region 1
		aD_b = 1.0026;
	}
	else if (r <= 6E6 + (ell.a + ell.b) / 2)
	{
		// region 2
		aD_b = 1.00092592;
	}
	else if (r <= 18E6 + (ell.a + ell.b) / 2)
	{
		// region 3
		aD_b = 0.999250297;
	}
	else // if (r <= 1E9 + (ell.a + ell.b) / 2)
	{
		// region 4
		aD_b = 0.997523508;
	}

	// (i = 0) geocentric to parametric
	auto sin_lat = z * aD_b;
	auto cos_lat = w;
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -3,
	/*.algo_author                 =*/ "B.R. Bowring, R. Toms",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.osti.gov/biblio/231228",
	/*.citation                    =*/ R"(Ralph M. Toms, An Improved Algorithm for Geocentric to Geodetic Coordinate Conversion)"
);

}
// }}}

namespace bowring_toms_1995_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);
	double aD_b;

	if (r <= 2E6 + (ell.a + ell.b)/2)
	{
		// region 1
		aD_b = 1.0026;
	}
	else if (r <= 6E6 + (ell.a + ell.b)/2)
	{
		// region 2
		aD_b = 1.00092592;
	}
	else if (r <= 18E6 + (ell.a + ell.b)/2)
	{
		// region 3
		aD_b = 0.999250297;
	}
	else // if (r <= 1E9 + (ell.a + ell.b)/2)
	{
		// region 4
		aD_b = 0.997523508;
	}

	// (i = 0) geocentric to parametric
	auto sin_lat = z * aD_b;
	auto cos_lat = w;
	normalize(cos_lat, sin_lat);

	// (i = 1) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	// geodetic to parametric
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);

	// (i = 2) parametric to geodetic
	sin_lat = z + ell.b * ell.ep2 * CB(sin_lat);
	cos_lat = w - ell.a * ell.e2 * CB(cos_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "B.R. Bowring, R. Toms",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.osti.gov/biblio/231228",
	/*.citation                    =*/ R"(Ralph M. Toms, An Improved Algorithm for Geocentric to Geodetic Coordinate Conversion)"
);

}
// }}}

namespace fukushima_1999_1
// {{{
{

template <std::floating_point T>
auto f(const T t, const T u, const T v, const T w)
{
	// w * t**4 + u * t**3 + v * t - w
	return w * t * t * t * t + u * t * t * t + v * t - w;
}

template <std::floating_point T>
auto fp(const T t, const T u, const T v, const T w)
{
	// 4 * w * t**3 + 3 * u * t**2 + v
	return 4 * w * t * t * t + 3 * u * t * t + v;
}

template <std::floating_point T>
auto fpp(const T t, const T u, [[gnu::unused]] const T v, const T w)
{
	// 12 * w * t**2 + 6 * u * t
	return 12 * w * t * t + 6 * u * t;
}

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	double sin_lat = 0;
	double cos_lat = 0;
	constexpr auto c = ell.a * ell.e2;
	constexpr auto ep = 1 - ell.f;
	const auto zp = ep * std::abs(z);
	const auto u = 2 * (zp - c);
	const auto v = 2 * (zp + c);

	const auto tM = (c - zp) / w;

	double t = 0;

	if (tM <= 0)
	{
		// Case 1
		t = (w - c + zp) / (w - c + 2 * zp);
	}
	else if (tM >= 1)
	{
		// Case 2
		t = w / (zp + c);
	}
	else
	{
		auto fM = f(tM, u, v, w);

		// Case 3
		if (fM >= 0)
		{
			// Case 3a
			// (same as Case 2)
			t = w / (zp + c);
		}
		else // fM < 0
		{
			// Case 3b
			// (same as Case 1)
			t = (w - c + zp) / (w - c + 2 * zp);
		}
	}

	/*
	// SDW: This does not work
	sin_lat = z;
	cos_lat = w * (1 - ell.f);

	t = std::tan(sin_lat / cos_lat);
	*/

	// i = 1
	t -= f(t, u, v, w) / fp(t, u, v, w);

	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
#if 0
	// SDW: This does not work
	sin_lat = 2 * t;
	cos_lat = (1 - t * t) * ep;

	// SDW: This does not work
	sin_lat = 2 * t * ep;
	cos_lat = (1 - t * t);
#else
	// SDW: This does not work
	sin_lat = (1 - t * t) * ep;
	cos_lat = 2 * t;

	// SDW: This does not work
	//sin_lat = (1 - t * t);
	//cos_lat = 2 * t * ep;
#endif

	if (z < 0)
		sin_lat = -sin_lat;

	// 2 * tan(x/2) / (1 - tan(x/2)^2) == tan(x)
	// https://www.wolframalpha.com/input/?i=2+*+tan(x%2F2)+%2F+(1+-+tan(x%2F2)%5E2),+tan(x)

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = (2 * w * ep * t + z * (1 - t * t) - ell.a * ep * (1 + t * t)) / std::sqrt(SQ(1 + t * t) - 4 * ell.e2 * t * t);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 20;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "Toshio Fukushima",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226311152_Fast_transform_from_geocentric_to_geodetic_coordinates",
	/*.citation                    =*/ R"(Fukushima, T. Journal of Geodesy (1999) 73: 603. https://doi.org/10.1007/s001900050271
https://link.springer.com/article/10.1007/s001900050271)"
);

}
// }}}

namespace fukushima_1999_customht_1
// {{{
{
#define USE_CUSTOM_HT

template <std::floating_point T>
auto f(const T t, const T u, const T v, const T w)
{
	// w * t**4 + u * t**3 + v * t - w
	return w * t * t * t * t + u * t * t * t + v * t - w;
}

template <std::floating_point T>
auto fp(const T t, const T u, const T v, const T w)
{
	// 4 * w * t**3 + 3 * u * t**2 + v
	return 4 * w * t * t * t + 3 * u * t * t + v;
}

template <std::floating_point T>
auto fpp(const T t, const T u, [[gnu::unused]] const T v, const T w)
{
	// 12 * w * t**2 + 6 * u * t
	return 12 * w * t * t + 6 * u * t;
}

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	double sin_lat = 0;
	double cos_lat = 0;
	constexpr auto c = ell.a * ell.e2;
	constexpr auto ep = 1 - ell.f;
	const auto zp = ep * std::abs(z);
	const auto u = 2 * (zp - c);
	const auto v = 2 * (zp + c);

	const auto tM = (c - zp) / w;

	double t = 0;

	if (tM <= 0)
	{
		// Case 1
		t = (w - c + zp) / (w - c + 2 * zp);
	}
	else if (tM >= 1)
	{
		// Case 2
		t = w / (zp + c);
	}
	else
	{
		auto fM = f(tM, u, v, w);

		// Case 3
		if (fM >= 0)
		{
			// Case 3a
			// (same as Case 2)
			t = w / (zp + c);
		}
		else // fM < 0
		{
			// Case 3b
			// (same as Case 1)
			t = (w - c + zp) / (w - c + 2 * zp);
		}
	}

	/*
	// SDW: This does not work
	sin_lat = z;
	cos_lat = w * (1 - ell.f);

	t = std::tan(sin_lat / cos_lat);
	*/

	// i = 1
	t -= f(t, u, v, w) / fp(t, u, v, w);

	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
#if 0
	// SDW: This does not work
	sin_lat = 2 * t;
	cos_lat = (1 - t * t) * ep;

	// SDW: This does not work
	sin_lat = 2 * t * ep;
	cos_lat = (1 - t * t);
#else
	// SDW: This does not work
	sin_lat = (1 - t * t) * ep;
	cos_lat = 2 * t;

	// SDW: This does not work
	//sin_lat = (1 - t * t);
	//cos_lat = 2 * t * ep;
#endif

	if (z < 0)
		sin_lat = -sin_lat;

	// 2 * tan(x/2) / (1 - tan(x/2)^2) == tan(x)
	// https://www.wolframalpha.com/input/?i=2+*+tan(x%2F2)+%2F+(1+-+tan(x%2F2)%5E2),+tan(x)

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = (2 * w * ep * t + z * (1 - t * t) - ell.a * ep * (1 + t * t)) / std::sqrt(SQ(1 + t * t) - 4 * ell.e2 * t * t);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 20;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "Toshio Fukushima",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226311152_Fast_transform_from_geocentric_to_geodetic_coordinates",
	/*.citation                    =*/ R"(Fukushima, T. Journal of Geodesy (1999) 73: 603. https://doi.org/10.1007/s001900050271
https://link.springer.com/article/10.1007/s001900050271)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace fukushima_2006_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto z_c = z * (1 - ell.f);
	const auto w_c = w * (1 - ell.f);
	const auto c = ell.a * ell.e2;

	double S_n;
	double C_n;

	double A_n;
	double B_n; // correction factor of Halley's method

	S_n = z;
	C_n = w_c;

	for (int i = 1; i <= max_iterations; ++i)
	{
		A_n = hypot(S_n, C_n);
		B_n = 1.5 * c * S_n * C_n * ((w * S_n - z_c * C_n) * A_n - c * S_n * C_n);

		// SDW: at great heights (e.g. 8000 km), there is overflow
		S_n = (z_c * CB(A_n) + c * CB(S_n)) * CB(A_n) - B_n * S_n;
		C_n = (w * CB(A_n) - c * CB(C_n)) * CB(A_n) - B_n * C_n;
	}

	auto sin_lat = S_n;
	auto cos_lat = C_n;
	cos_lat *= (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -5,
	/*.algo_author                 =*/ "Toshio Fukushima",
	/*.code_copyright              =*/ "Toshio Fukushima",
	/*.license                     =*/ "Unknown",
	/*.orig_impl_lang              =*/ "Fortran",
	/*.url                         =*/ "https://www.researchgate.net/publication/227215135_Transformation_from_Cartesian_to_Geodetic_Coordinates_Accelerated_by_Halley%27s_Method",
	/*.citation                    =*/ R"(Fukushima, Toshio. (2006). Transformation from Cartesian to Geodetic Coordinates Accelerated by Halley’s Method. Journal of Geodesy. 79. 689-693. 10.1007/s00190-006-0023-2.)"
);

}
// }}}

namespace fukushima_2006_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto z_c = z * (1 - ell.f);
	const auto w_c = w * (1 - ell.f);
	const auto c = ell.a * ell.e2;

	double S_n;
	double C_n;

	double A_n;
	double B_n; // correction factor of Halley's method

	S_n = z;
	C_n = w_c;

	for (int i = 1; i <= max_iterations; ++i)
	{
		A_n = hypot(S_n, C_n);
		B_n = 1.5 * c * S_n * C_n * ((w * S_n - z_c * C_n) * A_n - c * S_n * C_n);

		// SDW: at great heights (e.g. 8000 km), there is overflow
		S_n = (z_c * CB(A_n) + c * CB(S_n)) * CB(A_n) - B_n * S_n;
		C_n = (w * CB(A_n) - c * CB(C_n)) * CB(A_n) - B_n * C_n;
	}

	auto sin_lat = S_n;
	auto cos_lat = C_n;
	cos_lat *= (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "Toshio Fukushima",
	/*.code_copyright              =*/ "Toshio Fukushima",
	/*.license                     =*/ "Unknown",
	/*.orig_impl_lang              =*/ "Fortran",
	/*.url                         =*/ "https://www.researchgate.net/publication/227215135_Transformation_from_Cartesian_to_Geodetic_Coordinates_Accelerated_by_Halley%27s_Method",
	/*.citation                    =*/ R"(Fukushima, Toshio. (2006). Transformation from Cartesian to Geodetic Coordinates Accelerated by Halley’s Method. Journal of Geodesy. 79. 689-693. 10.1007/s00190-006-0023-2.)"
);

}
// }}}

namespace geographiclib
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	//const auto w = hypot(x, y);
	//const double sin_lambda = w != 0 ? y / w : 0;
	//const double cos_lambda = w != 0 ? x / w : 1;
	//ht = std::sqrt(w2 + z2);      // Distance to center of earth
	double sin_lat;
	double cos_lat;
	// Treat prolate spheroids by swapping w and z here and by switching
	// the arguments to phi = atan2(...) at the end.
	const auto p = w2 / ell.a2;
	// (1 - ell.e2) / ell.a2 == 1 / (Rp * Rp)
	const auto q = (z2 / ell.a2) * (1 - ell.e2);
	const auto r = (p + q - e4) / 6;
	if (!(e4 * q == 0 && r <= 0))
	{
		// Avoid possible division by zero when r = 0 by multiplying
		// equations for s and t by r^3 and r, resp.
		const auto S = e4 * p * q / 4; // S = r^3 * s
		const auto r3 = r * r * r;
		const auto disc = S * (2 * r3 + S);
		auto u = r;
		if (disc >= 0)
		{
			const auto r2 = r * r;

			auto T3 = S + r3;
			// Pick the sign on the sqrt to maximize abs(T3).  This minimizes
			// loss of precision due to cancellation.  The result is unchanged
			// because of the way the T_ is used in definition of u.
			T3 += T3 < 0 ? -std::sqrt(disc) : std::sqrt(disc); // T3 = (r * t)^3
			// N.B. cbrt always returns the real root.  cbrt(-8) = -2.
			const auto T_ = std::cbrt(T3); // T_ = r * t
			// T_ can be zero; but then r2 / T_ -> 0.
			u += T_ + (T_ != 0 ? r2 / T_ : 0);
		}
		else
		{
			// T_ is complex, but the way u is defined the result is real.
			const auto ang = std::atan2(std::sqrt(-disc), -(S + r3));
			// There are three possible cube roots.  We choose the root which
			// avoids cancellation.  Note that disc < 0 implies that r < 0.
			u += 2 * r * std::cos(ang / 3);
		}
		const auto v = std::sqrt(u * u + e4 * q); // guaranteed positive
		// Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
		// e2*e2 * q / (v - u) because u ~ e^4 when q is small and u < 0.
		const auto uv = u < 0 ? e4 * q / (v - u) : u + v; // u+v, guaranteed positive
		// Need to guard against w_ going negative due to roundoff in uv - q.
		const auto w_ = std::max(0.0, 0.5 * ell.e2 * (uv - q) / v);
		// Rearrange expression for k to avoid loss of accuracy due to
		// subtraction.  Division by 0 not possible because uv > 0, w_ >= 0.
		const auto k = uv / (std::sqrt(uv + w_ * w_) + w_);
		const auto k1 = ell.f >= 0 ? k : k - ell.e2;
		const auto k2 = ell.f >= 0 ? k + ell.e2 : k;
		const auto H = hypot(z/k1, w/k2);
		sin_lat = (z/k1) / H;
		cos_lat = (w/k2) / H;
#ifdef USE_CUSTOM_HT
		const auto d = k1 * w / k2;
		ht = (1 - (1 - ell.e2)/k1) * hypot(d, z);
#endif
	}
	else
	{
		// e2*e2 * q == 0 && r <= 0
		// This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
		// (prolate, rotation axis) and the generation of 0/0 in the general
		// formulas for phi and ht.  using the general formula and division by 0
		// in formula for ht.  So handle this case by taking the limits:
		// f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e2*e2 - p)
		// f < 0: w -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e2*e2 - p)
		const auto zz = std::sqrt((ell.f >= 0 ? e4 - p : p) / (1 - ell.e2));
		const auto xx = std::sqrt( ell.f <  0 ? e4 - p : p);
		const auto H = hypot(zz, xx);
		sin_lat = zz / H;
		cos_lat = xx / H;
		if (z < 0)
			sin_lat = -sin_lat; // for tiny negative z (not for prolate)

#ifdef USE_CUSTOM_HT
		ht = -ell.a * (ell.f >= 0 ? (1 - ell.e2) : 1) * H / ell.e2;
		// (1-e2)/e2 == 1/ep2
#endif
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	//lon_rad = std::atan2(sin_lambda, cos_lambda);

#ifdef USE_CUSTOM_HT
#else
	// sin and cos are sufficiently accurate (and normalized) to be used in the ell.get_ht function, but this is more accurate
	ht = ell.get_ht(w, z, lat_rad);
	// XXX: is this better?
	//normalize(cos_lat, sin_lat);
	//ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Charles F. F. Karney",
	/*.code_copyright              =*/ "Charles Karney",
	/*.license                     =*/ "MIT/X11",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "https://sourceforge.net/p/geographiclib/code/ci/release/tree/src/Geocentric.cpp",
	/*.citation                    =*/ R"(Geodesics on an ellipsoid of revolution
Charles F. F. Karney
https://arxiv.org/abs/1102.1215)"
);

}
// }}}

namespace geographiclib_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	//const auto w = hypot(x, y);
	//const double sin_lambda = w != 0 ? y / w : 0;
	//const double cos_lambda = w != 0 ? x / w : 1;
	//ht = std::sqrt(w2 + z2);      // Distance to center of earth
	double sin_lat;
	double cos_lat;
	// Treat prolate spheroids by swapping w and z here and by switching
	// the arguments to phi = atan2(...) at the end.
	const auto p = w2 / ell.a2;
	// (1 - ell.e2) / ell.a2 == 1 / (Rp * Rp)
	const auto q = (z2 / ell.a2) * (1 - ell.e2);
	const auto r = (p + q - e4) / 6;
	if (!(e4 * q == 0 && r <= 0))
	{
		// Avoid possible division by zero when r = 0 by multiplying
		// equations for s and t by r^3 and r, resp.
		const auto S = e4 * p * q / 4; // S = r^3 * s
		const auto r3 = r * r * r;
		const auto disc = S * (2 * r3 + S);
		auto u = r;
		if (disc >= 0)
		{
			const auto r2 = r * r;

			auto T3 = S + r3;
			// Pick the sign on the sqrt to maximize abs(T3).  This minimizes
			// loss of precision due to cancellation.  The result is unchanged
			// because of the way the T_ is used in definition of u.
			T3 += T3 < 0 ? -std::sqrt(disc) : std::sqrt(disc); // T3 = (r * t)^3
			// N.B. cbrt always returns the real root.  cbrt(-8) = -2.
			const auto T_ = std::cbrt(T3); // T_ = r * t
			// T_ can be zero; but then r2 / T_ -> 0.
			u += T_ + (T_ != 0 ? r2 / T_ : 0);
		}
		else
		{
			// T_ is complex, but the way u is defined the result is real.
			const auto ang = std::atan2(std::sqrt(-disc), -(S + r3));
			// There are three possible cube roots.  We choose the root which
			// avoids cancellation.  Note that disc < 0 implies that r < 0.
			u += 2 * r * std::cos(ang / 3);
		}
		const auto v = std::sqrt(u * u + e4 * q); // guaranteed positive
		// Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
		// e2*e2 * q / (v - u) because u ~ e^4 when q is small and u < 0.
		const auto uv = u < 0 ? e4 * q / (v - u) : u + v; // u+v, guaranteed positive
		// Need to guard against w_ going negative due to roundoff in uv - q.
		const auto w_ = std::max(0.0, 0.5 * ell.e2 * (uv - q) / v);
		// Rearrange expression for k to avoid loss of accuracy due to
		// subtraction.  Division by 0 not possible because uv > 0, w_ >= 0.
		const auto k = uv / (std::sqrt(uv + w_ * w_) + w_);
		const auto k1 = ell.f >= 0 ? k : k - ell.e2;
		const auto k2 = ell.f >= 0 ? k + ell.e2 : k;
		const auto H = hypot(z/k1, w/k2);
		sin_lat = (z/k1) / H;
		cos_lat = (w/k2) / H;
#ifdef USE_CUSTOM_HT
		const auto d = k1 * w / k2;
		ht = (1 - (1 - ell.e2)/k1) * hypot(d, z);
#endif
	}
	else
	{
		// e2*e2 * q == 0 && r <= 0
		// This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
		// (prolate, rotation axis) and the generation of 0/0 in the general
		// formulas for phi and ht.  using the general formula and division by 0
		// in formula for ht.  So handle this case by taking the limits:
		// f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e2*e2 - p)
		// f < 0: w -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e2*e2 - p)
		const auto zz = std::sqrt((ell.f >= 0 ? e4 - p : p) / (1 - ell.e2));
		const auto xx = std::sqrt( ell.f <  0 ? e4 - p : p);
		const auto H = hypot(zz, xx);
		sin_lat = zz / H;
		cos_lat = xx / H;
		if (z < 0)
			sin_lat = -sin_lat; // for tiny negative z (not for prolate)

#ifdef USE_CUSTOM_HT
		ht = -ell.a * (ell.f >= 0 ? (1 - ell.e2) : 1) * H / ell.e2;
		// (1-e2)/e2 == 1/ep2
#endif
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	//lon_rad = std::atan2(sin_lambda, cos_lambda);

#ifdef USE_CUSTOM_HT
#else
	// sin and cos are sufficiently accurate (and normalized) to be used in the ell.get_ht function, but this is more accurate
	ht = ell.get_ht(w, z, lat_rad);
	// XXX: is this better?
	//normalize(cos_lat, sin_lat);
	//ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Charles F. F. Karney",
	/*.code_copyright              =*/ "Charles Karney",
	/*.license                     =*/ "MIT/X11",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "https://sourceforge.net/p/geographiclib/code/ci/release/tree/src/Geocentric.cpp",
	/*.citation                    =*/ R"(Geodesics on an ellipsoid of revolution
Charles F. F. Karney
https://arxiv.org/abs/1102.1215)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace geotransformCpp
// {{{
{

constexpr int _line_begin = __LINE__;
//void Gcc_To_Gdc_Converter::Convert(int count, const Gcc_Coord_3d gcc[], Gdc_Coord_3d gdc[] )
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	// UPPER BOUNDS ON POINT

	//const auto ARat1  = std::pow((ell.a + 50005), 2);
	constexpr auto ARat1  = SQ(ell.a + 50005);
	//const auto ARat2  = (ARat1) / std::pow((ell.b + 50005), 2);
	constexpr auto ARat2  = (ARat1) / SQ(ell.b + 50005);

	// LOWER BOUNDS ON POINT

	//const auto BRat1  = std::pow((ell.a - 10005), 2);
	constexpr auto BRat1  = SQ(ell.a - 10005);
	//const auto BRat2  = (BRat1) / std::pow((ell.b - 10005), 2);
	constexpr auto BRat2  = (BRat1) / SQ(ell.b - 10005);

	constexpr auto B1 = 0.100225438677758E+01;
	constexpr auto B2 = -0.393246903633930E-04;
	constexpr auto B3 = 0.241216653453483E+12;
	constexpr auto B4 = 0.133733602228679E+14;
	constexpr auto B5 = 0.984537701867943E+00;

	[[gnu::unused]] double /*w2,w,z2,*/testu,testb,top,top2,rr,q,s12,rnn,s1,/*zp2,wp,wp2,*/cf,gee,alpha,cl,arg2,p,xarg,r2,r1,ro,
		   s,roe,arg,v,zo;

	//for (int i=0; i < count; i++)
	//{

	/* CHECK FOR SPECIAL CASES*/

	if (!(x == 0)) {;} /* null statement */
	else
	{
		if (y > 0)
		{
			lat_rad = M_PI_2;
		}
		else
		{
			if (y < 0)
			{
				//lon_rad = -M_PI_2;
			}
			else
			{
				if (z > 0)
				{
					lat_rad = M_PI_2;
					//lon_rad = 0;
					//ht = z;
					ht = z - ell.b;
					return;
				}
				else
				{
					if (z < 0)
					{
						lat_rad = -M_PI_2;
						//lon_rad = 0;
						//ht = z;
						ht = -(z + ell.b);
						return;
					}
					else
					{
						lat_rad = 0;
						//lon_rad = 0;
						ht = 0;
						return;
					}
				}
			}
		}
	}

	/* END OF SPECIAL CASES */

	testu = w2 + ARat2 * z2;
	testb = w2 + BRat2 * z2;

	if ((testb > BRat1) && (testu < ARat1))
	{
		/*POINT IS BETWEEN-10 KIL AND 50 KIL, SO COMPUTE TANGENT LATITUDE */

		top = z * (B1 + (B2 * w2 + B3) /
				(B4 + w2 * B5 + z2));

		top2 = top * top;

		rr = top2 + w2;

		q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
		/* ****************************************************************

		COMPUTE H IN LINE SQUARE ROOT OF 1-ell.e2*SIN*SIN.  USE SHORT BINOMIAL
		EXPANSION PLUS ONE ITERATION OF NEWTON'S METHOD FOR SQUARE ROOTS.
		*/

		s12 = top2 / rr;

		rnn = ell.a / ((0.25 - 0.25 * ell.e2 * s12 + .9999944354799 / 4) + (0.25 - 0.25 * ell.e2 * s12) / (0.25 - 0.25 * ell.e2 * s12 + .9999944354799 / 4));
		s1 = top / q;

		/******************************************************************/

		/* TEST FOR H NEAR POLE.  if SIN(_)**2 <= SIN(45.)**2 THEN NOT NEAR A POLE.*/

		if (s12 < 0.5)
		{
			ht = q - rnn;
		}
		else
		{
			//ht = z / s1 + ((ell.e2 - 1) * rnn);
			ht = z / s1 - (1 - ell.e2) * rnn;
		}
#endif
		//lat_rad = std::atan(top / w);
		lat_rad = std::atan2(top, w);
		//lon_rad = std::atan2(y, x);
	}
	/* POINT ABOVE 50 KILOMETERS OR BELOW -10 KILOMETERS  */
	else /* Do Exact Solution  ************ */
	{
		cf = 54 * ell.b2 * z2;

		//gee = w2 - ((ell.e2 - 1) * z2) - ell.e2 * (ell.a2 - ell.b2);
		gee = w2 + (1 - ell.e2) * z2 - ell.e2 * (ell.a2 - ell.b2);

		alpha = cf / (gee * gee);

		cl = e4 * w2 * alpha / gee;

		arg2 = cl * (cl + 2);

		s1 = 1 + cl + std::sqrt(arg2);

		//s = std::pow(s1, (1.0 / 3));
		s = std::cbrt(s1);

		//p = alpha / (3 * std::pow((s + (1 / s) + 1), 2));
		p = alpha / (3 * SQ(s + (1 / s) + 1));

		xarg = 1 + (2 * e4 * p);

		q = std::sqrt(xarg);

		r2= -p * (2 * (1 - ell.e2) * z2 / (q * (1 + q)) + w2);

		r1 = (1 + (1 / q));

		r2 /= ell.a2;

		/* DUE TO PRECISION ERRORS THE ARGUMENT MAY BECOME NEGATIVE IF SO SET THE ARGUMENT TO ZERO.*/

		if (r1 + r2 > 0)
		{
			ro = ell.a * std::sqrt(0.5 * (r1 + r2));
		}
		else
		{
			ro = 0;
		}

		ro = ro - p * ell.e2 * w / (1 + q);

		// unused expression arg0 = std::pow((w - ell.e2 * ro),2) + z2;

		roe = ell.e2 * ro;
		//arg = std::pow((w - roe), 2) + z2;
		arg = SQ(w - roe) + z2;
		v = std::sqrt(arg - ell.e2 * z2);

		zo = (ell.b2 / ell.a) * z / v;

#ifdef USE_CUSTOM_HT
		// b2/a == b*(1-f)
		ht = std::sqrt(arg) * (1 - (ell.b2 / ell.a) / v);
#endif

		// (a2-b2)/b2 == ep2
		//top = z + ((ell.a2 - ell.b2) / ell.b2) * zo;
		top = z + ell.ep2 * zo;

		//lat_rad = std::atan(top / w);
		lat_rad = std::atan2(top, w);
		//lon_rad = std::atan2(y, x);
	}  /* end of Exact solution */

	//} // end for

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "R. Toms",
	/*.code_copyright              =*/ "SEDRIS",
	/*.license                     =*/ "\"This software was developed for use by the United States Government with unlimited rights.\"",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "http://www.ai.sri.com/geotransform/download.shtml",
	/*.citation                    =*/ R"(geotransformCpp 1.5
file:Gcc_To_Gdc_Converter.cpp)"
);

}
// }}}

namespace geotransformCpp_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
//void Gcc_To_Gdc_Converter::Convert(int count, const Gcc_Coord_3d gcc[], Gdc_Coord_3d gdc[] )
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	// UPPER BOUNDS ON POINT

	//const auto ARat1  = std::pow((ell.a + 50005), 2);
	constexpr auto ARat1  = SQ(ell.a + 50005);
	//const auto ARat2  = (ARat1) / std::pow((ell.b + 50005), 2);
	constexpr auto ARat2  = (ARat1) / SQ(ell.b + 50005);

	// LOWER BOUNDS ON POINT

	//const auto BRat1  = std::pow((ell.a - 10005), 2);
	constexpr auto BRat1  = SQ(ell.a - 10005);
	//const auto BRat2  = (BRat1) / std::pow((ell.b - 10005), 2);
	constexpr auto BRat2  = (BRat1) / SQ(ell.b - 10005);

	constexpr auto B1 = 0.100225438677758E+01;
	constexpr auto B2 = -0.393246903633930E-04;
	constexpr auto B3 = 0.241216653453483E+12;
	constexpr auto B4 = 0.133733602228679E+14;
	constexpr auto B5 = 0.984537701867943E+00;

	[[gnu::unused]] double /*w2,w,z2,*/testu,testb,top,top2,rr,q,s12,rnn,s1,/*zp2,wp,wp2,*/cf,gee,alpha,cl,arg2,p,xarg,r2,r1,ro,
		   s,roe,arg,v,zo;

	//for (int i=0; i < count; i++)
	//{

	/* CHECK FOR SPECIAL CASES*/

	if (!(x == 0)) {;} /* null statement */
	else
	{
		if (y > 0)
		{
			lat_rad = M_PI_2;
		}
		else
		{
			if (y < 0)
			{
				//lon_rad = -M_PI_2;
			}
			else
			{
				if (z > 0)
				{
					lat_rad = M_PI_2;
					//lon_rad = 0;
					//ht = z;
					ht = z - ell.b;
					return;
				}
				else
				{
					if (z < 0)
					{
						lat_rad = -M_PI_2;
						//lon_rad = 0;
						//ht = z;
						ht = -(z + ell.b);
						return;
					}
					else
					{
						lat_rad = 0;
						//lon_rad = 0;
						ht = 0;
						return;
					}
				}
			}
		}
	}

	/* END OF SPECIAL CASES */

	testu = w2 + ARat2 * z2;
	testb = w2 + BRat2 * z2;

	if ((testb > BRat1) && (testu < ARat1))
	{
		/*POINT IS BETWEEN-10 KIL AND 50 KIL, SO COMPUTE TANGENT LATITUDE */

		top = z * (B1 + (B2 * w2 + B3) /
				(B4 + w2 * B5 + z2));

		top2 = top * top;

		rr = top2 + w2;

		q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
		/* ****************************************************************

		COMPUTE H IN LINE SQUARE ROOT OF 1-ell.e2*SIN*SIN.  USE SHORT BINOMIAL
		EXPANSION PLUS ONE ITERATION OF NEWTON'S METHOD FOR SQUARE ROOTS.
		*/

		s12 = top2 / rr;

		rnn = ell.a / ((0.25 - 0.25 * ell.e2 * s12 + .9999944354799 / 4) + (0.25 - 0.25 * ell.e2 * s12) / (0.25 - 0.25 * ell.e2 * s12 + .9999944354799 / 4));
		s1 = top / q;

		/******************************************************************/

		/* TEST FOR H NEAR POLE.  if SIN(_)**2 <= SIN(45.)**2 THEN NOT NEAR A POLE.*/

		if (s12 < 0.5)
		{
			ht = q - rnn;
		}
		else
		{
			//ht = z / s1 + ((ell.e2 - 1) * rnn);
			ht = z / s1 - (1 - ell.e2) * rnn;
		}
#endif
		//lat_rad = std::atan(top / w);
		lat_rad = std::atan2(top, w);
		//lon_rad = std::atan2(y, x);
	}
	/* POINT ABOVE 50 KILOMETERS OR BELOW -10 KILOMETERS  */
	else /* Do Exact Solution  ************ */
	{
		cf = 54 * ell.b2 * z2;

		//gee = w2 - ((ell.e2 - 1) * z2) - ell.e2 * (ell.a2 - ell.b2);
		gee = w2 + (1 - ell.e2) * z2 - ell.e2 * (ell.a2 - ell.b2);

		alpha = cf / (gee * gee);

		cl = e4 * w2 * alpha / gee;

		arg2 = cl * (cl + 2);

		s1 = 1 + cl + std::sqrt(arg2);

		//s = std::pow(s1, (1.0 / 3));
		s = std::cbrt(s1);

		//p = alpha / (3 * std::pow((s + (1 / s) + 1), 2));
		p = alpha / (3 * SQ(s + (1 / s) + 1));

		xarg = 1 + (2 * e4 * p);

		q = std::sqrt(xarg);

		r2= -p * (2 * (1 - ell.e2) * z2 / (q * (1 + q)) + w2);

		r1 = (1 + (1 / q));

		r2 /= ell.a2;

		/* DUE TO PRECISION ERRORS THE ARGUMENT MAY BECOME NEGATIVE IF SO SET THE ARGUMENT TO ZERO.*/

		if (r1 + r2 > 0)
		{
			ro = ell.a * std::sqrt(0.5 * (r1 + r2));
		}
		else
		{
			ro = 0;
		}

		ro = ro - p * ell.e2 * w / (1 + q);

		// unused expression arg0 = std::pow((w - ell.e2 * ro),2) + z2;

		roe = ell.e2 * ro;
		//arg = std::pow((w - roe), 2) + z2;
		arg = SQ(w - roe) + z2;
		v = std::sqrt(arg - ell.e2 * z2);

		zo = (ell.b2 / ell.a) * z / v;

#ifdef USE_CUSTOM_HT
		// b2/a == b*(1-f)
		ht = std::sqrt(arg) * (1 - (ell.b2 / ell.a) / v);
#endif

		// (a2-b2)/b2 == ep2
		//top = z + ((ell.a2 - ell.b2) / ell.b2) * zo;
		top = z + ell.ep2 * zo;

		//lat_rad = std::atan(top / w);
		lat_rad = std::atan2(top, w);
		//lon_rad = std::atan2(y, x);
	}  /* end of Exact solution */

	//} // end for

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "R. Toms",
	/*.code_copyright              =*/ "SEDRIS",
	/*.license                     =*/ "\"This software was developed for use by the United States Government with unlimited rights.\"",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "http://www.ai.sri.com/geotransform/download.shtml",
	/*.citation                    =*/ R"(geotransformCpp 1.5
file:Gcc_To_Gdc_Converter.cpp)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace gersten_1961
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);
	const auto e_ = ell.a * ell.e2 / r;
	const auto c = w / r;
	const auto s = z / r;

	const auto cos_lat_lat_ = 1 - 0.5 * e_ * e_ * s * s * c * c;
	const auto sin_lat_lat_ = s * c * (1 + e_ - (2 * e_ - 0.5 * ell.e2) * s * s);

	const auto s2 = s * cos_lat_lat_ + c * sin_lat_lat_;
	const auto c2 = c * cos_lat_lat_ - s * sin_lat_lat_;

	// SDW: taking the sqrt doesn't help
	auto sin_lat = s2;
	auto cos_lat = c2;

	if (z < 0)
		sin_lat = -sin_lat;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
/*
#ifdef USE_CUSTOM_HT
	ht = r - ell.a + 0.5 * ell.a * ell.e2 * s * s * (1 + e_ - (0.25 * ell.e2 - e_) * s * s);
#else
#endif
*/
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "R. H. Gersten",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "None",
	/*.citation                    =*/ R"(R. H. Gersten, "Geodetic Sub-Latitude and Altitude of a Space Vehicle", Journal of the Astronautical Sciences, Vol. 3, No. 1, Spring 1961, pp. 28-29)"
);

}
// }}}

namespace halley_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= halley_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_halley_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HalleysMethod.html",
	/*.citation                    =*/ R"(Halley's Method: Δφ = (f * f') / (f' * f' - 0.5 * f * f''))"
);

}
// }}}

namespace halley_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= halley_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_halley_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -10,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HalleysMethod.html",
	/*.citation                    =*/ R"(Halley's Method: Δφ = (f * f') / (f' * f' - 0.5 * f * f''))"
);

}
// }}}

namespace halley_quick_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = halley_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_halley_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HalleysMethod.html",
	/*.citation                    =*/ R"(Halley's Method: Δφ = (f * f') / (f' * f' - 0.5 * f * f''))"
);

}
// }}}

namespace halley_quick_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = halley_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_halley_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HalleysMethod.html",
	/*.citation                    =*/ R"(Halley's Method: Δφ = (f * f') / (f' * f' - 0.5 * f * f''))"
);

}
// }}}

namespace heikkinen_1982
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto F = 54 * ell.b2 * z2;
	//const auto G = w2 + (1 - ell.e2) * z2 - ell.a2 * ell.e2 * ell.e2;
	// (a2-b2)*e2 == a2*e4
	const auto G = w2 + (1 - ell.e2) * z2 - ell.a2 * e4;
	const auto C = e4 * F * w2 / CB(G);
	const auto S = std::cbrt(1 + C + std::sqrt(C * (C + 2)));
	const auto tmp1 = S + 1 / S + 1;
	const auto P = F / (3 * tmp1 * tmp1 * G * G);
	const auto Q = std::sqrt(1 + 2 * e4 * P);

	const auto r0_1 = -P * ell.e2 * w / (1 + Q);
	const auto r0_2 = 0.5 * ell.a2 * (1 + 1 / Q);
	const auto r0_3 = (P * (1 - ell.e2) * z2) / (Q * (1 + Q));
	const auto r0_4 = P * w2 / 2;

	// std::abs needed for large heights
	const auto r0 = r0_1 + std::sqrt(std::abs(r0_2 - r0_3 - r0_4));
	const auto V = hypot(w - ell.e2 * r0, (1 - ell.f) * z);
	//const auto z0 = b2 * z / (a * V);
	// b2/a == a*(1-e2)
	//const auto z0 = a * (1-e2) * (z / V);

	//A = z + ell.ep2 * z0;
	//A = z + ell.ep2 * ell.a * (1 - ell.e2) * (z / V);
	// ep2 == e2/(1-e2)
	auto sin_lat = z * (1 + ell.a * ell.e2 / V);
	// SDW: this is more accurate
	//auto sin_lat = z + ell.a * ell.e2 * z / V;
	auto cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	const auto U = hypot(w - ell.e2 * r0, z);
	//ht = U * (1 - ell.b2 / (ell.a * V));
	// b2/a == (1-e2)*a
	// SDW: this is more accurate
	ht = U * (1 - (1 - ell.e2) * ell.a / V);
	// SDW: accuracy is comparable to above, but no obvious benefit
	//ht = U - U * (1 - ell.e2) * (ell.a / V);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "M. Heikkinen",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "None",
	/*.citation                    =*/ R"(Heikkinen, M. "Closed formulas for the calculation of spatial geodetic coordinates from rectangular coordinates." Journal of Surveying 107.1982 (1982): 207-211.)"
);

}
// }}}

namespace heikkinen_1982_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto F = 54 * ell.b2 * z2;
	//const auto G = w2 + (1 - ell.e2) * z2 - ell.a2 * ell.e2 * ell.e2;
	// (a2-b2)*e2 == a2*e4
	const auto G = w2 + (1 - ell.e2) * z2 - ell.a2 * e4;
	const auto C = e4 * F * w2 / CB(G);
	const auto S = std::cbrt(1 + C + std::sqrt(C * (C + 2)));
	const auto tmp1 = S + 1 / S + 1;
	const auto P = F / (3 * tmp1 * tmp1 * G * G);
	const auto Q = std::sqrt(1 + 2 * e4 * P);

	const auto r0_1 = -P * ell.e2 * w / (1 + Q);
	const auto r0_2 = 0.5 * ell.a2 * (1 + 1 / Q);
	const auto r0_3 = (P * (1 - ell.e2) * z2) / (Q * (1 + Q));
	const auto r0_4 = P * w2 / 2;

	// std::abs needed for large heights
	const auto r0 = r0_1 + std::sqrt(std::abs(r0_2 - r0_3 - r0_4));
	const auto V = hypot(w - ell.e2 * r0, (1 - ell.f) * z);
	//const auto z0 = b2 * z / (a * V);
	// b2/a == a*(1-e2)
	//const auto z0 = a * (1-e2) * (z / V);

	//A = z + ell.ep2 * z0;
	//A = z + ell.ep2 * ell.a * (1 - ell.e2) * (z / V);
	// ep2 == e2/(1-e2)
	auto sin_lat = z * (1 + ell.a * ell.e2 / V);
	// SDW: this is more accurate
	//auto sin_lat = z + ell.a * ell.e2 * z / V;
	auto cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	const auto U = hypot(w - ell.e2 * r0, z);
	//ht = U * (1 - ell.b2 / (ell.a * V));
	// b2/a == (1-e2)*a
	// SDW: this is more accurate
	ht = U * (1 - (1 - ell.e2) * ell.a / V);
	// SDW: accuracy is comparable to above, but no obvious benefit
	//ht = U - U * (1 - ell.e2) * (ell.a / V);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "M. Heikkinen",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "None",
	/*.citation                    =*/ R"(Heikkinen, M. "Closed formulas for the calculation of spatial geodetic coordinates from rectangular coordinates." Journal of Surveying 107.1982 (1982): 207-211.)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace householder_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= householder_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_householder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HouseholdersMethod.html",
	/*.citation                    =*/ R"(Householder's Method: Δφ = (f / f') * (1 + 0.5 * f * f'' / (f' * f')))"
);

}
// }}}

namespace householder_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= householder_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_householder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -10,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HouseholdersMethod.html",
	/*.citation                    =*/ R"(Householder's Method: Δφ = (f / f') * (1 + 0.5 * f * f'' / (f' * f')))"
);

}
// }}}

namespace householder_quick_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = householder_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_householder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HouseholdersMethod.html",
	/*.citation                    =*/ R"(Householder's Method: Δφ = (f / f') * (1 + 0.5 * f * f'' / (f' * f')))"
);

}
// }}}

namespace householder_quick_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = householder_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_householder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/HouseholdersMethod.html",
	/*.citation                    =*/ R"(Householder's Method: Δφ = (f / f') * (1 + 0.5 * f * f'' / (f' * f')))"
);

}
// }}}

namespace jat_spacetime_geodetic
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

#if 0
	// SDW: This threshold is too large to be accurate.
	constexpr double MACHEPS = 1.1920929E-7;
	constexpr auto eps = 1000.0 * MACHEPS; // Convergence criterion
	constexpr auto eps_a = eps * ell.a;
#else
	constexpr double eps_a = 1E-4;
#endif

	ht = -ell.a;

	double sin_lat = 0;
	double ZdZ = 0;
	double Nh = 0;
	double Rn = 0;
	auto dZ = ell.e2 * z;
	double dZ_new = 0;

	while (std::abs(dZ - dZ_new) > eps_a)
	{
		ZdZ = z + dZ;
		Nh = hypot(w, ZdZ);
		sin_lat = ZdZ / Nh;
		Rn = ell.get_Rn(sin_lat);
		dZ = dZ_new;
		dZ_new = Rn * ell.e2 * sin_lat;
	}

	lat_rad = std::atan2(ZdZ, w);

	// SDW: This is not accurate.
	//ht = Nh - Rn;

	ht = ell.get_ht(w, z, lat_rad);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -3,
	/*.algo_author                 =*/ "O. Montenbruck, E. Gill",
	/*.code_copyright              =*/ "NASA",
	/*.license                     =*/ "NASA Open Source Agreement",
	/*.orig_impl_lang              =*/ "Java",
	/*.url                         =*/ "https://sourceforge.net/p/jat/code/HEAD/tree/trunk/src/jat/coreNOSA/timeRef/Geodetic.java",
	/*.citation                    =*/ R"(Java Astrodynamics Toolkit (JAT))" // citation
);

}
// }}}

namespace jones_2002_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto a2 = ell.a2;
	const auto b2 = ell.b2;
	const auto a_e2 = ell.a * ell.e2;

	double sin_lat;
	double cos_lat;
	//double u0;

	// geocentric to parametric

	if (w2 / a2 + z2 / b2 >= 1)
	{
		// region 1: the point lies on or outside the ellipse

		sin_lat = z;
		cos_lat = w * (1 - ell.f);
	}
	else // w2 / a2 + z2 / b2 < 1
	{
		if (w <= a_e2 + z / (1 - ell.f))
		{
			// region 2: the point lies inside the ellipse with high latitude and passes through the 2 vertices of the evolute

			sin_lat = z * (1 - ell.f) + a_e2;
			cos_lat = w;
		}
		else // w > a_e2 + z / (1 - ell.f)
		{
			// region 3: the point lies inside the ellipse with low latitude and passes through the 2 vertices of the evolute

			sin_lat = z * (1 - ell.f);
			cos_lat = w - a_e2;
		}
	}

	//u0 = std::atan2(sin_lat, cos_lat);

	// SDW: what is F'(u) ?

	// i = 1
	normalize(cos_lat, sin_lat);

	// SDW: this probably isn't right

	sin_lat = z * (1 - ell.f) + a_e2 * sin_lat;
	cos_lat = w;

	// parametric to geodetic
	cos_lat *= (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 1,
	/*.algo_author                 =*/ "G. C. Jones",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/225462506_New_solutions_for_the_geodetic_coordinate_transformation",
	/*.citation                    =*/ R"(C. Jones, G. (2002). New solutions for the geodetic coordinate transformation. Journal of Geodesy. 76. 437-446. 10.1007/s00190-002-0267-4.)"
);

}
// }}}

namespace ligas_2011_I_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double J[2][2];
	double J_inv[2][2];
	double f_X[2];
	double X_delta[2];

	auto we = w * ell.a / r;
	auto ze = z * ell.b / r;

	for (int i = 1; i <= max_iterations; ++i)
	{
		f_X[0] = ligas_f1(w, we, z, ze);
		f_X[1] = ligas_f2(we, ze);
		ligas_Jacobian(w, we, z, ze, J);
		inv(J, J_inv);
		mul(J_inv, f_X, X_delta);
		we -= X_delta[0];
		ze -= X_delta[1];
	}

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	// SDW: This is not accurate
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_ligas_util;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "M. Ligas, P. Banasik",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/258491472_Conversion_between_Cartesian_and_geodetic_coordinates_on_a_rotational_ellipsoid_by_solving_a_system_of_nonlinear_equations",
	/*.citation                    =*/ R"(Ligas, Marcin & Banasik, Piotr. (2011). Conversion between Cartesian and geodetic coordinates on a rotational ellipsoid by solving a system of nonlinear equations. Geodesy and Cartography. 60. 145-159. 10.2478/v10277-012-0013-x.)"
);

}
// }}}

namespace ligas_2011_I_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double J[2][2];
	double J_inv[2][2];
	double f_X[2];
	double X_delta[2];

	auto we = w * ell.a / r;
	auto ze = z * ell.b / r;

	for (int i = 1; i <= max_iterations; ++i)
	{
		f_X[0] = ligas_f1(w, we, z, ze);
		f_X[1] = ligas_f2(we, ze);
		ligas_Jacobian(w, we, z, ze, J);
		inv(J, J_inv);
		mul(J_inv, f_X, X_delta);
		we -= X_delta[0];
		ze -= X_delta[1];
	}

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	// SDW: This is not accurate
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_ligas_util;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -8,
	/*.algo_author                 =*/ "M. Ligas, P. Banasik",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/258491472_Conversion_between_Cartesian_and_geodetic_coordinates_on_a_rotational_ellipsoid_by_solving_a_system_of_nonlinear_equations",
	/*.citation                    =*/ R"(Ligas, Marcin & Banasik, Piotr. (2011). Conversion between Cartesian and geodetic coordinates on a rotational ellipsoid by solving a system of nonlinear equations. Geodesy and Cartography. 60. 145-159. 10.2478/v10277-012-0013-x.)"
);

}
// }}}

namespace lin_wang_1995_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

#if 0
	const auto tmp_a = ell.a2 * z2 + ell.b2 * w2;
	const auto tmp_b = 2 * (ell.a2 * ell.a2 * z2 + ell.b2 * ell.b2 * w2);

	auto m = (ell.a * ell.b * CB(std::sqrt(tmp_a)) - ell.a2 * ell.b2 * tmp_a) / tmp_b;
#else
	const auto tmp = z2 / ell.b2 + w2 / ell.a2;

	auto m = 0.5 * ell.a2 * ell.b2 * tmp * (std::sqrt(tmp) - 1) / (z2 / (1 - ell.e2) + w2 * (1 - ell.e2));
#endif

	for (int i = 1; i <= max_iterations; ++i)
	{
		m -= lin_wang_1995_delta_m(w2, z2, m);
	}

	const auto we = std::abs(w / (1 + 2 * m / ell.a2));
	const auto ze =         (z / (1 + 2 * m / ell.b2));

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_lin_wang_1995_delta_m;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "Kuo-Chi Lin, Jie Wang",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/BF00806742",
	/*.citation                    =*/ R"(Lin, KC. & Wang, J. Bulletin Géodésique (1995) 69: 300. https://doi.org/10.1007/BF00806742)"
);

}
// }}}

namespace lin_wang_1995_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

#if 0
	const auto tmp_a = ell.a2 * z2 + ell.b2 * w2;
	const auto tmp_b = 2 * (ell.a2 * ell.a2 * z2 + ell.b2 * ell.b2 * w2);

	auto m = (ell.a * ell.b * CB(std::sqrt(tmp_a)) - ell.a2 * ell.b2 * tmp_a) / tmp_b;
#else
	const auto tmp = z2 / ell.b2 + w2 / ell.a2;

	auto m = 0.5 * ell.a2 * ell.b2 * tmp * (std::sqrt(tmp) - 1) / (z2 / (1 - ell.e2) + w2 * (1 - ell.e2));
#endif

	for (int i = 1; i <= max_iterations; ++i)
	{
		m -= lin_wang_1995_delta_m(w2, z2, m);
	}

	const auto we = std::abs(w / (1 + 2 * m / ell.a2));
	const auto ze =         (z / (1 + 2 * m / ell.b2));

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_lin_wang_1995_delta_m;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Kuo-Chi Lin, Jie Wang",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/BF00806742",
	/*.citation                    =*/ R"(Lin, KC. & Wang, J. Bulletin Géodésique (1995) 69: 300. https://doi.org/10.1007/BF00806742)"
);

}
// }}}

namespace lin_wang_1995_customht_1
// {{{
{
#define USE_CUSTOM_HT

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

#if 0
	const auto tmp_a = ell.a2 * z2 + ell.b2 * w2;
	const auto tmp_b = 2 * (ell.a2 * ell.a2 * z2 + ell.b2 * ell.b2 * w2);

	auto m = (ell.a * ell.b * CB(std::sqrt(tmp_a)) - ell.a2 * ell.b2 * tmp_a) / tmp_b;
#else
	const auto tmp = z2 / ell.b2 + w2 / ell.a2;

	auto m = 0.5 * ell.a2 * ell.b2 * tmp * (std::sqrt(tmp) - 1) / (z2 / (1 - ell.e2) + w2 * (1 - ell.e2));
#endif

	for (int i = 1; i <= max_iterations; ++i)
	{
		m -= lin_wang_1995_delta_m(w2, z2, m);
	}

	const auto we = std::abs(w / (1 + 2 * m / ell.a2));
	const auto ze =         (z / (1 + 2 * m / ell.b2));

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_lin_wang_1995_delta_m;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "Kuo-Chi Lin, Jie Wang",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/BF00806742",
	/*.citation                    =*/ R"(Lin, KC. & Wang, J. Bulletin Géodésique (1995) 69: 300. https://doi.org/10.1007/BF00806742)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace lin_wang_1995_customht_2
// {{{
{
#define USE_CUSTOM_HT

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

#if 0
	const auto tmp_a = ell.a2 * z2 + ell.b2 * w2;
	const auto tmp_b = 2 * (ell.a2 * ell.a2 * z2 + ell.b2 * ell.b2 * w2);

	auto m = (ell.a * ell.b * CB(std::sqrt(tmp_a)) - ell.a2 * ell.b2 * tmp_a) / tmp_b;
#else
	const auto tmp = z2 / ell.b2 + w2 / ell.a2;

	auto m = 0.5 * ell.a2 * ell.b2 * tmp * (std::sqrt(tmp) - 1) / (z2 / (1 - ell.e2) + w2 * (1 - ell.e2));
#endif

	for (int i = 1; i <= max_iterations; ++i)
	{
		m -= lin_wang_1995_delta_m(w2, z2, m);
	}

	const auto we = std::abs(w / (1 + 2 * m / ell.a2));
	const auto ze =         (z / (1 + 2 * m / ell.b2));

	auto sin_lat = ze;
	auto cos_lat = we * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = hypot(w - we, z - ze);

	if (w + std::abs(z) < we + std::abs(ze))
	{
		ht = -ht;
	}
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_lin_wang_1995_delta_m;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Kuo-Chi Lin, Jie Wang",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/BF00806742",
	/*.citation                    =*/ R"(Lin, KC. & Wang, J. Bulletin Géodésique (1995) 69: 300. https://doi.org/10.1007/BF00806742)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace long_1974
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	const auto lat_c = std::atan2(z, w);

	// lat_c == geocentric latitude
	// r == geocentric distance, units of a

	lat_rad = lat_c + ell.f * std::sin(2 * lat_c) / r + ell.f * ell.f * ((1 / (w2 + z2) - 1 / (4 * r)) * std::sin(4 * lat_c));

	ht = ell.get_ht(w, z, lat_rad);
/*
#ifdef USE_CUSTOM_HT
	ht = (r - 1) + ell.f * 0.5 * (1 - std::cos(2 * lat_c)) + ell.f * ell.f * ((1 / (4 * r) - 1 / 16.0)) * (1 - std::cos(4 * lat_c));
#else
#endif
*/
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "S. Long",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19740021140.pdf",
	/*.citation                    =*/ R"(Derivation of Transformation Formulas Between Geocentric and Geodetic Coordinates for Nonzero Altitudes
By Sheila Ann T. Long
Langley Research Center
NASA Technical Note TN D-7522)"
);

}
// }}}

namespace naive_I_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	ht = 0;
	double Rn = 1;

	auto sin_lat = z * (Rn + ht);
	auto cos_lat = w * (Rn * (1 - ell.e2) + ht);

	for (int i = 1; i <= max_iterations; ++i)
	{
		normalize(cos_lat, sin_lat);
		Rn = ell.get_Rn(sin_lat);
		// SDW: this is the only call to the ell.get_ht with extra Rn arg
		ht = ell.get_ht(w, z, sin_lat, cos_lat, Rn);

		sin_lat = z * (Rn + ht);
		cos_lat = w * (Rn * (1 - ell.e2) + ht);
	}

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 0,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://clynchg3c.com/Technote/geodesy/coordcvt.pdf",
	/*.citation                    =*/ R"(Naive iteration: must guess initial lat and ht)"
);

}
// }}}

namespace naive_I_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	ht = 0;
	double Rn = 1;

	auto sin_lat = z * (Rn + ht);
	auto cos_lat = w * (Rn * (1 - ell.e2) + ht);

	for (int i = 1; i <= max_iterations; ++i)
	{
		normalize(cos_lat, sin_lat);
		Rn = ell.get_Rn(sin_lat);
		// SDW: this is the only call to the ell.get_ht with extra Rn arg
		ht = ell.get_ht(w, z, sin_lat, cos_lat, Rn);

		sin_lat = z * (Rn + ht);
		cos_lat = w * (Rn * (1 - ell.e2) + ht);
	}

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -3,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://clynchg3c.com/Technote/geodesy/coordcvt.pdf",
	/*.citation                    =*/ R"(Naive iteration: must guess initial lat and ht)"
);

}
// }}}

namespace naive_II_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);

	for (int i = 1; i <= max_iterations; ++i)
	{
		normalize(cos_lat, sin_lat);
		const auto Rn = ell.get_Rn(sin_lat);
		sin_lat = z + Rn * sin_lat * ell.e2 / 2;
		cos_lat = w - Rn * cos_lat * ell.e2 / 2;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 0,
	/*.algo_author                 =*/ "R. Hirvonen, H. Moritz",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://www.dtic.mil/docs/citations/AD0420541",
	/*.citation                    =*/ R"(Naive iteration: must guess initial lat
Hirvonen, R., and H. Moritz
Practical Computations of Gravity at High Altitudes
Report No. 27
Institude of Geodesy, Photogrammetry and Cartography
The Ohio State University
Columbus, 1963)"
);

}
// }}}

namespace naive_II_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);

	for (int i = 1; i <= max_iterations; ++i)
	{
		normalize(cos_lat, sin_lat);
		const auto Rn = ell.get_Rn(sin_lat);
		sin_lat = z + Rn * sin_lat * ell.e2 / 2;
		cos_lat = w - Rn * cos_lat * ell.e2 / 2;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -3,
	/*.algo_author                 =*/ "R. Hirvonen, H. Moritz",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://www.dtic.mil/docs/citations/AD0420541",
	/*.citation                    =*/ R"(Naive iteration: must guess initial lat
Hirvonen, R., and H. Moritz
Practical Computations of Gravity at High Altitudes
Report No. 27
Institude of Geodesy, Photogrammetry and Cartography
The Ohio State University
Columbus, 1963)"
);

}
// }}}

namespace newton_raphson_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= newton_raphson_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_newton_raphson_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -1,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/NewtonsMethod.html",
	/*.citation                    =*/ R"(Newton-Raphson Method: Δφ = f / f')"
);

}
// }}}

namespace newton_raphson_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= newton_raphson_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_newton_raphson_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -10,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/NewtonsMethod.html",
	/*.citation                    =*/ R"(Newton-Raphson Method: Δφ = f / f')"
);

}
// }}}

namespace newton_raphson_quick_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = newton_raphson_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_newton_raphson_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -1,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/NewtonsMethod.html",
	/*.citation                    =*/ R"(Newton-Raphson Method: Δφ = f / f')"
);

}
// }}}

namespace newton_raphson_quick_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = newton_raphson_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_newton_raphson_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/NewtonsMethod.html",
	/*.citation                    =*/ R"(Newton-Raphson Method: Δφ = f / f')"
);

}
// }}}

namespace olson_1996
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

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

	double s;
	double c;
	double ss;

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

#ifdef USE_CUSTOM_HT
	ht = f + m * p / 2;
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "D.K. Olson",
	/*.code_copyright              =*/ "\"U.S. Government work, U.S. copyright does not apply.\"",
	/*.license                     =*/ "Public Domain",
	/*.orig_impl_lang              =*/ "C",
	/*.url                         =*/ "https://www.researchgate.net/publication/3002552_Converting_Earth-Centered_Earth-Fixed_Coordinates_to_Geodetic_Coordinates",
	/*.citation                    =*/ R"(D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to geodetic coordinates," in IEEE Transactions on Aerospace and Electronic Systems, vol. 32, no. 1, pp. 473-476, Jan. 1996, doi: 10.1109/7.481290.
URL: https://ieeexplore.ieee.org/document/481290

Olson, D.K.. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. Aerospace and Electronic Systems, IEEE Transactions on. 32. 473 - 476. 10.1109/7.481290.

Converted to C++ and modified by Steven Ward.  No rights reserved.)"
);

}
// }}}

namespace olson_1996_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

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

	double s;
	double c;
	double ss;

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

#ifdef USE_CUSTOM_HT
	ht = f + m * p / 2;
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "D.K. Olson",
	/*.code_copyright              =*/ "\"U.S. Government work, U.S. copyright does not apply.\"",
	/*.license                     =*/ "Public Domain",
	/*.orig_impl_lang              =*/ "C",
	/*.url                         =*/ "https://www.researchgate.net/publication/3002552_Converting_Earth-Centered_Earth-Fixed_Coordinates_to_Geodetic_Coordinates",
	/*.citation                    =*/ R"(D. K. Olson, "Converting Earth-centered, Earth-fixed coordinates to geodetic coordinates," in IEEE Transactions on Aerospace and Electronic Systems, vol. 32, no. 1, pp. 473-476, Jan. 1996, doi: 10.1109/7.481290.
URL: https://ieeexplore.ieee.org/document/481290

Olson, D.K.. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. Aerospace and Electronic Systems, IEEE Transactions on. 32. 473 - 476. 10.1109/7.481290.

Converted to C++ and modified by Steven Ward.  No rights reserved.)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace openglobe
// {{{
{

constexpr int _line_begin = __LINE__;
template <std::floating_point T>
struct Vector3D
{
	T x{};
	T y{};
	T z{};

	Vector3D() = default;

	Vector3D(
		const T _x,
		const T _y,
		const T _z):
		x(_x),
		y(_y),
		z(_z)
	{}

	auto length_sq() const
	{
		return x*x + y*y + z*z;
	}

	auto length() const
	{
		return std::sqrt(length_sq());
	}

	auto normalize() const
	{
		const auto l = length();

		if (l == 0)
		{
			return Vector3D<T>{};
		}

		return this->scale(1 / l);
	}

	auto scale(const T s) const
	{
		return Vector3D<T>{x*s, y*s, z*s};
	}

	auto dot(const Vector3D<T>& that) const
	{
		return this->x * that.x + this->y * that.y + this->z * that.z;
	}

	auto subtract(const Vector3D<T>& that) const
	{
		return Vector3D<T>{this->x - that.x, this->y - that.y, this->z - that.z};
	}

	auto multiply(const Vector3D<T>& that) const
	{
		return Vector3D<T>{this->x * that.x, this->y * that.y, this->z * that.z};
	}
};

template <std::floating_point T>
struct Geodetic2D
{
	T lat_rad{};
	T lon_rad{};

	Geodetic2D() = default;

	Geodetic2D(
		const T _lat_rad,
		const T _lon_rad):
		lat_rad(_lat_rad),
		lon_rad(_lon_rad)
	{}

};

template <std::floating_point T>
struct Geodetic3D
{
	T lat_rad{};
	T lon_rad{};
	T ht{};

	Geodetic3D() = default;

	Geodetic3D(
		const T _lat_rad,
		const T _lon_rad,
		const T _ht):
		lat_rad(_lat_rad),
		lon_rad(_lon_rad),
		ht(_ht)
	{}

	Geodetic3D(const Geodetic2D<T>& g, T _ht = 0):
		lat_rad(g.lat_rad),
		lon_rad(g.lon_rad),
		ht(_ht)
	{}
};

template <std::floating_point T>
auto ScaleToGeodeticSurface(const Vector3D<T>& position)
{
	const auto x = position.x;
	const auto y = position.y;
	const auto z = position.z;

	const auto x2 = x * x;
	const auto y2 = y * y;
	const auto z2 = z * z;

	auto beta = 1 / std::sqrt(
			x2 / WGS84<T>.a2 +
			y2 / WGS84<T>.a2 +
			z2 / WGS84<T>.b2);
	auto n = Vector3D<T>{
			beta * x / WGS84<T>.a2,
			beta * y / WGS84<T>.a2,
			beta * z / WGS84<T>.b2}.length();
	auto alpha = (1 - beta) * (position.length() / n);

	T da = 0;
	T db = 0;

	T s = 0;
	T dSdA = 1;

	constexpr T dist_threshold = 1E-10;

	do
	{
		alpha -= s / dSdA;

		da = 1 + alpha / WGS84<T>.a2;
		db = 1 + alpha / WGS84<T>.b2;

		auto da2 = da * da;
		auto db2 = db * db;

		auto da3 = da * da2;
		auto db3 = db * db2;

		s = x2 / (WGS84<T>.a2 * da2) +
			y2 / (WGS84<T>.a2 * da2) +
			z2 / (WGS84<T>.b2 * db2) - 1;

		dSdA = -2 * (
				x2 / (WGS84<T>.a2 * WGS84<T>.a2 * da3) +
				y2 / (WGS84<T>.a2 * WGS84<T>.a2 * da3) +
				z2 / (WGS84<T>.b2 * WGS84<T>.b2 * db3));
	}
	while (std::abs(s) > dist_threshold);

	return Vector3D<T>{
			x / da,
			y / da,
			z / db};
}

template <std::floating_point T>
auto GeodeticSurfaceNormal(const Vector3D<T>& positionOnEllipsoid)
{
	return positionOnEllipsoid.multiply(Vector3D<T>{
				1 / WGS84<T>.a2,
				1 / WGS84<T>.a2,
				1 / WGS84<T>.b2}).normalize();
}

template <std::floating_point T>
auto ToGeodetic2D(const Vector3D<T>& positionOnEllipsoid)
{
	auto n = GeodeticSurfaceNormal(positionOnEllipsoid);
	return Geodetic2D<T>{
			std::asin(n.z / n.length()), // lat_rad
			/*std::atan2(n.y, n.x)*/0 // lon_rad is already calculated
		};
}

template <std::floating_point T>
auto ToGeodetic3D(const Vector3D<T>& position)
{
	const auto p = ScaleToGeodeticSurface(position);
	const auto h = position.subtract(p);
	// auto height = Math.sign(h.Dot(position)) * h.length();
	const auto height = std::copysign(h.length(), h.dot(position));
	const auto g2d = ToGeodetic2D(p);
	return Geodetic3D<T>{g2d, height};
}

void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
	const Vector3D<double> position{x, y, z};
	const auto result = ToGeodetic3D(position);
	lat_rad = result.lat_rad;
	lon_rad = result.lon_rad;
	ht = result.ht;
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

// this algorithm does not include the common declarations
const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra - _lines_common_first_decls,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -5,
	/*.algo_author                 =*/ "Cozzi and Ohlarik",
	/*.code_copyright              =*/ "Cozzi and Ohlarik",
	/*.license                     =*/ "MIT",
	/*.orig_impl_lang              =*/ "C#",
	/*.url                         =*/ "https://github.com/virtualglobebook/OpenGlobe/blob/master/Source/Core/Geometry/Ellipsoid.cs",
	/*.citation                    =*/ R"(OpenGlobe file:Source/Core/Geometry/Ellipsoid.cs)"
);

}
// }}}

namespace ozone_1985
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS_CHECKED

	// a2 - b2 = a2 * e2
	// SDW: NaNs at the equator
	//const auto M = 0.5 * (ell.a * w - ell.a2 * ell.e2) / (ell.b * std::abs(z));
	//const auto N = 0.5 * (ell.a * w + ell.a2 * ell.e2) / (ell.b * std::abs(z));
	//const auto M = 0.5 * ell.a * (w - ell.a * ell.e2) / (ell.b * std::abs(z));
	//const auto N = 0.5 * ell.a * (w + ell.a * ell.e2) / (ell.b * std::abs(z));
	const auto M = 0.5 * (w - ell.a * ell.e2) / ((1 - ell.f) * std::abs(z));
	const auto N = 0.5 * (w + ell.a * ell.e2) / ((1 - ell.f) * std::abs(z));

	const auto V = 4 * M * N + 1;
	const auto W = 2 * (N * N - M * M);
	const auto tmp1 = std::sqrt(CB(V / 3) + W * W / 4);
	const auto I = std::cbrt(tmp1 + W / 2) - std::cbrt(tmp1 - W / 2);
	const auto J = std::sqrt(2 * I + 4 * M * M);
	const auto K = 2 * (N - M * I) / J;
	const auto tmp2 = 2 * M + J;
	const auto G = tmp2 * tmp2 - 4 * (I - K);
	const auto u = (tmp2 + std::sqrt(G)) / 2;

	// tangent half-angle formula, but this is u**2 - 1 instead of 1 - u**2
	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
	auto sin_lat = 2 * u;
	auto cos_lat = (u * u - 1) * (1 - ell.f);

	if (z < 0)
		sin_lat = -sin_lat;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ true,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "M. Ozone",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/237622491_A_Comparative_Analysis_of_the_Performance_of_Iterative_and_Non-iterative_Solutions_to_the_Cartesian_to_Geodetic_Coordinate_Transformation",
	/*.citation                    =*/ R"(Ozone, M. I, (1985). Non-Iterative Solution of the φ Equation, Surveying and Mapping, Vol. 45, No. 2, pp.  169-171.
Copied from George P. Gerdan and Rodney E. Deakin)"
);

}
// }}}

namespace paul_1973
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto alpha = (w2 + ell.a2 * e4) / (1 - ell.e2);
	const auto alpha2 = alpha * alpha;

	// e2/(1-e2) == ep2
	//const auto beta = (w2 - ell.e2 * ell.a2) / (1 - ell.e2);
	const auto beta = (w2 - ell.a2 * e4) / (1 - ell.e2);
	const auto beta2 = beta * beta;

	const auto q = 1 + 13.5 * z2 * (alpha2 - beta2) / CB(beta + z2);
	//const auto q = 1 + 13.5 * z2 * (alpha2 - beta) / CB(beta + z2);
	//const auto q2 = q * q;
	const auto sqrt_q2m1 = std::sqrt(q * q - 1 < 0 ? 0 : q * q - 1);
	const auto u = std::cbrt(q + sqrt_q2m1) + std::cbrt(q - sqrt_q2m1);
	//const auto t = beta + z2 * (p + 1 / p) / 12 - beta / 6 + z2 / 12;
	const auto t = ((beta + z2) / 12) * u - beta / 6 + z2 / 12;

	auto sqrt_t = std::sqrt(std::abs(t));
#if 0
	if (z < 0)
		sqrt_t = -sqrt_t;

	auto root_2 = std::sqrt(std::abs(z2 / 4 - beta / 2 - t + 0.25 * alpha * z / sqrt_t));
	// SDW: this doesn't work either
	//auto root_2 = std::sqrt(std::abs(z2 / 4 - beta / 2 - t + 0.25 * ell.a * z / sqrt_t));
	if (z < 0)
		root_2 = -root_2;

	const auto sigma = root_1 + root_2;

	auto sin_lat = sigma + z / 2;
	auto cos_lat = w;
#else
	auto sin_lat = z / 2 + sqrt_t + std::sqrt(z2 / 4 - beta / 2 - t + 0.25 * alpha * z / sqrt_t);
	auto cos_lat = w;
#endif

	/*
	if (z < 0)
		sin_lat = -sin_lat;
	*/

	//tan_lat = (z / 2 + root_1 + std::sqrt(z2 / 4 - beta / 2 - t + alpha * z / (4 * root_1))) / w;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "M.K. Paul",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/BF02522075",
	/*.citation                    =*/ R"(Paul, M.K. Bull. Geodesique (1973) 108: 135. https://doi.org/10.1007/BF02522075)");

}
// }}}

namespace pollard_2002_ht_1
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double ze = 0;
	//double we = 0;
	double we_N = 0;
	double m;
	double n;
	double r_;
	double s;
	double t;
	double sin_lat = 0;
	double cos_lat = 0;

	ze = ell.b * z / r;

	we_N = hypot(w, z + ell.ep2 * ze);
	m = w / we_N;
	n = z / we_N;
	r_ = m * m + n * n / (1 - ell.e2);
	s = m * w + z / (1 - ell.e2);
	t = w2 + z2 / (1 - ell.e2) - ell.a2;
	ht = (s - std::sqrt(std::abs(s * s - r_ * t))) / r_;

	// i = 1
	ze = z - n * ht;

	we_N = hypot(w, z + ell.ep2 * ze);
	m = w / we_N;
	n = z / we_N;
	r_ = m * m + n * n / (1 - ell.e2);
	s = m * w + z / (1 - ell.e2);
	t = w2 + z2 / (1 - ell.e2) - ell.a2;
	ht = (s - std::sqrt(std::abs(s * s - r_ * t))) / r_;

	ze = z - n * ht;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;
	// SDW: is this necessary?
	sin_lat *= (1 - ell.f);
	//normalize(cos_lat, sin_lat);

	lat_rad = std::atan2(sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "J. Pollard",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226563131_Iterative_vector_methods_for_computing_geodetic_latitude_and_height_from_rectangular_coordinates",
	/*.citation                    =*/ R"(Pollard, J. (2002). Iterative vector methods for computing geodetic latitude and height from rectangular coordinates. Journal of Geodesy. 76. 36-40. 10.1007/s001900100220.
https://link.springer.com/article/10.1007%2Fs001900100220)"
);

}
// }}}

namespace pollard_2002_naive_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double ze = 0;
	double sin_lat = z / r;
	double cos_lat = 0;

	ze = ell.b * sin_lat;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	// i = 1
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);
	ze = ell.b * sin_lat;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -2,
	/*.algo_author                 =*/ "J. Pollard",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226563131_Iterative_vector_methods_for_computing_geodetic_latitude_and_height_from_rectangular_coordinates",
	/*.citation                    =*/ R"(Pollard, J. (2002). Iterative vector methods for computing geodetic latitude and height from rectangular coordinates. Journal of Geodesy. 76. 36-40. 10.1007/s001900100220.
https://link.springer.com/article/10.1007%2Fs001900100220)"
);

}
// }}}

namespace pollard_2002_naive_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double ze = 0;
	double sin_lat = z / r;
	double cos_lat = 0;

	ze = ell.b * sin_lat;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	// i = 1
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);
	ze = ell.b * sin_lat;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	// i = 2
	sin_lat *= (1 - ell.f);
	normalize(cos_lat, sin_lat);
	ze = ell.b * sin_lat;
	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "J. Pollard",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226563131_Iterative_vector_methods_for_computing_geodetic_latitude_and_height_from_rectangular_coordinates",
	/*.citation                    =*/ R"(Pollard, J. (2002). Iterative vector methods for computing geodetic latitude and height from rectangular coordinates. Journal of Geodesy. 76. 36-40. 10.1007/s001900100220.
https://link.springer.com/article/10.1007%2Fs001900100220)"
);

}
// }}}

namespace pollard_2002_newton_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double ze = 0;
	double sin_lat = z / r;
	double cos_lat = 0;
	double tmp = 0;

	ze = ell.b * sin_lat;

	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;
	sin_lat *= (1 - ell.f);

	// i = 1
	normalize(cos_lat, sin_lat);
#if 0
	tmp = ell.a * ell.e2 * CB(cos_lat) / w;
	ze = (ell.b * sin_lat - ze * tmp) / (1 - tmp);
#else
	tmp = ell.a * ell.e2 * CB(cos_lat);
	ze = (ell.b * sin_lat * w - ze * tmp) / (w - tmp);
#endif

	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "J. Pollard",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226563131_Iterative_vector_methods_for_computing_geodetic_latitude_and_height_from_rectangular_coordinates",
	/*.citation                    =*/ R"(Pollard, J. (2002). Iterative vector methods for computing geodetic latitude and height from rectangular coordinates. Journal of Geodesy. 76. 36-40. 10.1007/s001900100220.
https://link.springer.com/article/10.1007%2Fs001900100220)"
);

}
// }}}

namespace pollard_2002_newton_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto r = std::sqrt(w2 + z2);

	double ze = 0;
	double sin_lat = z / r;
	double cos_lat = 0;
	double tmp = 0;

	ze = ell.b * sin_lat;

	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;
	sin_lat *= (1 - ell.f);

	// i = 1
	normalize(cos_lat, sin_lat);
#if 0
	tmp = ell.a * ell.e2 * CB(cos_lat) / w;
	ze = (ell.b * sin_lat - ze * tmp) / (1 - tmp);
#else
	tmp = ell.a * ell.e2 * CB(cos_lat);
	ze = (ell.b * sin_lat * w - ze * tmp) / (w - tmp);
#endif

	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;
	sin_lat *= (1 - ell.f);

	// i = 2
	normalize(cos_lat, sin_lat);
#if 0
	tmp = ell.a * ell.e2 * CB(cos_lat) / w;
	ze = (ell.b * sin_lat - ze * tmp) / (1 - tmp);
#else
	tmp = ell.a * ell.e2 * CB(cos_lat);
	ze = (ell.b * sin_lat * w - ze * tmp) / (w - tmp);
#endif

	sin_lat = z + ell.ep2 * ze;
	cos_lat = w;

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "J. Pollard",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/226563131_Iterative_vector_methods_for_computing_geodetic_latitude_and_height_from_rectangular_coordinates",
	/*.citation                    =*/ R"(Pollard, J. (2002). Iterative vector methods for computing geodetic latitude and height from rectangular coordinates. Journal of Geodesy. 76. 36-40. 10.1007/s001900100220.
https://link.springer.com/article/10.1007%2Fs001900100220)"
);

}
// }}}

namespace schroder_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= schroder_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_schroder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -1,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/SchroedersMethod.html",
	/*.citation                    =*/ R"(Schroeder's Method: Δφ = (f * f') / (f' * f' - f * f''))"
);

}
// }}}

namespace schroder_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	lat_rad = std::atan2(sin_lat, cos_lat);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		lat_rad -= schroder_delta_lat(w, z, sin_lat, cos_lat);

		sin_lat = std::sin(lat_rad);
		cos_lat = std::cos(lat_rad);
	}

	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_schroder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -10,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/SchroedersMethod.html",
	/*.citation                    =*/ R"(Schroeder's Method: Δφ = (f * f') / (f' * f' - f * f''))"
);

}
// }}}

namespace schroder_quick_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = schroder_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_schroder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -1,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/SchroedersMethod.html",
	/*.citation                    =*/ R"(Schroeder's Method: Δφ = (f * f') / (f' * f' - f * f''))"
);

}
// }}}

namespace schroder_quick_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.e2);
	normalize(sin_lat, cos_lat);

	for (int i = 1; i <= max_iterations; ++i)
	{
		const auto delta_lat_rad = schroder_delta_lat(w, z, sin_lat, cos_lat);

		const auto sin_delta_lat = std::sin(delta_lat_rad);
		const auto cos_delta_lat = std::cos(delta_lat_rad);

		const auto next_sin_lat = sin_lat * cos_delta_lat - cos_lat * sin_delta_lat;
		const auto next_cos_lat = cos_lat * cos_delta_lat + sin_lat * sin_delta_lat;

		sin_lat = next_sin_lat;
		cos_lat = next_cos_lat;
	}

	lat_rad = std::atan2(sin_lat, cos_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_schroder_delta_lat;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Steven Ward",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://mathworld.wolfram.com/SchroedersMethod.html",
	/*.citation                    =*/ R"(Schroeder's Method: Δφ = (f * f') / (f' * f' - f * f''))"
);

}
// }}}

namespace sedris
// {{{
{

constexpr int _line_begin = __LINE__;
template <std::floating_point T>
auto gee(const T h, const T rn)
{
	return (rn + h) / ((1 - WGS84<T>.e2) * rn + h);
}

/*
 * STRUCT: SRM_GC_GD_Specific_Constants
 *
 *   Used to hold precomputed constants specific to a GC to GD transformation
 */
template <std::floating_point T>
struct SRM_GC_GD_Specific_Constants
{
	T b1[5];
	T b2[5];
	T b3[5];
	T b4[5];
	T b5[5];
	T u[5];
	T v[5];
};

template <std::floating_point T>
void set_gc_to_gd_constants(SRM_GC_GD_Specific_Constants<T>& gc_gd_spec)
{
	/*old function prototype
	void tf_set_gc_to_gd_constants
	(
	const SE_ERM_DATA                  *e_constants,
	const SE_SRF_PARAMETERS            *srf,
	TF_GC_GD_SPECIFIC_CONSTANTS *gc_gd_spec
	)
*/

	T del[5] = {NAN, NAN, NAN, NAN, NAN};
	T hmn = NAN;
	T hmx = NAN;
	T g1 = NAN;
	T g2 = NAN;
	T g3 = NAN;
	T g4 = NAN;
	T gm = NAN;
	T hm = NAN;
	T d1 = NAN;
	T d2 = NAN;
	T d3 = NAN;
	T d4 = NAN;
	T d5 = NAN;
	T d6 = NAN;
	T sm = NAN;
	T rnm = NAN;
	T zm = NAN;
	T wm = NAN;
	T z2 = NAN;
	T w2 = NAN;
	T d7 = NAN;
	T d8 = NAN;
	T d9 = NAN;
	T d10 = NAN;
	T a1 = NAN;
	T a2 = NAN;
	T a3 = NAN;
	T a4 = NAN;
	T a5 = NAN;
	T f1 = NAN;
	T f2 = NAN;

	//if (WGS84<T>.e > 1E-12)
	if (WGS84<T>.e2 > 1E-24)
	{
		del[0] = 3E4;
		del[1] = 5E4;
		del[2] = 2.2E7 - 5E4;
		del[3] = 4.5E8 - 2.2E7;
		del[4] = 1E10;

		hmn = -30000;

		for (int i = 0; i < 5; i++)
		{
			hmx = hmn + del[i];

			g1 = gee(hmn, WGS84<T>.a);  /* changed 0 to 1E-14 */
			g2 = gee(hmx, WGS84<T>.a);

			g3 = gee(hmx, WGS84<T>.a / (1 - WGS84<T>.f));
			g4 = gee(hmn, WGS84<T>.a / (1 - WGS84<T>.f));

			hm = (hmx - hmn) * 0.5;

			gm = gee(hm, WGS84<T>.a / std::sqrt(1 - WGS84<T>.e2 / 2));

			d1 = SQ(WGS84<T>.b + hmx) * g3 - SQ(WGS84<T>.b + hmn) * g4;
			d1 = d1 / (SQ(WGS84<T>.b + hmx) - SQ(WGS84<T>.b + hmn));

			d2 = -(g4 - g3) / (SQ(WGS84<T>.b + hmx) - SQ(WGS84<T>.b + hmn));

			d3 = SQ(WGS84<T>.b + hmx) * d2 - g3;

			d4 = SQ(WGS84<T>.b + hmx) * (g3 - d1);

			d5 = (g1 + d3) /  SQ(WGS84<T>.a + hmn) - (g2 + d3) / SQ(WGS84<T>.a + hmx);

			d5 = d5 / (g2 - g1);

			d6 = 1 / SQ(WGS84<T>.a + hmx) - 1 / SQ(WGS84<T>.a + hmn);

			d6 = d6 * d4 / (g2 - g1);

			sm = M_SQRT1_2;  /*sin(π/4) == 1/sqrt(2)*/

			//rnm = WGS84<T>.a / std::sqrt(1 - WGS84<T>.e2 * sm * sm);
			rnm = ell.get_Rn(sm);

			zm = ((1 - WGS84<T>.e2) * rnm + hm) * sm;
			wm = (rnm + hm) * M_SQRT1_2;  /*cos(π/4) == 1/sqrt(2)*/

			z2 = zm * zm;
			w2 = wm * wm;

			d7 = (z2 * d2 - d3 - gm - gm * w2 * d5) / w2;

			d8 = (z2 * gm - d4 - z2 * d1 + gm * w2 * d6) / w2;

			d9 = (g2 + d3) / SQ(WGS84<T>.a + hmx) + g2 * d5;

			d10 = -d4 / SQ(WGS84<T>.a + hmx) + d6 * g2;

			a4 = (d8 - d10) / (d9 + d7);
			a2 = d9 * a4 + d10;
			a5 = d5 * a4 + d6;
			a1 = d4 - d3 * a4;
			a3 = d1 + d2 * a4;

			gc_gd_spec.b1[i] = a3;
			gc_gd_spec.b2[i] = a1 - a3 * a4;
			gc_gd_spec.b3[i] = a2 - a3 * a5;
			gc_gd_spec.b4[i] = a4;
			gc_gd_spec.b5[i] = a5;

			f1 = SQ(WGS84<T>.a + hmn - 1E-8);
			f2 = SQ(WGS84<T>.b + hmn - 1E-8);

			gc_gd_spec.u[i] = f1 / f2;
			gc_gd_spec.v[i] = f1;
			hmn = hmx;
		}
	}
} /* end set_gc_to_gd_constants */

/*!This macro is from the AMIP codebase (CDK transplanted 8/19/2003)
It quickly computes the The Rn earth radius.  Implicit
in it are assumptions about the eccentricity of the
ellipsoids being less than that for Clark1880 ellipsoid.
The routine uses empirically optimized coefficients
whose derivation is not directly obvious from the code.
The routine also assumes the variable e_constants is
in scope in the caller of this macro and that the A
member contains the semi-major ellipsoid axis.
Algorithm derived by Ralph Toms, SRI.
*/
#define COMPUTE_RN_FAST(arg,answer) \
{\
  const auto _alpha = 1.004244;\
  const auto _d1 = -0.5 * ell.e2;\
  const auto _x = _d1 * arg;\
  const auto _ak = 0.5 + _x;\
  const auto _z = 1 - _alpha * _x;\
  answer = ell.a * _z * (1.5 - _ak * _z * _z);\
}

//unsigned int SRM_ChangeGC_GD(
	//void                   *constants,
	//const double          source_generic_coord[4],
	//double          dest_generic_coord[4],
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	constexpr short int REGION_1 = 0;
	constexpr short int REGION_2 = 1;
	constexpr short int REGION_3 = 2;
	constexpr short int REGION_4 = 3;
	constexpr short int REGION_5 = 4;
	constexpr short int REGION_UNDEFINED = 5;
	constexpr short int REGION_SPHERICAL = 6;

	/*
	The old amip prototype:

	SE_COORD_STATUS_CODE_ENUM
	tf_gc_gd_conversion(const TF_CONVERT_COORD_SYSTEM_PAIR  convert_params_ptr,
	const SE_COORD* src_coord_ptr,
	SE_COORD* dest_coord_ptr)
	*/

	double arg;
	double p;
	double s1;
	[[gnu::unused]] double s12;

	double alpha;
	double arg2;
	double cf;
	double cl;
	double q;
	double r_1;
	double r_2;
	double ro;
	double roe;
	double rr;
	double s;
	double top;
	double top2;
	double v;
	double xarg;
	double zo;

	double lowerBound;
	double upperBound;
	double ga;
	double bot;

	short int region = REGION_UNDEFINED;

	static SRM_GC_GD_Specific_Constants<double> gc_gd_spec;

	static bool first_run = true;

	if (first_run)
	{
		set_gc_to_gd_constants(gc_gd_spec);
		first_run = false;
	}

	//unsigned int status = 0;

	bool done_by_special_case = false;

	//if (std::abs(ell.e) <= 1E-12)
	if (ell.e2 <= 1E-24)
	{
		const auto r = std::sqrt(w2 + z2);

		/*This is the spherical special case*/
		/*lon_rad=std::atan2(y, x);*/
		/*This is done below at the end to avoid duplicating the code*/
		lat_rad = std::atan2(z /*z*/, w);
		//ht = std::sqrt(w2 + z2) - ell.a;
		ht = r - ell.a;
		region = REGION_SPHERICAL;
		goto END_REGION_CHECK;
	}
	else
	{
		/* Check for special cases */
		if (std::abs(x) <= 1E-12)
		{
			if (y > 0)
			{
				//lon_rad = M_PI_2;
			}
			else
			{
				if (y < 0)
				{
					//lon_rad = -M_PI_2;
				}
				else /* y == 0 */
				{
					if (z >= 0)
					{
						lat_rad = M_PI_2;
						//lon_rad = 0;
						ht = z - ell.b;

						done_by_special_case = true;

					} /* end if z> 0 */
					else if (z <= 0)
					{
						lat_rad = -M_PI_2;
						//lon_rad =  0;
						ht = -(z + ell.b);

						done_by_special_case = true;

					} /* end if z < 0 */
					else
					{
						/* at origin */

						//return 32;
						return;

					} /* end else at origin */
				} /* end else y == 0 */
			} /* end else y < 0 */
		} /* end if x == 0 */
	}
	/* end of special cases */

	if (!done_by_special_case)
	{
		/* region 2      0-50 kilometers*/

		lowerBound = w2 + gc_gd_spec.u[REGION_2] * z2;
		upperBound = w2 + gc_gd_spec.u[REGION_3] * z2;
		if ((lowerBound >= gc_gd_spec.v[REGION_2]) && (upperBound <= gc_gd_spec.v[REGION_3]))
		{
			region = REGION_2;
			goto END_REGION_CHECK;
		}
		else
		{
			/* region 3  50 - 23000 kilometers*/
			lowerBound = upperBound;
			upperBound = w2 + z2 * gc_gd_spec.u[REGION_4];
			if ((lowerBound >= gc_gd_spec.v[REGION_3]) && (upperBound <= gc_gd_spec.v[REGION_4]))
			{
				region = REGION_3;
				goto END_REGION_CHECK;
			}
			else
			{
				/* region 1 -30 to 0 kilometers */
				lowerBound = w2 + z2 * gc_gd_spec.u[REGION_1];
				upperBound = w2 + z2 * gc_gd_spec.u[REGION_2];
				if ((lowerBound >= gc_gd_spec.v[REGION_1]) && (upperBound <= gc_gd_spec.v[REGION_2]))
				{
					region = REGION_1;
					goto END_REGION_CHECK;
				}
				else
				{
					/* region 4  23000 to 10e6 kilometers */
					lowerBound = upperBound;
					upperBound = w2 + z2 * gc_gd_spec.u[REGION_5];
					if ((lowerBound >= gc_gd_spec.v[REGION_4]) && (upperBound <= gc_gd_spec.v[REGION_5]))
					{
						region = REGION_4;
						goto END_REGION_CHECK;
					}
					else
					{
						/*Declare region 5 unless the following test fails*/
						region=REGION_5;
						/* region 5 < -30 or > 10e6 kilometers*/

					}
				}
			}
		}

END_REGION_CHECK:

		/* Approximation to g function*/

		ga = gc_gd_spec.b1[region] +
			(gc_gd_spec.b2[region] + gc_gd_spec.b3[region] * w2) /
			(gc_gd_spec.b4[region] + gc_gd_spec.b5[region] * w2 + z2);

		/*
		GA=B1(II)+(B2(II)+B3(II)*W2)/(B4(II)+B5(II)*W2+Z2)
		*/

		switch (region)
		{
		case REGION_1:
		case REGION_2:
			top = z * ga;
			lat_rad = std::atan2(top, w);
			top2 = top * top;
			rr = top2 + w2;
			q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
			s12 = top2 / rr;

#if 1
			double Rn;

			/*uses a Newton-Raphson single iteration with
			excellent first guess usin only multiplications*/

			COMPUTE_RN_FAST(s12, Rn);
#if 0
			fprintf(stderr,"rn fast is %.15f\n",Rn);
#endif

#endif

			if (s12 <= 0.5)  /* If between +- 45 degrees lattitude */
			{
				ht = q - Rn;
			}
			else
			{
				//ht = q / ga + ell.a * (-(1 - ell.e2)) * Rn / ell.a;
				ht = q / ga - (1 - ell.e2) * Rn;
			}
#endif

			/******************************************************************
			Done below at end of function as optimization
			 lon_rad = std::atan2(src.source_y, src.source_x);
			****************************************************************/
			break;

		case REGION_3:
		case REGION_4:

			{

				/* correct g by using it as a first guess into the bowring formula*/
				top = z * ga * (1 - ell.f);
				top2 = top * top;

				rr = top2 + w2;
				q = std::sqrt(rr);

				{
					double sn;
					double cn;
					double s3;
					double c3;
					sn = top / q;
					cn = w / q;
					s3 = CB(sn);
					c3 = CB(cn);
					top = z + ell.b * ell.ep2 * s3;
					bot = w - ell.e2 * ell.a * c3;
				}
				top2 = top * top;

				rr = top2 + bot * bot;
				q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
				s12 = top2 / rr;

				/*Fast Rn computation*/
				{
					double Rn;
					COMPUTE_RN_FAST(s12, Rn);
					/*Fast Computation of Rn = a/std::sqrt(1-Eps2*sin_squared(latitude))*/

					if (s12 <= 0.5)
					{
						ht = w * q / bot - Rn;
					}
					else
					{
						ht = z * q / top - (1 - ell.e2) * Rn;
					}
				}
#endif

				lat_rad = std::atan2(top, bot);
			}
			break;
		case REGION_5:
			{
				double gee;

#if 0
				/*!\note This check has been replaced by a higer level check stipulating
				that we are at no lesser altitude than negative b.
				*/

				/*
				 * we are not allowed to be within 50 Kilometers of
				 * the center of the earth
				 */

				if (w2 + z2 < 50000
				   )
				{
#ifdef CONV_DEBUG
					fprintf(f_ptr, "[SEDRIS CONVERSIONS] Point within 50 "\
							"kilometers of the earth's center.  Outside "\
							"SEDRIS requirement range\n");
#endif
					return 32;
				}
#endif

				cf = 54 * ell.b2 * z2;
				gee = w2 + (1 - ell.e2) * z2 - ell.a2 * e4;
				{
					alpha = cf / (gee * gee);
					cl = e4 * w2 * alpha / gee;
				}
				arg2 = cl * (cl + 2);
				s1 = 1 + cl + std::sqrt(arg2);

				//s = pow(s1, .333333333333333333333333333333);
				s = std::cbrt(s1);

				{
					auto temp = s / (s * s + 1 + s);
					p = alpha * temp * temp / 3;
				}

				xarg = 1 + (2 * e4 * p);
				q = std::sqrt(xarg);

				{
					r_2 = -p * (2 * (1 - ell.e2) * z2 / (q * (1 + q)) + w2);
					r_1 = 1 + 1 / q;
					r_2 /= ell.a2;
					/*
					 * DUE TO PRECISION ERRORS THE ARGUMENT MAY BECOME
					 * NEGATIVE. IF SO SET THE ARGUMENT TO ZERO.
					 */
					if (r_1 + r_2 > 0)
					{
						ro = ell.a * std::sqrt(0.5 * (r_1 + r_2));
					}
					else
					{
						ro = 0;
					}
					ro -= p * ell.e2 * w / (1 + q);
				}

				roe = ell.e2 * ro;
				arg = SQ(w - roe) + z2;
				v   = std::sqrt(arg - ell.e2 * z2);
				{
					zo = ell.a * (1 - ell.e2) * z / v;
#ifdef USE_CUSTOM_HT
					ht = std::sqrt(arg) * (1 - ell.a * (1 - ell.e2) / v);
#endif
				}

				top = z + ell.ep2 * zo;
				lat_rad = std::atan2(top, w);

				/************************************************************

				As an optimization, this is done below
				 lon_rad = std::atan2(src.source_y, src.source_x);
				************************************************************/

				/* end of Exact solution */
			}
			break;
		case REGION_SPHERICAL:
			/*This case gets out of the switch statement and lets the longitude
			 *be computed at the end of the routine with a single return statement
			 *Conclusion:  Don't Remove this case or there won't be a good way
			 *to skip the region test and this won't work very well.
			 */
			break;

		default:
			//status = 32;
#if 0
			fprintf(stderr,"GC to GD: region error, region = %i line %d \n",region,__LINE__);
#endif
			break;
		}

		/*Since logitude calculation is common to all regions and requires lots of time,
		we'll just do it here to break interlock and hopefully have it occur
		during some of the function return bureaucracy*/

		//lon_rad = std::atan2(y, x);
	}

	//return status;

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
#undef COMPUTE_RN_FAST
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -5,
	/*.algo_author                 =*/ "Cameron Kellough?",
	/*.code_copyright              =*/ "SEDRIS",
	/*.license                     =*/ "\"This software was developed for use by the United States Government with unlimited rights.\"",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "http://www.sedris.org/sdk_4.0.1/sdk_c.htm",
	/*.citation                    =*/ R"(sedris_c_sdk_4.1.0
file:src/lib/srm/src/lib/cpp_impl/srf_impl/srm_gc.cpp)"
);

}
// }}}

namespace sedris_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
template <std::floating_point T>
auto gee(const T h, const T rn)
{
	return (rn + h) / ((1 - WGS84<T>.e2) * rn + h);
}

/*
 * STRUCT: SRM_GC_GD_Specific_Constants
 *
 *   Used to hold precomputed constants specific to a GC to GD transformation
 */
template <std::floating_point T>
struct SRM_GC_GD_Specific_Constants
{
	T b1[5];
	T b2[5];
	T b3[5];
	T b4[5];
	T b5[5];
	T u[5];
	T v[5];
};

template <std::floating_point T>
void set_gc_to_gd_constants(SRM_GC_GD_Specific_Constants<T>& gc_gd_spec)
{
	/*old function prototype
	void tf_set_gc_to_gd_constants
	(
	const SE_ERM_DATA                  *e_constants,
	const SE_SRF_PARAMETERS            *srf,
	TF_GC_GD_SPECIFIC_CONSTANTS *gc_gd_spec
	)
*/

	T del[5] = {NAN, NAN, NAN, NAN, NAN};
	T hmn = NAN;
	T hmx = NAN;
	T g1 = NAN;
	T g2 = NAN;
	T g3 = NAN;
	T g4 = NAN;
	T gm = NAN;
	T hm = NAN;
	T d1 = NAN;
	T d2 = NAN;
	T d3 = NAN;
	T d4 = NAN;
	T d5 = NAN;
	T d6 = NAN;
	T sm = NAN;
	T rnm = NAN;
	T zm = NAN;
	T wm = NAN;
	T z2 = NAN;
	T w2 = NAN;
	T d7 = NAN;
	T d8 = NAN;
	T d9 = NAN;
	T d10 = NAN;
	T a1 = NAN;
	T a2 = NAN;
	T a3 = NAN;
	T a4 = NAN;
	T a5 = NAN;
	T f1 = NAN;
	T f2 = NAN;

	//if (WGS84<T>.e > 1E-12)
	if (WGS84<T>.e2 > 1E-24)
	{
		del[0] = 3E4;
		del[1] = 5E4;
		del[2] = 2.2E7 - 5E4;
		del[3] = 4.5E8 - 2.2E7;
		del[4] = 1E10;

		hmn = -30000;

		for (int i = 0; i < 5; i++)
		{
			hmx = hmn + del[i];

			g1 = gee(hmn, WGS84<T>.a);  /* changed 0 to 1E-14 */
			g2 = gee(hmx, WGS84<T>.a);

			g3 = gee(hmx, WGS84<T>.a / (1 - WGS84<T>.f));
			g4 = gee(hmn, WGS84<T>.a / (1 - WGS84<T>.f));

			hm = (hmx - hmn) * 0.5;

			gm = gee(hm, WGS84<T>.a / std::sqrt(1 - WGS84<T>.e2 / 2));

			d1 = SQ(WGS84<T>.b + hmx) * g3 - SQ(WGS84<T>.b + hmn) * g4;
			d1 = d1 / (SQ(WGS84<T>.b + hmx) - SQ(WGS84<T>.b + hmn));

			d2 = -(g4 - g3) / (SQ(WGS84<T>.b + hmx) - SQ(WGS84<T>.b + hmn));

			d3 = SQ(WGS84<T>.b + hmx) * d2 - g3;

			d4 = SQ(WGS84<T>.b + hmx) * (g3 - d1);

			d5 = (g1 + d3) /  SQ(WGS84<T>.a + hmn) - (g2 + d3) / SQ(WGS84<T>.a + hmx);

			d5 = d5 / (g2 - g1);

			d6 = 1 / SQ(WGS84<T>.a + hmx) - 1 / SQ(WGS84<T>.a + hmn);

			d6 = d6 * d4 / (g2 - g1);

			sm = M_SQRT1_2;  /*sin(π/4) == 1/sqrt(2)*/

			//rnm = WGS84<T>.a / std::sqrt(1 - WGS84<T>.e2 * sm * sm);
			rnm = ell.get_Rn(sm);

			zm = ((1 - WGS84<T>.e2) * rnm + hm) * sm;
			wm = (rnm + hm) * M_SQRT1_2;  /*cos(π/4) == 1/sqrt(2)*/

			z2 = zm * zm;
			w2 = wm * wm;

			d7 = (z2 * d2 - d3 - gm - gm * w2 * d5) / w2;

			d8 = (z2 * gm - d4 - z2 * d1 + gm * w2 * d6) / w2;

			d9 = (g2 + d3) / SQ(WGS84<T>.a + hmx) + g2 * d5;

			d10 = -d4 / SQ(WGS84<T>.a + hmx) + d6 * g2;

			a4 = (d8 - d10) / (d9 + d7);
			a2 = d9 * a4 + d10;
			a5 = d5 * a4 + d6;
			a1 = d4 - d3 * a4;
			a3 = d1 + d2 * a4;

			gc_gd_spec.b1[i] = a3;
			gc_gd_spec.b2[i] = a1 - a3 * a4;
			gc_gd_spec.b3[i] = a2 - a3 * a5;
			gc_gd_spec.b4[i] = a4;
			gc_gd_spec.b5[i] = a5;

			f1 = SQ(WGS84<T>.a + hmn - 1E-8);
			f2 = SQ(WGS84<T>.b + hmn - 1E-8);

			gc_gd_spec.u[i] = f1 / f2;
			gc_gd_spec.v[i] = f1;
			hmn = hmx;
		}
	}
} /* end set_gc_to_gd_constants */

/*!This macro is from the AMIP codebase (CDK transplanted 8/19/2003)
It quickly computes the The Rn earth radius.  Implicit
in it are assumptions about the eccentricity of the
ellipsoids being less than that for Clark1880 ellipsoid.
The routine uses empirically optimized coefficients
whose derivation is not directly obvious from the code.
The routine also assumes the variable e_constants is
in scope in the caller of this macro and that the A
member contains the semi-major ellipsoid axis.
Algorithm derived by Ralph Toms, SRI.
*/
#define COMPUTE_RN_FAST(arg,answer) \
{\
  const auto _alpha = 1.004244;\
  const auto _d1 = -0.5 * ell.e2;\
  const auto _x = _d1 * arg;\
  const auto _ak = 0.5 + _x;\
  const auto _z = 1 - _alpha * _x;\
  answer = ell.a * _z * (1.5 - _ak * _z * _z);\
}

//unsigned int SRM_ChangeGC_GD(
	//void                   *constants,
	//const double          source_generic_coord[4],
	//double          dest_generic_coord[4],
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	constexpr short int REGION_1 = 0;
	constexpr short int REGION_2 = 1;
	constexpr short int REGION_3 = 2;
	constexpr short int REGION_4 = 3;
	constexpr short int REGION_5 = 4;
	constexpr short int REGION_UNDEFINED = 5;
	constexpr short int REGION_SPHERICAL = 6;

	/*
	The old amip prototype:

	SE_COORD_STATUS_CODE_ENUM
	tf_gc_gd_conversion(const TF_CONVERT_COORD_SYSTEM_PAIR  convert_params_ptr,
	const SE_COORD* src_coord_ptr,
	SE_COORD* dest_coord_ptr)
	*/

	double arg;
	double p;
	double s1;
	[[gnu::unused]] double s12;

	double alpha;
	double arg2;
	double cf;
	double cl;
	double q;
	double r_1;
	double r_2;
	double ro;
	double roe;
	double rr;
	double s;
	double top;
	double top2;
	double v;
	double xarg;
	double zo;

	double lowerBound;
	double upperBound;
	double ga;
	double bot;

	short int region = REGION_UNDEFINED;

	static SRM_GC_GD_Specific_Constants<double> gc_gd_spec;

	static bool first_run = true;

	if (first_run)
	{
		set_gc_to_gd_constants(gc_gd_spec);
		first_run = false;
	}

	//unsigned int status = 0;

	bool done_by_special_case = false;

	//if (std::abs(ell.e) <= 1E-12)
	if (ell.e2 <= 1E-24)
	{
		const auto r = std::sqrt(w2 + z2);

		/*This is the spherical special case*/
		/*lon_rad=std::atan2(y, x);*/
		/*This is done below at the end to avoid duplicating the code*/
		lat_rad = std::atan2(z /*z*/, w);
		//ht = std::sqrt(w2 + z2) - ell.a;
		ht = r - ell.a;
		region = REGION_SPHERICAL;
		goto END_REGION_CHECK;
	}
	else
	{
		/* Check for special cases */
		if (std::abs(x) <= 1E-12)
		{
			if (y > 0)
			{
				//lon_rad = M_PI_2;
			}
			else
			{
				if (y < 0)
				{
					//lon_rad = -M_PI_2;
				}
				else /* y == 0 */
				{
					if (z >= 0)
					{
						lat_rad = M_PI_2;
						//lon_rad = 0;
						ht = z - ell.b;

						done_by_special_case = true;

					} /* end if z> 0 */
					else if (z <= 0)
					{
						lat_rad = -M_PI_2;
						//lon_rad =  0;
						ht = -(z + ell.b);

						done_by_special_case = true;

					} /* end if z < 0 */
					else
					{
						/* at origin */

						//return 32;
						return;

					} /* end else at origin */
				} /* end else y == 0 */
			} /* end else y < 0 */
		} /* end if x == 0 */
	}
	/* end of special cases */

	if (!done_by_special_case)
	{
		/* region 2      0-50 kilometers*/

		lowerBound = w2 + gc_gd_spec.u[REGION_2] * z2;
		upperBound = w2 + gc_gd_spec.u[REGION_3] * z2;
		if ((lowerBound >= gc_gd_spec.v[REGION_2]) && (upperBound <= gc_gd_spec.v[REGION_3]))
		{
			region = REGION_2;
			goto END_REGION_CHECK;
		}
		else
		{
			/* region 3  50 - 23000 kilometers*/
			lowerBound = upperBound;
			upperBound = w2 + z2 * gc_gd_spec.u[REGION_4];
			if ((lowerBound >= gc_gd_spec.v[REGION_3]) && (upperBound <= gc_gd_spec.v[REGION_4]))
			{
				region = REGION_3;
				goto END_REGION_CHECK;
			}
			else
			{
				/* region 1 -30 to 0 kilometers */
				lowerBound = w2 + z2 * gc_gd_spec.u[REGION_1];
				upperBound = w2 + z2 * gc_gd_spec.u[REGION_2];
				if ((lowerBound >= gc_gd_spec.v[REGION_1]) && (upperBound <= gc_gd_spec.v[REGION_2]))
				{
					region = REGION_1;
					goto END_REGION_CHECK;
				}
				else
				{
					/* region 4  23000 to 10e6 kilometers */
					lowerBound = upperBound;
					upperBound = w2 + z2 * gc_gd_spec.u[REGION_5];
					if ((lowerBound >= gc_gd_spec.v[REGION_4]) && (upperBound <= gc_gd_spec.v[REGION_5]))
					{
						region = REGION_4;
						goto END_REGION_CHECK;
					}
					else
					{
						/*Declare region 5 unless the following test fails*/
						region=REGION_5;
						/* region 5 < -30 or > 10e6 kilometers*/

					}
				}
			}
		}

END_REGION_CHECK:

		/* Approximation to g function*/

		ga = gc_gd_spec.b1[region] +
			(gc_gd_spec.b2[region] + gc_gd_spec.b3[region] * w2) /
			(gc_gd_spec.b4[region] + gc_gd_spec.b5[region] * w2 + z2);

		/*
		GA=B1(II)+(B2(II)+B3(II)*W2)/(B4(II)+B5(II)*W2+Z2)
		*/

		switch (region)
		{
		case REGION_1:
		case REGION_2:
			{
			top = z * ga;
			lat_rad = std::atan2(top, w);
			top2 = top * top;
			rr = top2 + w2;
			q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
			s12 = top2 / rr;

#if 1
			double Rn;

			/*uses a Newton-Raphson single iteration with
			excellent first guess usin only multiplications*/

			COMPUTE_RN_FAST(s12, Rn);
#if 0
			fprintf(stderr,"rn fast is %.15f\n",Rn);
#endif

#endif

			if (s12 <= 0.5)  /* If between +- 45 degrees lattitude */
			{
				ht = q - Rn;
			}
			else
			{
				//ht = q / ga + ell.a * (-(1 - ell.e2)) * Rn / ell.a;
				ht = q / ga - (1 - ell.e2) * Rn;
			}
#endif

			/******************************************************************
			Done below at end of function as optimization
			 lon_rad = std::atan2(src.source_y, src.source_x);
			****************************************************************/
			}
			break;

		case REGION_3:
		case REGION_4:

			{

				/* correct g by using it as a first guess into the bowring formula*/
				top = z * ga * (1 - ell.f);
				top2 = top * top;

				rr = top2 + w2;
				q = std::sqrt(rr);

				{
					double sn;
					double cn;
					double s3;
					double c3;
					sn = top / q;
					cn = w / q;
					s3 = CB(sn);
					c3 = CB(cn);
					top = z + ell.b * ell.ep2 * s3;
					bot = w - ell.e2 * ell.a * c3;
				}
				top2 = top * top;

				rr = top2 + bot * bot;
				q = std::sqrt(rr);

#ifdef USE_CUSTOM_HT
				s12 = top2 / rr;

				/*Fast Rn computation*/
				{
					double Rn;
					COMPUTE_RN_FAST(s12, Rn);
					/*Fast Computation of Rn = a/std::sqrt(1-Eps2*sin_squared(latitude))*/

					if (s12 <= 0.5)
					{
						ht = w * q / bot - Rn;
					}
					else
					{
						ht = z * q / top - (1 - ell.e2) * Rn;
					}
				}
#endif

				lat_rad = std::atan2(top, bot);
			}
			break;
		case REGION_5:
			{
				double gee;

#if 0
				/*!\note This check has been replaced by a higer level check stipulating
				that we are at no lesser altitude than negative b.
				*/

				/*
				 * we are not allowed to be within 50 Kilometers of
				 * the center of the earth
				 */

				if (w2 + z2 < 50000
				   )
				{
#ifdef CONV_DEBUG
					fprintf(f_ptr, "[SEDRIS CONVERSIONS] Point within 50 "\
							"kilometers of the earth's center.  Outside "\
							"SEDRIS requirement range\n");
#endif
					return 32;
				}
#endif

				cf = 54 * ell.b2 * z2;
				gee = w2 + (1 - ell.e2) * z2 - ell.a2 * e4;
				{
					alpha = cf / (gee * gee);
					cl = e4 * w2 * alpha / gee;
				}
				arg2 = cl * (cl + 2);
				s1 = 1 + cl + std::sqrt(arg2);

				//s = pow(s1, .333333333333333333333333333333);
				s = std::cbrt(s1);

				{
					auto temp = s / (s * s + 1 + s);
					p = alpha * temp * temp / 3;
				}

				xarg = 1 + (2 * e4 * p);
				q = std::sqrt(xarg);

				{
					r_2 = -p * (2 * (1 - ell.e2) * z2 / (q * (1 + q)) + w2);
					r_1 = 1 + 1 / q;
					r_2 /= ell.a2;
					/*
					 * DUE TO PRECISION ERRORS THE ARGUMENT MAY BECOME
					 * NEGATIVE. IF SO SET THE ARGUMENT TO ZERO.
					 */
					if (r_1 + r_2 > 0)
					{
						ro = ell.a * std::sqrt(0.5 * (r_1 + r_2));
					}
					else
					{
						ro = 0;
					}
					ro -= p * ell.e2 * w / (1 + q);
				}

				roe = ell.e2 * ro;
				arg = SQ(w - roe) + z2;
				v   = std::sqrt(arg - ell.e2 * z2);
				{
					zo = ell.a * (1 - ell.e2) * z / v;
#ifdef USE_CUSTOM_HT
					ht = std::sqrt(arg) * (1 - ell.a * (1 - ell.e2) / v);
#endif
				}

				top = z + ell.ep2 * zo;
				lat_rad = std::atan2(top, w);

				/************************************************************

				As an optimization, this is done below
				 lon_rad = std::atan2(src.source_y, src.source_x);
				************************************************************/

				/* end of Exact solution */
			}
			break;
		case REGION_SPHERICAL:
			/*This case gets out of the switch statement and lets the longitude
			 *be computed at the end of the routine with a single return statement
			 *Conclusion:  Don't Remove this case or there won't be a good way
			 *to skip the region test and this won't work very well.
			 */
			break;

		default:
			//status = 32;
#if 0
			fprintf(stderr,"GC to GD: region error, region = %i line %d \n",region,__LINE__);
#endif
			break;
		}

		/*Since logitude calculation is common to all regions and requires lots of time,
		we'll just do it here to break interlock and hopefully have it occur
		during some of the function return bureaucracy*/

		//lon_rad = std::atan2(y, x);
	}

	//return status;

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
#undef COMPUTE_RN_FAST
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "Cameron Kellough?",
	/*.code_copyright              =*/ "SEDRIS",
	/*.license                     =*/ "\"This software was developed for use by the United States Government with unlimited rights.\"",
	/*.orig_impl_lang              =*/ "C++",
	/*.url                         =*/ "http://www.sedris.org/sdk_4.0.1/sdk_c.htm",
	/*.citation                    =*/ R"(sedris_c_sdk_4.1.0
file:src/lib/srm/src/lib/cpp_impl/srf_impl/srm_gc.cpp)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace shu_2010_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto tmp = ell.a2 * z2 + ell.b2 * w2;
	auto k = (w2 + z2) * (std::sqrt(tmp) - ell.a * ell.b) / tmp;

	for (int i = 1; i <= max_iterations; ++i)
	{
		k -= shu_2010_delta_k(w2, z2, k);
	}

	const auto p = ell.a + ell.b * k;
	const auto q = ell.b + ell.a * k;

	auto sin_lat = z * p;
	auto cos_lat = w * q * (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = k * hypot(w * ell.b / p, z * ell.a / q);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_shu_2010_delta_k;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -7,
	/*.algo_author                 =*/ "C. Shu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.academia.edu/10319914/An_iterative_algorithm_to_compute_geodetic_coordinates",
	/*.citation                    =*/ R"(Chanfang Shu, Fei Li, An iterative algorithm to compute geodetic coordinates, Computers & Geosciences, Volume 36, Issue 9, 2010, Pages 1145-1149, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2010.02.004.)"
);

}
// }}}

namespace shu_2010_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto tmp = ell.a2 * z2 + ell.b2 * w2;
	auto k = (w2 + z2) * (std::sqrt(tmp) - ell.a * ell.b) / tmp;

	for (int i = 1; i <= max_iterations; ++i)
	{
		k -= shu_2010_delta_k(w2, z2, k);
	}

	const auto p = ell.a + ell.b * k;
	const auto q = ell.b + ell.a * k;

	auto sin_lat = z * p;
	auto cos_lat = w * q * (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = k * hypot(w * ell.b / p, z * ell.a / q);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_shu_2010_delta_k;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "C. Shu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.academia.edu/10319914/An_iterative_algorithm_to_compute_geodetic_coordinates",
	/*.citation                    =*/ R"(Chanfang Shu, Fei Li, An iterative algorithm to compute geodetic coordinates, Computers & Geosciences, Volume 36, Issue 9, 2010, Pages 1145-1149, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2010.02.004.)"
);

}
// }}}

namespace shu_2010_customht_1
// {{{
{
#define USE_CUSTOM_HT

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto tmp = ell.a2 * z2 + ell.b2 * w2;
	auto k = (w2 + z2) * (std::sqrt(tmp) - ell.a * ell.b) / tmp;

	for (int i = 1; i <= max_iterations; ++i)
	{
		k -= shu_2010_delta_k(w2, z2, k);
	}

	const auto p = ell.a + ell.b * k;
	const auto q = ell.b + ell.a * k;

	auto sin_lat = z * p;
	auto cos_lat = w * q * (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = k * hypot(w * ell.b / p, z * ell.a / q);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_shu_2010_delta_k;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -5,
	/*.algo_author                 =*/ "C. Shu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.academia.edu/10319914/An_iterative_algorithm_to_compute_geodetic_coordinates",
	/*.citation                    =*/ R"(Chanfang Shu, Fei Li, An iterative algorithm to compute geodetic coordinates, Computers & Geosciences, Volume 36, Issue 9, 2010, Pages 1145-1149, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2010.02.004.)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace shu_2010_customht_2
// {{{
{
#define USE_CUSTOM_HT

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto tmp = ell.a2 * z2 + ell.b2 * w2;
	auto k = (w2 + z2) * (std::sqrt(tmp) - ell.a * ell.b) / tmp;

	for (int i = 1; i <= max_iterations; ++i)
	{
		k -= shu_2010_delta_k(w2, z2, k);
	}

	const auto p = ell.a + ell.b * k;
	const auto q = ell.b + ell.a * k;

	auto sin_lat = z * p;
	auto cos_lat = w * q * (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);

#ifdef USE_CUSTOM_HT
	ht = k * hypot(w * ell.b / p, z * ell.a / q);
#else
	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_shu_2010_delta_k;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "C. Shu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.academia.edu/10319914/An_iterative_algorithm_to_compute_geodetic_coordinates",
	/*.citation                    =*/ R"(Chanfang Shu, Fei Li, An iterative algorithm to compute geodetic coordinates, Computers & Geosciences, Volume 36, Issue 9, 2010, Pages 1145-1149, ISSN 0098-3004, https://doi.org/10.1016/j.cageo.2010.02.004.)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace sofair_1993
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// b2 / (a2 - b2) == 1 / ep2
	const auto p = 1 / ell.ep2;
	const auto q = ell.b2 * (ell.a2 * w2 + ell.b2 * z2 - SQ(ell.a2 - ell.b2)) / (SQ(ell.a2 - ell.b2) * z2);
	//const auto r = -SQ(ell.b2) / ((ell.a2 - ell.b2) * z2);
	// SDW: NaNs at the equator
	const auto s = -POW6(ell.b) / (SQ(ell.a2 - ell.b2) * z2);

	const auto f = -SQ(q) / 12;
	const auto h = -CB(q) / 108 - 0.5 * ell.a2 * POW4(ell.b2) * w2 / (POW4(ell.a2 - ell.b2) * POW4(z));

	const auto g = SQ(h) / 4 + CB(f) / 27;

	const auto u = std::cbrt(std::sqrt(g) - h / 2);
	const auto v = std::cbrt(-std::sqrt(g) - h / 2);

	const auto t_p = u + v;
	const auto t = t_p + q / 6;

	const auto c = std::sqrt(p * p - q + 2 * t);
	const auto d = std::sqrt(t * t - s);

	const auto k = (-(p - c) + std::sqrt((p - c) * (p - c) - 4 * (t - d))) / 2;

	auto z_ = z * k;

	auto w_ = ell.a * std::sqrt(1 - SQ(z_ / ell.b));

	const auto Ne = std::sqrt(POW4(ell.a) * SQ(z_) + POW4(ell.b) * SQ(w_)) / ell.b2;

	auto sin_lat = (ell.a2 / ell.b2) * z_ / Ne;

	// clamp to valid range
	if (sin_lat < -1) { sin_lat = -1; }
	else if (sin_lat > 1) { sin_lat = 1; }

	lat_rad = std::asin(sin_lat);

	//ht = ell.get_ht(w, z, lat_rad);
	const auto cos_lat = cos_from_sin(sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "I. Sofair",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "http://www.dtic.mil/dtic/tr/fulltext/u2/a265783.pdf",
	/*.citation                    =*/ R"(Sofair, I. (1985, revised 1993) A Method for Calculating Exact Geodetic Latitude and Altitude - NSWC TR 85-85 AD-A265 783)"
);

}
// }}}

namespace sofair_2000
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto p = std::abs(z) / ell.ep2;
	const auto s = w2 / (ell.e2 * ell.ep2);
	//const auto q = p * p - ell.b2 + w2 / (ell.e2 * ell.ep2);
	const auto q = p * p - ell.b2 + s;

	const auto u = p / std::sqrt(q);
	const auto v = ell.b2 * u * u / q;

	//const auto P = 27 * ell.a2 * (w2) * z2 / (ell.ep2 * ell.ep2 * ell.ep2 * ell.ep2);
	const auto P = 27 * v * s / q;
	const auto Q_ = std::sqrt(P + 1) + std::sqrt(P);
	//const auto Q = std::pow(std::sqrt(P + q * q * q) + std::sqrt(P), 2.0/3.0);
	//const auto Q = std::pow(Q_, 2.0/3.0);
	const auto Q = std::cbrt(Q_ * Q_);
	//const auto Qp = std::pow(std::sqrt(P + q * q * q) - std::sqrt(P), 2.0/3.0);
	//const auto t = (q + Q + Qp) / 6;
	const auto t = (1 + Q + 1 / Q) / 6;
	//const auto c = std::sqrt(p2 - q + 2 * t);
	const auto c = std::sqrt(u * u - 1 + 2 * t);
	const auto w_ = (c - u) / 2;

	//const auto Z = std::copysign((c - p + std::sqrt(2 * p2 - q - 2 * t - 2 * p * c + 4 * std::sqrt(t * t + p2 * b2))) / 2, z);
	//const auto Z = std::copysign(std::sqrt(q), z) * (w_ + std::sqrt(std::sqrt(t * t + v) - u * w_ - t / 2 - 0.25));
	auto Z = std::sqrt(q) * (w_ + std::sqrt(std::sqrt(t * t + v) - u * w_ - t / 2 - 0.25));
	if (z < 0)
		Z = -Z;

	const auto Ne = ell.a * std::sqrt(1 + Z * Z * ell.ep2 / ell.b2);

	auto sin_lat = (ell.ep2 + 1) * (Z / Ne);

	// clamp to valid range
	if (sin_lat < -1) { sin_lat = -1; }
	else if (sin_lat > 1) { sin_lat = 1; }

	lat_rad = std::asin(sin_lat);

	//ht = ell.get_ht(w, z, lat_rad);
	const auto cos_lat = cos_from_sin(sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 0,
	/*.algo_author                 =*/ "I. Sofair",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/246433285_Improved_Method_for_Calculating_Exact_Geodetic_Latitude_and_Altitude_Revisited",
	/*.citation                    =*/ R"(Sofair, I. (2000). Improved Method for Calculating Exact Geodetic Latitude and Altitude Revisited. Journal of Guidance Control and Dynamics - J GUID CONTROL DYNAM. 23. 369-369. 10.2514/2.4534.)"
);

}
// }}}

namespace sudano_1997
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	const auto w4 = w2 * w2;
	const auto z4 = z2 * z2;

	const auto a2 = ell.a2;
	const auto e2 = ell.e2;
	const auto a4 = a2 * a2;

	const auto a_0 = z4;
	const auto a_1 = -2 * ell.a2 * SQ(e2) * z2 - 2 * w2 * z2 - 2 * z4 - 2 * e2 * z4;
	const auto a_2 = a4 * POW4(e2) - 2 * a2 * e2 * w2 + w4 + 4 * a2 * SQ(e2) * z2 + 2 * a2 * CB(e2) * z2
		+ 2 * w2 * z2 + 4 * e2 * w2 * z2 + z4 + 4 * e2 * z4 + SQ(e2) * z4;
	const auto a_3 = -2 * a4 * POW4(e2) + 2 * a2 * SQ(e2) * w2 + 2 * a2 * CB(e2) * w2 - 2 * e2 * w4
		-2 * a2 * SQ(e2) * z2 - 4 * a2 * CB(e2) * z2 - 4 * e2 * w2 * z2 - 2 * SQ(e2) * w2 * z2
		-2 * e2 * z4 - 2 * SQ(e2) * z4;
	const auto a_4 = a4 * POW4(e2) - 2 * a2 * CB(e2) * w2 + SQ(e2) * w4 + 2 * a2 * CB(e2) * z2
		+ 2 * SQ(e2) * w2 * z2 + SQ(e2) * z4;

	const auto b_1 = 2 * CB(a_2) - 9 * a_1 * a_2 * a_3 + 27 * a_0 * SQ(a_3) + 27 * SQ(a_1) * a_4 - 72 * a_0 * a_2 * a_4;
	const auto b_2 = SQ(a_2) - 3 * a_1 * a_3 + 12 * a_0 * a_4;
	const auto b_3 = std::cbrt(b_1 + std::sqrt(SQ(b_1) - 4 * CB(b_2)));
	const auto b_4 = std::sqrt(
			SQ(a_3) / (4 * SQ(a_4))
			- 2 * a_2 / (3 * a_4)
			+ std::cbrt(2 * b_2) / (3 * a_4 * b_3)
			+ b_3 / (3 * std::cbrt(2 * a_4))
			);
	const auto b_5 = SQ(a_3) / (2 * SQ(a_4))
		- 4 * a_2 / (3 * a_4)
		+ std::cbrt(2 * b_2) / (3 * a_4 * b_3)
		- b_3 / (3 * std::cbrt(2 * a_4))
		+ (CB(a_3) - 4 * a_2 * a_3 * a_4 + 8 * a_1 * SQ(a_4)) / (4 * CB(a_4) * b_4)
		;

	const auto s2 = std::sqrt(b_5) / 2 - a_3 / (4 * a_4) - b_4 / 2;
	auto sin_lat = std::sqrt(s2);

	if (z < 0)
		sin_lat = -sin_lat;

	lat_rad = std::asin(sin_lat);

	//ht = ell.get_ht(w, z, lat_rad);
	const auto cos_lat = cos_from_sin(sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
/*
#ifdef USE_CUSTOM_HT
	const auto Rn = ell.get_Rn(sin_lat);
	ht = (std::abs(z) + w * Rn * std::cos(lat_rad + (1 - ell.e2)) * std::abs(sin_lat)) / (std::cos(lat_rad) + std::abs(sin_lat));
#else
#endif
*/
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "J. Sudano",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude",
	/*.citation                    =*/ R"(Sudano, John. (1997). An exact conversion from an Earth-centered coordinate system to latitude, longitude and altitude. 646 - 650 vol.2. 10.1109/NAECON.1997.622711.)"
);

}
// }}}

namespace turner_2013
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// p == a/b-1
	// a/b == 1/(1-f)
	const auto p = 1 / (1 - ell.f) - 1;
	ht = 0;
	double lat_c = 0;

	const auto r = std::sqrt(w2 + z2);
	const auto h0 = r - ell.b;
	// SDW: these do not work
	//const auto lat0 = std::atan2(z, r - w);
	const auto lat0 = 2 * std::atan2(r - w, z);

	ht += h0;
	lat_c += lat0;

	const auto cos_lat0 = std::cos(lat0);
	const auto cos2_lat0 = cos_lat0 * cos_lat0;
	//const auto cos3_lat0 = CB(cos_lat0);
	const auto cos4_lat0 = POW4(cos_lat0);
	//const auto cos5_lat0 = POW5(cos_lat0);
	const auto cos6_lat0 = POW6(cos_lat0);

	const auto h0pb = h0 + ell.b;
	const auto h0mb = h0 - ell.b;
	//const auto h0p3b = h0 + 3 * ell.b;
	const auto h0m3b = h0 - 3 * ell.b;

	const auto sin_2lat0 = std::sin(2 * lat0);
	const auto sin2_2lat0 = sin_2lat0 * sin_2lat0;

	const auto C1 = -4 * CB(h0) + 37 * CB(ell.b) - 66 * ell.b2 * h0 + 33 * ell.b * h0 * h0;
	const auto C2 = CB(h0) - 31 * CB(ell.b) + 75 * ell.b2 * h0 - 33 * ell.b * h0 * h0;
	const auto C3 = 3 * ell.b * (ell.b2 + 3 * h0 * h0 - 6 * ell.b * h0);
	const auto C4 = 139 * CB(ell.b) - 5 * CB(h0) + 49 * ell.b * h0 * h0 - 143 * ell.b2 * h0;
	const auto C5 = CB(h0) - 127 * CB(ell.b) + 163 * ell.b2 * h0 - 45 * ell.b * h0 * h0;
	const auto C6 = 4 * ell.b * (-10 * ell.b * h0 + h0 * h0 + 5 * ell.b2);
	const auto C7 = 4 * POW4(h0) + 118 * POW4(ell.b) + 198 * ell.b2 * h0 * h0 - 266 * CB(ell.b) * h0 - 54 * ell.b * CB(h0);
	const auto C8 = -155 * POW4(ell.b) - 2 * POW4(h0) + 67 * ell.b * CB(h0) + 421 * CB(ell.b) * h0 - 315 * ell.b2 * h0 * h0;
	const auto C9 = 49 * POW4(ell.b) - 15 * ell.b * CB(h0) - 185 * CB(ell.b) * h0 + 135 * CB(ell.b) * h0 * h0;
	const auto C10 = 2 * ell.b2 * (10 * ell.b * h0 - 5 * h0 * h0 - ell.b2);

	const auto h1 = -ell.b * cos2_lat0;
	const auto lat1 = (ell.b - h0) * sin_2lat0 / (2 * (h0pb));
	ht += p * h1;
	lat_c += p * lat1;

	const auto h2 = -ell.b * (h0m3b) * sin2_2lat0 / (8 * (h0pb));
	const auto lat2 = (h0 * h0 - 4 * ell.b * h0 + 3 * ell.b2) * std::sin(4 * lat0) / (8 * h0pb * h0pb) + sin_2lat0 / 4;
	ht += p * p * h2;
	lat_c += p * p * lat2;

	const auto h3 = ell.b * sin2_2lat0 * (h0m3b * h0m3b * cos2_lat0 + 4 * ell.b * (h0mb)) / (8 * h0pb * h0pb);
	const auto lat3 = sin_2lat0 * (C1 * cos4_lat0 + C2 * cos2_lat0 + C3) / (6 * CB(h0pb));
	ht += CB(p) * h3;
	lat_c += CB(p) * lat3;

	const auto h4 = ell.b * sin2_2lat0 * (C4 * cos4_lat0 + C5 * cos2_lat0 + C6) / (32 * CB(h0pb));
	const auto lat4 = sin_2lat0 * (C7 * cos6_lat0 + C8 * cos4_lat0 + C9 * cos2_lat0 + C10) / (4 * POW4(h0pb));
	ht += POW4(p) * h4;
	lat_c += POW4(p) * lat4;

	// SDW: this is wrong
	//lat_rad = std::atan2((1 - ell.f) * std::cos(lat_c), std::sin(lat_c));

	auto sin_lat = std::sin(lat_c);
	auto cos_lat = std::cos(lat_c);
	cos_lat *= (1 - ell.f);

	lat_rad = std::atan2(sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ 99,
	/*.algo_author                 =*/ "J. Turner, T. Elgohary",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/269038945_A_Simple_Perturbation_Algorithm_for_Inverting_the_Cartesian_to_Geodic_Transformation",
	/*.citation                    =*/ R"(Turner, James & Elgohary, Tarek. (2013). A Simple Perturbation Algorithm for Inverting the Cartesian to Geodic Transformation..)"
);

}
// }}}

namespace vermeille_2004
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto p = w2 / ell.a2;
	const auto q = (1 - ell.e2) * z2 / ell.a2;
	const auto r = (p + q - e4) / 6;
	const auto r3 = CB(r);

	const auto s = 0.25 * e4 * p * q / r3;
	const auto t = std::cbrt(1 + s + std::sqrt(s * (2 + s)));

	const auto u = r * (1 + t + 1 / t);
	const auto v = std::sqrt(u * u + e4 * q);
	const auto w_ = 0.5 * ell.e2 * (u + v - q) / v;
	const auto k = std::sqrt(u + v + w_ * w_) - w_;

	const auto D = k * w / (k + ell.e2);
	const auto tmp = hypot(D, z);

	lat_rad = 2 * std::atan2(z, D + tmp);

	ht = ell.get_ht(w, z, lat_rad);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "H. Vermeille",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/s00190-004-0375-4",
	/*.citation                    =*/ R"(Vermeille, H. Journal of Geodesy (2004) 78: 94. https://doi.org/10.1007/s00190-004-0375-4)"
);

}
// }}}

namespace vermeille_2011
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto p = w2 / ell.a2;
	const auto q = (z2 / ell.a2) * (1 - ell.e2);
	const auto r = (p + q - e4) / 6;
	const auto r3 = CB(r);

	const auto e4pq = e4 * p * q;
	//const auto sqrt_e4pq = ell.e2 * std::sqrt(p * q);
	const auto sqrt_e4pq = ell.e2 * (1 - ell.f) * std::abs(z) * w / ell.a2;

	const auto tmp1 = 8 * r3 + e4pq;

	// outside the evolute
	if (tmp1 > 0)
	{
		const auto tmp2 = std::sqrt(tmp1);

		const auto u = r + 0.5 * std::cbrt(SQ(tmp2 + sqrt_e4pq)) + 0.5 * std::cbrt(SQ(tmp2 - sqrt_e4pq));

		const auto v = std::sqrt(u * u + e4 * q);
		const auto w_ = 0.5 * ell.e2 * ((u + v) - q) / v;
		const auto k = (u + v) / (std::sqrt(w_ * w_ + u + v) + w_);
		const auto D = k * w / (k + ell.e2);

		const auto tmp3 = hypot(D, z);
		lat_rad = 2 * std::atan2(z, tmp3 + D);
#ifdef USE_CUSTOM_HT
		//ht = (k + ell.e2 - 1) * tmp3 / k;
		ht = (k - (1 - ell.e2)) * tmp3 / k;
#endif
	}
	// inside the evolute and not in the singular disc
	else // tmp1 <= 0
	{
		if (q != 0)
		{
			const auto tmp2_neg = std::sqrt(-tmp1);

			const auto tmp4 = 2 * std::atan2(sqrt_e4pq, tmp2_neg + std::sqrt(-8 * r3)) / 3;
			//const auto u = -4 * r * std::sin(tmp4) * std::cos(M_PI / 6 + tmp4);
			/*
			https://www.wolframalpha.com/input/?i=-4+*+sin(x)+*+cos(pi%2F6+%2B+x)
			-4 * sin(x) * cos(π/6 + x) == 1 - 2 * sin(2*x + π/6)
			*/
			const auto u = r * (1 - 2 * sin(2 * tmp4 + M_PI / 6));

			const auto v = std::sqrt(u * u + e4 * q);
			const auto w_ = 0.5 * ell.e2 * ((u + v) - q) / v;
			const auto k = (u + v) / (std::sqrt(w_ * w_ + (u + v)) + w_);
			const auto D = k * w / (k + ell.e2);

			const auto tmp3 = hypot(D, z);
			lat_rad = 2 * std::atan2(z, tmp3 + D);
#ifdef USE_CUSTOM_HT
			//ht = (k + ell.e2 - 1) * tmp3 / k;
			ht = (k - (1 - ell.e2)) * tmp3 / k;
#endif
		}
		// in the singular disc
		else // if (q == 0 && p <= e2*e2)
		{
			const auto tmp3 = std::sqrt(ell.e2 - p);
			lat_rad = 2 * std::atan2(std::sqrt(e4 - p), (ell.e * tmp3 + (1 - ell.f) * std::sqrt(p)));
#ifdef USE_CUSTOM_HT
			ht = -ell.a * (1 - ell.f) * tmp3 / ell.e;
#endif
		}
	}

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "H. Vermeille",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/s00190-010-0419-x",
	/*.citation                    =*/ R"(Vermeille, H. J Geod (2011) 85: 105. https://doi.org/10.1007/s00190-010-0419-x)"
);

}
// }}}

namespace vermeille_2011_customht
// {{{
{
#define USE_CUSTOM_HT

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto e4 = ell.e2 * ell.e2;

	const auto p = w2 / ell.a2;
	const auto q = (z2 / ell.a2) * (1 - ell.e2);
	const auto r = (p + q - e4) / 6;
	const auto r3 = CB(r);

	const auto e4pq = e4 * p * q;
	//const auto sqrt_e4pq = ell.e2 * std::sqrt(p * q);
	const auto sqrt_e4pq = ell.e2 * (1 - ell.f) * std::abs(z) * w / ell.a2;

	const auto tmp1 = 8 * r3 + e4pq;

	// outside the evolute
	if (tmp1 > 0)
	{
		const auto tmp2 = std::sqrt(tmp1);

		const auto u = r + 0.5 * std::cbrt(SQ(tmp2 + sqrt_e4pq)) + 0.5 * std::cbrt(SQ(tmp2 - sqrt_e4pq));

		const auto v = std::sqrt(u * u + e4 * q);
		const auto w_ = 0.5 * ell.e2 * ((u + v) - q) / v;
		const auto k = (u + v) / (std::sqrt(w_ * w_ + u + v) + w_);
		const auto D = k * w / (k + ell.e2);

		const auto tmp3 = hypot(D, z);
		lat_rad = 2 * std::atan2(z, tmp3 + D);
#ifdef USE_CUSTOM_HT
		//ht = (k + ell.e2 - 1) * tmp3 / k;
		ht = (k - (1 - ell.e2)) * tmp3 / k;
#endif
	}
	// inside the evolute and not in the singular disc
	else // tmp1 <= 0
	{
		if (q != 0)
		{
			const auto tmp2_neg = std::sqrt(-tmp1);

			const auto tmp4 = 2 * std::atan2(sqrt_e4pq, tmp2_neg + std::sqrt(-8 * r3)) / 3;
			//const auto u = -4 * r * std::sin(tmp4) * std::cos(M_PI / 6 + tmp4);
			/*
			https://www.wolframalpha.com/input/?i=-4+*+sin(x)+*+cos(pi%2F6+%2B+x)
			-4 * sin(x) * cos(π/6 + x) == 1 - 2 * sin(2*x + π/6)
			*/
			const auto u = r * (1 - 2 * sin(2 * tmp4 + M_PI / 6));

			const auto v = std::sqrt(u * u + e4 * q);
			const auto w_ = 0.5 * ell.e2 * ((u + v) - q) / v;
			const auto k = (u + v) / (std::sqrt(w_ * w_ + (u + v)) + w_);
			const auto D = k * w / (k + ell.e2);

			const auto tmp3 = hypot(D, z);
			lat_rad = 2 * std::atan2(z, tmp3 + D);
#ifdef USE_CUSTOM_HT
			//ht = (k + ell.e2 - 1) * tmp3 / k;
			ht = (k - (1 - ell.e2)) * tmp3 / k;
#endif
		}
		// in the singular disc
		else // if (q == 0 && p <= e2*e2)
		{
			const auto tmp3 = std::sqrt(ell.e2 - p);
			lat_rad = 2 * std::atan2(std::sqrt(e4 - p), (ell.e * tmp3 + (1 - ell.f) * std::sqrt(p)));
#ifdef USE_CUSTOM_HT
			ht = -ell.a * (1 - ell.f) * tmp3 / ell.e;
#endif
		}
	}

#ifdef USE_CUSTOM_HT
#else
	ht = ell.get_ht(w, z, lat_rad);
#endif
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "H. Vermeille",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://link.springer.com/article/10.1007/s00190-010-0419-x",
	/*.citation                    =*/ R"(Vermeille, H. J Geod (2011) 85: 105. https://doi.org/10.1007/s00190-010-0419-x)"
);

#undef USE_CUSTOM_HT
}
// }}}

namespace wu_2003_1
// {{{
{

constexpr int max_iterations = 1;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS_CHECKED

	const auto A = ell.b * z;
	const auto B = 2 * ell.a * w;
	constexpr auto C = 2 * ell.a2 * ell.e2;

	// SDW: NaNs at the poles
	//const auto tan_lat_0 = z / ((1 - ell.f) * w);

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.f);

	/*
	** +- sqrt(1+cot^2) == 1/sin
	** -1/tan + 1/sin = tan(x/2)
	** 2*atan(tan(x/2)) == x
	*/
	normalize(cos_lat, sin_lat);

	// SDW: NaNs at the equator
	//auto t = -1 / tan_lat_0 + 1 / sin_lat;
	auto t = (1 - cos_lat) / sin_lat;

	for (int i = 1; i <= max_iterations; ++i)
	{
		t -= wu_2003_delta_t(A, B, C, t);
	}

	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
	sin_lat = 2 * t;
	cos_lat = (1 - t * t) * (1 - ell.f);

	// 2 * tan(x/2) / (1 - tan(x/2)^2) == tan(x)
	// https://www.wolframalpha.com/input/?i=2+*+tan(x%2F2)+%2F+(1+-+tan(x%2F2)%5E2),+tan(x)

	lat_rad = std::atan2(sin_lat, cos_lat);

	// SDW: This is not accurate.
	//auto tan_lat = 2 * t / ((1 - ell.f) * (1 - t * t));
	// std::sqrt(1 + 1 / (tan_lat * tan_lat)) == 1 / std::abs(sin_lat);
	//ht = std::copysign(z - ell.b * 2 * t / (1 + t * t), lat_rad) * std::sqrt(1 + 1 / (tan_lat * tan_lat));
	// SDW: This is not accurate.
	//ht = std::copysign(z - ell.b * 2 * t / (1 + t * t), lat_rad) / std::abs(sin_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_wu_2003_delta_t;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ true,
	/*.ilog10_mean_dist_err        =*/ 2,
	/*.algo_author                 =*/ "Y. Wu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/3003672_Algorithm_of_Earth-Centered_Earth-Fixed_Coordinates_to_Geodetic_Coordinates",
	/*.citation                    =*/ R"(https://ieeexplore.ieee.org/document/1261144/
Wu, Yuanxin & Wang, Ping & Hu, Xiaoping. (2003). Algorithm of Earth-Centered Earth-Fixed Coordinates to Geodetic Coordinates. Aerospace and Electronic Systems, IEEE Transactions on. 39. 1457 - 1461. 10.1109/TAES.2003.1261144.)"
);

}
// }}}

namespace wu_2003_2
// {{{
{

constexpr int max_iterations = 2;

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS_CHECKED

	const auto A = ell.b * z;
	const auto B = 2 * ell.a * w;
	constexpr auto C = 2 * ell.a2 * ell.e2;

	// SDW: NaNs at the poles
	//const auto tan_lat_0 = z / ((1 - ell.f) * w);

	auto sin_lat = z;
	auto cos_lat = w * (1 - ell.f);

	/*
	** +- sqrt(1+cot^2) == 1/sin
	** -1/tan + 1/sin = tan(x/2)
	** 2*atan(tan(x/2)) == x
	*/
	normalize(cos_lat, sin_lat);

	// SDW: NaNs at the equator
	//auto t = -1 / tan_lat_0 + 1 / sin_lat;
	auto t = (1 - cos_lat) / sin_lat;

	for (int i = 1; i <= max_iterations; ++i)
	{
		t -= wu_2003_delta_t(A, B, C, t);
	}

	// https://en.wikipedia.org/wiki/Tangent_half-angle_formula
	sin_lat = 2 * t;
	cos_lat = (1 - t * t) * (1 - ell.f);

	// 2 * tan(x/2) / (1 - tan(x/2)^2) == tan(x)
	// https://www.wolframalpha.com/input/?i=2+*+tan(x%2F2)+%2F+(1+-+tan(x%2F2)%5E2),+tan(x)

	lat_rad = std::atan2(sin_lat, cos_lat);

	// SDW: This is not accurate.
	//auto tan_lat = 2 * t / ((1 - ell.f) * (1 - t * t));
	// std::sqrt(1 + 1 / (tan_lat * tan_lat)) == 1 / std::abs(sin_lat);
	//ht = std::copysign(z - ell.b * 2 * t / (1 + t * t), lat_rad) * std::sqrt(1 + 1 / (tan_lat * tan_lat));
	// SDW: This is not accurate.
	//ht = std::copysign(z - ell.b * 2 * t / (1 + t * t), lat_rad) / std::abs(sin_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = _lines_wu_2003_delta_t;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ true,
	/*.ilog10_mean_dist_err        =*/ -4,
	/*.algo_author                 =*/ "Y. Wu, et al.",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/3003672_Algorithm_of_Earth-Centered_Earth-Fixed_Coordinates_to_Geodetic_Coordinates",
	/*.citation                    =*/ R"(https://ieeexplore.ieee.org/document/1261144/
Wu, Yuanxin & Wang, Ping & Hu, Xiaoping. (2003). Algorithm of Earth-Centered Earth-Fixed Coordinates to Geodetic Coordinates. Aerospace and Electronic Systems, IEEE Transactions on. 39. 1457 - 1461. 10.1109/TAES.2003.1261144.)"
);

}
// }}}

namespace zhang_2005
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	// Eq.(9);
	const auto alpha = w2 / ell.a2 + z2 / ell.b2;
	const auto beta = w2 / ell.b2 + z2 / ell.a2;
	//const auto gamma = ell.b / ell.a + ell.a / ell.b;
	constexpr auto gamma = (1 - ell.f) + 1 / (1 - ell.f);
	constexpr auto gamma2 = gamma * gamma;
	constexpr auto gamma4 = gamma2 * gamma2;

	// Eq.(11);
	const auto p = 2 - beta - gamma2 / 2;
	const auto q = beta * gamma - (2 / gamma) * (alpha + beta);
	const auto r = 1 + beta - (gamma2 / 4) * (2 + beta) + gamma4 / 16;

	//Algebra prediction begain
	//Eq.(20)
	const auto G = (4 + beta - gamma2) / 12;
	const auto G2 = G * G;
	const auto G3 = G * G * G;

	const auto es = -((alpha + beta) * (alpha + beta) / gamma2 - alpha * beta) / 32;
	//Eq.(19)
	const auto E = -G2;
	const auto F = G3 + es;
	//Eq.(21)
	auto E3F2 = es * (es + 2 * G3);
	// theoretically, E3F2>=0; but sometimes E3F2 numerical value is negative
	// in multivalue region.
	E3F2 = std::abs(E3F2);
	// theoretically, F+std::sqrt(E3F2)>=0; but sometimes F+std::sqrt(E3F2) numerical value is negative
	auto m = std::cbrt(std::abs(F + std::sqrt(E3F2)));
	// theoretically, F-std::sqrt(E3F2)>=0; but sometimes F-std::sqrt(E3F2) numerical value is negative
	m = m - E/m;
	// theoretically, -m-p/6>=0; but sometimes -m-p/6 numerical value is negative
	const auto s0 = std::abs(-m - p / 6);
	// theoretically, u+v=std::sqrt((p/2+2*s0)^2-r)-p/2-s0 >=0;
	const auto uv = std::sqrt(std::abs(std::sqrt(std::abs((p / 2 + 2 * s0) * (p / 2 + 2 * s0) - r)) - p / 2 - s0));
	// theoretically, (p/2+2*s0)^2-r>=0;
	//auto t0 = -sign(q) * std::sqrt(s0) + uv;
	auto t0 = -std::copysign(std::sqrt(s0), q) + uv;
	//Algebra prediction end

	// sometimes the Algebra prediction t exists numerical value error, but we
	// can see this t as a good initial value and get a new solution from the
	// basic quartic equation untill the process converges satisfactorily.
	// The author find that only once iteration, the new solution is precision
	// and satisfy our application to 1 mm over geodetic height range from
	// -6x10^6 to 10^10 m.

	//Algebra correction begin
	t0 = std::sqrt(-p / 2 + std::sqrt(p * p - 4 * (q * t0 + r)) / 2);
	//Algebra correction end

	//eigenvalue or eigen-root that we publish
	const auto lmt = t0 - gamma / 2;

	const auto w0 = w * (ell.a / (ell.a + ell.b * lmt));
	const auto z0 = z * (ell.b / (ell.b + ell.a * lmt));

	auto sin_lat = z0;
	auto cos_lat = w0 * (1 - ell.e2);

	lat_rad = std::atan2(sin_lat, cos_lat); //[-π/2,+π/2]

	// SDW: this is not accurate enough
	//ht = sign(lmt) * hypot(x - x0, y - y0, z - z0);
	//ht = std::copysign(hypot(x - x0, y - y0, z - z0), lmt);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "CD. Zhang, et al.",
	/*.code_copyright              =*/ "CD. Zhang, et al.",
	/*.license                     =*/ "Unknown",
	/*.orig_impl_lang              =*/ "MATLAB",
	/*.url                         =*/ "https://link.springer.com/article/10.1007%2Fs00190-005-0487-5",
	/*.citation                    =*/ R"(Zhang, CD., Hsu, H., Wu, X. et al. J Geodesy (2005) 79: 413. https://doi.org/10.1007/s00190-005-0487-5)"
);

}
// }}}

namespace zhu_1993
// {{{
{

constexpr int _line_begin = __LINE__;
void ecef_to_geodetic(const double x, const double y, const double z,
                      double& lat_rad, double& lon_rad, double& ht)
{
COMMON_FIRST_DECLS

	constexpr auto l = ell.e2 / 2;
	constexpr auto l2 = l * l;
	const auto m = w2 / ell.a2;
	//const auto n = SQ(z * (1 - ell.e2) / ell.b);
	//const auto n = z2 * SQ((1 - ell.e2)) / ell.b2;
	//const auto n = (z2 / ell.a2) * (1 - ell.e2); // (1-e2) / a2 == 1/(Rp**2)
	const auto n = (1 - ell.e2) * z2 / ell.a2;

	const auto i = -(2 * l2 + m + n) / 2;
	const auto k = l2 * (l2 - m - n);

	const auto mnl2 = m * n * l2;
	const auto q = CB(m + n - 4 * l2) / 216 + mnl2;

	const auto D = std::sqrt((2 * q - mnl2) * mnl2);
	const auto Beta = i / 3 - std::cbrt(q + D) - std::cbrt(q - D);
	// (Beta - i) / 2 might be negative
	const auto t = std::sqrt(std::sqrt(Beta * Beta - k) - (Beta + i) / 2) - std::copysign(std::sqrt(std::max(0.0, (Beta - i) / 2)), m - n);

	auto sin_lat = z;
	auto cos_lat = w * (t - l) / (t + l);

	lat_rad = std::atan2(sin_lat, cos_lat);

	normalize(cos_lat, sin_lat);
	ht = ell.get_ht(w, z, sin_lat, cos_lat);
}
constexpr int _line_end = __LINE__;

constexpr int _lines_extra = 0;

const auto func_info = func_info_t(
	/*.func_ref                    =*/ ecef_to_geodetic,
	/*.num_lines                   =*/ _line_end - _line_begin + _lines_extra,
	/*.needs_code_for_corner_cases =*/ false,
	/*.ilog10_mean_dist_err        =*/ -9,
	/*.algo_author                 =*/ "Jijie Zhu",
	/*.code_copyright              =*/ "Steven Ward",
	/*.license                     =*/ "OSL-3.0",
	/*.orig_impl_lang              =*/ "None",
	/*.url                         =*/ "https://www.researchgate.net/publication/245432223_Exact_Conversion_of_Earth-Centered_Earth-Fixed_Coordinates_to_Geodetic_Coordinates",
	/*.citation                    =*/ R"(https://arc.aiaa.org/doi/10.2514/3.21016
https://ieeexplore.ieee.org/document/303772/
Zhu, Jijie. (1993). Exact Conversion of Earth-Centered Earth-Fixed Coordinates to Geodetic Coordinates. Journal of Guidance Control and Dynamics - J GUID CONTROL DYNAM. 16. 389-391. 10.2514/3.21016.)"
);

}
// }}}
