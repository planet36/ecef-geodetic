// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// circular angle class
/**
\file
\author Steven Ward
*/

#pragma once

#include "angle_unit.hpp"
#include "isclose.hpp"

/// an angle class where the scalar value and unit of measurement are preserved
template <angle_unit U, std::floating_point T>
class angle
{
    /*
    /// Every angle<U2, T2> is a friend of this.
    template <angle_unit U2, std::floating_point T2>
    friend class angle;
    */

private:
    /// the angle's value (in angle units)
    T value{};

public:
    /// default ctor
    angle() = default;

    /// ctor
    constexpr angle(const T x) : value(x) {}

    /// copy ctor
    angle(const angle&) = default;

    /// dtor
    ~angle() = default;

    /// conversion ctor
    template <angle_unit U2, std::floating_point T2>
    constexpr angle(const angle<U2, T2>& a) :
    value(convert_from<U2, U>(a.scalar()))
    {}

    /// default assignment
    angle& operator=(const angle&) = default;

    constexpr angle& operator+=(const angle& that)
    {
        value += that.value;
        return *this;
    }

    constexpr angle& operator-=(const angle& that)
    {
        value -= that.value;
        return *this;
    }

    constexpr angle& operator*=(const T& x)
    {
        value *= x;
        return *this;
    }

    constexpr angle& operator/=(const T& x)
    {
        value /= x;
        return *this;
    }

    [[nodiscard]] constexpr T to_mrad() const
    {
        return convert_from<U, angle_unit::milliradian>(value);
    }

    [[nodiscard]] constexpr T to_rad() const
    {
        return convert_from<U, angle_unit::radian>(value);
    }

    [[nodiscard]] constexpr T to_rev() const
    {
        return convert_from<U, angle_unit::revolution>(value);
    }

    [[nodiscard]] constexpr T to_deg() const
    {
        return convert_from<U, angle_unit::degree>(value);
    }

    [[nodiscard]] constexpr T to_arcmin() const
    {
        return convert_from<U, angle_unit::arcminute>(value);
    }

    [[nodiscard]] constexpr T to_arcsec() const
    {
        return convert_from<U, angle_unit::arcsecond>(value);
    }

    /// get the scalar value of the angle (in angle units)
    /**
    * This is useful for getting the raw internal value regardless of its unit
    * of measurement.
    */
    [[nodiscard]] constexpr T scalar() const { return value; }

    /// get the unit of measurement of the angle
    [[nodiscard]] constexpr angle_unit units() const { return U; }

    /// convert to a different data type
    template <std::floating_point T2>
    [[nodiscard]] constexpr auto to() const
    {
        return angle<U, T2>{*this};
    }

    /// convert to a different angle unit
    template <angle_unit U2>
    [[nodiscard]] constexpr auto to() const
    {
        return angle<U2, T>{*this};
    }

    /// convert to a different angle unit and data type
    template <angle_unit U2, std::floating_point T2>
    [[nodiscard]] constexpr auto to() const
    {
        return angle<U2, T2>{*this};
    }

    auto operator<=>(const angle&) const = default;
};

/// alias template
template <std::floating_point T>
using ang_mrad = angle<angle_unit::milliradian, T>;

/// alias template
template <std::floating_point T>
using ang_rad = angle<angle_unit::radian, T>;

/// alias template
template <std::floating_point T>
using ang_rev = angle<angle_unit::revolution, T>;

/// alias template
template <std::floating_point T>
using ang_deg = angle<angle_unit::degree, T>;

/// alias template
template <std::floating_point T>
using ang_arcmin = angle<angle_unit::arcminute, T>;

/// alias template
template <std::floating_point T>
using ang_arcsec = angle<angle_unit::arcsecond, T>;

// factory functions

/// create angle<angle_unit::milliradian, T>
template <std::floating_point T>
constexpr auto
make_ang_mrad(const T x)
{
    return ang_mrad<T>{x};
}

/// create angle<angle_unit::radian, T>
template <std::floating_point T>
constexpr auto
make_ang_rad(const T x)
{
    return ang_rad<T>{x};
}

/// create angle<angle_unit::revolution, T>
template <std::floating_point T>
constexpr auto
make_ang_rev(const T x)
{
    return ang_rev<T>{x};
}

/// create angle<angle_unit::degree, T>
template <std::floating_point T>
constexpr auto
make_ang_deg(const T x)
{
    return ang_deg<T>{x};
}

/// create angle<angle_unit::arcminute, T>
template <std::floating_point T>
constexpr auto
make_ang_arcmin(const T x)
{
    return ang_arcmin<T>{x};
}

/// create angle<angle_unit::arcsecond, T>
template <std::floating_point T>
constexpr auto
make_ang_arcsec(const T x)
{
    return ang_arcsec<T>{x};
}

// user-defined literals

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_mrad(const long double x)
{
    using T = long double;
    return ang_mrad<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_mrad(const unsigned long long int x)
{
    using T = long double;
    return ang_mrad<T>{static_cast<T>(x)};
}

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_rad(const long double x)
{
    using T = long double;
    return ang_rad<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_rad(const unsigned long long int x)
{
    using T = long double;
    return ang_rad<T>{static_cast<T>(x)};
}

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_rev(const long double x)
{
    using T = long double;
    return ang_rev<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_rev(const unsigned long long int x)
{
    using T = long double;
    return ang_rev<T>{static_cast<T>(x)};
}

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_deg(const long double x)
{
    using T = long double;
    return ang_deg<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_deg(const unsigned long long int x)
{
    using T = long double;
    return ang_deg<T>{static_cast<T>(x)};
}

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_arcmin(const long double x)
{
    using T = long double;
    return ang_arcmin<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_arcmin(const unsigned long long int x)
{
    using T = long double;
    return ang_arcmin<T>{static_cast<T>(x)};
}

// Note: The C++ Standard requires the input parameter to be 'long double'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_arcsec(const long double x)
{
    using T = long double;
    return ang_arcsec<T>{x};
}

// Note: The C++ Standard requires the input parameter to be 'unsigned long long int'.
/// User-defined literal to create an angle
/**
\sa https://en.cppreference.com/w/cpp/language/user_literal
\param x the given value
\return an angle of the given value
*/
constexpr auto
operator ""_arcsec(const unsigned long long int x)
{
    using T = long double;
    return ang_arcsec<T>{static_cast<T>(x)};
}

/// primary template
template <angle_unit U, std::floating_point T>
struct const_angle;

/// partial specialization for milliradians
template <std::floating_point T>
struct const_angle<angle_unit::milliradian, T>
{
    using angle = ang_mrad<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{mrad_per_rad * std::numbers::pi / 4};
    static constexpr angle quarter_turn{mrad_per_rad * std::numbers::pi / 2};
    static constexpr angle half_turn{mrad_per_rad * std::numbers::pi};
    static constexpr angle full_turn{mrad_per_rad * 2 * std::numbers::pi};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_mrad = const_angle<angle_unit::milliradian, T>;

/// partial specialization for radians
template <std::floating_point T>
struct const_angle<angle_unit::radian, T>
{
    using angle = ang_rad<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{std::numbers::pi / 4};
    static constexpr angle quarter_turn{std::numbers::pi / 2};
    static constexpr angle half_turn{std::numbers::pi};
    static constexpr angle full_turn{2 * std::numbers::pi};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_rad = const_angle<angle_unit::radian, T>;

/// partial specialization for revolutions
template <std::floating_point T>
struct const_angle<angle_unit::revolution, T>
{
    using angle = ang_rev<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{T{0.125}};
    static constexpr angle quarter_turn{T{0.25}};
    static constexpr angle half_turn{T{0.5}};
    static constexpr angle full_turn{T{1}};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_rev = const_angle<angle_unit::revolution, T>;

/// partial specialization for degrees
template <std::floating_point T>
struct const_angle<angle_unit::degree, T>
{
    using angle = ang_deg<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{T{45}};
    static constexpr angle quarter_turn{T{90}};
    static constexpr angle half_turn{T{180}};
    static constexpr angle full_turn{T{360}};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_deg = const_angle<angle_unit::degree, T>;

/// partial specialization for arcminutes
template <std::floating_point T>
struct const_angle<angle_unit::arcminute, T>
{
    using angle = ang_arcmin<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{T{2700}};
    static constexpr angle quarter_turn{T{5400}};
    static constexpr angle half_turn{T{10'800}};
    static constexpr angle full_turn{T{21'600}};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_arcmin = const_angle<angle_unit::arcminute, T>;

/// partial specialization for arcseconds
template <std::floating_point T>
struct const_angle<angle_unit::arcsecond, T>
{
    using angle = ang_arcsec<T>;
    static constexpr angle zero{};
    static constexpr angle eighth_turn{T{162'000}};
    static constexpr angle quarter_turn{T{324'000}};
    static constexpr angle half_turn{T{648'000}};
    static constexpr angle full_turn{T{1'296'000}};
    static constexpr angle inf{T{INFINITY}};
    static constexpr angle nan{T{NAN}};
};

/// alias template
template <std::floating_point T>
using const_ang_arcsec = const_angle<angle_unit::arcsecond, T>;

/// convert to angle from quadrants
template <std::floating_point T>
constexpr auto
convert_from_quadrant(const T x_quadrant)
{
    return ang_rev<T>{x_quadrant / quadrants_per_rev};
}

/// convert to quadrants from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_quadrant(const angle<U, T>& a)
{
    return quadrants_per_rev * a.to_rev();
}

/// convert to angle from sextants
template <std::floating_point T>
constexpr auto
convert_from_sextant(const T x_sextant)
{
    return ang_rev<T>{x_sextant / sextants_per_rev};
}

/// convert to sextants from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_sextant(const angle<U, T>& a)
{
    return sextants_per_rev * a.to_rev();
}

/// convert to angle from octants
template <std::floating_point T>
constexpr auto
convert_from_octant(const T x_octant)
{
    return ang_rev<T>{x_octant / octants_per_rev};
}

/// convert to octants from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_octant(const angle<U, T>& a)
{
    return octants_per_rev * a.to_rev();
}

/// convert to angle from hexacontades
template <std::floating_point T>
constexpr auto
convert_from_hexacontade(const T x_hexacontade)
{
    return ang_rev<T>{x_hexacontade / hexacontades_per_rev};
}

/// convert to hexacontades from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_hexacontade(const angle<U, T>& a)
{
    return hexacontades_per_rev * a.to_rev();
}

/// convert to angle from binary degrees
template <std::floating_point T>
constexpr auto
convert_from_binary_degree(const T x_binary_degree)
{
    return ang_rev<T>{x_binary_degree / binary_degrees_per_rev};
}

/// convert to binary degrees from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_binary_degree(const angle<U, T>& a)
{
    return binary_degrees_per_rev * a.to_rev();
}

/// convert to angle from gradians
template <std::floating_point T>
constexpr auto
convert_from_gradian(const T x_gradian)
{
    return ang_rev<T>{x_gradian / gradians_per_rev};
}

/// convert to gradians from angle
template <angle_unit U, std::floating_point T>
constexpr auto
convert_to_gradian(const angle<U, T>& a)
{
    return gradians_per_rev * a.to_rev();
}

/// unary plus (positive operator)
template <angle_unit U, std::floating_point T>
constexpr auto
operator+(const angle<U, T>& a)
{
    return angle<U, T>{+a.scalar()};
}

/// unary minus (negative operator)
template <angle_unit U, std::floating_point T>
constexpr auto
operator-(const angle<U, T>& a)
{
    return angle<U, T>{-a.scalar()};
}

/// angle<U, T> + angle<U, T>
template <angle_unit U, std::floating_point T>
constexpr auto
operator+(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return angle<U, T>{a1.scalar() + a2.scalar()};
}

/// angle<U, T> - angle<U, T>
template <angle_unit U, std::floating_point T>
constexpr auto
operator-(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return angle<U, T>{a1.scalar() - a2.scalar()};
}

/// angle<U, T> / angle<U, T>
template <angle_unit U, std::floating_point T>
constexpr auto
operator/(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return a1.scalar() / a2.scalar();
}

/// angle<U, T> + angle<U, T2>
template <angle_unit U, std::floating_point T, std::floating_point T2>
constexpr auto
operator+(const angle<U, T>& a1, const angle<U, T2>& a2)
{
    using result_type = typename std::common_type_t<T, T2>;
    return angle<U, result_type>{a1.scalar() + a2.scalar()};
}

/// angle<U, T> - angle<U, T2>
template <angle_unit U, std::floating_point T, std::floating_point T2>
constexpr auto
operator-(const angle<U, T>& a1, const angle<U, T2>& a2)
{
    using result_type = typename std::common_type_t<T, T2>;
    return angle<U, result_type>{a1.scalar() - a2.scalar()};
}

/// angle<U, T> * T2
template <angle_unit U, std::floating_point T, typename T2>
requires std::is_arithmetic_v<T2>
constexpr auto
operator*(const angle<U, T>& a, const T2& x)
{
    using result_type = typename std::common_type_t<T, T2>;
    return angle<U, result_type>{a.scalar() * x};
}

/// T2 * angle<U, T>
template <angle_unit U, std::floating_point T, typename T2>
requires std::is_arithmetic_v<T2>
constexpr auto
operator*(const T2& x, const angle<U, T>& a)
{
    // commutative property
    return a * x;
}

/// angle<U, T> / T2
template <angle_unit U, std::floating_point T, typename T2>
requires std::is_arithmetic_v<T2>
constexpr auto
operator/(const angle<U, T>& a, const T2& x)
{
    using result_type = typename std::common_type_t<T, T2>;
    return angle<U, result_type>{a.scalar() / x};
}

/// angle<U, T> / angle<U, T2>
template <angle_unit U, std::floating_point T, std::floating_point T2>
constexpr auto
operator/(const angle<U, T>& a1, const angle<U, T2>& a2)
{
    return a1.scalar() / a2.scalar();
}

/// angle<U, T> + angle<U2, T2>
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
operator+(const angle<U, T>& a1, const angle<U2, T2>& a2)
{
    return a1 + a2.template to<U>();
}

/// angle<U, T> - angle<U2, T2>
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
operator-(const angle<U, T>& a1, const angle<U2, T2>& a2)
{
    return a1 - a2.template to<U>();
}

/// angle<U, T> / angle<U2, T2>
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
operator/(const angle<U, T>& a1, const angle<U2, T2>& a2)
{
    return a1 / a2.template to<U>();
}

/// is acute?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/AcuteAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_acute(const angle<U, T>& a)
{
    return a < const_angle<U, T>::quarter_turn;
}

/// is right?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/RightAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_right(const angle<U, T>& a)
{
    return a == const_angle<U, T>::quarter_turn;
}

/// is (approximately) right?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/RightAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_right_approx(const angle<U, T>& a,
                const T allowed_rel_diff = 1E-12,
                const T allowed_abs_diff = 0)
{
    return isclose(a.scalar(), const_angle<U, T>::quarter_turn.scalar(),
                   allowed_rel_diff, allowed_abs_diff);
}

/// is obtuse?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/ObtuseAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_obtuse(const angle<U, T>& a)
{
    return (a > const_angle<U, T>::quarter_turn) &&
           (a < const_angle<U, T>::half_turn);
}

/// is straight?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/StraightAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_straight(const angle<U, T>& a)
{
    return a == const_angle<U, T>::half_turn;
}

/// is (approximately) straight?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/StraightAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_straight_approx(const angle<U, T>& a,
                   const T allowed_rel_diff = 1E-12,
                   const T allowed_abs_diff = 0)
{
    return isclose(a.scalar(), const_angle<U, T>::half_turn.scalar(),
                   allowed_rel_diff, allowed_abs_diff);
}

/// is reflex?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/ReflexAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_reflex(const angle<U, T>& a)
{
    return a > const_angle<U, T>::half_turn;
}

/// is full?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/FullAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_full(const angle<U, T>& a)
{
    return a == const_angle<U, T>::full_turn;
}

/// is (approximately) full?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/FullAngle.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
is_full_approx(const angle<U, T>& a,
               const T allowed_rel_diff = 1E-12,
               const T allowed_abs_diff = 0)
{
    return isclose(a.scalar(), const_angle<U, T>::full_turn.scalar(),
                   allowed_rel_diff, allowed_abs_diff);
}

/// are complementary?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/ComplementaryAngles.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
are_complementary(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return is_right(a1 + a2);
}

/// are (approximately) complementary?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/ComplementaryAngles.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
are_complementary_approx(const angle<U, T>& a1,
                         const angle<U, T>& a2,
                         const T allowed_rel_diff = 1E-12,
                         const T allowed_abs_diff = 0)
{
    return is_right_approx(a1 + a2, allowed_rel_diff, allowed_abs_diff);
}

/// are supplementary?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/SupplementaryAngles.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
are_supplementary(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return is_straight(a1 + a2);
}

/// are (approximately) supplementary?
/**
\pre \a a is non-negative
\sa https://mathworld.wolfram.com/SupplementaryAngles.html
*/
template <angle_unit U, std::floating_point T>
constexpr bool
are_supplementary_approx(const angle<U, T>& a1,
                         const angle<U, T>& a2,
                         const T allowed_rel_diff = 1E-12,
                         const T allowed_abs_diff = 0)
{
    return is_straight_approx(a1 + a2, allowed_rel_diff, allowed_abs_diff);
}

/// get the absolute value of the angle
template <angle_unit U, std::floating_point T>
constexpr auto
abs(const angle<U, T>& a)
{
    return angle<U, T>{std::abs(a.scalar())};
}

/// sin
template <angle_unit U, std::floating_point T>
constexpr T
sin(const angle<U, T>& a)
{
    return std::sin(a.to_rad());
}

/// cos
template <angle_unit U, std::floating_point T>
constexpr T
cos(const angle<U, T>& a)
{
    return std::cos(a.to_rad());
}

/// tan
template <angle_unit U, std::floating_point T>
constexpr T
tan(const angle<U, T>& a)
{
    return std::tan(a.to_rad());
}

/// sin_cos
template <angle_unit U, std::floating_point T>
constexpr void
sin_cos(const angle<U, T>& a, T& s, T& c)
{
    s = std::sin(a.to_rad());
    c = std::cos(a.to_rad());
}

/// inverse sin
template <std::floating_point T>
constexpr auto
a_asin(const T x)
{
    return ang_rad<T>{std::asin(x)};
}

/// inverse cos
template <std::floating_point T>
constexpr auto
a_acos(const T x)
{
    return ang_rad<T>{std::acos(x)};
}

/// inverse tan
template <std::floating_point T>
constexpr auto
a_atan(const T x)
{
    return ang_rad<T>{std::atan(x)};
}

/// inverse tan2
template <std::floating_point T>
constexpr auto
a_atan2(const T y, const T x)
{
    return ang_rad<T>{std::atan2(y, x)};
}

// XXX: hyperbolic functions don't really deal with circular angles

/// sinh
template <angle_unit U, std::floating_point T>
constexpr T
sinh(const angle<U, T>& a)
{
    return std::sinh(a.to_rad());
}

/// cosh
template <angle_unit U, std::floating_point T>
constexpr T
cosh(const angle<U, T>& a)
{
    return std::cosh(a.to_rad());
}

/// tanh
template <angle_unit U, std::floating_point T>
constexpr T
tanh(const angle<U, T>& a)
{
    return std::tanh(a.to_rad());
}

/// inverse sinh
template <std::floating_point T>
constexpr auto
a_asinh(const T x)
{
    return ang_rad<T>{std::asinh(x)};
}

/// inverse cosh
template <std::floating_point T>
constexpr auto
a_acosh(const T x)
{
    return ang_rad<T>{std::acosh(x)};
}

/// inverse tanh
template <std::floating_point T>
constexpr auto
a_atanh(const T x)
{
    return ang_rad<T>{std::atanh(x)};
}

/// get the IEEE remainder of dividing \a a1 by \a a2
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U, std::floating_point T>
constexpr auto
ieee_remainder(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return angle<U, T>{std::remainder(a1.scalar(), a2.scalar())};
}

/// get the IEEE remainder of dividing \a a1 by \a a2
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
ieee_remainder(const angle<U, T>& a1, const angle<U2, T2>& a2)
{
    using result_type = typename std::common_type_t<T, T2>;
    return ieee_remainder(a1.template to<result_type>(),
                          a2.template to<U, result_type>());
}

/// get the fmod remainder of dividing \a a1 by \a a2
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U, std::floating_point T>
constexpr auto
fmod_remainder(const angle<U, T>& a1, const angle<U, T>& a2)
{
    return angle<U, T>{std::fmod(a1.scalar(), a2.scalar())};
}

/// get the fmod remainder of dividing \a a1 by \a a2
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
fmod_remainder(const angle<U, T>& a1, const angle<U2, T2>& a2)
{
    using result_type = typename std::common_type_t<T, T2>;
    return ieee_remainder(a1.template to<result_type>(),
                          a2.template to<U, result_type>());
}

/// normalize the angle
/**
The result will be within the interval [-0.5, 0.5] revolutions.
\param[in,out] a the angle
*/
template <angle_unit U, std::floating_point T>
constexpr void
normalize_angle_signed(angle<U, T>& a)
{
    a = ieee_remainder(a, const_angle<U, T>::full_turn);
}

/// normalize the angle
/**
The result will be within the interval [0, 1] revolutions.
\param[in,out] a the angle
*/
template <angle_unit U, std::floating_point T>
constexpr void
normalize_angle_unsigned(angle<U, T>& a)
{
    normalize_angle_signed(a);

    if (a.scalar() < 0)
        a += const_angle<U, T>::full_turn;
}

/// normalize the geodetic latitude
/**
The result will be within the interval [-90, 90] degrees.
\param[in,out] lat the geodetic latitude
*/
template <angle_unit U, std::floating_point T>
constexpr void
normalize_latitude(angle<U, T>& lat)
{
    normalize_angle_signed(lat);

    if (lat < -const_angle<U, T>::quarter_turn)
        lat = -const_angle<U, T>::half_turn - lat;
    else if (lat > const_angle<U, T>::quarter_turn)
        lat = const_angle<U, T>::half_turn - lat;
}

/// normalize the geodetic longitude
/**
The result will be within the interval [-180, 180] degrees.
\param[in,out] lon the geodetic longitude
*/
template <angle_unit U, std::floating_point T>
constexpr void
normalize_longitude(angle<U, T>& lon)
{
    normalize_angle_signed(lon);
}

/// normalize the geodetic coordinate
/**
\param[in,out] lat the geodetic latitude
\param[in,out] lon the geodetic longitude
*/
template <angle_unit U, std::floating_point T>
constexpr void
normalize_geodetic(angle<U, T>& lat, angle<U, T>& lon)
{
    normalize_angle_signed(lat);

    if (lat < -const_angle<U, T>::quarter_turn)
    {
        lat = -const_angle<U, T>::half_turn - lat;
        lon += const_angle<U, T>::half_turn;
    }
    else if (lat > const_angle<U, T>::quarter_turn)
    {
        lat = const_angle<U, T>::half_turn - lat;
        lon += const_angle<U, T>::half_turn;
    }

    normalize_longitude(lon);
}

/// get the 0-based quadrant that the angle is in
/**
\param a the angle
\return the 0-based quadrant that the angle is in
\retval 0 for quadrant I
\retval 1 for quadrant II
\retval 2 for quadrant III
\retval 3 for quadrant IV
*/
template <angle_unit U, std::floating_point T>
constexpr auto
get_quadrant(angle<U, T> a)
{
    normalize_angle_unsigned(a);
    return static_cast<int>(a / const_angle<U, T>::quarter_turn);
}

/// get the difference between the angles going from the first angle to the second angle
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U, std::floating_point T>
constexpr auto
angle_diff(angle<U, T> a1, angle<U, T> a2)
{
    a2 -= a1;
    normalize_angle_signed(a2);
    return a2;
}

/// get the difference between the angles going from the first angle to the second angle
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U,
          std::floating_point T,
          angle_unit U2,
          std::floating_point T2>
constexpr auto
angle_diff(angle<U, T> a1, angle<U2, T2> a2)
{
    using result_type = typename std::common_type_t<T, T2>;
    return angle_diff(a1.template to<result_type>(),
                      a2.template to<U, result_type>());
}

/// get the remainder and part of the quotient upon division of \a a1 by \a a2
/**
\param a1 the first angle
\param a2 the second angle
*/
template <angle_unit U, std::floating_point T>
constexpr auto
remquo(const angle<U, T>& a1, const angle<U, T>& a2, int& quo)
{
    return angle<U, T>{std::remquo(a1.scalar(), a2.scalar(), &quo)};
}

// XXX: this seems like more trouble than it's worth

template <std::floating_point T>
constexpr bool
atan_boundary_case(const T y)
{
    switch (std::fpclassify(y))
    {
    case FP_INFINITE:
    case FP_NAN:
    case FP_ZERO:
        return true;
        break;

    case FP_NORMAL:
    case FP_SUBNORMAL:
    default:
        return false;
        break;
    }
}

template <std::floating_point T>
constexpr bool
atan2_boundary_case(const T y, const T x)
{
    return atan_boundary_case(y) || atan_boundary_case(x);
}

/// sin_cos with argument reduction
template <angle_unit U, std::floating_point T>
constexpr void
sin_cos_r(angle<U, T> a, T& s, T& c)
{
    int quo = 0;

    a = remquo(a, const_angle<U, T>::quarter_turn, quo);

    T _s{};
    T _c{};
    sin_cos(a, _s, _c);

    switch (quo % 4U) // 360/90 == 4
    {
    default: break; // prevent warning: switch-default
    case 0U: s =  _s; c =  _c; break; // shift by 0*pi/2 (0 deg)
    case 1U: s =  _c; c = -_s; break; // shift by 1*pi/2 (90 deg)
    case 2U: s = -_s; c = -_c; break; // shift by 2*pi/2 (180 deg)
    case 3U: s = -_c; c =  _s; break; // shift by 3*pi/2 (270 deg)
    }
}

/// sin with argument reduction
template <angle_unit U, std::floating_point T>
constexpr T
sin_r(angle<U, T> a)
{
    int quo = 0;

    a = remquo(a, const_angle<U, T>::quarter_turn, quo);

    T s{};

    switch (quo % 4U) // 360/90 == 4
    {
    default: break; // prevent warning: switch-default
#if 1
    case 0U: s =  sin(a); break; // shift by 0*pi/2 (0 deg)
    case 1U: s =  cos(a); break; // shift by 1*pi/2 (90 deg)
    case 2U: s = -sin(a); break; // shift by 2*pi/2 (180 deg)
    case 3U: s = -cos(a); break; // shift by 3*pi/2 (270 deg)
#else
    case 0U: s =  sinp(a.to_rad()); break; // shift by 0*pi/2 (0 deg)
    case 1U: s =  cosp(a.to_rad()); break; // shift by 1*pi/2 (90 deg)
    case 2U: s = -sinp(a.to_rad()); break; // shift by 2*pi/2 (180 deg)
    case 3U: s = -cosp(a.to_rad()); break; // shift by 3*pi/2 (270 deg)
#endif
    }

    return s;
}

/// cos with argument reduction
template <angle_unit U, std::floating_point T>
constexpr T
cos_r(angle<U, T> a)
{
    int quo = 0;

    a = remquo(a, const_angle<U, T>::quarter_turn, quo);

    T c{};

    switch (quo % 4U)
    {
    default: break; // prevent warning: switch-default
#if 1
    case 0U: c =  cos(a); break; // shift by 0*pi/2 (0 deg)
    case 1U: c = -sin(a); break; // shift by 1*pi/2 (90 deg)
    case 2U: c = -cos(a); break; // shift by 2*pi/2 (180 deg)
    case 3U: c =  sin(a); break; // shift by 3*pi/2 (270 deg)
#else
    case 0U: c =  cosp(a.to_rad()); break; // shift by 0*pi/2 (0 deg)
    case 1U: c = -sinp(a.to_rad()); break; // shift by 1*pi/2 (90 deg)
    case 2U: c = -cosp(a.to_rad()); break; // shift by 2*pi/2 (180 deg)
    case 3U: c =  sinp(a.to_rad()); break; // shift by 3*pi/2 (270 deg)
#endif
    }

    return c;
}

/// tan with argument reduction
template <angle_unit U, std::floating_point T>
constexpr T
tan_r(const angle<U, T>& a)
{
    T s{};
    T c{};
    sin_cos_r(a, s, c);

    return s / c;
}

/// atan2 with argument reduction
template <std::floating_point T>
constexpr auto
a_atan2_r(T y, T x)
{
    if (atan2_boundary_case(y, x))
        return a_atan2(y, x);

    int q = 0;
    //auto a = a_atan2_r_helper(y, x, q);

    q = 0;

    if (std::abs(y) > std::abs(x))
    {
        // x, y = y, -x
        x = -x;
        std::swap(x, y);
        q = 1; // shift by 1*pi/2 (90 deg)
    }

    if (x < 0)
    {
        // x, y = -x, -y
        x = -x;
        y = -y;
        q += 2; // shift by 2*pi/2 (180 deg)
    }

    auto a = a_atan2(y, x);

    switch (q)
    {
    default: break; // prevent warning: switch-default
    //case 0: a += 0; break; // shift by 0 deg
    case 1: a += const_angle<a.units(), T>::quarter_turn; break; // shift by 90 deg
    //case 2: a += (y < 0 ? const_angle<U, T>::half_turn : -const_angle<U, T>::half_turn); break; // shift by 180 or -180 deg
    case 2: a += (std::signbit(y) ? const_angle<a.units(), T>::half_turn : -const_angle<a.units(), T>::half_turn); break; // shift by 180 or -180 deg
    case 3: a += -const_angle<a.units(), T>::quarter_turn; break; // shift by -90 deg
    }

    return a;
}

/// atan with argument reduction
template <std::floating_point T>
constexpr auto
a_atan_r(const T y)
{
    if (atan_boundary_case(y))
        return a_atan(y);

    return a_atan2_r(y, T{1});
}
