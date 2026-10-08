// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// angle unit enum
/**
* \file
* \author Steven Ward
*/

#pragma once

#include "angle-conv-utils.hpp"

#include <concepts>
#include <string_view>

/// the fundamental unit of measurement of an angle
/**
* \sa https://en.wikipedia.org/wiki/Angle#Units
* \sa https://en.wikipedia.org/wiki/Angular_unit
*/
enum struct angle_unit : unsigned char
{
    milliradian,
    radian     ,
    revolution ,
    degree     ,
    arcminute  ,
    arcsecond  ,
};

/// convert the \c angle_unit to a string
template <angle_unit U>
[[nodiscard]] constexpr std::string_view
angle_unit_to_string() noexcept
{
    if      constexpr (U == angle_unit::milliradian) {return "mrad"  ;}
    else if constexpr (U == angle_unit::radian     ) {return "rad"   ;}
    else if constexpr (U == angle_unit::revolution ) {return "rev"   ;}
    else if constexpr (U == angle_unit::degree     ) {return "deg"   ;}
    else if constexpr (U == angle_unit::arcminute  ) {return "arcmin";}
    else if constexpr (U == angle_unit::arcsecond  ) {return "arcsec";}
    else static_assert(false, "unhandled angle_unit");
}

/// in the forward order of angle units, get the next \c angle_unit from \a U
template <angle_unit U>
constexpr angle_unit
angle_unit_forward() noexcept
{
    if      constexpr (U == angle_unit::milliradian) {return angle_unit::radian    ;}
    else if constexpr (U == angle_unit::radian     ) {return angle_unit::revolution;}
    else if constexpr (U == angle_unit::revolution ) {return angle_unit::degree    ;}
    else if constexpr (U == angle_unit::degree     ) {return angle_unit::arcminute ;}
    else if constexpr (U == angle_unit::arcminute  ) {return angle_unit::arcsecond ;}
    else if constexpr (U == angle_unit::arcsecond  ) {return angle_unit::arcsecond ;} // same
    else static_assert(false, "unhandled angle_unit");
}

/// in the reverse order of angle units, get the next \c angle_unit from \a U
template <angle_unit U>
constexpr angle_unit
angle_unit_backward() noexcept
{
    if      constexpr (U == angle_unit::milliradian) {return angle_unit::milliradian;} // same
    else if constexpr (U == angle_unit::radian     ) {return angle_unit::milliradian;}
    else if constexpr (U == angle_unit::revolution ) {return angle_unit::radian     ;}
    else if constexpr (U == angle_unit::degree     ) {return angle_unit::revolution ;}
    else if constexpr (U == angle_unit::arcminute  ) {return angle_unit::degree     ;}
    else if constexpr (U == angle_unit::arcsecond  ) {return angle_unit::arcminute  ;}
    else static_assert(false, "unhandled angle_unit");
}

/// get the next angle_unit from \a U to \a U2
template <angle_unit U, angle_unit U2>
constexpr angle_unit
angle_unit_toward() noexcept
{
    if constexpr (U == U2)
        return U;

    if constexpr (U < U2)
        return angle_unit_forward<U>();
    else
        return angle_unit_backward<U>();
}

/// convert \a x to \a To units from \a From units
/**
* \tparam To the angle unit to convert to
* \tparam From the angle unit to convert from
* \param x the given angle value
* \return \a x converted to \a To units from \a From units
*/
template <angle_unit To, angle_unit From, std::floating_point T>
constexpr T
convert_to(const T& x) noexcept
{
    // convert same angle unit (i.e. do nothing)

    if constexpr (To == From) { return x; }

    // convert between adjacent units

    if constexpr (To == angle_unit::radian && From == angle_unit::milliradian) {return rad_from_mrad(x);}
    if constexpr (To == angle_unit::milliradian && From == angle_unit::radian) {return mrad_from_rad(x);}

    if constexpr (To == angle_unit::revolution && From == angle_unit::radian) {return rev_from_rad(x);}
    if constexpr (To == angle_unit::radian && From == angle_unit::revolution) {return rad_from_rev(x);}

    if constexpr (To == angle_unit::degree && From == angle_unit::revolution) {return deg_from_rev(x);}
    if constexpr (To == angle_unit::revolution && From == angle_unit::degree) {return rev_from_deg(x);}

    if constexpr (To == angle_unit::arcminute && From == angle_unit::degree) {return arcmin_from_deg(x);}
    if constexpr (To == angle_unit::degree && From == angle_unit::arcminute) {return deg_from_arcmin(x);}

    if constexpr (To == angle_unit::arcsecond && From == angle_unit::arcminute) {return arcsec_from_arcmin(x);}
    if constexpr (To == angle_unit::arcminute && From == angle_unit::arcsecond) {return arcmin_from_arcsec(x);}

    // TODO: it might be possible to do conversion of the irregular units here

    // convert between non-adjacent units

    constexpr angle_unit inter = angle_unit_toward<To, From>();

    return convert_to<To, inter>(convert_to<inter, From>(x));
}

/// convert \a x from \a From units to \a To units
/**
* \tparam From the angle unit to convert from
* \tparam To the angle unit to convert to
* \param x the given angle value
* \return \a x converted from \a From units to \a To units
*/
template <angle_unit From, angle_unit To, std::floating_point T>
constexpr T
convert_from(const T& x) noexcept
{
    return convert_to<To, From>(x);
}
