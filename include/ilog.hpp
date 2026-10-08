// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: MPL-2.0

/// functions to get integral part of logarithm
/**
* \file
* \author Steven Ward
* \sa https://en.cppreference.com/w/cpp/numeric/math/ilogb
* \sa https://en.cppreference.com/w/cpp/numeric/math/log2
* \sa https://en.cppreference.com/w/cpp/numeric/math/log10
*/

#pragma once

#include <cmath>
#include <concepts>
#include <limits>

/// return the base-2 logarithm of \a x as a signed integer
/**
* \return the floor of the base-2 logarithm of |\a x|, or what \c std::ilogb returns when
* \a x is zero, infinite, or NaN
*/
template <std::floating_point T>
constexpr int
ilog2(const T x) noexcept
{
    if constexpr (std::numeric_limits<T>::radix == 2)
        return std::ilogb(x);
    else
    {
        // Match std::ilogb where the logarithm is not finite, as ilog10 does.
        if (x == 0)
            return FP_ILOGB0;
        if (std::isnan(x))
            return FP_ILOGBNAN;
        if (std::isinf(x))
            return std::numeric_limits<int>::max();
        return std::floor(std::log2(std::abs(x)));
    }
}

/// return the base-10 logarithm of \a x as a signed integer
/**
* \return the floor of the base-10 logarithm of |\a x|, or what \c std::ilogb returns when
* \a x is zero, infinite, or NaN
*/
template <std::floating_point T>
constexpr int
ilog10(const T x) noexcept
{
    if constexpr (std::numeric_limits<T>::radix == 10)
        return std::ilogb(x);
    else
    {
        // Match std::ilogb where the logarithm is not finite.  The logarithm of zero is
        // -infinity, and of a negative value or NaN is NaN, and converting either to int
        // is undefined behavior.
        if (x == 0)
            return FP_ILOGB0;
        if (std::isnan(x))
            return FP_ILOGBNAN;
        if (std::isinf(x))
            return std::numeric_limits<int>::max();
        return std::floor(std::log10(std::abs(x)));
    }
}
