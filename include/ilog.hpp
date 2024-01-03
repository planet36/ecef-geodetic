// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// functions to get integral part of logarithm
/**
\file
\author Steven Ward
\sa https://en.cppreference.com/w/cpp/numeric/math/ilogb
\sa https://en.cppreference.com/w/cpp/numeric/math/log2
\sa https://en.cppreference.com/w/cpp/numeric/math/log10
*/

#pragma once

#include <cmath>
#include <concepts>
#include <limits>

/// return the base-2 logarithm of \a x as a signed integer
template <std::floating_point T>
int
ilog2(const T x)
{
	if constexpr (std::numeric_limits<T>::radix == 2)
		return std::ilogb(x);
	else
		return std::floor(std::log2(x));
}

/// return the base-10 logarithm of \a x as a signed integer
template <std::floating_point T>
int
ilog10(const T x)
{
	if constexpr (std::numeric_limits<T>::radix == 10)
		return std::ilogb(x);
	else
		return std::floor(std::log10(x));
}
