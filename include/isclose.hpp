// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// function to determine if two floating-point numbers are close to each other
/**
\file
\author Steven Ward
*/

#pragma once

#include <cmath>
#include <concepts>

/// are two floating-point numbers are close to each other?
/**
This is similar to the Python function \c math.isclose.
\sa https://docs.python.org/3/library/math.html#math.isclose
\sa https://peps.python.org/pep-0485/
\sa https://github.com/python/cpython/blob/main/Modules/mathmodule.c#L2997
\pre \a allowed_rel_diff is positive
\pre \a allowed_abs_diff is non-negative
\param allowed_rel_diff maximum relative difference for being considered "close", relative to the magnitude of the input values
\param allowed_abs_diff maximum absolute difference for being considered "close", regardless of the magnitude of the input values
*/
template <std::floating_point T>
bool
isclose(const T a,
        const T b,
        const T allowed_rel_diff = 1E-9,
        const T allowed_abs_diff = 0)
{
	if (a == b)
		return true;

	if (!std::isfinite(a) || !std::isfinite(b))
		return false;

	const auto actual_abs_diff = std::abs(a - b);
	const auto max_abs = std::max(std::abs(a), std::abs(b));
	return actual_abs_diff <= allowed_abs_diff ||
	       actual_abs_diff <= allowed_rel_diff * max_abs;
}
