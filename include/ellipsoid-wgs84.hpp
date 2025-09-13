// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// WGS 84 Reference Ellipsoid
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

#include "ellipsoid.hpp"

#include <concepts>

template <std::floating_point T>
inline constexpr Ellipsoid<T> WGS84{6'378'137.0L, 298.257223563L};
