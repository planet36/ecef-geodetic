// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// A map of function names to function info objects
/**
\file
\author Steven Ward
*/

#pragma once

#include "ecef_to_geodetic-funcs.hpp"

#include <map>
#include <string>

/*
# command to generate this map:
grep namespace ecef_to_geodetic-funcs.hpp | awk '{print "{\"" $2 "\", " $2 "::func_info},"}'
*/

const std::map<std::string, func_info_t> map_func_name_to_func_info {

{"borkowski_1989", borkowski_1989::func_info},
{"bowring_1976_x1", bowring_1976_x1::func_info},
{"bowring_1976_x2", bowring_1976_x2::func_info},
{"bowring_1985_x1", bowring_1985_x1::func_info},
{"bowring_1985_x2", bowring_1985_x2::func_info},
{"bowring_toms_1995_x1", bowring_toms_1995_x1::func_info},
{"bowring_toms_1995_x2", bowring_toms_1995_x2::func_info},
{"fukushima_1999_x1", fukushima_1999_x1::func_info},
{"fukushima_1999_customht_x1", fukushima_1999_customht_x1::func_info},
{"fukushima_2006_x1", fukushima_2006_x1::func_info},
{"fukushima_2006_x2", fukushima_2006_x2::func_info},
{"geographiclib", geographiclib::func_info},
{"geographiclib_customht", geographiclib_customht::func_info},
{"geotransformCpp", geotransformCpp::func_info},
{"geotransformCpp_customht", geotransformCpp_customht::func_info},
{"gersten_1961", gersten_1961::func_info},
{"halley_x1", halley_x1::func_info},
{"halley_x2", halley_x2::func_info},
{"halley_quick_x1", halley_quick_x1::func_info},
{"halley_quick_x2", halley_quick_x2::func_info},
{"heikkinen_1982", heikkinen_1982::func_info},
{"heikkinen_1982_customht", heikkinen_1982_customht::func_info},
{"householder_x1", householder_x1::func_info},
{"householder_x2", householder_x2::func_info},
{"householder_quick_x1", householder_quick_x1::func_info},
{"householder_quick_x2", householder_quick_x2::func_info},
{"jat_geodetic", jat_geodetic::func_info},
{"jones_2002_x1", jones_2002_x1::func_info},
{"ligas_2011_I_x1", ligas_2011_I_x1::func_info},
{"ligas_2011_I_x2", ligas_2011_I_x2::func_info},
{"lin_wang_1995_x1", lin_wang_1995_x1::func_info},
{"lin_wang_1995_x2", lin_wang_1995_x2::func_info},
{"lin_wang_1995_customht_x1", lin_wang_1995_customht_x1::func_info},
{"lin_wang_1995_customht_x2", lin_wang_1995_customht_x2::func_info},
{"long_1974", long_1974::func_info},
{"naive_I_x1", naive_I_x1::func_info},
{"naive_I_x2", naive_I_x2::func_info},
{"naive_II_x1", naive_II_x1::func_info},
{"naive_II_x2", naive_II_x2::func_info},
{"newton_raphson_x1", newton_raphson_x1::func_info},
{"newton_raphson_x2", newton_raphson_x2::func_info},
{"newton_raphson_quick_x1", newton_raphson_quick_x1::func_info},
{"newton_raphson_quick_x2", newton_raphson_quick_x2::func_info},
{"olson_1996", olson_1996::func_info},
{"olson_1996_customht", olson_1996_customht::func_info},
{"openglobe", openglobe::func_info},
{"ozone_1985", ozone_1985::func_info},
{"paul_1973", paul_1973::func_info},
{"pollard_2002_ht_x1", pollard_2002_ht_x1::func_info},
{"pollard_2002_naive_x1", pollard_2002_naive_x1::func_info},
{"pollard_2002_naive_x2", pollard_2002_naive_x2::func_info},
{"pollard_2002_newton_x1", pollard_2002_newton_x1::func_info},
{"pollard_2002_newton_x2", pollard_2002_newton_x2::func_info},
{"schroder_x1", schroder_x1::func_info},
{"schroder_x2", schroder_x2::func_info},
{"schroder_quick_x1", schroder_quick_x1::func_info},
{"schroder_quick_x2", schroder_quick_x2::func_info},
{"sedris", sedris::func_info},
{"sedris_customht", sedris_customht::func_info},
{"shu_2010_x1", shu_2010_x1::func_info},
{"shu_2010_x2", shu_2010_x2::func_info},
{"shu_2010_customht_x1", shu_2010_customht_x1::func_info},
{"shu_2010_customht_x2", shu_2010_customht_x2::func_info},
{"sofair_1993", sofair_1993::func_info},
{"sofair_2000", sofair_2000::func_info},
{"sudano_1997", sudano_1997::func_info},
{"turner_2013", turner_2013::func_info},
{"vermeille_2004", vermeille_2004::func_info},
{"vermeille_2004_customht", vermeille_2004_customht::func_info},
{"vermeille_2011", vermeille_2011::func_info},
{"vermeille_2011_customht", vermeille_2011_customht::func_info},
{"wu_2003_x1", wu_2003_x1::func_info},
{"wu_2003_x2", wu_2003_x2::func_info},
{"zhang_2005", zhang_2005::func_info},
{"zhu_1993", zhu_1993::func_info},

};
