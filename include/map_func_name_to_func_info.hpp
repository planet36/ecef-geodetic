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
{"bowring_1976_1", bowring_1976_1::func_info},
{"bowring_1976_2", bowring_1976_2::func_info},
{"bowring_1985_1", bowring_1985_1::func_info},
{"bowring_1985_2", bowring_1985_2::func_info},
{"bowring_toms_1995_1", bowring_toms_1995_1::func_info},
{"bowring_toms_1995_2", bowring_toms_1995_2::func_info},
{"fukushima_1999_1", fukushima_1999_1::func_info},
{"fukushima_1999_customht_1", fukushima_1999_customht_1::func_info},
{"fukushima_2006_1", fukushima_2006_1::func_info},
{"fukushima_2006_2", fukushima_2006_2::func_info},
{"geographiclib", geographiclib::func_info},
{"geographiclib_customht", geographiclib_customht::func_info},
{"geotransformCpp", geotransformCpp::func_info},
{"geotransformCpp_customht", geotransformCpp_customht::func_info},
{"gersten_1961", gersten_1961::func_info},
{"halley_1", halley_1::func_info},
{"halley_2", halley_2::func_info},
{"halley_quick_1", halley_quick_1::func_info},
{"halley_quick_2", halley_quick_2::func_info},
{"heikkinen_1982", heikkinen_1982::func_info},
{"heikkinen_1982_customht", heikkinen_1982_customht::func_info},
{"householder_1", householder_1::func_info},
{"householder_2", householder_2::func_info},
{"householder_quick_1", householder_quick_1::func_info},
{"householder_quick_2", householder_quick_2::func_info},
{"jat_spacetime_geodetic", jat_spacetime_geodetic::func_info},
{"jones_2002_1", jones_2002_1::func_info},
{"ligas_2011_I_1", ligas_2011_I_1::func_info},
{"ligas_2011_I_2", ligas_2011_I_2::func_info},
{"lin_wang_1995_1", lin_wang_1995_1::func_info},
{"lin_wang_1995_2", lin_wang_1995_2::func_info},
{"lin_wang_1995_customht_1", lin_wang_1995_customht_1::func_info},
{"lin_wang_1995_customht_2", lin_wang_1995_customht_2::func_info},
{"long_1974", long_1974::func_info},
{"naive_I_1", naive_I_1::func_info},
{"naive_I_2", naive_I_2::func_info},
{"naive_II_1", naive_II_1::func_info},
{"naive_II_2", naive_II_2::func_info},
{"newton_raphson_1", newton_raphson_1::func_info},
{"newton_raphson_2", newton_raphson_2::func_info},
{"newton_raphson_quick_1", newton_raphson_quick_1::func_info},
{"newton_raphson_quick_2", newton_raphson_quick_2::func_info},
{"olson_1996", olson_1996::func_info},
{"olson_1996_customht", olson_1996_customht::func_info},
{"openglobe", openglobe::func_info},
{"ozone_1985", ozone_1985::func_info},
{"paul_1973", paul_1973::func_info},
{"pollard_2002_ht_1", pollard_2002_ht_1::func_info},
{"pollard_2002_naive_1", pollard_2002_naive_1::func_info},
{"pollard_2002_naive_2", pollard_2002_naive_2::func_info},
{"pollard_2002_newton_1", pollard_2002_newton_1::func_info},
{"pollard_2002_newton_2", pollard_2002_newton_2::func_info},
{"schroder_1", schroder_1::func_info},
{"schroder_2", schroder_2::func_info},
{"schroder_quick_1", schroder_quick_1::func_info},
{"schroder_quick_2", schroder_quick_2::func_info},
{"sedris", sedris::func_info},
{"sedris_customht", sedris_customht::func_info},
{"shu_2010_1", shu_2010_1::func_info},
{"shu_2010_2", shu_2010_2::func_info},
{"shu_2010_customht_1", shu_2010_customht_1::func_info},
{"shu_2010_customht_2", shu_2010_customht_2::func_info},
{"sofair_1993", sofair_1993::func_info},
{"sofair_2000", sofair_2000::func_info},
{"sudano_1997", sudano_1997::func_info},
{"turner_2013", turner_2013::func_info},
{"vermeille_2004", vermeille_2004::func_info},
{"vermeille_2004_customht", vermeille_2004_customht::func_info},
{"vermeille_2011", vermeille_2011::func_info},
{"vermeille_2011_customht", vermeille_2011_customht::func_info},
{"wu_2003_1", wu_2003_1::func_info},
{"wu_2003_2", wu_2003_2::func_info},
{"zhang_2005", zhang_2005::func_info},
{"zhu_1993", zhu_1993::func_info},

};
