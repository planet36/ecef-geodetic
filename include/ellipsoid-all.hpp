// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Reference Ellipsoid Names and Constants
/**
\file
\author Steven Ward
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08

\sa https://nsgreg.nga.mil/doc/view?i=4085
\sa https://web.archive.org/web/20181220230431/https://earth-info.nga.mil/GandG/publications/NGA_STND_0036_1_0_0_WGS84/NGA.STND.0036_1.0.0_WGS84.pdf

\verbatim
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

#include <array>
#include <concepts>
#include <string_view>

/// Reference Ellipsoid ID Codes enum
/**
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08

Page C-2
\verbatim
To facilitate correct referencing, a standardized two letter code is also
included to identify the different ellipsoids and/or their “versions”
pertaining to the different values of the semi-major axis (a).
\endverbatim

Page C-3
\verbatim
Appendix C.1
Reference Ellipsoid Names and Constants
Used for Datum Transformations*

* Refer to Appendices D, E, and F.
** As accepted by NGA.
*** Through adoption of a new yard to meter conversion factor in the referenced
country.
\endverbatim
*/
enum struct ELLIPSOID_ID_CODE : unsigned short
{
    AA,
    AN,
    BR,
    BN,
    CC,
    CD, // **
    CG,
    EB,
    EA,
    EC, // ***
    EF, // ***
    EE,
    ED, // ***
    RF,
    HE,
    HO,
    ID,
    IN,
    KA,
    AM,
    FA,
    SA,
    WO,
    WD,
    WE,
};

/// Convert the Reference Ellipsoid ID Code to a string
constexpr std::string_view
ellipsoid_id_code_to_str(const ELLIPSOID_ID_CODE ell_id_code)
{
    switch (ell_id_code)
    {
    case ELLIPSOID_ID_CODE::AA: return "AA"; break;
    case ELLIPSOID_ID_CODE::AN: return "AN"; break;
    case ELLIPSOID_ID_CODE::BR: return "BR"; break;
    case ELLIPSOID_ID_CODE::BN: return "BN"; break;
    case ELLIPSOID_ID_CODE::CC: return "CC"; break;
    case ELLIPSOID_ID_CODE::CD: return "CD"; break;
    case ELLIPSOID_ID_CODE::CG: return "CG"; break;
    case ELLIPSOID_ID_CODE::EB: return "EB"; break;
    case ELLIPSOID_ID_CODE::EA: return "EA"; break;
    case ELLIPSOID_ID_CODE::EC: return "EC"; break;
    case ELLIPSOID_ID_CODE::EF: return "EF"; break;
    case ELLIPSOID_ID_CODE::EE: return "EE"; break;
    case ELLIPSOID_ID_CODE::ED: return "ED"; break;
    case ELLIPSOID_ID_CODE::RF: return "RF"; break;
    case ELLIPSOID_ID_CODE::HE: return "HE"; break;
    case ELLIPSOID_ID_CODE::HO: return "HO"; break;
    case ELLIPSOID_ID_CODE::ID: return "ID"; break;
    case ELLIPSOID_ID_CODE::IN: return "IN"; break;
    case ELLIPSOID_ID_CODE::KA: return "KA"; break;
    case ELLIPSOID_ID_CODE::AM: return "AM"; break;
    case ELLIPSOID_ID_CODE::FA: return "FA"; break;
    case ELLIPSOID_ID_CODE::SA: return "SA"; break;
    case ELLIPSOID_ID_CODE::WO: return "WO"; break;
    case ELLIPSOID_ID_CODE::WD: return "WD"; break;
    case ELLIPSOID_ID_CODE::WE: return "WE"; break;
    default: return "???"; break;
    }
}

/// Reference Ellipsoid Constants array
/**
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08

Page C-3
\verbatim
Appendix C.1
Reference Ellipsoid Names and Constants
Used for Datum Transformations*

* Refer to Appendices D, E, and F.
** As accepted by NGA.
*** Through adoption of a new yard to meter conversion factor in the referenced
country.
\endverbatim
*/
template <std::floating_point T>
constexpr std::array<Ellipsoid<T>, 25> reference_ellipsoids{
    /*AA*/ Ellipsoid<T>{6'377'563.396L, 299.3249646L  },
    /*AN*/ Ellipsoid<T>{6'378'160.0L  , 298.25L       },
    /*BR*/ Ellipsoid<T>{6'377'397.155L, 299.1528128L  },
    /*BN*/ Ellipsoid<T>{6'377'483.865L, 299.1528128L  },
    /*CC*/ Ellipsoid<T>{6'378'206.4L  , 294.9786982L  },
    /*CD*/ Ellipsoid<T>{6'378'249.145L, 293.465L      }, // **
    /*CG*/ Ellipsoid<T>{6'378'249.2L  , 293.4660208L  },
    /*EB*/ Ellipsoid<T>{6'377'298.556L, 300.8017L     },
    /*EA*/ Ellipsoid<T>{6'377'276.345L, 300.8017L     },
    /*EC*/ Ellipsoid<T>{6'377'301.243L, 300.8017L     }, // ***
    /*EF*/ Ellipsoid<T>{6'377'309.613L, 300.8017L     }, // ***
    /*EE*/ Ellipsoid<T>{6'377'304.063L, 300.8017L     },
    /*ED*/ Ellipsoid<T>{6'377'295.664L, 300.8017L     }, // ***
    /*RF*/ Ellipsoid<T>{6'378'137.0L  , 298.257222101L},
    /*HE*/ Ellipsoid<T>{6'378'200.0L  , 298.3L        },
    /*HO*/ Ellipsoid<T>{6'378'270.0L  , 297.0L        },
    /*ID*/ Ellipsoid<T>{6'378'160.0L  , 298.247L      },
    /*IN*/ Ellipsoid<T>{6'378'388.0L  , 297.0L        },
    /*KA*/ Ellipsoid<T>{6'378'245.0L  , 298.3L        },
    /*AM*/ Ellipsoid<T>{6'377'340.189L, 299.3249646L  },
    /*FA*/ Ellipsoid<T>{6'378'155.0L  , 298.3L        },
    /*SA*/ Ellipsoid<T>{6'378'160.0L  , 298.25L       },
    /*WO*/ Ellipsoid<T>{6'378'300.58L , 296.0L        },
    /*WD*/ Ellipsoid<T>{6'378'135.0L  , 298.26L       },
    /*WE*/ Ellipsoid<T>{6'378'137.0L  , 298.257223563L},
};

/// Reference Ellipsoid Names array
/**
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08

Page C-3
\verbatim
Appendix C.1
Reference Ellipsoid Names and Constants
Used for Datum Transformations*

* Refer to Appendices D, E, and F.
** As accepted by NGA.
*** Through adoption of a new yard to meter conversion factor in the referenced
country.
\endverbatim
*/
constexpr std::array<std::string_view, 25> reference_ellipsoid_names{
    /*AA*/ "Airy 1830",
    /*AN*/ "Australian National 1966",
    /*BR*/ "Bessel 1841: Ethiopia, Indonesia, Japan and Korea",
    /*BN*/ "Bessel 1841: Namibia",
    /*CC*/ "Clarke 1866",
    /*CD*/ "Clarke 1880", // **
    /*CG*/ "Clarke 1880 (IGN)",
    /*EB*/ "Everest: Brunei and E. Malaysia (Sabah and Sarawak)",
    /*EA*/ "Everest: India 1830",
    /*EC*/ "Everest: India 1956", // ***
    /*EF*/ "Everest: Pakistan", // ***
    /*EE*/ "Everest: W. Malaysia and Singapore 1948",
    /*ED*/ "Everest: W. Malaysia 1969", // ***
    /*RF*/ "Geodetic Reference System 1980",
    /*HE*/ "Helmert 1906",
    /*HO*/ "Hough 1960",
    /*ID*/ "Indonesian 1974",
    /*IN*/ "International 1924",
    /*KA*/ "Krassovsky 1940",
    /*AM*/ "Modified Airy",
    /*FA*/ "Modified Fischer 1960",
    /*SA*/ "South American 1969",
    /*WO*/ "War Office 1924",
    /*WD*/ "WGS 1972",
    /*WE*/ "WGS 1984",
};
