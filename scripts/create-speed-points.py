# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''
Create 2D ECEF points to be used in the speed tests.

The points are typical inputs at geodetic latitudes 5 to 85 degrees and heights -1000 to
10,000 meters, all at longitude 0.
'''

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'
__version__ = '2026-10-10'

from ellipsoid import WGS84

# Latitudes 0 and 90 and height 0 are left out on purpose.  Some algorithms take a shortcut
# on the polar axis, on the equatorial plane, or exactly on the ellipsoid, and the speed test
# would credit them with it.
all_lat_deg = range(5, 90, 5)
all_ht = [ht for ht in range(-1000, 11000, 1000) if ht != 0]

for lat_deg in all_lat_deg:
    for ht in all_ht:
        (w, z) = WGS84.geodetic_2d_to_ecef(lat_deg, ht)
        print(f'{w} {z}')
