# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''
Create 2D geodetic points near the equatorial plane.

The latitudes are plus and minus 10**-k degrees for k = 5 to 13.  At the surface they put the
points from about 1.1 m down to 11 nm from the plane.  The heights are -10 km to 100 km.
'''

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'
__version__ = '2026-10-11'

from decimal import Decimal

# The latitudes are exact in decimal, so the text holds every digit.  Latitude 0 is left
# out, because the region files already put points on the plane.
all_lat_deg = sorted(sign * Decimal(10) ** -k for k in range(5, 14) for sign in (-1, 1))

# Height 0 is left out, because some algorithms take a shortcut exactly on the ellipsoid.
all_ht = sorted((set(range(-100, 1100, 100)) | set(range(-10_000, 101_000, 1_000))) - {0})

for lat_deg in all_lat_deg:
    for ht in all_ht:
        print(f'{lat_deg:f} {ht}')
