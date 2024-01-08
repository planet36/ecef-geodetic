# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name

'''
Read ECEF coordinates from stdin and convert them to Geodetic coordinates.

If 2 fields are given, they are interpreted as
    W (meters)
    Z (meters)

    The output coordinate represents
        geodetic latitude (degrees)
        ellipsoid height (meters)

If 3 fields are given, they are interpreted as
    X (meters)
    Y (meters)
    Z (meters)

    The output coordinate represents
        geodetic latitude (degrees)
        geodetic longitude (degrees)
        ellipsoid height (meters)
'''

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

import sys

from ellipsoid import WGS84

for line in sys.stdin:
    fields = line.strip().split()

    if len(fields) == 2:
        w = float(fields[0])
        z = float(fields[1])
        (lat_deg, ht) = WGS84.ecef_2d_to_geodetic(w, z)
        print(f"{lat_deg} {ht}")

    elif len(fields) == 3:
        x = float(fields[0])
        y = float(fields[1])
        z = float(fields[2])
        (lat_deg, lon_deg, ht) = WGS84.ecef_to_geodetic(x, y, z)
        print(f"{lat_deg} {lon_deg} {ht}")
