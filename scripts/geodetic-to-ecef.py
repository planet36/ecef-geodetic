# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''
Read Geodetic coordinates from stdin and convert them to ECEF coordinates.

If 2 fields are given, they are interpreted as
    geodetic latitude (degrees)
    ellipsoid height (meters)

    The output coordinate represents
        W (meters)
        Z (meters)

If 3 fields are given, they are interpreted as
    geodetic latitude (degrees)
    geodetic longitude (degrees)
    ellipsoid height (meters)

    The output coordinate represents
        X (meters)
        Y (meters)
        Z (meters)
'''

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'

import sys

from ellipsoid import WGS84

for (line_num, line) in enumerate(sys.stdin, start=1):
    match line.split():
        case []:
            continue

        case [lat_deg, ht]:
            (w, z) = WGS84.geodetic_2d_to_ecef(float(lat_deg), float(ht))
            print(f"{w} {z}")

        case [lat_deg, lon_deg, ht]:
            (x, y, z) = WGS84.geodetic_to_ecef(float(lat_deg), float(lon_deg), float(ht))
            print(f"{x} {y} {z}")

        case fields:
            sys.exit(f'Error: line {line_num}: expected 2 or 3 fields, got {len(fields)}')
