# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name

'''Convert Cartesian coordinates (x, y) to Polar coordinates (r, theta (degrees)).'''

from cmath import polar
from math import degrees
import sys

for line in sys.stdin:
    (x, y) = map(float, line.split())
    (r, theta_rad) = polar(complex(x, y))
    print(r, degrees(theta_rad))
