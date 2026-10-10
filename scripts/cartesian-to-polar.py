# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''Convert Cartesian coordinates (x, y) to Polar coordinates (r, theta (degrees)).'''

import sys
from cmath import polar
from math import degrees

for line in sys.stdin:
    if not line.strip():
        continue
    (x, y) = map(float, line.split())
    (r, theta_rad) = polar(complex(x, y))
    print(r, degrees(theta_rad))
