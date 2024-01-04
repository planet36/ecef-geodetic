# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name

'''Convert Polar coordinates (r, theta (degrees)) to Cartesian coordinates (x, y).'''

from cmath import rect
from math import radians
import sys

for line in sys.stdin:
    (r, theta_deg) = map(float, line.split())
    z = rect(r, radians(theta_deg))
    print(z.real, z.imag)
