# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''Convert Polar coordinates (r, theta (degrees)) to Cartesian coordinates (x, y).'''

import sys
from cmath import rect
from math import radians

for line in sys.stdin:
    if not line.strip():
        continue
    (r, theta_deg) = map(float, line.split())
    z = rect(r, radians(theta_deg))
    print(z.real, z.imag)
