# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name
# pylint: disable=no-member

'''
Create 2D ECEF points to be used in the speed tests.
'''

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'
__version__ = '2024-01-03'

import gmpy2

# https://gmpy2.readthedocs.io/en/latest/mpfr.html

gmpy2.set_context(gmpy2.ieee(256))

log10_of_2 = gmpy2.log10(2)
pi = gmpy2.const_pi()
cos = gmpy2.cos
sin = gmpy2.sin

# WGS 84
a = gmpy2.mpfr('6378137.0')
b = gmpy2.mpfr('6356752.31424517949756396659963365515679817131108549733884857165128320852')


binary256_digits10 = int(gmpy2.ceil(log10_of_2 * gmpy2.ieee(256).precision))
# 72
binary192_digits10 = int(gmpy2.ceil(log10_of_2 * gmpy2.ieee(192).precision))
# 53
binary128_digits10 = int(gmpy2.ceil(log10_of_2 * gmpy2.ieee(128).precision))
# 35
binary64_digits10 = int(gmpy2.ceil(log10_of_2 * gmpy2.ieee(64).precision))
# 16

zero_threshold = gmpy2.mpfr(f'1E-{binary128_digits10}')

# pylint: disable=redefined-outer-name
def fix_zero(x: gmpy2.mpfr) -> gmpy2.mpfr:
    '''Treat -0.0 and very small numbers (e.g. 1.3E-65) as 0.0'''
    if x.is_zero() or (abs(x) < zero_threshold):
        return gmpy2.zero()
    return x

# graph of ellipsoid and evolute
# https://www.desmos.com/calculator/0kv3gs1lzg
# https://www.desmos.com/calculator/vgwsyhnjvm

all_t = (gmpy2.zero(),
         pi/6,
         pi/3,
         pi/2,
         )

all_r = ((a/500, b/500), # inside the evolute
         (a/2, b/2), # inside the ellipsoid
         (a, b), # on the ellipsoid
         (a*2, b*2), # outside the ellipsoid
         )

points = []

points.append((gmpy2.zero(), gmpy2.zero()))

for t in all_t:
    for r in all_r:
        (w, z) = (r[0] * cos(t), r[1] * sin(t))
        w = fix_zero(w)
        z = fix_zero(z)
        points.append((w, z))
        if not z.is_zero():
            points.append((w, -z))

    # https://mathworld.wolfram.com/EllipseEvolute.html
    # points on the evolute
    (w, z) = ((a**2 - b**2) / a * cos(t)**3, (b**2 - a**2) / b * sin(t)**3)
    w = fix_zero(w)
    z = fix_zero(z)
    points.append((w, z))
    if not z.is_zero():
        points.append((w, -z))

for p in points:
    print(f'{p[0]:.{binary128_digits10}NG} {p[1]:.{binary128_digits10}NG}')
    # Note: read_coords_ecef only supports decimal float
    #print(f'{p[0]:NA} {p[1]:NA}')
