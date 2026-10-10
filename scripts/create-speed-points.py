# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: MPL-2.0

# pylint: disable=invalid-name

'''
Create 2D ECEF points to be used in the speed tests.
'''

__author__ = 'Steven Ward'
__license__ = 'MPL-2.0'
__version__ = '2026-10-10'

from decimal import Context, Decimal, getcontext, localcontext

# Significant decimal digits of the IEEE binary256 and binary128 formats
binary256_digits10 = 72
binary128_digits10 = 35

getcontext().prec = binary256_digits10

# The pi, cos, and sin functions are copied from the recipes at
# <https://docs.python.org/3/library/decimal.html#recipes>.


def pi() -> Decimal:
    '''Compute Pi to the current precision.'''
    with localcontext() as ctx:
        ctx.prec += 2 # extra digits for intermediate steps
        three = Decimal(3)
        # pylint: disable-next=redefined-outer-name
        lasts, t, s, n, na, d, da = Decimal(0), three, three, 1, 0, 0, 24
        while s != lasts:
            lasts = s
            n, na = n + na, na + 8
            d, da = d + da, da + 32
            t = (t * n) / d
            s += t
    return +s # unary plus applies the new precision


def cos(x: Decimal) -> Decimal:
    '''Return the cosine of x as measured in radians.'''
    with localcontext() as ctx:
        ctx.prec += 2
        i, lasts, s, fact, num, sign = 0, Decimal(0), Decimal(1), 1, Decimal(1), 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i - 1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
    return +s


def sin(x: Decimal) -> Decimal:
    '''Return the sine of x as measured in radians.'''
    with localcontext() as ctx:
        ctx.prec += 2
        i, lasts, s, fact, num, sign = 1, Decimal(0), x, 1, x, 1
        while s != lasts:
            lasts = s
            i += 2
            fact *= i * (i - 1)
            num *= x * x
            sign *= -1
            s += num / fact * sign
    return +s


# WGS 84
a = Decimal('6378137.0')
b = Decimal('6356752.31424517949756396659963365515679817131108549733884857165128320852')

zero_threshold = Decimal(f'1E-{binary128_digits10}')


def fix_zero(x: Decimal) -> Decimal:
    '''Treat -0.0 and very small numbers (e.g. 1.3E-65) as 0.0'''
    if x.is_zero() or (abs(x) < zero_threshold):
        return Decimal(0)
    return x


output_context = Context(prec=binary128_digits10)


def to_str(x: Decimal) -> str:
    '''Format x to binary128 precision without trailing zeros.'''
    # Note: read_coords_ecef only supports decimal float
    s = f'{x.normalize(output_context):f}'
    # Keep the decimal point on whole numbers, as in "0.0" and "6378137.0".
    return s if '.' in s else s + '.0'


# graph of ellipsoid and evolute
# https://www.desmos.com/calculator/0kv3gs1lzg
# https://www.desmos.com/calculator/vgwsyhnjvm

all_t = (Decimal(0),
         pi()/6,
         pi()/3,
         pi()/2,
         )

all_r = ((a/500, b/500), # inside the evolute
         (a/2, b/2), # inside the ellipsoid
         (a, b), # on the ellipsoid
         (a*2, b*2), # outside the ellipsoid
         )

points = []

points.append((Decimal(0), Decimal(0)))

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
    print(f'{to_str(p[0])} {to_str(p[1])}')
