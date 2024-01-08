# SPDX-FileCopyrightText: Steven Ward
# SPDX-License-Identifier: OSL-3.0

# pylint: disable=invalid-name
# pylint: disable=line-too-long
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring

__author__ = 'Steven Ward'
__license__ = 'OSL-3.0'

import math
import numpy as np

# pylint: disable=too-many-instance-attributes
class Ellipsoid:

    # pylint: disable=too-many-locals
    def __init__(self,
                 _a: float,
                 _f_recip: float,
                 _GM: float = 3.986004418E14,
                 _omega: float = 7.292115E-5):

        # defining parameters
        a = _a # semi-major axis (equatorial radius of the earth) (meters)
        f = 1 / _f_recip # (a-b)/a # flattening factor of the earth
        GM = _GM # geocentric gravitational constant (m³/s²)
        omega = _omega # nominal mean angular velocity of the earth (rad/s)

        # derived geometric constants
        b = a*(1-f) # semi-minor axis (meters)
        a2 = a*a # a squared
        b2 = b*b # b squared
        fp = f/(1-f) # (a-b)/b # second flattening
        n = f/(2-f) # (a-b)/(a+b) # third flattening
        e2 = f*(2-f) # (a2-b2)/a2 # first eccentricity squared
        e = math.sqrt(e2) # first eccentricity
        ep2 = e2/(1-e2) # (a2-b2)/b2 # second eccentricity squared
        ep = math.sqrt(ep2) # second eccentricity
        epp2 = e2/(2-e2) # (a2-b2)/(a2+b2) # third eccentricity squared
        epp = math.sqrt(epp2) # third eccentricity
        c2 = a2 - b2 # linear eccentricity squared
        c = math.sqrt(c2) # linear eccentricity
        alpha = math.asin(e) # angular eccentricity # acos(b/a)

        # derived physical constants
        gamma_e = 9.7803253359 # normal gravity at the equator (on the ellipsoid) (m/s²)
        gamma_p = 9.8321849379 # normal gravity at the poles (on the ellipsoid) (m/s²)
        k = (1 - f) * gamma_p / gamma_e - 1 # Somigliana's Formula - normal gravity formula constant
        m = omega * omega * a2 * b / GM # normal gravity formula constant

        self.a       = a
        self.f       = f
        self.GM      = GM
        self.omega   = omega
        self.b       = b
        self.a2      = a2
        self.b2      = b2
        self.fp      = fp
        self.n       = n
        self.e2      = e2
        self.e       = e
        self.ep2     = ep2
        self.ep      = ep
        self.epp2    = epp2
        self.epp     = epp
        self.c2      = c2
        self.c       = c
        self.alpha   = alpha
        self.gamma_e = gamma_e
        self.gamma_p = gamma_p
        self.k       = k
        self.m       = m

    def get_Rn(self, sin_lat: np.array) -> np.array:
        d2 = 1 - self.e2 * sin_lat * sin_lat
        d = np.sqrt(d2)
        return self.a / d

    def get_R(self, sin_lat: np.array) -> np.array:
        return self.get_Rn(sin_lat) * np.sqrt(1 - self.e2 * sin_lat * sin_lat * (2 - self.e2))

    def get_Rm(self, sin_lat: np.array) -> np.array:
        d2 = 1 - self.e2 * sin_lat * sin_lat
        d = np.sqrt(d2)
        return self.a * (1 - self.e2) / (d2 * d)

    def get_gamma(self, sin_lat: np.array) -> np.array:
        d2 = 1 - self.e2 * sin_lat * sin_lat
        d = np.sqrt(d2)
        return self.gamma_e * (1 + self.k * sin_lat * sin_lat) / d

    def get_gamma_h(self, sin_lat: np.array, ht: np.array) -> np.array:
        return self.get_gamma(sin_lat) * (1
                - 2 * ht * (1 + self.f + self.m - 2 * self.f * sin_lat * sin_lat) / self.a
                + 3 * ht * ht / self.a2)

    # pylint: disable=too-many-arguments
    def get_ht(self, w: np.array, z: np.array, sin_lat: np.array, cos_lat: np.array, Rn: np.array) -> np.array:
        # pylint: disable=no-else-return
        # https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
        # cos(45°) == 1/√(2)
        if cos_lat > 1 / np.sqrt(2): # Equatorial
            return w / cos_lat - Rn
        else: # Polar
            return z / sin_lat - Rn * (1 - self.e2)

    def geodetic_to_ecef(self, lat_deg: np.array, lon_deg: np.array, ht: np.array) -> np.array:

        lat_rad = np.radians(lat_deg)
        lon_rad = np.radians(lon_deg)

        sin_lat = np.sin(lat_rad)
        cos_lat = np.cos(lat_rad)

        sin_lon = np.sin(lon_rad)
        cos_lon = np.cos(lon_rad)

        Rn = self.get_Rn(sin_lat)

        x = (Rn + ht) * cos_lat * cos_lon
        y = (Rn + ht) * cos_lat * sin_lon
        z = (Rn * (1 - self.e2) + ht) * sin_lat

        return np.stack((x, y, z), axis=1)

    def geodetic_2d_to_ecef(self, lat_deg: np.array, ht: np.array) -> np.array:

        lat_rad = np.radians(lat_deg)

        sin_lat = np.sin(lat_rad)
        cos_lat = np.cos(lat_rad)

        Rn = self.get_Rn(sin_lat)

        w = (Rn + ht) * cos_lat
        z = (Rn * (1 - self.e2) + ht) * sin_lat

        return np.stack((w, z), axis=1)

    def ecef_to_geodetic(self, x: np.array, y: np.array, z: np.array) -> np.array:
        '''Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476. https://doi.org/10.1109/7.481290

Converted to Python and modified by Steven Ward.  No rights reserved.
'''

        w2 = x * x + y * y
        w = np.sqrt(w2)
        z2 = z * z
        lon_rad = np.atan2(y, x)

        a1 = self.a * self.e2
        a2 = a1 * a1
        a3 = a1 * self.e2 / 2
        a4 = 2.5 * a2
        a5 = a1 + a3
        #a6 = (1 - self.e2)

        r2 = w2 + z2
        r = np.sqrt(r2)

        s2 = z2 / r2
        c2 = w2 / r2
        u = a2 / r
        v = a3 - a4 / r

        s = 0
        c = 0
        ss = 0

        # cos(45°)² == ½
        if c2 > 0.5: # Equatorial
            s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r)
            lat_rad = np.asin(s)
            ss = s * s
            c = np.sqrt(1 - ss)
        else: # Polar
            c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r)
            lat_rad = np.acos(c)
            ss = 1 - c * c
            s = np.sqrt(ss)

            if z < 0:
                lat_rad = -lat_rad
                s = -s

        d2 = 1 - self.e2 * ss
        Rn = self.a / np.sqrt(d2)
        Rm = (1 - self.e2) * Rn / d2
        rf = (1 - self.e2) * Rn
        u = w - Rn * c
        v = z - rf * s
        f = c * u + s * v
        m = c * v - s * u
        p = m / (Rm + f)

        lat_rad += p

        ht = f + m * p / 2

        return np.stack((np.degrees(lat_rad), np.degrees(lon_rad), ht), axis=1)

    def ecef_2d_to_geodetic(self, w: np.array, z: np.array) -> np.array:
        '''Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476. https://doi.org/10.1109/7.481290

Converted to Python and modified by Steven Ward.  No rights reserved.
'''

        w2 = w * w
        z2 = z * z

        a1 = self.a * self.e2
        a2 = a1 * a1
        a3 = a1 * self.e2 / 2
        a4 = 2.5 * a2
        a5 = a1 + a3
        #a6 = (1 - self.e2)

        r2 = w2 + z2
        r = np.sqrt(r2)

        s2 = z2 / r2
        c2 = w2 / r2
        u = a2 / r
        v = a3 - a4 / r

        s = 0
        c = 0
        ss = 0

        # cos(45°)² == ½
        if c2 > 0.5: # Equatorial
            s = (z / r) * (1 + c2 * (a1 + u + s2 * v) / r)
            lat_rad = np.asin(s)
            ss = s * s
            c = np.sqrt(1 - ss)
        else: # Polar
            c = (w / r) * (1 - s2 * (a5 - u - c2 * v) / r)
            lat_rad = np.acos(c)
            ss = 1 - c * c
            s = np.sqrt(ss)

            if z < 0:
                lat_rad = -lat_rad
                s = -s

        d2 = 1 - self.e2 * ss
        Rn = self.a / np.sqrt(d2)
        Rm = (1 - self.e2) * Rn / d2
        rf = (1 - self.e2) * Rn
        u = w - Rn * c
        v = z - rf * s
        f = c * u + s * v
        m = c * v - s * u
        p = m / (Rm + f)

        lat_rad += p

        ht = f + m * p / 2

        return np.stack((np.degrees(lat_rad), ht), axis=1)

WGS84 = Ellipsoid(6_378_137.0, 298.257223563)
'''
Source:
NGA.STND.0036_1.0.0_WGS84 2014-07-08

Appendix C.1
Reference Ellipsoid Names and Constants
Used for Datum Transformations*

* Refer to Appendices D, E, and F.
'''
