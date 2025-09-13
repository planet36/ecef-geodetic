// SPDX-FileCopyrightText: Donald K. Olson
// SPDX-License-Identifier: CC-PDDC
// SPDX-License-Identifier: CC0-1.0

#include <math.h>

/// convert from ECEF to geodetic
/**
* Olson, D. K. (1996). Converting Earth-Centered, Earth-Fixed Coordinates to Geodetic Coordinates. IEEE Transactions on Aerospace and Electronic Systems, 32(1), 473–476. https://doi.org/10.1109/7.481290
*
* U.S. Government work, U.S. copyright does not apply.
*
* \sa https://ieeexplore.ieee.org/document/481290
* \param x X coordinate (meters)
* \param y Y coordinate (meters)
* \param z Z coordinate (meters)
* \param[out] lat geodetic latitude (radians)
* \param[out] lon geodetic longitude (radians)
* \param[out] ht ellipsoid height (meters)
*/
void
latlon (x,y,z,lat,lon,ht)/*c code*/
double x,y,z,*lat,*lon,*ht;
{
    double a = 6378137.0; /*wgs-84*/
    double e2 = 6.6943799901377997e-3;
    double a1 = 4.2697672707157535e+4;
    double a2 = 1.8230912546075455e+9;
    double a3 = 1.4291722289812413e+2;
    double a4 = 4.5577281365188637e+9;
    double a5 = 4.2840589930055659e+4;
    double a6 = 9.9330562000986220e-1;
    /* a1 = a*e2, a2 = a1*a1, a3 = a1*e2/2, */
    /* a4 = (5/2)*a2, a5 = a1 + a3, a6 = 1 - e2 */
    double zp,w2,w,z2,r2,r,s2,c2,s,c,ss;
    double g,rg,rf,u,v,m,f,p;
    zp = fabs(z);
    w2 = x*x + y*y;
    w = sqrt(w2);
    z2 = z*z;
    r2 = w2 + z2;
    r = sqrt(r2);
    if (r < 100000.)
    {
        *lat = 0.;
        *lon = 0.;
        *ht = -1.e7;
        return;
    }
    *lon = atan2(y,x);
    s2 = z2/r2;
    c2 = w2/r2;
    u = a2/r;
    v = a3 - a4/r;
    if (c2 > .3)
    {
        s = (zp/r)*(1. + c2*(a1 + u + s2*v)/r);
        *lat = asin(s);
        ss = s*s;
        c = sqrt(1. - ss);
    }
    else
    {
        c = (w/r)*(1. - s2*(a5 - u - c2*v)/r);
        *lat = acos(c);
        ss = 1. - c*c;
        s = sqrt(ss);
    }
    g = 1. - e2*ss;
    rg = a/sqrt(g);
    rf = a6*rg;
    u = w - rg*c;
    v = zp - rf*s;
    f = c*u + s*v;
    m = c*v - s*u;
    p = m/(rf/g + f);
    *lat = *lat + p;
    *ht = f + m*p/2.;
    if (z < 0.)
        *lat = -*lat;
    return;
}
