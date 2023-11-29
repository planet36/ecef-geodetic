# ecef-geodetic

A collection of ECEF-to-geodetic coordinate conversion functions

Here are the abstract and recommendations sections from my paper _Comparison of Earth-Centered, Earth-Fixed to Geodetic Coordinate Conversion Algorithms_ published in the [2020 SIW](https://web.archive.org/web/20200813024039/https://www.sisostds.org/2020SIW.aspx) conference.

### Abstract

> The conversion of Earth-Centered, Earth-Fixed (ECEF) coordinates to geodetic latitude is inexact, and many algorithms for it have been developed in recent decades.
This paper compares the performance, as measured by accuracy and speed, of approximately 40 such algorithms using input locations ranging from the centroid of the World Geodetic System 1984 (WGS 84) ellipsoid to beyond lunar orbit.
Other factors such as code size, copyright, and licensing are presented, and suggestions of specific algorithms are made for various simulation use cases.


### Recommendations

>For the given input data, over a dozen algorithms yielded nanometer accuracy, which should be enough for everyone.™

>The **olson_1996** _(custom height)_ algorithm is recommended for ECEF-to-geodetic coordinate conversion. It was the fastest algorithm tested, of an acceptable size, and, being in the public domain, may be copied freely.

>When using an iterative algorithm, doing 2 iterations is recommended.

>If converting ECEF coordinates near the center of the Earth is required, the [GeographicLib](https://github.com/geographiclib/geographiclib) software library by Karney is recommended. It lists features such as conversion between ECEF, geodetic, UTM/UPS, and MGRS coordinates, as well as geoid height, geomagnetic, and geodesic calculations.

