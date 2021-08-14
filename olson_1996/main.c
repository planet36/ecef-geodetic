// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

#include <math.h>
#include <stdio.h>

#define DEG_PER_RAD (180 / M_PI)

extern void latlon(
	double x, double y, double z,
	double* lat, double* lon, double* ht);

int main()
{
	while (1)
	{
		double x = 0, y = 0, z = 0, lat = 0, lon = 0, ht = 0;
		if (fscanf(stdin, "%lf %lf %lf", &x, &y, &z) == EOF)
			break;
		latlon(x, y, z, &lat, &lon, &ht);
		printf("%.12f %.12f %.12f\n", lat * DEG_PER_RAD, lon * DEG_PER_RAD, ht);
	}
	return 0;
}
