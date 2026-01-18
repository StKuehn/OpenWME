/*
Copyright (c) 2023 by Steffen Kühn, steffen.kuehn@aurinovo.de

This file is part of OpenWME, an electromagnetic field solver based on
Weber-Maxwell electrodynamics.

OpenWME is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.

OpenWME is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.
*/


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "vector.h"

TVector::TVector(void)
{
	x = 0;
	y = 0;
	z = 0;
}

TVector::TVector(sim_double x, sim_double y, sim_double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

TVector::TVector(const TVector& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
}

TVector::~TVector()
{

}

TVector operator+(const TVector &v1, const TVector &v2)
{
	TVector o = v1;

	o.x += v2.x;
	o.y += v2.y;
	o.z += v2.z;

	return o;
}

TVector operator-(const TVector &v1, const TVector &v2)
{
	TVector o = v1;

	o.x -= v2.x;
	o.y -= v2.y;
	o.z -= v2.z;

	return o;
}

sim_double operator*(const TVector &v1, const TVector &v2)
{
	sim_double o = 0;

	o += v1.x * v2.x;
	o += v1.y * v2.y;
	o += v1.z * v2.z;

	return o;
}

TVector operator*(sim_double v1, const TVector &v2)
{
	TVector o;

	o.x = v1 * v2.x;
	o.y = v1 * v2.y;
	o.z = v1 * v2.z;

	return o;
}

TVector operator*(const TVector &v1, sim_double v2)
{
	TVector o;

	o.x = v1.x * v2;
	o.y = v1.y * v2;
	o.z = v1.z * v2;

	return o;
}

TVector operator^(const TVector &v1, const TVector &v2)
{
	TVector o;

	o.x = -v1.z * v2.y + v1.y * v2.z;
	o.y = v1.z * v2.x - v1.x * v2.z;
	o.z = -v1.y * v2.x + v1.x * v2.y;

	return o;
}

TVector operator/(const TVector &v1, sim_double v2)
{
	TVector o;

	o.x = v1.x / v2;
	o.y = v1.y / v2;
	o.z = v1.z / v2;

	return o;
}

TVector operator-(const TVector &v)
{
	TVector o = v;

	o.x = -v.x;
	o.y = -v.y;
	o.z = -v.z;

	return o;
}

sim_double nrm(const TVector &v)
{
	return sqrt(v * v);
}

sim_double nrm(sim_double v)
{
	if (v > 0) return v;
	else return -v;
}

