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

#ifndef OPENWME_VECTOR_H
#define OPENWME_VECTOR_H

#include "openwme_types.h"

/*
A class for convenient representation of equations with vectors.
*/

class TVector
{
public:
	TVector(void);
	TVector(sim_double x, sim_double y, sim_double z);
	TVector(const TVector& v);
	~TVector();

	friend TVector operator+(const TVector &v1, const TVector &v2);
	friend TVector operator-(const TVector &v1, const TVector &v2);
	friend TVector operator-(const TVector &v);
	friend sim_double operator*(const TVector &v1, const TVector &v2);
	friend TVector operator*(sim_double v1, const TVector &v2);
	friend TVector operator*(const TVector &v1, sim_double v2);
	friend TVector operator/(const TVector &v1, sim_double v2);

	// Cross product
	friend TVector operator^(const TVector &v1, const TVector &v2);

	// Vector norm
	friend sim_double nrm(const TVector &v);

	sim_double x;
	sim_double y;
	sim_double z;
};

sim_double nrm(sim_double v);

#endif
