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

#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "scene.h"
#include "display.h"

const TVector black = TVector(0.25, 0.25, 0.25);
const TVector red = TVector(1, 0, 0);
const TVector blue = TVector(0, 0, 1);
// the time step of the simulation
const sim_double dt = 0.01 * s;
// current strength
const sim_double current = 1 * A;
// number of current elements
const int segments = 300;
// length of the wires
const sim_double wire_length = 30 * m;
// distance between the two wires
const sim_double wire_distance = 1 * m;

int main()
{
	TScene scene;
	std::deque<TParticle*> Wire1;
	std::deque<TParticle*> Wire2;

	// create the Wires
	for (int i = -segments / 2; i <= segments / 2; i++)
	{
		TParticle* CurrentElement1 = scene.Add_Particle();
		CurrentElement1->ToDCCurrentElement(0.01 * g, current, wire_length / segments, TVector(1, 0, 0), NULL);
		Wire1.push_back(CurrentElement1);

		TParticle* CurrentElement2 = scene.Add_Particle();
		CurrentElement2->ToDCCurrentElement(0.01 * g, current, wire_length / segments, TVector(-1, 0, 0), NULL);
		Wire2.push_back(CurrentElement2);

		CurrentElement1->SetLinearTrajectory(TVector((wire_length / segments)*i, 0, -wire_distance / 2), TVector(0, 0, 0));
		CurrentElement2->SetLinearTrajectory(TVector((wire_length / segments)*i, 0, wire_distance / 2), TVector(0, 0, 0));
	}

	// setup of the electromagnetic forces between the current elements
	for (std::size_t i = 0; i < Wire1.size(); i++)
	{
		for (std::size_t j = 0; j < Wire2.size(); j++)
		{
			scene.Add_WeberMaxwellForce(Wire1[i], Wire2[j]);
		}
	}

	// check if the strength of the force is equal to 2e-7 N
	scene.TimeStep(dt);
	TVector ft;
	for (std::size_t j = 0; j < Wire1.size(); j++)
	{
		ft = ft + Wire1[j]->force;
	}
	ft = ft / wire_length;
	printf("%.3e N\n", (float)nrm(ft));

	return 0;
}
