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
const TVector blue = TVector(0, 0, 1);
const TVector red = TVector(1, 0, 0);
// the time step of the simulation
const sim_double dt = 1 * ns;
// current strength
const sim_double current = 10 * A;
// number of current elements
const int segments = 100;
// radius of the Helmholtz coil
const sim_double radius = 40 * cm;

int main()
{
	TScene scene;
	std::deque<TParticle*> HelmholtzCoil;

	// create the Helmholtz coil
	for (int i = 0; i < segments; i++)
	{
		TParticle* MetallAtoms = scene.Add_Particle();
		TVector r(radius * sin(2 * pi / segments * i), 0, radius * cos(2 * pi / segments * i));
		TVector d(cos(2 * pi / segments * i), 0, -sin(2 * pi / segments * i));
		MetallAtoms->SetLinearTrajectory(r, TVector(0, 0, 0));
		MetallAtoms->ToDCCurrentElement(1 * g, current, (2 * pi * radius) / segments, d, NULL);
		HelmholtzCoil.push_back(MetallAtoms);
	}

	// create a negative point charge with an elementary charge that moves with 0.1 percent
	// of the speed of light
	TParticle* PointCharge = scene.Add_Particle();
	PointCharge->ToPointCharge(me, -e);
	PointCharge->SetFreeTrajectory(TVector(0, 0, 10 * cm), TVector(0.001 * c, 0, 0), 1000, false);

	// define an electromagnetic force between each current element of the Helmholtz coil
	// and the point charge
	for (std::size_t i = 0; i < HelmholtzCoil.size(); i++)
	{
		scene.Add_WeberMaxwellForce(HelmholtzCoil[i], PointCharge);
	}

	// open a window to display the the setup
	TDisplay display(700, 700, -50 * cm, -50 * cm, 50 * cm, 50 * cm, "Point Charge in a Helmholtz Coil (Calculated without Lorentz force)");

	for (int i = 0; ; i++)
	{
		// time step
		scene.TimeStep(dt);

		if (i % 10 == 0)
		{
			// display simulation
			display.Clear();

			// show Helmholtz Coil with current direction
			for (std::size_t j = 0; j < HelmholtzCoil.size(); j++)
			{
				display.DrawParticle(HelmholtzCoil[j], i * dt, 7, black);
			}
			for (std::size_t j = 0; j < HelmholtzCoil.size(); j++)
			{
				display.PlotVector(HelmholtzCoil[j]->dvec * 10, HelmholtzCoil[j]->GetPosition(i * dt), 10, 10, blue);
			}

			display.DrawParticle(PointCharge, i * dt, 5, red);
			display.DrawTime(i * dt / us, "us", 1, black);

			display.Draw(true);
		}
	}

	return 0;
}
