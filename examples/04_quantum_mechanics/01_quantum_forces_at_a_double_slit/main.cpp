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
const sim_double dt = 0.05 * ns;

sim_double AmplModFunc(sim_double t)
{
	if (t < 0) return 0;
	return 1;
}

int main()
{
	TScene scene;
	std::deque<TParticle*> Field;
	// for the representation of one of the many double slits
	std::deque<TParticle*> FirstDoubleSlit;
	bool first = true;

	// create a grid of transmitters and double slits
	for (int ix = -25; ix <= 25; ix++)
	{
		for (int iz = -25; iz <= 25; iz++)
		{
			// create the transmitter (a copper bullet with a weight of 1g, in which 0.1 percent of the conduction
			// electrons are oscillating)
			TParticle* Transmitter = scene.Add_Particle();
			// convert the point charge into a Hertzian dipole
			Transmitter->ToHertzianDipole(1 * g, 9.47682e21 * e * 0.001, 1 * nm, 2.5 * GHz, 0, TVector(0, 0, 1), AmplModFunc);
			// set position
			Transmitter->SetLinearTrajectory(TVector(ix * 2 * cm, 0, iz * 2 * cm), TVector(0, 0, 0));
			// we want to know the force that would be generated on the transmitter by the reflected wave
			Field.push_back(Transmitter);

			// create a double slit exclusively for this transmitter (each transmitter gets its own double slit)
			std::deque<TParticle*> DoubleSlit;
			for (int iz = -60; iz <= 60; iz++)
			{
				if ((iz > -12) && (iz < -4)) continue;
				if ((iz <  12) && (iz >  4)) continue;
				TParticle* Reflector = scene.Add_Particle();
				// convert the point charge into a Hertzian dipole
				Reflector->ToHertzianDipole(1 * g, 9.47682e21 * e, 0, 0, 0, TVector(0, 0, 1), NULL);
				// make reflective
				Reflector->MakeReflective(1e-1, 100);
				// set position
				Reflector->SetLinearTrajectory(TVector(0 * cm, 0, iz * 1.5 * cm), TVector(0, 0, 0));
				DoubleSlit.push_back(Reflector);
				if (first) FirstDoubleSlit.push_back(Reflector);
			}
			first = false;

			// create electromagnetic forces between transmitter and double slit
			for (std::size_t j = 0; j < DoubleSlit.size(); j++)
			{
				scene.Add_WeberMaxwellForce(DoubleSlit[j], Transmitter);
			}
		}
	}

	// open a window to display the field and motion of the Hertzian dipole
	TDisplay display(1000, 1000, -50 * cm, -50 * cm, 50 * cm, 50 * cm, "\"Quantum\" forces at a double slit");

	for (int i = 0; ; i++)
	{
		// time step
		scene.TimeStep(dt);

		display.Clear();
		for (std::size_t j = 0; j < Field.size(); j++)
		{
			Field[j]->force = Field[j]->force_cumulative / (i + 1);
		}
		display.PlotForceFieldXZ(Field, 0, 3e3, 8, 8, black);
		for (std::size_t j = 0; j < FirstDoubleSlit.size(); j++)
		{
			display.DrawParticle(FirstDoubleSlit[j], 0, 5, blue);
		}
		display.DrawTime(i * dt / ns, "ns", 1, blue);
		display.Draw(true);
	}

	return 0;
}
