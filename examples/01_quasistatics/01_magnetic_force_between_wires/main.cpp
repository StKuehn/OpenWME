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
const sim_double dt = 0.5 * ms;
// current strength
const sim_double current = 10 * A;
// number of current elements
const int segments = 30;
// length of the wires
const sim_double wire_length = 1 * m;
// distance between the two wires
const sim_double wire_distance = 0.2 * m;
// spring constant of the harmonic force between the current elements
const sim_double spring_constant = 1;
// friction parameter of the harmonic force
const sim_double friction = 0.9;

sim_double AmplModFunc1(sim_double t)
{
	sim_double freq = 0.5 * Hz;
	if (sin(2 * pi * freq * t) < 0) return -1;
	return 1;
}

sim_double AmplModFunc2(sim_double t)
{
	sim_double freq = 0.25 * Hz;
	if (sin(2 * pi * freq * t) < 0) return -1;
	return 1;
}

int main()
{
	TScene scene;
	std::deque<TParticle*> Wire1;
	std::deque<TParticle*> Wire2;

	// create the Wires
	for (int i = -segments / 2; i <= segments / 2; i++)
	{
		TParticle* CurrentElement1 = scene.Add_Particle();
		CurrentElement1->ToDCCurrentElement(0.01 * g, current, wire_length / segments, TVector(1, 0, 0), AmplModFunc1);
		Wire1.push_back(CurrentElement1);

		TParticle* CurrentElement2 = scene.Add_Particle();
		CurrentElement2->ToDCCurrentElement(0.01 * g, current, wire_length / segments, TVector(-1, 0, 0), AmplModFunc2);
		Wire2.push_back(CurrentElement2);

		if ((i == -segments / 2) || (i == segments / 2))
		{
			CurrentElement1->SetLinearTrajectory(TVector((wire_length / segments)*i, 0, -wire_distance / 2), TVector(0, 0, 0));
			CurrentElement2->SetLinearTrajectory(TVector((wire_length / segments)*i, 0, wire_distance / 2), TVector(0, 0, 0));
		}
		else
		{
			CurrentElement1->SetFreeTrajectory(TVector((wire_length / segments)*i, 0, -wire_distance / 2), TVector(0, 0, 0), 200, false);
			CurrentElement2->SetFreeTrajectory(TVector((wire_length / segments)*i, 0, wire_distance / 2), TVector(0, 0, 0), 200, false);
		}
	}

	// setup of the electromagnetic forces between the current elements
	for (std::size_t i = 0; i < Wire1.size(); i++)
	{
		for (std::size_t j = 0; j < Wire2.size(); j++)
		{
			scene.Add_WeberMaxwellForce(Wire1[i], Wire2[j]);
		}
	}

	// add harmonic forces between neighboring current elements
	for (std::size_t i = 0; i < Wire1.size() - 1; i++)
	{
		scene.Add_HarmonicForce(Wire1[i], Wire1[i + 1], spring_constant, friction);
	}
	for (std::size_t i = 0; i < Wire2.size() - 1; i++)
	{
		scene.Add_HarmonicForce(Wire2[i], Wire2[i + 1], spring_constant, friction);
	}

	// open a window to display the the setup
	TDisplay display(700, 350, -50 * cm, -25 * cm, 50 * cm, 25 * cm, "Magnetic force between current-carrying wires");

	for (int i = 0; ; i++)
	{
		// time step
		scene.TimeStep(dt);

		if (i % 10 == 0)
		{
			// display simulation
			display.Clear();

			// show the wires
			for (std::size_t j = 0; j < Wire1.size(); j++)
			{
				display.DrawParticle(Wire1[j], i * dt, 5, black);
			}
			for (std::size_t j = 0; j < Wire2.size(); j++)
			{
				display.DrawParticle(Wire2[j], i * dt, 5, black);
			}

			// show the current directions
			TVector color1 = red;
			if (Wire1[0]->GetVelocityCurElem(i * dt).x < 0) color1 = blue;
			for (std::size_t j = 0; j < Wire1.size(); j++)
			{
				display.PlotVector(Wire1[j]->GetVelocityCurElem(i * dt), Wire1[j]->GetPosition(i * dt), 10, 10, color1);
			}
			TVector color2 = red;
			if (Wire2[0]->GetVelocityCurElem(i * dt).x < 0) color2 = blue;
			for (std::size_t j = 0; j < Wire2.size(); j++)
			{
				display.PlotVector(Wire2[j]->GetVelocityCurElem(i * dt), Wire2[j]->GetPosition(i * dt), 10, 10, color2);
			}

			display.DrawTime(i * dt / s, "s", 1, black);

			display.Draw(true);
		}
	}

	return 0;
}
