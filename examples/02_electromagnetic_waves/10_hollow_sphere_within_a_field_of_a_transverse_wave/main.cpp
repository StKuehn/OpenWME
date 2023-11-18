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
const TVector green = TVector(0, 1, 0);
const TVector red = TVector(1, 0, 0);
const TVector blue = TVector(0, 0, 1);

const sim_double dt = 0.05 * ns;
const sim_double freq = 1 * GHz;
const sim_double amp = 1 * nm;
const sim_double charge = 1e14 * e;
const sim_double mass = 1 * g;
const sim_double maxtime = 20 * ns;
const TVector dir = TVector(0, 0, 1);
const sim_double radius = 40*cm;
const int shield_points = 100;
const sim_double reflpara = -3.5e18;
const sim_double refldelay = 0.35 * ns;

sim_double AmplModFunc(sim_double t)
{
	if (t < -0.33*ns) return 0;
	return 1;
}

int main()
{
	TDisplay display(600, 600, -1 * m, -1 * m, 1 * m, 1 * m, "Hollow sphere within a field of a transverse wave");

	TScene scene;
	std::deque<TParticle*> Field;

	// the transmitter
	std::deque<TParticle*> Transmitter;
	for (int i = -15; i <= 15; i++)
	{
		TParticle* HertzianDipole = scene.Add_Particle();
		HertzianDipole->ToHertzianDipole(mass, charge, amp, freq, 0, dir, AmplModFunc);
		HertzianDipole->SetLinearTrajectory(TVector(-1.1*m, 0, i * 10 * cm), TVector(0, 0, 0));
		Transmitter.push_back(HertzianDipole);
	}

	// the shield
	std::deque<TParticle*> Shield;
	for (int i = 0; i <= shield_points; i++)
	{
		TParticle* Reflector = scene.Add_Particle();
		Reflector->ToHertzianDipole(mass, charge/shield_points, 0, 0, 0, TVector(0, 0, 1), NULL);
		Reflector->MakeReflective(reflpara, refldelay, 1000);
		Reflector->SetLinearTrajectory(TVector(radius*sin((2*pi*i)/shield_points), 0, radius*cos((2*pi*i)/shield_points)), TVector(0, 0, 0));
		Shield.push_back(Reflector);
	}

	// create electromagnetic forces between transmitter and shield
	for (std::size_t j = 0; j < Shield.size(); j++)
	{
		for (std::size_t k = 0; k < Transmitter.size(); k++)
		{
			scene.Add_WeberMaxwellForce(Shield[j], Transmitter[k]);
		}
	}
	// create electromagnetic forces between the parts of the shield
	for (std::size_t j = 0; j < Shield.size(); j++)
	{
		for (std::size_t k = 0; k < Shield.size(); k++)
		{
			// avoid nearfield interaction between shield parts
			if (nrm(Shield[j]->GetPosition(0) - Shield[k]->GetPosition(0)) < 1.9*radius) continue;
			scene.Add_WeberMaxwellForce(Shield[j], Shield[k]);
		}
	}

	// create a grid of test charges to represent the field
	for (int ix = -50; ix <= 50; ix++)
	{
		for (int iz = -50; iz <= 50; iz++)
		{
			TParticle* Probe = scene.Add_Probe(TVector(ix * 2 * cm, 0, iz * 2 * cm), TVector(0, 0, 0));
			Field.push_back(Probe);
			for (std::size_t j = 0; j < Transmitter.size(); j++)
			{
				scene.Add_WeberMaxwellForce(Transmitter[j], Probe);
			}
			for (std::size_t j = 0; j < Shield.size(); j++)
			{
				scene.Add_WeberMaxwellForce(Shield[j], Probe);
			}
		}
	}

	// simulation
	for (int i = 0; i * dt <= maxtime; i++)
	{
		// time step
		scene.TimeStep(dt);

		// display simulation
		display.Clear();
		display.PlotForceFieldXZ(Field, i * dt, 1e1, 10, 10, black);
		display.DrawTime(i * dt / ns, "ns", 1, blue);
		for (std::size_t j = 0; j < Transmitter.size(); j++)
		{
			display.DrawParticle(Transmitter[j], i * dt, 5, green);
		}
		for (std::size_t j = 0; j < Shield.size(); j++)
		{
			display.DrawParticle(Shield[j], i * dt, 7, blue);
		}
		display.Draw(true);
	}

	return 0;
}
