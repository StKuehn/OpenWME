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
const sim_double radius = 40 * cm;
const sim_double wg_length = 2 * m;
const sim_double wg_diameter = 50 * cm;
const int wg_points = 40;
const sim_double reflpara = -5e14;
const sim_double refldelay = 0.35 * ns;
const int wgh_points = 10;
const sim_double wg_amp = 5 * cm;

sim_double AmplModFunc(sim_double t)
{
	if (t < -0.33 * ns) return 0;
	return 1;
}

int main()
{
	TDisplay display(600, 600, -1 * m, -1 * m, 1 * m, 1 * m, "Waveguide");

	TScene scene;
	std::deque<TParticle*> Field;

	// the transmitter
	std::deque<TParticle*> Transmitter;
	TParticle* HertzianDipole = scene.Add_Particle();
	HertzianDipole->ToHertzianDipole(mass, charge, amp, freq, 0, dir, AmplModFunc);
	HertzianDipole->SetLinearTrajectory(TVector(-wg_length / 2, 0, 0), TVector(0, 0, 0));
	Transmitter.push_back(HertzianDipole);

	// the waveguide
	std::deque<TParticle*> Waveguide;
	std::deque<TParticle*> Contour;
	for (int k = 0; k < wgh_points; k++)
	{
		for (int i = 0; i <= wg_points; i++)
		{
			TParticle* Reflector = scene.Add_Particle();
			Reflector->ToHertzianDipole(mass, charge, 0, 0, 0, TVector(0, 0, 1), NULL);
			Reflector->MakeReflective(reflpara, refldelay, 1000);
			sim_double x = i * wg_length / wg_points - wg_length / 2;
			sim_double c = wg_amp * sin(pi * i * wg_length / wg_points);
			sim_double y = c + wg_diameter / 2 * sin(2 * pi / wgh_points * k);
			sim_double z = c + wg_diameter / 2 * cos(2 * pi / wgh_points * k);
			Reflector->SetLinearTrajectory(TVector(x, y, z), TVector(0, 0, 0));
			Waveguide.push_back(Reflector);
			if ((k == 0) || (k == wgh_points / 2))
			{
				Contour.push_back(Reflector);
			}
		}
	}

	// create electromagnetic forces between transmitter and waveguide
	for (std::size_t j = 0; j < Waveguide.size(); j++)
	{
		for (std::size_t k = 0; k < Transmitter.size(); k++)
		{
			scene.Add_WeberMaxwellForce(Waveguide[j], Transmitter[k]);
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
			for (std::size_t j = 0; j < Waveguide.size(); j++)
			{
				scene.Add_WeberMaxwellForce(Waveguide[j], Probe);
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
		for (std::size_t j = 0; j < Contour.size(); j++)
		{
			display.DrawParticle(Contour[j], i * dt, 7, blue);
		}
		if (i * dt >= 10 * ns)
		{
			display.Draw(true);
		}
		else
		{
			display.Draw(false);
		}
	}

	return 0;
}
