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

const sim_double maxtime = 10 * ns;
const sim_double dt = 0.05 * ns;
const sim_double freq = 3 * GHz;
const sim_double amp = 100 * nm;
const sim_double charge = 1e14 * e;
const TVector direction = TVector(0, 0, 1);
const sim_double reflpara = 10000;
const sim_double mass = 1 * g;
const int n = 30;

sim_double AmplModFunc(sim_double t)
{
	if (t < 0) return 0;
	return 1;
}

int main()
{
	TScene scene;
	std::deque<TParticle*> Mirror;

	// create a transmitter at location x = -0.5 m
	TParticle* Transmitter1 = scene.Add_Particle();
	Transmitter1->ToHertzianDipole(mass, charge, amp, freq, 0, direction, AmplModFunc);
	Transmitter1->SetLinearTrajectory(TVector(-0.5 * m, 0, 0), TVector(0, 0, 0));
	// create a probe that is one meter distant from the transmitter
	TParticle* RefProbe1 = scene.Add_Probe(TVector(0.5 * m, 0, 0), TVector(0, 0, 0));
	// we want to know the RMS value of the field strength at a distance of one meter
	scene.Add_WeberMaxwellForce(RefProbe1, Transmitter1);

	// create a mirror at location x = 0.0 m
	TParticle* Transmitter2 = scene.Add_Particle();
	Transmitter2->ToHertzianDipole(mass, charge, amp, freq, 0, direction, AmplModFunc);
	Transmitter2->SetLinearTrajectory(TVector(-0.5 * m, 0, 0), TVector(0, 0, 0));
	for (int imz = -n; imz <= n; imz++)
	{
		for (int imy = -n; imy <= n; imy++)
		{
			TParticle* Reflector = scene.Add_Particle();
			Reflector->ToHertzianDipole(1, 1, 0, 0, 0, TVector(0, 0, 1), NULL);
			Reflector->SetLinearTrajectory(TVector(0, imy * 2 * cm, imz * 2 * cm), TVector(0, 0, 0));
			// make reflective
			Reflector->MakeReflective(reflpara, 1000);
			Mirror.push_back(Reflector);
		}
	}
	// we want to compare the RMS value of the reflected field
	TParticle* RefProbe2 = scene.Add_Probe(TVector(-0.5 * m, 0, 0), TVector(0, 0, 0));
	for (std::size_t j = 0; j < Mirror.size(); j++)
	{
		scene.Add_WeberMaxwellForce(Mirror[j], Transmitter2);
		scene.Add_WeberMaxwellForce(Mirror[j], RefProbe2);
	}

	int steps = 0;
	sim_double rms1 = 0;
	sim_double rms2 = 0;
	for (; steps * dt <= maxtime; steps++)
	{
		// time step
		scene.TimeStep(dt);
		rms1 += RefProbe1->force * RefProbe1->force;
		rms2 += RefProbe2->force * RefProbe2->force;
	}
	rms1 = sqrt(rms1 / steps);
	rms2 = sqrt(rms2 / steps);

	printf("RMS (transmitter): %.2f V/m, ERP: %.2f W, RMS (mirror) %.2f V/m, ratio: %.3f\n",
		   (double)rms1, (double)(pow((rms1 * 1 * m), 2) / (30 * 1.64)), (double)rms2, (double)(rms2 / rms1));

	return 0;
}
