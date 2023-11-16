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

// these parameters represent a transmitter with an ERP of approximately 100W (see calibration)
const sim_double maxtime = 10 * ns;
const sim_double dt = 0.05 * ns;
const sim_double freq = 3 * GHz;
const sim_double amp = 100 * nm;
const sim_double charge = 1e14 * e;
const TVector direction = TVector(0, 0, 1);
const sim_double reflpara = 10000;
const sim_double mass = 1 * g;

sim_double AmplModFunc(sim_double t)
{
	if (t < 0) return 0;
	return 1;
}

int main()
{
	TScene scene;
	std::deque<TParticle*> Field;
	std::deque<TParticle*> FirstDoubleSlit;
	bool first = true;

	// create a grid of transmitters and double slits
	for (int ix = -25; ix <= 25; ix++)
	{
		for (int iz = -25; iz <= 25; iz++)
		{
			TParticle* Transmitter = scene.Add_Particle();
			Transmitter->ToHertzianDipole(mass, charge, amp, freq, 0, direction, AmplModFunc);
			Transmitter->SetLinearTrajectory(TVector(ix * 2 * cm, 0, iz * 2 * cm), TVector(0, 0, 0));
			Field.push_back(Transmitter);

			// create a double slit exclusively for this transmitter (each transmitter gets its own double slit)
			std::deque<TParticle*> DoubleSlit;
			for (int iz = -60; iz <= 60; iz++)
			{
				sim_double z = iz * 1.5 * cm;

				if ((z > -18 * cm) && (z < -6 * cm)) continue;
				if ((z <  18 * cm) && (z >  6 * cm)) continue;

				TParticle* Reflector = scene.Add_Particle();
				Reflector->ToHertzianDipole(1, 1, 0, 0, 0, TVector(0, 0, 1), NULL);
				Reflector->SetLinearTrajectory(TVector(0, 0,  z), TVector(0, 0, 0));
				Reflector->MakeReflective(reflpara, 0, 1000);
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
		display.PlotForceFieldXZ(Field, 0, 1e10, 8, 8, black);
		for (std::size_t j = 0; j < FirstDoubleSlit.size(); j++)
		{
			display.DrawParticle(FirstDoubleSlit[j], 0, 5, blue);
		}
		display.DrawTime(i * dt / ns, "ns", 1, blue);
		if (i * dt >= maxtime)
		{
			display.Draw(true);
			break;
		}
		else
		{
			display.Draw(false);
		}
	}

	// output the result in a form that can be read by Mathematica
	FILE* fp = fopen("result.txt", "w");
	fprintf(fp, "Field = {");
	for (std::size_t j = 0; j < Field.size(); j++)
	{
		if (j != 0) fprintf(fp, ",");
		fprintf(fp, "{{%.5f,%.5f}, {%.5f,%.5f}}",
				(double)Field[j]->GetPosition(0).x,
				(double)Field[j]->GetPosition(0).z,
				(double)(Field[j]->force.x / nN),
				(double)(Field[j]->force.z / nN));
	}
	fprintf(fp, "};\n");
	fclose(fp);

	return 0;
}
