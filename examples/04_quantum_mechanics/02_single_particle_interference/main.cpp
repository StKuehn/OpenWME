/*
Copyright (c) 2026 by Steffen Kühn, steffen.kuehn@aurinovo.de

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

const sim_double h = 6.626e-34; // Planck constant

const TVector black = TVector(0.25, 0.25, 0.25);
const TVector green = TVector(0, 1, 0);
const TVector red = TVector(1, 0, 0);
const TVector blue = TVector(0, 0, 1);

const TVector ex = TVector(1, 0, 0);
const TVector ez = TVector(0, 0, 1);
const TVector zv = TVector(0, 0, 0);

// time step size
const sim_double dt = 5e-17 * s;
// parameter that specifies how many time steps of a particle's history are stored
const int his = 10000;
// determines how many particles are simulated
const int nz = 100;
// defines after how many time steps the display is updated (set to zero if nothing should be displayed)
const int plotsw = 3;
// determines whether the graphics are saved as PNG files
const bool produce_pngs = false;
// defines the size of the plotted area
const sim_double plotrng = 5 * nm;
// experimental
const bool spherical_wave = false;

// center-to-center distance between the two openings
const sim_double od = 1.5 * nm;
// size of an opening
const sim_double os = 1.0 * nm;
// distance of the detector from the barrier
const sim_double L = 0.1 * m;
// starting position of a particle
const sim_double x0 = -5.0 * nm;
// size parameters of the barrier
const int barrier_nbr_z = 50;
const int barrier_nbr_y = 10;

// set this value TRUE if only a plot of the ponderomotive field is needed
const bool plot_field_only = false;
const bool plot_field = plot_field_only || ((plotsw != 0) & true);
// gain factor of the force so that the representation of the force arrows becomes visible
const sim_double plotgain = 3e11;

typedef struct
{
	sim_double mass;
	sim_double speed;
	sim_double em_freq;
	sim_double charge;
	sim_double reflp;
	TVector amp;
} TSimPara;

std::deque<TVector> barrier;

static void calc_pond_field(const TSimPara& simpara, const std::deque<TVector>& barrier, std::deque<TParticle*>& field)
{
	const int n = 75;
	static TScene plot_scene;

	for (int ix = -n; ix <= n; ix++)
	{
		for (int iz = -n; iz <= n; iz++)
		{
			TParticle* probe = plot_scene.Add_Particle();
			probe->ToPointCharge(simpara.mass, simpara.charge);
			probe->SetLinearTrajectory(TVector(ix * plotrng / n, 0, iz * plotrng / n), zv);
			TPonderomotiveForce* pondf = plot_scene.Add_PonderomotiveForce(probe, simpara.em_freq, 0, simpara.amp, spherical_wave);
			for (std::size_t j = 0; j < barrier.size(); j++)
			{
				pondf->AddReflector(barrier[j], simpara.reflp);
			}
			field.push_back(probe);
		}
	}

	plot_scene.TimeStep(0, true);
}

static void create_barrier(void)
{
	for (int ir = -barrier_nbr_z / 2; ir <= barrier_nbr_z / 2; ir++)
	{
		for (int jr = -barrier_nbr_y / 2; jr <= barrier_nbr_y / 2; jr++)
		{
			sim_double zp = ir * 0.2 * nm;
			sim_double yp = jr * 0.2 * nm;
			bool slit1 = ((zp >= -od / 2 - os / 2) && (zp <= -od / 2 + os / 2));
			bool slit2 = ((zp >= od / 2 - os / 2) && (zp <=  od / 2 + os / 2));
			if (slit1 || slit2) continue;
			barrier.push_back(TVector(0, yp, zp));
		}
	}
}

static void mainloop(const TSimPara& simpara)
{
	bool sc = true;
	FILE* fpd = fopen("detector.txt", "w");
	FILE* fpt = fopen("trajectories.txt", "w");
	TScene* scene;

	// calculate the starting positions
	std::deque<TVector> startpos;
	for (int spi = 0; spi < nz; spi ++)
	{
		sim_double A = 1.5 * nm;
		sim_double B = -1.5 * nm;
		sim_double z0 = A + ((B - A) * spi) / (nz - 1);
		startpos.push_back(TVector(x0, 0, z0));
	}

	printf("speed = %f c\n", (float)(simpara.speed / c));
	printf("frequency = %f 10^17 Hz\n", (float)(2 * simpara.em_freq) / 1e17);
	printf("wavelength = %f nm\n", (float)(c / (2 * simpara.em_freq)) / nm);
	printf("barrier N = %i\n", (int)barrier.size());

	TDisplay display(1200, 1200, -plotrng, -plotrng, plotrng, plotrng, "single partice interference/diffraction");

	std::deque<TParticle*> field;
	if (plot_field) calc_pond_field(simpara, barrier, field);

	if (plot_field_only)
	{
		display.Clear();
		if (plot_field) display.PlotForceFieldXZ(field, 0, plotgain, 7, 7, black);
		for (std::size_t j = 0; j < barrier.size(); j++)
		{
			display.DrawDisc(barrier[j], 8, blue);
		}
		display.Draw(true);
		return;
	}

	printf("%i trajectories\n", (int)startpos.size());
	for (std::size_t ip = 0; ip < startpos.size(); ip++)
	{
		printf("\rprocess trajectory %i", (int)(ip + 1));
		fflush(stdout);

		scene =  new TScene;
		TParticle* particle = scene->Add_Particle();
		particle->ToPointCharge(simpara.mass, simpara.charge);
		particle->SetFreeTrajectory(startpos[ip], ex * simpara.speed, his, false, false);
		TPonderomotiveForce* pondf = scene->Add_PonderomotiveForce(particle, simpara.em_freq, 0, simpara.amp, spherical_wave);

		for (std::size_t ib = 0; ib < barrier.size(); ib++)
		{
			pondf->AddReflector(barrier[ib], simpara.reflp);
		}

		if (ip == 0)
		{
			fprintf(fpd, "detector={");
		}

		// simulation
		for (int i = 0;; i++)
		{
			bool evaluate = (plotsw != 0) ? (i % plotsw == 0) : false;
			// time step
			scene->TimeStep(dt, evaluate);

			// display simulation
			if (evaluate)
			{
				display.Clear();
				if (plot_field) display.PlotForceFieldXZ(field, 0, plotgain, 7, 7, black);
				display.DrawTime(i * dt * 1e15, "fs", 1, blue);
				display.DrawParticle(particle, i * dt, 8, green);
				for (std::size_t j = 0; j < (std::size_t)i; j++)
				{
					display.DrawParticle(particle, j * dt, 2, green);
				}
				for (std::size_t j = 0; j < barrier.size(); j++)
				{
					display.DrawDisc(barrier[j], 8, blue);
				}
				display.Draw(produce_pngs);
			}

			TVector pos = particle->GetPosition(i * dt);
			TVector vel = particle->GetVelocity(i * dt);

			if (i == 0)
			{
				if (ip != 0) fprintf(fpt, "};\n");
				fprintf(fpt, "trajectory[%i]={", (int)ip);
			}
			else
			{
				fprintf(fpt, ",");
			}
			fprintf(fpt, "{%f,%f,%f}", (float)(i * dt / 1e-15), (float)(pos.x / nm), (float)(pos.z / nm));

			if (pos.x > -x0)
			{
				if (!sc) fprintf(fpd, ",");
				sc = false;
				float v = (float)((vel.z * (L - pos.x)) / (vel.x));
				fprintf(fpd, "{%f,%f}", (float)startpos[ip].z / nm, v);
				fflush(fpd);
				break;
			}
			if (pos.x < x0) break;
			if (i * dt >  4 * nrm(x0) / simpara.speed) break;
		}
		delete scene;
	}

	fprintf(fpd, "};");
	fclose(fpd);

	fprintf(fpt, "};");
	fclose(fpt);

	printf("\ndone\n");
}

int main(void)
{
	TSimPara simpara;

	create_barrier();

	// physical parameters
	simpara.mass = me;
	simpara.charge = 1.65e12 * e;
	simpara.speed = 0.004 * c;
	simpara.em_freq = 0.5 * c / h * simpara.mass * simpara.speed;
	simpara.amp = TVector(1 * m / (s * s), 0, 0);
	simpara.reflp = simpara.charge * simpara.charge / me;

	mainloop(simpara);

	return 0;
}
