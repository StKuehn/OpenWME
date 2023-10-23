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
#include "particle.h"

TParticle::TParticle(void)
{
	this->mass = 1 * kg;
	this->charge = 1 * C;
	this->reflection_parameter = 0;
	this->probe = false;
	this->v_dc = 0;
	this->dvec = TVector(0, 0, 0);
	this->ampl = 0;
	this->freq = 0;
	this->phase = 0;
	this->current_element = false;
	this->electro_dynamics = true;
	this->ac_current_element = false;
	this->ac_ampmod = NULL;
	this->dc_ampmod = NULL;
	this->trj = new TLinearTrajectory(TVector(0, 0, 0), TVector(0, 0, 0));
	this->force = TVector(0, 0, 0);
	this->force_cur = TVector(0, 0, 0);
	this->force_cumulative = TVector(0, 0, 0);
}

TParticle::~TParticle()
{
	delete this->trj;
}

void TParticle::ToPointCharge(sim_double mass, sim_double charge)
{
	this->mass = mass;
	this->charge = charge;
}

void TParticle::ToHertzianDipole(sim_double mass, sim_double charge, sim_double ampl, sim_double freq, sim_double phase, TVector dvec, TAmplModFunc ampmod)
{
	this->mass = mass;
	this->charge = charge;
	this->ampl = ampl;
	this->freq = freq;
	this->phase = phase;
	this->ac_ampmod = ampmod;
	this->electro_dynamics = true;
	this->current_element = true;
	this->ac_current_element = false;
	this->dvec = dvec / nrm(dvec);
}

void TParticle::ToDCCurrentElement(sim_double mass, sim_double current, sim_double len, TVector dvec, TAmplModFunc ampmod)
{
	// too small values causes numerical problems
	this->v_dc = 1000 * m / s;
	this->mass = mass;
	this->charge = (current * nrm(len)) / (2 * this->v_dc);
	this->dc_ampmod = ampmod;
	this->electro_dynamics = false;
	this->current_element = true;
	this->ac_current_element = false;
	this->dvec = dvec / nrm(dvec);
}

void TParticle::ToACCurrentElement(sim_double mass, sim_double current, sim_double len, TVector dvec, sim_double freq, sim_double phase, TAmplModFunc ampmod)
{
	sim_double v_dc = 1000 * m / s;
	this->mass = mass;
	this->ampl = v_dc / (2 * pi * freq);
	this->freq = freq;
	this->phase = phase;
	this->charge = (current * nrm(len)) / (2 * v_dc);
	this->ac_ampmod = ampmod;
	this->electro_dynamics = true;
	this->current_element = true;
	this->ac_current_element = true;
	this->dvec = dvec / nrm(dvec);
}

void TParticle::ClearForces(void)
{
	force = TVector(0, 0, 0);
	force_cur = TVector(0, 0, 0);
}

void TParticle::SetFreeTrajectory(TVector r0, TVector v0, int max_history, bool cons_kin_energy)
{
	delete this->trj;
	this->trj = new TFreeTrajectory(r0, v0, this->force_cur_history.max_history, cons_kin_energy);
}

void TParticle::SetLinearTrajectory(TVector r0, TVector v0)
{
	delete this->trj;
	this->trj = new TLinearTrajectory(r0, v0);
}

void TParticle::SetFixedTrajectory(TVector r0, TTrajectoryFunc pf, TTrajectoryFunc vf, TTrajectoryFunc af)
{
	delete this->trj;
	this->trj = new TFixedTrajectory(r0, pf, vf, af);
}

void TParticle::MakeReflective(sim_double reflection_parameter, int max_history)
{
	this->reflection_parameter = reflection_parameter;
	force_cur_history.max_history = max_history;
}

TVector TParticle::GetPosition(sim_double t)
{
	return this->trj->GetPosition(t);
}

TVector TParticle::GetVelocity(sim_double t)
{
	return this->trj->GetVelocity(t);
}

TVector TParticle::GetPositionCurElem(sim_double t)
{
	TVector r;
	if ((current_element) && (electro_dynamics) && (!ac_current_element))
	{
		// only for Hertzian Dipols
		sim_double ac_amp = ampl;
		if (ac_ampmod) ac_amp *= ac_ampmod(t);
		r = dvec * ac_amp * sin(2 * pi * freq * t + phase);
	}
	return r;
}

TVector TParticle::GetVelocityCurElem(sim_double t)
{
	TVector v;
	if (current_element)
	{
		if (electro_dynamics)
		{
			sim_double ac_amp = ampl;
			if (ac_ampmod) ac_amp *= ac_ampmod(t);
			v = dvec * 2 * ac_amp * freq * pi * cos(2 * pi * freq * t + phase);
		}
		else
		{
			v = dvec * v_dc;
			if (dc_ampmod) v = v * dc_ampmod(t);
		}
	}
	return v;
}

TVector TParticle::GetAccelerationCurElem(sim_double t)
{
	TVector a;
	if (current_element)
	{
		if (electro_dynamics)
		{
			sim_double ac_amp = ampl;
			if (ac_ampmod) ac_amp *= ac_ampmod(t);
			a = -4 * ac_amp * dvec * freq * freq * pi * pi * sin(2 * pi * freq * t + phase);
			if (reflection_parameter != 0)
			{
				// this is a heuristic and should be further improved
				TVector fch = force_cur_history.GetValue(t);
				a = a + reflection_parameter * fch;
			}
		}
	}
	return a;
}

TVector TParticle::GetAcceleration(sim_double t)
{
	return this->trj->GetAcceleration(t);
}

void TParticle::TimeStep(sim_double dt)
{
	if (reflection_parameter != 0)
	{
		force_cur_history.TimeStep(dt, force_cur);
	}
	force_cumulative = force_cumulative + force;
	trj->TimeStep(dt, force / this->mass);
}

