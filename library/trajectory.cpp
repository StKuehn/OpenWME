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

#include <math.h>
#include <stdio.h>
#include "trajectory.h"
#include "constants.h"

static std::size_t SearchT(sim_double t, std::deque<sim_double>& tb)
{
	std::size_t start = 0;
	std::size_t end = tb.size();

	for (;;)
	{
		std::size_t pos = (end - start) / 2 + start;

		if ((end - start) < 4)
		{
			break;
		}

		if (t < tb[pos])
		{
			end = pos;
		}
		else
		{
			start = pos;
		}
	}

	for (std::size_t i = start; i < end; i++)
	{
		if (t < tb[i])
		{
			return i;
		}
	}

	return end;
}

static TVector Interpolate(sim_double t, std::deque<sim_double>& tb, std::deque<TVector>* r, std::deque<TVector>* v)
{
	std::size_t i = SearchT(t, tb);

	if (i == 0)
	{
		// extrapolate
		TVector rc = TVector(0, 0, 0);

		if (v != NULL)
		{
			rc = (*v)[0] * (tb[i] - t);
		}

		return (*r)[0] - rc;
	}
	else if (i == tb.size())
	{
		// extrapolate
		TVector rc = TVector(0, 0, 0);

		if (v != NULL)
		{
			rc = (*v)[i - 1] * (t - tb[i - 1]);
		}

		return (*r)[i - 1] + rc;
	}
	else
	{
		// interpolate
		sim_double t0 = tb[i - 1];
		sim_double t1 = tb[i];
		sim_double rel = (t1 - t) / (t1 - t0);
		return rel * (*r)[i - 1] + (1 - rel) * (*r)[i];
	}
}

TTrajectory::~TTrajectory()
{
}

TFreeTrajectory::TFreeTrajectory(TVector r0, TVector v0, int max_history, bool cons_kin_energy)
{
	this->cons_kin_energy = cons_kin_energy;
	this->t.push_back(0);
	this->r.push_back(r0);
	this->v.push_back(v0);
	this->a.push_back(TVector(0, 0, 0));
	this->max_history = max_history;
}

TFreeTrajectory::~TFreeTrajectory()
{
}

TVector TFreeTrajectory::GetPosition(sim_double t)
{
	return Interpolate(t, this->t, &this->r, &this->v);
}

TVector TFreeTrajectory::GetVelocity(sim_double t)
{
	return Interpolate(t, this->t, &this->v, &this->a);
}

TVector TFreeTrajectory::GetAcceleration(sim_double t)
{
	return Interpolate(t, this->t, &this->a, NULL);
}

void TFreeTrajectory::TimeStep(sim_double dt, TVector a)
{
	// see: https://www.compadre.org/PICUP/resources/Numerical-Integration/
	// velocity Verlet algorithm
	sim_double tn = this->t.back();
	TVector rn = this->r.back();
	TVector vn = this->v.back();
	TVector an = this->a.back();

	sim_double tnp = tn + dt;
	TVector anp = a;
	TVector rnp = rn + vn * dt + 0.5 * an * dt * dt;
	TVector vnp = vn + 0.5 * (anp + an) * dt;

	if (cons_kin_energy)
	{
		vnp = nrm(vn) / nrm(vnp) * vnp;
	}

	this->t.push_back(tnp);
	this->r.push_back(rnp);
	this->v.push_back(vnp);
	this->a.push_back(anp);

	if (this->t.size() > (std::size_t)max_history)
	{
		this->t.pop_front();
		this->r.pop_front();
		this->v.pop_front();
		this->a.pop_front();
	}
}

TLinearTrajectory::TLinearTrajectory(TVector r, TVector v)
{
	this->r = r;
	this->v = v;
}

TLinearTrajectory::~TLinearTrajectory()
{
}

TVector TLinearTrajectory::GetPosition(sim_double t)
{
	return r + v * t;
}

TVector TLinearTrajectory::GetVelocity(sim_double t)
{
	return v;
}

TVector TLinearTrajectory::GetAcceleration(sim_double t)
{
	return TVector(0, 0, 0);
}

void TLinearTrajectory::TimeStep(sim_double dt, TVector a)
{
}

TFixedTrajectory::TFixedTrajectory(TVector r0, TTrajectoryFunc pf, TTrajectoryFunc vf, TTrajectoryFunc af)
{
	this->r0 = r0;
	this->pf = pf;
	this->vf = vf;
	this->af = af;
}

TFixedTrajectory::~TFixedTrajectory()
{
}

TVector TFixedTrajectory::GetPosition(sim_double t)
{
	return pf(t, r0);
}

TVector TFixedTrajectory::GetVelocity(sim_double t)
{
	return vf(t, r0);
}

TVector TFixedTrajectory::GetAcceleration(sim_double t)
{
	return af(t, r0);
}

void TFixedTrajectory::TimeStep(sim_double dt, TVector a)
{
}

TTimeBuffer::TTimeBuffer(void)
{
	t.push_back(0);
	v.push_back(TVector(0, 0, 0));
}

TVector TTimeBuffer::GetValue(sim_double t)
{
	return Interpolate(t, this->t, &v, NULL);
}

void TTimeBuffer::TimeStep(sim_double dt, TVector v)
{
	sim_double tn = this->t.back() + dt;
	this->t.push_back(tn);
	this->v.push_back(v);
	if (this->t.size() > (std::size_t)max_history)
	{
		this->t.pop_front();
		this->v.pop_front();
	}
}



