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
#include "forces.h"
#include "constants.h"

TForce::TForce()
{
	this->p1 = NULL;
	this->p2 = NULL;
}

TForce::~TForce()
{
}

TWeberMaxwellForce::TWeberMaxwellForce(TParticle* p1, TParticle* p2)
{
	this->p1 = p1;
	this->p2 = p2;
}

sim_double TWeberMaxwellForce::CalcTs(TParticle* src, TParticle* dst, sim_double t)
{
	// note that the fixed point iteration might not converge if relative velocities occur which
	// exceed the speed of light

	sim_double ts = t - nrm(dst->GetPosition(t) - src->GetPosition(t)) / c;
	sim_double tso;
	sim_double td = 1e30;
	sim_double tdo;

	do
	{
		tdo = td;
		tso = ts;
		TVector rd = dst->GetPosition(ts);
		TVector rs = src->GetPosition(ts);
		ts = t - nrm(rd - rs) / c;
		td = nrm(tso - ts);
	}
	while (tdo > td);

	return ts;
}

void TWeberMaxwellForce::Calculate(sim_double t)
{
	if ((p1 == NULL) || (p2 == NULL)) return;

	TVector r1, r2, v1, a1, v2, a2, rc, vc, ac, f11, f12, f21, f22;

	if ((p1->electro_dynamics) && (p2->electro_dynamics))
	{
		sim_double tc = CalcTs(p1, p2, t);
		rc = p2->GetPosition(tc) - p1->GetPosition(tc);
		vc = p2->GetVelocity(tc) - p1->GetVelocity(tc);
		ac = p2->GetAcceleration(tc) - p1->GetAcceleration(tc);

		r1 = p1->GetPositionCurElem(tc);
		v1 = p1->GetVelocityCurElem(tc);
		a1 = p1->GetAccelerationCurElem(tc);
		r2 = p2->GetPositionCurElem(tc);
		v2 = p2->GetVelocityCurElem(tc);
		a2 = p2->GetAccelerationCurElem(tc);

		f11 = CalcWeberMaxwellForce(p1->charge, p2->charge, t, tc, rc + r2 - r1, vc + v2 - v1, ac + a2 - a1);
		if (p1->current_element)
		{
			f21 = CalcWeberMaxwellForce(-p1->charge, p2->charge, t, tc, rc + r2 + r1, vc + v2 + v1, ac + a2 + a1);
		}
		if (p2->current_element)
		{
			f12 = CalcWeberMaxwellForce(p1->charge, -p2->charge, t, tc, rc - r2 - r1, vc - v2 - v1, ac - a2 - a1);
		}
		if (p1->current_element && p2->current_element)
		{
			f22 = CalcWeberMaxwellForce(-p1->charge, -p2->charge, t, tc, rc - r2 + r1, vc - v2 + v1, ac - a2 + a1);
		}
	}
	else
	{
		rc = p2->GetPosition(t) - p1->GetPosition(t);
		vc = p2->GetVelocity(t) - p1->GetVelocity(t);

		v1 = p1->GetVelocityCurElem(t);
		v2 = p2->GetVelocityCurElem(t);

		f11 = CalcWeberForce(p1->charge, p2->charge, rc, vc + v2 - v1);
		if (p1->current_element)
		{
			f21 = CalcWeberForce(-p1->charge, p2->charge, rc, vc + v2 + v1);
		}
		if (p2->current_element)
		{
			f12 = CalcWeberForce(p1->charge, -p2->charge, rc, vc - v2 - v1);
		}
		if (p1->current_element && p2->current_element)
		{
			f22 = CalcWeberForce(-p1->charge, -p2->charge, rc, vc - v2 + v1);
		}
	}

	if (!p2->probe)
	{
		p1->force = p1->force - (f11 + f12 + f21 + f22);
		p1->force_cur = p1->force_cur - (f11 + f12);
	}

	if (!p1->probe)
	{
		p2->force = p2->force + (f11 + f12 + f21 + f22);
		p2->force_cur = p2->force_cur + (f11 + f21);
	}
}

TVector TWeberMaxwellForce::CalcClassicalWeberForce(sim_double q1, sim_double q2, TVector r, TVector v)
{
	if ((q1 == 0) || (q2 == 0))
	{
		return TVector(0, 0, 0);

	}

	return (q1 * q2) / (4 * pi * eps0) * (1 + (v * v) / (c * c) - 3 / 2 * pow((r * v) / (c * nrm(r)), 2)) * r / pow(nrm(r), 3);
}

TVector TWeberMaxwellForce::CalcModernWeberForce(sim_double q1, sim_double q2, TVector r, TVector v)
{
	if ((q1 == 0) || (q2 == 0))
	{
		return TVector(0, 0, 0);
	}

	// The classical Weber force is an approximation of this formula for small velocities. This formula here
	// is in turn an approximation of the Weber-Maxwell force.
	return (q1 * q2 * r * sqrt(1 - v * v / (c * c))) / (4 * pi * eps0 * pow(r * r - pow(nrm(r ^ v) / c, 2), 1.5));
}

TVector TWeberMaxwellForce::CalcWeberForce(sim_double q1, sim_double q2, TVector r, TVector v)
{
	return CalcModernWeberForce(q1, q2, r, v);
}

TVector TWeberMaxwellForce::CalcRt(TVector rc, TVector vc, sim_double t, sim_double tc)
{
	return rc + vc * (t - tc);
}

sim_double TWeberMaxwellForce::CalcRh(TVector rt, TVector rc, TVector vc)
{
	return sqrt(rt * rt - ((rc ^ vc) * (rc ^ vc)) / (c * c));
}

TVector TWeberMaxwellForce::CalcWeberMaxwellForce(sim_double q1, sim_double q2, sim_double t, sim_double tc, TVector rc, TVector vc, TVector ac)
{
	if ((q1 == 0) || (q2 == 0))
	{
		return TVector(0, 0, 0);
	}

	TVector rt = CalcRt(rc, vc, t, tc);
	sim_double rh = CalcRh(rt, rc, vc);
	sim_double h1 = c * c - vc * vc;
	sim_double h2 = c * c * (t - tc) + rc * vc;
	return (q1 * q2 * (ac * (c - (rt * vc) / rh) + rt * (h1 * (h1 - rc * ac)) / (h2 * rh))) / (4 * pi * eps0 * c * c * h2 * pow((1 - vc * vc / (c * c)), 3 / 2));
}

THarmonicForce::THarmonicForce(TParticle* p1, TParticle* p2, sim_double spring_constant, sim_double friction)
{
	this->p1 = p1;
	this->p2 = p2;
	this->spring_constant = spring_constant;
	this->friction = friction;
	this->nsp = nrm(p1->GetPosition(0) - p2->GetPosition(0));
}

void THarmonicForce::Calculate(sim_double t)
{
	TVector dist = p1->GetPosition(t) - p2->GetPosition(t);
	TVector dvel = p1->GetVelocity(t) - p2->GetVelocity(t);

	sim_double nrmd = nrm(dist);
	if (nrmd < 0.01 * nm) nrmd = 0.01 * nm;
	TVector f = spring_constant * (dist - nsp * dist / nrmd);

	sim_double nrmv = nrm(dvel);
	if (nrmv > nm / s)
	{
		f = f + friction * nrm(f) * dvel / nrmv;
	}

	if (!p2->probe)
	{
		p1->force = p1->force - f;
		p1->force_cur = p1->force_cur - f / 2;
	}

	if (!p1->probe)
	{
		p2->force = p2->force + f;
		p2->force_cur = p2->force_cur + f / 2;
	}
}



