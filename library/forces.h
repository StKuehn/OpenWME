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

#ifndef OPENWME_FORCE_H
#define OPENWME_FORCE_H

#include <deque>
#include "particle.h"

/*
In this file we define objects that represent forces. Currently, there are only two types of forces.
More will be added gradually. For example, gravity.
*/

class TForce
{
public:
	TForce();
	virtual ~TForce();
	virtual void Calculate(sim_double t, sim_double dt, bool with_probes) = 0;
	virtual void Update(void) = 0;

protected:
	TParticle* p1;
	TParticle* p2;
};

/*
The Weber-Maxwell force is the electromagnetic force between two point charges or current elements. It describes
all electromagnetic effects, such as magnetism, induction, Lorentz force and also contains all effects related to
electromagnetic waves. Note that the Weber-Maxwell force is symmetric and obeys Newton's third law.

The deduction of the formula from Maxwell's equations is described here:
https://www.techrxiv.org/articles/preprint/24087840
*/

class TWeberMaxwellForce: public TForce
{
public:
	TWeberMaxwellForce(TParticle* p1, TParticle* p2);
	virtual void Calculate(sim_double t, sim_double dt, bool with_probes);
	virtual void Update(void);

private:
	sim_double CalcTs(TParticle* src, TParticle* dst, sim_double t, sim_double dt);
	TVector CalcWeberMaxwellForce(sim_double q1, sim_double q2, sim_double t, sim_double tc, const TVector& rc, const TVector& vc, const TVector& ac);
	TVector CalcClassicalWeberForce(sim_double q1, sim_double q2, const TVector& r, const TVector& v);
	TVector CalcModernWeberForce(sim_double q1, sim_double q2, const TVector& r, const TVector& v);
	TVector CalcWeberForce(sim_double q1, sim_double q2, const TVector& r, const TVector& v);

	TVector f11, f12, f21, f22;
};

/*
This object describes the force of a spring. The force acts in such a way that the distance between the two
particles tries to remain the same during the simulation.
*/

class THarmonicForce: public TForce
{
public:
	THarmonicForce(TParticle* p1, TParticle* p2, sim_double spring_constant, sim_double friction);
	virtual void Calculate(sim_double t, sim_double dt, bool with_probes);
	virtual void Update(void);

private:
	sim_double spring_constant;
	sim_double friction;
	sim_double nsp;
	TVector f;
};

class TPonderomotiveForce: public TForce
{
public:
	TPonderomotiveForce(TParticle* p, sim_double freq, sim_double freq_res, TVector a, bool spherical_wave);
	void AddReflector(const TVector& pos, sim_double refl_para);
	virtual void Calculate(sim_double t, sim_double dt, bool with_probes);
	virtual void Update(void);

private:
	TVector AkThi(int k, int i);
	TVector GradhkThi(int k, int i);
	void Update_h(int k);
	void UpdateR(sim_double t);

	bool spherical_wave;
	std::deque<TVector> positions;
	std::deque<TVector> r;
	std::deque<sim_double> rn;
	std::deque<sim_double> rn2;
	std::deque<sim_double> rn4;
	std::deque<sim_double> rp;
	std::deque<TVector> h;
	TVector f;
	sim_double w;
	sim_double wr;
	TVector av;
	sim_double q0;
	TVector r0;
};

#endif
