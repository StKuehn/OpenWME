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
	virtual void Calculate(sim_double t, sim_double dt) = 0;

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
	virtual void Calculate(sim_double t, sim_double dt);

private:
	sim_double CalcTs(TParticle* src, TParticle* dst, sim_double t, sim_double dt);
	TVector CalcWeberMaxwellForce(sim_double q1, sim_double q2, sim_double t, sim_double tc, TVector rc, TVector vc, TVector ac);
	TVector CalcClassicalWeberForce(sim_double q1, sim_double q2, TVector r, TVector v);
	TVector CalcModernWeberForce(sim_double q1, sim_double q2, TVector r, TVector v);
	TVector CalcWeberForce(sim_double q1, sim_double q2, TVector r, TVector v);
};

/*
This object describes the force of a spring. The force acts in such a way that the distance between the two
particles tries to remain the same during the simulation.
*/

class THarmonicForce: public TForce
{
public:
	THarmonicForce(TParticle* p1, TParticle* p2, sim_double spring_constant, sim_double friction);
	virtual void Calculate(sim_double t, sim_double dt);

private:
	sim_double spring_constant;
	sim_double friction;
	sim_double nsp;
};

#endif
