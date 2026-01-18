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

#ifndef OPENWME_TRAJECTORY_H
#define OPENWME_TRAJECTORY_H

#include <deque>
#include "vector.h"

/*
In this file we define objects that can be used to describe trajectories. These are
characterized by the fact that they describe location, velocity and acceleration of
an object at a time t. Furthermore, the passing of time can be described by calling
the function TimeStep.
*/

class TTrajectory
{
public:
	virtual ~TTrajectory();

	virtual TVector GetPosition(sim_double t) = 0;
	virtual TVector GetVelocity(sim_double t) = 0;
	virtual TVector GetAcceleration(sim_double t) = 0;

	virtual void TimeStep(sim_double dt, TVector a) = 0;
};

class TTimeBuffer
{
public:
	TTimeBuffer(void);
	TVector GetValue(sim_double t);
	void TimeStep(sim_double dt, TVector v);
	int max_history;

private:
	std::deque<sim_double> t;
	std::deque<TVector> v;
};

// TFreeTrajectory describes a freely moving object that changes speed and
// direction under the influence of forces.
class TFreeTrajectory: public TTrajectory
{
public:
	// r0 - starting postition at time 0
	// r0 - starting velocity at time 0
	// max_history - the number of time steps that the trajectory reaches in the past.
	//               the value should be greater than the largest distance in the simulation
	//               divided by the speed of light and the time step.
	// cons_kin_energy - if True, the kinetic energy of the particle is kept constant.
	//                   the parameter should usually be False.
	TFreeTrajectory(TVector r0, TVector v0, int max_history, bool cons_kin_energy);
	virtual ~TFreeTrajectory();

	virtual TVector GetPosition(sim_double t);
	virtual TVector GetVelocity(sim_double t);
	virtual TVector GetAcceleration(sim_double t);

	virtual void TimeStep(sim_double dt, TVector a);

private:
	bool cons_kin_energy;
	std::deque<sim_double> t;
	std::deque<TVector> r;
	std::deque<TVector> v;
	std::deque<TVector> a;
	int max_history;
};

// TLinearTrajectory represents the trajectory of a particle that always moves
// on a straight line
class TLinearTrajectory: public TTrajectory
{
public:
	// r - starting postition at time 0
	// v - constant velocity
	TLinearTrajectory(TVector r, TVector v);
	virtual ~TLinearTrajectory();

	virtual TVector GetPosition(sim_double t);
	virtual TVector GetVelocity(sim_double t);
	virtual TVector GetAcceleration(sim_double t);

	virtual void TimeStep(sim_double dt, TVector a);

private:
	TVector r;
	TVector v;
};

// TFixedTrajectory is intended to describe arbitrary but fixed trajectories
typedef TVector(*TTrajectoryFunc)(sim_double t, TVector r0);

class TFixedTrajectory: public TTrajectory
{
public:
	// r0 - starting postition at time 0
	// pf - function pointer to the trajectory r(t)
	// vf - function pointer to r'(t)
	// af - function pointer to r''(t)
	TFixedTrajectory(TVector r0, TTrajectoryFunc pf, TTrajectoryFunc vf, TTrajectoryFunc af);
	virtual ~TFixedTrajectory();

	virtual TVector GetPosition(sim_double t);
	virtual TVector GetVelocity(sim_double t);
	virtual TVector GetAcceleration(sim_double t);

	virtual void TimeStep(sim_double dt, TVector a);

private:
	TVector r0;
	TTrajectoryFunc pf;
	TTrajectoryFunc vf;
	TTrajectoryFunc af;
};

#endif
