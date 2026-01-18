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

#ifndef OPENWME_PARTICLE_H
#define OPENWME_PARTICLE_H

#include <stdio.h>
#include "trajectory.h"

/*
In this file an object is defined, which can flexibly describe different types of point-like particles.
The simplest particle of this kind is a point charge. More complex particles are the Hertzian dipole and the
DC or AC current element. The correct usage of the different types is illustrated by the example programs.
*/

typedef sim_double(*TAmplModFunc)(sim_double t);

class TParticle
{
public:
	TParticle(void);
	~TParticle();

	// Turns the particle into a simple point charge.
	// mass - mass in kg
	// charge - electric charge in C
	void ToPointCharge(sim_double mass, sim_double charge);

	// Turns the particle into a Hertzian dipole. A Hertzian dipole is an electrically neutral
	// particle consisting of two inversely equal electric charges at the same location. Both charges
	// move sinusoidally with respect to each other and therefore radiate an electromagnetic wave.
	// Unlike the similar AC current element, in the case of the Hertzian dipole the spatial displacement
	// of the two charges with respect to each other is not assumed to be zero.
	// mass - mass in kg
	// charge - electric charge of the first charge in C
	// ampl - maximum spatial displacement of the first charge with respect to the center of mass in m
	// freq - frequency in Hz
	// phase - phase
	// dvec - direction in which the first point charge oscillates
	// ampmod - optional amplitude modulation, set NULL when not needed
	void ToHertzianDipole(sim_double mass, sim_double charge, sim_double ampl, sim_double freq,
						  sim_double phase, TVector dvec, TAmplModFunc ampmod);

	// Similar to ToHertzianDipole, but the charge quantities do not have to be inversely equal
	void ToChargedDipole(sim_double mass, sim_double charge1, sim_double charge2, sim_double ampl,
						 sim_double freq, sim_double phase, TVector dvec, TAmplModFunc ampmod);

	// Transforms the particle into a DC current element. A DC current element can be thought of as a
	// short segment of a wire in which a direct current flows. In reality, many DC current elements must
	// always be connected to a closed conductor loop.
	// mass - mass in kg
	// current - electric current in A
	// len - length of the wire segment
	// dvec - direction in which the first point charge moves
	// ampmod - optional amplitude modulation, set NULL when not needed
	void ToDCCurrentElement(sim_double mass, sim_double current, sim_double len, TVector dvec,
							TAmplModFunc ampmod);

	// Transforms the particle into a AC current element. An AC current element can be thought of as a
	// short segment of a wire in which a alternating current flows. In reality, low frequency AC current
	// elements must be connected to a closed conductor loop. Interrupted wires should be best modeled
	// with Hertzian dipoles.
	void ToACCurrentElement(sim_double mass, sim_double current, sim_double len, TVector dvec,
							sim_double freq, sim_double phase, TAmplModFunc ampmod);

	// Defines that the particle can move freely under the influence of forces (see trajectory.h).
	void SetFreeTrajectory(TVector r0, TVector v0, int max_history, bool cons_kin_energy, bool cons_location = false);

	// Defines that the particle moves only along a given line (see trajectory.h).
	void SetLinearTrajectory(TVector r0, TVector v0);

	// Defines that the particle moves along an arbitrarily shaped trajectory, but does not react
	// to forces (see trajectory.h).
	void SetFixedTrajectory(TVector r0, TTrajectoryFunc pf, TTrajectoryFunc vf, TTrajectoryFunc af);

	// Makes the particle reflective. The usage is explained in the examples.
	void MakeReflective(sim_double amplitude, sim_double delay, int max_history);

	// Performs a time step.
	void TimeStep(sim_double dt);

	// returns the position of the center of gravity at time t
	TVector GetPosition(sim_double t);

	// returns the velocity of the center of gravity at time t
	TVector GetVelocity(sim_double t);

	// returns the acceleration of the center of gravity at time t
	TVector GetAcceleration(sim_double t);

	// returns the position of the first charge of the current element relative to the center of gravity at time t
	TVector GetPositionCurElem(sim_double t);

	// returns the velocity of the first charge of the current element relative to the center of gravity at time t
	TVector GetVelocityCurElem(sim_double t);

	// returns the acceleration of the first charge of the current element relative to the center of gravity at time t
	TVector GetAccelerationCurElem(sim_double t);

	// set all forces to zero
	void ClearForces(void);

	// the net force
	TVector force;
	// force on the first charge for current elements
	TVector force_cur;
	// the cumulative net force
	TVector force_cumulative;

	// inertial mass in kg
	sim_double mass;

	// electric charge (in C). for the Hertzian dipole, this parameter corresponds to the amount of
	// positive charge.
	sim_double charge;

	// if this flag is set, the point charge becomes a current element, i.e. it consists of two point
	// charges of the same size at the same location but with different velocities and accelerations.
	bool current_element;

	// if this flag is True, the finite propagation speed of the electromagnetic force is fully taken
	// into account
	bool electro_dynamics;

	// this flag indicates that it is an AC current element
	bool ac_current_element;

	// constant speed of the moving charge (in m/s)
	sim_double v_dc;

	// frequency for modelling AC (in Hz)
	sim_double freq;

	// phase for modelling AC
	sim_double phase;

	// oscillation amplitude for modelling AC (in m)
	sim_double ampl;

	// direction vector of the motion of the moving charge
	TVector dvec;

	// if this flag is set, then this particle does not exert any force on other particles, but
	// it is able to receive some. such particles do not exist in reality and serve only for the
	// visualization of fields.
	bool probe;

	// optional amplitude modulations
	TAmplModFunc dc_ampmod;
	TAmplModFunc ac_ampmod;

private:
	// reflection parameters. define the strength and phase of the particle's reflection of incident waves
	sim_double reflection_amplitude;
	sim_double reflection_delay;

	// contains all information about the motion of the center of gravity
	TTrajectory* trj;

	TTimeBuffer force_cur_history;
};

#endif
