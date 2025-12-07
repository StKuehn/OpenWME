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

#ifndef OPENWME_SCENE_H
#define OPENWME_SCENE_H

#include "forces.h"

/*
In this file an object is defined that contains all particles, current elements and forces of a setup.
Usually, a single object of this type is the central object of a simulation.
*/

class TScene
{
public:
	TScene(void);
	~TScene();

	// Adds an empty particle whose properties can be modified later. There are different types,
	// such as point particles, Hertzian dipoles, DC current elements or AC current elements.
	TParticle* Add_Particle(void);

	// This function adds a special point particle that has charge 1 and mass 1 and can receive
	// forces but not emit them. Such objects do not exist in nature, one needs them only for the
	// representation of fields. If only the Weber-Maxwell force acts on a probe, then the force corresponds
	// to the numerical value of the electromagnetic field strength with the unit V/m.
	// Note that probes can move. They measure then the fields in the corresponding inertial frame.
	// r - starting position
	// v - constant velocity (usually zero)
	TParticle* Add_Probe(TVector r, TVector v);

	// This function connects the two given particles with an electromagnetic force. Note that the
	// electromagnetic force is always symmetric, so the order of the two particles is irrelevant. This
	// symmetry corresponds to Newton's third law.
	TWeberMaxwellForce* Add_WeberMaxwellForce(TParticle* p1, TParticle* p2);

	// This function adds a spring force. The two given particles therefore try to maintain their initial
	// distance in the course of the simulation. Here, too, the order of the particles does not matter,
	// because a spring force, like the electromagnetic force, obeys Newton's third law.
	THarmonicForce* Add_HarmonicForce(TParticle* p1, TParticle* p2, sim_double spring_constant, sim_double friction);

	// Calcuates all forces and performs a time step.
	void TimeStep(sim_double dt, bool with_probes = true);

	// Returns the current time of the simulation.
	sim_double GetCurrentTime(void);

private:
	sim_double t;
	std::deque<TForce*> Forces;
	std::deque<TParticle*> Particles;
};

#endif

