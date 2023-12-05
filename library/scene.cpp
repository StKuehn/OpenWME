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

#include "scene.h"

TScene::TScene(void)
{
	this->t = 0;
}

TScene::~TScene()
{
	for (std::size_t i = 0; i < Forces.size(); i++)
	{
		delete Forces[i];
	}

	for (std::size_t i = 0; i < Particles.size(); i++)
	{
		delete Particles[i];
	}
}

void TScene::TimeStep(sim_double dt)
{
	for (std::size_t i = 0; i < Particles.size(); i++)
	{
		Particles[i]->ClearForces();
	}

	for (std::size_t i = 0; i < Forces.size(); i++)
	{
		Forces[i]->Calculate(this->t, dt);
	}

	for (std::size_t i = 0; i < Particles.size(); i++)
	{
		Particles[i]->TimeStep(dt);
	}

	this->t += dt;
}

TWeberMaxwellForce* TScene::Add_WeberMaxwellForce(TParticle* p1, TParticle* p2)
{
	TWeberMaxwellForce* WeberMaxwellForce = new TWeberMaxwellForce(p1, p2);
	Forces.push_back(WeberMaxwellForce);
	return WeberMaxwellForce;
}

THarmonicForce* TScene::Add_HarmonicForce(TParticle* p1, TParticle* p2, sim_double spring_constant, sim_double friction)
{
	THarmonicForce* HarmonicForce = new THarmonicForce(p1, p2, spring_constant, friction);
	Forces.push_back(HarmonicForce);
	return HarmonicForce;
}

TParticle* TScene::Add_Particle(void)
{
	Particles.push_back(new TParticle);
	return Particles.back();
}

TParticle* TScene::Add_Probe(TVector r, TVector v)
{
	TParticle* particle = Add_Particle();
	particle->SetLinearTrajectory(r, v);
	particle->probe = true;
	particle->current_element = false;
	return particle;
}

sim_double TScene::GetCurrentTime(void)
{
	return this->t;
}
