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

#ifndef OPENWME_CONSTANTS_H
#define OPENWME_CONSTANTS_H

/*
In this file, physical constants and units are defined.
*/

#include "openwme_types.h"

// Meter
const sim_double m  = 1;
const sim_double km = 1e+3 * m;
const sim_double cm = 0.01 * m;
const sim_double mm = 1e-3 * m;
const sim_double um = 1e-6 * m;
const sim_double nm = 1e-9 * m;

// Second
const sim_double s  = 1;
const sim_double ms = 1e-3 * s;
const sim_double us = 1e-6 * s;
const sim_double ns = 1e-9 * s;

// Hertz
const sim_double Hz  = 1 / s;
const sim_double kHz = 1e3 * Hz;
const sim_double MHz = 1e6 * Hz;
const sim_double GHz = 1e9 * Hz;

// Kilogramm
const sim_double kg = 1;
const sim_double g = 1e-3 * kg;
const sim_double mg = 1e-6 * kg;
const sim_double ug = 1e-9 * kg;

// Ampere
const sim_double A = 1;

// Newton
const sim_double N = kg * m / (s*s);
const sim_double mN = 1e-3 * N;
const sim_double uN = 1e-6 * N;
const sim_double nN = 1e-9 * N;
const sim_double pN = 1e-12 * N;

// Coulomb
const sim_double C = (A*s);

// Volt
const sim_double V = (kg*m*m) / (A*s*s*s);

// speed of light in vacuum
const sim_double c = 299792458 * m / s;

// vacuum permittivity
const sim_double eps0 = 8.854187817e-12 * (A*s) / (V*m);

// vacuum permeability
const sim_double mu0 = 1 / (eps0*c*c);

// elementary charge
const sim_double e = 1.602176565e-19 * C;
// mass electron
const sim_double me = 9.10938356e-31 * kg;

// mass proton
const sim_double mp = 1.672621898e-27 * kg;

const sim_double pi = 3.14159265358979;

#endif
