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

#ifndef OPENWME_DISPLAY_H
#define OPENWME_DISPLAY_H

#include <cairo.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

#include "scene.h"

/*
In this file, a simple object is provided for graphical representation.
Internally the object is based on Cairo (https://www.cairographics.org/).
*/

class TDisplay
{
public:
	TDisplay(int SizeX, int SizeY, sim_double x1,  sim_double y1, sim_double x2,  sim_double y2, const char* label);
	~TDisplay();

	void Clear(void);
	void PlotVector(TVector v, TVector r, sim_double mscf, sim_double maxh, TVector color);
	void PlotForceFieldXZ(const std::deque<TParticle*>& field, sim_double t, sim_double amp, sim_double mscf, sim_double maxh, TVector color);
	void DrawTime(sim_double t, const char* unit, sim_double size, TVector color);
	void DrawParticle(TParticle* p, sim_double t, sim_double size, TVector color);
	void Draw(bool write_to_file);

private:
	void mirror(cairo_t* cr);
	void draw_arrow(cairo_t* cr, TVector v, TVector r, sim_double scf, sim_double mscf, sim_double maxh);
	sim_double x_to_win(sim_double x);
	sim_double z_to_win(sim_double y);
	void copy_buffer_to_screen(cairo_surface_t* cs_buffer);
	void handle_events(void);
	void WriteToFile(void);

	int SizeX;
	int SizeY;
	sim_double x1;
	sim_double y1;
	sim_double x2;
	sim_double y2;
	Display* dis;
	Window win;
	cairo_surface_t* cs_buffer;
	cairo_t* cr_buffer;
	int png_nbr;
};

#endif

