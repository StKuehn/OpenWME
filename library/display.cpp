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
#include <stdlib.h>
#include <math.h>
#include <cairo-xlib.h>
#include <cairo-svg.h>
#include <sys/stat.h>
#include <stdio.h>
#include <termios.h>
#include <sys/resource.h>

#include "display.h"

TDisplay::TDisplay(int SizeX, int SizeY, sim_double x1,  sim_double y1, sim_double x2,  sim_double y2, const char* label)
{
	this->SizeX = SizeX;
	this->SizeY = SizeY;
	this->x1 = x1;
	this->x2 = x2;
	this->y1 = y1;
	this->y2 = y2;

	cs_buffer = NULL;
	cr_buffer = NULL;

	png_nbr = 0;

	dis = XOpenDisplay(NULL);
	win = XCreateSimpleWindow(dis, RootWindow(dis, 0), 1, 1, this->SizeX, this->SizeY, 0, WhitePixel(dis, 0), WhitePixel(dis, 0));
	XMapWindow(dis, win);
	XStoreName(dis, win, label);
}

TDisplay::~TDisplay()
{
	if (cs_buffer != NULL)
	{
		cairo_destroy(cr_buffer);
		cr_buffer = NULL;
	}

	if (cr_buffer != NULL)
	{
		cairo_surface_destroy(cs_buffer);
		cs_buffer = NULL;
	}

	XCloseDisplay(dis);
}

void TDisplay::Clear(void)
{
	if ((cs_buffer != NULL) || (cr_buffer != NULL))
	{
		cairo_destroy(cr_buffer);
		cairo_surface_destroy(cs_buffer);
		handle_events();
		cs_buffer = NULL;
		cr_buffer = NULL;
	}

	//cs_buffer = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, this->SizeX, this->SizeY);
	cs_buffer = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, this->SizeX, this->SizeY);
	cr_buffer = cairo_create(cs_buffer);
	cairo_set_source_rgb(cr_buffer, 1, 1, 1);
	cairo_paint(cr_buffer);
	mirror(cr_buffer);
}

void TDisplay::Draw(bool write_to_file)
{
	if (cs_buffer == NULL)
	{
		return;
	}

	cairo_surface_t* cs_screen = cairo_xlib_surface_create(dis, win, DefaultVisual(dis, 0),  this->SizeX,  this->SizeY);
	cairo_t* cr_screen = cairo_create(cs_screen);
	cairo_set_source_surface(cr_screen, cs_buffer, 0, 0);
	cairo_paint(cr_screen);
	cairo_destroy(cr_screen);
	cairo_surface_destroy(cs_screen);

	if (write_to_file)
	{
		WriteToFile();
	}
}

void TDisplay::WriteToFile(void)
{
	if (cs_buffer == NULL)
	{
		return;
	}

	// to JPG: find . -name "*.png" -print0 | xargs -0 mogrify -format jpg -quality 50
	// to MP4: ffmpeg -r 50 -i %05d_output.png test.mp4
	// to WEBM: ffmpeg -r 50 -i %05d_output.png -c:v libvpx -crf 10 -b:v 5M -c:a libvorbis test.webm
	static int png_nbr = 1;
	char filename [100];

	mkdir("./png", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	snprintf(filename, sizeof(filename), "./png/%5.5i_output.png", png_nbr++);
	cairo_surface_write_to_png(cs_buffer, filename);
}

void TDisplay::PlotVector(TVector v, TVector r, sim_double mscf, sim_double maxh, TVector color)
{
	TVector ri(x_to_win(r.x), 0, z_to_win(r.z));
	cairo_set_source_rgb(cr_buffer, color.x, color.y, color.z);
	cairo_set_line_width(cr_buffer, 0.5);
	draw_arrow(cr_buffer, v, ri, 1, mscf, maxh);
}

void TDisplay::PlotForceFieldXZ(const std::deque<TParticle*>& field, sim_double t, sim_double amp, sim_double mscf, sim_double maxh, TVector color)
{
	cairo_set_source_rgb(cr_buffer, color.x, color.y, color.z);
	cairo_set_line_width(cr_buffer, 0.5);

	for (std::size_t i = 0; i < field.size(); i++)
	{
		TVector r, ri, f;
		r = field[i]->GetPosition(t);
		f = field[i]->force;
		ri.x = x_to_win(r.x);
		ri.y = 0;
		ri.z = z_to_win(r.z);
		draw_arrow(cr_buffer, f, ri, amp, mscf, maxh);
	}
}

void TDisplay::DrawTime(sim_double t, const char* unit, sim_double size, TVector color)
{
	char str[100];
	cairo_set_source_rgb(cr_buffer, color.x, color.y, color.z);
	cairo_select_font_face(cr_buffer, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr_buffer, 20 * size);
	cairo_move_to(cr_buffer, 10, SizeY - 30);
	snprintf(str, sizeof(str), "%.2f %s", (double)t, unit);
	mirror(cr_buffer);
	cairo_show_text(cr_buffer, str);
	mirror(cr_buffer);
}

void TDisplay::DrawParticle(TParticle* p, sim_double t, sim_double size, TVector color)
{
	TVector r = p->GetPosition(t);
	DrawDisc(r, size, color);
}

void TDisplay::DrawDisc(TVector pos, sim_double size, TVector color)
{
	cairo_set_source_rgb(cr_buffer, color.x, color.y, color.z);
	cairo_arc(cr_buffer, x_to_win(pos.x), z_to_win(pos.z), size, 0, 2 * M_PI);
	cairo_fill(cr_buffer);
}

void TDisplay::mirror(cairo_t* cr)
{
	cairo_translate(cr, 0, this->SizeY / 2);
	cairo_scale(cr, 1, -1);
	cairo_translate(cr, 0, -this->SizeY / 2);
}

void TDisplay::draw_arrow(cairo_t* cr, TVector v, TVector r, sim_double scf,
						  sim_double mscf, sim_double maxh)
{
	TVector s, e, f;

	f.x = scf * v.x;
	f.y = 0;
	f.z = scf * v.z;

	sim_double nrm = sqrt(f * f);

	if (nrm > mscf)
	{
		scf = (mscf * scf) / nrm;
		f.x = scf * v.x;
		f.z = scf * v.z;
		nrm = mscf;
	}

	sim_double headlen = nrm;

	if (headlen > maxh) headlen = maxh;

	sim_double rad = 0.3;

	s = r - f;
	e = r + f;

	sim_double angle = atan2(e.z - s.z, e.x - s.x) + M_PI;

	sim_double xl = e.x + headlen * cos(angle - rad);
	sim_double yl = e.z + headlen * sin(angle - rad);
	sim_double xr = e.x + headlen * cos(angle + rad);
	sim_double yr = e.z + headlen * sin(angle + rad);

	cairo_move_to(cr, s.x + 0.5, s.z + 0.5);
	cairo_line_to(cr, e.x + 0.5, e.z + 0.5);
	cairo_stroke(cr);

	cairo_line_to(cr, xl + 0.5, yl + 0.5);
	cairo_line_to(cr, xr + 0.5, yr + 0.5);
	cairo_line_to(cr, e.x + 0.5, e.z + 0.5);
	cairo_stroke_preserve(cr);
	cairo_fill(cr);
}

sim_double TDisplay::x_to_win(sim_double x)
{
	return (SizeX * x1) / (x1 - x2) - (SizeX / (x1 - x2)) * x;
}

sim_double TDisplay::z_to_win(sim_double y)
{
	return (SizeY * y1) / (y1 - y2) - (SizeY / (y1 - y2)) * y;
}

void TDisplay::handle_events(void)
{
	XEvent event;
	XNextEvent(dis, &event);
}
