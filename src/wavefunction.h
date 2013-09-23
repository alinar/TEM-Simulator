/*
 * Copyright 2008-2010, Hans Rullgard, Stockholm University and 
 * Lars-Goran Ofverstedt, Karolinska Institute
 *
 * This file is part of TEM Simulator.
 *
 * TEM Simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TEM Simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TEM Simulator.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef WAVEFUNCTION_HEADER
#define WAVEFUNCTION_HEADER
#include "structs.h"


wavefunction *new_wavefunction(electronbeam *ed, optics *o);

void delete_wavefunction(wavefunction *wf);

int wavefunction_init(wavefunction *wf, double det_pixel_size, double det_range[4]);

void wavefunction_reset(wavefunction *wf);

int wavefunction_prop_el_opt(wavefunction *wf, long tilt, double z);

int wavefunction_prop_el_vac(wavefunction *wf, double dist);

int wavefunction_set_incoming(wavefunction *wf);

int wavefunction_adj_phase(wavefunction *wf);

int wavefunction_prop_inel(wavefunction *wf, long tilt, double z);

int wavefunction_apply_bg_blur(wavefunction *wf, long tilt);

int wavefunction_set_inel_ctf(wavefunction *wf);

int wavefunction_propagate(simulation *sim, wavefunction *wf, double slice_th, long tilt);

#endif
