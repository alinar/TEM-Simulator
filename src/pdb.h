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

#ifndef PDB_HEADER
#define PDB_HEADER
#include "macros.h"
#include "structs.h"

#define NUM_SCATTERING_PARAMETERS 5

int get_pot_from_pdb(particle *p, electronbeam *ed);

int get_scattering_parameters(double a[NUM_SCATTERING_PARAMETERS], double b[NUM_SCATTERING_PARAMETERS], int atomic_num);

#endif
