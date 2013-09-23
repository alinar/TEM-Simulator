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

#ifndef PARTICLESET_HEADER
#define PARTICLESET_HEADER
#include "structs.h"


/***********************************************************************
 * Function:  new_particleset
 * Purpose:   Allocate and return a pointer to a new particleset struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new particleset struct.
 */

particleset *new_particleset(const char *name);


/***********************************************************************
 * Function:  delete_particleset
 * Purpose:   Free memory allocated by particleset struct.
 * Arguments: ps - Pointer to particleset to be deleted.
 */

void delete_particleset(particleset *ps);


/***********************************************************************
 * Function:  particleset_init
 * Purpose:   Do initializations necessary to make particleset struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling particleset_init on a particleset object 
 *            which is already initialized has no effect. Once the 
 *            particleset object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: ps - Pointer to particleset struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int particleset_init(particleset *ps, 
                     simulation *sim);


/***********************************************************************
 * Function:  particleset_reset
 * Purpose:   Reset particleset object to the uninitialized state and 
 *            release memory allocated by particleset_init. Has no effect
 *            if particleset object is already in the uninitialized state.
 *            After particleset_reset has been called, particleset_init can 
 *            be used to initialize particleset object again, possibly with
 *            different input parameters.
 * Arguments: ps - Pointer to particleset struct to be reset.
 */

void particleset_reset(particleset *ps);


/***********************************************************************
 * Function:  particleset_write_log
 * Purpose:   Write input parameters of particleset object to the log file.
 * Arguments: ps - Pointer to particleset struct.
 */

int particleset_write_log(particleset *ps);

int rotate_particle(matrix *pm, 
                    particleset *ps, 
                    long i);

int get_particle_pos(double pos[3], 
                     particleset *ps, 
                     long i);

int get_particle_coord(matrix *pm, 
                       double pos[3], 
                       particleset *ps, 
                       long i);

#endif
