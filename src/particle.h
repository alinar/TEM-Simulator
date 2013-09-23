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

#ifndef PARTICLE_HEADER
#define PARTICLE_HEADER
#include "array.h"
#include "matrix.h"
#include "structs.h"


/***********************************************************************
 * Function:  new_particle
 * Purpose:   Allocate and return a pointer to a new particle struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new particle struct.
 */

particle *new_particle(const char *name);


/***********************************************************************
 * Function:  delete_particle
 * Purpose:   Free memory allocated by particle struct.
 * Arguments: p - Pointer to particle to be deleted.
 */

void delete_particle(particle *p);


/***********************************************************************
 * Function:  particle_init
 * Purpose:   Do initializations necessary to make particle struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling particle_init on a particle object 
 *            which is already initialized has no effect. Once the 
 *            particle object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: p - Pointer to particle struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int particle_init(particle *p, 
                  simulation *sim);


/***********************************************************************
 * Function:  particle_reset
 * Purpose:   Reset particle object to the uninitialized state and 
 *            release memory allocated by particle_init. Has no effect
 *            if particle object is already in the uninitialized state.
 *            After particle_reset has been called, particle_init can 
 *            be used to initialize particle object again, possibly with
 *            different input parameters.
 * Arguments: p - Pointer to particle struct to be reset.
 */

void particle_reset(particle *p);


/***********************************************************************
 * Function:  particle_write_log
 * Purpose:   Write input parameters of particle object to the log file.
 * Arguments: p - Pointer to particle struct.
 */

int particle_write_log(particle *p);


/***********************************************************************
 * Function:  particle_project
 * Purpose:   Project potential map of one particle onto phase shift 
 *            plane in wavefunction struct.
 * Arguments: p - Pointer to particle struct.
 *            wf - Pointer to wavefunction struct.
 *            proj_matrix - 3x3 matrix mapping particle coordinate system
 *            to microscope coordinate system.
 *            pos - coordinates of particle origin.
 * Return:    0 on success, nonzero on failure.
 */

int particle_project(particle *p, 
                     wavefunction *wf, 
                     matrix *proj_matrix, 
                     double pos[3]);


/***********************************************************************
 * Function:  particle_hits_wavefunction
 * Purpose:   Check if projection of particle hits the area covered by
 *            wavefunction struct.
 * Arguments: p - Pointer to particle struct.
 *            wf - Pointer to wavefunction struct.
 *            pm - 3x3 matrix mapping particle coordinate system
 *            to microscope coordinate system.
 *            pos - coordinates of particle origin.
 * Return:    Nonzero if projection of particle hits wavefunction area, 
 *            0 otherwise.
 */

int particle_hits_wavefunction(particle *p, 
                               wavefunction *wf, 
                               matrix *pm, 
                               double pos[3]);

#endif
