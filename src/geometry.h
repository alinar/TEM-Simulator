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

#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include "input.h"
#include "structs.h"


/***********************************************************************
 * Function:  new_geometry
 * Purpose:   Allocate and return a pointer to a new geometry struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new geometry struct.
 */

geometry *new_geometry(const char *name);


/***********************************************************************
 * Function:  delete_geometry
 * Purpose:   Free memory allocated by geometry struct.
 * Arguments: g - Pointer to geometry to be deleted.
 */

void delete_geometry(geometry *g);


/***********************************************************************
 * Function:  geometry_init
 * Purpose:   Do initializations necessary to make geometry struct ready
 *            for use in computations once parameters have been read from 
 *            input file. Calling geometry_init on a geometry object 
 *            which is already initialized has no effect. Once the 
 *            geometry object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: g - Pointer to geometry struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int geometry_init(geometry *g, 
                  simulation *sim);


/***********************************************************************
 * Function:  geometry_reset
 * Purpose:   Reset geometry object to the uninitialized state and 
 *            release memory allocated by geometry_init. Has no effect
 *            if geometry object is already in the uninitialized state.
 *            After geometry_reset has been called, geometry_init can 
 *            be used to initialize geometry object again, possibly with
 *            different input parameters.
 * Arguments: g - Pointer to geometry struct to be reset.
 */

void geometry_reset(geometry *g);


/***********************************************************************
 * Function:  geometry_write_log
 * Purpose:   Write input parameters of geometry object to the log file.
 * Arguments: g - Pointer to geometry struct.
 * Return:    0 on success, nonzero on failure.
 */

int geometry_write_log(geometry *g);


/***********************************************************************
 * Function:  set_proj_matrix
 * Purpose:   Set entries in projection matrix mapping the sample coordinate
 *            system to the microscope coordinate system for a given tilt.
 * Arguments: pm - Pointer to 3x3 matrix. On successful return is the 
 *            mapping matrix from sample coordinate system to microscope
 *            coordinate system.
 *            g - Pointer to geometry struct.
 *            tilt - Number of the image in tilt series.
 * Return:    0 on success, nonzero on failure.
 */

int set_proj_matrix(matrix *pm, 
                    geometry *g, 
                    long tilt);


/***********************************************************************
 * Function:  get_particle_geom
 * Purpose:   Get the mapping from particle coordinate system to microscope 
 *            coordinate system for a given particle and a given tilt.
 * Arguments: pm - Pointer to 3x3 matrix. On successful return is the 
 *            mapping matrix from particle coordinate system to 
 *            microscope coordinate system.
 *            pos - Array of 3 doubles. On successful return, contains the
 *            coordinates of the particle origin in the microscope 
 *            coordinate system.
 *            ps - Pointer to particleset struct to which the particle
 *            belongs.
 *            i - number of the particle in the particleset.
 *            g - Pointer to geometry struct.
 *            tilt - Number of the image in tilt series.
 * Return:    0 on success, nonzero on failure.
 */

int get_particle_geom(matrix *pm, 
                      double pos[3], 
                      particleset *ps, 
                      long i, 
                      geometry *g, 
                      long tilt);


/***********************************************************************
 * Function:  get_sample_geom
 * Purpose:   Get the mapping from sample coordinate system to microscope 
 *            coordinate system for a given tilt.
 * Arguments: pm - Pointer to 3x3 matrix. On successful return is the 
 *            mapping matrix from sample coordinate system to 
 *            microscope coordinate system.
 *            pos - Array of 3 doubles. On successful return, contains the
 *            coordinates of the sample origin in the microscope 
 *            coordinate system.
 *            s - Pointer to sample.
 *            i - number of the particle in the particleset.
 *            g - Pointer to geometry struct.
 *            tilt - Number of the image in tilt series.
 * Return:    0 on success, nonzero on failure.
 */

int get_sample_geom(matrix *pm, 
                    double pos[3], 
                    sample *s, 
                    geometry *g, 
                    long tilt);

#endif
