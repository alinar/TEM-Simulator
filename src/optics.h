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

#ifndef OPTICS_HEADER
#define OPTICS_HEADER
#include "input.h"
#include "structs.h"


/***********************************************************************
 * Function:  new_optics
 * Purpose:   Allocate and return a pointer to a new optics struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new optics struct.
 */

optics *new_optics(const char *name);


/***********************************************************************
 * Function:  delete_optics
 * Purpose:   Free memory allocated by optics struct.
 * Arguments: o - Pointer to optics to be deleted.
 */

void delete_optics(optics *o);


/***********************************************************************
 * Function:  optics_init
 * Purpose:   Do initializations necessary to make optics struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling optics_init on a optics object 
 *            which is already initialized has no effect. Once the 
 *            optics object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: o - Pointer to optics struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int optics_init(optics *o, 
                simulation *sim);


/***********************************************************************
 * Function:  optics_reset
 * Purpose:   Reset optics object to the uninitialized state and 
 *            release memory allocated by optics_init. Has no effect
 *            if optics object is already in the uninitialized state.
 *            After optics_reset has been called, optics_init can 
 *            be used to initialize optics object again, possibly with
 *            different input parameters.
 * Arguments: o - Pointer to optics struct to be reset.
 */

void optics_reset(optics *o);


/***********************************************************************
 * Function:  optics_write_log
 * Purpose:   Write input parameters of optics object to the log file.
 * Arguments: o - Pointer to optics struct.
 */

int optics_write_log(optics *o);


/***********************************************************************
 * Function:  optics_get_aperture
 * Purpose:   Get the diameter of the aperture in the focal plane, 
 *            converted to standard length unit.
 * Arguments: o - Pointer to optics struct.
 * Return:    Diameter of the aperture in the focal plane, converted to 
 *            standard length unit.
 */

double optics_get_aperture(optics *o);


/***********************************************************************
 * Function:  optics_get_focal_length
 * Purpose:   Get the focal length of the primary lens, converted to 
 *            standard length unit.
 * Arguments: o - Pointer to optics struct.
 * Return:    Focal length of primary lens, converted to standard length 
 *            unit.
 */

double optics_get_focal_length(optics *o);


/***********************************************************************
 * Function:  optics_get_cond_ap_angle
 * Purpose:   Get the aperture angle of the condenser lens converted to
 *            radians.
 * Arguments: o - Pointer to optics struct.
 * Return:    Aperture angle of the condenser lens, converted to radians.
 */

double optics_get_cond_ap_angle(optics *o);


/***********************************************************************
 * Function:  optics_get_cs
 * Purpose:   Get the spherical aberration converted to standard length 
 *            unit.
 * Arguments: o - Pointer to optics struct.
 * Return:    Spherical aberration converted to standard length unit.
 */

double optics_get_cs(optics *o);


/***********************************************************************
 * Function:  optics_get_cc
 * Purpose:   Get the chromatic aberration converted to standard length 
 *            unit.
 * Arguments: o - Pointer to optics struct.
 * Return:    Chromatic aberration converted to standard length unit.
 */

double optics_get_cc(optics *o);


/***********************************************************************
 * Function:  optics_get_defocus
 * Purpose:   Get the defocus in image number tilt, converted to standard
 *            length unit.
 * Arguments: o - Pointer to optics struct.
 *            tilt - Number of image in tilt series, zero-based.
 * Return:    Defocus converted to standard length unit.
 */

double optics_get_defocus(optics *o, 
                          long tilt);

#endif
