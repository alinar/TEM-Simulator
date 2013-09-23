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


#ifndef VOLUME_HEADER
#define VOLUME_HEADER
#include "input.h"
#include "matrix.h"
#include "structs.h"


/***********************************************************************
 * Function:  new_volume
 * Purpose:   Allocate and return a pointer to a new volume struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new volume struct.
 */

volume *new_volume(const char *name);


/***********************************************************************
 * Function:  delete_volume
 * Purpose:   Free memory allocated by volume struct.
 * Arguments: v - Pointer to volume to be deleted.
 */

void delete_volume(volume *v);


/***********************************************************************
 * Function:  volume_init
 * Purpose:   Do initializations necessary to make volume struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling volume_init on a volume object 
 *            which is already initialized has no effect. Once the 
 *            volume object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: v - Pointer to volume struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int volume_init(volume *v, simulation *sim);


/***********************************************************************
 * Function:  volume_reset
 * Purpose:   Reset volume object to the uninitialized state and 
 *            release memory allocated by volume_init. Has no effect
 *            if volume object is already in the uninitialized state.
 *            After volume_reset has been called, volume_init can 
 *            be used to initialize volume object again, possibly with
 *            different input parameters.
 * Arguments: v - Pointer to volume struct to be reset.
 */

void volume_reset(volume *v);


/***********************************************************************
 * Function:  volume_write_log
 * Purpose:   Write input parameters of volume object to the log file.
 * Arguments: v - Pointer to volume struct.
 */

int volume_write_log(volume *v);

int volume_write_maps(volume *v);

int volume_get_potential(volume *v, simulation *sim);

#endif
