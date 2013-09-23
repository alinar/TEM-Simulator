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


#ifndef SAMPLE_HEADER
#define SAMPLE_HEADER
#include "structs.h"


/***********************************************************************
 * Function:  new_sample
 * Purpose:   Allocate and return a pointer to a new sample struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new sample struct.
 */

sample *new_sample(const char *name);


/***********************************************************************
 * Function:  delete_sample
 * Purpose:   Free memory allocated by sample struct.
 * Arguments: s - Pointer to sample to be deleted.
 */

void delete_sample(sample *s);


/***********************************************************************
 * Function:  sample_init
 * Purpose:   Do initializations necessary to make sample struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling sample_init on a sample object 
 *            which is already initialized has no effect. Once the 
 *            sample object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: s - Pointer to sample struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int sample_init(sample *s, simulation *sim);


/***********************************************************************
 * Function:  sample_reset
 * Purpose:   Reset sample object to the uninitialized state and 
 *            release memory allocated by sample_init. Has no effect
 *            if sample object is already in the uninitialized state.
 *            After sample_reset has been called, sample_init can 
 *            be used to initialize sample object again, possibly with
 *            different input parameters.
 * Arguments: s - Pointer to sample struct to be reset.
 */

void sample_reset(sample *s);


/***********************************************************************
 * Function:  sample_write_log
 * Purpose:   Write input parameters of sample object to the log file.
 * Arguments: s - Pointer to sample struct.
 */

int sample_write_log(sample *s);

int background_project(sample *sam, wavefunction *wf, matrix *pm, double pos[3]);

double sample_ice_pot(void);

double sample_ice_abs_pot(double acc_en);

double sample_support_pot(void);

double sample_support_abs_pot(double acc_en);

#endif
