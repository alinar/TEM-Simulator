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

#ifndef DETECTOR_HEADER
#define DETECTOR_HEADER
#include "input.h"
#include "structs.h"

/***********************************************************************
 * Function:  new_detector
 * Purpose:   Allocate and return a pointer to a new detector struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new detector struct.
 */

detector *new_detector(const char *name);


/***********************************************************************
 * Function:  delete_detector
 * Purpose:   Free memory allocated by detector struct.
 * Arguments: d - Pointer to detector to be deleted.
 */

void delete_detector(detector *d);


/***********************************************************************
 * Function:  detector_init
 * Purpose:   Do initializations necessary to make detector struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling detector_init on a detector object 
 *            which is already initialized has no effect. Once the 
 *            detector object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: d - Pointer to detector struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int detector_init(detector *d, 
                  simulation *sim);


/***********************************************************************
 * Function:  detector_init_share
 * Purpose:   Do initializations necessary to make detector struct ready
 *            for use in computations once parameters have been from 
 *            input file. Calling detector_init on a detector object 
 *            which is already initialized has no effect. Once the 
 *            detector object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: d - Pointer to detector struct to be initialized.
 *            count, count_ft - Pointers to arrays, which the detector
 *            arrays will share memory with. Can be null pointers, in
 *            which case the memory will not be shared.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int detector_init_share(detector *d, 
                        simulation *sim, 
                        array *count, 
                        array *count_ft);


/***********************************************************************
 * Function:  detector_init_all
 * Purpose:   Initialize all uninitialized detectors in simulation, so 
 *            that they share the same memory.
 * Arguments: sim - Pointer to simulation object. 
 * Return:    0 on success, nonzero on failure.
 */

int detector_init_all(simulation *sim);


/***********************************************************************
 * Function:  detector_reset
 * Purpose:   Reset detector object to the uninitialized state and 
 *            release memory allocated by detector_init. Has no effect
 *            if detector object is already in the uninitialized state.
 *            After detector_reset has been called, detector_init can 
 *            be used to initialize detector object again, possibly with
 *            different input parameters.
 * Arguments: d - Pointer to detector struct to be reset.
 */

void detector_reset(detector *d);


/***********************************************************************
 * Function:  detector_write_log
 * Purpose:   Write input parameters of detector object to the log file.
 * Arguments: d - Pointer to detector struct.
 */

int detector_write_log(detector *d);


/***********************************************************************
 * Function:  detector_write_image
 * Purpose:   Add an output image to the image file and update file 
 *            header if necessary.
 * Arguments: d - Pointer to detector struct, assumed to be initialized.
 * Return:    0 on success, nonzero on failure.
 */

int detector_write_image(detector *d);


/***********************************************************************
 * Function:  detector_apply_quantization
 * Purpose:   Apply quantization noise to the detector image being 
 *            processed.
 * Arguments: d - Pointer to detector struct, assumed to be initialized.
 * Return:    0 on success, nonzero on failure.
 */

int detector_apply_quantization(detector *d);


/***********************************************************************
 * Function:  detector_apply_mtf
 * Purpose:   Apply modulation transfer function (detector blurring) to
 *            the detector image being processed.
 * Arguments: d - Pointer to detector struct, assumed to be initialized.
 * Return:    0 on success, nonzero on failure.
 */

int detector_apply_mtf(detector *d);


/***********************************************************************
 * Function:  detector_get_intensity
 * Purpose:   Get electron wave intensity in the detector plane from
 *            wave function object.
 * Arguments: d - Pointer to detector struct, assumed to be initialized.
 *            wf - Pointer to wave function object, assumed to be 
 *                 initialized.
 *            tilt - number of image in tilt series.
 * Return:    0 on success, nonzero on failure.
 */

int detector_get_intensity(detector *d, 
                           wavefunction *wf, 
                           long tilt);


/***********************************************************************
 * Function:  detector_get_pixel_size
 * Purpose:   Get size of detector pixels in standard units.
 * Arguments: d - Pointer to detector struct, assumed to be initialized.
 * Return:    Size of detector pixels in standard units.
 */

double detector_get_pixel_size(detector *d);

#endif
