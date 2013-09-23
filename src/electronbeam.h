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

#ifndef ELECTRONBEAM_HEADER
#define ELECTRONBEAM_HEADER
#include "input.h"
#include "structs.h"

/***********************************************************************
 * Function:  new_electronbeam
 * Purpose:   Allocate and return a pointer to a new electronbeam struct.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to a new detector struct.
 */

electronbeam *new_electronbeam(const char *name);


/***********************************************************************
 * Function:  delete_electronbeam
 * Purpose:   Free memory allocated by electronbeam struct.
 * Arguments: d - Pointer to detector to be deleted.
 */

void delete_electronbeam(electronbeam *ed);


/***********************************************************************
 * Function:  electronbeam_init
 * Purpose:   Do initializations necessary to make electronbeam struct 
 *            ready for use in computations once parameters have been
 *            read from input file. Calling electronbeam_init on an 
 *            electronbeam object which is already initialized has no 
 *            effect. Once the electronbeam object has been initialized, 
 *            input parameters can not be changed.
 * Arguments: ed - Pointer to electronbeam struct to be initialized.
 *            sim - Pointer to simulation object. Used to access input 
 *                  of other objects if necessary.
 * Return:    0 on success, nonzero on failure.
 */

int electronbeam_init(electronbeam *ed, 
                      simulation *sim);


/***********************************************************************
 * Function:  electronbeam_reset
 * Purpose:   Reset electronbeam object to the uninitialized state and 
 *            release memory allocated by electronbeam_init. Has no effect
 *            if electronbeam object is already in the uninitialized state.
 *            After electronbeam_reset has been called, detector_init can 
 *            be used to initialize electronbeam object again, possibly with
 *            different input parameters.
 * Arguments: ed - Pointer to electronbeam struct to be reset.
 */

void electronbeam_reset(electronbeam *ed);


/***********************************************************************
 * Function:  electronbeam_write_log
 * Purpose:   Write input parameters of electronbeam object to the log file.
 * Arguments: ed - Pointer to electronbeam struct.
 */

int electronbeam_write_log(electronbeam *ed);


/***********************************************************************
 * Function:  wave_number
 * Purpose:   Compute wave number of electron wave for given acceleration
 *            energy.
 * Arguments: acc_en - Acceleration energy of electrons.
 * Return:    Wave number is standard units.
 */

double wave_number(double acc_en);


/***********************************************************************
 * Function:  potential_conv_factor
 * Purpose:   Compute factor by which projected potential should be 
 *            multiplied to obtain phase shift of electron wave.
 * Arguments: acc_en - Acceleration energy of electrons.
 * Return:    Multiplicative factor.
 */

double potential_conv_factor(double acc_en);


/***********************************************************************
 * Function:  diff_cross_sec
 * Purpose:   Compute differential cross section for elastic and inelastic
 *            scattering of electrons from a single neutral atom according 
 *            to model in Ludwig Reimer, Transmission electron microscopy: 
 *            physics of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Scattering angle.
 * Return:    Differential cross section.
 */

double diff_cross_sec(double acc_en, 
                      double Z, 
                      double theta);


/***********************************************************************
 * Function:  cross_sec
 * Purpose:   Compute total cross section for elastic and inelastic
 *            scattering of electrons from a single neutral atom according 
 *            to model in Ludwig Reimer, Transmission electron microscopy: 
 *            physics of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 * Return:    Cross section.
 */

double cross_sec(double acc_en, 
                 double Z);


/***********************************************************************
 * Function:  cross_sec_thr
 * Purpose:   Compute cross section for elastic and inelastic scattering 
 *            to angles greater than theta of electrons from a single 
 *            neutral atom according to model in Ludwig Reimer, 
 *            Transmission electron microscopy: physics of image formation 
 *            and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Smallest scattering angle.
 * Return:    Cross section
 */

double cross_sec_thr(double acc_en, 
                     double Z, 
                     double theta);


/***********************************************************************
 * Function:  diff_el_cross_sec
 * Purpose:   Compute differential cross section for elastic scattering 
 *            of electrons from a single neutral atom according to model 
 *            in Ludwig Reimer, Transmission electron microscopy: physics 
 *            of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Scattering angle.
 * Return:    Differential cross section.
 */

double diff_el_cross_sec(double acc_en, 
                         double Z, 
                         double theta);


/***********************************************************************
 * Function:  el_cross_sec
 * Purpose:   Compute total cross section for elastic scattering of 
 *            electrons from a single neutral atom according to model in 
 *            Ludwig Reimer, Transmission electron microscopy: physics 
 *            of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 * Return:    Cross section.
 */

double el_cross_sec(double acc_en, 
                    double Z);


/***********************************************************************
 * Function:  el_cross_sec_thr
 * Purpose:   Compute cross section for elastic scattering to angles 
 *            greater than theta of electrons from a single neutral atom 
 *            according to model in Ludwig Reimer, Transmission electron 
 *            microscopy: physics of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Smallest scattering angle.
 * Return:    Cross section
 */

double el_cross_sec_thr(double acc_en, 
                        double Z, 
                        double theta);


/***********************************************************************
 * Function:  diff_inel_cross_sec
 * Purpose:   Compute differential cross section for inelastic scattering 
 *            of electrons from a single neutral atom according to model 
 *            in Ludwig Reimer, Transmission electron microscopy: physics 
 *            of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Scattering angle.
 * Return:    Differential cross section.
 */

double diff_inel_cross_sec(double acc_en, 
                           double Z, 
                           double theta);


/***********************************************************************
 * Function:  inel_cross_sec
 * Purpose:   Compute total cross section for inelastic scattering of 
 *            electrons from a single neutral atom according to model in 
 *            Ludwig Reimer, Transmission electron microscopy: physics 
 *            of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 * Return:    Cross section.
 */

double inel_cross_sec(double acc_en, 
                      double Z);


/***********************************************************************
 * Function:  inel_cross_sec_thr
 * Purpose:   Compute cross section for inelastic scattering to angles 
 *            greater than theta of electrons from a single neutral atom 
 *            according to model in Ludwig Reimer, Transmission electron 
 *            microscopy: physics of image formation and microanalysis.
 * Arguments: acc_en - Acceleration energy of electrons.
 *            Z - Atomic number.
 *            theta - Smallest scattering angle.
 * Return:    Cross section
 */

double inel_cross_sec_thr(double acc_en, 
                          double Z, 
                          double alpha);


/***********************************************************************
 * Function:  electronbeam_get_acc_energy
 * Purpose:   Get acceleration energy of electrons as determined by 
 *            acceleration voltage read from input.
 * Arguments: ed - Pointer to electronbeam struct.
 * Return:    Acceleration energy of electrons.
 */

double electronbeam_get_acc_energy(electronbeam *ed);


/***********************************************************************
 * Function:  electronbeam_get_energy_spread
 * Purpose:   Get energy spread of electorns as determined by input 
 *            parameters.
 * Arguments: ed - Pointer to electronbeam struct.
 * Return:    Energy spread of electrons.
 */

double electronbeam_get_energy_spread(electronbeam *ed);

#endif
