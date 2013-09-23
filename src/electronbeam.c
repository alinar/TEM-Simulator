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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "electronbeam.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "macros.h"
#include "matrix.h"
#include "misc.h"
#include "random.h"
#include "simulation.h"


/***********************************************************************
 * Function:  electronbeam_param_table
 * Purpose:   Initialize table of input parameters used by electronbeam.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to parameter table.
 */

param_table *electronbeam_param_table(const char *name);


/***********************************************************************
 * Function:  electronbeam_check_input
 * Purpose:   Check validity of input parameters for electronbeam object.
 * Arguments: ed - Pointer to electronbeam struct.
 * Return:    0 if input is OK, nonzero otherwise.
 */

int electronbeam_check_input(electronbeam *ed);


/***********************************************************************
 * Function:  electronbeam_gen_dose
 * Purpose:   Generate electron dose for each micrograph in tilt series
 *            from total or average dose and variation. Called during
 *            initialization.
 * Arguments: ed - Pointer to electronbeam struct.
 *            ntilts - Number of images in tilt series.
 * Return:    0 if input is OK, nonzero otherwise.
 */

int electronbeam_gen_dose(electronbeam *ed, long ntilts);


/***********************************************************************
 * Function:  electronbeam_read_dose_from_file
 * Purpose:   Read electron dose for each micrograph from a text file.
 *            Called during initialization.
 * Arguments: ed - Pointer to electronbeam struct.
 *            ntilts - Number of images in tilt series.
 * Return:    0 if input is OK, nonzero otherwise.
 */

int electronbeam_read_dose_from_file(electronbeam *ed, long ntilts);


/***********************************************************************
 * Function:  electronbeam_write_dose_to_file
 * Purpose:   Write electron dose for each micrograph to a text file.
 *            Called during initialization.
 * Arguments: ed - Pointer to electronbeam struct.
 *            ntilts - Number of images in tilt series.
 * Return:    0 if input is OK, nonzero otherwise.
 */

int electronbeam_write_dose_to_file(electronbeam *ed);


/****************************************************************************/

electronbeam *new_electronbeam(const char *name){
  electronbeam *ed = malloc(sizeof(electronbeam));
  ed->param = electronbeam_param_table(name);
  ed->dose.data = NULL;
  ed->init = 0;
  return ed;
}

/****************************************************************************/

void delete_electronbeam(electronbeam *ed){
  electronbeam_reset(ed);
  delete_param_table(ed->param);
  free(ed);
}

/****************************************************************************/

param_table *electronbeam_param_table(const char *name){
  param_table *pt = new_param_table(8, TYPE_ELECTRONBEAM, name);
  add_param_req_constr(pt, PAR_ACC_VOLTAGE, "d", 1, 1e4);
  add_param_req_constr(pt, PAR_ENERGY_SPREAD, "d", 0, HUGE_VAL);
  add_param_req(pt, PAR_GEN_DOSE, "b");
  add_param_opt_constr(pt, PAR_TOTAL_DOSE, "d", 0, HUGE_VAL);
  add_param_opt_constr(pt, PAR_DOSE_PER_IM, "d", 0, HUGE_VAL);
  add_param_def_constr(pt, PAR_DOSE_SD, "d", "0", 0, 1);
  add_param_opt(pt, PAR_DOSE_FILE_IN, "s");
  add_param_opt(pt, PAR_DOSE_FILE_OUT, "s");
  set_comp_descr(pt, "The electronbeam component controls the properties of \
the electron beam. An electronbeam component is required for simulation of \
micrographs.");
  set_param_descr(pt, PAR_ACC_VOLTAGE, "Acceleration voltage of the \
microscope in kV.");
  set_param_descr(pt, PAR_ENERGY_SPREAD, "Energy spread of the beam in eV.");
  set_param_descr(pt, PAR_GEN_DOSE, "Controls whether the electron dose \
for each image should be automatically generated or read from a file.");
  set_param_descr(pt, PAR_TOTAL_DOSE, "If the dose in each image is to be \
automatically generated, either total_dose or dose_per_im must be provided. \
total_dose sets the total dose in all the images, measured in electrons per nm^2.");
  set_param_descr(pt, PAR_DOSE_PER_IM, "If the dose in each image is to be \
automatically generated, either total_dose or dose_per_im must be provided. \
dose_per_im sets the average dose in each image, measured in electrons per nm^2.");
  set_param_descr(pt, PAR_DOSE_SD, "Standard deviation of the average dose \
used in each image, relative to the average dose per image, if the dose in \
each image is automaticlly generated.");
  set_param_descr(pt, PAR_DOSE_FILE_IN, "A text file from which to read the \
dose of each image, if it is not to be automatically generated.");
  set_param_descr(pt, PAR_DOSE_FILE_OUT, "A text file to which the dose in \
each image is written if it has been automatically generated.");
  return pt;
}

/****************************************************************************/

int electronbeam_init(electronbeam *ed, simulation *sim){
  long ntilts;
  geometry *g;
  if(ed->init) return 0;
  if(electronbeam_check_input(ed)){
    WARNING("Error initializing electronbeam: incomplete input data.\n");
    return 1;
  }
  g = get_geometry(sim, "");
  if(g == NULL){
    WARNING("Geometry required for initialization of electronbeam.\n");
    return 1;
  }
  if(geometry_init(g, sim)) return 1;
  ntilts = get_param_long(g->param, PAR_NTILTS);
  ed->init = 1;
  if(get_param_boolean(ed->param, PAR_GEN_DOSE)){
    if(electronbeam_gen_dose(ed, ntilts)){
      electronbeam_reset(ed);
      return 1;
    }
    if(param_isset(ed->param, PAR_DOSE_FILE_OUT)){
      electronbeam_write_dose_to_file(ed);
    }
  }
  else {
    if(electronbeam_read_dose_from_file(ed, ntilts)){
      electronbeam_reset(ed);
      return 1;
    }
  }
  param_table_set_lock(ed->param, 1);
  write_log_comment("Electronbeam component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void electronbeam_reset(electronbeam *ed){
  if(0 == ed->init) return;
  free_matrix(&ed->dose);
  param_table_set_lock(ed->param, 0);
  ed->init = 0;
}

/****************************************************************************/

int electronbeam_check_input(electronbeam *ed){
  int total, per_im;
  if(check_params(ed->param)) return 1;
  if(get_param_boolean(ed->param, PAR_GEN_DOSE)){
    total = param_isset(ed->param, PAR_TOTAL_DOSE);
    per_im = param_isset(ed->param, PAR_DOSE_PER_IM);
    if(!(total || per_im)){
      WARNING("Electron dose must be specified as total dose or dose per image.");
      return 1;
    }
    if(total && per_im){
      WARNING("Both total dose and dose per image specified.\nTotal dose value will be used in simulation.");
    }
  }
  else {
    if(require_param(ed->param, PAR_DOSE_FILE_IN)) return 1;
  }
  return 0;
}

/****************************************************************************/

int electronbeam_write_log(electronbeam *ed){
  return write_parameters_to_log(ed->param);
}

/****************************************************************************/

int electronbeam_gen_dose(electronbeam *ed, long ntilts){
  long i;
  double avg_dose, dose_sd;
  if(0 == ed->init){
    WARNING("Error generating random dose data: Electronbeam component has not been initialized.\n");
    return 1;
  }
  if(ntilts <= 0){
    WARNING("Error generating random dose data: ntilts must be positive.\n");
    return 1;
  }
  init_matrix(&(ed->dose), ntilts, 1);
  if(param_isset(ed->param, PAR_TOTAL_DOSE)){
    avg_dose = get_param_double(ed->param, PAR_TOTAL_DOSE) / ed->dose.m;
  }
  else {
    avg_dose = get_param_double(ed->param, PAR_DOSE_PER_IM);
  }
  dose_sd = get_param_double(ed->param, PAR_DOSE_SD);
  for(i = 0; i < ed->dose.m; i++){
    set_matrix_entry(&(ed->dose), i, 0, max_d(0, avg_dose * rand_gauss(1, dose_sd)));
  }
  return 0;
}

/****************************************************************************/

int electronbeam_read_dose_from_file(electronbeam *ed, long ntilts){
  if(0 == ed->init){
    WARNING("Error reading dose data from file: Electronbeam component has not been initialized.\n");
    return 1;
  }
  if(ntilts <= 0){
    WARNING("Error reading random dose data: ntilts must be positive.\n");
    return 1;
  }
  if(read_matrix_text(&(ed->dose), get_param_string(ed->param, PAR_DOSE_FILE_IN))) return 1;
  if(shrink_matrix(&(ed->dose), ntilts, 1)){
    WARNING("Too few dose data found in file %s to match ntilts = %i.\n", 
	    get_param_string(ed->param, PAR_DOSE_FILE_IN), (int)ntilts);
    return 1;
  }
  return 0;
}

/****************************************************************************/

int electronbeam_write_dose_to_file(electronbeam *ed){
  if(0 == ed->init){
    WARNING("Error reading dose data from file: Electronbeam component has not been initialized.\n");
    return 1;
  }
  return write_matrix_text(&(ed->dose), get_param_string(ed->param, PAR_DOSE_FILE_OUT), 1, "dose");
}

/****************************************************************************/

double wave_number(double acc_en){
  return 2*M_PI/(SPEED_OF_LIGHT*PLANCKS_CONSTANT) * sqrt(acc_en * (acc_en + 2*ELEC_REST_ENERGY));
}

/****************************************************************************/

double potential_conv_factor(double acc_en){
  return 4*M_PI*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE/PLANCKS_CONSTANT_SQ
    * POTENTIAL_UNIT * (1 + acc_en/ELEC_REST_ENERGY)/wave_number(acc_en); 
}

/****************************************************************************/

double diff_cross_sec(double acc_en, double Z, double theta){
  return diff_el_cross_sec(acc_en, Z, theta) + diff_inel_cross_sec(acc_en, Z, theta);
}

/****************************************************************************/

double cross_sec(double acc_en, double Z){
  return el_cross_sec(acc_en, Z) + inel_cross_sec(acc_en, Z);
}

/****************************************************************************/

double cross_sec_thr(double acc_en, double Z, double theta){
  return el_cross_sec_thr(acc_en, Z, theta) + inel_cross_sec_thr(acc_en, Z, theta);
}

/****************************************************************************/

double diff_el_cross_sec(double acc_en, double Z, double theta){
  double relang2, a;
  double k = wave_number(acc_en);
  double theta02 = pow(Z, 1/3.0)/(k * BOHR_RADIUS);
  theta02 *= theta02;
  relang2 = theta*theta/theta02;
  a = BOHR_RADIUS * (1 + acc_en/ELEC_REST_ENERGY) / (1 + relang2);
  return 4 * pow(Z, 2/3.0) * a * a;
}

/****************************************************************************/

double el_cross_sec(double acc_en, double Z){
  double k = wave_number(acc_en);
  double a = (1 + acc_en/ELEC_REST_ENERGY)/k;
  return 4 * M_PI * pow(Z, 4/3.0) * a * a;
}

/****************************************************************************/

double el_cross_sec_thr(double acc_en, double Z, double theta){
  double relang2, a;
  double k = wave_number(acc_en);
  double theta02 = pow(Z, 1/3.0)/(k * BOHR_RADIUS);
  theta02 *= theta02;
  relang2 = theta*theta/theta02;
  a = (1 + acc_en/ELEC_REST_ENERGY)/k;
  return 4 * M_PI * pow(Z, 4/3.0) * a * a / (1 + relang2);
}

/****************************************************************************/

double diff_inel_cross_sec(double acc_en, double Z, double theta){
  double k, thetaE2, theta02, theta2, a, b;
  k = wave_number(acc_en);
  thetaE2 = 13.5 * ONE_ELECTRONVOLT * Z / (4 * acc_en);
  thetaE2 *= thetaE2;
  theta02 = pow(Z, 1/3.0)/(k * BOHR_RADIUS);
  theta02 *= theta02;
  theta2 = theta*theta;
  a = 1 + theta2/theta02;
  b = (1 + acc_en / ELEC_REST_ENERGY) / (k * k * BOHR_RADIUS * (theta2 + thetaE2));
  return 4 * Z * b * b * (1 - 1/(a*a));
}

/****************************************************************************/

double inel_cross_sec(double acc_en, double Z){
  double k, thetaE2, theta02, a, b;
  k = wave_number(acc_en);
  thetaE2 = 13.5 * ONE_ELECTRONVOLT * Z / (4 * acc_en);
  thetaE2 *= thetaE2;
  theta02 = pow(Z, 1/3.0)/(k * BOHR_RADIUS);
  theta02 *= theta02;
  a = theta02 - thetaE2;
  b = (1 + acc_en / ELEC_REST_ENERGY) / (k * k * BOHR_RADIUS * a);
  return 4 * M_PI * Z * b * b / a
           * (-3*theta02*theta02 + 4*theta02*thetaE2 - thetaE2*thetaE2
	      + 2*theta02*theta02*(log(theta02)-log(thetaE2)));
}

/****************************************************************************/

double inel_cross_sec_thr(double acc_en, double Z, double angle){
  double relang2, a;
  double k = wave_number(acc_en);
  double theta02 = pow(Z, 1/3.0)/(k * BOHR_RADIUS);
  theta02 *= theta02;
  relang2 = angle*angle/theta02;
  a = (1 + acc_en/ELEC_REST_ENERGY)/k;
  return 4 * M_PI * pow(Z, 1/3.0) * a * a * (2*log(1 + 1/relang2) - 1/(1 + relang2));
}

/****************************************************************************/

double electronbeam_get_acc_energy(electronbeam *ed){
  return ACC_ENERGY_UNIT * get_param_double(ed->param, PAR_ACC_VOLTAGE);
}

/****************************************************************************/

double electronbeam_get_energy_spread(electronbeam *ed){
  return ENERGY_SPREAD_UNIT * get_param_double(ed->param, PAR_ENERGY_SPREAD);
}
