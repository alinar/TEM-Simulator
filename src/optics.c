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
#include "optics.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "macros.h"
#include "matrix.h"
#include "misc.h"
#include "random.h"
#include "simulation.h"

/****************************************************************************/

param_table *optics_param_table(const char *name); 

int optics_check_input(optics *o);

/*
 * optics_set_random_defocus
 * Sets random values of the defocus for each tilt.
 * Assumes that the matrix storing defocus values has been initialized
 * and that defocus_nominal, defocus_syst_error, and defocus_nonsyst_error
 * have been initialized.
 */
int optics_generate_random_defocus(optics *o, long ntilts);

/*
 * read_defocus_from_file
 * Read defocus values from file.
 * Assumes that the matrix storing defocus values has been initialized
 * and that the file has at least as many entries as the matrix
 */
int optics_read_defocus_from_file(optics *o, long ntilts);

/*
 * write_defocus_to_file
 * Write defocus values to file.
 * Assumes that the matrix of defocus values has been initialized
 */
int optics_write_defocus_to_file(optics *o);

/****************************************************************************/

optics *new_optics(const char *name){
  optics *o = malloc(sizeof(optics));
  o->param = optics_param_table(name);
  o->defocus.data = NULL;
  o->init = 0;
  return o;
}

/****************************************************************************/

void delete_optics(optics *o){
  optics_reset(o);
  delete_param_table(o->param);
  free(o);
}

/****************************************************************************/

param_table *optics_param_table(const char *name){
  param_table *pt = new_param_table(12, TYPE_OPTICS, name);
  add_param_req_constr(pt, PAR_MAGNIFICATION, "d", 1, HUGE_VAL);
  add_param_req(pt, PAR_GEN_RAND_DEFOCUS, "b");
  add_param_opt(pt, PAR_DEFOCUS_NOMINAL, "d");
  add_param_def_constr(pt, PAR_DEFOCUS_SYST_ERROR, "d", "0", 0, HUGE_VAL);
  add_param_def_constr(pt, PAR_DEFOCUS_NONSYST_ERROR, "d", "0", 0, HUGE_VAL);
  add_param_opt(pt, PAR_DEFOCUS_FILE_IN, "s");
  add_param_opt(pt, PAR_DEFOCUS_FILE_OUT, "s");
  add_param_req(pt, PAR_CS, "d");
  add_param_req(pt, PAR_CC, "d");
  add_param_req_constr(pt, PAR_APERTURE, "d", 0, HUGE_VAL);
  add_param_req_constr(pt, PAR_FOCAL_LENGTH, "d", 0.01, HUGE_VAL);
  add_param_req_constr(pt, PAR_COND_AP_ANGLE, "d", 0, HUGE_VAL);
  set_comp_descr(pt, "The optics component specifies the characteristics of \
the optical system in the microscope. Most of the optical parameters stay the \
same in all the images of a tilt series, but the defocus can be varied from \
image to image. An optics component is required for simulation of micrographs.");
  set_param_descr(pt, PAR_MAGNIFICATION, "The magnification of the microscope.");
  set_param_descr(pt, PAR_GEN_RAND_DEFOCUS, "Controls if the defocus in each \
image should be automatically generated as random variations around a nominal \
value, or if an arbitrary sequence of defocus values should be read from a file.");
  set_param_descr(pt, PAR_DEFOCUS_NOMINAL, "The nominal defocus value used if \
gen_defocus is set to yes. Defocus is measured in micrometer and a positive \
value is intrepreted as underfocus.");
  set_param_descr(pt, PAR_DEFOCUS_SYST_ERROR, "Standard deviation of a random \
error applied to the defocus values if gen_defocus is set to yes. The same \
error is used for all images in the tilt series. Measured in micrometer.");
  set_param_descr(pt, PAR_DEFOCUS_SYST_ERROR, "Standard deviation of a random \
error applied to the defocus values if gen_defocus is set to yes. A new random \
error is computed for each image. Measured in micrometer.");
  set_param_descr(pt, PAR_DEFOCUS_FILE_IN, "A text file which the defocus of \
each image in the tilt series is read from if gen_defocus is set to no.");
  set_param_descr(pt, PAR_DEFOCUS_FILE_OUT, "A text file which the defocus of \
each image in the tilt series is written to if gen_defocus is set to yes. If \
the parameter is not specified, defocus values are not written to file.");
  set_param_descr(pt, PAR_CS, "The spherical aberration measured in mm.");
  set_param_descr(pt, PAR_CC, "The chromatic aberration measured in mm.");
  set_param_descr(pt, PAR_APERTURE, "The diameter of the aperture in the focal \
plane, measured in micrometer.");
  set_param_descr(pt, PAR_FOCAL_LENGTH, "The focal length of the primary lens \
measured in mm.");
  set_param_descr(pt, PAR_COND_AP_ANGLE, "The aperture angle of the condenser \
lens, measured in milliradian.");
  return pt;
}

/****************************************************************************/

int optics_init(optics *o, simulation *sim){
  long ntilts;
  geometry *g;
  if(o->init) return 0;
  if(optics_check_input(o)){
    WARNING("Error initializing optics: incomplete input data.\n");
    return 1;
  }
  g = get_geometry(sim, "");
  if(g == NULL){
    WARNING("Geometry required for initialization of optics.\n");
    return 1;
  }
  if(geometry_init(g, sim)) return 1;
  ntilts = get_param_long(g->param, PAR_NTILTS);
  o->init = 1;
  if(get_param_boolean(o->param, PAR_GEN_RAND_DEFOCUS)){
    if(optics_generate_random_defocus(o, ntilts)){
      optics_reset(o);
      return 1;
    }
    if(param_isset(o->param, PAR_DEFOCUS_FILE_OUT)){
      optics_write_defocus_to_file(o);
    }
  }
  else {
    if(optics_read_defocus_from_file(o, ntilts)){
      optics_reset(o);
      return 1;
    }
  }
  param_table_set_lock(o->param, 1);
  write_log_comment("Optics component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void optics_reset(optics *o){
  if(0 == o->init) return;
  free_matrix(&(o->defocus));
  param_table_set_lock(o->param, 0);
  o->init = 0;
}

/****************************************************************************/

int optics_check_input(optics *o){
  if(check_params(o->param)) return 1;
  if(get_param_boolean(o->param, PAR_GEN_RAND_DEFOCUS)){
    if(require_param(o->param, PAR_DEFOCUS_NOMINAL) || require_param(o->param, PAR_DEFOCUS_SYST_ERROR)
       || require_param(o->param, PAR_DEFOCUS_NONSYST_ERROR)) return 1;
  }
  else {
    if(require_param(o->param, PAR_DEFOCUS_FILE_IN)) return 1;
  }
  return 0;
}

/****************************************************************************/

int optics_write_log(optics *o){
  return write_parameters_to_log(o->param);
}

/****************************************************************************/

int optics_generate_random_defocus(optics *o, long ntilts){
  long i;
  double def_avg, def_var;
  if(0 == o->init){
    WARNING("Error generating random defocus: Optics component has not been initialized.\n");
    return 1;
  }
  if(ntilts <= 0){
    WARNING("Error generating random defocus: ntilts must be positive.\n");
    return 1;
  }
  init_matrix(&(o->defocus), ntilts, 1);
  def_avg = rand_gauss(get_param_double(o->param, PAR_DEFOCUS_NOMINAL), get_param_double(o->param, PAR_DEFOCUS_SYST_ERROR));
  def_var = get_param_double(o->param, PAR_DEFOCUS_NONSYST_ERROR);
  for(i = 0; i < o->defocus.m; i++){
    set_matrix_entry(&(o->defocus), i, 0, rand_gauss(def_avg, def_var));
  }
  return 0;
}

/****************************************************************************/

int optics_read_defocus_from_file(optics *o, long ntilts){
  if(0 == o->init){
    WARNING("Error reading defocus data: Optics component has not been initialized.\n");
    return 1;
  }
  if(ntilts <= 0){
    WARNING("Error reading defocus data: ntilts must be positive.\n");
    return 1;
  }
  if(read_matrix_text(&(o->defocus), get_param_string(o->param, PAR_DEFOCUS_FILE_IN))) return 1;
  if(shrink_matrix(&(o->defocus), ntilts, 1)){
    WARNING("Too few defocus data found in file %s to match ntilts = %i.\n", 
	    get_param_string(o->param, PAR_DEFOCUS_FILE_IN), (int)ntilts);
    return 1;
  }
  return 0;
}

/****************************************************************************/

int optics_write_defocus_to_file(optics *o){
  if(0 == o->init){
    WARNING("Error writing defocus file: Optics component has not been initialized.\n");
    return 1;
  }
  return write_matrix_text(&(o->defocus), get_param_string(o->param, PAR_DEFOCUS_FILE_OUT), 1, "defocus");
}

/****************************************************************************/

double optics_get_aperture(optics *o){
  return APERTURE_UNIT * get_param_double(o->param, PAR_APERTURE);
}

/****************************************************************************/

double optics_get_focal_length(optics *o){
  return FOCAL_LENGTH_UNIT * get_param_double(o->param, PAR_FOCAL_LENGTH);
}

/****************************************************************************/

double optics_get_cond_ap_angle(optics *o){
  return COND_AP_ANGLE_UNIT * get_param_double(o->param, PAR_COND_AP_ANGLE);
}

/****************************************************************************/

double optics_get_cs(optics *o){
  return CS_UNIT * get_param_double(o->param, PAR_CS);
}

/****************************************************************************/

double optics_get_cc(optics *o){
  return CC_UNIT * get_param_double(o->param, PAR_CC);
}

/****************************************************************************/

double optics_get_defocus(optics *o, long tilt){
  return DEFOCUS_UNIT * get_matrix_entry(&(o->defocus), tilt, 0);
}
