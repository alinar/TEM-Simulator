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

#include "macros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "sample.h"
#include "array.h"
#include "electronbeam.h"
#include "functions.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "matrix.h"
#include "misc.h"
#include "particle.h"
#include "particleset.h"
#include "pdb.h"
#include "simulation.h"
#include "wavefunction.h"

/****************************************************************************/

param_table *sample_param_table(const char *name);

int sample_check_input(sample *s);

/****************************************************************************/

sample *new_sample(const char *name){
  sample  *s = malloc(sizeof(sample));
  s->param = sample_param_table(name);
  s->init = 0;
  return s;
}

/****************************************************************************/

void delete_sample(sample *s){
  sample_reset(s);
  delete_param_table(s->param);
  free(s);
}

/****************************************************************************/

param_table *sample_param_table(const char *name){
  param_table *pt = new_param_table(5, TYPE_SAMPLE, name);
  add_param_req_constr(pt, PAR_DIAMETER, "d", 0, HUGE_VAL);
  add_param_req_constr(pt, PAR_THICKNESS_CENTER, "d", 0, HUGE_VAL);
  add_param_req_constr(pt, PAR_THICKNESS_EDGE, "d", 0, HUGE_VAL);
  add_param_def(pt, PAR_OFFSET_X, "d", "0");
  add_param_def(pt, PAR_OFFSET_Y, "d", "0");
  set_comp_descr(pt, "The sample component defines the shape of the sample, \
which is a round hole in a supporting slab filled with ice. A sample component \
is required for simulation of micrographs and for generating a volume map.");
  set_param_descr(pt, PAR_DIAMETER, "The diameter of the hole measured in nm.");
  set_param_descr(pt, PAR_THICKNESS_CENTER, "The thickness of the ice at the \
center of the hole, measured in nm.");
  set_param_descr(pt, PAR_THICKNESS_EDGE, "The thickness of the ice at the \
egde of the hole, and also the thickness of the supporting slab, measured in nm.");
  set_param_descr(pt, PAR_OFFSET_X, "The x coordinate of the center of the \
hole, measured in nm.");
  set_param_descr(pt, PAR_OFFSET_Y, "The y coordinate of the center of the \
hole, measured in nm.");
  return pt;
}

/****************************************************************************/

int sample_init(sample *s, simulation *sim){
  if(s->init) return 0;
  if(sample_check_input(s)){
    WARNING("Error initializing sample: incomplete input data.\n");
    return 1;
  }
  /* Nothing to be done */
  s->init = 1;
  param_table_set_lock(s->param, 1);
  write_log_comment("Sample component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void sample_reset(sample *s){
  s->init = 0;
  param_table_set_lock(s->param, 1);
  /* Nothing to be done */
}

/****************************************************************************/

int sample_check_input(sample *s){
  if(check_params(s->param)) return 1;
  return 0;
}

/****************************************************************************/

int sample_write_log(sample *s){
  return write_parameters_to_log(s->param);
}

/****************************************************************************/

double find_root(double a, double b, double c){
  /* Finds the root closest to 0 of the equation a*x^2 + 2*b*x + c = 0.
   * Assumes that the equation has real roots.
   * The coefficient a can be 0, but b must be nonzero. */
  double b2 = b*b;
  if(fabs(a*c) < 1e-2*b2){
    return -0.5*(1 + 0.25*a*c/(b2))*c/b;
  }
  else {
    return (b/a)*(sqrt(1-a*c/b2) - 1);
  }
}

/****************************************************************************/

int background_project(sample *sam, wavefunction *wf, matrix *pm, double pos[3]){
  double k[2], m[3], x[2], ic, jc, k2, sk2, r2, a, b, c, d, s, t, z1, z2, z3, z4, zb, zc,
    path_ice, path_support, diameter;
  double *pl1, *pl2;
  long i, j, l, imax, jmax;
  const double e1 = 1e-6;
  const double e3 = 1e-3;

  if(0 == sam->init){
    WARNING("Error in backgound_project: Sample component has not been initialized.\n");
    return 1;
  }
  if(0 == wf->init){
    WARNING("Error in backgound_project: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((3 != pm->m) || (3 != pm->n)){
    WARNING("Error in background_project: Wrong size of projection matrix. Matrix must be 3 x 3.\n");
    return 1;
  }
  /* If projection direction is parallel to sample plane, do nothing.
   * This case is not relevant to simulate and would lead to infinite pathlengths. */
  if(fabs(get_matrix_entry(pm, 2, 2)) < e3){
    fill_array(&(wf->pathlength_ice), 0);
    fill_array(&(wf->pathlength_supp), 0);
    WARNING("Ignoring background in projection parallel to sample plane.\n");
    return 0;
  }
  for(j = 0; j < 2; j++){
    k[j] = get_matrix_entry(pm, 2, j)/get_matrix_entry(pm, 2, 2);
  }
  k2 = k[0]*k[0] + k[1]*k[1];
  sk2 = sqrt(1 + k2);
  diameter = get_param_double(sam->param, PAR_DIAMETER);
  r2 = 0.25 * diameter * diameter;
  zb = 0.5 * get_param_double(sam->param, PAR_THICKNESS_EDGE);
  zc = 0.5 * get_param_double(sam->param, PAR_THICKNESS_CENTER);
  a = 2*(zb - zc)/r2;
  imax = wf->phase.values.size[1];
  jmax = wf->phase.values.size[2];
  ic = 0.5*(imax-1);
  jc = 0.5*(jmax-1);
  pl1 = wf->pathlength_ice.data;
  pl2 = wf->pathlength_supp.data;
  for(j = 0; j < jmax; j++){
    for(i = 0; i < imax; i++){
      x[0] = (i-ic)*wf->phase.basis[0] + (j-jc)*wf->phase.basis[2] + wf->phase.offset[0] - pos[0];
      x[1] = (i-ic)*wf->phase.basis[1] + (j-jc)*wf->phase.basis[3] + wf->phase.offset[1] - pos[1];
      for(l = 0; l < 3; l++){
	m[l] = x[0] * get_matrix_entry(pm, 0, l) + x[1] * get_matrix_entry(pm, 1, l);
      }
      m[0] -= k[0]*m[2];
      m[1] -= k[1]*m[2];
      s = k[0]*m[0] + k[1]*m[1];
      t = m[0]*m[0] + m[1]*m[1];
      if(k2 < e1){
	/* Treat rays as perpendicular to sample */
	if(t < r2){
	  path_support = 0;
	  path_ice = 2*zc + a*t;
	}
	else {
	  path_support = 2*zb;
	  path_ice = 0;
	}
      }
      else {
	/* Ray is not perpendicular to sample */
	/*First find intersections of ray with boundary between ice and support */
	b = s/k2;
	c = (t - r2)/k2;
	d = b*b - c;
	if(d > 0){
	  z1 = -b - sqrt(d);
	  z2 = -b + sqrt(d);
	}
	else {
	  z1 = z2 = zb + 1; /* Arbitrary interval which does not intersect [-zb, zb] */
	}
	if((z2 <= -zb)||(z1 >= zb)){
	  /* Ray does not intersect the ice */
	  path_support = 2 * sk2 * zb;
	  path_ice = 0;
	}
	else {
	  /* Ray intersects the ice */
	  path_support = 0;
	  if(z1 < -zb){
	    /* Ray enters through surface of ice. Find intersection point */
	    z3 = find_root(a*k2, a*s + 1, a*t + 2*zc);
	  }
	  else {
	    z3 = z1;
	    path_support += sk2 * (z1 + zb);
	  }
	  if(z2 > zb){
	    /* Ray exits through surface of ice. Find intersection point */
	    z4 = find_root(a*k2, a*s - 1, a*t + 2*zc);
	  }
	  else {
	    z4 = z2;
	    path_support += sk2 * (zb - z2);
	  }
	  /* In some cases, for high tilt angles which should not occur in practice, 
	   * the computation of path_ice is not correct. However, this is only for rays 
	   * passing through the support which would not be used for reconstruction in 
	   * any case. */
	  path_ice = sk2 * (z4 - z3);
	}
      }
      *pl1 = path_ice;
      *pl2 = path_support;
      pl1++; pl2++;
    } /* End for loop */
  } /* End for loop */
  return 0;
}

/****************************************************************************/

double sample_ice_pot(void){
  double a[NUM_SCATTERING_PARAMETERS], b[NUM_SCATTERING_PARAMETERS], s = 0;
  int i;
  get_scattering_parameters(a, b, 8);
  for(i = 0; i < NUM_SCATTERING_PARAMETERS; ++i){
    s += a[i];
  }
  get_scattering_parameters(a, b, 1);
  for(i = 0; i < NUM_SCATTERING_PARAMETERS; ++i){
    s += 2*a[i];
  }
  s *= ONE_ANGSTROM; /* Correct for units of scattering parameters */
  s *= PLANCKS_CONSTANT_SQ/(2*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE)*ICE_DENS;
  return s;
}

/****************************************************************************/

double sample_ice_abs_pot(double acc_en){
  return PLANCKS_CONSTANT_SQ/(8*M_PI*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE*POTENTIAL_UNIT)*wave_number(acc_en)
         /(1+acc_en/ELEC_REST_ENERGY)*ICE_DENS*(cross_sec(acc_en,8) + 2*cross_sec(acc_en,1));
}

/****************************************************************************/

double sample_support_pot(void){
  double a[NUM_SCATTERING_PARAMETERS], b[NUM_SCATTERING_PARAMETERS], s = 0;
  int i;
  get_scattering_parameters(a, b, 6);
  for(i = 0; i < NUM_SCATTERING_PARAMETERS; ++i){
    s += a[i];
  }
  s *= ONE_ANGSTROM; /* Correct for units of scattering parameters */
  s *= PLANCKS_CONSTANT_SQ/(2*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE)*CARBON_DENS;
  return s;
}

/****************************************************************************/

double sample_support_abs_pot(double acc_en){
  return PLANCKS_CONSTANT_SQ/(8*M_PI*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE*POTENTIAL_UNIT)*wave_number(acc_en)
         /(1+acc_en/ELEC_REST_ENERGY)*CARBON_DENS*cross_sec(acc_en,6);
}
