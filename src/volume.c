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
#include <string.h>
#include "volume.h"
#include "array.h"
#include "electronbeam.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "misc.h"
#include "mrc.h"
#include "particle.h"
#include "particleset.h"
#include "sample.h"
#include "simulation.h"

/****************************************************************************/

param_table *volume_param_table(const char *name);

int volume_check_input(volume *v);

int volume_write_map(volume *v, array *a, const char *fn);

int volume_add_particle(volume *v, particle *p, matrix *rm, double pos[3]);

int volume_add_background(volume *v, simulation *sim);

/****************************************************************************/

volume *new_volume(const char *name){
  volume *v = malloc(sizeof(volume));
  v->param = volume_param_table(name);
  v->pot_re.data = NULL;
  v->pot_im.data = NULL;
  v->init = 0;
  return v;
}

/****************************************************************************/

void delete_volume(volume *v){
  volume_reset(v);
  delete_param_table(v->param);
  free(v);
}

/****************************************************************************/

param_table *volume_param_table(const char *name){
  param_table *pt = new_param_table(12, TYPE_VOLUME, name);
  add_param_req_constr(pt, PAR_VOXEL_SIZE, "d", 0.01, 10);
  add_param_req_constr(pt, PAR_NX, "l", 1, 1e3);
  add_param_req_constr(pt, PAR_NY, "l", 1, 1e3);
  add_param_req_constr(pt, PAR_NZ, "l", 1, 1e3);
  add_param_def(pt, PAR_OFFSET_X, "d", "0");
  add_param_def(pt, PAR_OFFSET_Y, "d", "0");
  add_param_def(pt, PAR_OFFSET_Z, "d", "0");
  add_param_def(pt, PAR_MAP_FILE_FORMAT, "s," PAR_FILE_FORMAT__MRC "," PAR_FILE_FORMAT__MRC_INT "," PAR_FILE_FORMAT__RAW,
                PAR_FILE_FORMAT__MRC);
  add_param_def(pt, PAR_MAP_AXIS_ORDER, "s,xyz,xzy,yxz,yzx,zxy,zyx", "xyz");
  add_param_def(pt, PAR_MAP_FILE_BYTE_ORDER, "s," PAR_BYTE_ORDER__BE "," PAR_BYTE_ORDER__LE "," PAR_BYTE_ORDER__NATIVE, 
                PAR_BYTE_ORDER__NATIVE);
  add_param_req(pt, PAR_MAP_FILE_RE_OUT, "s");
  add_param_opt(pt, PAR_MAP_FILE_IM_OUT, "s");
  set_comp_descr(pt, "The volume component plays no role in the actual \
simulation of micrographs. It is used to export to a file a map of the entire \
sample or a subvolume. This map can be used as a reference, for example when \
evaluating reconstructions made from the tilt series.");
  set_param_descr(pt, PAR_VOXEL_SIZE, "The size of voxels in the volume map, \
measured in nm.");
  set_param_descr(pt, PAR_NX, "The number of voxels along the x axis of the \
volume map.");
  set_param_descr(pt, PAR_NY, "The number of voxels along the y axis of the \
volume map.");
  set_param_descr(pt, PAR_NZ, "The number of voxels along the z axis of the \
volume map.");
  set_param_descr(pt, PAR_OFFSET_X, "The x coordinate of the center of the \
volume map, measured in nm.");
  set_param_descr(pt, PAR_OFFSET_Y, "The y coordinate of the center of the \
volume map, measured in nm.");
  set_param_descr(pt, PAR_OFFSET_Z, "The z coordinate of the center of the \
volume map, measured in nm.");
  set_param_descr(pt, PAR_MAP_FILE_FORMAT, "The format of exported map files. \
\"" PAR_FILE_FORMAT__MRC "\" means MRC file with 4-byte floating point data, \"" PAR_FILE_FORMAT__MRC_INT  "\" means MRC \
file with 2-byte integer data. \"" PAR_FILE_FORMAT__RAW "\" is the same format as \"" PAR_FILE_FORMAT__MRC "\" but \
without the file header.");
  set_param_descr(pt, PAR_MAP_AXIS_ORDER, "Controls the order in which voxel \
values appear in exported map files. The first letter is the fastest varying \
index, and the last letter is the slowest varying. \"xyz\" means that the \
simulator's use of x, y, and z coordinates is consistent with the MRC file \
convention.");
  set_param_descr(pt, PAR_MAP_FILE_BYTE_ORDER, "Controls if exported map \
files are in little endian or big endian format.");
  set_param_descr(pt, PAR_MAP_FILE_RE_OUT, "The name of a file to which the \
real part of the volume map is written.");
  set_param_descr(pt, PAR_MAP_FILE_IM_OUT, "The name of a file to which the \
imaginary part of the volume map is written. If not defined, the imaginary \
part is not written to file.");
  return pt;
}

/****************************************************************************/

int volume_init(volume *v, simulation *sim){
  long nx, ny, nz;
  if(v->init) return 0;
  if(volume_check_input(v)){
    WARNING("Error initializing volume: Incomplete input data.\n");
    return 1;
  }
  nx = get_param_long(v->param, PAR_NX);
  ny = get_param_long(v->param, PAR_NY);
  nz = get_param_long(v->param, PAR_NZ);
  if((nx <= 0) || (ny <= 0) || (nz <= 0)){
    WARNING("Error initializing volume: Dimensions nx, ny and nz must all be positive.\n");
    return 1;
  }
  init_array(&(v->pot_re), nx, ny, nz);
  if(param_isset(v->param, PAR_MAP_FILE_IM_OUT)){
    init_array(&(v->pot_im), nx, ny, nz);
  }
  v->init = 1;
  param_table_set_lock(v->param, 1);
  write_log_comment("Volume component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void volume_reset(volume *v){
  if(0 == v->init) return;
  free_array(&(v->pot_re));
  free_array(&(v->pot_im));
  param_table_set_lock(v->param, 0);
  v->init = 0;
}

/****************************************************************************/

int volume_check_input(volume *v){
  return check_params(v->param);
}

/****************************************************************************/

int volume_write_log(volume *v){
  return write_parameters_to_log(v->param);
}

/****************************************************************************/

int volume_write_maps(volume *v){
  int ret = 0;
  if(0 == v->init){
    WARNING("Error writing volume maps: Volume component has not been initialized.\n");
    return 1;
  }
  if(volume_write_map(v, &(v->pot_re), get_param_string(v->param, PAR_MAP_FILE_RE_OUT))) ret = 1;
  if(param_isset(v->param, PAR_MAP_FILE_IM_OUT)){
    if(volume_write_map(v, &(v->pot_im), get_param_string(v->param, PAR_MAP_FILE_IM_OUT))) ret = 1;
  }
  return ret;
}

/****************************************************************************/

int volume_write_map(volume *v, array *a, const char *fn){
  int ret, rev;
  double voxel_size = get_param_double(v->param, PAR_VOXEL_SIZE);
  const char *map_axis_order = get_param_string(v->param, PAR_MAP_AXIS_ORDER);
  const char *format = get_param_string(v->param, PAR_MAP_FILE_FORMAT);
  rev = reverse_byte_order(get_param_string(v->param, PAR_MAP_FILE_BYTE_ORDER));
  if(0 == strcmp(format, PAR_FILE_FORMAT__RAW)){
    ret = write_array_float4b(a, fn, no_header, map_axis_order, rev, voxel_size); 
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC)){
    ret = write_array_float4b(a, fn, mrc_header, map_axis_order, rev, voxel_size);
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT)){
    ret = write_array_int2b(a, fn, mrc_header, map_axis_order, rev, voxel_size, POTENTIAL_UNIT/INT_FILE_POT_UNIT);
  }
  else {
    WARNING("Unknown file format %s.\n", format);
    return 1;
  }
  if(ret){
    WARNING("Error writing data to file %s.\n", fn);
    return 1;
  }
  write_log_comment("Volume map written to file %s.\n\n", fn);
  return 0;
}

/****************************************************************************/

int volume_add_particle(volume *v, particle *p, matrix *rm, double pos[3]){
  long imin[3], imax[3], i, j, k, l, m, ip, jp, kp, d[3], c, n1, n1n2;
  double x[3], y[3], ic[3], d0[3], famp, fph, r, re, im, w[3][2], wp, voxel_size_v, voxel_size_p;
  double *vdr, *vdi, *wd, *pdr, *pdi;
  boolean use_imag_pot_v;
  array pot_re, pot_im, weight;
  if(0 == v->init){
    WARNING("Error in volume_add_particle: Volume component has not been initialized.\n");
    return 1;
  }
  if(0 == p->init){
    WARNING("Error in volume_add_particle: Particle component has not been initialized.\n");
    return 1;
  }
  if((3 != rm->m) || (3 != rm->n)){
    WARNING("Error in volume_add_particle: Wrong size of rotation matrix. Matrix must be 3 x 3.\n");
    return 0;
  }
  ic[0] = 0.5*(p->pot_re.size[0] - 1);
  ic[1] = 0.5*(p->pot_re.size[1] - 1);
  ic[2] = 0.5*(p->pot_re.size[2] - 1);
  x[0] = -get_param_double(v->param, PAR_OFFSET_X);
  x[1] = -get_param_double(v->param, PAR_OFFSET_Y);
  x[2] = -get_param_double(v->param, PAR_OFFSET_Z);
  voxel_size_v = get_param_double(v->param, PAR_VOXEL_SIZE);
  voxel_size_p = get_param_double(p->param, PAR_VOXEL_SIZE);
  use_imag_pot_v = param_isset(v->param, PAR_MAP_FILE_IM_OUT);
  for(l = 0; l < 3; l++){
    d0[l] = x[l]/voxel_size_v;
    x[l] += pos[l];
    r = 0;
    for(m = 0; m < 3; m++){
      r += ic[m] * fabs(get_matrix_entry(rm, l, m));
    }
    y[l] = x[l] + r;
    x[l] -= r;
  }
  for(l = 0; l < 3; l++){
    imin[l] = max_l(-1, (long)floor(x[l]/voxel_size_v + 0.5*(v->pot_re.size[l] - 1)));
    imax[l] = min_l(1 + v->pot_re.size[l], (long)ceil(y[l]/voxel_size_v + 0.5*(v->pot_re.size[l] - 1)));
  }
  if((imin[0] >= imax[0])||(imin[1] >= imax[1])||(imin[2] >= imax[2])){
    return 0; /* Intersection between particle and volume is empty */
  }
  for(l = 0; l < 3; l++){
    d0[l] += 0.5*(v->pot_re.size[l] - 1) - imin[l];
  }
  init_array(&pot_re, imax[0]-imin[0], imax[1]-imin[1], imax[2]-imin[2]);
  init_array(&weight, imax[0]-imin[0], imax[1]-imin[1], imax[2]-imin[2]);
  fill_array(&pot_re, 0);
  fill_array(&weight, 0);
  if(use_imag_pot_v){
    init_array(&pot_im, imax[0]-imin[0], imax[1]-imin[1], imax[2]-imin[2]);
    fill_array(&pot_im, 0);
  }
  pdr = p->pot_re.data;
  if(p->use_imag_pot){
    pdi = p->pot_im.data;
    famp = 1;
    fph = 1;
  }
  else {
    pdi = p->pot_re.data;
    famp = get_param_double(p->param, PAR_FAMP);
    fph = sqrt(1 - famp * famp);
  }
  n1 = pot_re.size[0];
  n1n2 = n1 * pot_re.size[1];
  for(kp = 0; kp < p->pot_re.size[2]; kp++){
    y[2] = voxel_size_p * (kp - ic[2]);
    for(jp = 0; jp < p->pot_re.size[1]; jp++){
      y[1] = voxel_size_p * (jp - ic[1]);
      for(ip = 0; ip < p->pot_re.size[0]; ip++){
	y[0] = voxel_size_p * (ip - ic[0]);
	for(l = 0; l < 3; l++){
	  x[l] = pos[l];
	  for(m = 0; m < 3; m++){
	    x[l] += y[m] * get_matrix_entry(rm, l, m);
	  }
	  x[l] /= voxel_size_v;
	  x[l] += d0[l];
	  d[l] = (long)floor(x[l]);
	  w[l][1] = x[l] - d[l];
	  w[l][0] = 1 - w[l][1];
	}
	if(array_index_in_range(&pot_re, d[0], d[1], d[2]) && array_index_in_range(&pot_re, d[0]+1, d[1]+1, d[2]+1)){
	  if(use_imag_pot_v){
	    vdr = get_array_entry_ptr(&pot_re, d[0], d[1], d[2]);
	    vdi = get_array_entry_ptr(&pot_im, d[0], d[1], d[2]);
	    wd = get_array_entry_ptr(&weight, d[0], d[1], d[2]);
	    re = fph * (*pdr);
	    im = famp * (*pdi);
	    for(k = 0; k <= 1; k++){
	      for(j = 0; j <= 1; j++){
		for(i = 0; i <= 1; i++){
		  wp = w[0][i]*w[1][j]*w[2][k];
		  c = i + n1*j + n1n2*k;
		  vdr[c] += wp*re;
		  vdi[c] += wp*im;
		  wd[c] += wp;
		}
	      }
	    }
	  }
	  else {
	    vdr = get_array_entry_ptr(&pot_re, d[0], d[1], d[2]);
	    wd = get_array_entry_ptr(&weight, d[0], d[1], d[2]);
	    re = fph * (*pdr);
	    for(k = 0; k <= 1; k++){
	      for(j = 0; j <= 1; j++){
		for(i = 0; i <= 1; i++){
		  wp = w[0][i]*w[1][j]*w[2][k];
		  c = i + n1*j + n1n2*k;
		  vdr[c] += wp*re;
		  wd[c] += wp;
		}
	      }
	    }
	  }
	}
	pdr++;
	pdi++;
      }
    }
  }
  wd = weight.data;
  m = nele_array(&weight);
  for(i = 0; i < m; i++){
    if(*wd > 0){
      *wd = 1/(*wd);
    }
    wd++;
  }
  mult_array(&weight, &pot_re);
  add_array_offset(&pot_re, &(v->pot_re), imin[0], imin[1], imin[2]);
  free_array(&pot_re);

  if(use_imag_pot_v){
    mult_array(&weight, &pot_im);
    add_array_offset(&pot_im, &(v->pot_im), imin[0], imin[1], imin[2]);
    free_array(&pot_im);
  }
  free_array(&weight);
  return 0;
}

/****************************************************************************/

int volume_add_background(volume *v, simulation *sim){
  double a, b, c, r, r2, rad, sm, u, w, x, y, z, zb, voxel_size, offset[3], 
    support_pot, support_abs = 0, ice_pot, ice_abs = 0, acc_en;
  long i, j, k;
  boolean use_imag_pot;
  sample *s = get_sample(sim, "");
  electronbeam *ed;
  if(volume_init(v, sim)){
    WARNING("Error in volume_add_background: Volume component could not be initialized.\n");
    return 1;
  }
  if(s == NULL || sample_init(s, sim)){
    WARNING("Error in volume_add_background: Sample component could not be initialized.\n");
    return 1;
  }
  support_pot = sample_support_pot();
  ice_pot = sample_ice_pot();
  use_imag_pot = param_isset(v->param, PAR_MAP_FILE_IM_OUT);
  if(use_imag_pot){
    ed = get_electronbeam(sim, "");
    if(ed == NULL || require_param(ed->param, PAR_ACC_VOLTAGE)){
      WARNING("Error in volume_add_background: Electronbeam component could not be initialized.\n");
      return 1;
    }
    acc_en = electronbeam_get_acc_energy(ed);
    support_abs = sample_support_abs_pot(acc_en);
    ice_abs = sample_ice_abs_pot(acc_en);
  }
  voxel_size = get_param_double(v->param, PAR_VOXEL_SIZE);
  offset[0] = get_param_double(v->param, PAR_OFFSET_X) - get_param_double(s->param, PAR_OFFSET_X)
    - voxel_size * 0.5*(v->pot_re.size[0] - 1);
  offset[1] = get_param_double(v->param, PAR_OFFSET_Y) - get_param_double(s->param, PAR_OFFSET_Y)
    - voxel_size * 0.5*(v->pot_re.size[1] - 1);
  offset[2] = get_param_double(v->param, PAR_OFFSET_Z) - voxel_size * 0.5*(v->pot_re.size[2] - 1);
  rad = 0.5 * get_param_double(s->param, PAR_DIAMETER);
  b = 0.5 * get_param_double(s->param, PAR_THICKNESS_CENTER);
  c = 0.5 * get_param_double(s->param, PAR_THICKNESS_EDGE);
  a = (c-b)/(rad*rad);
  sm = voxel_size;  /* Amount of smoothing at boundaries */
  for(k = 0; k < v->pot_re.size[2]; k++){
    z = fabs(offset[2] + k * voxel_size);
    for(j = 0; j < v->pot_re.size[1]; j++){
      y = offset[1] + j * voxel_size;
      for(i = 0; i < v->pot_re.size[0]; i++){
	x = offset[0] + i * voxel_size;
	r2 = x*x + y*y;
	r = sqrt(r2);
	zb = (r < rad)?a*r2+b:c;
	if(r < rad-sm){
	  u = 1;
	}
	else if(r < rad+sm){
	  u = 0.5*(rad + sm - r)/sm;
	}
	else {
	  u = 0;
	}
	if(z < zb-sm){
	  w = 1;
	}
	else if(z < zb+sm){
	  w = 0.5*(zb + sm - z)/sm;
	}
	else {
	  w = 0;
	}
	add_to_array_entry(&(v->pot_re), i, j, k, w*(u*ice_pot + (1-u)*support_pot));
	if(use_imag_pot){
	  add_to_array_entry(&(v->pot_im), i, j, k, w*(u*ice_abs + (1-u)*support_abs));
	}
      }
    }
  }
  return 0;
}

/****************************************************************************/

int volume_get_potential(volume *v, simulation *sim){
  int j;
  long i;
  double pos[3];
  matrix rm;
  particle *p;
  particleset *ps;
  if(volume_init(v, sim)){
    WARNING("Error in volume_get_potential: Volume component could not be initialized.\n");
    return 1;
  }
  for(j = 0; j < sim->num_particleset; j++){
    if(particleset_init(sim->particleset[j], sim)){
      WARNING("Error in volume_get_potential: Particleset component number %i could not be initialized.\n", j+1);
      return 1;
    }
  }
  init_matrix(&rm, 3, 3);
  fill_array(&(v->pot_re), 1e-10); /* Some software seems to have problems with MRC files containing 0 */
  if(param_isset(v->param, PAR_MAP_FILE_IM_OUT)){
    fill_array(&(v->pot_im), 1e-10);
  }
  for(j = 0; j < sim->num_particleset; j++){
    ps = sim->particleset[j];
    p = get_particle(sim, get_param_string(ps->param, PAR_PARTICLE_TYPE));
    if(p == NULL || particle_init(p, sim)){
      WARNING("Particle %s not found.\n", get_param_string(ps->param, PAR_PARTICLE_TYPE));
      return 1;
    }
    for(i = 0; i < ps->coordinates.m; i++){
      if(get_particle_coord(&rm, pos, ps, i)) return 1;
      if(volume_add_particle(v, p, &rm, pos)) return 1;
    }
  }
  if(volume_add_background(v, sim)) return 1;
  return 0;
}

