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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "particle.h"
#include "array.h"
#include "electronbeam.h"
#include "/opt/local/include/fftw3.h"
#include "functions.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "matrix.h"
#include "misc.h"
#include "mrc.h"
#include "pdb.h"
#include "random.h"
#include "sample.h"
#include "simulation.h"
#include "wavefunction.h"
/****************************************************************************/

param_table *particle_param_table(const char *name);

int particle_check_input(particle *p);

int particle_read_maps(particle *p);

int particle_read_map(particle *p, array *a, const char *fn);

int particle_write_maps(particle *p);

int particle_write_map(particle *p, array *a, const char *fn);

int generate_random_particle(array *a, double contrast, double smoothness, boolean positive);

/****************************************************************************/

particle *new_particle(const char *name){
  particle *p = malloc(sizeof(particle));
  p->param = particle_param_table(name);
  p->pot_re.data = NULL;
  p->pot_im.data = NULL;
  p->lap_pot_re.data = NULL;
  p->lap_pot_im.data = NULL;
  p->init = 0;
  return p;
}

/****************************************************************************/

void delete_particle(particle *p){
  particle_reset(p);
  delete_param_table(p->param);
  free(p);
}

/****************************************************************************/

param_table *particle_param_table(const char *name){
  param_table *pt = new_param_table(22, TYPE_PARTICLE, name);
  add_param_opt_constr(pt, PAR_VOXEL_SIZE, "d", 0.01, 10);
  add_param_opt_constr(pt, PAR_NX, "l", 1, 1e3);
  add_param_opt_constr(pt, PAR_NY, "l", 1, 1e3);
  add_param_opt_constr(pt, PAR_NZ, "l", 1, 1e3);
  add_param_def(pt, PAR_USE_IMAG_POT, "b", YES_STRING);
  add_param_opt_constr(pt, PAR_FAMP, "d", 0, 1);
  add_param_def(pt, PAR_USE_DEFOCUS_CORR, "b", NO_STRING);
  add_param_req(pt, PAR_SOURCE, "s," PAR_SOURCE__PDB "," PAR_SOURCE__MAP "," PAR_SOURCE__RANDOM);
  add_param_opt(pt, PAR_PDB_FILE_IN, "s");
  add_param_opt(pt, PAR_PDB_TRANSF_FILE_IN, "s");
  add_param_def(pt, PAR_ADD_HYDROGEN, "b", YES_STRING);
  add_param_def(pt, PAR_MAP_FILE_FORMAT, "s," PAR_FILE_FORMAT__MRC "," PAR_FILE_FORMAT__MRC_INT "," PAR_FILE_FORMAT__RAW, 
                PAR_FILE_FORMAT__MRC);
  add_param_def(pt, PAR_MAP_AXIS_ORDER, "s,xyz,xzy,yxz,yzx,zxy,zyx", "xyz");
  add_param_def(pt, PAR_MAP_FILE_BYTE_ORDER, "s," PAR_BYTE_ORDER__BE "," PAR_BYTE_ORDER__LE "," PAR_BYTE_ORDER__NATIVE, 
                PAR_BYTE_ORDER__NATIVE);
  add_param_opt(pt, PAR_MAP_FILE_RE_IN, "s");
  add_param_opt(pt, PAR_MAP_FILE_IM_IN, "s");
  add_param_opt(pt, PAR_MAP_FILE_RE_OUT, "s");
  add_param_opt(pt, PAR_MAP_FILE_IM_OUT, "s");
  add_param_opt(pt, PAR_CONTRAST_RE, "d");
  add_param_opt(pt, PAR_CONTRAST_IM, "d");
  add_param_opt_constr(pt, PAR_SMOOTHNESS, "d", 0, 10);
  add_param_opt(pt, PAR_MAKE_POSITIVE, "b");
  set_comp_descr(pt, "A particle component defines a single particle, for example \
a macromolecule. Once a particle is defined, one or several copies of the particle \
can be placed in the sample using a particleset component. A simulation of \
micrographs can use one or several particle components, or possibly none at all \
if an empty sample is to be simulated.");
  set_param_descr(pt, PAR_VOXEL_SIZE, "The voxel size in nm of the map of the \
particle potential used in the simulations.");
  set_param_descr(pt, PAR_NX, "The number of voxels along the x-axis in the map \
of the particle potential. The size of the particle map must be specified if the \
map is to be read from a raw data file (without header) or randomly generated. It \
can optionally be specified if the map is generated from a PDB file. If the map is \
read from an MRC file the size is determined by the file header.");
  set_param_descr(pt, PAR_NY, "The number of voxels along the y-axis in the map of \
the particle potential. The size of the particle map must be specified if the map \
is to be read from a raw data file (without header) or randomly generated. It can \
optionally be specified if the map is generated from a PDB file. If the map is read \
from an MRC file the size is determined by the file header.");
  set_param_descr(pt, PAR_NZ, "The number of voxels along the z-axis in the map of \
the particle potential. The size of the particle map must be specified if the map \
is to be read from a raw data file (without header) or randomly generated. It can \
optionally be specified if the map is generated from a PDB file. If the map is read \
from an MRC file the size is determined by the file header.");
  set_param_descr(pt, PAR_USE_IMAG_POT, "Controls whether a separate map of an \
\"imaginary\" absorption potential should be used in the simulation along with the \
real valued electrostatic potential. If a separate map is not used, the imaginary \
potential is taken to be a multiple of the real potential, given by the parameter \
famp. When using a PDB file as source, " PAR_USE_IMAG_POT " should be set to yes in order \
to make the simulation as correct as possible.");
  set_param_descr(pt, PAR_FAMP, "If use_imag_pot is set to no, the same map is used \
to define both the real and the imaginary potential. The imaginary potential is \
taken as the map multiplied by famp, and the real potential as the map multiplied \
by sqrt(1-famp^2).");
  set_param_descr(pt, PAR_USE_DEFOCUS_CORR, "A thick sample is always split into a \
number of slices to account for the varying defocus within the sample. This parameter \
tells the simulator to also correct for the varying defocus within each slice when \
simulating the interaction with this particle type. The computational work increases \
and the effect on the end result is usually negligible.");
  set_param_descr(pt, PAR_SOURCE, "Defines what kind of source of information is used \
to obtain a potential map. Possible sources of the potential map are data in a PDB \
file, reading a map directly from an MRC file or a raw data file, or generating a \
random map.");
  set_param_descr(pt, PAR_PDB_FILE_IN, "The name of a PDB file used to generate a \
potential map if source is set to \"" PAR_SOURCE__PDB "\".");
  set_param_descr(pt, PAR_PDB_TRANSF_FILE_IN, "The name of a text file which defines \
a number of geometric transformations applied to the atomic coordinates in the PDB \
file. This can be used to define particles built up of several identical subunits. \
If the parameter is not given, the atomic coordinates are used just as they appear \
in the PDB file without any transformation applied. See the manual for details on \
the file format.");
  set_param_descr(pt, PAR_ADD_HYDROGEN, "Tells the simulator to add hydrogen atoms, \
which are usually not specified in PDB files, at the appropriate locations in the \
amino acid residuals. The hydrogen atoms are added at the coordinates of the atom \
to which they are bonded. This is of course not really correct, but it is a good \
enough approximation for many purposes.");
  set_param_descr(pt, PAR_MAP_FILE_FORMAT, "The format of files used for input and \
output of potential maps. \"" PAR_FILE_FORMAT__MRC "\" means MRC file with 4-byte floating point data, \
\"" PAR_FILE_FORMAT__MRC_INT "\" means MRC file with 2-byte integer data. \"" PAR_FILE_FORMAT__RAW "\" is the same format as \
\"" PAR_FILE_FORMAT__MRC "\" but without the file header. For input files, \
\"" PAR_FILE_FORMAT__MRC "\" and \"" PAR_FILE_FORMAT__MRC_INT "\" are \
equivalent, the type of MRC file is determined by the file header. The potential \
map is assumed to be in V for floating point file formats, and mV for integer file \
formats.");
  set_param_descr(pt, PAR_MAP_AXIS_ORDER, "Controls the order in which voxel values \
appear in input and output map files. The first letter is the fastest varying index, \
and the last letter is the slowest varying. \"xyz\" means that the simulator's use \
of x, y, and z coordinates is consistent with the MRC file convention.");
  set_param_descr(pt, PAR_MAP_FILE_BYTE_ORDER, "Controls if input and output binary \
files are in little endian or big endian format. For MRC input files the parameter \
is ignored, instead the format is determined by the file header.");
  set_param_descr(pt, PAR_MAP_FILE_RE_IN, "The name of a file from which the \
potential map is read if source is set to \"" PAR_SOURCE__MAP "\".");
  set_param_descr(pt, PAR_MAP_FILE_IM_IN, "The name of a file from which the \
imaginary absorption potential map is read if source is set to \"map\" and \
use_imag_pot is set to yes.");
  set_param_descr(pt, PAR_MAP_FILE_RE_OUT, "The name of a file to which the \
potential map is written if " PAR_SOURCE " is set to \"" PAR_SOURCE__PDB "\" or \"" PAR_SOURCE__RANDOM "\". If the \
parameter is not given, the potential map is not written to file.");
  set_param_descr(pt, PAR_MAP_FILE_IM_OUT, "The name of a file to which the \
imaginary absorption potential map is written if source is set to \"" PAR_SOURCE__PDB "\" or \
\"" PAR_SOURCE__RANDOM "\" and use_imag_pot is set to yes. If the parameter is not given, the \
imaginary potential map is not written to file.");
  set_param_descr(pt, PAR_CONTRAST_RE, "The overall contrast level of a randomly \
generated potential map.");
  set_param_descr(pt, PAR_CONTRAST_IM, "The overall contrast level of a randomly \
generated imaginary absorption potential map.");
  set_param_descr(pt, PAR_SMOOTHNESS, "The smoothness of randomly generated maps. \
Higher values means smoother maps. The numerical value does not have any obvious \
physical interpretation.");
  set_param_descr(pt, PAR_MAKE_POSITIVE, "Tells the simulator to make random maps \
always have values above the background level, rather than randomly going above and \
below the background level. (Values always below the background value can be \
accomplished by setting " PAR_MAKE_POSITIVE " to " YES_STRING " and making " PAR_CONTRAST_RE " \
and " PAR_CONTRAST_IM " negative.)");
  return pt;
}

/****************************************************************************/

int particle_init(particle *p, simulation *sim){
  const char *source;
  long nx, ny, nz;
  double voxel_size;
  electronbeam *ed;
  if(p->init) return 0;
  if(particle_check_input(p)){
    WARNING("Error initializing particle: incomplete input data.\n");
    return 1;
  }
  p->init = 1;
  p->use_imag_pot = get_param_boolean(p->param, PAR_USE_IMAG_POT);
  p->use_defocus_corr = get_param_boolean(p->param, PAR_USE_DEFOCUS_CORR);
  p->rev_byte_order = reverse_byte_order(get_param_string(p->param, PAR_MAP_FILE_BYTE_ORDER));
  source = get_param_string(p->param, PAR_SOURCE);
  if(0 == strcmp(source, PAR_SOURCE__MAP)){
    if(particle_read_maps(p)){
      particle_reset(p);
      return 1;
    }
  }
  else if(0 == strcmp(source, PAR_SOURCE__RANDOM)){
    nx =  get_param_long(p->param, PAR_NX);
    ny =  get_param_long(p->param, PAR_NY);
    nz =  get_param_long(p->param, PAR_NZ);
    init_array(&(p->pot_re), nx, ny, nz);
    generate_random_particle(&(p->pot_re), get_param_double(p->param, PAR_CONTRAST_RE), 
			     get_param_double(p->param, PAR_SMOOTHNESS), get_param_boolean(p->param, PAR_MAKE_POSITIVE));
    if(p->use_imag_pot){
      init_array(&(p->pot_im), nx, ny, nz);
      generate_random_particle(&(p->pot_im), get_param_double(p->param, PAR_CONTRAST_IM), 
			       get_param_double(p->param, PAR_SMOOTHNESS), get_param_boolean(p->param, PAR_MAKE_POSITIVE));
    }
    particle_write_maps(p);
  }
  else if(0 == strcmp(source, PAR_SOURCE__PDB)){
    ed = get_electronbeam(sim, "");
    if(get_pot_from_pdb(p, ed)){
      WARNING("Could not generate potential map from pdb file.\n");
      p->init = 0;
      return 1;
    }
    write_log_comment("Particle map generated from pbd file.\n");
    particle_write_maps(p);
  }
  add_array_const(&(p->pot_re), -boundary_mean_array(&(p->pot_re)));
  if(p->use_imag_pot){
    add_array_const(&(p->pot_im), -boundary_mean_array(&(p->pot_im)));
  }
  if(p->use_defocus_corr){
    voxel_size = get_param_double(p->param, PAR_VOXEL_SIZE);
    init_array(&(p->lap_pot_re), p->pot_re.size[0], p->pot_re.size[1], p->pot_re.size[2]);
    laplace_array(&(p->pot_re), &(p->lap_pot_re));
    mult_array_const(&(p->lap_pot_re), 1/(voxel_size * voxel_size));
    if(p->use_imag_pot){
      init_array(&(p->lap_pot_im), p->pot_im.size[0], p->pot_im.size[1], p->pot_im.size[2]);
      laplace_array(&(p->pot_im), &(p->lap_pot_im));
      mult_array_const(&(p->lap_pot_im), 1/(voxel_size * voxel_size));
    }
  }
  param_table_set_lock(p->param, 1);
  write_log_comment("Particle component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void particle_reset(particle *p){
  if(0 == p->init) return;
  free_array(&(p->pot_re));
  free_array(&(p->pot_im));
  free_array(&(p->lap_pot_re));
  free_array(&(p->lap_pot_im));
  param_table_set_lock(p->param, 0);
  p->init = 0;
}

/****************************************************************************/

int particle_write_log(particle *p){
  return write_parameters_to_log(p->param);
}

/****************************************************************************/

int particle_check_input(particle *p){
  const char *source;
  boolean use_imag_pot;
  if(check_params(p->param)) return 1;
  source = get_param_string(p->param, PAR_SOURCE);
  use_imag_pot = get_param_boolean(p->param, PAR_USE_IMAG_POT);

  if(0 == strcmp(source, PAR_SOURCE__MAP)){
    if(0 == strcmp(PAR_FILE_FORMAT__RAW, get_param_string(p->param, PAR_MAP_FILE_FORMAT))){
      if(require_param(p->param, PAR_NX) || require_param(p->param, PAR_NY) || require_param(p->param, PAR_NZ)
	 || require_param(p->param, PAR_VOXEL_SIZE) || require_param(p->param, PAR_MAP_FILE_RE_IN)) return 1;
    }
    if(use_imag_pot){
      if(require_param(p->param, PAR_MAP_FILE_IM_IN)) return 1;
    }
  }
  else if(0 == strcmp(source, PAR_SOURCE__RANDOM)){
    if(require_param(p->param, PAR_NX) || require_param(p->param, PAR_NY) || require_param(p->param, PAR_NZ)
       || require_param(p->param, PAR_VOXEL_SIZE) || require_param(p->param, PAR_CONTRAST_RE) 
       || require_param(p->param, PAR_SMOOTHNESS) || require_param(p->param, PAR_MAKE_POSITIVE))
      return 1;
    if(use_imag_pot){
      if(require_param(p->param, PAR_CONTRAST_IM)) return 1;
    }
  }
  else if(0 == strcmp(source, PAR_SOURCE__PDB)){
    if(require_param(p->param, PAR_VOXEL_SIZE) || require_param(p->param, PAR_PDB_FILE_IN)) return 1;
  }
  else {
    WARNING("Undefined particle source %s.\n", source);
    return 1;
  }
  if(0 == use_imag_pot){
    if(require_param(p->param, PAR_FAMP)) return 1;
  }
  return 0;
}

/****************************************************************************/

int particle_read_maps(particle *p){
  if(0 == p->init){
    WARNING("Error reading particle maps: Particle component has not been initialized.\n");
    return 1;
  }
  if(particle_read_map(p, &(p->pot_re), get_param_string(p->param, PAR_MAP_FILE_RE_IN))) return 1;
  if(p->use_imag_pot){
    if(particle_read_map(p, &(p->pot_im), get_param_string(p->param, PAR_MAP_FILE_IM_IN))) return 1;
    if(!same_size_array(&(p->pot_re), &(p->pot_im))){
      WARNING("Real and imaginary maps must have the same size.\n");
      return 1;
    }
  }
  return 0;
}

/****************************************************************************/

int particle_read_map(particle *p, array *a, const char *fn){
  long nx, ny, nz, n;
  int ret, i;
  mrcheaderdata mrchead;
  double voxel_size;
  const char *map_axis_order = get_param_string(p->param, PAR_MAP_AXIS_ORDER);
  const char *format = get_param_string(p->param, PAR_MAP_FILE_FORMAT);
  if(0 == strcmp(format, PAR_FILE_FORMAT__RAW)){
    if(require_param(p->param, PAR_NX) || require_param(p->param, PAR_NY) || require_param(p->param, PAR_NZ)){
      WARNING("Failed to read particle map from raw file %s. The size of the map must be known.\n", fn);
      return 1;
    }
    nx =  get_param_long(p->param, PAR_NX);
    ny =  get_param_long(p->param, PAR_NY);
    nz =  get_param_long(p->param, PAR_NZ);
    init_array(a, nx, ny, nz);
    ret = read_array_float4b_raw(a, fn, map_axis_order, p->rev_byte_order);
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC) || 0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT)){
    ret = read_array_mrc(a, fn, map_axis_order, &mrchead);
    if(!ret){
      if(mrchead.mode == 1){
        mult_array_const(a, INT_FILE_POT_UNIT/POTENTIAL_UNIT);
      }
      voxel_size = 0;
      n = 0;
      for(i = 0; i < 3; i++){
	voxel_size += mrchead.cell[i];
	n += mrchead.size[i];
      }
      voxel_size /= n;
      if(voxel_size != 0.0){
	for(i = 0; i < 3; i++){
	  if(fabs(1.0 - mrchead.cell[i]/(voxel_size*mrchead.size[i])) > 0.01){
	    WARNING("Map read from MRC file %s has non-cubic voxels, which will not be handled correctly.\n", fn);
	    break;
	  }
	}
      }
      voxel_size *= ONE_ANGSTROM;
      if(param_isset(p->param, PAR_VOXEL_SIZE)){
        if(voxel_size != 0.0 && fabs(1.0 - get_param_double(p->param, PAR_VOXEL_SIZE)/voxel_size) > 0.01){
	  WARNING("Voxel size set in input does not agree with voxel size found in MRC file %s.\n", fn);
	}
      }
      else {
	if(voxel_size > 0.0){
	  set_param_double(p->param, PAR_VOXEL_SIZE, voxel_size);
	}
	else {
	  WARNING("Voxel size of MRC file %s could not be determined. Set voxel size of particle %s manually.\n", fn, p->param->name);
	  ret = 1;
	}
      }
    }
  }
  else {
    WARNING("Unknown file format %s.\n", format);
    return 1;
  }
  if(ret){
    WARNING("Error reading data from file %s.\n", fn);
    return 1;
  }
  write_log_comment("Particle map read from file %s.\n", fn);
  return ret;
}

/****************************************************************************/

int particle_write_maps(particle *p){
  int ret = 0;
  if(0 == p->init){
    WARNING("Error writing particle maps: Particle component has not been initialized.\n");
    return 1;
  }
  if(param_isset(p->param, PAR_MAP_FILE_RE_OUT)){
    if(particle_write_map(p, &(p->pot_re), get_param_string(p->param, PAR_MAP_FILE_RE_OUT))) ret = 1;
  }
  if(p->use_imag_pot && param_isset(p->param, PAR_MAP_FILE_IM_OUT)){
    if(particle_write_map(p, &(p->pot_im), get_param_string(p->param, PAR_MAP_FILE_IM_OUT))) ret = 1;
  }
  return ret;
}

/****************************************************************************/

int particle_write_map(particle *p, array *a, const char *fn){
  int ret;
  double voxel_size = get_param_double(p->param, PAR_VOXEL_SIZE);
  const char *map_axis_order = get_param_string(p->param, PAR_MAP_AXIS_ORDER);
  const char *format = get_param_string(p->param, PAR_MAP_FILE_FORMAT);
  if(0 == strcmp(format, PAR_FILE_FORMAT__RAW)){
    ret = write_array_float4b(a, fn, no_header, map_axis_order, p->rev_byte_order, voxel_size); 
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC)){
    ret = write_array_float4b(a, fn, mrc_header, map_axis_order, p->rev_byte_order, voxel_size);
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT)){
    ret = write_array_int2b(a, fn, mrc_header, map_axis_order, p->rev_byte_order, voxel_size, POTENTIAL_UNIT/INT_FILE_POT_UNIT);
  }
  else {
    WARNING("Unknown file format %s.\n", format);
    return 1;
  }
  if(ret){
    WARNING("Error writing data to file %s.\n", fn);
    return 1;
  }
  write_log_comment("Particle map written to file %s.\n", fn);
  return 0;
}

/****************************************************************************/

int particle_project(particle *p, wavefunction *wf, matrix *proj_matrix, double pos[3]){
  long dim[3], n[3], i, steps[3], k0, k1, k2, step0, step1, step2, nextcol, j0, j1;
  double w, w0, w1, w2, w3, s0, s1, ds0, ds1, r, t0, t1, z, dz0 = 0, dz1 = 0, dz2 = 0, ar, ai, famp, fph, acc_en, voxel_size;
  double *pdr, *pdi = NULL, *plr = NULL, *pli = NULL, *vd;
  vecf2d pf;
  if(0 == p->init){
    WARNING("Error in particle_project: Particle component has not been initialized.\n");
    return 1;
  }
  if(0 == wf->init){
    WARNING("Error in particle_project: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((3 != proj_matrix->m) || (3 != proj_matrix->n)){
    WARNING("Error in particle_hits_wavefunction: Wrong size of projection matrix. Matrix must be 3 x 3.\n");
    return 0;
  }
  voxel_size = get_param_double(p->param, PAR_VOXEL_SIZE);
  acc_en = electronbeam_get_acc_energy(wf->ed);

  if(fabs(get_matrix_entry(proj_matrix, 2, 2)) > fabs(get_matrix_entry(proj_matrix, 2, 0))){
    if(fabs(get_matrix_entry(proj_matrix, 2, 2)) > fabs(get_matrix_entry(proj_matrix, 2, 1))){
      /* Projection line closest to z axis */
      dim[0] = 0;
      dim[1] = 1;
      dim[2] = 2;
    }
    else {
      /* Projection line closest to y axis */
      dim[0] = 0;
      dim[1] = 2;
      dim[2] = 1;
    }
  }
  else {
    if(fabs(get_matrix_entry(proj_matrix, 2, 1)) > fabs(get_matrix_entry(proj_matrix, 2, 0))){
      /* Projection line closest to y axis */
      dim[0] = 0;
      dim[1] = 2;
      dim[2] = 1;
    }
    else {
      /* Projection line closest to x axis */
      dim[0] = 1;
      dim[1] = 2;
      dim[2] = 0;
    }
  }
  n[0] = p->pot_re.size[dim[0]];
  n[1] = p->pot_re.size[dim[1]];
  n[2] = p->pot_re.size[dim[2]];
  steps[0] = 1; steps[1] = p->pot_re.size[0]; steps[2] = (p->pot_re.size[0])*(p->pot_re.size[1]);
  step0 = steps[dim[0]];
  step1 = steps[dim[1]] - n[0]*step0;
  step2 = steps[dim[2]] - n[1]*(step1 + n[0]*step0);

  ds0 = -get_matrix_entry(proj_matrix, 2, dim[0])/get_matrix_entry(proj_matrix, 2, dim[2]);
  ds1 = -get_matrix_entry(proj_matrix, 2, dim[1])/get_matrix_entry(proj_matrix, 2, dim[2]);
  w = voxel_size * sqrt(1 + ds0*ds0 + ds1*ds1) * potential_conv_factor(acc_en);

  init_array(&(pf.values), 2, (long)ceil(n[0] + n[2]*fabs(ds0) + 1), (long)ceil(n[1] + n[2]*fabs(ds1) + 1));
  fill_array(&(pf.values), 0);
  pf.basis[0] = voxel_size*get_matrix_entry(proj_matrix, 0, dim[0]);
  pf.basis[1] = voxel_size*get_matrix_entry(proj_matrix, 1, dim[0]);
  pf.basis[2] = voxel_size*get_matrix_entry(proj_matrix, 0, dim[1]);
  pf.basis[3] = voxel_size*get_matrix_entry(proj_matrix, 1, dim[1]);
  pf.offset[0] = pos[0];
  pf.offset[1] = pos[1];

  nextcol = 2*pf.values.size[1];
  s0 = 0.5*(pf.values.size[1] - n[0]) - 0.5*(n[2]-1)*ds0;
  s1 = 0.5*(pf.values.size[2] - n[1]) - 0.5*(n[2]-1)*ds1;

  pdr = p->pot_re.data;
  if(p->use_imag_pot){
    pdi = p->pot_im.data;
  }
  if(p->use_defocus_corr){
    dz0 = -(0.5 * voxel_size / wave_number(acc_en)) * get_matrix_entry(proj_matrix, 2, dim[0]);
    dz1 = -(0.5 * voxel_size / wave_number(acc_en)) * get_matrix_entry(proj_matrix, 2, dim[1]);
    dz2 = -(0.5 * voxel_size / wave_number(acc_en)) * get_matrix_entry(proj_matrix, 2, dim[2]);
    plr = p->lap_pot_re.data;
    if(p->use_imag_pot){
      pli = p->lap_pot_im.data;
    }
  }

  for(k2 = 0; k2 < n[2]; k2++){
    j0 = (long)floor(s0); j1 = (long)floor(s1);
    t0 = s0 - j0; t1 = s1 - j1;
    w0 = w*(1-t0)*(1-t1); w1 = w*t0*(1-t1); w2 = w*(1-t0)*t1; w3 = w*t0*t1;
    if(p->use_imag_pot){
      if(p->use_defocus_corr){
	for(k1 = 0; k1 < n[1]; k1++){
	  vd = get_array_entry_ptr(&(pf.values), 0, j0, j1+k1);
	  z = -0.5*(n[0]-1)*dz0 + (k1 - 0.5*(n[1]-1))*dz1 + (k2 - 0.5*(n[2]-1))*dz2;
	  for(k0 = 0; k0 < n[0]; k0++){
	    ar = *pdr - z*(*pli);
	    ai = *pdi + z*(*plr);
	    vd[0] += w0*ar;
	    vd[2] += w1*ar;
	    vd[nextcol] += w2*ar;
	    vd[nextcol+2] += w3*ar;
	    vd[1] += w0*ai;
	    vd[3] += w1*ai;
	    vd[nextcol+1] += w2*ai;
	    vd[nextcol+3] += w3*ai;
	    pdr += step0;
	    pdi += step0;
	    plr += step0;
	    pli += step0;
	    z += dz0;
	    vd += 2;
	  }
	  pdr += step1;
	  pdi += step1;
	  plr += step1;
	  pli += step1;
	}
	pdr += step2;
	pdi += step2;
	plr += step2;
	pli += step2;
      }
      else {
	for(k1 = 0; k1 < n[1]; k1++){
	  vd = get_array_entry_ptr(&(pf.values), 0, j0, j1+k1);
	  for(k0 = 0; k0 < n[0]; k0++){
	    vd[0] += w0*(*pdr);
	    vd[2] += w1*(*pdr);
	    vd[nextcol] += w2*(*pdr);
	    vd[nextcol+2] += w3*(*pdr);
	    vd[1] += w0*(*pdi);
	    vd[3] += w1*(*pdi);
	    vd[nextcol+1] += w2*(*pdi);
	    vd[nextcol+3] += w3*(*pdi);
	    pdr += step0;
	    pdi += step0;
	    vd += 2;
	  }
	  pdr += step1;
	  pdi += step1;
	}
	pdr += step2;
	pdi += step2;
      }
    }
    else {
      if(p->use_defocus_corr){
	for(k1 = 0; k1 < n[1]; k1++){
	  vd = get_array_entry_ptr(&(pf.values), 0, j0, j1+k1);
	  z = -0.5*(n[0]-1)*dz0 + (k1 - 0.5*(n[1]-1))*dz1 + (k2 - 0.5*(n[2]-1))*dz2;
	  for(k0 = 0; k0 < n[0]; k0++){
	    vd[0] += w0*(*pdr);
	    vd[2] += w1*(*pdr);
	    vd[nextcol] += w2*(*pdr);
	    vd[nextcol+2] += w3*(*pdr);
	    ai = z*(*plr);
	    vd[1] += w0*ai;
	    vd[3] += w1*ai;
	    vd[nextcol+1] += w2*ai;
	    vd[nextcol+3] += w3*ai;
	    pdr += step0;
	    plr += step0;
	    z += dz0;
	    vd += 2;
	  }
	  pdr += step1;
	  plr += step1;
	}
	pdr += step2;
	plr += step2;
      }
      else {
	for(k1 = 0; k1 < n[1]; k1++){
	  vd = get_array_entry_ptr(&(pf.values), 0, j0, j1+k1);
	  for(k0 = 0; k0 < n[0]; k0++){
	    vd[0] += w0*(*pdr);
	    vd[2] += w1*(*pdr);
	    vd[nextcol] += w2*(*pdr);
	    vd[nextcol+2] += w3*(*pdr);
	    pdr += step0;
	    vd += 2;
	  }
	  pdr += step1;
	}
	pdr += step2;
      }
    }
    s0 += ds0; s1 += ds1;
  }

  if(0 == p->use_imag_pot){
    famp = get_param_double(p->param, PAR_FAMP);
    fph = sqrt(1 - famp * famp);
    vd = pf.values.data;
    /* Multiply by complex number (fph + i*famp) */
    for(i = 0; i < pf.values.size[1]*pf.values.size[2]; i++){
      r = vd[0]*fph - vd[1]*famp;
      vd[1] = vd[0]*famp + vd[1]*fph;
      vd[0] = r;
      vd += 2;
    }
  }
  vecf2d_add(&pf, &(wf->phase));
  free_array(&(pf.values));
  return 0;
}

/****************************************************************************/

int particle_hits_wavefunction(particle *p, wavefunction *wf, matrix *pm, double pos[3]){
  long i, j;
  double a, b, c, voxel_size;
  if(0 == p->init){
    WARNING("Error in particle_hits_wavefunction: Particle component has not been initialized.\n");
    return 0;
  }
  if(0 == wf->init){
    WARNING("Error in particle_hits_wavefunction: Wavefunction object has not been initialized.\n");
    return 0;
  }
  if((3 != pm->m) || (3 != pm->n)){
    WARNING("Error in particle_hits_wavefunction: Wrong size of projection matrix. Matrix must be 3 x 3.\n");
    return 0;
  }
  voxel_size = get_param_double(p->param, PAR_VOXEL_SIZE);
  for(j = 0; j < 2; j++){
    a = 0;
    for(i = 0; i < 3; i++){
      a += 0.5 * p->pot_re.size[i] * voxel_size 
	* fabs(wf->phase.basis[1+2*j]*get_matrix_entry(pm, 0, i) - wf->phase.basis[2*j]*get_matrix_entry(pm, 1, i));
    }
    b = fabs(wf->phase.basis[1+2*j]*(pos[0] - wf->phase.offset[0]) - wf->phase.basis[2*j]*(pos[1] - wf->phase.offset[1]));
    c = 0.5 * wf->phase.values.size[2-j] * fabs(wf->phase.basis[0] * wf->phase.basis[3] - wf->phase.basis[1] * wf->phase.basis[2]);
    if(b > a+c){
      return 0; /* particle projects outside wavefunction area */
    }
  }
  return 1; /* Particle likely to hit wavefunction area */
}

/****************************************************************************/

int generate_random_particle(array *a, double contrast, double smoothness, boolean positive){
  double xi, xj, xk, dxi, dxj, dxk, r, s, ic, jc, kc, c, c2;
  double *x, *transform;
  long n[3], i, j, k;
  fftw_plan fft;
  const long nel = nele_array(a);
  const double b0 = 1 + 2/(1 + smoothness);
  const double b1 = 0.1;
  n[0] = a->size[0]/2 + 1;
  n[1] = a->size[1];
  n[2] = a->size[2];
  dxi = (a->size[0] + a->size[1] + a->size[2])/(3.0 * a->size[0]);
  dxj = (a->size[0] + a->size[1] + a->size[2])/(3.0 * a->size[1]);
  dxk = (a->size[0] + a->size[1] + a->size[2])/(3.0 * a->size[2]);
  transform = malloc(2*n[0]*n[1]*n[2]*sizeof(double));
  fft = fftw_plan_dft_c2r_3d(a->size[2], a->size[1], a->size[0], (fftw_complex*)transform, a->data, FFTW_ESTIMATE);
  x = transform;
  for(k = 0; k < n[2]; k++){
    xk = dxk * min_l(k, n[2]-k);
    for(j = 0; j < n[1]; j++){
      xj = dxj * min_l(j, n[1]-j);
      for(i = 0; i < n[0]; i++){
	xi = dxi * i;
	s = pow(xi*xi + xj*xj + xk*xk + b0, -1-smoothness);
	*x = rand_gauss(0, s);
	x++;
	*x = rand_gauss(0, s);
	x++;
      }
    }
  }
  if(positive){
    transform[0] = pow(b0, -1-smoothness);
  }
  fftw_execute(fft);
  fftw_destroy_plan(fft);
  free(transform);
  if(positive){
    c = b1 * sqrt(norm_sq_array(a)/nel);
    c2 = c*c;
    x = a->data;
    for(i = 0; i < nel; i++){
      *x = sqrt((*x)*(*x) + c2) - c;
      x++;
    }
  }
  ic = 0.5*(a->size[0] - 1);
  jc = 0.5*(a->size[1] - 1);
  kc = 0.5*(a->size[2] - 1);
  s = 0;
  x = a->data;
  for(k = 0; k < a->size[2]; k++){
    xk = k/kc - 1;
    for(j = 0; j < a->size[1]; j++){
      xj = j/jc - 1;
      for(i = 0; i < a->size[0]; i++){
	xi = i/ic - 1;
	r = xi*xi + xj*xj + xk*xk;
	if(r < 1){
	  *x *= 1 + cos(r*r * M_PI);
	  s += fabs(*x);
	}
	else {
	  *x = 0;
	}
	x++;
      }
    }
  }
  mult_array_const(a, contrast * nel * M_PI / (6 * s));
  return 0;
}

/****************************************************************************/

int init_blank_similar_particle(particle *particle_org, particle *new_particle){
	/*with sharing parameter tables*/
	array_index_type m	=	particle_org->pot_re.size[0];
	array_index_type n	=	particle_org->pot_re.size[1];
	array_index_type o	=	particle_org->pot_re.size[2];
	/*initiating new_particle:*/
		new_particle->param	=	particle_org->param;
		init_array(&new_particle->pot_re,m,n,o);
		init_array(&new_particle->pot_im,m,n,o);
		init_array(&new_particle->lap_pot_re,m,n,o);
		init_array(&new_particle->lap_pot_im,m,n,o);
		new_particle->use_defocus_corr	=	particle_org->use_defocus_corr;
		new_particle->use_imag_pot		=	particle_org->use_imag_pot;
		new_particle->rev_byte_order	=	particle_org->rev_byte_order;
		new_particle->init				=	1;
		return 0;
}

