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
#include "particleset.h"
#include "array.h"
#include "input.h"
#include "log.h"
#include "matrix.h"
#include "misc.h"
#include "particle.h"
#include "random.h"
#include "sample.h"
#include "simulation.h"

/****************************************************************************/

param_table *particleset_param_table(const char *name);

int particleset_check_input(particleset *ps);

int particleset_read_coord(particleset *ps);

int particleset_write_coord(particleset *ps);

int generate_random_coordinates(particleset *ps, 
                                simulation *sim);

/****************************************************************************/

particleset *new_particleset(const char *name){
  particleset *ps = malloc(sizeof(particleset));
  ps->init = 0;
  ps->param = particleset_param_table(name);
  return ps;
}

/****************************************************************************/

void delete_particleset(particleset *ps){
  particleset_reset(ps);
  free(ps);
}

/****************************************************************************/

param_table *particleset_param_table(const char *name){
  param_table *pt = new_param_table(16, TYPE_PARTICLESET, name);
  add_param_req(pt, PAR_PARTICLE_TYPE, "s");
  add_param_opt(pt, PAR_PARTICLE_COORDS, "s," PAR_PARTICLE_COORDS__FILE "," PAR_PARTICLE_COORDS__RANDOM "," PAR_PARTICLE_COORDS__GRID);
  add_param_opt(pt, PAR_GEN_RAND_POSITIONS, "b");
  add_param_opt(pt, PAR_WHERE, "s," PAR_WHERE__VOLUME "," PAR_WHERE__SURFACE "," PAR_WHERE__TOP "," PAR_WHERE__BOTTOM);
  add_param_opt_constr(pt, PAR_NUM_PARTICLES, "l", 0, 1e6);
  add_param_opt_constr(pt, PAR_PARTICLE_CONC, "d", 0, 1e6);
  add_param_opt_constr(pt, PAR_OCCUPANCY, "d", 0, 100);
  add_param_def_constr(pt, PAR_NX, "l", "1", 1, 1000);
  add_param_def_constr(pt, PAR_NY, "l", "1", 1, 1000);
  add_param_def_constr(pt, PAR_NZ, "l", "1", 1, 1000);
  add_param_def(pt, PAR_OFFSET_X, "d", "0");
  add_param_def(pt, PAR_OFFSET_Y, "d", "0");
  add_param_def(pt, PAR_OFFSET_Z, "d", "0");
  add_param_opt_constr(pt, PAR_PARTICLE_DIST, "d", 0, 1e4);
  add_param_opt(pt, PAR_COORD_FILE_IN, "s");
  add_param_opt(pt, PAR_COORD_FILE_OUT, "s");
  set_comp_descr(pt, "A particleset component specifies the locations of one \
or several identical particles in the sample. The particle itself is defined \
by a particle component. A simulation of micrographs can use one or several \
particleset components, or possibly none at all if an empty sample is to be \
simulated.");
  set_param_descr(pt, PAR_PARTICLE_TYPE, "The name of the particle component \
used in the particleset.");
  set_param_descr(pt, PAR_PARTICLE_COORDS, "Specifies how the particle coordinates \
should be generated. \"" PAR_PARTICLE_COORDS__FILE "\" means that positions and orientations are read \
from a text file. \"" PAR_PARTICLE_COORDS__RANDOM "\" means that positions and orientations are randomly \
generated. \"" PAR_PARTICLE_COORDS__GRID "\" means that the positions are on a regular grid, and the \
orientations are random.");
  set_param_descr(pt, PAR_GEN_RAND_POSITIONS, "This parameter is included only for \
backward compatibility and will be removed. Use " PAR_PARTICLE_COORDS " instead.");
  set_param_descr(pt, PAR_WHERE, "Possible locations of randomly placed \
particles. Available choices are the entire sample volume, near both surfaces \
of the sample, near the top surface, or near the bottom surface.");
  set_param_descr(pt, PAR_NUM_PARTICLES, "When generating random particle \
positions, one of the parameters num_particles, particle_conc, or occupancy \
must be defined. num_particles defines the number of particles to place in the \
sample. The parameter is optional if positions are read from a file, and defines \
the number of particle positions in the file to use.");
  set_param_descr(pt, PAR_PARTICLE_CONC, "When generating random particle \
positions, one of the parameters num_particles, particle_conc, or occupancy must \
be defined. particle_conc defines the number of particles per cubic micrometer \
of the sample if " PAR_WHERE " = \"" PAR_WHERE__VOLUME "\", or the number of particles per square \
micrometer on the surface otherwise.");
  set_param_descr(pt, PAR_OCCUPANCY, "When generating random particle positions, \
one of the parameters num_particles, particle_conc, or occupancy must be defined. \
occupancy defines the fraction of the sample volume occupied by particles if \
where = \"volume\", or the fraction of the sample surface occupied in other cases.");
  set_param_descr(pt, PAR_NX, "The size of the grid in the x axis direction when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_NY, "The size of the grid in the y axis direction when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_NZ, "The size of the grid in the z axis direction when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_OFFSET_X, "The x coordinate of the grid center when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_OFFSET_Y, "The y coordinate of the grid center when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_OFFSET_Z, "The z coordinate of the grid center when \
particle positions are on a regular grid.");
  set_param_descr(pt, PAR_PARTICLE_DIST, "The distance in nanometer between the \
centers of adjacent particles when particle positions are on a regular grid.");
  set_param_descr(pt, PAR_COORD_FILE_IN, "The name of a text file from which the \
positions of particles are read if they are not randomly generated. See the manual \
for details about the file format.");
  set_param_descr(pt, PAR_COORD_FILE_OUT, "The name of a text file which randomly \
generated particle positions are written to. See the manual for details about the \
file format.");
  return pt;
}

/****************************************************************************/

int particleset_init(particleset *ps, simulation *sim){
  if(ps->init) return 0;
  if(particleset_check_input(ps)){
    WARNING("Error initializing particleset: incomplete input data.\n");
    return 1;
  }
  ps->init = 1;
  if(0 == strcmp(get_param_string(ps->param, PAR_PARTICLE_COORDS), PAR_PARTICLE_COORDS__FILE)){
    if(particleset_read_coord(ps)){
      particleset_reset(ps);
      return 1;
    }
  }
  else {
    if(generate_random_coordinates(ps, sim)){
      particleset_reset(ps);
      return 1;
    }
    if(param_isset(ps->param, PAR_COORD_FILE_OUT)){
      particleset_write_coord(ps);
    }
  }
  param_table_set_lock(ps->param, 1);
  write_log_comment("Particleset component initialized.\n\n");
  return 0;
}

/****************************************************************************/

void particleset_reset(particleset *ps){
  if(0 == ps->init) return;
  free_matrix(&(ps->coordinates));
  param_table_set_lock(ps->param, 0);
  ps->init = 0;
}

/****************************************************************************/

int particleset_check_input(particleset *ps){
  if(check_params(ps->param)) return 1;
  if(param_isset(ps->param, PAR_GEN_RAND_POSITIONS) && !param_isset(ps->param, PAR_PARTICLE_COORDS)){
    if(get_param_boolean(ps->param, PAR_GEN_RAND_POSITIONS)){
      set_param(ps->param, PAR_PARTICLE_COORDS, PAR_PARTICLE_COORDS__RANDOM);
    }
    else {
      set_param(ps->param, PAR_PARTICLE_COORDS, PAR_PARTICLE_COORDS__FILE);
    }
    WARNING("Parameter " TYPE_PARTICLESET ":" PAR_GEN_RAND_POSITIONS " is obsolete and will "
	    "be removed. Use " PAR_PARTICLE_COORDS " instead.\n");
  }
  if(0 == strcmp(get_param_string(ps->param, PAR_PARTICLE_COORDS), PAR_PARTICLE_COORDS__RANDOM)){
    if(require_param(ps->param, PAR_WHERE)) return 1;
  }
  else if (0 == strcmp(get_param_string(ps->param, PAR_PARTICLE_COORDS), PAR_PARTICLE_COORDS__GRID)){
    if(require_param(ps->param, PAR_PARTICLE_DIST)) return 1;
  }
  else {
    if(require_param(ps->param, PAR_COORD_FILE_IN)) return 1;
  }
  return 0;
}

/****************************************************************************/

int particleset_write_log(particleset *ps){
  return write_parameters_to_log(ps->param);
}

/****************************************************************************/

int particleset_read_coord(particleset *ps){
  long npart;
  double conv[6];
  if(0 == ps->init){
    WARNING("Error reading particle coordinates: Particleset component has not been initialized.\n");
    return 1;
  }
  conv[0] = conv[1] = conv[2] = 1;
  conv[3] = conv[4] = conv[5] = ANGLE_UNIT;
  if(read_matrix_text_conv(&(ps->coordinates), get_param_string(ps->param, PAR_COORD_FILE_IN), 6, conv)) return 1;
  if(param_isset(ps->param, PAR_NUM_PARTICLES)){
    npart = get_param_long(ps->param, PAR_NUM_PARTICLES);
  }
  else {
    npart = ps->coordinates.m;
    set_param_long(ps->param, PAR_NUM_PARTICLES, npart);
    write_log_comment("Number of particles read from particle coordinate file: %i.\n", (int)npart);
  }
  if(shrink_matrix(&(ps->coordinates), npart, 6)){
    WARNING("Too few data in particle coordinate input file %s.\n", get_param_string(ps->param, PAR_COORD_FILE_IN));
    return 1;
  }
  return 0;
}

/****************************************************************************/

int particleset_write_coord(particleset *ps){
  double conv[6];
  if(0 == ps->init){
    WARNING("Error writing particle coordinates: Particleset component has not been initialized.\n");
    return 1;
  }
  conv[0] = conv[1] = conv[2] = 1;
  conv[3] = conv[4] = conv[5] = ANGLE_UNIT;
  return write_matrix_text_conv(&(ps->coordinates), get_param_string(ps->param, PAR_COORD_FILE_OUT),
				6, conv, 6, "x", "y", "z", "phi", "theta", "psi");
}

/****************************************************************************/

int rotate_particle(matrix *pm, particleset *ps, long i){
  if(0 == ps->init){
    WARNING("Error in rotate_particle: Particleset component has not been initialized.\n");
    return 1;
  }
  if((3 != pm->m) || (3 != pm->n)){
    WARNING("Error in rotate_particle: Projection matrix has wrong size. Matrix must be 3 x 3.\n");
    return 1;
  }
  if((i < 0) || (i >= ps->coordinates.m)){
    WARNING("Error in rotate_particle: Index i is out of range.\n");
    return 1;
  }
  rotate3drows(pm, 0, 0, 1, -get_matrix_entry(&(ps->coordinates), i, 5));
  rotate3drows(pm, 1, 0, 0, -get_matrix_entry(&(ps->coordinates), i, 4));
  rotate3drows(pm, 0, 0, 1, -get_matrix_entry(&(ps->coordinates), i, 3));
  return 0;
}

/****************************************************************************/

int get_particle_pos(double pos[3], particleset *ps, long i){
  long j;
  if(0 == ps->init){
    WARNING("Error in get_particle_pos: Particleset component has not been initialized.\n");
    return 1;
  }
  if((i < 0) || (i >= ps->coordinates.m)){
    WARNING("Error in get_particle_pos: Index i is out of range.\n");
    return 1;
  }
  for(j = 0; j < 3; j++){
    pos[j] = get_matrix_entry(&(ps->coordinates), i, j);
  }
  return 0;
}

/****************************************************************************/

int get_particle_coord(matrix *pm, double pos[3], particleset *ps, long i){
  if(get_particle_pos(pos, ps, i)) return 1;
  fill_matrix_diag(pm, 1);
  if(rotate_particle(pm, ps, i)) return 1;
  return 0;
}

/****************************************************************************/

#define IN_VOLUME 0
#define ON_SURFACES 1
#define ON_TOP_SURFACE 2
#define ON_BOTTOM_SURFACE 3

int generate_random_coordinates(particleset *ps, simulation *sim){
  long i, j, num, nx=0, ny=0, nz=0;
  double r, rad=0, rad2=0, a=0, b=0, x, y, z, zs, zmax=0, u, ori[3], dist=0,
    thickness_center=0, thickness_edge=0, diameter=0, offset_x=0, offset_y=0, offset_z=0, vol, area, vox_s;
  int grid = 0;
  int where = IN_VOLUME;
  const char *where_str;
  particle *p;
  sample *s;
  if(0 == ps->init){
    WARNING("Error in generate_random_coordinates: Particleset component has not been initialized.\n");
    return 1;
  }
  if(0 == strcmp(get_param_string(ps->param, PAR_PARTICLE_COORDS), PAR_PARTICLE_COORDS__RANDOM)){
    p = get_particle(sim, get_param_string(ps->param, PAR_PARTICLE_TYPE));
    if(p == NULL){
      WARNING("Particle %s not found.\n", get_param_string(ps->param, PAR_PARTICLE_TYPE));
      return 1;
    }
    if(particle_init(p, sim)) return 1;
    s = get_sample(sim, "");
    if(s == NULL){
      WARNING("Sample required for initialization of particleset.\n");
      return 1;
    }
    if(sample_init(s, sim)) return 1;
    where_str = get_param_string(ps->param, PAR_WHERE);
    if(0 == strcmp(where_str, PAR_WHERE__VOLUME)){
      where = IN_VOLUME;
    }
    else if(0 == strcmp(where_str, PAR_WHERE__SURFACE)){
      where = ON_SURFACES;
    }
    else if(0 == strcmp(where_str, PAR_WHERE__TOP)){
      where = ON_TOP_SURFACE;
    }
    else if(0 == strcmp(where_str, PAR_WHERE__BOTTOM)){
      where = ON_BOTTOM_SURFACE;
    }
    thickness_center = get_param_double(s->param, PAR_THICKNESS_CENTER);
    thickness_edge = get_param_double(s->param, PAR_THICKNESS_EDGE);
    diameter = get_param_double(s->param, PAR_DIAMETER);
    offset_x = get_param_double(s->param, PAR_OFFSET_X);
    offset_y = get_param_double(s->param, PAR_OFFSET_Y);
    vox_s = get_param_double(p->param, PAR_VOXEL_SIZE);
    r = (p->pot_re.size[0] + p->pot_re.size[1] + p->pot_re.size[2]) * vox_s / 6.0;
    if(param_isset(ps->param, PAR_NUM_PARTICLES)){
      num = get_param_long(ps->param, PAR_NUM_PARTICLES);
    }
    else if(param_isset(ps->param, PAR_PARTICLE_CONC)){
      if(where == IN_VOLUME){
	vol = 0.125 * M_PI * diameter * diameter * (thickness_center + thickness_edge);
	num = round_to_int(vol * VOL_CONC_UNIT * get_param_double(ps->param, PAR_PARTICLE_CONC));
      }
      else {
	area = 0.25 * M_PI * diameter * diameter;
	if(where == ON_SURFACES) area *= 2;
	num = round_to_int(area * SURF_CONC_UNIT * get_param_double(ps->param, PAR_PARTICLE_CONC));
      }
    }
    else if(param_isset(ps->param, PAR_OCCUPANCY)){
      if(where == IN_VOLUME){
	vol = 0.125 * M_PI * diameter * diameter * (thickness_center + thickness_edge);
	num = round_to_int(vol * 3 / (4 * M_PI * r * r * r) * get_param_double(ps->param, PAR_OCCUPANCY));
      }
      else {
	area = 0.25 * M_PI * diameter * diameter;
	if(where == ON_SURFACES) area *= 2;
	num = round_to_int(area / (M_PI * r * r) * get_param_double(ps->param, PAR_OCCUPANCY));
      }
    }
    else {
      WARNING("To generate random particle coordinates, one of the parameters " PAR_NUM_PARTICLES ",\n "
PAR_PARTICLE_CONC ", or " PAR_OCCUPANCY " must be defined.\n");
      return 1;
    }
    r = (p->pot_re.size[0] + p->pot_re.size[1] + p->pot_re.size[2]) * vox_s / 6.0; /*Estimate of particle radius */
    if(2*r > min_d(thickness_center, diameter)){
      WARNING("Error in generate_random_coordinates: particle to big to fit in sample");
    }
    rad = 0.5 * diameter - r;
    rad2 = rad*rad;
    a = 2*(thickness_edge - thickness_center)/(diameter * diameter);
    b = 0.5 * thickness_center - r;
    zmax = a*rad2 + b;
  }
  else {
    nx = get_param_long(ps->param, PAR_NX);
    ny = get_param_long(ps->param, PAR_NY);
    nz = get_param_long(ps->param, PAR_NZ);
    dist = get_param_double(ps->param, PAR_PARTICLE_DIST);
    offset_x = get_param_double(ps->param, PAR_OFFSET_X);
    offset_y = get_param_double(ps->param, PAR_OFFSET_Y);
    offset_z = get_param_double(ps->param, PAR_OFFSET_Z);
    num = nx*ny*nz;
    grid = 1;
  }
  init_matrix(&(ps->coordinates), num, 6);
  for(i = 0; i < num; /* empty */){
    if(grid){
      j = i/nx;
      x = (i%nx - 0.5*(nx-1))*dist + offset_x;
      y = (j%ny - 0.5*(ny-1))*dist + offset_y;
      z = (j/ny - 0.5*(nz-1))*dist + offset_z;
    }
    else {
      x = rand_uniform(-rad, rad);
      y = rand_uniform(-rad, rad);
      u = x*x + y*y;
      if(u >= rad2) continue;
      zs = a*u + b;
      z = 0;
      switch(where){
      case IN_VOLUME: /* put particles in entire volume */
	z = rand_uniform(-zmax, zmax);
	if(fabs(z) > zs){
	  continue;
	}
	break;
      case ON_SURFACES: /* put particles on both surfaces */
	z = (rand_uniform(-1,1) > 0)?zs:-zs;
	break;
      case ON_TOP_SURFACE: /* put particles on top surface */
	z = -zs;
	break;
      case ON_BOTTOM_SURFACE: /* put particles on bottom surface */
	z = zs;
      }
      x += offset_x;
      y += offset_y;
    }
    rand_orientation(ori);
    set_matrix_entry(&(ps->coordinates), i, 0, x);
    set_matrix_entry(&(ps->coordinates), i, 1, y);
    set_matrix_entry(&(ps->coordinates), i, 2, z);
    set_matrix_entry(&(ps->coordinates), i, 3, ori[0]);
    set_matrix_entry(&(ps->coordinates), i, 4, ori[1]);
    set_matrix_entry(&(ps->coordinates), i, 5, ori[2]);
    i++;
  }
  write_log_comment("Generated %i particle coordinates.\n", (int)ps->coordinates.m);
  return 0;
}

