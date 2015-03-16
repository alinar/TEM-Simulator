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


#ifndef STRUCTS_HEADER
#define STRUCTS_HEADER
#include "array.h"
#include "fftw3.h"
#include "functions.h"
#include "input.h"
#include "macros.h"
#include "matrix.h"
#include "mrc.h"


typedef struct {
  vecf2d count;
  array count_ft;
  fftw_plan fftplan_f;
  fftw_plan fftplan_b;
  mrcheaderdata file_header;
  double pixel;
  param_table *param;
  char init;
} detector;

typedef struct {
  matrix dose;
  param_table *param;
  char init;
} electronbeam;

enum tilt_modes {mode_tiltseries, mode_single_particle};

typedef struct {
  matrix data;
  matrix errors;
  enum tilt_modes tilt_mode;
  param_table *param;
  char init;
} geometry;

typedef struct {
  matrix defocus;
  param_table *param;
  char init;
} optics;

typedef struct {
  array pot_re;
  array pot_im;
  array lap_pot_re;
  array lap_pot_im;
  boolean use_imag_pot;
  boolean use_defocus_corr;
  int rev_byte_order;
  param_table *param;
  char init;
} particle;

typedef struct {
  matrix coordinates;
  param_table *param;
  char init;
} particleset;

typedef struct {
  param_table *param;
  char init;
} sample;

typedef struct {
  array pot_re;
  array pot_im;
  param_table *param;
  char init;
} volume;

typedef struct {
  vecf2d phase;
  array wf;
  array wf_ft;
  array inel;
  array inel_ft;
  array inel_ft_copy;
  array pathlength_ice;
  array pathlength_supp;
  vecf2d intens;
  matrix inel_ctf;
  double inel_ctf_dx;
  double pixel_size;
  fftw_plan fftplan_wf_f;
  fftw_plan fftplan_wf_b;
  fftw_plan fftplan_inel_f;
  fftw_plan fftplan_inel_b;
  electronbeam *ed;
  optics *opt;
  char init;
} wavefunction;


#define MAX_NUM_DETECTOR 10
#define MAX_NUM_ELECTRONBEAM 1
#define MAX_NUM_GEOMETRY 1
#define MAX_NUM_OPTICS 1
#define MAX_NUM_PARTICLE 20
#define MAX_NUM_PARTICLESET 20
#define MAX_NUM_SAMPLE 1
#define MAX_NUM_VOLUME 20

typedef struct {
  detector *detector[MAX_NUM_DETECTOR];
  electronbeam *electronbeam[MAX_NUM_ELECTRONBEAM];
  geometry *geometry[MAX_NUM_GEOMETRY];
  optics *optics[MAX_NUM_OPTICS];
  particle *particle[MAX_NUM_PARTICLE];
  particleset *particleset[MAX_NUM_PARTICLESET];
  sample *sample[MAX_NUM_SAMPLE];
  volume *volume[MAX_NUM_VOLUME];
  int num_detector;
  int num_electronbeam;
  int num_geometry;
  int num_optics;
  int num_particle;
  int num_particleset;
  int num_sample;
  int num_volume;
  param_table *param;
  char init;
} simulation;


#endif
