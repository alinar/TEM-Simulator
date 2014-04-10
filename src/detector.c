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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "detector.h"
#include "array.h"
#include "fftw3.h"
#include "functions.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "misc.h"
#include "optics.h"
#include "random.h"
#include "simulation.h"


/***********************************************************************
 * Function:  detector_param_table
 * Purpose:   Initialize table of input parameters used by detector.
 * Arguments: name - Name used in log file and error messages.
 * Return:    Pointer to parameter table.
 */

param_table *detector_param_table(const char *name);


/***********************************************************************
 * Function:  detector_check_input
 * Purpose:   Check validity of input parameters for detector object.
 * Arguments: d - Pointer to detector struct.
 * Return:    0 if input is OK, nonzero otherwise.
 */

int detector_check_input(detector *d);


/***********************************************************************
 * Function:  detector_create_image_file
 * Purpose:   Write header part of image output file. Called during 
 *            initialization. Images are added to the file (and the 
 *            header is updated if necessary) by subsequent calls to 
 *            detector_write_image.
 * Arguments: d - Pointer to detector struct.
 *            sim - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure.
 */

int detector_create_image_file(detector *d,
			       simulation *sim);


/****************************************************************************/

detector *new_detector(const char *name){
  detector *d = malloc(sizeof(detector));
  d->param = detector_param_table(name);
  d->count.values.data = NULL;
  d->count_ft.data = NULL;
  d->init = 0;
  return d;
}

/****************************************************************************/

void delete_detector(detector *d){
  detector_reset(d);
  delete_param_table(d->param);
  free(d);
}

/****************************************************************************/

param_table *detector_param_table(const char *name){
  param_table *pt = new_param_table(18, TYPE_DETECTOR, name);
  add_param_req_constr(pt, PAR_DET_PIX_X, "l", 1, 1e4);
  add_param_req_constr(pt, PAR_DET_PIX_Y, "l", 1, 1e4);
  add_param_def_constr(pt, PAR_PADDING, "l", "20", 0, 1e3);
  add_param_req_constr(pt, PAR_PIXEL_SIZE, "d", 1e-2, 1e3);
  add_param_req_constr(pt, PAR_GAIN, "d", 0, HUGE_VAL);
  add_param_req_constr(pt, PAR_DQE, "d", 0.01, 1);
  add_param_req(pt, PAR_USE_QUANTIZATION, "b");
  add_param_req_constr(pt, PAR_MTF_A, "d", 0, 1);
  add_param_req_constr(pt, PAR_MTF_B, "d", 0, 1);
  add_param_req_constr(pt, PAR_MTF_C, "d", 0, 1);
  add_param_req_constr(pt, PAR_MTF_ALPHA, "d", 0, HUGE_VAL);
  add_param_req_constr(pt, PAR_MTF_BETA, "d", 0, HUGE_VAL);
  add_param_def_constr(pt, PAR_MTF_P, "d", "1", 0, 10);
  add_param_def_constr(pt, PAR_MTF_Q, "d", "1", 0, 10);
  add_param_req(pt, PAR_IMAGE_FILE_OUT, "s");
  add_param_def(pt, PAR_IMAGE_FILE_FORMAT, "s," PAR_FILE_FORMAT__MRC "," PAR_FILE_FORMAT__MRC_INT "," PAR_FILE_FORMAT__RAW, 
                PAR_FILE_FORMAT__MRC);
  add_param_def(pt, PAR_IMAGE_AXIS_ORDER, "s,xy,yx", "xy");
  add_param_def(pt, PAR_IMAGE_FILE_BYTE_ORDER, "s," PAR_BYTE_ORDER__BE "," PAR_BYTE_ORDER__LE "," PAR_BYTE_ORDER__NATIVE, 
                PAR_BYTE_ORDER__NATIVE);
  set_comp_descr(pt, "The detector component simulates the formation of \
an image once the electron wave reaches the scintillator and writes the \
images to a file. One or several detectors with different properties can \
be used in the same simulation.");
  set_param_descr(pt, PAR_DET_PIX_X, "The number of detector pixels in \
the x direction.");
  set_param_descr(pt, PAR_DET_PIX_Y, "The number of detector pixels in \
the y direction.");
  set_param_descr(pt, PAR_PADDING, "Number of pixels of padding added \
on each side of the detector in the internal computations.");
  set_param_descr(pt, PAR_PIXEL_SIZE, "Physical size of the detector \
pixels in micrometer."); 
  set_param_descr(pt, PAR_GAIN, "The average number of counts produced \
by each electron striking the detector.");
  set_param_descr(pt, PAR_DQE, "Detector quantum efficiency");
  set_param_descr(pt, PAR_USE_QUANTIZATION, "Controls whether the electron \
wave is quantized into discrete electrons. This is the primary source of \
noise in a real detector.");
  set_param_descr(pt, PAR_MTF_A, "One of the parameters controlling the \
modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_B, "One of the parameters controlling the \
modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_C, "One of the parameters controlling the \
modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_ALPHA, "One of the parameters controlling \
the modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_BETA, "One of the parameters controlling \
the modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_P, "One of the parameters controlling the \
modulation transfer function of the detector.");
  set_param_descr(pt, PAR_MTF_Q, "One of the parameters controlling the \
modulation transfer function of the detector.");
  set_param_descr(pt, PAR_IMAGE_FILE_OUT, "Name of the file to which \
images are written.");
  set_param_descr(pt, PAR_IMAGE_FILE_FORMAT, "File format to use for \
output images.");
  set_param_descr(pt, PAR_IMAGE_AXIS_ORDER, "Controls the order in which \
pixel values are written in the output file. \"xy\" means that x is the \
fastest varying index.");
  set_param_descr(pt, PAR_IMAGE_FILE_BYTE_ORDER, "Controls if the output \
file should be big endian or little endian.");
  return pt;
}

/****************************************************************************/

int detector_init(detector *d, simulation *sim){
  return detector_init_share(d, sim, NULL, NULL);
}

/****************************************************************************/

int detector_init_share(detector *d, simulation *sim, array *count, array *count_ft){
  long nx, ny, padding;
  if(d->init) return 0;
  if(detector_check_input(d)){
    WARNING("Error initializing detector: incomplete input data.\n");
    return 1;
  }
  padding = get_param_long(d->param, PAR_PADDING);
  nx = get_param_long(d->param, PAR_DET_PIX_X) + 2*padding;
  ny = get_param_long(d->param, PAR_DET_PIX_Y) + 2*padding;
  if(count == NULL || init_array_shared(&(d->count.values), 1, nx, ny, count)){
    init_array(&(d->count.values), 1, nx, ny);
  }
  d->count.basis[0] = d->count.basis[3] = detector_get_pixel_size(d);
  d->count.basis[1] = d->count.basis[2] = 0;
  d->count.offset[0] = d->count.offset[1] = 0;
  if(count_ft == NULL || init_array_shared(&(d->count_ft), 2, 1 + nx/2, ny, count_ft)){
    init_array(&(d->count_ft), 2, 1 + nx/2, ny);
  }
  d->fftplan_f = fftw_plan_dft_r2c_2d(d->count.values.size[2], d->count.values.size[1], 
				      d->count.values.data, (fftw_complex*)d->count_ft.data,
				      FFTW_ESTIMATE);
  d->fftplan_b = fftw_plan_dft_c2r_2d(d->count.values.size[2], d->count.values.size[1], 
				      (fftw_complex*)d->count_ft.data, d->count.values.data,
				      FFTW_ESTIMATE);
  d->init = 1;
  if(detector_create_image_file(d, sim)){
    WARNING("Error creating output file %s.\n", get_param_string(d->param, PAR_IMAGE_FILE_OUT));
    detector_reset(d);
    return 1;
  }
  param_table_set_lock(d->param, 1);
  write_log_comment("Detector component initialized.\n\n");
  return 0;
}

/****************************************************************************/

int detector_init_all(simulation *sim){
  array count, count_ft;
  long ncount = 0, ncount_ft = 0, nx, ny, padding;
  int i, c = 0;
  for(i = 0; i < sim->num_detector; i++){
    if(!sim->detector[i]->init){
      padding = get_param_long(sim->detector[i]->param, PAR_PADDING);
      nx = get_param_long(sim->detector[i]->param, PAR_DET_PIX_X) + 2*padding;
      ny = get_param_long(sim->detector[i]->param, PAR_DET_PIX_Y) + 2*padding;
      ncount = max_l(ncount, nx*ny);
      ncount_ft = max_l(ncount_ft, 2*(1+nx/2)*ny);
      c++;
    }
  }
  if(c == 0){
    return 0;
  }
  if(init_array_alloc(&count, 1, 1, 1, ncount) || init_array_alloc(&count_ft, 1, 1, 1, ncount_ft)){
    return 1;
  }
  for(i = 0; i < sim->num_detector; i++){
    if(detector_init_share(sim->detector[i], sim, &count, &count_ft)){
      free_array(&count);
      free_array(&count_ft);
      return 1;
    }		   
  }
  free_array(&count);
  free_array(&count_ft);
  return 0;
}

/****************************************************************************/

void detector_reset(detector *d){
  if(0 == d->init) return;
  fftw_destroy_plan(d->fftplan_f);
  fftw_destroy_plan(d->fftplan_b);
  free_array(&(d->count.values));
  free_array(&(d->count_ft));
  d->init = 0;
  param_table_set_lock(d->param, 0);
}

/****************************************************************************/

int detector_check_input(detector *d){
  return check_params(d->param);
}

/****************************************************************************/

int detector_write_log(detector *d){
  return write_parameters_to_log(d->param);
}

/****************************************************************************/

int detector_create_image_file(detector *d, simulation *sim){
  FILE *fp;
  const char *fn;
  const char *format;
  int i, fmt = 0;
  optics *o;
  geometry *g = NULL;
  if(0 == d->init){
    WARNING("Error creating image file: Detector component has not been initialized.\n");
    return 1;
  }
  fn = get_param_string(d->param, PAR_IMAGE_FILE_OUT);
  format = get_param_string(d->param, PAR_IMAGE_FILE_FORMAT);
  if(0 == strcmp(format, PAR_FILE_FORMAT__MRC)){
    fmt = 1;
  }
  else if(0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT)){
    fmt = 2;
  }
  if(fmt > 0){
    o = get_optics(sim, "");
    if(o == NULL){
      WARNING("Optics required for initialization of detector.\n");
      return 1;
    }
    g = get_geometry(sim, "");
    if(g == NULL){
      WARNING("Geometry required for initialization of detector.\n");
      return 1;
    }
    if(optics_init(o, sim) || geometry_init(g, sim)) return 1;
    d->pixel = detector_get_pixel_size(d) / get_param_double(o->param, PAR_MAGNIFICATION) / ONE_ANGSTROM;
  }
  else {
    d->pixel = 0;
  }
  i = 0;
  if(0 == strcmp("yx", get_param_string(d->param, PAR_IMAGE_AXIS_ORDER))) i = 1;
  d->file_header.size[i] = get_param_long(d->param, PAR_DET_PIX_X);
  d->file_header.size[1-i] = get_param_long(d->param, PAR_DET_PIX_Y);
  d->file_header.size[2] = 0;
  d->file_header.mode = 2;
  d->file_header.cell[0] = d->file_header.size[0] * d->pixel;
  d->file_header.cell[1] = d->file_header.size[1] * d->pixel;
  d->file_header.cell[2] = 0;
  d->file_header.dims[i] = 1;
  d->file_header.dims[1-i] = 2;
  d->file_header.dims[2] = 3;
  d->file_header.amin = 0;
  d->file_header.amax = 0;
  d->file_header.amean = 0;
  d->file_header.next = 0;
  d->file_header.rev = reverse_byte_order(get_param_string(d->param, PAR_IMAGE_FILE_BYTE_ORDER));
  if(fmt == 2) d->file_header.mode = 1;
  FOPEN(fp, fn, "wb");
  if(fp == NULL){
    WARNING("Could not open file %s for writing.\n", fn);
    return 1;
  }
  if(fmt > 0){
    if(get_param_boolean(g->param, PAR_GEN_TILT_DATA)){
      if(write_mrc_header_ext(fp, &(d->file_header), &(g->data))){
	fclose(fp);
	return 1;
      }
    }
    else {
      if(write_mrc_header(fp, &(d->file_header))){
	fclose(fp);
	return 1;
      }
    }
  }
  fclose(fp);
  return 0;
}

/****************************************************************************/

int detector_write_image(detector *d){
  long size[3], steps[3];
  double min_max_mean[3];
  int ret;
  FILE *fp;
  const char *image_axis_order;
  const char *fn;
  const char *format;
  const double factor = 1;
  long padding;
  if(0 == d->init){
    WARNING("Error writing image file: Detector component has not been initialized.\n");
    return 1;
  }
  fn = get_param_string(d->param, PAR_IMAGE_FILE_OUT);
  image_axis_order = get_param_string(d->param, PAR_IMAGE_AXIS_ORDER);
  format = get_param_string(d->param, PAR_IMAGE_FILE_FORMAT);
  padding = get_param_long(d->param, PAR_PADDING);
  size[0] = d->file_header.size[0];
  size[1] = d->file_header.size[1];
  size[2] = 1;
  if(0 == strcmp(image_axis_order, "xy")){
    steps[0] = 1;
    steps[1] = 2*padding;
    steps[2] = 0; /* This value is unimportant since size[2] = 1. */
  }
  else {
    steps[0] = size[1] + 2*padding;
    steps[1] = 1 - steps[0] * size[0];
    steps[2] = 0; /* This value is unimportant since size[2] = 1. */
  }
  FOPEN(fp, fn, "ab");
  if(fp == NULL){
    WARNING("Could not open file %s for writing.\n", fn);
    return 1;
  }
  if(0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT)){
    ret = write_int2b_data(d->count.values.data + padding*(1 + d->count.values.size[1]), fp, size, steps, 
			   d->file_header.rev, factor, min_max_mean);
  }
  else {
    ret = write_float4b_data(d->count.values.data + padding*(1 + d->count.values.size[1]), fp, size, steps, 
			     d->file_header.rev, min_max_mean);
  }
  fclose(fp);
  if(ret){
    WARNING("Error writing to file %s.\n", fn);
    return 1;
  }
  d->file_header.size[2]++;
  d->file_header.cell[2] += d->pixel;
  if(1 == d->file_header.size[2]){ /* min and max are slightly different for first image. */
    d->file_header.amin = min_max_mean[0];
    d->file_header.amax = min_max_mean[1];
    d->file_header.amean = min_max_mean[2];
  }
  else {
    d->file_header.amin = min_d(d->file_header.amin, min_max_mean[0]);
    d->file_header.amax = max_d(d->file_header.amax, min_max_mean[1]);
    d->file_header.amean = ((d->file_header.size[2] - 1) * d->file_header.amean + min_max_mean[2])/d->file_header.size[2];
  }
  if((0 == strcmp(format, PAR_FILE_FORMAT__MRC)) || (0 == strcmp(format, PAR_FILE_FORMAT__MRC_INT))){
    FOPEN(fp, fn, "rb+");
    if(update_mrc_header(fp, &(d->file_header))){
      WARNING("Error writing to file %s.\n", fn);
      fclose(fp);
      return 1;
    }
    fclose(fp);
  }
  write_log_comment("Detector image number %i written to file %s.\n", d->file_header.size[2], fn);
  return 0;
}

/****************************************************************************/

int detector_apply_quantization(detector *d){
  long i, n;
  double *x, a, c;
  if(0 == d->init){
    WARNING("Error applying quantization noise: Detector component has not been initialized.\n");
    return 1;
  }
  a = 1/get_param_double(d->param, PAR_DQE) - 1;
  n = nele_array(&(d->count.values));
  x = d->count.values.data;
  for(i = 0; i < n; i++){
    c = (double)rand_poisson(*x);
    *x = rand_wald(c, a*c);
    x++;
  }
  return 0;
}

/****************************************************************************/

int detector_apply_mtf(detector *d){
  long m, n, i, j;
  double g, dxi, dxj, xi, xj, x2, a, b, c, alpha, beta, p, q, mtf, two_p, two_q, y;
  double *ft;
  if(0 == d->init){
    WARNING("Error applying detector MTF: Detector component has not been initialized.\n");
    return 1;
  }
  m = d->count.values.size[1];
  n = d->count.values.size[2];
  dxi = 1.0/m;
  dxj = 1.0/n;
  a = get_param_double(d->param, PAR_MTF_A);
  b = get_param_double(d->param, PAR_MTF_B);
  c = get_param_double(d->param, PAR_MTF_C);
  alpha = get_param_double(d->param, PAR_MTF_ALPHA);
  beta = get_param_double(d->param, PAR_MTF_BETA);
  p = get_param_double(d->param, PAR_MTF_P);
  q = get_param_double(d->param, PAR_MTF_Q);
  g = get_param_double(d->param, PAR_GAIN)/(m * n * (a + b + c));
  two_p = pow(2.0, p);
  two_q = pow(2.0, q);
  fftw_execute(d->fftplan_f);
  ft = d->count_ft.data;
  for(j = 0; j < d->count_ft.size[2]; j++){
    xj = dxj*min_l(j, n-j);
    for(i = 0; i < d->count_ft.size[1]; i++){
      xi = dxi*min_l(i, m-i);
      x2 = xi*xi + xj*xj;
      mtf = c;
      y = 1 + alpha*x2;
      mtf += a * two_p * pow(y, -p) / (2 + (two_p - 2) * pow(y, -p));
      y = 1 + beta*x2;
      mtf += b * two_q * pow(y, -q) / (2 + (two_q - 2) * pow(y, -q));
      mtf *= g;
      *ft *= mtf;
      ft++;
      *ft *= mtf;
      ft++;
    }
  }
  fftw_execute(d->fftplan_b);
  return 0;
}

/****************************************************************************/

int detector_get_intensity(detector *d, wavefunction *wf, long tilt){
  double dose_per_pix, magnification;
  if(0 == d->init){
    WARNING("Error computing intensity: Detector component has not been initialized.\n");
    return 1;
  }
  if(0 == wf->init){
    WARNING("Error computing intensity: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((tilt < 0)||(tilt >= wf->ed->dose.m)){
    WARNING("Error computing intensity: Tilt number out of range.\n");
    return 1;
  }
  magnification = get_param_double(wf->opt->param, PAR_MAGNIFICATION);
  dose_per_pix = detector_get_pixel_size(d) / magnification;
  dose_per_pix *= dose_per_pix;                               /* Area of detector pixels */
  dose_per_pix *= get_matrix_entry(&(wf->ed->dose), tilt, 0); /* Area of pixel times dose per unit area */
  fill_array(&(d->count.values), 0);
  vecf2d_add(&(wf->intens), &(d->count));
  mult_array_const(&(d->count.values), dose_per_pix);
  return 0;
}

/****************************************************************************/

double detector_get_pixel_size(detector *d){
  return DET_PIXEL_UNIT * get_param_double(d->param, PAR_PIXEL_SIZE);
}
