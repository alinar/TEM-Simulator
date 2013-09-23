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
#include <string.h>
#include "mrc.h"
#include "log.h"
#include "misc.h"

/****************************************************************************/

int check_types(void){
  if((sizeof(int4b) == 4) && (sizeof(float4b) == 4) && (sizeof(int2b) == 2)) return 0;
  WARNING("Wrong size of types used for reading and writing data files.\nThe program will not be able to read and write files correctly.\nTry changing the the typedefs in mrc.h and recompiling.\n");
  return 1;
}

/****************************************************************************/

int read_array_float4b_raw(array *a, const char *fn, const char *axis_order, int rev){
  long size[3], steps[3];
  int ret;
  FILE *fp;
  if(check_types()) return 6;
  if(set_steps(size, steps, NULL, a, axis_order)) return 3;
  FOPEN(fp, fn, "rb");
  if(fp == NULL) return 1;
  ret = read_float4b_data(a->data, fp, size, steps, rev);
  fclose(fp);
  return ret;  
}

/****************************************************************************/

int read_array_mrc(array *a, const char *fn, const char *axis_order, mrcheaderdata *m){
  long size[3], steps[3];
  int dim[3], ret;
  FILE *fp;
  if(check_types()) return 6;
  if(set_dims(dim, axis_order)) return 3;
  FOPEN(fp, fn, "rb");
  if(fp == NULL){
    return 1;
  }
  if(read_mrc_header(fp, m)){
    fclose(fp);
    return 4;
  }
  init_array(a, m->size[dim[0]], m->size[dim[1]], m->size[dim[2]]);
  if(set_steps(size, steps, NULL, a, axis_order)){
    fclose(fp);
    return 3;
  }
  if(1 == m->mode){
    ret = read_int2b_data(a->data, fp, size, steps, m->rev);
  }
  else if(2 == m->mode){
    ret = read_float4b_data(a->data, fp, size, steps, m->rev);
  }
  else {
    WARNING("No support for mrc mode %i found in file %s.\n", m->mode, fn);
    fclose(fp);
    return 5;
  }
  fclose(fp);
  return ret;
}

/****************************************************************************/

int read_float4b_data(double *a, FILE *fp, long size[3], long steps[3], int rev){
  float4b *x, *xd;
  double *ad;
  long i, j, k, n;
  size_t count;
  if(check_types()) return 6;
  n = size[0] * size[1] * size[2];
  x = malloc(n * sizeof(float4b));
  count = fread(x, sizeof(float4b), n, fp);
  if(rev){
    xd = x;
    for(i = 0; i < n; i++){
      rev_byte_order_4b(xd);
      xd++;
    }
  }
  if(count == n){
    xd = x;
    ad = a;
    for(k = 0; k < size[2]; k++){
      for(j = 0; j < size[1]; j++){
	for(i = 0; i < size[0]; i++){
	  *ad = (double)(*xd);
	  xd++;
	  ad += steps[0];
	}
	ad += steps[1];
      }
      ad += steps[2];
    }
    free(x);
    return 0;
  }
  else {
    free(x);
    return 2;
  }
}

/****************************************************************************/

int read_int2b_data(double *a, FILE *fp, long size[3], long steps[3], int rev){
  int2b *x, *xd;
  double *ad;
  long i, j, k, n;
  size_t count;
  if(check_types()) return 6;
  n = size[0] * size[1] * size[2];
  x = malloc(n * sizeof(int2b));
  count = fread(x, sizeof(int2b), n, fp);
  if(rev){
    xd = x;
    for(i = 0; i < n; i++){
      rev_byte_order_2b(xd);
      xd++;
    }
  }
  if(count == n){
    xd = x;
    ad = a;
    for(k = 0; k < size[2]; k++){
      for(j = 0; j < size[1]; j++){
	for(i = 0; i < size[0]; i++){
	  *ad = (double)(*xd);
	  xd++;
	  ad += steps[0];
	}
	ad += steps[1];
      }
      ad += steps[2];
    }
    free(x);
    return 0;
  }
  else {
    free(x);
    return 2;
  }
}

/****************************************************************************/

int write_array_float4b(const array *a, const char *fn, enum file_header_type header_type, const char *axis_order, int rev, double voxel_size){
  long size[3], steps[3];
  int ret;
  FILE *fp;
  mrcheaderdata m;
  if(check_types()) return 6;
  if(set_steps(size, steps, m.dims, a, axis_order)) return 3;
  FOPEN(fp, fn, "wb");
  if(fp == NULL){
    return 1;
  }
  if(header_type == mrc_header){ /* Write mrc header */
    voxel_size /= ONE_ANGSTROM;
    m.size[0] = (int)size[0];
    m.size[1] = (int)size[1];
    m.size[2] = (int)size[2];
    m.cell[0] = voxel_size*size[0];
    m.cell[1] = voxel_size*size[1];
    m.cell[2] = voxel_size*size[2];
    m.mode = 2;
    m.amin = min_array(a);
    m.amax = max_array(a);
    m.amean = mean_array(a);
    m.next = 0;
    m.rev = rev;
    if(write_mrc_header(fp, &m)){
      fclose(fp);
      return 4;
    }
  }
  ret = write_float4b_data(a->data, fp, size, steps, rev, NULL);
  fclose(fp);
  return ret;
}

/****************************************************************************/

int write_array_int2b(const array *a, const char *fn, enum file_header_type header_type, const char *axis_order, int rev, double voxel_size, double conv_factor){
  long size[3], steps[3];
  int ret;
  FILE *fp;
  mrcheaderdata m;
  if(check_types()) return 6;
  if(set_steps(size, steps, m.dims, a, axis_order)) return 3;
  FOPEN(fp, fn, "wb");
  if(fp == NULL){
    return 1;
  }
  if(header_type == mrc_header){ /* Write mrc header */
    voxel_size /= ONE_ANGSTROM;
    m.size[0] = (int)size[0];
    m.size[1] = (int)size[1];
    m.size[2] = (int)size[2];
    m.cell[0] = voxel_size*size[0];
    m.cell[1] = voxel_size*size[1];
    m.cell[2] = voxel_size*size[2];
    m.mode = 1;
    m.amin = conv_factor * min_array(a);
    m.amax = conv_factor * max_array(a);
    m.amean = conv_factor * mean_array(a);
    m.next = 0;
    m.rev = rev;
    if(write_mrc_header(fp, &m)){
      fclose(fp);
      return 4;
    }
  }
  ret = write_int2b_data(a->data, fp, size, steps, rev, conv_factor, NULL);
  fclose(fp);
  return ret;
}

/****************************************************************************/

int write_float4b_data(double *a, FILE *fp, long size[3], long steps[3], int rev, double min_max_mean[3]){
  float4b *x, *xd;
  double *ad;
  long i, j, k, n;
  size_t count;
  double s, t;
  if(check_types()) return 6;
  n = size[0] * size[1] * size[2];
  x = malloc(n * sizeof(float4b));
  xd = x;
  ad = a;
  for(k = 0; k < size[2]; k++){
    for(j = 0; j < size[1]; j++){
      for(i = 0; i < size[0]; i++){
	*xd = (float4b)(*ad);
	xd++;
	ad += steps[0];
      }
      ad += steps[1];
    }
    ad += steps[2];
  }
  if(NULL != min_max_mean){
    xd = x;
    min_max_mean[0] = *xd;
    min_max_mean[1] = *xd;
    min_max_mean[2] = 0;
    for(k = 0; k < size[2]; k++){
      s = 0;
      for(j = 0; j < size[1]; j++){
	t = 0;
	for(i = 0; i < size[0]; i++){
	  if(*xd < min_max_mean[0]) min_max_mean[0] = *xd;
	  if(*xd > min_max_mean[1]) min_max_mean[1] = *xd;
	  t += *xd;
	  xd++;
	}
	s += t;
      }
      min_max_mean[2] += s;
    }
    min_max_mean[2] /= n;
  }
  if(rev){
    xd = x;
    for(i = 0; i < n; i++){
      rev_byte_order_4b(xd);
      xd++;
    }
  }
  count = fwrite(x, sizeof(float4b), n, fp);
  free(x);
  if(count == n){
    return 0;
  }
  else {
    return 2;
  }
}

/****************************************************************************/

int write_int2b_data(double *a, FILE *fp, long size[3], long steps[3], int rev, double factor, double min_max_mean[3]){
  int2b *x, *xd;
  double *ad;
  long i, j, k, n;
  size_t count;
  double s, t;
  if(check_types()) return 6;
  n = size[0] * size[1] * size[2];
  x = malloc(n * sizeof(int2b));
  xd = x;
  ad = a;
  for(k = 0; k < size[2]; k++){
    for(j = 0; j < size[1]; j++){
      for(i = 0; i < size[0]; i++){
	*xd = (int2b)(*ad * factor);
	xd++;
	ad += steps[0];
      }
      ad += steps[1];
    }
    ad += steps[2];
  }
  if(NULL != min_max_mean){
    xd = x;
    min_max_mean[0] = *xd;
    min_max_mean[1] = *xd;
    min_max_mean[2] = 0;
    for(k = 0; k < size[2]; k++){
      s = 0;
      for(j = 0; j < size[1]; j++){
	t = 0;
	for(i = 0; i < size[0]; i++){
	  if(*xd < min_max_mean[0]) min_max_mean[0] = *xd;
	  if(*xd > min_max_mean[1]) min_max_mean[1] = *xd;
	  t += *xd;
	  xd++;
	}
	s += t;
      }
      min_max_mean[2] += s;
    }
    min_max_mean[2] /= n;
  }
  if(rev){
    xd = x;
    for(i = 0; i < n; i++){
      rev_byte_order_2b(xd);
      xd++;
    }
  }
  count = fwrite(x, sizeof(int2b), n, fp);
  free(x);
  if(count == n){
    return 0;
  }
  else {
    return 2;
  }
}

/****************************************************************************/

int set_dims(int dim[3], const char *axis_order){
  if(0 == strncmp(axis_order, "xyz", 3)){
    dim[0] = 0; dim[1] = 1; dim[2] = 2;
  }
  else if(0 == strncmp(axis_order, "xzy", 3)){
    dim[0] = 0; dim[2] = 1; dim[1] = 2;
  }
  else if(0 == strncmp(axis_order, "yxz", 3)){
    dim[1] = 0; dim[0] = 1; dim[2] = 2;
  }
  else if(0 == strncmp(axis_order, "yzx", 3)){
    dim[1] = 0; dim[2] = 1; dim[0] = 2;
  }
  else if(0 == strncmp(axis_order, "zxy", 3)){
    dim[2] = 0; dim[0] = 1; dim[1] = 2;
  }
  else if(0 == strncmp(axis_order, "zyx", 3)){
    dim[2] = 0; dim[1] = 1; dim[0] = 2;
  }
  else {
    WARNING("Axis order must be a permutation of xyz.\n");
    return 1;
  }
  return 0;
}

/****************************************************************************/

int set_steps(long size[3], long steps[3], int dims[3], const array *a, const char *axis_order){
  int dim[3];
  if(set_dims(dim, axis_order)) return 1;
  size[dim[0]] = a->size[0];
  size[dim[1]] = a->size[1];
  size[dim[2]] = a->size[2];
  steps[dim[0]] = 1;
  steps[dim[1]] = a->size[0];
  steps[dim[2]] = a->size[0] * a->size[1];
  steps[1] -= size[0] * steps[0];
  steps[2] -= size[1] * (size[0] * steps[0] + steps[1]);
  if(dims){
    dims[dim[0]] = 1;
    dims[dim[1]] = 2;
    dims[dim[2]] = 3;
  }
  return 0;
}

/****************************************************************************/

int write_mrc_header(FILE *fp, mrcheaderdata *m){
  int i;
  if(check_types()) return 1;
  rewind(fp);
  for(i = 0; i < 1024 + m->next; i++){
    fputc(0, fp);
  }
  if(write_int4b_field(fp, m->size[0], 0, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[1], 4, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[2], 8, m->rev)) return 1;
  if(write_int4b_field(fp, m->mode, 12, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[0], 28, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[1], 32, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[2], 36, m->rev)) return 1;
  if(write_float4b_field(fp, m->cell[0], 40, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, m->cell[1], 44, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, m->cell[2], 48, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, 90, 52, m->rev)) return 1;     /* cell angles */
  if(write_float4b_field(fp, 90, 56, m->rev)) return 1;     /* cell angles */
  if(write_float4b_field(fp, 90, 60, m->rev)) return 1;     /* cell angles */
  if(write_int4b_field(fp, m->dims[0], 64, m->rev)) return 1;  /* mapc */
  if(write_int4b_field(fp, m->dims[1], 68, m->rev)) return 1;  /* mapr */
  if(write_int4b_field(fp, m->dims[2], 72, m->rev)) return 1;  /* maps */
  if(write_float4b_field(fp, m->amin, 76, m->rev)) return 1;
  if(write_float4b_field(fp, m->amax, 80, m->rev)) return 1;
  if(write_float4b_field(fp, m->amean, 84, m->rev)) return 1;
  if(write_int4b_field(fp, m->next, 92, m->rev)) return 1;
  if(fseek(fp, 208, SEEK_SET)) return 1;
  if(EOF == fputs("MAP ", fp)) return 1;
  if(write_int2b_field(fp, 0x1144, 212, m->rev)) return 1;  /* flag for big- or little endian */
  if(write_int4b_field(fp, 1, 220, m->rev)) return 1;       /* number of labels */
  if(fseek(fp, 224, SEEK_SET)) return 1;
  if(EOF == fputs("Created by TEM-simulator, version                                               ", fp)) return 1;
  if(fseek(fp, 258, SEEK_SET)) return 1;
  if(EOF == fputs(VERSION_NUMBER, fp)) return 1;
  if(fseek(fp, 1024 + m->next, SEEK_SET)) return 1;
  return 0;
}

/****************************************************************************/

int update_mrc_header(FILE *fp, mrcheaderdata *m){
  if(check_types()) return 1;
  rewind(fp);
  if(write_int4b_field(fp, m->size[0], 0, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[1], 4, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[2], 8, m->rev)) return 1;
  if(write_int4b_field(fp, m->mode, 12, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[0], 28, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[1], 32, m->rev)) return 1;
  if(write_int4b_field(fp, m->size[2], 36, m->rev)) return 1;
  if(write_float4b_field(fp, m->cell[0], 40, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, m->cell[1], 44, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, m->cell[2], 48, m->rev)) return 1;     /* cell lengths */
  if(write_float4b_field(fp, m->amin, 76, m->rev)) return 1;
  if(write_float4b_field(fp, m->amax, 80, m->rev)) return 1;
  if(write_float4b_field(fp, m->amean, 84, m->rev)) return 1;
  return 0;
}

/****************************************************************************/

int write_mrc_header_ext(FILE *fp, mrcheaderdata *m, matrix *tiltangles){
  /* Tilt axis in degrees measured from y axis. Assumed to be the same for all images */
  double tilt_axis = -90 + get_matrix_entry(tiltangles, 0, 2) / ONE_DEGREE;
  long i;
  int ang;
  m->next = (int)(2*tiltangles->m);
  if(write_mrc_header(fp, m)) return 1;
  if(write_int2b_field(fp, 2, 128, m->rev)) return 1;       /* flag for extended header */
  if(write_int2b_field(fp, 1, 130, m->rev)) return 1;       /* flag for extended header */
  if(write_int4b_field(fp, 2, 220, m->rev)) return 1;       /* number of labels */
  if(fseek(fp, 304, SEEK_SET)) return 1;
  if(EOF == fputs("          Tilt axis rotation angle =                                            ", fp)) return 1;
  if(fseek(fp, 341, SEEK_SET)) return 1;
  if(EOF == fprintf(fp, "%5.2f", tilt_axis)) return 1;
  for(i = 0; i < tiltangles->m; i++){
    ang = round_to_int(100/ONE_DEGREE*get_matrix_entry(tiltangles, i, 1));
    write_int2b_field(fp, ang, 1024+2*i, m->rev);
  }
  return 0;
}

/****************************************************************************/

int read_mrc_header(FILE *fp, mrcheaderdata *m){
  int nx;
  const int nxmax = 10000;
  if(check_types()) return 1;
  m->rev = 0;
  if(read_int4b_field(fp, &nx, 0, m->rev)) return 1;
  if((nx < 0)||(nx > nxmax)){
    m->rev = 1;
  }
  if(read_int4b_field(fp, &(m->size[0]), 0, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->size[1]), 4, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->size[2]), 8, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->mode), 12, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->cell[0]), 40, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->cell[1]), 44, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->cell[2]), 48, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->dims[0]), 64, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->dims[1]), 68, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->dims[2]), 72, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->amin), 76, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->amax), 80, m->rev)) return 1;
  if(read_float4b_field(fp, &(m->amean), 84, m->rev)) return 1;
  if(read_int4b_field(fp, &(m->next), 92, m->rev)) return 1;
  if(fseek(fp, 1024 + m->next, SEEK_SET)) return 1;
  if((m->size[0] < 0)||(m->size[0] > nxmax)||(m->size[1] < 0)||(m->size[1] > nxmax)
     ||(m->size[2] < 0)||(m->size[2] > nxmax)){
    WARNING("Unable to determine byte order of MRC file.\n");
    return 1;
  }
  return 0;
}

/****************************************************************************/

int write_int2b_field(FILE *fp, int value, int pos, int rev){
  int2b v = (int2b)value;
  if(rev) rev_byte_order_2b(&v);
  if(fseek(fp, pos, SEEK_SET)) return 1; /* fseek failed */
  if(1 != fwrite(&v, 2, 1, fp)) return 1; /* fwrite failed */
  return 0;
}

/****************************************************************************/

int write_int4b_field(FILE *fp, int value, int pos, int rev){
  int4b v = (int4b)value;
  if(rev) rev_byte_order_4b(&v);
  if(fseek(fp, pos, SEEK_SET)) return 1; /* fseek failed */
  if(1 != fwrite(&v, 4, 1, fp)) return 1; /* fwrite failed */
  return 0;
}

/****************************************************************************/

int write_float4b_field(FILE *fp, double value, int pos, int rev){
  float4b v = (float4b)value;
  if(rev) rev_byte_order_4b(&v);
  if(fseek(fp, pos, SEEK_SET)) return 1; /* fseek failed */
  if(1 != fwrite(&v, 4, 1, fp)) return 1; /* fwrite failed */
  return 0;
}

/****************************************************************************/

int read_int4b_field(FILE *fp, int *value, int pos, int rev){
  int4b v;
  if(fseek(fp, pos, SEEK_SET)) return 1; /* fseek failed */
  if(1 != fread(&v, 4, 1, fp)) return 1; /* fread failed */
  if(rev) rev_byte_order_4b(&v);
  *value = (int)v;
  return 0;
}

/****************************************************************************/

int read_float4b_field(FILE *fp, double *value, int pos, int rev){
  float4b v;
  if(fseek(fp, pos, SEEK_SET)) return 1; /* fseek failed */
  if(1 != fread(&v, 4, 1, fp)) return 1; /* fread failed */
  if(rev) rev_byte_order_4b(&v);
  *value = (double)v;
  return 0;
}

/****************************************************************************/

void rev_byte_order_4b(void *a){
  char *p = (char*)a;
  char tmp;
  tmp = p[0];
  p[0] = p[3];
  p[3] = tmp;
  tmp = p[1];
  p[1] = p[2];
  p[2] = tmp;
}

/****************************************************************************/

void rev_byte_order_2b(void *a){
  char *p = (char*)a;
  char tmp;
  tmp = p[0];
  p[0] = p[1];
  p[1] = tmp;
}

/****************************************************************************/

int reverse_byte_order(const char *s){
  if(0 == strcmp(s, PAR_BYTE_ORDER__BE)) return machine_is_le();
  if(0 == strcmp(s, PAR_BYTE_ORDER__LE)) return !machine_is_le();
  return 0;
}

/****************************************************************************/

int machine_is_le(void){
  int4b a = 1;
  char *b = (char*)(void*)&a;
  return (int)*b;
}
