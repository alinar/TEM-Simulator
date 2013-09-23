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

#ifndef MRC_HEADER
#define MRC_HEADER
#include "array.h"
#include "matrix.h"

typedef struct {
  int size[3];
  int mode;
  double cell[3];
  int dims[3];
  double amin;
  double amax;
  double amean;
  int next;
  int rev;
} mrcheaderdata;

typedef int int4b;     /* This must be a 4 byte integer type */
typedef short int2b;   /* This must be a 2 byte integer type */
typedef float float4b; /* This must be a 4 byte floating point type */

enum file_header_type {no_header, mrc_header};

int read_array_float4b_raw(array *a, const char *fn, const char *axis_order, int rev);

int read_array_mrc(array *a, const char *fn, const char *axis_order, mrcheaderdata *m);

int read_float4b_data(double *a, FILE *fp, long size[3], long steps[3], int rev);

int read_int2b_data(double *a, FILE *fp, long size[3], long steps[3], int rev);

int write_array_float4b(const array *a, const char *fn, enum file_header_type header_type, const char *axis_order, int rev, double voxel_size);

int write_array_int2b(const array *a, const char *fn, enum file_header_type header_type, const char *axis_order, int rev, double voxel_size, double conv_factor);

int write_float4b_data(double *a, FILE *fp, long size[3], long steps[3], int rev, double min_max_mean[3]);

int write_int2b_data(double *a, FILE *fp, long size[3], long steps[3], int rev, double factor, double min_max_mean[3]);

int set_dims(int dim[3], const char *axis_order);

int set_steps(long size[3], long steps[3], int dims[3], const array *a, const char *axis_order);

int write_mrc_header(FILE *fp, mrcheaderdata *m);

int update_mrc_header(FILE *fp, mrcheaderdata *m);

int write_mrc_header_ext(FILE *fp, mrcheaderdata *m, matrix *tiltangles);

int read_mrc_header(FILE *fp, mrcheaderdata *m);

int write_int2b_field(FILE *fp, int value, int pos, int rev);

int write_int4b_field(FILE *fp, int value, int pos, int rev);

int write_float4b_field(FILE *fp, double value, int pos, int rev);

int read_int4b_field(FILE *fp, int *value, int pos, int rev);

int read_float4b_field(FILE *fp, double *value, int pos, int rev);

void rev_byte_order_4b(void *a);

void rev_byte_order_2b(void *a);

int reverse_byte_order(const char *s);

int machine_is_le(void);

#endif
