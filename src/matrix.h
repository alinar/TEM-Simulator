/*
 * Copyright 2008, Hans Rullgard, Stockholm University and 
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

#ifndef MATRIX_HEADER
#define MATRIX_HEADER
#include <stdarg.h>

typedef struct {
  long m, n;
  double *data;
} matrix;

void init_matrix(matrix *a, long m, long n);

void free_matrix(matrix *a);

long nele_matrix(const matrix *a);

void set_matrix_entry(matrix *a, long i, long j, double x);

void add_to_matrix_entry(matrix *a, long i, long j, double x);

double get_matrix_entry(const matrix *a, long i, long j);

double *get_matrix_entry_ptr(matrix *a, long i, long j);

void fill_matrix(matrix *a, double x);

void fill_matrix_diag(matrix *a, double x);

int shrink_matrix(matrix *a, long m, long n);

void swap_matrix_data(matrix *a, matrix *b);

void copy_matrix(const matrix *a, matrix *b);

void add_matrix(matrix *a, matrix *b, double x);

double sc_prod_matrix(matrix *a, matrix *b);

double dist_sq_matrix(matrix *a, matrix *b);

double norm_sq_matrix(matrix *a);

void rotate3drows(matrix *a, double x, double y, double z, double v);

void matrix_mult(const matrix *a, const matrix *b, matrix *c);

void solve_2x2_matrix(const matrix *a, const matrix *b, matrix *c);

int read_matrix_text(matrix *a, const char *fn);

int read_matrix_text_conv(matrix *a, const char *fn, long nconv, const double *conv);

int write_matrix_text_valist(const matrix *a, const char *fn, long nconv, const double *conv, long nhead, va_list header);

int write_matrix_text(const matrix *a, const char *fn, long nhead, ...);

int write_matrix_text_conv(const matrix *a, const char *fn, long nconv, const double *conv, long nhead, ...);

#endif
