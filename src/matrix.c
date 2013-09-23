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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "log.h"
#include "misc.h"

/****************************************************************************/

void init_matrix(matrix *a, long m, long n){
  a->m = m; a->n = n;
  if((a->m > 0) && (a->n > 0)){
    a->data = malloc(a->m*a->n*sizeof(double));
  }
  else {
    a->data = NULL;
  }
}

/****************************************************************************/

void free_matrix(matrix *a){
  a->m = 0; a->n = 0;
  if(a->data != NULL){
    free(a->data);
    a->data = NULL;
  }
}

/****************************************************************************/

long nele_matrix(const matrix *a){
  return (a->m)*(a->n);
}

/****************************************************************************/

void set_matrix_entry(matrix *a, long i, long j, double x){
  double *d;
  if((i < a->m)&&(j < a->n)&&(i >= 0)&&(j >= 0)){  
    d = a->data;
    d += i + j*a->m;
    *d = x;
  }
}

/****************************************************************************/

void add_to_matrix_entry(matrix *a, long i, long j, double x){
  double *d;
  if((i < a->m)&&(j < a->n)&&(i >= 0)&&(j >= 0)){  
    d = a->data;
    d += i + j*a->m;
    *d += x;
  }
}

/****************************************************************************/

double get_matrix_entry(const matrix *a, long i, long j){
  double *d;
  if((i >= a->m)||(j >= a->n)||(i < 0)||(j < 0)){
    return 0;
  }
  else {
    d = a->data;
    d += i + j*a->m;
    return *d;
  }
}

/****************************************************************************/

double *get_matrix_entry_ptr(matrix *a, long i, long j){
  return a->data + i + j*a->m;
}

/****************************************************************************/

void fill_matrix(matrix *a, double x){
  long i;
  double *d = a->data;
  for(i = (a->m)*(a->n); i > 0; i--){
    *d = x;
    d++;
  }
}

/****************************************************************************/

void fill_matrix_diag(matrix *a, double x){
  long i, n;
  fill_matrix(a, 0);
  n = min_l(a->m, a->n);
  for(i = 0; i < n; i++){
    set_matrix_entry(a, i, i, x);
  }
}

/****************************************************************************/

int shrink_matrix(matrix *a, long m, long n){
  if(m > a->m || n > a->n){
    return 1;
  }
  if(m != a->m || n != a->n){
    matrix b;
    long i, j;
    init_matrix(&b, m, n);
    for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
        set_matrix_entry(&b, i, j, get_matrix_entry(a, i, j));
      }
    }
    free_matrix(a);
    a->data = b.data;
    a->m = b.m;
    a->n = b.n;
  }
  return 0;
}

/****************************************************************************/

void swap_matrix_data(matrix *a, matrix *b){
  double *tmp;
  if((a->m)*(a->n) != (b->m)*(b->n)){
    err_message("Error in swap_matrix_data: arguments have different size");
  };
  tmp = a->data; a->data = b->data; b->data = tmp;
}

/****************************************************************************/

void copy_matrix(const matrix *a, matrix *b){
  double *da, *db;
  long k;
  if(nele_matrix(a) == nele_matrix(b)){
    b->m = a->m; b->n = a->n;
  }
  else {
    free_matrix(b);
    init_matrix(b, a->m, a->n);
  }
  da = a->data; db = b->data;
  for(k = nele_matrix(a); k > 0; k--){
    *db = *da;
    da++; db++;
  }
}

/****************************************************************************/

void add_matrix(matrix *a, matrix *b, double x){
  double *da, *db;
  long k;
  if((a->m)*(a->n) != (b->m)*(b->n)){
    err_message("Error in add_matrix: arguments have different size");
  }
  da = a->data; db = b->data;
  for(k = (a->m)*(a->n); k > 0; k--){
    *db += x*(*da);
    da++; db++;
  }
}

/****************************************************************************/

void mult_matrix_by_sc(matrix *a, double x){
  double *d;
  long k;
  d = a->data;
  for(k = (a->m)*(a->n); k > 0; k--){
    *d *= x;
    d++;
  }
}

/****************************************************************************/

double sc_prod_matrix(matrix *a, matrix *b){
  long k;
  double *da, *db;
  double s = 0;
  if((a->m)*(a->n) != (b->m)*(b->n)){
    err_message("Error in sc_prod_matrix: arguments have different size");
  };
  da = a->data; db = b->data;
  for(k = (a->m)*(a->n); k > 0; k--){
    s += (*da)*(*db);
    da++; db++;
  };
  return s;
}

/****************************************************************************/

double dist_sq_matrix(matrix *a, matrix *b){
  long k;
  double *da, *db;
  double s = 0;
  if((a->m)*(a->n) != (b->m)*(b->n)){
    err_message("Error in dist_sq_matrix: arguments have different size");
  };
  da = a->data; db = b->data;
  for(k = (a->m)*(a->n); k > 0; k--){
    s += pow((*da)-(*db),2);
    da++; db++;
  };
  return s;
}

/****************************************************************************/

double norm_sq_matrix(matrix *a){
  long k;
  double *d;
  double s = 0;
  d = a->data;
  for(k = (a->m)*(a->n); k > 0; k--){
    s += (*d)*(*d);
    d++;
  };
  return s;
}

/****************************************************************************/

void rotate3drows(matrix *a, double x, double y, double z, double v){
  double k, sn, cs;
  matrix b, c;
  if(a->n != 3){
    err_message("Error in rotate3drows: matrix must have three columns");
  }
  init_matrix(&b, 3, 3);
  init_matrix(&c, a->m, a->n);
  k = 1/sqrt(x*x + y*y + z*z);
  x *= k; y*= k; z*= k;
  cs = cos(v); sn = sin(v);
  b.data[0] = cs + x*x*(1-cs);
  b.data[1] = x*y*(1-cs)-z*sn;
  b.data[2] = x*z*(1-cs)+y*sn;
  b.data[3] = x*y*(1-cs)+z*sn;
  b.data[4] = cs + y*y*(1-cs);
  b.data[5] = y*z*(1-cs)-x*sn;
  b.data[6] = x*z*(1-cs)-y*sn;
  b.data[7] = y*z*(1-cs)+x*sn;
  b.data[8] = cs + z*z*(1-cs);
  matrix_mult(a, &b, &c);
  copy_matrix(&c, a);
  free_matrix(&b);
  free_matrix(&c);
}

/****************************************************************************/

void matrix_mult(const matrix *a, const matrix *b, matrix *c){
  long p, q, r, i, j, k;
  double *ad, *bd, *cd, *cd2;
  if((a->n != b->m)||(a->m != c->m)||(b->n != c->n)){
    err_message("Error in matrix_mult: incompatible sizes");
  }
  p = a->m; q = a->n; r = b->n;
  fill_matrix(c, 0);
  bd = b->data;
  for(k = 0; k < r; k++){
    ad = a->data;
    cd = c->data + k*p;
    for(j = 0; j < q; j++){
      cd2 = cd;
      for(i = 0; i < p; i++){
	*cd2 += (*ad)*(*bd);
	cd2++;
	ad++;
      }
      bd++;
    }
  }
}

/****************************************************************************/

void solve_2x2_matrix(const matrix *a, const matrix *b, matrix *c){
  matrix ainv;
  if((a->m != 2)||(a->n != 2)){
    err_message("Error in solve_2x2_matrix: Matrix a must be 2x2.\n");
  }
  if((b->m != 2)||(c->m != 2)||(b->n != c->n)){
    err_message("Error in solve_2x2_matrix: dimension mismatch.\n");
  }
  init_matrix(&ainv, 2, 2);
  ainv.data[0] = a->data[3];
  ainv.data[1] = -a->data[1];
  ainv.data[2] = -a->data[2];
  ainv.data[3] = a->data[0];
  mult_matrix_by_sc(&ainv, 1/(ainv.data[0]*ainv.data[3] - ainv.data[1]*ainv.data[2]));
  matrix_mult(&ainv, b, c);
  free_matrix(&ainv);
}

/****************************************************************************/

int read_matrix_text(matrix *a, const char *fn){
  long i, j, k, m, n;
  char *l1, *l2;
  char line[MAX_LINE_LENGTH];
  double x;
  FILE *fp;
  FOPEN(fp, fn, "r");
  if(fp == NULL){
    WARNING("Could not open file %s for reading.\n", fn);
    return 1;
  }
  j = 0;
  line[0] = COMMENT_CHAR;
  while(line[0] == COMMENT_CHAR){
    j++;
    if(NULL == fgetline(line, MAX_LINE_LENGTH, fp)){
      WARNING("Error reading file %s on line %i.\n", fn, j);
      fclose(fp);
      return 1;
    }
  }
  m = strtol(line, &l1, 10);
  n = strtol(l1, &l2, 10);
  if((line == l1)||(l1 == l2)){
    WARNING("Error reading file %s on line %i.\n", fn, j);
    fclose(fp);
    return 1;
  }
  if((m < 0)||(n < 0)){
    WARNING("Error reading file %s: negative matrix dimension found.\n", fn);
    fclose(fp);
    return 1;
  }
  init_matrix(a, m, n);
  i = 0;
  while(i < a->m){
    j++;
    if(NULL == fgetline(line, MAX_LINE_LENGTH, fp)){
      WARNING("Error reading file %s on line %i.\n", fn, j);
      fclose(fp);
      return 1;
    }
    if(line[0] != COMMENT_CHAR){
      l1 = line;
      for(k = 0; k < a->n; k++){
	x = strtod(l1, &l2);
	if(l1 == l2){
	  WARNING("Error reading file %s on line %i.\n", fn, j);
	  fclose(fp);
	  return 1;
	}
	l1 = l2;
	set_matrix_entry(a, i, k, x);
      }
      i++;
    }  
  }
  fclose(fp);
  return 0;
}

/****************************************************************************/

int read_matrix_text_conv(matrix *a, const char *fn, long nconv, const double *conv){
  long i, j;
  double *ad;
  if(read_matrix_text(a, fn)) return 1;
  ad = a->data;
  for(j = 0; j < a->n && j < nconv; j++){
    for(i = 0; i < a->m; i++){
      *ad *= conv[j];
      ad++;
    }
  }
  return 0;
}

/****************************************************************************/

int write_matrix_text_valist(const matrix *a, const char *fn, long nconv, const double *conv, long nhead, va_list header){
  long i, j;
  double aij;
  FILE *fp;
  FOPEN(fp, fn, "w");
  if(fp == NULL){
    WARNING("Could not open file %s for writing.\n", fn);
    return 1;
  }
  fprintf(fp, "%c File created by TEM-simulator, version %s.\n", COMMENT_CHAR, VERSION_NUMBER);
  fprintf(fp, "%i  %i\n", (int)a->m, (int)a->n);
  if(header != NULL){
    fprintf(fp, "%c ", COMMENT_CHAR);
    for(j = 0; j < nhead; j++){
      fprintf(fp, "%12.12s  ", va_arg(header, char*));
    }
    fprintf(fp, "\n");
  }
  for(i = 0; i < a->m; i++){
    for(j = 0; j < a->n; j++){
      aij = get_matrix_entry(a, i, j);
      fprintf(fp, "  %12.4f", (float)(j < nconv ? aij/conv[j] : aij));
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/****************************************************************************/

int write_matrix_text(const matrix *a, const char *fn, long nhead, ...){
  int ret;
  va_list header;
  va_start(header, nhead);
  ret = write_matrix_text_valist(a, fn, 0, NULL, nhead, header);
  va_end(header);
  return ret;
}

/****************************************************************************/

int write_matrix_text_conv(const matrix *a, const char *fn, long nconv, const double *conv, long nhead, ...){
  va_list header;
  int ret;
  va_start(header, nhead);
  ret = write_matrix_text_valist(a, fn, nconv, conv, nhead, header);
  va_end(header);
  return ret;
}
