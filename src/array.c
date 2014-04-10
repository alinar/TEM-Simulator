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
#include "array.h"
#include "fftw3.h"
#include "log.h"
#include "misc.h"

/* #define DEBUGGING */

/****************************************************************************/

int init_array(array *a, array_index_type n1, array_index_type n2, array_index_type n3){
  array_index_type n = n1*n2*n3*sizeof(array_data_type);
  double c = ((double)n)/(((double)n1) * ((double)n2) * ((double)n3) * sizeof(array_data_type));
  if(!a) {
    WARNING("Null pointer\n");
    return 1;
  }
  if(n1 > 0 && n2 > 0 && n3 > 0){
    if(c < 0.99 || c > 1.01){
      WARNING("Integer overflow\n");
      a->data = NULL;
      return 1;
    }
    if(!(a->data = fftw_malloc(n))){
      WARNING("Allocation failed\n");
      a->ref_count = NULL;
      return 1;
    }
    if(!(a->ref_count = malloc(sizeof(int)))){
      WARNING("Allocation failed\n");
      fftw_free(a->data);
      a->data = NULL;
      return 1;
    }
    a->size[0] = n1; a->size[1] = n2; a->size[2] = n3;
    a->nalloc = n1*n2*n3;
    a->prealloc = 0;
    *(a->ref_count) = 1;
  }
  else {
    a->data = NULL;
    a->nalloc = 0;
    a->ref_count = NULL;
    a->prealloc = 0;
  }
  return 0;
}

/****************************************************************************/

int init_array_alloc(array *a, array_index_type n1, array_index_type n2, array_index_type n3, array_index_type nalloc){
  array_index_type n = n1*n2*n3*sizeof(array_data_type);
  array_index_type nb = nalloc*sizeof(array_data_type);
  double c = ((double)n)/(((double)n1) * ((double)n2) * ((double)n3) * sizeof(array_data_type));
  if(!a) {
    WARNING("Null pointer\n");
    return 1;
  }
  if(n1 > 0 && n2 > 0 && n3 > 0){
    if(c < 0.99 || c > 1.01){
      WARNING("Integer overflow\n");
      a->data = NULL;
      return 1;
    }
    c = ((double)nb)/(((double)nalloc)*sizeof(array_data_type));
    if(c < 0.99 || c > 1.01){
      WARNING("Integer overflow\n");
      a->data = NULL;
      return 1;
    }
    if(nb < n){
      nb = n;
      nalloc = n1*n2*n3;
    }
    if(!(a->data = fftw_malloc(nb))){
      WARNING("Allocation failed\n");
      a->ref_count = NULL;
      return 1;
    }
    if(!(a->ref_count = malloc(sizeof(int)))){
      WARNING("Allocation failed\n");
      fftw_free(a->data);
      a->data = NULL;
      return 1;
    }
    a->size[0] = n1; a->size[1] = n2; a->size[2] = n3;
    a->nalloc = nalloc;
    a->prealloc = 0;
    *(a->ref_count) = 1;
  }
  else {
    a->data = NULL;
    a->nalloc = 0;
    a->ref_count = NULL;
    a->prealloc = 0;
  }
  return 0;
}

/****************************************************************************/

int init_array_shared(array *a, array_index_type n1, array_index_type n2, array_index_type n3, array *b){
  array_index_type n = n1*n2*n3;
  double c = ((double)n)/(((double)n1) * ((double)n2) * ((double)n3));
  if(!a || !b) {
    WARNING("Null pointer\n");
    return 1;
  }
  if(!array_initialized(b)){
    return 1;
  }
  if(n1 > 0 && n2 > 0 && n3 > 0){
    if(c < 0.99 || c > 1.01){
      WARNING("Integer overflow\n");
      a->data = NULL;
      return 1;
    }
    if(n > b->nalloc){
      WARNING("Error in init_array_shared: Not enough allocated memory\n");
      a->data = NULL;
      return 1;
    }
    a->size[0] = n1;
    a->size[1] = n2;
    a->size[2] = n3;
    a->nalloc = b->nalloc;
    a->data = b->data;
    a->prealloc = b->prealloc;
    a->ref_count = b->ref_count;
    if (a->ref_count != NULL) {
      (*(a->ref_count))++;
    }
  }
  else {
    a->data = NULL;
    a->nalloc = 0;
    a->ref_count = NULL;
    a->prealloc = 0;
  }
  return 0;
}

/****************************************************************************/

int init_array_prealloc(array *a, array_index_type n1, array_index_type n2, array_index_type n3, void *x){
  array_index_type n = n1*n2*n3*sizeof(array_data_type);
  double c = ((double)n)/(((double)n1) * ((double)n2) * ((double)n3) * sizeof(array_data_type));
  if(!a){
    WARNING("Null pointer\n");
    return 1;
  }
  if(n1 > 0 && n2 > 0 && n3 > 0){
    if(c < 0.99 || c > 1.01){
      WARNING("Integer overflow\n");
      a->data = NULL;
      return 1;
    }
    if(!x){
      WARNING("Null pointer\n");
      return 1;
    }
    a->size[0] = n1;
    a->size[1] = n2;
    a->size[2] = n3;
    a->nalloc = n1*n2*n3;
    a->data = x;
    a->prealloc = 1;
    a->ref_count = NULL;
  }
  else{
    a->data = NULL;
    a->ref_count = NULL;
    a->prealloc = 1;
  }
  return 0;
}

/****************************************************************************/

void free_array(array *a){
  if(a){
    a->size[0] = 0; a->size[1] = 0; a->size[2] = 0;
    if((a->prealloc == 0) && (a->ref_count != NULL)){
      (*(a->ref_count))--;
      if (*(a->ref_count) <= 0) {
	free(a->ref_count);
        if(a->data != NULL){
          fftw_free(a->data);
	}
      }
    }
    a->data = NULL;
    a->ref_count = NULL;
    a->nalloc = 0;
  }
}

/****************************************************************************/

int array_initialized(const array *a){
  if(!a || !a->data || !a->size){
    WARNING("Uninitialized array\n");
    return 0;
  }
  return 1;
}

/****************************************************************************/

array_index_type nele_array(const array *a){
  if(a){
    return (a->size[0])*(a->size[1])*(a->size[2]);
  }
  else {
    return 0;
  }
}

/****************************************************************************/

int array_index_in_range(const array *a, array_index_type i, array_index_type j, array_index_type k){
#ifdef DEBUGGING
  if(!array_initialized(a)) return 0;
#endif
  return (i >= 0)&&(i < a->size[0])&&(j >= 0)&&(j < a->size[1])&&(k >= 0)&&(k < a->size[2]);
}

/****************************************************************************/

void set_array_entry(array *a, array_index_type i, array_index_type j, array_index_type k, array_data_type x){
  array_data_type *d;
  if(array_index_in_range(a,i,j,k)){
    d = a->data;
    d += i + (j + k*a->size[1])*a->size[0];
    *d = x;
  }
}

/****************************************************************************/

void add_to_array_entry(array *a, array_index_type i, array_index_type j, array_index_type k, array_data_type x){
  array_data_type *d;
  if(array_index_in_range(a,i,j,k)){
    d = a->data;
    d += i + (j + k*a->size[1])*a->size[0];
    *d += x;
  }
}

/****************************************************************************/

array_data_type get_array_entry(const array *a, array_index_type i, array_index_type j, array_index_type k){
  if(array_index_in_range(a,i,j,k)){
    return *(a->data + i + (j + k*a->size[1])*a->size[0]);
  }
  else {
    return 0;
  }
}

/****************************************************************************/

array_data_type *get_array_entry_ptr(array *a, array_index_type i, array_index_type j, array_index_type k){
  if(array_index_in_range(a,i,j,k)){
    return a->data + i + (j + k*a->size[1])*a->size[0];
  }
  else {
    return NULL;
  }
}  

/****************************************************************************/

int fill_array(array *a, array_data_type x){
  array_index_type i;
  array_data_type *d;
  if(!array_initialized(a)) return 1;
  d = a->data;
  for(i = (a->size[0])*(a->size[1])*(a->size[2]); i > 0; i--){
    *d = x;
    d++;
  }
  return 0;
}

/****************************************************************************/

int same_size_array(const array *a, const array *b){
  if(!a || !b || !a->size || !b->size){
    WARNING("Uninitialized array\n");
    return 0;
  }
  return (a->size[0] == b->size[0] && a->size[1] && b->size[1] && a->size[2] == b->size[2]);
}

/****************************************************************************/

int copy_array(const array *a, array *b){
  array_data_type *da, *db;
  array_index_type k;
  if(!array_initialized(a) || !array_initialized(b)) return 1;
  if(!same_size_array(a,b)){
    WARNING("Error in copy_array: arguments have different size\n");
    return 1;
  }
  da = a->data; db = b->data;
  for(k = nele_array(a); k > 0; k--){
    *db = *da;
    da++; db++;
  }
  return 0;
}

/****************************************************************************/

int add_array_const(array *a, array_data_type x){
  array_data_type *da;
  array_index_type k;
  if(!array_initialized(a)) return 1;
  if(x != 0.0){
    da = a->data;
    for(k = nele_array(a); k > 0; k--){
      *da += x;
      da++;
    }
  }
  return 0;
}

/****************************************************************************/

int add_array(const array *a, array *b, array_data_type x){
  array_data_type *da, *db;
  array_index_type k;
  if(!array_initialized(a) || !array_initialized(b)) return 1;
  if(!same_size_array(a,b)){
    WARNING("Error in add_array: arguments have different size\n");
  }
  da = a->data; db = b->data;
  if(x == 1.0){
    for(k = nele_array(a); k > 0; k--){
      *db += *da;
      da++; db++;
    }
  }
  else if(x == -1.0){
    for(k = nele_array(a); k > 0; k--){
      *db -= *da;
      da++; db++;
    }
  }
  else if(x != 0.0){
    for(k = nele_array(a); k > 0; k--){
      *db += x*(*da);
      da++; db++;
    }
  }
  return 0;
}

/****************************************************************************/

int add_array_offset(array *a, array *b, array_index_type i, array_index_type j, array_index_type k){
  array_index_type imin, imax, ia, jmin, jmax, ja, kmin, kmax, ka;
  array_data_type *da, *db;
  if(!array_initialized(a) || !array_initialized(b)) return 1;
  imin = max_l(0, -i);
  imax = min_l(a->size[0], b->size[0]-i);
  jmin = max_l(0, -j);
  jmax = min_l(a->size[1], b->size[1]-j);
  kmin = max_l(0, -k);
  kmax = min_l(a->size[2], b->size[2]-k);
  for(ka = kmin; ka < kmax; ka++){
    for(ja = jmin; ja < jmax; ja++){
      da = get_array_entry_ptr(a, imin, ja, ka);
      db = get_array_entry_ptr(b, i+imin, j+ja, k+ka);
      for(ia = imin; ia < imax; ia++){
	*db += *da;
	da++; db++;
      }
    }
  }
  return 0;
}

/****************************************************************************/

int mult_array_const(array *a, array_data_type x){
  array_data_type *da;
  array_index_type k;
  if(!array_initialized(a)) return 1;
  da = a->data;
  if(x == 0.0){
    fill_array(a, 0.0);
  }
  else if(x != 1.0){
    for(k = nele_array(a); k > 0; k--){
      *da *= x;
      da++;
    }
  }
  return 0;
}

/****************************************************************************/

int mult_array(const array *a, array *b){
  array_index_type i, n;
  array_data_type *da, *db;
  if(!array_initialized(a) || !array_initialized(b)) return 1;
  if(!same_size_array(a,b)){
    WARNING("Error in mult_array_elem: arguments have different size\n");
    return 1;
  }
  da = a->data;
  db = b->data;
  n = nele_array(a);
  for(i = 0; i < n; i++){
    *db *= *da;
    da++; db++;
  }
  return 0;
}

/****************************************************************************/

array_data_type norm_sq_array(const array *a){
  array_index_type i, j, k;
  array_data_type *d;
  array_data_type si, sj, sk;
  if(!array_initialized(a)){
    return 0;
  }
  d = a->data;
  sk = 0;
  for(k = 0; k < a->size[2]; k++){
    sj = 0;
    for(j = 0; j < a->size[1]; j++){
      si = 0;
      for(i = 0; i < a->size[0]; i++){
	si += (*d)*(*d);
	d++;
      }
      sj += si;
    }
    sk += sj;
  }
  return sk;
}

/****************************************************************************/

array_data_type max_array(const array *a){
  array_index_type k, n;
  array_data_type *da;
  array_data_type x;
  if(!array_initialized(a)){
    return 0;
  }
  da = a->data;
  x = *da;
  n = nele_array(a);
  for(k = 0; k < n; k++){
    if(*da > x) x = *da;
    da++;
  }
  return x;
}

/****************************************************************************/

array_data_type min_array(const array *a){
  array_index_type k, n;
  array_data_type *da;
  array_data_type x;
  if(!array_initialized(a)){
    return 0;
  }
  da = a->data;
  x = *da;
  n = nele_array(a);
  for(k = 0; k < n; k++){
    if(*da < x) x = *da;
    da++;
  }
  return x;
}

/****************************************************************************/

array_data_type mean_array(const array *a){
  array_index_type k, n;
  array_data_type x = 0;
  array_data_type *da;
  if(!array_initialized(a)){
    return 0;
  }
  da = a->data;
  n = nele_array(a);
  for(k = 0; k < n; k++){
    x += *da;
    da++;
  }
  if(n > 0){
    return x/n;
  }
  else {
    return 0;
  }
}

/****************************************************************************/

array_data_type boundary_mean_array(const array *a){
  array_index_type i, j, k, n;
  array_data_type b;
  array_data_type *ad;
  if(!array_initialized(a)){
    return 0;
  }
  ad = a->data;
  b = 0;
  n = 0;
  for(k = 1; k <= a->size[2]; k++){
    for(j = 1; j <= a->size[1]; j++){
      if((j > 1)&&(j < a->size[1])&&(k > 1)&&(k < a->size[2])){
	b += *ad;
	ad += a->size[0] - 1;
	b += *ad;
	ad++;
	n += 2;
      }
      else {
	for(i = 1; i <= a->size[0]; i++){
	  b += *ad;
	  ad++;
	  n++;
	}
      }   
    }
  }
  if(n > 0){
    return b/n;
  }
  else {
    return 0;
  }
}

/****************************************************************************/

int laplace_array(const array *a, array *b){
  array_data_type x[6];
  array_data_type *ad, *bd;
  array_index_type i, j, k, n1, n1n2;
  if(!same_size_array(a,b)){
    WARNING("Error in laplace_array: dimension mismatch\n");
    return 1;
  }
  ad = a->data;
  bd = b->data;
  n1 = a->size[0];
  n1n2 = a->size[0]*a->size[1];
  for(k = 0; k < a->size[2]; k++){
    for(j = 0; j < a->size[1]; j++){
      for(i = 0; i < a->size[0]; i++){
	if(i == 0){x[0] = 0;}
	else {x[0] = ad[-1];}
	if(i == a->size[0]-1){x[1] = 0;}
	else {x[1] = ad[1];}
	if(j == 0){x[2] = 0;}
	else {x[2] = ad[-n1];}
	if(j == a->size[1]-1){x[3] = 0;}
	else {x[3] = ad[-n1];}
	if(k == 0){x[4] = 0;}
	else {x[4] = ad[-n1n2];}
	if(k == a->size[2]-1){x[5] = 0;}
	else {x[5] = ad[n1n2];}
	*bd = x[0]+x[1]+x[2]+x[3]+x[4]+x[5] - 6*(*ad);
	ad++; bd++;
      }
    }
  }
  return 0;
}
