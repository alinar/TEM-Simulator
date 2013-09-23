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
#include "functions.h"
#include "matrix.h"
#include "array.h"
#include "misc.h"

/****************************************************************************/

void interp_add(double *x, long nx, long mx, double dx, double x0, double *y, long ny, long my, double dy, double y0);

void double_shear_coeff(const vecf2d *a, const vecf2d *b, const int dimsa[2], const int dimsb[2], double coeff[2][3]);

/****************************************************************************/

void vecf2d_xyrange(const vecf2d *a, double range[4]){
  double s;
  s = 0.5 * (a->values.size[1] * fabs(a->basis[0]) + a->values.size[2] * fabs(a->basis[2]));
  range[0] = a->offset[0] - s;
  range[1] = a->offset[0] + s;
  s = 0.5 * (a->values.size[1] * fabs(a->basis[1]) + a->values.size[2] * fabs(a->basis[3]));
  range[2] = a->offset[1] - s;
  range[3] = a->offset[1] + s;
}

/****************************************************************************/

void vecf2d_add(const vecf2d *a, vecf2d *b){
  int dimsa[2], dimsb[2];
  double coeff[2][3], c[4], *x, *px, *pd, x0, y0;
  long i, j, k, nx, imin, imax, jmin, jmax, stepsa[2], stepsb[2];

  for(i = 0; i < 4; i++){
    dimsa[0] = i%2; dimsa[1] = 1 - dimsa[0];
    dimsb[0]=  i/2; dimsb[1] = 1 - dimsb[0];
    double_shear_coeff(a, b, dimsa, dimsb, coeff);
    c[i] = max_d(max_d(fabs(coeff[0][0]), fabs(coeff[1][0])), max_d(fabs(1/coeff[0][0]), fabs(1/coeff[1][0])));
  }
  i = min_index(c, 4);
  dimsa[0] = i%2; dimsa[1] = 1 - dimsa[0];
  dimsb[0]=  i/2; dimsb[1] = 1 - dimsb[0];
  double_shear_coeff(a, b, dimsa, dimsb, coeff);
  stepsa[0] = a->values.size[0];
  stepsa[1] = stepsa[0] * a->values.size[1];
  stepsb[0] = b->values.size[0];
  stepsb[1] = stepsb[0] * b->values.size[1];

  imin = max_l(0, (long)floor(0.5*(a->values.size[1+dimsa[1]] - 1)
			      - fabs(coeff[1][1]) * 0.5*(b->values.size[1+dimsb[0]] - 1)
			      - fabs(coeff[1][0]) * 0.5*(b->values.size[1+dimsb[1]] - 1)
			      + coeff[1][2]));
  imax = min_l(a->values.size[1+dimsa[1]], 
	       1 + (long)ceil(0.5*(a->values.size[1+dimsa[1]] - 1)
			      + fabs(coeff[1][1]) * 0.5*(b->values.size[1+dimsb[0]] - 1)
			      + fabs(coeff[1][0]) * 0.5*(b->values.size[1+dimsb[1]] - 1)
			      + coeff[1][2]));
  jmin = max_l(0, (long)floor(0.5*(b->values.size[1+dimsb[0]] - 1)
			      - fabs(1/coeff[0][0]) * 0.5*(a->values.size[1+dimsa[0]] - 1)
			      - fabs(coeff[0][1]/coeff[0][0]) * 0.5*(a->values.size[1+dimsa[1]] - 1)
			      - coeff[0][2]/coeff[0][0]));
  jmax = min_l(b->values.size[1+dimsb[0]], 
	       1 + (long)ceil(0.5*(b->values.size[1+dimsb[0]] - 1)
			      + fabs(1/coeff[0][0]) * 0.5*(a->values.size[1+dimsa[0]] - 1)
			      + fabs(coeff[0][1]/coeff[0][0]) * 0.5*(a->values.size[1+dimsa[1]] - 1)
			      - coeff[0][2]/coeff[0][0]));
  if((imin >= imax)||(jmin >= jmax)){
    return;
  }
  nx = (imax-imin)*(jmax-jmin);
  x = malloc(nx*sizeof(double));
  for(k = 0; k < a->values.size[0]; k++){
    for(i = 0; i < nx; i++){
      x[i] = 0;
    }
    px = x;
    pd = a->values.data + k + imin * stepsa[dimsa[1]];
    x0 = 0;
    y0 = coeff[0][2] + 0.5*(jmin+jmax - b->values.size[1+dimsb[0]])*coeff[0][0] + (imin - 0.5*a->values.size[1+dimsa[1]])*coeff[0][1];
    for(i = imin; i < imax; i++){
      interp_add(pd, a->values.size[1+dimsa[0]], stepsa[dimsa[0]], 1, x0, px, jmax-jmin, imax-imin, coeff[0][0], y0);
      y0 += coeff[0][1];
      pd += stepsa[dimsa[1]];
      px++;
    }
    px = x;
    pd = b->values.data + k + jmin * stepsb[dimsb[0]];
    x0 = 0.5*(imin+imax - a->values.size[1+dimsa[1]]);
    y0 = coeff[1][2] + (jmin - 0.5*b->values.size[1+dimsb[0]])*coeff[1][1];
    for(j = jmin; j < jmax; j++){
      interp_add(px, imax-imin, 1, 1, x0, pd, b->values.size[1+dimsb[1]], stepsb[dimsb[1]], coeff[1][0], y0);
      y0 += coeff[1][1];
      pd += stepsb[dimsb[0]];
      px += imax - imin;
    }
  }
  free(x);
}

/****************************************************************************/

void interp_add(double *x, long nx, long mx, double dx, double x0, double *y, long ny, long my, double dy, double y0){
  long kx, ky;
  double x2, y2, z1, z2;
  if(dx < 0){
    x += mx*(nx-1);
    dx = -dx;
    mx = -mx;
  }
  if(dy < 0){
    y += my*(ny-1);
    dy = -dy;
    my = -my;
  }
  x2 = x0 - 0.5*dx*nx;
  y2 = y0 - 0.5*dy*ny;
  z2 = max_d(x2, y2);
  if(x2 < y2){
    kx = (long)floor((y2 - x2)/dx);
    ky = 0;
  }
  else {
    ky = (long)floor((x2 - y2)/dy);
    kx = 0;
  }
  x2 += (1+kx)*dx;
  y2 += (1+ky)*dy;
  x += mx*kx;
  y += my*ky;
  while((kx < nx)&&(ky < ny)){
    z1 = z2;
    z2 = min_d(x2, y2);
    *y += (z2-z1)/dy*(*x);
    if(x2 < y2){
      kx++;
      x2 += dx;
      x += mx;
    }
    else {
      ky++;
      y2 += dy;
      y += my;
    }
  }
}

/****************************************************************************/

void double_shear_coeff(const vecf2d *a, const vecf2d *b, const int dimsa[2], const int dimsb[2], double coeff[2][3]){
  long i, j;
  matrix u, v, w;
  init_matrix(&u, 2, 2);
  init_matrix(&v, 2, 3);
  init_matrix(&w, 2, 3);
  for(i = 0; i <= 1; i++){
    for(j = 0; j <= 1; j++){
      set_matrix_entry(&u, i, j, a->basis[i+2*dimsa[j]]);
      set_matrix_entry(&v, i, j, b->basis[i+2*dimsb[j]]);
    }
    set_matrix_entry(&v, i, 2, b->offset[i] - a->offset[i]);
  }
  solve_2x2_matrix(&u, &v, &w);
  coeff[1][0] = get_matrix_entry(&w, 1, 1);
  coeff[1][1] = get_matrix_entry(&w, 1, 0);
  coeff[1][2] = get_matrix_entry(&w, 1, 2);
  coeff[0][1] = get_matrix_entry(&w, 0, 1)/coeff[1][0];
  coeff[0][0] = get_matrix_entry(&w, 0, 0) - coeff[0][1]*coeff[1][1];
  coeff[0][2] = get_matrix_entry(&w, 0, 2) - coeff[0][1]*coeff[1][2];
  free_matrix(&u);
  free_matrix(&v);
  free_matrix(&w);
}


