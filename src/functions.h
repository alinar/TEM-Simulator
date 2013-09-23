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

#ifndef FUNCTIONS
#define FUNCTIONS
#include "array.h"

/***********************************************************************
 * The vecf2d struct represents a vector valued (for example complex
 * valued) function of two variables. The function is represented by its
 * values on a lattice of points in the plane: If values is an k x m x n 
 * array, values(:,i,j) is the value of the function at the point 
 * ((i-i0)*basis[0]+(j-j0)*basis[2]+offset[0], (i-i0)*basis[1]+(j-j0)*basis[3]+offset[1])
 * where i0 = 0.5*(m-1), j0 = 0.5*(n-1).
 ***********************************************************************/

typedef struct {
  array values;
  double basis[4];
  double offset[2];
} vecf2d;


/***********************************************************************
 * Function:  vecf2d_xyrange
 * Purpose:   Compute the range of x and y coordinates in which the 
 *            function is sampled.
 * Arguments: a - pointer to the vecf2d object.
 *            range - after return, contains {xmin, xmax, ymin, ymax}.
 */

void vecf2d_xyrange(const vecf2d *a, 
                    double range[4]);


/***********************************************************************
 * Function:  vecf2d_add
 * Purpose:   Add the function a to b, using linear interpolation to 
 *            determine the values of a at the sampling points of b.
 *            a is treated as 0 outside the parallelogram where it is
 *            sampled.
 * Arguments: a, b - pointers to the vecf2d objects.
 */

void vecf2d_add(const vecf2d *a, 
                vecf2d *b);

#endif
