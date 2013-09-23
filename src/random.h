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

#ifndef RANDOM_HEADER
#define RANDOM_HEADER

/*
 * rand_uniform(a, b)
 * uniform random variable in a <= x < b.
 */
double rand_uniform(double a, double b);

/*
 * rand_poisson(a)
 * Poisson random variable with mean a.
 */
long rand_poisson(double a);

/*
 * rand_gauss(a, b)
 * Gaussian random variable with mean a and standard deviation b.
 */
double rand_gauss(double a, double b);

double rand_wald(double a, double b);

/*
 * rand_exp(a)
 * Exponential random variable with mean a.
 */
double rand_exp(double a);

/*
 * rand_orientation(double angles[3])
 * Generate Euler angles of random rotation
 */
void rand_orientation(double angles[3]);


#endif
