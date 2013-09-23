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
#include "random.h"

/****************************************************************************/

double rand_uniform(double a, double b){
  return a + (b-a)*rand()/(RAND_MAX + 1.0);
}

/****************************************************************************
 *
 * rand_poisson
 * This implementation of poisson random variable follows 
 * the algorithm developed in 
 * J. H. Ahrens, U. Dieter, Computer generation of poisson deviates from
 * modified normal distributions, ACM Transactions on Mathematical
 * Software, vol. 8, no. 2, June 1982, 163--179.
 *--------------------------------------------------------------------------*/

long rand_poisson(double a){
  long j, l, k;
  double s, d, u, omega, b1, b2, c, c0, c1, c2, c3, px, py, fx, fy, e = 0, x, x2, t, v, delta;
  if(a < 10){
    j = 0;
    while(a >= 0){
      j++;
      s = (1.0 + rand())/(RAND_MAX + 2.0);
      a += log(s);
    }
    return j-1;
  }
  else {
    s = sqrt(a);
    d = 6*a*a;
    l = (long)floor(a - 1.1484);
    k = (long)floor(rand_gauss(a, s));
    if(k >= l){
      return k;
    }
    else if(k >= 0){
      u = rand_uniform(0,1);
      if(d*u >= pow(a-k, 3)){
	return k;
      }
      omega = 1/(s*sqrt(2*M_PI));
      b1 = 1/(24*a);
      b2 = 0.3*b1*b1;
      c3 = b1*b2/7;
      c2 = b2 - 15*c3;
      c1 = b1 - 6*b2 + 45*c3;
      c0 = 1 - b1 + 3*b2 - 15*c3;
      c = 0.1069/a;
      if(k < 10){
	px = -a;
	py = 1;
	for(j = 1; j <= k; j++){
	  py *= a/j;
	}
      }
      else {
	delta = 1.0/(12.0*k); delta -= 4.8*pow(delta, 3);
	v = (a - k)/k;
	if(fabs(v) <= 0.25){
	  px = (((((((0.12500596*v - 0.13847944)*v + 0.14218783)*v - 0.16612694)*v + 0.20001178)*v - 0.25000678)*v + 0.33333328)*v - 0.49999999)*k*v*v - delta;
	}
	else {
	  px = k*log(1 + v) - a + k - delta;
	}
	py = 1/sqrt(2*M_PI*k);
      }
      x = (k - a + 0.5)/s;
      x2 = x*x;
      fx = -0.5*x2;
      fy = omega * (((c3*x2 + c2)*x2 + c1)*x2 + c0);
      if(fy*(1-u) <= py*exp(px - fx)){
	return k;
      }
    }
    else {
      omega = 1/(s*sqrt(2*M_PI));
      b1 = 1/(24*a);
      b2 = 0.3*b1*b1;
      c3 = b1*b2/7;
      c2 = b2 - 15*c3;
      c1 = b1 - 6*b2 + 45*c3;
      c0 = 1 - b1 + 3*b2 - 15*c3;
      c = 0.1069/a;
    }
    while(1){
      t = -1;
      while(t < -0.6744){
	e = rand_exp(1);
	t = (rand_uniform(-1,1) > 0)?(1.8+e):(1.8-e);
      }
      k = (long)floor(a + s*t);
      if(k < 10){
	px = -a;
	py = 1;
	for(j = 1; j <= k; j++){
	  py *= a/j;
	}
      }
      else {
	delta = 1.0/(12.0*k); delta -= 4.8*pow(delta, 3);
	v = (a - k)/k;
	if(fabs(v) <= 0.25){
	  px = (((((((0.12500596*v - 0.13847944)*v + 0.14218783)*v - 0.16612694)*v + 0.20001178)*v - 0.25000678)*v + 0.33333328)*v - 0.49999999)*k*v*v - delta;
	}
	else {
	  px = k*log(1 + v) - a + k - delta;
	}
	py = 1/sqrt(2*M_PI*k);
      }
      x = (k - a + 0.5)/s;
      x2 = x*x;
      fx = -0.5*x2;
      fy = omega * (((c3*x2 + c2)*x2 + c1)*x2 + c0);
      if(c*rand_uniform(0,1) <= py * exp(px + e) - fy * exp(fx + e)){
	return k;
      }
    }
  }
}

/****************************************************************************/

double rand_gauss(double a, double b){
  double s, t;
  s = (1.0 + rand())/(RAND_MAX + 2.0);
  t = rand()*2.0*M_PI/(RAND_MAX + 1.0);
  return a + b*sqrt(-2.0*log(s))*cos(t);
}

/****************************************************************************/

/* Wald or inverse Gaussian distribution with mean a and variance b */
double rand_wald(double a, double b){
  double c, x, u;
  if((a <= 0) || (b < 0)) return 0; /* Illegal values */
  if(b < 1e-6) return a;
  c = (a*a*a)/b;
  x = rand_gauss(0, 1);
  x *= x;
  x *= a;
  x = a + a/(2*c) * (x - sqrt(x*(4*c + x)));
  u = rand_uniform(0,a+x);
  if(u <= a) return x;
  else return a*a/x;
}

/****************************************************************************/

double rand_exp(double a){
  double s;
  s = (1.0 + rand())/(RAND_MAX + 2.0);
  return -a*log(s);
}

/****************************************************************************/

void rand_orientation(double angles[3]){
  angles[0] = rand_uniform(0, 2*M_PI);
  angles[1] = acos(rand_uniform(-1,1));
  angles[2] = rand_uniform(0, 2*M_PI);
}
