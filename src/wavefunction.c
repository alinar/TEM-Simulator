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
#include "wavefunction.h"
#include "array.h"
#include "electronbeam.h"
#include "/opt/local/include/fftw3.h"
#include "functions.h"
#include "geometry.h"
#include "log.h"
#include "misc.h"
#include "optics.h"
#include "particle.h"
#include "particleset.h"
#include "sample.h"
#include "simulation.h"

/****************************************************************************/

wavefunction *new_wavefunction(electronbeam *ed, optics *o){
  wavefunction *wf = malloc(sizeof(wavefunction));
  wf->ed = ed;
  wf->opt = o;
  wf->phase.values.data = NULL;
  wf->wf.data = NULL;
  wf->wf_ft.data = NULL;
  wf->inel.data = NULL;
  wf->inel_ft.data = NULL;
  wf->inel_ft_copy.data = NULL;
  wf->pathlength_ice.data = NULL;
  wf->pathlength_supp.data = NULL;
  wf->intens.values.data = NULL;
  wf->init = 0;
  return wf;
}

/****************************************************************************/

void delete_wavefunction(wavefunction *wf){
  wavefunction_reset(wf);
  free(wf);
}

/****************************************************************************/

/* How the arrays are layed out in memory
 * Arrays in the same vertical column share the same memory
 *   -----------------------------------------------------------------------------
 *  |               wf                 |          phase           |     intens    |
 *   -----------------------------------------------------------------------------
 *  |              wf_ft               |   wf_inel   | wf_inel_ft |
 *   -------------------------------------------------------------
 *  | pathlength_ice | pathlength_supp |
 *   ----------------------------------
 *                   | wf_inel_ft_copy |
 *                    -----------------
 */

int wavefunction_init(wavefunction *wf, double det_pixel_size, double det_range[4]){
  long pixels_x, pixels_y;
  double magn;
  int i;
  if(wf->init) return 0;
  if(0 == wf->ed->init){
    WARNING("Error initializing wavefunction: Electronbeam member has not been initialized.\n");
    return 1;
  }
  if(0 == wf->opt->init){
    WARNING("Error initializing wavefunction: Optics member has not been initialized.\n");
    return 1;
  }
  wf->init = 1;
  magn = get_param_double(wf->opt->param, PAR_MAGNIFICATION);
  wf->pixel_size = det_pixel_size / magn;
  pixels_x = (long)ceil((det_range[1] - det_range[0])/det_pixel_size);
  pixels_y = (long)ceil((det_range[3] - det_range[2])/det_pixel_size);
  wf->phase.basis[0] = wf->pixel_size;
  wf->phase.basis[1] = 0;
  wf->phase.basis[2] = 0;
  wf->phase.basis[3] = wf->pixel_size;
  wf->phase.offset[0] = 0.5 * (det_range[0] + det_range[1]) / magn;
  wf->phase.offset[1] = 0.5 * (det_range[2] + det_range[3]) / magn;
  for(i = 0; i < 4; i++){
    wf->intens.basis[i] = magn * wf->phase.basis[i];
  }
  for(i = 0; i < 2; i++){
    wf->intens.offset[i] = magn * wf->phase.offset[i];
  }
  init_array_alloc(&(wf->wf), 2, pixels_x, pixels_y, pixels_x*pixels_y + 2*(1+pixels_x/2)*pixels_y);
  init_array_shared(&(wf->wf_ft), 2, pixels_x, pixels_y, &(wf->wf));
  init_array_shared(&(wf->pathlength_ice), 1, pixels_x, pixels_y, &(wf->wf));
  init_array_prealloc(&(wf->pathlength_supp), 1, pixels_x, pixels_y, wf->wf.data + pixels_x*pixels_y);
  init_array_prealloc(&(wf->inel_ft_copy), 2, 1+pixels_x/2, pixels_y, wf->wf.data + pixels_x*pixels_y);
  init_array_alloc(&(wf->phase.values), 2, pixels_x, pixels_y, pixels_x*pixels_y + 2*(1+pixels_x/2)*pixels_y);
  init_array_shared(&(wf->inel), 1, pixels_x, pixels_y, &(wf->phase.values));
  init_array_prealloc(&(wf->inel_ft), 2, 1+pixels_x/2, pixels_y, wf->phase.values.data + pixels_x*pixels_y);
  init_array(&(wf->intens.values), 1, pixels_x, pixels_y);
  wf->fftplan_wf_f = fftw_plan_dft_2d(wf->wf.size[2], wf->wf.size[1], 
                                      (fftw_complex*)wf->wf.data, (fftw_complex*)wf->wf_ft.data,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
  wf->fftplan_wf_b = fftw_plan_dft_2d(wf->wf.size[2], wf->wf.size[1], 
                                      (fftw_complex*)wf->wf_ft.data, (fftw_complex*)wf->wf.data,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
  wf->fftplan_inel_f = fftw_plan_dft_r2c_2d(wf->inel.size[2], wf->inel.size[1], 
                                            wf->inel.data, (fftw_complex*)wf->inel_ft.data,
                                            FFTW_ESTIMATE);
  wf->fftplan_inel_b = fftw_plan_dft_c2r_2d(wf->inel.size[2], wf->inel.size[1], 
                                            (fftw_complex*)wf->inel_ft.data, wf->inel.data,
                                            FFTW_ESTIMATE);
  wavefunction_set_inel_ctf(wf);
  write_log_comment("Wavefunction object initialized.\n\n");
  return 0;
}

/****************************************************************************/

void wavefunction_reset(wavefunction *wf){
  fftw_destroy_plan(wf->fftplan_wf_f);
  fftw_destroy_plan(wf->fftplan_wf_b);
  fftw_destroy_plan(wf->fftplan_inel_f);
  fftw_destroy_plan(wf->fftplan_inel_b);
  free_array(&(wf->phase.values));
  free_array(&(wf->wf));
  free_array(&(wf->wf_ft));
  free_array(&(wf->inel));
  free_array(&(wf->inel_ft));
  free_array(&(wf->inel_ft_copy));
  free_array(&(wf->pathlength_ice));
  free_array(&(wf->pathlength_supp));
  free_array(&(wf->intens.values));
  free_matrix(&(wf->inel_ctf));
}

/****************************************************************************/

int wavefunction_prop_el_opt(wavefunction *wf, long tilt, double z){
  double a, b, c, d, e, dxi, dxj, xi, xj, x, x2, x4, x6, env, chi, acc_en, acc_en_spr, df, 
    ccp, ap_ang, cs, k, B, ctfr, ctfi, re;
  const double epsilon = 0.5/ELEC_REST_ENERGY;
  long n, m, i, j;
  double *ft, *wd, *id;
  if(0 == wf->init){
    WARNING("Error in wavefunction_prop_el_opt: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((tilt < 0) || (tilt >= wf->opt->defocus.m)){
    WARNING("Error in wavefunction_prop_el_opt: Tilt number out of range.\n");
    return 1;
  }
  m = wf->wf_ft.size[1];
  n = wf->wf_ft.size[2];
  dxi = 2*M_PI/(m*wf->pixel_size);
  dxj = 2*M_PI/(n*wf->pixel_size);
  acc_en = electronbeam_get_acc_energy(wf->ed);
  k = wave_number(acc_en);
  acc_en_spr = electronbeam_get_energy_spread(wf->ed);
  df = z + optics_get_defocus(wf->opt, tilt); 
  cs = optics_get_cs(wf->opt);
  ccp = (1 + 2*epsilon*acc_en)/(acc_en*(1 + epsilon*acc_en)) * optics_get_cc(wf->opt);
  ap_ang = optics_get_cond_ap_angle(wf->opt);
  a = 0.5*df/k;
  b = -0.25*cs/(k*k*k);
  c = -0.25*df*df*ap_ang*ap_ang;
  d = (0.5*df*cs*ap_ang*ap_ang - 0.0625*acc_en_spr*acc_en_spr*ccp*ccp)/(k*k);
  e = -0.25*cs*cs*ap_ang*ap_ang/(k*k*k*k);
  B = min_d(M_PI/wf->pixel_size, 0.5 * k * optics_get_aperture(wf->opt) / optics_get_focal_length(wf->opt));
  fftw_execute(wf->fftplan_wf_f);
  ft = wf->wf_ft.data;
  for(j = 0; j < n; j++){
    xj = dxj * min_l(j, n - j);
    for(i = 0; i < m; i++){
      xi = dxi * min_l(i, m - i);
      x2 = xi*xi + xj*xj;
      x = sqrt(x2);
      if(x < B){
	x4 = x2*x2;
	x6 = x2*x4;
	env = 1.0/(n*m)*exp(c*x2 + d*x4 + e*x6);
	chi = a*x2 + b*x4;
	ctfr = env*cos(chi);
	ctfi = env*sin(chi);
      }
      else {
	ctfr = 0;
	ctfi = 0;
      }
      re = ctfr * ft[0] - ctfi * ft[1];
      ft[1] = ctfi * ft[0] + ctfr * ft[1];
      ft[0] = re;
      ft += 2;
    }
  }
  fftw_execute(wf->fftplan_wf_b);
  n = wf->wf.size[1]*wf->wf.size[2];
  wd = wf->wf.data;
  id = wf->intens.values.data;
  for (i = 0; i < n; i++){
    *id += wd[0]*wd[0] + wd[1]*wd[1];
    wd += 2;
    id++;
  }
  return 0;
}

/****************************************************************************/

int wavefunction_prop_el_vac(wavefunction *wf, double dist){
  double a, dxi, dxj, xi, xj, x2, env, chi, ctfr, ctfi, re, wavenum;
  long n, m, i, j;
  double *ft;
  if(0 == wf->init){
    WARNING("Error in wavefunction_prop_el_vac: Wavefunction object has not been initialized.\n");
    return 1;
  }
  m = wf->wf_ft.size[1];
  n = wf->wf_ft.size[2];
  dxi = 2*M_PI/(m*wf->pixel_size);
  dxj = 2*M_PI/(n*wf->pixel_size);
  wavenum = wave_number(electronbeam_get_acc_energy(wf->ed));
  a = -0.5*dist/wavenum;
  env = 1.0/(n*m);
  fftw_execute(wf->fftplan_wf_f);
  ft = wf->wf_ft.data;
  for(j = 0; j < n; j++){
    xj = dxj * min_l(j, n - j);
    for(i = 0; i < m; i++){
      xi = dxi * min_l(i, m - i);
      x2 = xi*xi + xj*xj;
      chi = a*x2;
      ctfr = env*cos(chi);
      ctfi = env*sin(chi);
      re = ctfr * ft[0] - ctfi * ft[1];
      ft[1] = ctfi * ft[0] + ctfr * ft[1];
      ft[0] = re;
      ft += 2;
    }
  }
  fftw_execute(wf->fftplan_wf_b);
  return 0;
}

/****************************************************************************/

int wavefunction_set_incoming(wavefunction *wf){
  long i, n;
  double *x;
  if(0 == wf->init){
    WARNING("Error in wavefunction_set_incoming: Wavefunction object has not been initialized.\n");
    return 1;
  }
  n = wf->wf.size[1] * wf->wf.size[2];
  x = wf->wf.data;
  for(i = 0; i < n; i++){
    x[0] = 1;
    x[1] = 0;
    x += 2;
  }
  fill_array(&(wf->intens.values), 0);
  return 0;
}

/****************************************************************************/

int wavefunction_adj_phase(wavefunction *wf){
  double *w, *ph, *inel, x, y, r, in;
  long n, i;
  if(0 == wf->init){
    WARNING("Error in wavefunction_adj_phase: Wavefunction object has not been initialized.\n");
    return 1;
  }
  n = wf->wf.size[1] * wf->wf.size[2];
  w = wf->wf.data;
  ph = wf->phase.values.data;
  inel = wf->inel.data;
  for(i = 0; i < n; i++){
    x = exp(-ph[1]);
    in = (1-x*x)*(w[0]*w[0]+w[1]*w[1]);
    y = x*sin(ph[0]);
    x *= cos(ph[0]);
    r = x*w[0] - y*w[1];
    w[1] = x*w[1] + y*w[0];
    w[0] = r;
    *inel = in;
    ph += 2;
    w += 2;
    inel++;
  }
  return 0;
}

/****************************************************************************/

int wavefunction_prop_inel(wavefunction *wf, long tilt, double z){
  long i, j, k;
  double df, dxi, dxj, xi, xj, s, ctf, c;
  double *id;
  if(0 == wf->init){
    WARNING("Error in wavefunction_prop_inel: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((tilt < 0) || (tilt >= wf->opt->defocus.m)){
    WARNING("Error in wavefunction_prop_inel: Tilt number out of range.\n");
    return 1;
  }
  df = fabs(z + optics_get_defocus(wf->opt, tilt));
  dxi = 2*M_PI*df/(wf->pixel_size * wf->inel.size[1]);
  dxj = 2*M_PI*df/(wf->pixel_size * wf->inel.size[2]);
  c = 1.0/(wf->inel.size[1]*wf->inel.size[2]);
  id = wf->inel_ft.data;
  fftw_execute(wf->fftplan_inel_f);
  for(j = 0; j < wf->inel_ft.size[2]; j++){
    xj = dxj * min_l(j, wf->inel.size[2]-j);
    for(i = 0; i < wf->inel_ft.size[1]; i++){
      xi = dxi * min_l(i, wf->inel.size[1]-i);
      s = sqrt(xi*xi + xj*xj)/wf->inel_ctf_dx;
      k = (long)floor(s);
      if (k < wf->inel_ctf.m-1) {
	s -= k;
        ctf = (1-s)*wf->inel_ctf.data[k] + s*wf->inel_ctf.data[k+1];
      }
      else {
        ctf = wf->inel_ctf.data[wf->inel_ctf.m-1];
      }
      ctf *= c;
      *id *= ctf;
      id++;
      *id *= ctf;
      id++;
    }
  }
  fftw_execute(wf->fftplan_inel_b);
  add_array(&(wf->inel), &(wf->intens.values), 1);
  return 0;
}

/****************************************************************************/

int wavefunction_apply_bg_blur(wavefunction *wf, long tilt){
  long i, j, k, l, n, m;
  double pl_min, pl_max, pl_step, pli, ctf, dxi, dxj, xi, xj, c, df, s, acc_en;
  double *pl1, *pl2, *ctfv, *id;
  if(0 == wf->init){
    WARNING("Error in wavefunction_apply_bg_blur: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if((tilt < 0) || (tilt >= wf->opt->defocus.m)){
    WARNING("Error in wavefunction_apply_inel_blur: Tilt number out of range.\n");
    return 1;
  }

  /* Find maximum and minimum pathlength */
  pl_min = HUGE_VAL; /* minimum pathlengt */
  pl_max = 0;        /* maximum pathlength of rays not passing through the support */
  pl1 = wf->pathlength_ice.data;
  pl2 = wf->pathlength_supp.data;
  n = wf->pathlength_ice.size[1] * wf->pathlength_ice.size[2];
  acc_en = electronbeam_get_acc_energy(wf->ed);
  c = sample_support_abs_pot(acc_en)/sample_ice_abs_pot(acc_en);
  for(i = 0; i < n; i++, pl1++, pl2++){
    *pl1 += c*(*pl2);
    pl_min = min_d(pl_min, *pl1);
    if(*pl2 < ONE_NANOMETER){
      pl_max = max_d(pl_max, *pl1);
    }
  }
  pl_max = max_d(pl_max, pl_min);
  if(pl_max < ONE_NANOMETER){
    return 0;
  }

  /* Determine number of thickness values to use */
  m = 20; /* maximum number */
  c = 10;
  if(pl_max - pl_min < ONE_NANOMETER){
    m = 0; /* all pathlengths can be considered equal */
  }
  else if(c*(pl_max-pl_min) < m*pl_min){
    m = (long)ceil(c*(pl_max - pl_min)/pl_min);
  }
  if(m > 0){
    pl_step = (pl_max - pl_min)/m;
  }
  else {
    pl_step = 1; /* just to avoid division by 0 */
  }

  /* Compute blurring for finite set of thicknesses and interpolate */
  df = fabs(optics_get_defocus(wf->opt, tilt));
  c = 1.0/(wf->inel.size[1]*wf->inel.size[2]);
  dxi = 2*M_PI*df/(wf->pixel_size * wf->inel.size[1]);
  dxj = 2*M_PI*df/(wf->pixel_size * wf->inel.size[2]);
  ctfv = wf->inel_ctf.data + wf->inel_ctf.m;
  copy_array(&(wf->intens.values), &(wf->inel));
  fftw_execute(wf->fftplan_inel_f);
  fill_array(&(wf->intens.values), 0);
  if(m > 0){
    /* fftw complex-to-real transform destroys its input, so a copy must be made */
    copy_array(&(wf->inel_ft), &(wf->inel_ft_copy));
  }
  for(l = 0, pli = pl_min; l <= m; l++, pli += pl_step){
    if(l > 0){
      copy_array(&(wf->inel_ft_copy), &(wf->inel_ft));
    }
    id = wf->inel_ft.data;
    for(j = 0; j < wf->inel_ft.size[2]; j++){
      xj = dxj * min_l(j, wf->inel.size[2]-j);
      for(i = 0; i < wf->inel_ft.size[1]; i++){
	xi = dxi * min_l(i, wf->inel.size[1]-i);
	s = sqrt(xi*xi + xj*xj)/wf->inel_ctf_dx;
	k = (long)floor(s);
	if (k < wf->inel_ctf.m-1) {
	  s -= k;
	  ctf = (1-s)*ctfv[k] + s*ctfv[k+1];
	}
	else {
	  ctf = ctfv[wf->inel_ctf.m-1];
	}
	ctf = c*exp(-ctf*pli);
	*id *= ctf;
	id++;
	*id *= ctf;
	id++;
      }
    }
    fftw_execute(wf->fftplan_inel_b);
    pl1 = wf->pathlength_ice.data;
    id = wf->inel.data;
    if(l < m){
      for(i = 0; i < n; i++, pl1++, id++){
	s = fabs(pli - *pl1);
	if(s < pl_step){
	  *id *= 1 - s/pl_step;
	}
	else {
	  *id = 0;
	}
      }
    }
    else {
      for(i = 0; i < n; i++, pl1++, id++){
        s = pli - *pl1;
	if(s >= 0){
          if(s < pl_step){
	    *id *= 1 - s/pl_step;
	  }
	  else {
	    *id = 0;
	  }
	}
	else {
	  *id *= exp(s*ctfv[0]);
	}
      }
    }
    add_array(&(wf->inel), &(wf->intens.values), 1);
  }
  return 0;
}

/****************************************************************************/

int wavefunction_set_inel_ctf(wavefunction *wf){
  double def_max = 0, dt, t_max, x_max, x, t, f, f0, acc_en;
  long nt, i, j;
  if(0 == wf->init){
    WARNING("Error in wavefunction_set_inel_ctf: Wavefunction object has not been initialized.\n");
    return 1;
  }
  for(i = 0; i < wf->opt->defocus.m; i++){
    def_max = max_d(def_max, fabs(optics_get_defocus(wf->opt, i)));
  }
  t_max = 0.5 * optics_get_aperture(wf->opt) / optics_get_focal_length(wf->opt);
  wf->inel_ctf_dx = 1/t_max;
  x_max = M_PI*(def_max+ONE_MICROMETER)/wf->pixel_size;
  dt = 1/x_max;
  nt = (long)floor(t_max/dt);
  init_matrix(&(wf->inel_ctf), (long)ceil(x_max/wf->inel_ctf_dx), 2);
  fill_matrix(&(wf->inel_ctf), 0);
  acc_en = electronbeam_get_acc_energy(wf->ed);
  f0 = cross_sec(acc_en, 6);
  for(i = 1; i <= nt; i++){
    t = i*dt;
    /* Differential inelastic cross section of carbon */
    f = 2*M_PI*t*dt/f0*diff_cross_sec(acc_en, 6, t);
    for(j = 0; j < wf->inel_ctf.m; j++){
      x = j*wf->inel_ctf_dx;
      wf->inel_ctf.data[j] += f*BESSEL0(x*t);
    }
    /* Differential cross section of water */
    f = 2*M_PI*t*dt*(diff_cross_sec(acc_en, 8, t) + 2*diff_cross_sec(acc_en, 1, t));
    for(j = 0; j < wf->inel_ctf.m; j++){
      x = j*wf->inel_ctf_dx;
      wf->inel_ctf.data[wf->inel_ctf.m + j] += f*BESSEL0(x*t);
    }
  }
  f0 = cross_sec(acc_en, 8) + 2*cross_sec(acc_en, 1);
  for(j = 0; j < wf->inel_ctf.m; j++){
    wf->inel_ctf.data[wf->inel_ctf.m + j] = ICE_DENS*(f0 - wf->inel_ctf.data[wf->inel_ctf.m + j]);
  }
  return 0;
}

/****************************************************************************/

int wavefunction_propagate(simulation *sim, wavefunction *wf, double slice_th, long tilt){
  int kmax, kmin, k, count = 0, j, c, first_hit = 1;
  int *kvec, *hit;
  long i, np;
  double pos[3];
  matrix pm;
  particle *p,*masked_particle;
  particleset *ps;
  sample *s = get_sample(sim, "");
  geometry *g = get_geometry(sim, "");
  if(s == NULL || sample_init(s, sim)){
    WARNING("Error in wavefunction_propagate: Sample component required.\n");
    return 1;
  }
  if(0 == wf->init){
    WARNING("Error in wavefunction_propagate: Wavefunction object has not been initialized.\n");
    return 1;
  }
  if(g == NULL || geometry_init(g, sim)){
    WARNING("Error in wavefunction_propagate: Geometry component required.\n");
    return 1;
  }
  for(j = 0; j < sim->num_particleset; j++){
    if(particleset_init(sim->particleset[j], sim)){
      WARNING("Error in wavefunction_propagate: Particleset component number %i could not be initialized.\n", j+1);
      return 1;
    }
  }
  if((tilt < 0) || (tilt >= g->data.m)){
    WARNING("Error in wavefunction_propagate: tilt number out of range.\n");
    return 1;
  }
  write_log_comment("\nComputing projection number %i.\n", tilt+1);
  wavefunction_set_incoming(wf); /* Set incoming wave to plane wave */
  init_matrix(&pm, 3, 3);
  /* Find range of slices */
  np = 0;
  for(j = 0; j < sim->num_particleset; j++){
    np += sim->particleset[j]->coordinates.m;
  }
  kmin = kmax = 0;
  if (np > 0) {
    kvec = malloc(np*sizeof(int));
    hit = malloc(np*sizeof(int));
    c = 0;
    for(j = 0; j < sim->num_particleset; j++){
      ps = sim->particleset[j];
      p = get_particle(sim, get_param_string(ps->param, PAR_PARTICLE_TYPE));
      if(p == NULL || particle_init(p, sim)){
	WARNING("Particle %s not found.\n", get_param_string(ps->param, PAR_PARTICLE_TYPE));
	return 1;
      }
      for(i = 0; i < sim->particleset[j]->coordinates.m; i++, c++){
	if(get_particle_geom(&pm, pos, sim->particleset[j], i, g, tilt)) return 1;
	if(particle_hits_wavefunction(p, wf, &pm, pos)){
	  k = round_to_int(pos[2]/slice_th);
          kvec[c] = k;
          hit[c] = 1;
	  if(first_hit || k > kmax){
	    kmax = k;
	  }
	  if(first_hit || k < kmin){
	    kmin = k;
	  }
          first_hit = 0;
	}
        else {
	  hit[c] = 0;
	}
      }
    }
    masked_particle	=	malloc(sizeof(particle));
    /* Project particles slice by slice */
    for(k = kmin; k <= kmax; k++){
      c = 0;
      fill_array(&(wf->phase.values), 0);
      for(j = 0; j < sim->num_particleset; j++){
        ps = sim->particleset[j];
        p = get_particle(sim, get_param_string(ps->param, PAR_PARTICLE_TYPE));
        if(p == NULL || particle_init(p, sim)){
	  WARNING("Particle %s not found.\n", get_param_string(ps->param, PAR_PARTICLE_TYPE));
	  return 1;
	}
    init_blank_similar_particle(p,masked_particle);
	for(i = 0; i < ps->coordinates.m; i++, c++){
	  if (hit[c] && k == kvec[c]){
	    if(get_particle_geom(&pm, pos, ps, i, g, tilt)) return 1;
	    pos[2] -= k*slice_th;
		mask_particle(p,masked_particle,s,ps,i);
	    if(particle_project(masked_particle, wf, &pm, pos)) return 1;
	    count++;
	  }
	}
	particle_reset(masked_particle);
      }
      wavefunction_adj_phase(wf);
      wavefunction_prop_inel(wf, tilt, k*slice_th);
      if(k < kmax){
	wavefunction_prop_el_vac(wf, slice_th);
      }
    }
    free(kvec);
    free(hit);
  }
  write_log_comment("Projected %i particles.\n", count);
  /* Propagate elastic wave to detector plane */
  if(wavefunction_prop_el_opt(wf, tilt, kmax*slice_th)) return 1;
  /* Correct for inelastic scattering in background */
  if(get_sample_geom(&pm, pos, s, g, tilt) || background_project(s, wf, &pm, pos)
     || wavefunction_apply_bg_blur(wf, tilt)) return 1;
  free_matrix(&pm);
  free(masked_particle);
  return 0;
}
