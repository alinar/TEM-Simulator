/*
 *
 * Copyright 2012 Ali Narangifard, Royal Institute of Technology and
 * Hans Rullgard, Stockholm University and 
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

#include "structs.h"
#include "array.h"

#ifndef MULTISLICE_HEADER
#define MULTISLICE_HEADER


/***********************************************************************
 * Function:  elec_spec_model_param_table
 * Purpose:   Initialize table of input parameters used for electron specimen
              interaction model.
 * Arguments: N/A
 * Return:    Pointer to parameter table.
 */
param_table *elec_spec_model_param_table();

/***********************************************************************
 * Function:  trilinear_interpol
 * Purpose:   calculates the 3-D interpolated value of a point.
 * Arguments: f - 2 by 2 by 2 array (***double). known values of corners of a cubic lattice.
			  x,y,z - arrays with length 2 determining positions of sides of lattice.
			  pos - array with length 3 determining position of the point interpolated.
 * Return:    interpolated value of the point.
 */
double trilinear_interpol(double f[][2][2], double x[], double y[], double z[], double pos[]);

/***********************************************************************
 * Function:  slice_background
 * Purpose:   adding the potential of ice to the slice
 * Arguments: slice_re - pointer to real part of slice
              slice_im - pointer to imaginary part of slice
			  sam - pointer to the sample
			  wavefunction - pointer to wavefunction (energy of wavefunction determines 
			  the imaginary part of the ice's potential
			  pixel_size
			  slice_pos_z - z value of centeral position of the slice.

 * Return:    0 if succeed. 
 */
int slice_background(array *slice_re,array *slice_im, sample *sam , wavefunction *wf, double pixel_size, double slice_pos_z);

/***********************************************************************
 * Function:  voxels_inside_particle
 * Purpose:   determine index margin of voxels belong to a particle inside a slice.
 * Arguments: in - pointer to the particle array.
              out - pointer to the slice array.
			  obj_pos_out - *double to position of the particle.
			  pixel_size_in - particle's array pixel size
			  pixel_size_out - slice's array pixel size
			  marg - **long int - 3 by 2- returned voxel indice of lower and upper margin of the
			  particle inside slice's grid. 
 * Return:    0 if succeed. 
 */
int voxels_inside_particle(array *in,array *out,double obj_pos_out[],double pixel_size_in,double pixel_size_out,long int marg[][2]);

/***********************************************************************
 * Function:  array_voxel_value
 * Purpose:   determine the value of a voxel of the slice which is inside a particle
 * Arguments: cio-cjo-cko x,y,z posiiton of the voxel.
              p - pointer to particle's array
			  psize_p - pixel size of particle
			  pos - pointer to position of the particle
			  rm - pointer to the rotation matrix of the particle.
 * Return:    determined value of the particle. 
 */

double array_voxel_value(double cio,double cjo,double cko, array *p,double psize_p,double *pos,matrix *rm);

/***********************************************************************
 * Function:  interpol_array
 * Purpose:   interpolation of slice voxels which are inside a particle.
 * Arguments: in - pointer to particle's array
			  out - pointer to slice's array
			  rm - pointer to rotation matrix of the particle
			  obj_pos_out - poisition of the particle inide the slice
			  pixel_size_in
			  pixel_size_out
 * Return:    0 if succeed.
 */

int interpol_array(array *in, array *out, matrix *rm, double obj_pos_out[] , double pixel_size_in, double pixel_size_out);

/***********************************************************************
 * Function:  m_initializations
 * Purpose:   initializations needed.
 * Arguments: sim - pointer to the simulation
			  tilt - tilt's number
			  wf - pointer to wavefunction
 * Return:    0 if succeed.
 */

int m_initializations(simulation *sim, long tilt, wavefunction *wf);

/***********************************************************************
 * Function:  particle_in_slice
 * Purpose:   check if a particle is located inside a slice 
 * Arguments: p - pointer to particle
			  pos - pointer to position of the particle
			  low_z - lower z value of the slice
			  high_z - upper z value of the slice
 * Return:    1 if the particle is located inside of the slice
			  0 otherwise
 */

int particle_in_slice(particle *p, double pos[3] , double low_z, double high_z );

/***********************************************************************
 * Function:  multislice_volume_create
 * Purpose:   creating a slice of specimen
 * Arguments: sim - pointer to simulation
			  wf - pointer to wavefuntion
			  slice_re - pointer to real part of slice
			  slice_im - pointer to imaginary part of slice
			  pixel_size
			  low_z - lower z value of the slice
			  high_z - upper z value of the slice
			  tilt - tilt number
 * Return:    1 if the particle is located inside of the slice
			  0 otherwise
 */

int multislice_volume_create (simulation *sim, wavefunction *wf, array *slice_re , array *slice_im,  double pixel_size, double low_z, double high_z , long tilt);

/***********************************************************************
 * Function:  elec_relative_mass
 * Purpose:   calculates electron's relative mass
 * Arguments: acc_energy - electron's energy
 * Return:    relative mass
 */

double elec_relative_mass(double acc_energy);

/***********************************************************************
 * Function:  proj_slice_to_phase
 * Purpose:   calculates the slice's project along z axis and phase-shift 
			  the wavefunction respecting to the projection.
 * Arguments: slice_re - pointer to real part of slice
			  slice_im - pointer to imaginary part of slice
			  w - pointer to wavefunction
			  pixel_size
 * Return:    zero on success
 */

int proj_slice_to_phase(array *slice_re, array *slice_im, wavefunction *w, double pixel_size);

/***********************************************************************
 * Function:  propagator_creator
 * Purpose:   returns the phase of sampled of the fourier transform of 
			  Fresnel propagator and rearrange the output to conform 2D 
			  DFT pattern.
			  (low freq in the corners, high freq in the middle)
 * Arguments: propagator_ph - pointer to the output phase of the propagator
			  wavelength
			  pixel_size
			  propagation_dist - propagation disctance
 * Return:    zero on success
 */

int propagator_creator (array *propagator_ph, double wavelength, double pixel_size, double propagation_dist);

/***********************************************************************
 * Function:  slice_propagation
 * Purpose:   applies Fresnel propagation into the wavefunction
 * Arguments: wf - pointer to wavefunction
			  propagation_dist - distance of propagation
 * Return:    zero on success
 */

int slice_propagation (wavefunction *wf,double propagation_dist);

/***********************************************************************
 * Function:  write_slice_on_file
 * Purpose:   saves the slice on an MRC file
 * Arguments: i - index of the slice to be used in the name of file
			  pixel_size
			  slice_re - pointer to real part of slice
			  slice_im - pointer to imaginary part of slice
 * Return:    zero on success
 */

int write_slice_on_file(long i,double pixel_size ,array *slice_re, array *slice_im);

/***********************************************************************
 * Function:  multislice
 * Purpose:   performs multislice/projection method in electron-specimen
		      interaction and then applies optic effects and calculates
			  outgoing wavefunction.
 * Arguments: sim - pointer to the simulation
			  wf - pointer to wavefunction
			  tilt - tilt number
 * Return:    zero on success
 */

int multislice(simulation *sim, wavefunction *wf, long tilt);

double get_max_wf(wavefunction *wf);

#endif
