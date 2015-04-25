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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "array.h"
#include "matrix.h"
#include "structs.h"
#include "log.h"
#include "fftw3.h"
#include "sample.h"
#include "simulation.h"
#include "wavefunction.h"
#include "input.h"
#include "particle.h"
#include "geometry.h"
#include "particleset.h"
#include "misc.h"
#include "electronbeam.h"
#include "sample.h"
#include "multislice.h"

param_table *elec_spec_model_param_table(){
	param_table *pt = new_param_table(3, TYPE_ELECTRON_SPECIMEN_MODEL,"");
	//delete_param_table in delete_simulation.
	add_param_def(pt, PAR_INTERACTION_MODEL, "s", PAR_INTERACTION_MODEL__LINEAR);
	add_param_def(pt, PAR_SLICE_THICKNESS, "d", "10");
	add_param_def(pt, PAR_SLICE_WRITE_VOLUME, "b", NO_STRING);
	set_comp_descr(pt, "The electron_specimen_model determines which electron-specimen\
 interaction model will be used in the simulation.");
	set_param_descr(pt, PAR_INTERACTION_MODEL, "The interaction approximation model used\
 in the simulator. The value is a string from \"linear\"/\"projection\"/\"multislice\".");
	set_param_descr(pt, PAR_SLICE_THICKNESS, "The slice thickness in the multislice approximation\
 in nm. It is ignored if interaction_model is set to \"linear\".");
	set_param_descr(pt,PAR_SLICE_WRITE_VOLUME,"If it is set to yes, each slice is written on two MRC\
 files for real and imaginary parts of scattering potential of the slice. the name of files are\
 determined automatically.It is ignored if interaction_model is set to \"linear\".");
	return pt;
}

double trilinear_interpol(double f[][2][2], double x[], double y[], double z[], double pos[]){
	double x_d=(pos[0]-x[0])/(x[1]-x[0]);
	double y_d=(pos[1]-y[0])/(y[1]-y[0]);
	double z_d=(pos[2]-z[0])/(z[1]-z[0]);
	double c[2][2];
	double d[2];

	c[0][0]=f[0][0][0]*(1-x_d)+f[1][0][0]*x_d;
	c[1][0]=f[0][1][0]*(1-x_d)+f[1][1][0]*x_d;
	c[0][1]=f[0][0][1]*(1-x_d)+f[1][0][1]*x_d;
	c[1][1]=f[0][1][1]*(1-x_d)+f[1][1][1]*x_d;

	d[0]=c[0][0]*(1-y_d)+c[1][0]*y_d;
	d[1]=c[0][1]*(1-y_d)+c[1][1]*y_d;

	return (d[0]*(1-z_d)+d[1]*z_d);
}

int slice_background(array *slice_re,array *slice_im, sample *sam , wavefunction *wf, double pixel_size, double slice_pos_z){
	double m=slice_re->size[0] , n=slice_re->size[1] , o=slice_re->size[2];
	double mo=0.5*(m-1) , no=0.5*(n-1) , oo=0.5*(o-1) ;
	double zb = 0.5 * get_param_double(sam->param, PAR_THICKNESS_EDGE);
	double zc = 0.5 * get_param_double(sam->param, PAR_THICKNESS_CENTER);
	double D = get_param_double(sam->param, PAR_DIAMETER);
	double a=zb-zc;
	double l1,l2;
	double R=(a*a + D*D/4)/(2*a);
	double zo=R+zc;
	/*double acc_en = electronbeam_get_acc_energy(wf->ed);*/
	double p_re = sample_ice_pot();
	double p_im = sample_ice_abs_pot(electronbeam_get_acc_energy(wf->ed));
	long i,j,k;
	/*printf("\ntest: ice potential:%f+i%f",p_re,p_im);*/
	fill_array(slice_re, p_re);
	fill_array(slice_im, p_im);
	if (a!=0){
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				for (k=0;k<o;k++){
					l1=(i-mo)*pixel_size*(i-mo)*pixel_size + (j-no)*pixel_size*(j-no)*pixel_size + (slice_pos_z+(k-oo)*pixel_size-zo)*(slice_pos_z+(k-oo)*pixel_size-zo);
					l2=(i-mo)*pixel_size*(i-mo)*pixel_size + (j-no)*pixel_size*(j-no)*pixel_size + (slice_pos_z+(k-oo)*pixel_size+zo)*(slice_pos_z+(k-oo)*pixel_size+zo);
					if (l1<(R*R) || l2<(R*R)) {
						set_array_entry(slice_re,i,j,k, 0);
						set_array_entry(slice_im,i,j,k, 0);
					}
				}
			}
		}
	}
	return 0;
}

int voxels_inside_particle(array *in,array *out,double obj_pos_out[],double pixel_size_in,double pixel_size_out,long int marg[][2]){
	double c1=((double)out->size[0]-1)/2 , c2=((double)out->size[1]-1)/2 , c3=((double)out->size[2]-1)/2 ;
	double r=(double)sqrt((double)(in->size[0]*in->size[0]+in->size[1]*in->size[1]+in->size[2]*in->size[2]))*pixel_size_in/2;
	marg[0][0]=(long)floor((c1*pixel_size_out + obj_pos_out[0]-r)/pixel_size_out);
	marg[0][1]=(long)ceil((c1*pixel_size_out + obj_pos_out[0]+r)/pixel_size_out);
	marg[1][0]=(long)floor((c2*pixel_size_out + obj_pos_out[1]-r)/pixel_size_out);
	marg[1][1]=(long)ceil((c2*pixel_size_out + obj_pos_out[1]+r)/pixel_size_out);
	marg[2][0]=(long)floor((c3*pixel_size_out + obj_pos_out[2]-r)/pixel_size_out);
	marg[2][1]=(long)ceil((c3*pixel_size_out + obj_pos_out[2]+r)/pixel_size_out);

	//if ((r>abs(cio-obj_pos_out[0])) && (r>abs(cjo-obj_pos_out[1])) && (r>abs(cko-obj_pos_out[2]))) return 1;
	return 0;
}

double array_voxel_value(double cio,double cjo,double cko, array *p,double psize_p,double *pos,matrix *rm){
	matrix x,y,or;
	long i,j,k,m=p->size[0],n=p->size[1],o=p->size[2];
	double output;
	double f[2][2][2],X[2],Y[2],Z[2];
	double cube_pos[3];
	init_matrix(&x,1,3);
	init_matrix(&y,1,3);
	init_matrix(&or,1,3);

	set_matrix_entry(&x,0,0,cio-pos[0]);
	set_matrix_entry(&x,0,1,cjo-pos[1]);
	set_matrix_entry(&x,0,2,cko-pos[2]);

	matrix_mult(&x,rm,&y);

	set_matrix_entry(&or,0,0,((double)m-1)*psize_p/2);
	set_matrix_entry(&or,0,1,((double)n-1)*psize_p/2);
	set_matrix_entry(&or,0,2,((double)o-1)*psize_p/2);

	add_matrix(&y,&or,+1);

	i=(long)floor(get_matrix_entry(&or,0,0)/psize_p);
	j=(long)floor(get_matrix_entry(&or,0,1)/psize_p);
	k=(long)floor(get_matrix_entry(&or,0,2)/psize_p);

	f[0][0][0]=get_array_entry(p,i,j,k);
	f[1][0][0]=get_array_entry(p,i+1,j,k);
	f[0][1][0]=get_array_entry(p,i,j+1,k);
	f[1][1][0]=get_array_entry(p,i+1,j+1,k);
	f[0][0][1]=get_array_entry(p,i,j,k+1);
	f[1][0][1]=get_array_entry(p,i+1,j,k+1);
	f[0][1][1]=get_array_entry(p,i,j+1,k+1);
	f[1][1][1]=get_array_entry(p,i+1,j+1,k+1);

	X[0]=(double)i;
	X[1]=(double)i+1;

	Y[0]=(double)j;
	Y[1]=(double)j+1;

	Z[0]=(double)k;
	Z[1]=(double)k+1;

	cube_pos[0]=(get_matrix_entry(&or,0,0)/psize_p);
	cube_pos[1]=(get_matrix_entry(&or,0,1)/psize_p);
	cube_pos[2]=(get_matrix_entry(&or,0,2)/psize_p);


	output=trilinear_interpol(f, X, Y, Z, cube_pos);

	//-----Check for NaN entry.
	if (!(output==output)){
		WARNING("\nNaN was produced in the slice!");
		return(0);
	}

	free_matrix(&x);
	free_matrix(&y);
	free_matrix(&or);

	return output;
}

int interpol_array(array *in, array *out, matrix *rm, double obj_pos_out[] , double pixel_size_in, double pixel_size_out){
	long int i,j,k;
	long int marg[3][2];
	double cio,cjo,cko;
	double intr;
	long  int m=out->size[0],n=out->size[1],o=out->size[2];	
	voxels_inside_particle(in,out,obj_pos_out,pixel_size_in,pixel_size_out,marg);
	for (i=marg[0][0];i<=marg[0][1];i++){
		for (j=marg[1][0];j<=marg[1][1];j++){
			for (k=marg[2][0];k<=marg[2][1];k++){
				cio=(i-(m-1)*0.5)*pixel_size_out;
				cjo=(j-(n-1)*0.5)*pixel_size_out;
				cko=(k-(o-1)*0.5)*pixel_size_out;

				intr= array_voxel_value(cio,cjo,cko, in, pixel_size_in,obj_pos_out,rm);
				add_to_array_entry(out,i,j,k, intr);				
			}
		}
	}

	return 0;
}

int m_initializations(simulation *sim, long tilt, wavefunction *wf){
	int j;
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
	return 0;
}

int particle_in_slice(particle *p, double pos[3] , double low_z, double high_z ){
	double v_size = get_param_double(p->param,PAR_VOXEL_SIZE);
	long a=p->pot_re.size[0];
	long b=p->pot_re.size[1];
	long c=p->pot_re.size[2];
	double l=v_size*sqrt((double)(a*a+b*b+c*c))/2;
	if ( ((pos[2]+l) > low_z)  &&  ((pos[2]-l) < high_z) ) 
		return 1;
	else 
		return 0;
}

int multislice_volume_create (simulation *sim, wavefunction *wf, array *slice_re , array *slice_im,  double pixel_size, double low_z, double high_z , long tilt){
	particleset *ps ; 
	particle *p ;
	int i,j;
	double pos[3];
	double p_v_size;
	matrix *rm ;
	geometry *g = get_geometry(sim, "");
	rm = malloc(sizeof(matrix));	
	init_matrix(rm,3,3);
	for(j = 0; j < sim->num_particleset; j++){
		ps = sim->particleset[j];
		p = get_particle(sim, get_param_string(ps->param, PAR_PARTICLE_TYPE));
		if(p == NULL || particle_init(p, sim)){
			WARNING("\nParticle %s not found.\n", get_param_string(ps->param, PAR_PARTICLE_TYPE));
			return 1;
		}
		//printf("\ntest:max and min entry of the particle (%s):%f,%f",get_param_string(ps->param, PAR_PARTICLE_TYPE),max_array(&p->pot_re),min_array(&p->pot_re));
		for(i = 0; i < sim->particleset[j]->coordinates.m; i++){
			if(get_particle_geom(rm, pos, sim->particleset[j], i, g, tilt)) return 1;
			if (particle_in_slice(p ,pos , low_z, high_z)==1){
				//relative position of the particle :
				pos[2]=pos[2] - ((low_z+high_z))/2;
				p_v_size = get_param_double(p->param,PAR_VOXEL_SIZE);
				//---test
				//printf("\np_v_size=%f, particle size(%d,%d,%d), pos=(%f,%f,%f)",p_v_size,p->pot_re.size[0],p->pot_re.size[1],p->pot_re.size[2],pos[0],pos[1],pos[2] );
				//getch();

				if (interpol_array(&p->pot_re,slice_re,rm,pos,p_v_size,pixel_size) || interpol_array(&p->pot_im,slice_im,rm,pos,p_v_size,pixel_size)){
					WARNING("\nan Error occurred in interpolation of particle %ld from particleset %ld in slice %f-%f",i,j,low_z,high_z);
					return 1;
				}
			}
		}
	}
	free_matrix(rm);
	free(rm);
	return 0;
}

double elec_relative_mass(double acc_energy){
	return ELECTRON_MASS * ( 1 + acc_energy / (ELECTRON_MASS * SPEED_OF_LIGHT* SPEED_OF_LIGHT));
}

int proj_slice_to_phase(array *slice_re, array *slice_im, wavefunction *w, double pixel_size){
	long m=slice_re->size[0], n=slice_re->size[1], l=slice_re->size[2];
	long i,j,k;
	double *ph;
	double sigma , wavelength, elec_rel_mass;
	double max=-1000000,max_i=-1000000,min=1000000,min_i=1000000;
	fill_array(&(w->phase.values), 0);	
	for (k=0 ; k<l ; k++){
		ph = w->phase.values.data;
		/*phase has the type array(2, pixels_x, pixels_y) which has column-major structure. ->wavefunction_adj_phase */
		for (j=0 ; j<n ; j++){
			for (i=0 ; i<m ; i++){
				//real potential
				ph[0]+=get_array_entry(slice_re,i,j,k);
				//imaginary potential
				ph[1]+=get_array_entry(slice_im,i,j,k);			
				ph += 2;
			}
		}
	}
	// multiplication by sigma * pixel_size
	elec_rel_mass = elec_relative_mass (electronbeam_get_acc_energy(w->ed));
	wavelength = 2 * M_PI * ONE_NANOMETER / wave_number(electronbeam_get_acc_energy(w->ed));
	sigma = (2 * M_PI * elec_rel_mass * ELEMENTARY_CHARGE * wavelength) / PLANCKS_CONSTANT_SQ ;
	ph = w->phase.values.data;
	for (i=0 ; i < (m*n) ; i++){
		ph[0] = ph[0] * sigma * pixel_size  ;
		ph[1] = ph[1] * sigma * pixel_size  ;
		if (ph[0] > max) max = ph[0] ;
		if (ph[1] > max_i) max_i = ph[1] ;
		if (ph[0] < min) min = ph[0] ;
		if (ph[1] < min_i) min_i = ph[1] ;
		ph += 2;
	}
	//printf("\ntest:sigma:%f pixel_size:%f wavelenght:%f  max of phase values: %f elec_mass=%e\n",sigma,pixel_size,wavelength, max_array(&w->phase.values),elec_rel_mass);

	return 0;
}

int propagator_creator (array *propagator_ph, double wavelength, double pixel_size, double propagation_dist){
	long m = propagator_ph->size[1];
	long n = propagator_ph->size[2];
	double i,j,ip,jp;
	double m2=(double)m*0.5 , n2=(double)n*0.5 ;
	for (i=0 ; i < m ; i++){
		for (j=0 ; j<n ; j++){
			if (i>m2) ip=m-i; 
			else ip=i;
			if (j>n2) jp=n-j;
			else jp=j;
			set_array_entry(propagator_ph, 0,(long)i,(long)j, - M_PI * ( (ip/m) * (ip/m) + (jp/n) * (jp/n) )  * propagation_dist * wavelength / (pixel_size*pixel_size) ) ;

		}
	}
	//printf("\npropagator created. min entry:%e or %e\n",max_array(propagator_ph),min_array(propagator_ph) );
	return 0;
}
/*//------------------------------------
int Write_wave_on_file_mult(array *a){
	FILE *f=NULL;
	char *file_name="wave-array.txt";
	size_t len = 0;
	char buffer[50];
	long size= a->size[0] * a->size[1] * a->size[2];
	long int m=a->size[0], n=a->size[1], o=a->size[2];
	long i,j,k ;
	f=fopen(file_name,"w");
    if (f== NULL) {
    printf("Error in opening a file..", file_name);
    return 1;
    }
	for (k=0; k<o; k++){
		for (j=0; j<n; j++){
			for (i=0; i<m; i++){
				sprintf(buffer,"%f\n\0",get_array_entry(a,i,j,k));
                len = strlen(buffer);
                fwrite(buffer, sizeof(char), 8 , f);
				fwrite("\n", sizeof(char), 1 , f);
				//printf("%s",buffer);
			}
		}
	}
	fclose(f);
	printf("\nwriting file is completetd!\n");
	getch();
	return 0;
}
-------------------------------------------------------------------------------------
 */
int slice_propagation (wavefunction *wf,double propagation_dist){
	long i,n;
	double *ph,*prop_ph;
	array *propagator_ph;
	long m =  2*wf->wf.size[1] * wf->wf.size[2];

	propagator_ph = malloc(sizeof(array));
	n = wf->phase.values.size[1] * wf->phase.values.size[2];
	if ( init_array (propagator_ph, 1, wf->phase.values.size[1], wf->phase.values.size[2]) ) return 1;
	if (propagator_creator (propagator_ph, 2*M_PI / wave_number(electronbeam_get_acc_energy(wf->ed)), wf->pixel_size,  propagation_dist)) return 1;

	fftw_execute(wf->fftplan_wf_f);
	mult_array_const(&wf->wf_ft,2.0/(double)m); //compensating the scale factor of fftw

	//printf("\nmax entry of fft of wf before phase adj : %e or %e \n ",max_array(&wf->wf_ft), min_array(&wf->wf_ft) );	
	fill_array(&(wf->phase.values), 0);

	ph = wf->phase.values.data;
	prop_ph = propagator_ph->data;

	for (i=0 ; i<n ; ++i){
		ph[0]=prop_ph[i];
		ph[1]=0;
		ph+=2;
	}

	if ( wavefunction_adj_phase(wf) ) return 1;
	//printf("\nmax entry of fft of wf after phase adj:%e or %e \n",max_array(&wf->wf), min_array(&wf->wf) );
	fftw_execute(wf->fftplan_wf_b);
	free_array(propagator_ph);
	free(propagator_ph);
	//printf("\nmax entry of wf after propagation : %e\n",max_array(&wf->wf) );	
	return 0;
}

int write_slice_on_file(long i,double pixel_size ,array *slice_re, array *slice_im){
	char *map_axis_order="xyz";
	char *mrc = ".mrc";
	char index[4] ;
	char fn_re[20] = "slice_re_";
	char fn_im[20] = "slice_im_";
	sprintf (index,"%d",(int)i);
	strcat (fn_re , index);
	strcat (fn_re,mrc);
	strcat (fn_im , index);
	strcat (fn_im,mrc);
	write_array_float4b(slice_re, fn_re , mrc_header,map_axis_order,0,pixel_size);
	write_array_float4b(slice_im, fn_im , mrc_header,map_axis_order,0,pixel_size);
	write_log_comment("\nThe slice number %d is written on %s and %s.",i,fn_re,fn_im);
	return 0;
}

double get_max_wf(wavefunction *wf){
	array_data_type *wd;
	long i,n;
	double max=-1000,min=1000;
	double max_ph=-1000,min_ph=1000;
	double l,ph;
	wd = wf->wf.data;
	n = wf->wf.size[1]*wf->wf.size[2];
	for (i = 0; i < n; i++){
		l = sqrt(wd[0]*wd[0] + wd[1]*wd[1]);
		ph = wd[0]/wd[1];
		wd += 2;
		if (l>max) max = l;
		if (l<min) min = l;
		if (ph>max) max_ph = ph;
		if (ph<min) min_ph = ph;
	}
	printf("max amp and phase of wf = %f , %f,\nmin amp and phase of wf = %f , %f \n",max,max_ph,min,min_ph);
	return l;
}

int multislice(simulation *sim, wavefunction *wf,  long tilt){

	long slices , pixels_x=wf->wf.size[1] , pixels_y=wf->wf.size[2] , pixels_th , pixels_in_slice , i;
	double pixel_size = wf->pixel_size;
	array *slice_re , *slice_im;
	double low_z,high_z;
	double propagation_dist ;
	double slice_pos;
	//calculating the number of slices:
	double sample_th = get_param_double(sim->sample[0]->param , PAR_THICKNESS_EDGE);
	double m_slice_th=10;
	int multislice_mode=1; //multislice/projection

	m_initializations(sim, tilt,  wf);
	if (0==strcmp(get_param_string(sim->elec_spec_model_param,PAR_INTERACTION_MODEL),PAR_INTERACTION_MODEL__MULTISLICE)){
		m_slice_th=get_param_double(sim->elec_spec_model_param,PAR_SLICE_THICKNESS);
		write_log_comment("\nmultislice method started.\n");
	}
	else if (0==strcmp(get_param_string(sim->elec_spec_model_param,PAR_INTERACTION_MODEL),PAR_INTERACTION_MODEL__PROJECTION)){
		m_slice_th=get_param_double(sim->elec_spec_model_param,PAR_SLICE_THICKNESS);
		multislice_mode=0;
		write_log_comment("\nProjection method started.\n");
	}

	pixels_th = (long)ceil(sample_th/pixel_size);
	pixels_in_slice = (long)ceil(m_slice_th/pixel_size);

	//quantized slice thickness 
	m_slice_th = pixels_in_slice * pixel_size;
	slices = (long)ceil(sample_th/m_slice_th);	

	//printf("\ntest:-%s-,-%d-,-%f-",get_param_string(sim->elec_spec_model_param,PAR_INTERACTION_MODEL),get_param_boolean(sim->elec_spec_model_param,PAR_SLICE_WRITE_VOLUME),get_param_double(sim->elec_spec_model_param,PAR_SLICE_THICKNESS));
	write_log_comment("\nsample_thickness:%f\tnumber of slices:%d\tpixlesize:%f\tslice thickness:%f\n",sample_th,slices,pixel_size,m_slice_th);
	write_log_comment("slice dimensions: %d by %d by %d pixels\n",pixels_x,pixels_y,pixels_in_slice);
	write_log_comment("Propagation Started...\n");

	slice_re=malloc(sizeof(array));
	slice_im=malloc(sizeof(array));

	if(init_array(slice_re,pixels_x,pixels_y,pixels_in_slice)){
		WARNING ("\nRunning out of memory!\nError occurred while initilizing the slices for multislice methode!\ntry smaller slice thickness or size of detector(s).");
		return 1;
	}
	if (init_array(slice_im,pixels_x,pixels_y,pixels_in_slice)){
		WARNING ("\nRunning out of memory!\nError occurred while initilizing the slices for multislice methode!\ntry smaller slice thickness or size of detector(s).");
		return 1;
	}


	wavefunction_set_incoming(wf);
	for ( i=0 ; i < slices ; i++){
		low_z = -1 * sample_th/2 + i * m_slice_th ;
		high_z = low_z + m_slice_th;
		propagation_dist = high_z - max_d ( -sample_th/2 , low_z) ;
		slice_pos = high_z - 0.5*m_slice_th;

		if (slice_background(slice_re,slice_im,sim->sample[0],wf,pixel_size,slice_pos)) return 1;
		if (multislice_volume_create (sim, wf, slice_re , slice_im, pixel_size, low_z, high_z ,0)){
			WARNING("\nError occurred in multislice_volume_create at %dth slice.",i);
			return 1;
		}
		//printf("\ntest:slice %d created. max entry of slice %d : %f",i,i,max_array(slice_re));

		if (multislice_mode){
			if (slice_propagation(wf, propagation_dist)){
				WARNING("\nError occurred in slice_popagation at slice %d.",i);
				return 1;
			}
		}
		if (proj_slice_to_phase(slice_re, slice_im, wf, pixel_size)){
			WARNING("\nError occured in proj_slice_to_phase at slice %d ." ,i);
			return 1;
		}
		wavefunction_adj_phase(wf);
		//printf("\nSlice %d was successfully projected",i);
		// write slice on file.
		if(get_param_boolean(sim->elec_spec_model_param,PAR_SLICE_WRITE_VOLUME))
			if (write_slice_on_file(i,pixel_size ,slice_re, slice_im)){
				WARNING("\nError occurred in writing the slice number %d on file",i);
			}
		//log
		write_log_comment("Calculations on slice %d is finished.\n",i);
	}
	//printf("\nDone!\n");
	free_array(slice_re);
	free_array(slice_im);
	free(slice_re);
	free(slice_im);
	get_max_wf(wf);
	//propagate to detector through optics
	if (multislice_mode){
		if(wavefunction_prop_el_opt(wf, tilt, sample_th/2 )) return 1;
	}
	else
		if(wavefunction_prop_el_opt(wf, tilt, 0 )) return 1;
	printf("intensity range after optics : %f - %f\n",min_array(&wf->intens.values),max_array(&wf->intens.values));
	return 0;
}
