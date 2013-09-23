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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pdb.h"
#include "electronbeam.h"
#include "input.h"
#include "log.h"
#include "misc.h"
#include "sample.h"


/****************************************************************************/

#define PDB_LINE_LENGTH 100
#define NUM_ELEM 98

typedef struct {
  char atom_name[5];
  char element_sym[3];
  int atomic_num;
  char residue_name[4];
  double coord[3];
  double occ;
  double temp_fact;
} pdbrecord;

int add_atomic_pot(array *vol, double xc, double yc, double zc, int atomic_num, double temp_fact, double voxel_size, double num);

int add_atomic_abs_pot(array *vol, double xc, double yc, double zc, int atomic_num, double temp_fact, 
		       double voxel_size, double num, double acc_en);

int add_gaussians(array *vol, double xc, double yc, double zc, double *a, double *b, int n, double voxel_size);

int pdb_find_box(double box[6], const char *fn, array *transf, double padding);

int parse_pdb_line(char *line, pdbrecord *p);

int read_pdb_transf_file(array *a, particle *p);

void transf_coord(double coord[3], pdbrecord *p, array *transf, long k);

int atomic_number_lookup(char *s);

void atom_name_lookup(int number, char name[3]);

int nhydro(char *seq, char *pos);

/****************************************************************************/

int get_pot_from_pdb(particle *p, electronbeam *ed){
  double box[6] = {0,0,0,0,0,0}, coord[3];
  const double padding = 5;
  double vox_sz_angst;
  const char *pdb_file;
  int num_atoms[NUM_ELEM+1];
  int nh = 0;
  long n, i, nx, ny, nz, j;
  double acc_en = 0, solvent_pot;
  double *vd;
  boolean add_hydrogen;
  char line[PDB_LINE_LENGTH+1], name[3];
  pdbrecord pdbrec;
  array transf;
  FILE *fp;
  if(0 == p->init){
    WARNING("Error in get_pot_from_pdb: Particle object has not been initialized.\n");
    return 1;
  }
  vox_sz_angst = get_param_double(p->param, PAR_VOXEL_SIZE)/ONE_ANGSTROM;
  pdb_file = get_param_string(p->param, PAR_PDB_FILE_IN);
  write_log_comment("\nGenerating potential from pdb file %s.\n", pdb_file);
  add_hydrogen = get_param_boolean(p->param, PAR_ADD_HYDROGEN);
  for(i = 0; i <= NUM_ELEM; i++){
    num_atoms[i] = 0;
  }
  read_pdb_transf_file(&transf, p);
  if(pdb_find_box(box, pdb_file, &transf, padding)) return 1;
  if(param_isset(p->param, PAR_NX)){
    nx = get_param_long(p->param, PAR_NX);
  }
  else {
    nx = (long)ceil((box[1] - box[0])/vox_sz_angst);
    set_param_long(p->param, PAR_NX, nx);
  }
  if(param_isset(p->param, PAR_NY)){
    ny = get_param_long(p->param, PAR_NY);
  }
  else {
    ny = (long)ceil((box[3] - box[2])/vox_sz_angst);
    set_param_long(p->param, PAR_NY, ny);
  }
  if(param_isset(p->param, PAR_NZ)){
    nz = get_param_long(p->param, PAR_NZ);
  }
  else {
    nz = (long)ceil((box[5] - box[4])/vox_sz_angst);
    set_param_long(p->param, PAR_NZ, nz);
  }
  box[0] += 0.5*(box[1] - box[0] - (nx - 1) * vox_sz_angst);
  box[2] += 0.5*(box[3] - box[2] - (ny - 1) * vox_sz_angst);
  box[4] += 0.5*(box[5] - box[4] - (nz - 1) * vox_sz_angst);
  if(p->use_imag_pot){
    if((NULL == ed) || require_param(ed->param, PAR_ACC_VOLTAGE)){
      WARNING("Can't determine acceleration energy needed to compute absorption potential");
      return 1;
    }
    acc_en = electronbeam_get_acc_energy(ed);
    init_array(&(p->pot_im), nx, ny, nz);
    fill_array(&(p->pot_im), 0);
  }
  init_array(&(p->pot_re), nx, ny, nz);
  fill_array(&(p->pot_re), 0);
  FOPEN(fp, pdb_file, "r");
  if(fp == NULL){
    WARNING("Could not open file %s for reading.\n", pdb_file);
    return 1;
  }
  while(!feof(fp)){
    if(NULL == fgets(line, PDB_LINE_LENGTH+1, fp)){
      break;
    }
    if(parse_pdb_line(line, &pdbrec)){
      if((pdbrec.atomic_num >= 0)&&(pdbrec.atomic_num <= NUM_ELEM)){
	num_atoms[pdbrec.atomic_num]++;
      }
      if(0 != pdbrec.atomic_num){
        if(add_hydrogen) {
          nh = nhydro(pdbrec.residue_name, pdbrec.atom_name);
          num_atoms[1] += nh;
	}
	else {
          nh = 0;
	}
	for(j = 0; j < transf.size[2]; j++){
	  transf_coord(coord, &pdbrec, &transf, j);
	  coord[0] -= box[0];
	  coord[1] -= box[2];
	  coord[2] -= box[4];
	  add_atomic_pot(&(p->pot_re), coord[0], coord[1], coord[2], pdbrec.atomic_num, pdbrec.temp_fact, vox_sz_angst, pdbrec.occ);
	  if(nh){
	    add_atomic_pot(&(p->pot_re), coord[0], coord[1], coord[2], 1, pdbrec.temp_fact, vox_sz_angst, nh*pdbrec.occ);
	  }
	  if(p->use_imag_pot){
	    add_atomic_abs_pot(&(p->pot_im), coord[0], coord[1], coord[2], pdbrec.atomic_num, 
			       pdbrec.temp_fact, vox_sz_angst, pdbrec.occ, acc_en);
	    if(nh){
	      add_atomic_abs_pot(&(p->pot_im), coord[0], coord[1], coord[2], 1, 
				 pdbrec.temp_fact, vox_sz_angst, nh*pdbrec.occ, acc_en);
	    }
	  }
	}
      }
    }
  }
  fclose(fp);
  n = nele_array(&(p->pot_re));
  vd = p->pot_re.data;
  solvent_pot = sample_ice_pot();
  for(i = 0; i < n; i++){
    *vd = max_d(solvent_pot, *vd);
    vd++;
  }
  if(p->use_imag_pot){
    vd = p->pot_im.data;
    solvent_pot = sample_ice_abs_pot(acc_en);
    for(i = 0; i < n; i++){
      *vd = max_d(solvent_pot, *vd);
      vd++;
    }
  }
  write_log_comment("Size of particle map: %i x %i x %i\n", (int)nx, (int)ny, (int)nz);
  write_log_comment("PDB coordinates of map center: (%f, %f, %f)\n", (float)(box[0] + 0.5*(nx-1)*vox_sz_angst), 
		    (float)(box[2] + 0.5*(ny-1)*vox_sz_angst), (float)(box[4] + 0.5*(nz-1)*vox_sz_angst));
  write_log_comment("Number of atoms used to compute potential:\n");
  for(i = 1; i <= NUM_ELEM; i++){
    if(num_atoms[i] > 0){
      atom_name_lookup(i, name);
      write_log_comment("%s : %i\n", name, num_atoms[i]);
    }
  }
  if(num_atoms[0] > 0){
    write_log_comment("Unknown entries in PDB file: %i\n", num_atoms[0]);
  }
  write_log_comment("Number of transformed copies of PDB model: %i\n\n", (int)transf.size[2]);
  free_array(&transf);
  return 0;
}

/****************************************************************************/

int add_atomic_pot(array *vol, double xc, double yc, double zc, int atomic_num, double temp_fact, double voxel_size, double num){
  int i;
  static double a[NUM_SCATTERING_PARAMETERS], b[NUM_SCATTERING_PARAMETERS];
  if(get_scattering_parameters(a, b, atomic_num)) return 1;
  for(i = 0; i < NUM_SCATTERING_PARAMETERS; i++){
    a[i] *= num;
    b[i] += temp_fact;
    b[i] = max_d(b[i], 200*voxel_size*voxel_size);
    a[i] *= 4*sqrt(M_PI)*pow(b[i], -1.5)*PLANCKS_CONSTANT_SQ/(ELECTRON_MASS*ELEMENTARY_CHARGE)/(ONE_ANGSTROM_SQ*POTENTIAL_UNIT);
    b[i] = 4*M_PI*M_PI/b[i];

  }
  return add_gaussians(vol, xc, yc, zc, a, b, NUM_SCATTERING_PARAMETERS, voxel_size);
}

/****************************************************************************/

int add_atomic_abs_pot(array *vol, double xc, double yc, double zc, int atomic_num, double temp_fact, 
		       double voxel_size, double num, double acc_en){
  double wn = wave_number(acc_en);
  double b = 4*M_PI*M_PI/max_d(temp_fact + 20, 200*voxel_size*voxel_size);
  double a = PLANCKS_CONSTANT_SQ/(8*M_PI*M_PI*ELECTRON_MASS*ELEMENTARY_CHARGE*ONE_ANGSTROM_CU*POTENTIAL_UNIT)
    *num*wn/(1 + acc_en/ELEC_REST_ENERGY)*pow(b/M_PI, 1.5)*cross_sec(acc_en, atomic_num);
  return add_gaussians(vol, xc, yc, zc, &a, &b, 1, voxel_size);
}

/****************************************************************************/

int add_gaussians(array *vol, double xc, double yc, double zc, double *a, double *b, int n, double voxel_size){
  long i, j, k, l, imin, imax, jmin, jmax, kmin, kmax;
  double x, y, z, u2, v2, z2, r, r2 = 0, s, t;
  double *vd;
  for(l = 0; l < n; l++){
    b[l] *= voxel_size*voxel_size;
    r2 = max_d(r2, 10/b[l]);
  }
  r = sqrt(r2);
  xc /= voxel_size;
  yc /= voxel_size;
  zc /= voxel_size;
  kmin = max_l(0, (long)ceil(zc - r));
  kmax = min_l(vol->size[2]-1, (long)floor(zc + r));
  for(k = kmin; k <= kmax; k++){
    z = zc - k;
    z2 = z*z;
    s = sqrt(r2 - z2);
    jmin = max_l(0, (long)ceil(yc - s));
    jmax = min_l(vol->size[1]-1, (long)floor(yc + s));
    for(j = jmin; j <= jmax; j++){
      y = yc - j;
      v2 = z2 + y*y;
      t = sqrt(r2-v2);
      imin = max_l(0, (long)ceil(xc - t));
      imax = min_l(vol->size[0]-1, (long)floor(xc + t));
      vd = get_array_entry_ptr(vol, imin, j, k);
      for(i = imin; i <= imax; i++){
	x = xc - i;
	u2 = v2 + x*x;
	for(l = 0; l < n; l++){
	  *vd += a[l] * exp(-b[l]*u2);
	}
	vd++;
      }
    }
  }
  return 0;
}

/****************************************************************************/

/* The following numerical values are from L.-M. Peng, G. Ren, S. L. Dudarev, 
 * M. J. Whelan, Robust Parameterization of Elastic and Absorptive Electron 
 * Atomic Scattering Factors, Acta Cryst. (1996). A52, 257-276, Table 1.
 */

int get_scattering_parameters(double a[NUM_SCATTERING_PARAMETERS], double b[NUM_SCATTERING_PARAMETERS], int atomic_num){
  const static double scfact[NUM_ELEM][2*NUM_SCATTERING_PARAMETERS] =
    {{0.0349,0.1201,0.1970, 0.0573, 0.1195,0.5347,3.5867,12.3471, 18.9525, 38.6269},
     {0.0317,0.0838,0.1526, 0.1334, 0.0164,0.2507,1.4751, 4.4938, 12.6646, 31.1653},
     {0.0750,0.2249,0.5548, 1.4954, 0.9354,0.3864,2.9383,15.3829, 53.5545,138.7337},
     {0.0780,0.2210,0.6740, 1.3867, 0.6925,0.3131,2.2381,10.1517, 30.9061, 78.3273},
     {0.0909,0.2551,0.7738, 1.2136, 0.4606,0.2995,2.1155, 8.3816, 24.1292, 63.1314},
     {0.0893,0.2563,0.7570, 1.0487, 0.3575,0.2465,1.7100, 6.4094, 18.6113, 50.2523},
     {0.1022,0.3219,0.7982, 0.8197, 0.1715,0.2451,1.7481, 6.1925, 17.3894, 48.1431},
     {0.0974,0.2921,0.6910, 0.6990, 0.2039,0.2067,1.3815, 4.6943, 12.7105, 32.4726},
     {0.1083,0.3175,0.6487, 0.5846, 0.1421,0.2057,1.3439, 4.2788, 11.3932, 28.7881},
     {0.1269,0.3535,0.5582, 0.4674, 0.1460,0.2200,1.3779, 4.0203,  9.4934, 23.1278},
     {0.2142,0.6853,0.7692, 1.6589, 1.4482,0.3334,2.3446,10.0830, 48.3037,138.2700},
     {0.2314,0.6866,0.9677, 2.1882, 1.1339,0.3278,2.2720,10.9241, 39.2898,101.9748},
     {0.2390,0.6573,1.2011, 2.5586, 1.2312,0.3138,2.1063,10.4163, 34.4552, 98.5344},
     {0.2519,0.6372,1.3795, 2.5082, 1.0500,0.3075,2.0174, 9.6746, 29.3744, 80.4732},
     {0.2548,0.6106,1.4541, 2.3204, 0.8477,0.2908,1.8740, 8.5176, 24.3434, 63.2996},
     {0.2497,0.5628,1.3899, 2.1865, 0.7715,0.2681,1.6711, 7.0267, 19.5377, 50.3888},
     {0.2443,0.5397,1.3919, 2.0197, 0.6621,0.2468,1.5242, 6.1537, 16.6687, 42.3086},
     {0.2385,0.5017,1.3428, 1.8899, 0.6079,0.2289,1.3694, 5.2561, 14.0928, 35.5361},
     {0.4115,-1.4031,2.2784,2.6742, 2.2162,0.3703,3.3874,13.1029, 68.9592,194.4329},
     {0.4054,1.3880,2.1602, 3.7532, 2.2063,0.3499,3.0991,11.9608, 53.9353,142.3892},
     {0.3787,1.2181,2.0594, 3.2618, 2.3870,0.3133,2.5856, 9.5813, 41.7688,116.7282},
     {0.3825,1.2598,2.0008, 3.0617, 2.0694,0.3040,2.4863, 9.2783, 39.0751,109.4583},
     {0.3876,1.2750,1.9109, 2.8314, 1.8979,0.2967,2.3780, 8.7981, 35.9528,101.7201},
     {0.4046,1.3696,1.8941, 2.0800, 1.2196,0.2986,2.3958, 9.1406, 37.4701,113.7121},
     {0.3796,1.2094,1.7815, 2.5420, 1.5937,0.2699,2.0455, 7.4726, 31.0604, 91.5622},
     {0.3946,1.2725,1.7031, 2.3140, 1.4795,0.2717,2.0443, 7.6007, 29.9714, 86.2265},
     {0.4118,1.3161,1.6493, 2.1930, 1.2830,0.2742,2.0372, 7.7205, 29.9680, 84.9383},
     {0.3860,1.1765,1.5451, 2.0730, 1.3814,0.2478,1.7660, 6.3107, 25.2204, 74.3146},
     {0.4314,1.3208,1.5236, 1.4671, 0.8562,0.2694,1.9223, 7.3474, 28.9892, 90.6246},
     {0.4288,1.2646,1.4472, 1.8294, 1.0934,0.2593,1.7998, 6.7500, 25.5860, 73.5284},
     {0.4818,1.4032,1.6561, 2.4605, 1.1054,0.2825,1.9785, 8.7546, 32.5238, 98.5523},
     {0.4655,1.3014,1.6088, 2.6998, 1.3003,0.2647,1.7926, 7.6071, 26.5541, 77.5238},
     {0.4517,1.2229,1.5852, 2.7958, 1.2638,0.2493,1.6436, 6.8154, 22.3681, 62.0390},
     {0.4477,1.1678,1.5843, 2.8087, 1.1956,0.2405,1.5442, 6.3231, 19.4610, 52.0233},
     {0.4798,1.1948,1.8695, 2.6953, 0.8203,0.2504,1.5963, 6.9653, 19.8492, 50.3233},
     {0.4546,1.0993,1.7696, 2.7068, 0.8672,0.2309,1.4279, 5.9449, 16.6752, 42.2243},
     {1.0160,2.8528,3.5466,-7.7804,12.1148,0.4853,5.0925,25.7851,130.4515,138.6775},
     {0.6703,1.4926,3.3368, 4.4600, 3.1501,0.3190,2.2287,10.3504, 52.3291,151.2216},
     {0.6894,1.5474,3.2450, 4.2126, 2.9764,0.3189,2.2904,10.0062, 44.0771,125.0120},
     {0.6719,1.4684,3.1668, 3.9557, 2.8920,0.3036,2.1249, 8.9236, 36.8458,108.2049},
     {0.6123,1.2677,3.0348, 3.3841, 2.3683,0.2709,1.7683, 7.2489, 27.9465, 98.5624},
     {0.6773,1.4798,3.1788, 3.0824, 1.8384,0.2920,2.0606, 8.1129, 30.5336,100.0658},
     {0.7082,1.6392,3.1993, 3.4327, 1.8711,0.2976,2.2106, 8.5246, 33.1456, 96.6377},
     {0.6735,1.4934,3.0966, 2.7254, 1.5597,0.2773,1.9716, 7.3249, 26.6891, 90.5581},
     {0.6413,1.3690,2.9854, 2.6952, 1.5433,0.2580,1.7721, 6.3854, 23.2549, 85.1517},
     {0.5904,1.1775,2.6519, 2.2875, 0.8689,0.2324,1.5019, 5.1591, 15.5428, 46.8213},
     {0.6377,1.3790,2.8294, 2.3631, 1.4553,0.2466,1.6974, 5.7656, 20.0943, 76.7372},
     {0.6364,1.4247,2.7802, 2.5973, 1.7886,0.2407,1.6823, 5.6588, 20.7219, 69.1109},
     {0.6768,1.6589,2.7740, 3.1835, 2.1326,0.2522,1.8545, 6.2936, 25.1457, 84.5448},
     {0.7224,1.9610,2.7161, 3.5603, 1.8972,0.2651,2.0604, 7.3011, 27.5493, 81.3349},
     {0.7106,1.9247,2.6149, 3.8322, 1.8899,0.2562,1.9646, 6.8852, 24.7648, 68.9168},
     {0.6947,1.8690,2.5356, 4.0013, 1.8955,0.2459,1.8542, 6.4411, 22.1730, 59.2206},
     {0.7047,1.9484,2.5940, 4.1526, 1.5057,0.2455,1.8638, 6.7639, 21.8007, 56.4395},
     {0.6737,1.7908,2.4129, 4.2100, 1.7058,0.2305,1.6890, 5.8218, 18.3928, 47.2496},
     {1.2704,3.8018,5.6618, 0.9205, 4.8105,0.4356,4.2058,23.4342,136.7783,171.7561},
     {0.9049,2.6076,4.8498, 5.1603, 4.7388,0.3066,2.4363,12.1821, 54.6135,161.9978},
     {0.8405,2.3863,4.6139, 5.1514, 4.7949,0.2791,2.1410,10.3400, 41.9148,132.0204},
     {0.8551,2.3915,4.5772, 5.0278, 4.5118,0.2805,2.1200,10.1808, 42.0633,130.9893},
     {0.9096,2.5313,4.5266, 4.6376, 4.3690,0.2939,2.2471,10.8266, 48.8842,147.6020},
     {0.8807,2.4183,4.4448, 4.6858, 4.1725,0.2802,2.0836,10.0357, 47.4506,146.9976},
     {0.9471,2.5463,4.3523, 4.4789, 3.9080,0.2977,2.2276,10.5762, 49.3619,145.3580},
     {0.9699,2.5837,4.2778, 4.4575, 3.5985,0.3003,2.2447,10.6487, 50.7994,146.4179},
     {0.8694,2.2413,3.9196, 3.9694, 4.5498,0.2653,1.8590, 8.3998, 36.7397,125.7089},
     {0.9673,2.4702,4.1148, 4.4972, 3.2099,0.2909,2.1014, 9.7067, 43.4270,125.9474},
     {0.9325,2.3673,3.8791, 3.9674, 3.7996,0.2761,1.9511, 8.9296, 41.5937,131.0122},
     {0.9505,2.3705,3.8218, 4.0471, 3.4451,0.2773,1.9469, 8.8862, 43.0938,133.1396},
     {0.9248,2.2428,3.6182, 3.7910, 3.7912,0.2660,1.8183, 7.9655, 33.1129,101.8139},
     {1.0373,2.4824,3.6558, 3.8925, 3.0056,0.2944,2.0797, 9.4156, 45.8056,132.7720},
     {1.0075,2.3787,3.5440, 3.6932, 3.1759,0.2816,1.9486, 8.7162, 41.8420,125.0320},
     {1.0347,2.3911,3.4619, 3.6556, 3.0052,0.2855,1.9679, 8.7619, 42.3304,125.6499},
     {0.9927,2.2436,3.3554, 3.7813, 3.0994,0.2701,1.8073, 7.8112, 34.4849,103.3526},
     {1.0295,2.2911,3.4110, 3.9497, 2.4925,0.2761,1.8625, 8.0961, 34.2712, 98.5295},
     {1.0190,2.2291,3.4097, 3.9252, 2.2679,0.2694,1.7962, 7.6944, 31.0942, 91.1089},
     {0.9853,2.1167,3.3570, 3.7981, 2.2798,0.2569,1.6745, 7.0098, 26.9234, 81.3910},
     {0.9914,2.0858,3.4531, 3.8812, 1.8526,0.2548,1.6518, 6.8845, 26.7234, 81.7215},
     {0.9813,2.0322,3.3665, 3.6235, 1.9741,0.2487,1.5973, 6.4737, 23.2817, 70.9254},
     {1.0194,2.0645,3.4425, 3.4914, 1.6976,0.2554,1.6475, 6.5966, 23.2269, 70.0272},
     {0.9148,1.8096,3.2134, 3.2953, 1.5754,0.2263,1.3813, 5.3243, 17.5987, 60.0171},
     {0.9674,1.8916,3.3993, 3.0524, 1.2607,0.2358,1.4712, 5.6758, 18.7119, 61.5286},
     {1.0033,1.9469,3.4396, 3.1548, 1.4180,0.2413,1.5298, 5.8009, 19.4520, 60.5753},
     {1.0689,2.1038,3.6039, 3.4927, 1.8283,0.2540,1.6715, 6.3509, 23.1531, 78.7099},
     {1.0891,2.1867,3.6160, 3.8031, 1.8994,0.2552,1.7174, 6.5131, 23.9170, 74.7039},
     {1.1007,2.2306,3.5689, 4.1549, 2.0382,0.2546,1.7351, 6.4948, 23.6464, 70.3780},
     {1.1568,2.4353,3.6459, 4.4064, 1.7179,0.2648,1.8786, 7.1749, 25.1766, 69.2821},
     {1.0909,2.1976,3.3831, 4.6700, 2.1277,0.2466,1.6707, 6.0197, 20.7657, 57.2663},
     {1.0756,2.1630,3.3178, 4.8852, 2.0489,0.2402,1.6169, 5.7644, 19.4568, 52.5009},
     {1.4282,3.5081,5.6767, 4.1964, 3.8946,0.3183,2.6889,13.4816, 54.3866,200.8321},
     {1.3127,3.1243,5.2988, 5.3891, 5.4133,0.2887,2.2897,10.8276, 43.5389,145.6109},
     {1.3128,3.1021,5.3385, 5.9611, 4.7562,0.2861,2.2509,10.5287, 41.7796,128.2973},
     {1.2553,2.9178,5.0862, 6.1206, 4.7122,0.2701,2.0636, 9.3051, 34.5977,107.9200},
     {1.3218,3.1444,5.4371, 5.6444, 4.0107,0.2827,2.2250,10.2454, 41.1162,124.4449},
     {1.3382,3.2043,5.4558, 5.4839, 3.6342,0.2838,2.2452,10.2519, 41.7251,124.9023},
     {1.5193,4.0053,6.5327,-0.1402, 6.7489,0.3213,2.8206,14.8878, 68.9103, 81.7257},
     {1.3517,3.2937,5.3213, 4.6466, 3.5714,0.2813,2.2418, 9.9952, 42.7939,132.1739},
     {1.2135,2.7962,4.7545, 4.5731, 4.4786,0.2483,1.8437, 7.5421, 29.3841,112.4579},
     {1.2937,3.1100,5.0393, 4.7546, 3.5031,0.2638,2.0341, 8.7101, 35.2992,109.4972},
     {1.2915,3.1023,4.9309, 4.6009, 3.4661,0.2611,2.0023, 8.4377, 34.1559,105.8911},
     {1.2089,2.7391,4.3482, 4.0047, 4.6497,0.2421,1.7487, 6.7262, 23.2153, 80.3108}};
  int i;
  if((atomic_num <= 0)||(atomic_num > NUM_ELEM)){
    return 1;
  }
  for(i = 0; i < NUM_SCATTERING_PARAMETERS; i++){
    a[i] = scfact[atomic_num-1][i];
    b[i] = scfact[atomic_num-1][i+NUM_SCATTERING_PARAMETERS];
  }
  return 0;
}

/****************************************************************************/

int pdb_find_box(double box[6], const char *fn, array *transf, double padding){
  char line[PDB_LINE_LENGTH+1];
  FILE *fp;
  pdbrecord p;
  double coord[3];
  int first_line = 1;
  long i;
  FOPEN(fp, fn, "r");
  if(fp == NULL){
    WARNING("Could not open file %s for reading.\n", fn);
    return 1;
  }
  while(!feof(fp)){
    if(NULL == fgets(line, PDB_LINE_LENGTH+1, fp)){
      break;
    }
    if(parse_pdb_line(line, &p)){
      if(first_line){
	transf_coord(coord, &p, transf, 0);
	box[0] = box[1] = coord[0];
	box[2] = box[3] = coord[1];
	box[4] = box[5] = coord[2];
	first_line = 0;
      }
      for(i = 0; i < transf->size[2]; i++){
	transf_coord(coord, &p, transf, i);
	box[0] = min_d(box[0], coord[0]); box[1] = max_d(box[1], coord[0]);
	box[2] = min_d(box[2], coord[1]); box[3] = max_d(box[3], coord[1]);
	box[4] = min_d(box[4], coord[2]); box[5] = max_d(box[5], coord[2]);
      }
    }
  }
  fclose(fp);
  if(first_line){
    WARNING("No atoms found in pdb file %s.\n", fn);
    return 2;
  }
  else {
    box[0] -= padding; box[2] -= padding; box[4] -= padding;
    box[1] += padding; box[3] += padding; box[5] += padding;
    return 0;
  }
}

/****************************************************************************/

int parse_pdb_line(char *line, pdbrecord *p){
#define MAX_FIELD_WIDTH 9
  char field[MAX_FIELD_WIDTH];
  int ret;

  if(strlen(line) < 78){
    return 0; /* Bad line */
  }
  if(0 == strncmp(line, "ATOM", 4)){
    ret = 1;
  }
  else if(0 == strncmp(line, "HETATM", 6)){
    ret = 2;
  }
  else {
    return 0; /* Other entries are unimportant */
  }

  /* Atom name, columns 13-16 */
  strncpy_term(p->atom_name, 5, line+12, 5);

  /* Residue name, columns 18-20 */
  strncpy_term(p->residue_name, 4, line+17, 4);

  /* x coordinate, columns 31-38 */
  strncpy_term(field, MAX_FIELD_WIDTH, line+30, 9);
  p->coord[0] = strtod(field, NULL);

  /* y coordinate, columns 39-46 */
  strncpy_term(field, MAX_FIELD_WIDTH, line+38, 9);
  p->coord[1] = strtod(field, NULL);

  /* z coordinate, columns 47-54 */
  strncpy_term(field, MAX_FIELD_WIDTH, line+46, 9);
  p->coord[2] = strtod(field, NULL);

  /* Occupancy, columns 55-60 */
  strncpy_term(field, MAX_FIELD_WIDTH, line+54, 7);
  p->occ = strtod(field, NULL);

  /* Temperature factor, columns 61-66 */
  strncpy_term(field, MAX_FIELD_WIDTH, line+60, 7);
  p->temp_fact = strtod(field, NULL);

  /* Element symbol, columns 77-78 */
  strncpy_term(p->element_sym, 3, line+76, 3);

  /* Look up atomic number */
  p->atomic_num = atomic_number_lookup(p->element_sym);

  return ret;

}

/****************************************************************************/

int read_pdb_transf_file(array *a, particle *p){
  matrix b, c;
  long i, j, k;
  if(param_isset(p->param, PAR_PDB_TRANSF_FILE_IN)){
    if(read_matrix_text(&b, get_param_string(p->param, PAR_PDB_TRANSF_FILE_IN))) return 1;
    if(6 == b.n){
      init_array(a, 3, 4, b.m);
      c.m = c.n = 3;
      for(i = 0; i < b.m; i++){
	c.data = get_array_entry_ptr(a, 0, 0, i);
	fill_matrix_diag(&c, 1);
	rotate3drows(&c, 0, 0, 1, -ANGLE_UNIT*get_matrix_entry(&b, i, 5));
	rotate3drows(&c, 1, 0, 0, -ANGLE_UNIT*get_matrix_entry(&b, i, 4));
	rotate3drows(&c, 0, 0, 1, -ANGLE_UNIT*get_matrix_entry(&b, i, 3));
	for(j = 0; j < 3; j++){
	  set_array_entry(a, j, 3, i, get_matrix_entry(&b, i, j));
	}
      }
    }
    else if(4 == b.n){
      init_array(a, 3, 4, b.m/3);
      for(k = 0; k < a->size[2]; k++){
	for(i = 0; i < 3; i++){
	  for(j = 0; j < 4; j++){
	    set_array_entry(a, i, j, k, get_matrix_entry(&b, 3*k+i, j));
	  }
	}
      }
    }
    else {
      free_matrix(&b);
      WARNING("Error reading pdb transform file %s: %i columns found, should be 4 or 6.\n", 
	      get_param_string(p->param, PAR_PDB_TRANSF_FILE_IN), (int)b.n);
      return 1;
    }
    free_matrix(&b);
  }
  else {
    init_array(a, 3, 4, 1);
    c.m = 3;
    c.n = 4;
    c.data = a->data;
    fill_matrix_diag(&c, 1);
  }
  return 0;
}

/****************************************************************************/

void transf_coord(double coord[3], pdbrecord *p, array *transf, long k){
  int i, j;
  for(i = 0; i < 3; i++){
    coord[i] = get_array_entry(transf, i, 3, k);
    for(j = 0; j < 3; j++){
      coord[i] += p->coord[j] * get_array_entry(transf, i, j, k);
    }
  }
}

/****************************************************************************/

const static char element_symbols[NUM_ELEM][3] = 
  {" H","HE","LI","BE"," B"," C"," N"," O"," F","NE",
   "NA","MG","AL","SI"," P"," S","CL","AR"," K","CA",
   "SC","TI"," V","CR","MN","FE","CO","NI","CU","ZN",
   "GA","GE","AS","SE","BR","KR","RB","SR"," Y","ZR",
   "NB","MO","TC","RU","RH","PD","AG","CD","IN","SN",
   "SB","TE"," I","XE","CS","BA","LA","CE","PR","ND",
   "PM","SM","EU","GD","TB","DY","HO","ER","TM","YB",
   "LU","HF","TA"," W","RE","OS","IR","PT","AU","HG",
   "TL","PB","BI","PO","AT","RN","FR","RA","AC","TH",
   "PA"," U","NP","PU","AM","CM","BK","CF"};

/****************************************************************************/

int atomic_number_lookup(char *s){
  const char first = 'B';
  const int numchar = 1 + ('Y' - first);
  static int *one_letter_elem = NULL;
  int i;
  if(' ' == s[0]){ /* If first character of s is space, we have a one letter element. Use a quick lookup. */
    if(NULL == one_letter_elem){ /* Perform initialization the first time this happens. */
      one_letter_elem = malloc(numchar * sizeof(int));
      for(i = 0; i < numchar; i++){
	one_letter_elem[i] = 0;
      }
      one_letter_elem['B' - first] = 5;
      one_letter_elem['C' - first] = 6;
      one_letter_elem['F' - first] = 9;
      one_letter_elem['H' - first] = 1;
      one_letter_elem['I' - first] = 53;
      one_letter_elem['K' - first] = 19;
      one_letter_elem['N' - first] = 7;
      one_letter_elem['O' - first] = 8;
      one_letter_elem['P' - first] = 15;
      one_letter_elem['S' - first] = 16;
      one_letter_elem['U' - first] = 92;
      one_letter_elem['V' - first] = 23;
      one_letter_elem['W' - first] = 74;
      one_letter_elem['Y' - first] = 39;
    }
    i = s[1] - first;
    if((i >= 0)&&(i < numchar)){
      return one_letter_elem[i];
    }
    else {
      return 0;
    }
  }
  else {
    i = 0;
    while((i < NUM_ELEM) && strncmp(s, element_symbols[i], 2)){
      i++;
    }
    if(i < NUM_ELEM){
      return i+1;
    }
    else {
      return 0;
    }
  }
}

/****************************************************************************/

void atom_name_lookup(int number, char name[3]){
  if((number > 0)&&(number <= NUM_ELEM)){
    name[0] = element_symbols[number-1][0];
    name[1] = element_symbols[number-1][1];
    name[2] = '\0';
  }
}

/****************************************************************************/

int nhydro(char *seq, char *pos){

  /* Water molecule*/
  if (strcmp(seq, "HOH") == 0)
    {
      return 2;
    }

  /*Peptide bond C and O*/
  if (strcmp(pos, " C  ") == 0 || strcmp(pos, " O  ") == 0)
    {
      return 0;
    }

  /*Peptide bond N*/
  if (strcmp(pos, " N  ") == 0)
    {
      if (strcmp(seq, "PRO") == 0)
        {
	  return 0;
        }
      else
        {
	  return 1;
        }
    }

  /*Alpha carbon CA*/
  if (strcmp(pos, " CA ") == 0)
    {
      if (strcmp(seq, "GLY") == 0)
        {
	  return 2;
        }
      else
        {
	  return 1;
        }
    }

  /*Side chain position nomenclature depends on amino acid or base*/
  if (strcmp(seq, "ALA") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 3;
    }
  else if (strcmp(seq, "VAL") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 1;
      if (strcmp(pos, " CG1") == 0) return 3;
      if (strcmp(pos, " CG2") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, "LEU") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 1;
      if (strcmp(pos, " CD1") == 0) return 3;
      if (strcmp(pos, " CD2") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, "ILE") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 1;
      if (strcmp(pos, " CG1") == 0) return 3;
      if (strcmp(pos, " CG2") == 0) return 2;
      if (strcmp(pos, " CD1") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, "SER") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " OG ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "MET") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      if (strcmp(pos, " CE ") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, "THR") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 1;
      if (strcmp(pos, " OG1") == 0) return 1;
      if (strcmp(pos, " CG2") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, "PHE") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CD1") == 0) return 1;
      if (strcmp(pos, " CD2") == 0) return 1;
      if (strcmp(pos, " CE1") == 0) return 1;
      if (strcmp(pos, " CE2") == 0) return 1;
      if (strcmp(pos, " CZ ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "TYR") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CD1") == 0) return 1;
      if (strcmp(pos, " CD2") == 0) return 1;
      if (strcmp(pos, " CE1") == 0) return 1;
      if (strcmp(pos, " CE2") == 0) return 1;
      if (strcmp(pos, " OH ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "TRP") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CD1") == 0) return 1;
      if (strcmp(pos, " NE1") == 0) return 1;
      if (strcmp(pos, " CE3") == 0) return 1;
      if (strcmp(pos, " CZ2") == 0) return 1;
      if (strcmp(pos, " CZ3") == 0) return 1;
      if (strcmp(pos, " CH2") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "CYS") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " SG ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "PRO") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      if (strcmp(pos, " CD ") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "ASP") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "ASN") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " ND2") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "GLU") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "GLN") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      if (strcmp(pos, " NE2") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "HIS") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " ND1") == 0) return 1;
      if (strcmp(pos, " CD2") == 0) return 1;
      if (strcmp(pos, " CE1") == 0) return 1;
      if (strcmp(pos, " NE2") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "ARG") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      if (strcmp(pos, " CD ") == 0) return 2;
      if (strcmp(pos, " NE ") == 0) return 1;
      if (strcmp(pos, " NH1") == 0) return 2;
      if (strcmp(pos, " NH2") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, "LYS") == 0)
    {
      if (strcmp(pos, " CB ") == 0) return 2;
      if (strcmp(pos, " CG ") == 0) return 2;
      if (strcmp(pos, " CD ") == 0) return 2;
      if (strcmp(pos, " CE ") == 0) return 2;
      if (strcmp(pos, " NZ ") == 0) return 3;
      return 0;
    }
  else if (strcmp(seq, " U") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 1;
      if (strcmp(pos, " O2'") == 0) return 1;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " N3 ") == 0) return 1;
      if (strcmp(pos, " C5 ") == 0) return 1;
      if (strcmp(pos, " C6 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "  T") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 1;
      if (strcmp(pos, " O2'") == 0) return 1;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " N3 ") == 0) return 1;
      if (strcmp(pos, " C7 ") == 0) return 3;
      if (strcmp(pos, " C6 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "  C") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 1;
      if (strcmp(pos, " O2'") == 0) return 1;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " N4 ") == 0) return 2;
      if (strcmp(pos, " C5 ") == 0) return 1;
      if (strcmp(pos, " C6 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "  A") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 1;
      if (strcmp(pos, " O2'") == 0) return 1;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " C8 ") == 0) return 1;
      if (strcmp(pos, " N6 ") == 0) return 2;
      if (strcmp(pos, " C2 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, "  G") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 1;
      if (strcmp(pos, " O2'") == 0) return 1;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " C8 ") == 0) return 1;
      if (strcmp(pos, " N1 ") == 0) return 1;
      if (strcmp(pos, " N2 ") == 0) return 2;
      return 0;
    }
  else if (strcmp(seq, " DT") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 2;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " N3 ") == 0) return 1;
      if (strcmp(pos, " C7 ") == 0) return 3;
      if (strcmp(pos, " C6 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, " DC") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 2;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " N4 ") == 0) return 2;
      if (strcmp(pos, " C5 ") == 0) return 1;
      if (strcmp(pos, " C6 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, " DA") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 2;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " C8 ") == 0) return 1;
      if (strcmp(pos, " N6 ") == 0) return 2;
      if (strcmp(pos, " C2 ") == 0) return 1;
      return 0;
    }
  else if (strcmp(seq, " DG") == 0)
    {
      if (strcmp(pos, " C5'") == 0) return 2;
      if (strcmp(pos, " C4'") == 0) return 1;
      if (strcmp(pos, " C3'") == 0) return 1;
      if (strcmp(pos, " C2'") == 0) return 2;
      if (strcmp(pos, " C1'") == 0) return 1;
      if (strcmp(pos, " C8 ") == 0) return 1;
      if (strcmp(pos, " N1 ") == 0) return 1;
      if (strcmp(pos, " N2 ") == 0) return 2;
      return 0;
    }

  return 0;
}
