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
#include <time.h>
#include "simulation.h"
#include "detector.h"
#include "electronbeam.h"
#include "geometry.h"
#include "input.h"
#include "log.h"
#include "misc.h"
#include "optics.h"
#include "particle.h"
#include "particleset.h"
#include "sample.h"
#include "volume.h"
#include "wavefunction.h"

/****************************************************************************/

param_table *simulation_param_table();

int simulation_check_input(simulation *s);

enum if_anonymous {no_match, match_unique, match_last};

int find_detector(simulation *s, const char *name, enum if_anonymous an);

int find_electronbeam(simulation *s, const char *name, enum if_anonymous an);

int find_geometry(simulation *s, const char *name, enum if_anonymous an);

int find_optics(simulation *s, const char *name, enum if_anonymous an);

int find_particle(simulation *s, const char *name, enum if_anonymous an);

int find_particleset(simulation *s, const char *name, enum if_anonymous an);

int find_sample(simulation *s, const char *name, enum if_anonymous an);

int find_volume(simulation *s, const char *name, enum if_anonymous an);

/****************************************************************************/

simulation *new_simulation(){
  simulation *s = malloc(sizeof(simulation));
  s->param = simulation_param_table();
  s->init = 0;
  s->num_detector = 0;
  s->num_electronbeam = 0;
  s->num_geometry = 0;
  s->num_optics = 0;
  s->num_particle = 0;
  s->num_particleset = 0;
  s->num_sample = 0;
  s->num_volume = 0;
  return s;
}

/****************************************************************************/

void delete_simulation(simulation *s){
  int i;
  for(i = 0; i < s->num_detector; i++){
    delete_detector(s->detector[i]);
  }
  for(i = 0; i < s->num_electronbeam; i++){
    delete_electronbeam(s->electronbeam[i]);
  }
  for(i = 0; i < s->num_geometry; i++){
    delete_geometry(s->geometry[i]);
  }
  for(i = 0; i < s->num_optics; i++){
    delete_optics(s->optics[i]);
  }
  for(i = 0; i < s->num_particleset; i++){
    delete_particleset(s->particleset[i]);
  }
  for(i = 0; i < s->num_sample; i++){
    delete_sample(s->sample[i]);
  }
  for(i = 0; i < s->num_volume; i++){
    delete_volume(s->volume[i]);
  }
  delete_param_table(s->param);
  free(s);
}

/****************************************************************************/

param_table *simulation_param_table(){
  param_table *pt = new_param_table(5, TYPE_SIMULATION, "");
  add_param_def(pt, PAR_GENERATE_MICROGRAPHS, "b", NO_STRING);
  add_param_def(pt, PAR_GENERATE_VOLUMES, "b", NO_STRING);
  add_param_def(pt, PAR_GENERATE_PARTICLE_MAPS, "b", NO_STRING);
  add_param_opt(pt, PAR_LOG_FILE, "s");
  add_param_opt(pt, PAR_RAND_SEED, "i");
  set_comp_descr(pt, "The simulation component defines what should happen when \
the simulator is run. Three kinds of actions are possible: generate a series of \
simulated micrographs (the main purpose of the program), generate one or several \
volume maps of the sample used in the simulation, and generate maps of one or \
several individual particles.");
  set_param_descr(pt, PAR_GENERATE_MICROGRAPHS, "Tells the simulator to generate \
a micrograph or a tilt series. Requires one each of the components sample, \
electronbeam, geometry, and optics, at least one detector, and for a meaningful \
simulation at least one particle and one particleset component.");
  set_param_descr(pt, PAR_GENERATE_VOLUMES, "Tells the simulator to generate one \
or several volume maps of all or part of the sample, as specified by volume \
components. Requires one sample component and (to be meaningful) at least one \
particle and one particleset component. Might also require an electronbeam \
component, since the imaginary or absorption potential can depend on the \
accelertion energy.");
  set_param_descr(pt, PAR_GENERATE_PARTICLE_MAPS, "Tells the simulator to \
generate map files for all particle components. This also happens as a side \
effect whenever one of the other two simulation tasks are run, at least for \
particles used in a particleset.");
  set_param_descr(pt, PAR_LOG_FILE, "The name of a log file which is continuously \
updated during the simulation. If the parameter is not defined, a unique name is \
automatically generated, based on the time when the simulation starts.");
  set_param_descr(pt, PAR_RAND_SEED, "A number which is used to seed the random \
number generator at the start of a simulation. If the parameter is not defined, \
the random number generator is instead seeded with the current time.");
  return pt;
}

/****************************************************************************/

int simulation_init(simulation *s){
  unsigned int seed;
  if(s->init) return 0;
  if(simulation_check_input(s)){
    WARNING("Invalid input of simulation parameters.\n");
    return 1;
  }
  if(param_isset(s->param, PAR_RAND_SEED)){
    srand((unsigned int)get_param_int(s->param, PAR_RAND_SEED));
  }
  else {
    seed = (unsigned int)time(NULL);
    srand(seed);
    write_log_comment("Random seed taken from clock: %i\n", seed);
  }
  param_table_set_lock(s->param, 1);
  s->init = 1;
  write_log_comment("Simulation parameters object initialized.\n\n");
  return 0;
}

/****************************************************************************/

int simulation_check_input(simulation *s){
  return check_params(s->param);
}

/****************************************************************************/

void simulation_write_log(simulation *s){
  int i;
  write_parameters_to_log(s->param);
  for(i = 0; i < s->num_detector; i++){
    detector_write_log(s->detector[i]);
  }
  for(i = 0; i < s->num_electronbeam; i++){
    electronbeam_write_log(s->electronbeam[i]);
  }
  for(i = 0; i < s->num_geometry; i++){
    geometry_write_log(s->geometry[i]);
  }
  for(i = 0; i < s->num_optics; i++){
    optics_write_log(s->optics[i]);
  }
  for(i = 0; i < s->num_particle; i++){
    particle_write_log(s->particle[i]);
  }
  for(i = 0; i < s->num_particleset; i++){
    particleset_write_log(s->particleset[i]);
  }
  for(i = 0; i < s->num_sample; i++){
    sample_write_log(s->sample[i]);
  }
  for(i = 0; i < s->num_volume; i++){
    volume_write_log(s->volume[i]);
  }
  write_block_sep_log("end of input");
}

/****************************************************************************/

#define LOGFILE_MAX_NAME_LENGTH 40
int simulation_create_log(simulation *s){
  unsigned int id, i;
  int ret;
  if(param_isset(s->param, PAR_LOG_FILE)){
    ret = create_log(get_param_string(s->param, PAR_LOG_FILE), 1);
  }
  else {
    char logfile[LOGFILE_MAX_NAME_LENGTH];
    id = (unsigned int)time(NULL);
    ret = 1;
    for(i = 0; i < 10; i++){
      SNPRINTF(logfile, LOGFILE_MAX_NAME_LENGTH, "TEM-simulator_%u.log", id+i);
      if(0 == create_log(logfile, 0)){
	ret = 0;
	break;
      }
    }
  }
  if(ret){
    WARNING("Failed to create log file.\n");
  }
  else {
    write_log_comment("Log file created by TEM-simulator, version %s.\n", VERSION_NUMBER);
  }
  return ret;
}

/****************************************************************************/

int read_input_data(simulation *s, char *fn){
  const char *class, *name;
  param_table *pt;
  if(s == NULL){
    return 1;
  }
  if(NULL == open_input_file(fn)){
    WARNING("Could not open input file %s.\n", fn);
    return 1;
  }
  while(1){
    class = current_input_block_class();
    name = current_input_block_name();
    if(0 == strcmp(class, TYPE_SIMULATION)){
      read_input_block(s->param);
    }
    else if(0 == strcmp(class, "end")){
      break;
    }
    else {
      pt = add_simcomp(s, class, name, 0);
      if(pt != NULL){
        read_input_block(pt);
      }
      else {
	WARNING("Skipping input block %s %s.\n", class, name);
	skip_input_block();
      }
    }
  }
  close_input_file();
  return 0;
}

/****************************************************************************/

param_table *add_simcomp(simulation *s, const char *class, const char *name, int reuse){
  int i, n;
  if(0 == strcmp(class, TYPE_DETECTOR)){
    i = find_detector(s, name, no_match);
    if(i < 0){
      n = s->num_detector;
      if(n < MAX_NUM_DETECTOR){
        s->detector[n] = new_detector(name);
	s->num_detector++;
        return s->detector[n]->param;
      }
    }
    else if(reuse){      
      return s->detector[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_ELECTRONBEAM)){
    i = find_electronbeam(s, name, no_match);
    if(i < 0){
      n = s->num_electronbeam;
      if(n < MAX_NUM_ELECTRONBEAM){
        s->electronbeam[n] = new_electronbeam(name);
	s->num_electronbeam++;
        return s->electronbeam[n]->param;
      }
    }
    else if(reuse){      
      return s->electronbeam[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_GEOMETRY)){
    i = find_geometry(s, name, no_match);
    if(i < 0){
      n = s->num_geometry;
      if(n < MAX_NUM_GEOMETRY){
        s->geometry[n] = new_geometry(name);
	s->num_geometry++;
        return s->geometry[n]->param;
      }
    }
    else if(reuse){      
      return s->geometry[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_OPTICS)){
    i = find_optics(s, name, no_match);
    if(i < 0){
      n = s->num_optics;
      if(n < MAX_NUM_OPTICS){
        s->optics[n] = new_optics(name);
	s->num_optics++;
        return s->optics[n]->param;
      }
    }
    else if(reuse){      
      return s->optics[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_PARTICLE)){
    i = find_particle(s, name, no_match);
    if(i < 0){
      n = s->num_particle;
      if(n < MAX_NUM_PARTICLE){
        s->particle[n] = new_particle(name);
	s->num_particle++;
        return s->particle[n]->param;
      }
    }
    else if(reuse){      
      return s->particle[i]->param;
    }
  }  
  else if(0 == strcmp(class, TYPE_PARTICLESET)){
    i = find_particleset(s, name, no_match);
    if(i < 0){
      n = s->num_particleset;
      if(n < MAX_NUM_PARTICLESET){
        s->particleset[n] = new_particleset(name);
	s->num_particleset++;
        return s->particleset[n]->param;
      }
    }
    else if(reuse){      
      return s->particleset[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_SAMPLE)){
    i = find_sample(s, name, no_match);
    if(i < 0){
      n = s->num_sample;
      if(n < MAX_NUM_SAMPLE){
        s->sample[n] = new_sample(name);
	s->num_sample++;
        return s->sample[n]->param;
      }
    }
    else if(reuse){      
      return s->sample[i]->param;
    }
  }
  else if(0 == strcmp(class, TYPE_VOLUME)){
    i = find_volume(s, name, no_match);
    if(i < 0){
      n = s->num_volume;
      if(n < MAX_NUM_VOLUME){
        s->volume[n] = new_volume(name);
	s->num_volume++;
        return s->volume[n]->param;
      }
    }
    else if(reuse){      
      return s->volume[i]->param;
    }
  }
  else {
    WARNING("Unknown component %s\n", class);
    return NULL;
  }
  if(i < 0){
    WARNING("Maximum number of %s exceeded.\n", class);
  }
  else {
    WARNING("A %s with the name %s already exists.\n", class, name);
  }
  return NULL;
}

/****************************************************************************/

int remove_simcomp(simulation *s, const char *class, const char *name){
  void **list = NULL;
  int num, i, j;
  if(0 == strcmp(class, TYPE_DETECTOR)){
    i = find_detector(s, name, match_last);
    if(i < 0) return 0;
    delete_detector(s->detector[i]);
    list = (void**)s->detector;
    num = --(s->num_detector);
  }
  else if(0 == strcmp(class, TYPE_ELECTRONBEAM)){
    i = find_electronbeam(s, name, match_last);
    if(i < 0) return 0;
    delete_electronbeam(s->electronbeam[i]);
    list = (void**)s->electronbeam;
    num = --(s->num_electronbeam);
  }
  else if(0 == strcmp(class, TYPE_GEOMETRY)){
    i = find_geometry(s, name, match_last);
    if(i < 0) return 0;
    delete_geometry(s->geometry[i]);
    list = (void**)s->geometry;
    num = --(s->num_geometry);
  }
  else if(0 == strcmp(class, TYPE_OPTICS)){
    i = find_optics(s, name, match_last);
    if(i < 0) return 0;
    delete_optics(s->optics[i]);
    list = (void**)s->optics;
    num = --(s->num_optics);
  }
  else if(0 == strcmp(class, TYPE_PARTICLE)){
    i = find_particle(s, name, match_last);
    if(i < 0) return 0;
    delete_particle(s->particle[i]);
    list = (void**)s->particle;
    num = --(s->num_particle);
  }
  else if(0 == strcmp(class, TYPE_PARTICLESET)){
    i = find_particleset(s, name, match_last);
    if(i < 0) return 0;
    delete_particleset(s->particleset[i]);
    list = (void**)s->particleset;
    num = --(s->num_particleset);
  }
  else if(0 == strcmp(class, TYPE_SAMPLE)){
    i = find_sample(s, name, match_last);
    if(i < 0) return 0;
    delete_sample(s->sample[i]);
    list = (void**)s->sample;
    num = --(s->num_sample);
  }
  else if(0 == strcmp(class, TYPE_VOLUME)){
    i = find_volume(s, name, match_last);
    if(i < 0) return 0;
    delete_volume(s->volume[i]);
    list = (void**)s->volume;
    num = --(s->num_volume);
  }
  else {
    return 0;
  }
  for(j = i; j < num; j++){
    list[j] = list[j+1];
  }
  return 1;
}

/****************************************************************************/

detector *get_detector(simulation *s, const char *name){
  int i = find_detector(s, name, match_unique);
  if(i >= 0){
    return s->detector[i];
  }
  return NULL;
}

/****************************************************************************/

electronbeam *get_electronbeam(simulation *s, const char *name){
  int i = find_electronbeam(s, name, match_unique);
  if(i >= 0){
    return s->electronbeam[i];
  }
  return NULL;
}

/****************************************************************************/

geometry *get_geometry(simulation *s, const char *name){
  int i = find_geometry(s, name, match_unique);
  if(i >= 0){
    return s->geometry[i];
  }
  return NULL;
}

/****************************************************************************/

optics *get_optics(simulation *s, const char *name){
  int i = find_optics(s, name, match_unique);
  if(i >= 0){
    return s->optics[i];
  }
  return NULL;
}

/****************************************************************************/

particle *get_particle(simulation *s, const char *name){
  int i = find_particle(s, name, match_unique);
  if(i >= 0){
    return s->particle[i];
  }
  return NULL;
}

/****************************************************************************/

particleset *get_particleset(simulation *s, const char *name){
  int i = find_particleset(s, name, match_unique);
  if(i >= 0){
    return s->particleset[i];
  }
  return NULL;
}

/****************************************************************************/

sample *get_sample(simulation *s, const char *name){
  int i = find_sample(s, name, match_unique);
  if(i >= 0){
    return s->sample[i];
  }
  return NULL;
}

/****************************************************************************/

volume *get_volume(simulation *s, const char *name){
  int i = find_volume(s, name, match_unique);
  if(i >= 0){
    return s->volume[i];
  }
  return NULL;
}

/****************************************************************************/

int find_detector(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_detector == 1){
      return 0;
    }
    if(an == match_last && s->num_detector > 0){
      return s->num_detector - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_detector; i++){
    if(0 == strcmp(name, s->detector[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_electronbeam(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_electronbeam == 1){
      return 0;
    }
    if(an == match_last && s->num_electronbeam > 0){
      return s->num_electronbeam - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_electronbeam; i++){
    if(0 == strcmp(name, s->electronbeam[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_geometry(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_geometry == 1){
      return 0;
    }
    if(an == match_last && s->num_geometry > 0){
      return s->num_geometry - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_geometry; i++){
    if(0 == strcmp(name, s->geometry[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_optics(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_optics == 1){
      return 0;
    }
    if(an == match_last && s->num_optics > 0){
      return s->num_optics - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_optics; i++){
    if(0 == strcmp(name, s->optics[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_particle(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_particle == 1){
      return 0;
    }
    if(an == match_last && s->num_particle > 0){
      return s->num_particle - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_particle; i++){
    if(0 == strcmp(name, s->particle[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_particleset(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_particleset == 1){
      return 0;
    }
    if(an == match_last && s->num_particleset > 0){
      return s->num_particleset - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_particleset; i++){
    if(0 == strcmp(name, s->particleset[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_sample(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_sample == 1){
      return 0;
    }
    if(an == match_last && s->num_sample > 0){
      return s->num_sample - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_sample; i++){
    if(0 == strcmp(name, s->sample[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int find_volume(simulation *s, const char *name, enum if_anonymous an){
  int n = -1, i;
  if(0 == strcmp(name, "")){
    if(an == match_unique && s->num_volume == 1){
      return 0;
    }
    if(an == match_last && s->num_volume > 0){
      return s->num_volume - 1;
    }
    return -1;
  }
  for(i = 0; i < s->num_volume; i++){
    if(0 == strcmp(name, s->volume[i]->param->name)){
      n = i;
      break;
    }
  }
  return n;
}

/****************************************************************************/

int generate_micrographs(simulation *s){
  double slice_th;
  long tilt, ntilts;
  wavefunction *wf;
  double pix_size, det_area[4], det_area_all[4];
  int i;

  if((0 == s->num_detector) || (0 == s->num_electronbeam) || (0 == s->num_geometry)
     || (0 == s->num_optics) || (0 == s->num_sample)){
    WARNING("Incomplete input data for generating micrographs.\n");
    return 1;
  }
  if(simulation_init(s)) return 1;
  if(geometry_init(s->geometry[0], s)) return 1;
  if(electronbeam_init(s->electronbeam[0], s)) return 1;
  if(optics_init(s->optics[0], s)) return 1;
  if(sample_init(s->sample[0], s)) return 1;
  for(i = 0; i < s->num_particleset; i++){
    if(particleset_init(s->particleset[i], s)) return 1;
  }
  if(detector_init_all(s)) return 1;
  ntilts = get_param_long(s->geometry[0]->param, PAR_NTILTS);
  /* Find detector area which needs to be covered by wavefunction */
  pix_size = detector_get_pixel_size(s->detector[0]);
  vecf2d_xyrange(&(s->detector[0]->count), det_area_all);
  for(i = 1; i < s->num_detector; i++){
    pix_size = min_d(pix_size, detector_get_pixel_size(s->detector[i]));
    vecf2d_xyrange(&(s->detector[i]->count), det_area);
    det_area_all[0] = min_d(det_area_all[0], det_area[0]);
    det_area_all[1] = max_d(det_area_all[1], det_area[1]);
    det_area_all[2] = min_d(det_area_all[2], det_area[2]);
    det_area_all[3] = max_d(det_area_all[3], det_area[3]);
  }
  wf = new_wavefunction(s->electronbeam[0], s->optics[0]);
  if(wavefunction_init(wf, pix_size, det_area_all)){
    delete_wavefunction(wf);
    return 1;
  }
  slice_th = wf->pixel_size * wf->pixel_size / (4*M_PI)
    * wave_number(electronbeam_get_acc_energy(s->electronbeam[0]));
  write_log_comment("\nGenerating micrographs.\n");
  for(tilt = 0; tilt < ntilts; tilt++){
    if(wavefunction_propagate(s, wf, slice_th, tilt)){
      WARNING("Error computing projection.\n");
      break;
    }
    for(i = 0; i < s->num_detector; i++){
      if(detector_get_intensity(s->detector[i], wf, tilt)){
	WARNING("Error computing intensity.\n");
	break;
      }
      if(get_param_boolean(s->detector[i]->param, PAR_USE_QUANTIZATION)){
	if(detector_apply_quantization(s->detector[i])){
	  WARNING("Error computing quantization noise.\n");
	  break;
	}
      }
      if(detector_apply_mtf(s->detector[i])){
	WARNING("Error applying MTF.\n");
	break;
      }
      detector_write_image(s->detector[i]);
    }
  }
  delete_wavefunction(wf);
  write_log_comment("\nSimulation complete.\n");
  return 0;
}

/****************************************************************************/

int generate_volumes(simulation *s){
  electronbeam *ed = NULL;
  int i;
  if(0 == s->num_sample){
    WARNING("Incomplete input data for generating volumes.\n");
    return 1;
  }
  if(simulation_init(s)) return 1;
  if(0 < s->num_electronbeam){
    ed = s->electronbeam[0];
  }
  if(sample_init(s->sample[0], s)) return 1;
  for(i = 0; i < s->num_particleset; i++){
    if(particleset_init(s->particleset[i], s)) return 1;
  }

  write_log_comment("\nGenerating volumes.\n");

  for(i = 0; i < s->num_volume; i++){
    volume_init(s->volume[i], s);
    volume_get_potential(s->volume[i], s);
    volume_write_maps(s->volume[i]);
    volume_reset(s->volume[i]);
  }
  return 0;
}

/****************************************************************************/

int generate_particle_maps(simulation *s){
  int i;
  particle *p;
  const char *source;
  if(simulation_init(s)) return 1;
  for(i = 0; i < s->num_particle; i++){
    p = s->particle[i];
    source = get_param_string(p->param, PAR_SOURCE);
    if((0 == strcmp(source, "pdb") || 0 == strcmp(source, "random"))
       && (param_isset(p->param, PAR_MAP_FILE_RE_OUT) || param_isset(p->param, PAR_MAP_FILE_IM_OUT))){
      particle_init(s->particle[i], s);
    }
  }
  return 0;
}

/****************************************************************************/

int simulation_run(simulation *s){
  int ret = 0;
  if(get_param_boolean(s->param, PAR_GENERATE_PARTICLE_MAPS)){
    if(generate_particle_maps(s)) ret = 1;
  }

  if(get_param_boolean(s->param, PAR_GENERATE_VOLUMES)){
    if(generate_volumes(s)) ret = 1;
  }

  if(get_param_boolean(s->param, PAR_GENERATE_MICROGRAPHS)){
    if(generate_micrographs(s)) ret = 1;
  }
  return ret;
}
