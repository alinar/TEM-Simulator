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


#ifndef SIMULATION_HEADER
#define SIMULATION_HEADER
#include "input.h"
#include "macros.h"
#include "structs.h"


/***********************************************************************
 * Function:  new_simulation
 * Purpose:   Allocate and return a pointer to a new simulation struct.
 * Return:    Pointer to a new simulation struct.
 */

simulation *new_simulation();


/***********************************************************************
 * Function:  delete_simulation
 * Purpose:   Free memory allocated by simulation struct.
 * Arguments: s - Pointer to simulation to be deleted.
 */

void delete_simulation(simulation *s);


/***********************************************************************
 * Function:  simulation_init
 * Purpose:   Do initializations necessary to make simulation struct ready
 *            for use in computations once parameters have been read from 
 *            input file. Calling simulation_init on a geometry object 
 *            which is already initialized has no effect. Once the 
 *            simulation object has been initialized, input parameters can 
 *            not be changed.
 * Arguments: s - Pointer to simulation struct to be initialized.
 * Return:    0 on success, nonzero on failure.
 */

int simulation_init(simulation *s);


/***********************************************************************
 * Function:  simulation_write_log
 * Purpose:   Write input parameters of simulation object and all simulation
 *            components to the log file.
 * Arguments: s - Pointer to simulation struct.
 */

void simulation_write_log(simulation *s);


/***********************************************************************
 * Function:  simulation_create_log
 * Purpose:   Create a log file to be used in the simulation. If the name
 *            of a log file is specified in the input, that file is used,
 *            overwriting the file if it already exists. If no file is 
 *            specified, a unique name is generated, and it is checked
 *            that the file does not already exist.
 * Arguments: s - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure to create log file.
 */

int simulation_create_log(simulation *s);


/***********************************************************************
 * Function:  read_input_data
 * Purpose:   Read input data from a file, creating simulation components
 *            and assigning parameters values.
 * Arguments: s - Pointer to simulation struct in which the simulation
 *            components are stored.
 *            fn - Name of input file.
 * Return:    0 on success, nonzero on failure. Fails if the input file
 *            can not be opened. If the file contains incorrect syntax
 *            a warning is printed for each line that can not be parsed.
 */

int read_input_data(simulation *s, char *fn);


/***********************************************************************
 * Function:  add_simcomp
 * Purpose:   Add a simulation component to the simulation.
 * Arguments: s - Pointer to simulation struct to which a component is 
 *            added.
 *            class - String identifying the type of component.
 *            name - Name of the component.
 *            reuse - If reuse is nonzero, and a simulation component
 *            with the same name and type already exists, it is reused.
 *            If reuse is zero, and a simulation component with the 
 *            same name and type already exists, no component is added
 *            and a null pointer is returned. The names of two components
 *            are only considered to be the same if the name is a 
 *            nonempty string.
 * Return:    Pointer to the parameter table of the added component
 */

param_table *add_simcomp(simulation *s, const char *class, const char *name, int reuse);


/***********************************************************************
 * Function:  remove_simcomp
 * Purpose:   Remove a simulation component from the simulation.
 * Arguments: s - Pointer to simulation struct from which a component is 
 *            removed.
 *            class - String identifying the type of component.
 *            name - Name of the component. If the name is an empty
 *            the last component of that type which was added is removed
 *            regardless of its name. If no component with the given
 *            name and type exists, nothing is changed.
 * Return:    Nonzero if a component was removed.
 */

int remove_simcomp(simulation *s, const char *class, const char *name);


/***********************************************************************
 * Function:  get_detector
 * Purpose:   Get a pointer to the detector component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of detector component. If name is an empty 
 *            string and there is just one detector component, a pointer
 *            to it is returned. Otherwise, if there is a detector 
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

detector *get_detector(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_electronbeam
 * Purpose:   Get a pointer to the electronbeam component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of electronbeam component. If name is an empty 
 *            string and there is just one electronbeam component, a pointer
 *            to it is returned. Otherwise, if there is an electronbeam
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

electronbeam *get_electronbeam(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_geometry
 * Purpose:   Get a pointer to the geometry component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of geometry component. If name is an empty 
 *            string and there is just one geometry component, a pointer
 *            to it is returned. Otherwise, if there is a geometry
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

geometry *get_geometry(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_optics
 * Purpose:   Get a pointer to the optics component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of optics component. If name is an empty 
 *            string and there is just one optics component, a pointer
 *            to it is returned. Otherwise, if there is an optics
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

optics *get_optics(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_particle
 * Purpose:   Get a pointer to the particle component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of particle component. If name is an empty 
 *            string and there is just one particle component, a pointer
 *            to it is returned. Otherwise, if there is a particle
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

particle *get_particle(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_particleset
 * Purpose:   Get a pointer to the particleset component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of particle component. If name is an empty 
 *            string and there is just one particleset component, a pointer
 *            to it is returned. Otherwise, if there is a particleset
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

particleset *get_particleset(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_sample
 * Purpose:   Get a pointer to the sample component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of sample component. If name is an empty 
 *            string and there is just one sample component, a pointer
 *            to it is returned. Otherwise, if there is a sample
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

sample *get_sample(simulation *s, const char *name);


/***********************************************************************
 * Function:  get_volume
 * Purpose:   Get a pointer to the volume component with a given name.
 * Arguments: s - Pointer to simulation struct.
 *            name - name of volume component. If name is an empty 
 *            string and there is just one volume component, a pointer
 *            to it is returned. Otherwise, if there is a volume
 *            component with the given name, a pointer to it is returned.
 *            In all other cases a null pointer is returned.
 */

volume *get_volume(simulation *s, const char *name);


/***********************************************************************
 * Function:  generate_micrographs
 * Purpose:   Run simulation of micrographs with the simulation components
 *            currently in the simulation struct.
 * Arguments: s - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure.
 */

int generate_micrographs(simulation *s);


/***********************************************************************
 * Function:  generate_volumes
 * Purpose:   Generate volume maps for all volume components in the 
 *            simulation struct.
 * Arguments: s - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure.
 */

int generate_volumes(simulation *s);


/***********************************************************************
 * Function:  generate_particle_maps
 * Purpose:   Generate particle maps for all particle components in the 
 *            simulation struct.
 * Arguments: s - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure.
 */

int generate_particle_maps(simulation *s);


/***********************************************************************
 * Function:  simulation_run
 * Purpose:   Run all simulations specified by input parameters 
 *            PAR_GENERATE_MICROGRAPHS, PAR_GENERATE_PARTICLE_MAPS, and
 *            PAR_GENERATE_VOLUMES.
 * Arguments: s - Pointer to simulation struct.
 * Return:    0 on success, nonzero on failure.
 */

int simulation_run(simulation *s);

#endif
