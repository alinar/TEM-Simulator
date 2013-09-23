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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "log.h"
#include "macros.h"
#include "misc.h"
#include "simulation.h"

/****************************************************************************/

void print_info(void){
  printf("\nTEM-simulator version %s\n", VERSION_NUMBER);
  printf("compiled on %s at %s\n\n", __DATE__, __TIME__);
  print_with_line_breaks("TEM-simulator is an open source program for simulation \
of transmission electron microscope images and tilt series.\n", 80, 0);
  printf("Usage: TEM-simulator <input file>\n\n");
  printf("For more information, type\n");
  printf("TEM-simulator -help\n");
  printf("or see <http://TEM-simulator.sourceforge.net>.\n\n");
  print_with_line_breaks("Copyright 2008-2010, Hans Rullgard, Stockholm University \
and Lars-Goran Ofverstedt, Karolinska Institute, Sweden.\n", 80, 0);
  print_with_line_breaks("TEM-simulator is free software: you can redistribute it \
and/or modify it under the terms of the GNU General Public License as published by \
the Free Software Foundation, either version 3 of the License, or (at your option) \
any later version.\n", 80, 0);
  print_with_line_breaks("TEM-simulator is distributed in the hope that it will be \
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY \
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more \
details.\n", 80, 0);
  print_with_line_breaks("You should have received a copy of the GNU General Public \
License along with TEM-simulator.  If not, see <http://www.gnu.org/licenses/>.\n", 80, 0);
}

/****************************************************************************/

void print_help(int argc, char **argv){
  simulation *sim;
  param_table *pt;

  if(argc <= 2){
    printf("\n");
    print_with_line_breaks("TEM-simulator is run from the command line with \
one argument, which is the name of an input file:", 80, 0);
    printf("TEM-simulator <input file>\n\n");
    print_with_line_breaks("The input file is a text file which specifies a \
number of components to be used in the simulation, and for each component \
assigns values to a number of parameters. The input file has the following \
structure:", 80, 0);
    printf("=== <component type 1> <component name 1> ===\n");
    printf("<parameter 1a> = <value 1a>\n<parameter 1b> = <value 1b>\n");
    printf("=== <component type 2> <component name 2> ===\n");
    printf("<parameter 2a> = <value 2a>\n<parameter 2b> = <value 2b>\n\n");
    print_with_line_breaks("Lines starting with an = sign, indicate the \
beginning of a component, and the lines following it assign values to \
parameters in that component. Component names are optional, and have a use \
only for some component types. Lines starting with a # sign are treated as \
comments and ignored.\n", 80, 0);
    printf("The following component types can be used:\n  ");
    printf(TYPE_SIMULATION); printf("\n  ");
    printf(TYPE_SAMPLE); printf("\n  ");
    printf(TYPE_PARTICLE); printf("\n  ");
    printf(TYPE_PARTICLESET); printf("\n  ");
    printf(TYPE_GEOMETRY); printf("\n  ");
    printf(TYPE_ELECTRONBEAM); printf("\n  ");
    printf(TYPE_OPTICS); printf("\n  ");
    printf(TYPE_DETECTOR); printf("\n  ");
    printf(TYPE_VOLUME); printf("\n\n");
    printf("For help on a particular component, type\n");
    printf("TEM-simulator -help <component type>\n\n");
    printf("For help on a particular parameter, type\n");
    printf("TEM-simulator -help <component type> <parameter>\n\n");
    printf("See the manual for more details.\n\n");
    return;
  }

  sim = new_simulation();
  if(0 == strcmp("simulation", argv[2])){
    pt = sim->param;
  }
  else {
    pt = add_simcomp(sim, argv[2], "", 0);
  }

  if(pt == NULL){
    return;
  }

  if(argc == 3){
    print_comp_help(pt);
  }
  else {
    print_param_help(pt, argv[3]);
  }
  
  delete_simulation(sim);
}

/****************************************************************************/

int main(int argc, char **argv){
  simulation *sim;

  if(argc < 2){
    print_info();
    return 0;
  }

  if(0 == strcmp(argv[1], "-help")){
    print_help(argc, argv);
    return 0;
  }

  sim = new_simulation();
  if(read_input_data(sim, argv[1])){
    return 1;
  }

  simulation_create_log(sim);
  write_log_comment("Input data read from file %s.\n", argv[1]);
  simulation_write_log(sim);
  simulation_run(sim);
  delete_simulation(sim);

  return 0;

}
