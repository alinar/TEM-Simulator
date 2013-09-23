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

#ifndef INPUT_HEADER
#define INPUT_HEADER
#include <stdio.h>
#include "macros.h"

/***********************************************************************
 * The param_table struct stores values of parameters read from a file,
 * for later use. When a parameter is read, it can optionally be checked
 * if the value belongs to an allowed range of values for that parameter.
 ***********************************************************************/

typedef char boolean;

typedef struct {
  char *name;     /* name of the parameter */
  char *descr;    /* description, used for help messages */
  char type;      /* type of the parameter: boolean, int, long, float, double, string */
  char set;       /* has parameter been set? */
  char req;       /* is parameter required? */
  char *def;      /* optional default value */
  char *allowed;  /* optional list of allowed values for string parameter */
  double minval;  /* optional minimum value of numerical parameter */
  double maxval;  /* optional maximum value of numerical parameter */
  void *addr;     /* pointer to address where value is stored */
} param;

typedef struct {
  param *params;  /* pointer to list of parameters */
  int length;     /* length of list (including unused parameters) */
  int used;       /* number of used parameters in list */
  char lock;      /* if the table is locked values can not be changed */
  char *name;     /* name of table (for log files and error messages) */
  char *class;    /* class of table (for log files and error messages) */
  char *descr;    /* description (for help messages) */
} param_table;


/***********************************************************************
 * Function:  new_param_table
 * Purpose:   Allocate memory for a new param_table with space for n
 *            parameters and return a pointer to it. Initially, none of 
 *            the parameters is used. The number of parameters can not 
 *            be changed once the table has been created.
 * Arguments: n - Maximum number of parameters that can be stored in table.
 *            class - The class of the table, used in log files and error 
 *                    messages. The intention is that all tables of the 
 *                    same class should contain the same parameters
 *                    (but they can be assigned different values).
 *            name - Name of the table, used for log files and error 
 *                   messages. The intention is that class and name 
 *                   together should be unique for each table.
 * Return:    Pointer to allocated table.
 */

param_table *new_param_table(int n, 
                             const char *class, 
                             const char *name);


/***********************************************************************
 * Function:  delete_param_table
 * Purpose:   Free memory allocated by new_param_table and for parameters
 *            inserted in the table.
 * Arguments: P - pointer to param_table to be deleted.
 */

void delete_param_table(param_table *P);


/***********************************************************************
 * Function:  param_table_set_lock
 * Purpose:   Set lock status of param_table. When table is locked, 
 *            functions which would otherwise change the value of 
 *            parameters have no effect and produce an error message.
 * Arguments: P - pointer to param_table.
 *            locked - new status of lock. If locked==0, the status is 
 *                     set to unlocked, otherwise it is set to locked.
 */

void param_table_set_lock(param_table *P, 
                          char locked);


/***********************************************************************
 * Function:  add_param
 * Purpose:   Add parameter to param_table. 
 * Arguments: P - pointer to param_table.
 *            name - name of parameter.
 *            type - type of parameter: "b" for boolean, "i" for int,
 *                   "l" for long, "f" for float, "d" for double, "s"
 *                   for string. Additionally, "s,<val1>,<val2>,<val3>"
 *                   means a string parameter with the allowed values
 *                   <val1>, <val2> or <val3>. The ',' can be replaced
 *                   by any other character.
 *            req - nonzero if parameter is required from input.
 *            min - minimum allowed value of numerical parameter.
 *            max - maximum allowed value of numerical parameter.
 * Return:    Pointer to added parameter.
 */

param *add_param(param_table *P, 
                 const char *name, 
                 const char *type, 
                 char req, 
                 double min, 
                 double max);


/***********************************************************************
 * Function:  add_param_req
 * Purpose:   Add required parameter without constraints to param_table. 
 * Arguments: P - pointer to param_table.
 *            name - name of parameter.
 *            type - type of parameter. See add_param.
 */

param *add_param_req(param_table *P, 
                     const char *name, 
                     const char *type);


/***********************************************************************
 * Function:  add_param_opt
 * Purpose:   Add optional parameter without constraints to param_table. 
 * Arguments: P - pointer to param_table.
 *            name - name of parameter.
 *            type - type of parameter. See add_param.
 */

param *add_param_opt(param_table *P, 
                     const char *name, 
                     const char *type);


/***********************************************************************
 * Function:  add_param_def
 * Purpose:   Add parameter with default value without constraints to 
 *            param_table. 
 * Arguments: P - pointer to param_table.
 *            name - name of parameter.
 *            type - type of parameter. See add_param.
 *            value - string representation of default value.
 */

param *add_param_def(param_table *P, 
                     const char *name, 
                     const char *type, 
                     const char *value);

param *add_param_req_constr(param_table *P, const char *name, const char *type, double minval, double maxval);

param *add_param_opt_constr(param_table *P, const char *name, const char *type, double minval, double maxval);

param *add_param_def_constr(param_table *P, const char *name, const char *type, const char *value, double minval, double maxval);

int set_param_char(param_table *pt, const char *name, char value);

int set_param_int(param_table *pt, const char *name, int value);

int set_param_long(param_table *pt, const char *name, long value);

int set_param_float(param_table *pt, const char *name, float value);

int set_param_double(param_table *pt, const char *name, double value);

int set_param_boolean(param_table *pt, const char *name, boolean value);

int set_param(param_table *pt, const char *name, const char *value);

int unset_param(param_table *pt, const char *name);

char get_param_char(param_table *pt, const char *name);

int get_param_int(param_table *pt, const char *name);

long get_param_long(param_table *pt, const char *name);

float get_param_float(param_table *pt, const char *name);

double get_param_double(param_table *pt, const char *name);

boolean get_param_boolean(param_table *pt, const char *name);

const char *get_param_string(param_table *pt, const char *name);

int set_comp_descr(param_table *pt, const char *descr);

int set_param_descr(param_table *pt, const char *name, const char *descr);

int print_comp_help(param_table *pt);

int print_param_help(param_table *pt, const char *name);

int check_params(param_table *P);

int require_param(param_table *P, const char *name);

int param_isset(param_table *P, const char *name);

FILE *open_input_file(const char *fn);

void close_input_file(void);

int read_input_block(param_table *P);

int skip_input_block(void);

const char *current_input_block_class(void);

const char *current_input_block_name(void);

int write_parameters_to_file(FILE *fp, param_table *P);

int write_parameters_to_log(param_table *P);

void write_block_sep_log(const char *name);

#endif
