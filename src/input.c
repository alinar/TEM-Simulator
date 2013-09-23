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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "input.h"
#include "log.h"
#include "macros.h"
#include "misc.h"


static FILE *current_input_file = NULL;
static char input_file_name[MAX_LINE_LENGTH+1];
static char current_block_class[MAX_NAME_LENGTH+1] = "none";
static char current_block_name[MAX_NAME_LENGTH+1] = "";
static int current_line_number;

/****************************************************************************/

param *add_param(param_table *P, const char *name, const char *type, char req, double min, double max);

void set_block_name(char *line);

int set_this_param(param *p, const char *value);

int set_param_line(param_table *pt, char *s);

param *find_param(param_table *pt, const char *name);

int check_param(param *p, const char *ptclass, const char *ptname);

/****************************************************************************/

param_table *new_param_table(int n, const char *class, const char *name){
  param_table *p;
  p = malloc(sizeof(param_table));
  p->params = malloc(n*sizeof(param));
  p->used = 0;
  p->length = n;
  p->lock = 0;
  p->name = strcpy_new(name);
  p->class = strcpy_new(class);
  p->descr = NULL;
  return p;
}

/****************************************************************************/

void delete_param_table(param_table *P){
  int i;
  for(i = 0; i < P->used; i++){
    free(P->params[i].name);
    if(P->params[i].descr != NULL) free(P->params[i].descr);
    if(P->params[i].addr != NULL) free(P->params[i].addr);
    if(P->params[i].allowed != NULL) free(P->params[i].allowed);
  }
  free(P->name);
  free(P->class);
  if(P->descr != NULL) free(P->descr);
  free(P->params);
  free(P);
}

/****************************************************************************/

void param_table_set_lock(param_table *P, char locked){
  P->lock = locked;
}

/****************************************************************************/

param *add_param(param_table *P, const char *name, const char *type, char req, double minval, double maxval){
  int n;
  if(P->used >= P->length){
    WARNING("Parameter table is full\n");
    return NULL;
  }
  if(0 == strlen(type)){
    WARNING("No type specified for parameter %s:%s.\n", P->class, name);
    return NULL;
  }
  else {
    n = P->used;
    P->used++;
    P->params[n].name = strcpy_new(name);
    P->params[n].descr = NULL;
    P->params[n].type = type[0];
    P->params[n].set = 0;
    P->params[n].req = req;
    P->params[n].def = NULL;
    P->params[n].allowed = NULL;
    P->params[n].minval = minval;
    P->params[n].maxval = maxval;
    P->params[n].addr = NULL;
    if(1 < strlen(type)) P->params[n].allowed = strcpy_new(type+1);
    return &(P->params[n]);
  }
}

/****************************************************************************/

param *add_param_req(param_table *P, const char *name, const char *type){
  return add_param(P, name, type, 1, -HUGE_VAL, HUGE_VAL);
}

/****************************************************************************/

param *add_param_opt(param_table *P, const char *name, const char *type){
  return add_param(P, name, type, 0, -HUGE_VAL, HUGE_VAL);
}

/****************************************************************************/

param *add_param_def(param_table *P, const char *name, const char *type, const char *value){
  param *par = add_param(P, name, type, 1, -HUGE_VAL, HUGE_VAL);
  if(par != NULL){
    set_this_param(par, value);
    par->def = strcpy_new(value);
    par->set = 2; /* Indicates that parameter is set to default value */
  }
  return par;
}

/****************************************************************************/

param *add_param_req_constr(param_table *P, const char *name, const char *type, double minval, double maxval){
  return add_param(P, name, type, 1, minval, maxval);
}

/****************************************************************************/

param *add_param_opt_constr(param_table *P, const char *name, const char *type, double minval, double maxval){
  return add_param(P, name, type, 0, minval, maxval);
}

/****************************************************************************/

param *add_param_def_constr(param_table *P, const char *name, const char *type, const char *value, double minval, double maxval){
  param *par = add_param(P, name, type, 1, minval, maxval);
  if(par != NULL){
    set_this_param(par, value);
    par->def = strcpy_new(value);
    par->set = 2; /* Indicates that parameter is set to default value */
  }
  return par;
}

/****************************************************************************/

int set_this_param(param *p, const char *value){
  int set = 1;
  switch(p->type){
  case 'c':
    if(p->addr == NULL) p->addr = malloc(sizeof(char));
    *(char*)p->addr = value[0]; 
    break;
  case 'i':
    if(p->addr == NULL) p->addr = malloc(sizeof(int));
    *(int*)p->addr = atoi(value);
    break;
  case 'l':
    if(p->addr == NULL) p->addr = malloc(sizeof(long));
    *(long*)p->addr = atol(value);
    break;
  case 'f':
    if(p->addr == NULL) p->addr = malloc(sizeof(float));
    *(float*)p->addr = (float)strtod(value, NULL);
    break;
  case 'd':
    if(p->addr == NULL) p->addr = malloc(sizeof(double));
    *(double*)p->addr = strtod(value, NULL);
    break;
  case 's':
    if(p->addr != NULL) free(p->addr);
    p->addr =  strcpy_new(value);
    break;
  case 'b':
    if(0 == strcmp(value, YES_STRING)){
      if(p->addr == NULL) p->addr = malloc(sizeof(boolean));
      *(boolean*)p->addr = 1;
    }
    else if(0 == strcmp(value, NO_STRING)){
      if(p->addr == NULL) p->addr = malloc(sizeof(boolean));
      *(boolean*)p->addr = 0;
    }
    else {
      WARNING("Warning: Boolean variable %s must be %s or %s.\n", p->name, YES_STRING, NO_STRING);
      set = 0;
    }
    break;
  default:
    WARNING("Warning: Undefined type %c\n", p->type);
    set = 0;
  }
  if(set){
    p->set = 1;
    return 0;
  }
  return 1;
}

/****************************************************************************/

int set_param_char(param_table *pt, const char *name, char value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('c' != p->type){
    WARNING("Parameter %s:%s is not of type char\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(char));
  *(char*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param_int(param_table *pt, const char *name, int value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('i' != p->type){
    WARNING("Parameter %s:%s is not of type int\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(int));
  *(int*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param_long(param_table *pt, const char *name, long value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('l' != p->type){
    WARNING("Parameter %s:%s is not of type long\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(long));
  *(long*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param_float(param_table *pt, const char *name, float value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('f' != p->type){
    WARNING("Parameter %s:%s is not of type float\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(float));
  *(float*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param_double(param_table *pt, const char *name, double value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('d' != p->type){
    WARNING("Parameter %s:%s is not of type double\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(double));
  *(double*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param_boolean(param_table *pt, const char *name, boolean value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if('b' != p->type){
    WARNING("Parameter %s:%s is not of type boolean\n", pt->class, p->name);
    return 1;
  }
  if((value != 0)&&(value != 1)){
    WARNING("Boolean parameter %s:%s can only be set to 0 or 1\n", pt->class, p->name);
    return 1;
  }
  if(p->addr == NULL) p->addr = malloc(sizeof(double));
  *(boolean*)p->addr = value;
  p->set = 1;
  return 0;
}

/****************************************************************************/

int set_param(param_table *pt, const char *name, const char *value){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  return set_this_param(p, value);
}

/****************************************************************************/

int unset_param(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(pt->lock){
    WARNING("Parameter table %s %s is locked\n", pt->class, pt->name);
    return 1;
  }
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if(p->def != NULL){
    set_this_param(p, p->def);
    p->set = 2;
  }
  else {
    p->set = 0;
    if(p->addr != NULL){
      free(p->addr);
    }
  }
  return 0;
}

/****************************************************************************/

char get_param_char(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('c' != p->type){
    WARNING("Parameter %s:%s is not of type char\n", pt->class, p->name);
    return 0;
  }
  return *(char*)p->addr;
}

/****************************************************************************/

int get_param_int(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('i' != p->type){
    WARNING("Parameter %s:%s is not of type int\n", pt->class, p->name);
    return 0;
  }
  return *(int*)p->addr;
}

/****************************************************************************/

long get_param_long(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('l' != p->type){
    WARNING("Parameter %s:%s is not of type long\n", pt->class, p->name);
    return 0;
  }
  return *(long*)p->addr;
}

/****************************************************************************/

float get_param_float(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('f' != p->type){
    WARNING("Parameter %s:%s is not of type float\n", pt->class, p->name);
    return 0;
  }
  return *(float*)p->addr;
}

/****************************************************************************/

double get_param_double(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('d' != p->type){
    WARNING("Parameter %s:%s is not of type double\n", pt->class, p->name);
    return 0;
  }
  return *(double*)p->addr;
}

/****************************************************************************/

boolean get_param_boolean(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('b' != p->type){
    WARNING("Parameter %s:%s is not of type boolean\n", pt->class, p->name);
    return 0;
  }
  return *(boolean*)p->addr;
}

/****************************************************************************/

const char *get_param_string(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  if('s' != p->type){
    WARNING("Parameter %s:%s is not of type string\n", pt->class, p->name);
    return 0;
  }
  return (char*)p->addr;
}

/****************************************************************************/

int set_comp_descr(param_table *pt, const char *descr){
  if(pt->descr != NULL) free(pt->descr);
  pt->descr = strcpy_new(descr);
  return 0;
}

/****************************************************************************/

int set_param_descr(param_table *pt, const char *name, const char *descr){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  if(p->descr != NULL) free(p->descr);
  p->descr = strcpy_new(descr);
  return 0;
}

/****************************************************************************/

int print_comp_help(param_table *pt){
  int i;
  printf("\nComponent: %s\n\n", pt->class);
  if(pt->descr != NULL){
    printf("Description: ");
    print_with_line_breaks(pt->descr, 80, 13);
    printf("\n");
  }
  printf("For help on a specific parameter, type\n");
  printf("TEM-simulator -help %s <parameter>\n\n", pt->class);
  printf("Supported parameters:\n");
  for(i = 0; i < pt->used; i++){
    printf("  %s\n", pt->params[i].name);
  }
  printf("\n");
  return 0;
}

/****************************************************************************/

int print_param_help(param_table *pt, const char *name){
  char *allowed, *substr;
  char delim[2];
  int first;
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 1;
  }
  printf("\nComponent: %s\n\nParameter: %s\n\n", pt->class, p->name);
  if(p->descr != NULL){
    printf("Description: ");
    print_with_line_breaks(p->descr, 80, 13);
    printf("\n");
  }
  printf("Allowed values: ");
  switch(p->type){
  case 'c':
    printf("any character\n\n");
    break;
  case 'i':
  case 'l':
    if(p->minval == -HUGE_VAL && p->maxval == HUGE_VAL){
      printf("any integer\n\n");
    }
    else {
      printf("any integer between %1.0f and %1.0f\n\n", (float)p->minval, (float)p->maxval);
    }
    break;
  case 'f':
  case 'd':
    if(p->minval == -HUGE_VAL && p->maxval == HUGE_VAL){
      printf("any real number\n\n");
    }
    else {
      printf("any real number between %1.2f and %1.2f\n\n", (float)p->minval, (float)p->maxval);
    }
    break;
  case 'b':
    printf("%s or %s\n\n", YES_STRING, NO_STRING);
    break;
  case 's':
    if(NULL == p->allowed){
      printf("any string\n\n");
    }
    else {
      printf("one of the strings ");
      delim[0] = p->allowed[0]; delim[1] = '\0'; /* First character in string is delimiter */
      allowed = strcpy_new(p->allowed+1);
      substr = STRTOK(allowed, delim);
      first = 1;
      while(NULL != substr){
        if(!first) printf(", ");
        printf("\"%s\"", substr);
        first = 0;
	substr = STRTOK(NULL, delim);
      }
      printf("\n\n");
      free(allowed);
    }
  }
  if(p->def != NULL){
    printf("Default: ");
    if('s' == p->type){
      printf("\"%s\"\n\n", p->def);
    }
    else {
      printf("%s\n\n", p->def);
    }
  }
  return 0;
}

/****************************************************************************/

int set_param_line(param_table *pt, char *s){
  char *name, *value, *c;
  int k;
  name = s;
  c = strchr(s, '=');
  if(NULL == c){
    WARNING("Unable to parse input line %i\n", current_line_number);
    return 1;
  }
  *c = '\0'; /* Terminate the parameter name part of the line */
  value = c + 1;
  name = STRTOK(name, " \t");
  if(NULL == name){
    WARNING("Unable to parse input line %i\n", current_line_number);
    return 1;
  }
  k = strspn(value, " \t"); /* Skip initial blank spaces after the equals sign */
  value += k;
  if('\0' == value[0]){
    WARNING("Unable to parse input line %i\n", current_line_number);
    return 1;
  }
  if('"' == value[0]){ /* quoted string */
    value++;
    c = strchr(value, '"'); /* Find the matching end quote */
    if(NULL == c){
      WARNING("Unable to parse input line %i\n", current_line_number);
      return 1;
    }
    *c = '\0';
  }
  else {
    value = STRTOK(value, " \t");
    if(NULL == c){
      WARNING("Unable to parse input line %i\n", current_line_number);
      return 1;
    }
  }
  return set_param(pt, name, value);
}

/****************************************************************************/

param *find_param(param_table *pt, const char *name){
  int i = 0;
  while((i < pt->used) && strcmp(name, pt->params[i].name)){i++;}
  if(i < pt->used){
    return &(pt->params[i]);
  }
  else {
    return NULL;
  }
}

/****************************************************************************/

int check_params(param_table *P){
  int i, j = 0;
  for(i = 0; i < P->used; i++){
    if(P->params[i].req || P->params[i].set){
      if(check_param(&(P->params[i]), P->class, P->name)) j = 1;
    }
  }
  return j;
}

/****************************************************************************/

int check_param(param *p, const char *ptclass, const char *ptname){
  char delim[2];
  char *substr;
  char *allowed;
  int check_bounds;
  double val;
  int ret = 0;
  if(!p->set){
    WARNING("Parameter %s undefined in %s %s\n", p->name, ptclass, ptname);
    return 1;
  }
  check_bounds = 1;
  val = 0;
  switch(p->type){
  case 'i':
    val = (double)*(int*)p->addr; break;
  case 'l':
    val = (double)*(long*)p->addr; break;
  case 'f':
    val = (double)*(float*)p->addr; break;
  case 'd':
    val = *(double*)p->addr; break;
  default:
    check_bounds = 0;
  }
  if(check_bounds && ((val < p->minval) || (val > p->maxval))){
    WARNING("Parameter %s in %s %s has illegal value %f. Value must be in the interval [%f,%f].\n", 
	    p->name, ptclass, ptname, (float)val, (float)p->minval, (float)p->maxval);
    ret = 1;
  }
  if('s' == p->type && NULL != p->allowed){
    delim[0] = p->allowed[0]; delim[1] = '\0'; /* First character in string is delimiter */
    allowed = strcpy_new(p->allowed+1);
    substr = STRTOK(allowed, delim);
    while((NULL != substr) && strcmp(substr, (char*)p->addr)){
      substr = STRTOK(NULL, delim);
    }
    free(allowed);
    if(NULL == substr){
      WARNING("Parameter %s in %s %s has illegal value %s. Permissible values are %s.\n", 
              p->name, ptclass, ptname, (char*)p->addr, p->allowed+1);
      ret = 1;
    }
  }
  return ret;
}

/****************************************************************************/

int require_param(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s required\n", pt->class, name);
    return 1;
  }
  return check_param(p, pt->class, pt->name);
}

/****************************************************************************/

int param_isset(param_table *pt, const char *name){
  param *p = find_param(pt, name);
  if(p == NULL){
    WARNING("Unknown parameter %s:%s\n", pt->class, name);
    return 0;
  }
  else {
    return p->set;
  }
}

/****************************************************************************/

FILE *open_input_file(const char *fn){
  char line[MAX_LINE_LENGTH+1];
  close_input_file();
  strncpy_term(input_file_name, MAX_LINE_LENGTH, fn, MAX_LINE_LENGTH);
  FOPEN(current_input_file, fn, "r");
  current_line_number = 0;
  if(current_input_file != NULL){
    while(!feof(current_input_file)){
      current_line_number++;
      fgetline(line, MAX_LINE_LENGTH+1, current_input_file);
      if(line[0] == '='){
        break;
      }
    }
    if(!feof(current_input_file)){
      set_block_name(line);
    }
    else {
      close_input_file();
    }
  }
  return current_input_file;
}

/****************************************************************************/

void close_input_file(void){
  if(current_input_file != NULL){
    fclose(current_input_file);
  }
  current_input_file = NULL;
  strncpy_term(current_block_class, MAX_NAME_LENGTH, "end", 4);
  strncpy_term(current_block_name, MAX_NAME_LENGTH, "", 1);
}

/****************************************************************************/

int read_input_block(param_table *P){
  char line[MAX_LINE_LENGTH+1];
  if(P->lock){
    WARNING("Parameter table %s %s is locked.\n", P->class, P->name);
    return 1;
  }
  if(NULL == current_input_file){
    WARNING("No input file open.\n");
    return 1;
  }
  while(!feof(current_input_file)){
    current_line_number++;
    if(NULL == fgetline(line, MAX_LINE_LENGTH+1, current_input_file)){
      break;
    }
    if((1 < strlen(line))&&(line[0] != COMMENT_CHAR)){
      if(line[0] == '='){
	break;
      }
      set_param_line(P, line);
    }
  }
  if(line[0] == '='){
    set_block_name(line);
  }
  else {
    close_input_file();
  }
  return 0;
}

/****************************************************************************/

int skip_input_block(void){
  char line[MAX_LINE_LENGTH+1];
  if(NULL == current_input_file){
    WARNING("No input file open.\n");
    return 1;
  }
  while(!feof(current_input_file)){
    current_line_number++;
    if(NULL == fgetline(line, MAX_LINE_LENGTH+1, current_input_file)){
      break;
    }
    if(line[0] == '='){
      break;
    }
  }
  if(line[0] == '='){
    set_block_name(line);
  }
  else {
    close_input_file();
  }
  return 0;
}

/****************************************************************************/

const char *current_input_block_class(void){
  return current_block_class;
}

/****************************************************************************/

const char *current_input_block_name(void){
  return current_block_name;
}

/****************************************************************************/

void set_block_name(char *line){
  char *tok = STRTOK(line, "= \n");
  if(tok != NULL){
    strncpy_term(current_block_class, MAX_NAME_LENGTH, tok, MAX_NAME_LENGTH);
    tok = STRTOK(NULL, "= \n");
    if(tok != NULL){
      strncpy_term(current_block_name, MAX_NAME_LENGTH, tok, MAX_NAME_LENGTH);
    }
    else {
      strncpy_term(current_block_name, MAX_NAME_LENGTH, "", 1);
    }
  }
  else {
    strncpy_term(current_block_class, MAX_NAME_LENGTH, "none", 5);
    strncpy_term(current_block_name, MAX_NAME_LENGTH, "", 1);
  }
}

/****************************************************************************/

int write_parameters_to_file(FILE *fp, param_table *P){
  int i;
  size_t k;
  char *val;
  if(fp == NULL){
    return 1;
  }
  for(i = 0; i < P->used; i++){
    if(1 != P->params[i].set){
      fprintf(fp, "%c ", COMMENT_CHAR); /* Parameter is undefined or has default value */
    }
    else {
      fprintf(fp, "  ");
    }
    if(0 == P->params[i].set){
      fprintf(fp, "%-25s   undefined\n", P->params[i].name);
    }
    else {
      switch(P->params[i].type){
      case 'c':
	fprintf(fp, "%-25s = %c", P->params[i].name, *(char*)(P->params[i].addr)); break;
      case 'i':
	fprintf(fp, "%-25s = %i", P->params[i].name, *(int*)(P->params[i].addr)); break;
      case 'l':
	fprintf(fp, "%-25s = %i", P->params[i].name, (int)*(long*)(P->params[i].addr)); break;
      case 'f':
	fprintf(fp, "%-25s = %f", P->params[i].name, *(float*)(P->params[i].addr)); break;
      case 'd':
	fprintf(fp, "%-25s = %f", P->params[i].name, (float)*(double*)(P->params[i].addr)); break;
      case 's':
	val = (char*)(P->params[i].addr);
	k = strcspn(val, " \t"); /* Check if value contains spaces */
	if(k < strlen(val)){
	  fprintf(fp, "%-25s = \"%s\"", P->params[i].name, val);
	}
	else {
	  fprintf(fp, "%-25s = %s", P->params[i].name, val);
	}
	break;
      case 'b':
	if(*(char*)P->params[i].addr == 1){
	  fprintf(fp, "%-25s = %s", P->params[i].name, YES_STRING);
	}
	else if(*(char*)P->params[i].addr == 0){
	  fprintf(fp, "%-25s = %s", P->params[i].name, NO_STRING);
	}
	else {
	  fprintf(fp, "%c Boolean %s has undefined value", COMMENT_CHAR, P->params[i].name);
	}
	break;
      }
      if(2 == P->params[i].set){
	fprintf(fp, "  (default)\n");
      }
      else {
	fprintf(fp, "\n");
      }
    }
  }
  return 0;
}

/****************************************************************************/

int write_parameters_to_log(param_table *P){
  int i;
  FILE *fp;
  write_log("%c\n=== %s %s ===\n%c\n", COMMENT_CHAR, P->class, P->name, COMMENT_CHAR);
  fp = open_log();
  i = write_parameters_to_file(fp, P);
  close_log();
  return i;
}

/****************************************************************************/

void write_block_sep_log(const char *name){
  write_log("%c\n=== %s ===\n%c\n", COMMENT_CHAR, name, COMMENT_CHAR);
}
