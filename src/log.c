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
#include <stdarg.h>
#include <string.h>
#include "log.h"
#include "macros.h"
#include "misc.h"

/****************************************************************************/

static char *logfile = NULL;
static int logfile_active = 0;
static FILE *logfile_ptr = NULL;

/****************************************************************************/

int write_log_comment_valist(int length, char *fmt, va_list ap);

/****************************************************************************/

int create_log(const char *fn, int overwrite){
  FILE *fp;
  finish_log();
  if(!overwrite){
    FOPEN(fp, fn, "a");
    if(fp != NULL && ftell(fp) > 0){
      fclose(fp); /* File exists and is nonempty */
      return 1;
    }
  }
  else {
    FOPEN(fp, fn, "w");
  }
  if(fp != NULL){
    fclose(fp);
    logfile = strcpy_new(fn);
    logfile_active = 1;
    return 0;
  }
  return 2;
}

/****************************************************************************/

void finish_log(void){
  close_log();
  if(logfile != NULL){
    free(logfile);
  }
  logfile = NULL;
  logfile_active = 0;
  logfile_ptr = NULL;
}

/****************************************************************************/

int log_exists(void){
  return logfile_active;
}

/****************************************************************************/

FILE *open_log(void){
  if(logfile_active && (logfile_ptr == NULL)){
    FOPEN(logfile_ptr, logfile, "a");
  }
  return logfile_ptr;
}

/****************************************************************************/

void close_log(void){
  if(logfile_ptr != NULL){
    fclose(logfile_ptr);
    logfile_ptr = NULL;
  }
}

/****************************************************************************/

void write_log(char *fmt, ...){
  FILE *fp;
  va_list ap;
  fp = open_log();
  if(fp != NULL){
    va_start(ap, fmt);
    vfprintf(fp, fmt, ap);
    va_end(ap);
    close_log();
  }
}

/****************************************************************************/

int write_log_comment(char *fmt, ...){
  int n, ret;
  va_list ap;
  va_start(ap, fmt);
  n = VSNPRINTFCOUNT(fmt, ap);
  va_end(ap);
  va_start(ap, fmt);
  ret = write_log_comment_valist(n, fmt, ap);
  va_end(ap);
  return ret;
}

/****************************************************************************/

int write_log_comment_valist(int length, char *fmt, va_list ap){
  int i, ret = 0;
  char *str;
  FILE *fp = open_log();
  if(fp == NULL) return 1;
  length++;
  str = malloc(length * sizeof(char));
  i = VSNPRINTF(str, length, fmt, ap);
  if((length <= i) || (i < 0)) ret = 1; /* Both conditions are needed to handle Unix and Windows versions */
  length = strlen(str);
  fputc(COMMENT_CHAR, fp);
  fputc(' ', fp);
  for(i = 0; i < length-1; i++){
    if(fputc(str[i], fp) == '\n'){
      fputc(COMMENT_CHAR, fp);
      fputc(' ', fp);
    }
  }
  fputc(str[length-1], fp);
  close_log();
  free(str);
  return ret;
}

/****************************************************************************/

void write_log_disp(char *fmt, ...){
  FILE *fp;
  va_list ap;
  fp = open_log();
  if(fp != NULL){
    va_start(ap, fmt);
    vfprintf(fp, fmt, ap);
    va_end(ap);
    close_log();
  }
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
}

/****************************************************************************/

void write_log_comment_disp(char *fmt, ...){
  va_list ap;
  int n;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  if(log_exists()){
    va_start(ap, fmt);
    n = VSNPRINTFCOUNT(fmt, ap);
    va_end(ap);
    va_start(ap, fmt);
    write_log_comment_valist(n, fmt, ap);
    va_end(ap);
  }
}
