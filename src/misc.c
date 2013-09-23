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
#include <stdarg.h>
#include <string.h>
#include "misc.h"
#include "log.h"

/****************************************************************************/

void err_message(char *fmt, ...){
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
  exit(1);
}

/****************************************************************************/

char *strcpy_new(const char *s){
  size_t n = 1 + strlen(s);
  char *p = malloc(n*sizeof(char));
  STRNCPY(p, n, s, n);
  return p;
}

/****************************************************************************/

char *strncpy_term(char *dest, size_t destsize, const char *source, size_t n){
  if(n > destsize) n = destsize;
  STRNCPY(dest, destsize, source, n-1);
  dest[n-1] = '\0';
  return dest;
}

/****************************************************************************/

char *fgetline(char *s, int num, FILE *fp){
  char *end;
  s = fgets(s, num, fp);
  if(NULL == s) return NULL;
  if(feof(fp)) return s;
  /*  end = strchr(s, '\n'); */
  end = strpbrk(s, "\n\r");
  if(NULL == end){
    while('\n' != fgetc(fp));
  }
  else {
    *end = '\0';
  }
  return s;
}

/****************************************************************************/

void print_with_line_breaks(const char *s, int linelength, int pos){
  char *line;
  const char *curr;
  int i, n;
  if(linelength <= 0) return;
  line = malloc((linelength+1)*sizeof(char));
  curr = s;
  n = pos;
  while(curr[0] != '\0'){
    while(n < linelength && curr[0] != '\0'){
      i = 1 + strcspn(curr+1, " \0");
      if(n == 0 && i > linelength){
        i = linelength;
      }
      n += i;
      if(n <= linelength){
        strncpy_term(line, linelength+1, curr, i+1);
        printf("%s", line);
        curr += i;
      }
    }
    printf("\n");
    curr += strspn(curr, " ");
    n = 0;
  }
  free(line);
}

/****************************************************************************/

long min_l(long x, long y){
  return  x>y?y:x;
}

/****************************************************************************/

long max_l(long x, long y){
  return  x>y?x:y;
}

/****************************************************************************/

long round_to_int(double x){
  return (long)floor(x + 0.5);
}

/****************************************************************************/

double min_d(double x, double y){
  return  x>y?y:x;
}

/****************************************************************************/

double max_d(double x, double y){
  return  x>y?x:y;
}

/****************************************************************************/

double max3(double x1, double x2, double x3){
  double x = x1>x2?x1:x2;
  return x>x3?x:x3;
}

/****************************************************************************/

double min3(double x1, double x2, double x3){
  double x = x1<x2?x1:x2;
  return x<x3?x:x3;
}

/****************************************************************************/

long max_index(double *x, long n){
  long i = 0, j;
  if(n <= 0){
    return 0;
  }
  for(j = 1; j < n; j++){
    if(x[j] > x[i]){
      i = j;
    }
  }
  return i;
}

/****************************************************************************/

long min_index(double *x, long n){
  long i = 0, j;
  if(n <= 0){
    return 0;
  }
  for(j = 1; j < n; j++){
    if(x[j] < x[i]){
      i = j;
    }
  }
  return i;
}

/****************************************************************************/

#ifdef USE_WINDOWS_MACROS

int vsnprintf_win(char *str, size_t size, const char *format, va_list ap){
  int r = vsnprintf_s(str, size, _TRUNCATE, format, ap);
  str[size-1] = '\0';
  return r;
}

char *strtok_win(char *str, const char *delim){
  static char *p = NULL;
  return strtok_s(str, delim, &p);
}

#endif
