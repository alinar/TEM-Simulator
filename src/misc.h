/*
 * Copyright 2008, Hans Rullgard, Stockholm University and 
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

#ifndef MISC_HEADER
#define MISC_HEADER
#include <stdio.h>

void err_message(char *fmt, ...);

char *strcpy_new(const char *s);

char *strncpy_term(char *dest, size_t destsize, const char *source, size_t n);

char *fgetline(char *s, int num, FILE *fp);

void print_with_line_breaks(const char *s, int linelength, int pos);

long min_l(long x, long y);

long max_l(long x, long y);

long round_to_int(double x);

double min_d(double x, double y);

double max_d(double x, double y);

double max3(double x1, double x2, double x3);

double min3(double x1, double x2, double x3);

long max_index(double *x, long n);

long min_index(double *x, long n);

int vsnprintf_win(char *str, size_t size, const char *format, va_list ap);

char *strtok_win(char *str, const char *delim);

#endif
