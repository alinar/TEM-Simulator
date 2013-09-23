/*
 * Copyright 2008-2009, Hans Rullgard, Stockholm University and 
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

#ifndef LOG
#define LOG
#include <stdarg.h>

int create_log(const char *fn, int overwrite);

void finish_log(void);

int log_exists(void);

FILE *open_log(void);

void close_log(void);

void write_log(char *fmt, ...);

int write_log_comment(char *fmt, ...);

void write_log_disp(char *fmt, ...);

void write_log_comment_disp(char *fmt, ...);

#endif
