
/**************************************************************************
*
* Cinderella -- COMPASS ONLINE FILTER
*
* Copyright (C) 2002-2004
*   Technische Universititaet Muenchen
*   Physik-Department
*   James Frank Strasse
*   D-85748 Garching
*   Germany
*
* Author(s) of this file: RK, TN
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
**************************************************************endofheader*/

/** \file
    This is the compatibility layer to allow CORAL to use Cinderella's
    silicon_timing.  To separate Silicon Time Reconstruction from Cinderella,
    copy silicon_timing.c, silicon_timing_types.h and silicon_timing_coral.h
    to the desired destination, renaming silicon_timing_coral.h to
    silicon_timing.h in the process.  To use the time reconstruction
    functions, silicon_timing.h needs to be included.  High-level functions
    that may be re-used are silicon_calculate_timings() and
    silicon_clusterize_plane() whereas interesting low-level functions are
    ratio_error(), single_timing() and join_timing(). */

#ifndef _SILICON_TIMING_H
#define _SILICON_TIMING_H

#include <cmath>
#include <cstdlib>
#include <iostream>


/* from messages.h */
#define _ERROR 0
#define _DEBUG3 0

#define MSG(LVL, SYS, ...) do {} while (0)
#define MSG_T(LVL, SYS, ...) do {} while (0)


/* from filter_malloc.h */
#define test_malloc_failed(ptr) do {            \
  if(!(ptr)) {                                  \
    printf("\nmalloc failed! aborting\n");      \
    exit(EXIT_FAILURE);                         \
  }                                             \
} while(0)


/** calloc() replacement: aborts if out of memory */
#define F_CALLOC(var, type, count) do {    \
  var = (type*)calloc(count, sizeof(type));     \
  test_malloc_failed(var);                      \
} while(0)


/* from histogram{,2}.h */
#define histogram_fill(...) do {} while (0)
#define histogram2_fill(...) do {} while (0)
#define histogram_fill_raw(...) do {} while (0)


/* from filter_util.h */
#define SQR(a) ((a)*(a))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MEAN(a,b) (((a)+(b))/2)

int static          inline
cind_iroundf(float f)
{
  return lroundf(ceil(f - 0.5));
}

static inline int 
cind_iround(double f)
{
  return lround(ceil(f - .5));
}

#ifdef __GNUC__
  #if ___GNUC__ >= 3
    #define __A_CONST__ __attribute__((__const__))
    #define __A_PURE__ __attribute__((__pure__))
  #else
    #define __A_CONST__
    #define __A_PURE__
  #endif

#else
  #define __A_CONST__
  #define __A_PURE__
#endif

#include "silicon_timing_types.h"

#undef __A_CONST__
#undef __A_PURE__

#endif
