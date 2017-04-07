/**************************************************************************
*
* Cinderella -- COMPASS ONLINE FILTER
*
* Copyright (C) 2002-2006
*   Technische Universititaet Muenchen
*   Physik-Department
*   James Frank Strasse
*   D-85748 Garching
*   Germany
*
* Author(s) of this file: RoK
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


#ifndef _CHIPSADC_HUFFMAN_H
#define _CHIPSADC_HUFFMAN_H

#include "config.h"

namespace CS {
  
#define WORD_LEN 32u

/** huffman tree */
struct hufftree_s {
  int                 val;
  uint64              freq;
  struct hufftree_s  *next[2];
};
typedef struct hufftree_s hufftree_t;

struct elt_s {
  int                 val;
  uint64              freq;
};
typedef struct elt_s elt_t;

/** structure for huffman statistics */
struct huffstat_s {
  elt_t              *elt;
#if 0
  int                *val;
                        /**< array of values */
  uint64             *freq;
                        /**< frequency of the values */
#endif
  unsigned int        size;     /**< size of both arrays */
};
typedef struct huffstat_s huffstat_t;

/** lookup table for huffman coding */
struct hufflookup_s {
  unsigned int        val;              /**< the original value */
  uint32              huf;
                        /**< the huffman code */
  unsigned int        length;           /**< the length of the huffman code */
};
typedef struct hufflookup_s hufflookup_t;


int                 sadc_huffman_encode(const hufflookup_t *lookup,
                                        const unsigned int lookup_size,
                                        const unsigned int *input,
                                        const unsigned int input_size,
                                        const unsigned int samplenum,
                                        uint32 **output);

int                 sadc_huffman_decode(const hufftree_t *hufftree,
                                        const unsigned int *input,
                                        const unsigned int input_size,
                                        unsigned int *sample_num,
                                        unsigned int **output);

int                 sanity_check(const hufftree_t *hufftree,
                                 const hufflookup_t *hufflookup,
                                 const unsigned lookup_size);

uint32              make_int_from_code(const char *code);

int                 create_lookup_from_string(const char *string,
                                              hufflookup_t **hufflookup);

int                 create_lookup_from_file(const char *path,
                                            hufflookup_t **hufflookup);

int                 val_code(const int val, unsigned int level,
                             const hufftree_t *hufftree, char *code);

int                 copy_tree_members(const hufftree_t *orig_hufftree,
                                      hufftree_t **new_hufftree);

hufftree_t         *copy_tree(const hufftree_t *hufftree);

hufftree_t         *hufftree_read_from_string(const char *string);

hufftree_t         *hufftree_read(const char *path);

int                 hufftree_dump(const char *path,
                                  const hufftree_t *hufftree,
                                  const huffstat_t *huffstat);

int                 huffstat_dump(const char *path, huffstat_t *huffstat);

huffstat_t         *huffstat_read(const char *path);

void                make_allstat(huffstat_t *huffstat, const int range_begin,
                                 const int range_end);

void                insert_value(huffstat_t *huffstat, int value);

void                sort_huffstat(huffstat_t *huffstat);

void                quick_sort_huffstat(huffstat_t *huffstat, unsigned al,
                                        unsigned ar);

hufftree_t         *make_hufftree(huffstat_t *huffstat);

void                calc_bit_per_val(hufftree_t *hufftree,
                                     const unsigned int level, float *bits,
                                     uint64 *freq);

void                free_tree(hufftree_t *hufftree);

int                 compar_huffstat(const void *e1, const void *e2);

} //namespace CS
#endif
