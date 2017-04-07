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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ChipSADC_huffman.h"


namespace CS {

inline void
free_not_null(void *mem)
{
  if (mem)
    free(mem);
}
  
/** 
 * encodes a set of values from the 'input' array by using the huffman lookup
 * table 'lookup' and returns a pointer to the encoded data array to 'output' array. 
 * The first 9 bit contain the nr. of samples-1 per channel.
 * They are followed by multiple blocks of 6 bit channel number trailed by 
 * 'samplenum' huffman encoded samples.
 * The size of the array is the return value of the function. 
 * Negative return values indicate errors 
 **/
int
sadc_huffman_encode(const hufflookup_t *lookup,
                    const unsigned int lookup_size, const unsigned int *input,
                    const unsigned int input_size,
                    const unsigned int samplenum, uint32 **output)
{
  unsigned int        i;
  unsigned int        size = 0;
  unsigned int        cur_size = (input_size / 6) + 1;  /* size of output
                                                           buffer */

  /* first 9 bit will be 'samplenum - 1' */
  /* the next 4 bit will be left free for channel nr */
  unsigned int        pos = 15;  /* the actual bit position inside the current word */
  uint32             *coded = NULL;

  coded = (uint32 *)calloc(cur_size, sizeof(unsigned int));

  memset(coded, 0, sizeof(uint32) * cur_size);

  /* first 9 bit must be 'samplenum - 1' and must not be encoded! */
  coded[0] = (samplenum - 1) & 511;
  /* next 6 bit must be channel nr and must not be encoded! */
  coded[0] += input[0] << 9;
  pos = 15;

  for (i = 1; i < input_size; i++) {

    if (i % (samplenum + 1) == 0) {
      if (pos + 6 > WORD_LEN - 1) {
        coded[size] += input[i] << pos;
        size++;
        if (size >= cur_size) {

          cur_size += 2;
          {
            uint32             *tmp;
            tmp = (uint32 *)realloc(coded, cur_size * sizeof(uint32));
            if (!tmp) {
              free(coded);
              return -2;
            }
            coded = tmp;
          }
          memset(&coded[size], 0, (cur_size - size) * sizeof(uint32));
        }
        coded[size] += input[i] >> (WORD_LEN - pos);
        pos = pos + 6 - WORD_LEN;
      } else {
        coded[size] += input[i] << pos;
        pos += 6;
      }
      continue;
    }

    if (input[i] >= lookup_size) {
      free_not_null(coded);
      return -1;
    }

    /* check if enough bit are available in current word */
    if (pos + lookup[input[i]].length > WORD_LEN - 1) {
      coded[size] += (lookup[input[i]].huf << pos);
      size++;
      if (size >= cur_size) {
        cur_size += 2;
        {
          uint32             *tmp;
          tmp = (uint32 *)realloc(coded, cur_size * sizeof(uint32));
          if (!tmp) {
            free(coded);
            return -2;
          }
          coded = tmp;
        }
        memset(&coded[size], 0, (cur_size - size) * sizeof(uint32));
      }
      coded[size] += (lookup[input[i]].huf >> (WORD_LEN - pos));
      pos = pos + lookup[input[i]].length - WORD_LEN;
    } else {
      coded[size] += (lookup[input[i]].huf << pos);
      pos += lookup[input[i]].length;
    }

  }

  *output = coded;
  return (int)size + 1;
}

/** 
 * decodes huffman encoded data from the 'input' buffer using huffman tree 
 * 'hufftree'.  The function returns the size of the 'output' array on
 * success, and <0 on error.  The decoded data (multiple blocks of channel
 * number followed by samples) is returned in the array 'output'.  In
 * 'sample_num' the number of samples per channel is returned. 
**/
int
sadc_huffman_decode(const hufftree_t *hufftree, const unsigned int *input,
                    const unsigned int input_size, unsigned int *sample_num,
                    unsigned int **output)
{
  int                 size = 0, cur_size = input_size * 2 + 2;
  unsigned int       *val = NULL;
  unsigned int        i, pos = 0;
  hufftree_t         *ht_pointer = (hufftree_t *)hufftree;
  unsigned int        samplenum = 0;


  val = (unsigned int *)calloc(cur_size, sizeof(unsigned int));

  /* first 9 bit are always 'samplenum - 1' */
  samplenum = (input[0] & 511) + 1;
  if (sample_num == NULL)
    return -1;

  *sample_num = samplenum;
  /* next 6 bit are always channel nr */
  val[size] = (input[0] >> 9) & 63;
  size++;
  pos = 15;

  for (i = 0; i < input_size; i++) {
    while (1) {
      /* check if actual element is a leaf */
      if (!ht_pointer->next[0] && !ht_pointer->next[1]) {
        val[size] = ht_pointer->val;
        size++;
        /* realloc new memory if running out of it */
        if (size >= cur_size) {
          cur_size += 2;
          {
            unsigned int       *tmp;
            tmp =
                (unsigned int *)realloc(val, cur_size * sizeof(unsigned int));
            if (!tmp) {
              free(val);
              return -1;
            }
            val = tmp;
          }
        }

        if (size % (samplenum + 1) == 0) {
          /* check if this was the last channel */
          if (i == input_size - 1) {
            /* check whether unused bits are padded to zero (this is an
               important sanity check on the input) */
            if (input[i] >> pos != 0) {
              free_not_null(val);
              return -2;
            }
            break;
          }

          /* pos--; */
          if (pos + 6 > WORD_LEN) {

            val[size] = (input[i] >> pos) & ((1 << (WORD_LEN - pos)) - 1);
            i++;
            val[size] +=
                (input[i] & ((1 << (pos + 6 - WORD_LEN)) - 1)) << (WORD_LEN -
                                                                   pos);
            pos = pos + 6 - WORD_LEN;
            size++;
            if (size >= cur_size) {
              cur_size += 2;
              {
                unsigned int       *tmp;
                tmp =
                    (unsigned int *)realloc(val,
                                            cur_size * sizeof(unsigned int));
                if (!tmp) {
                  free(val);
                  return -1;
                }
                val = tmp;
              }
            }

          } else {
            val[size] = (input[i] >> pos) & 63;
            pos += 6;
            size++;
            if (size >= cur_size) {
              cur_size += 2;
              {
                unsigned int       *tmp;
                tmp =
                    (unsigned int *)realloc(val,
                                            cur_size * sizeof(unsigned int));
                if (!tmp) {
                  free(val);
                  return -1;
                }
                val = tmp;
              }
            }
          }

        }
        ht_pointer = (hufftree_t *)hufftree;
        if (pos > WORD_LEN - 1)
          break;
      }

      ht_pointer = ht_pointer->next[(input[i] >> pos) & 1];
      pos++;

      if (pos > WORD_LEN - 1)
        break;
    }
    pos = 0;
  }

  *output = val;
  return size;
}

/**
 * checks if tree and lookuptable are ok. Return value is <0 on error.
**/
int
sanity_check(const hufftree_t *hufftree, const hufflookup_t *hufflookup,
             const unsigned int lookup_size)
{
  int                 i;

  for (i = 0; i < lookup_size; i++) {
    int                 pos = 0;
    hufftree_t         *ht_pointer = (hufftree_t *)hufftree;

    while (1) {
      ht_pointer = ht_pointer->next[(hufflookup[i].huf >> pos) & 1];

      if (!ht_pointer)
        return -1;

      if (!ht_pointer->next[0] && !ht_pointer->next[1]) {
        if (ht_pointer->val != hufflookup[i].val)
          return -1;
        else
          break;
      }

      pos++;
      if (pos > hufflookup[i].length)
        return -1;
    }
  }

  return 0;
}

/** makes a code in binary format from a code in string format and reverts it 
*/
uint32
make_int_from_code(const char *code)
{
  unsigned int        intcode = 0;
  int                 i;

  for (i = strlen(code) - 1; i >= 0; i--) {
    intcode = intcode << 1;
    if (code[i] == '1')
      intcode += 1;
  }

  return intcode;
}

/** 
 * creates a lookup table as an array 'hufflookup' for all values in the tree 
 * stored in c-string 'string' and returns its size. The values in the array are 
 * sorted by quantity beginning from "0". If the return value is < 0 an
 * error occured 
**/
int
create_lookup_from_string(const char *string, hufflookup_t **hufflookup)
{
  int                 cur_size = 4096, val = 0;
  unsigned int        str_pos = 0;
  char                code[1024];
  hufflookup_t       *lookup = NULL;

  lookup = (hufflookup_t *)calloc(cur_size, sizeof(hufflookup_t));
  memset(lookup, 0, cur_size * sizeof(hufflookup_t));

  for (;;) {
    int                 ret;

    if (string[str_pos] == '\n') {
      str_pos++;
      continue;
    }

    ret = sscanf(&string[str_pos], "%i\t%s\n", &val, code);
    if (ret == EOF)
      break;

    if (ret != 2) {
      free_not_null(lookup);
      return -2;
    }

    if (val >= cur_size) {
      cur_size = val + 1;
      {
        hufflookup_t       *tmp;
        tmp =
            (hufflookup_t *)realloc(lookup, cur_size * sizeof(hufflookup_t));
        if (!tmp) {
          free(lookup);
          return -2;
        }
        lookup = tmp;
      }
    }

    /* check if element already exists */
    if (lookup[val].val != 0 || lookup[val].huf != 0
        || lookup[val].length != 0) {
      free_not_null(lookup);
      return -1;
    }
    lookup[val].val = val;
    lookup[val].huf = make_int_from_code(code);
    lookup[val].length = strlen(code);

    for (; str_pos < strlen(string); str_pos++) {
      if (string[str_pos] == '\n') {
        str_pos++;
        break;
      }
    }
    if (str_pos >= strlen(string) - 1)
      break;

  }

  *hufflookup = lookup;
  return cur_size;
}

/** 
 * creates a lookup table as an array 'hufflookup' for all values in the tree 
 * stored in file 'path' and returns its size. The values in the array are 
 * sorted by quantity beginning from "0". If the return value is < 0 an
 * error occured 
**/
int
create_lookup_from_file(const char *path, hufflookup_t **hufflookup)
{
  int                 cur_size = 1024, val = 0;
  char                code[1024];
  FILE               *Ftree = NULL;
  hufflookup_t       *lookup = NULL;

  Ftree = fopen(path, "r");
  if (!Ftree)
    return -1;

  lookup = (hufflookup_t *)calloc(cur_size, sizeof(hufflookup_t));

  for (;;) {
    int                 ret;

    ret = fscanf(Ftree, "%i\t%s\n", &val, code);
    if (ret == EOF)
      break;
    if (ret != 2) {
      free_not_null(lookup);
      return -2;
    }

    if (val >= cur_size) {
      cur_size = val + 1;
      {
        hufflookup_t       *tmp;
        tmp =
            (hufflookup_t *)realloc(lookup, cur_size * sizeof(hufflookup_t));
        if (!tmp) {
          free(lookup);
          return -2;
        }
        lookup = tmp;
      }
    }

    lookup[val].val = val;
    lookup[val].huf = make_int_from_code(code);
    lookup[val].length = strlen(code);

  }


  fclose(Ftree);
  *hufflookup = lookup;
  return cur_size;
}

/** 
 * returns the code-length of a value 'val' in the tree.
 * if 'val' is not found the return value is -1.
 * The code of the value 'val' is stored as a nullterminated string in 'code'.
**/
int
val_code(const int val, unsigned int level, const hufftree_t *hufftree,
         char *code)
{
  int                 the_level;

  if (hufftree->next[0] || hufftree->next[1]) {
    if (hufftree->next[0]) {
      code[level] = '0';
      code[level + 1] = '\0';
      the_level = val_code(val, level + 1, hufftree->next[0], code);
      if (the_level > 0)
        return the_level;
    }

    if (hufftree->next[1]) {
      code[level] = '1';
      code[level + 1] = '\0';
      the_level = val_code(val, level + 1, hufftree->next[1], code);
      if (the_level > 0)
        return the_level;
    }

  } else {
    if (hufftree->val == val)
      return level;
  }

  return -1;
}

/**
 * looks if 'orig_hufftree' is a leaf. If it is a leaf the values are copied
 * to 'new_hufftree' otherwise the function calls itself with the next tree
 * member as argument. The return value is "0" on success, otherwise the return
 * value differs from "0".
**/
int
copy_tree_members(const hufftree_t *orig_hufftree, hufftree_t **new_hufftree)
{
  *new_hufftree = (hufftree_t *)calloc(1, sizeof(hufftree_t));
  if (!*new_hufftree)
    return -1;

  if ((!(orig_hufftree->next[0])) && (!(orig_hufftree->next[1]))) {
    (*new_hufftree)->val = orig_hufftree->val;
    (*new_hufftree)->freq = orig_hufftree->freq;
    return 0;
  }

  if (orig_hufftree->next[0]) {
    int                 ret;

    ret =
        copy_tree_members(orig_hufftree->next[0],
                          &((*new_hufftree)->next[0]));

    if (ret != 0) {
      free_not_null(*new_hufftree);
      return ret;
    }
  }

  if (orig_hufftree->next[1]) {
    int                 ret;

    ret =
        copy_tree_members(orig_hufftree->next[1],
                          &((*new_hufftree)->next[1]));

    if (ret != 0) {
      free_not_null(*new_hufftree);
      return ret;
    }
  }

  return 0;
}

/**
 * this function duplicates a given 'hufftree' and returns it.
 * In case of an error NULL is the return value
**/
hufftree_t         *
copy_tree(const hufftree_t *hufftree)
{
  hufftree_t         *newTree = NULL;

  if (copy_tree_members(hufftree, &newTree) != 0)
    return NULL;

  return newTree;
}

/** 
 * reads a huffman tree from a c-string.
 * If return value is NULL an error occured 
**/
hufftree_t         *
hufftree_read_from_string(const char *string)
{
  hufftree_t         *hufftree = NULL;
  unsigned int        str_pos = 0;

  hufftree = (hufftree_t *)calloc(1, sizeof(hufftree_t));

  memset(hufftree, 0, sizeof(hufftree_t));

  while (1) {
    int                 i, ret, val;
    char                code[1024];
    hufftree_t         *ht_pointer = NULL;

    memset(code, 0, sizeof(code));

    ret = sscanf(&string[str_pos], "%i\t%s", &val, code);
    if (ret == EOF)
      break;

    if (ret != 2) {
      return NULL;
    }

    ht_pointer = hufftree;

    for (i = 0; (unsigned int)i < sizeof(code); i++) {
      int                 tmp;

      if (code[i] == '\0') {
        ht_pointer->val = (unsigned int)val;
        break;
      }

      if (code[i] == '0')
        tmp = 0;
      else
        tmp = 1;

      if ((tmp != 0) && (tmp != 1)) {
        free_tree(hufftree);
        return NULL;
      }

      if (!ht_pointer->next[tmp]) {
        ht_pointer->next[tmp] = (hufftree_t *)calloc(1, sizeof(hufftree_t));
        memset(ht_pointer->next[tmp], 0, sizeof(hufftree_t));
        ht_pointer = ht_pointer->next[tmp];
        continue;
      } else {
        ht_pointer = ht_pointer->next[tmp];
        continue;
      }

    }

    for (; str_pos < strlen(string); str_pos++) {
      if (string[str_pos] == '\n') {
        str_pos++;
        break;
      }
    }
    if (str_pos >= strlen(string) - 1)
      break;

  }

  return hufftree;
}

/** 
 * reads a huffman tree from file 'path'. 
 * If return value is NULL an error occured 
**/

hufftree_t         *
hufftree_read(const char *path)
{
  hufftree_t         *hufftree = NULL;
  FILE               *Ftree = NULL;

  Ftree = fopen(path, "r");
  if (!Ftree)
    return NULL;

  hufftree = (hufftree_t *)calloc(1, sizeof(hufftree_t));
  memset(hufftree, 0, sizeof(hufftree_t));

  while (1) {
    int                 i, ret, val;
    char                code[1024];
    hufftree_t         *ht_pointer = NULL;

    memset(code, 0, sizeof(code));

    ret = fscanf(Ftree, "%i\t%s\n", &val, code);

    if (ret == EOF)
      break;
    if (ret != 2) {
      return NULL;
    }

    ht_pointer = hufftree;

    for (i = 0; i < (int)sizeof(code); i++) {
      int                 tmp;

      if (code[i] == '\0') {
        ht_pointer->val = (unsigned int)val;
        break;
      }

      if (code[i] == '0')
        tmp = 0;
      else
        tmp = 1;

      if ((tmp != 0) && (tmp != 1)) {
        free_tree(hufftree);
        return NULL;
      }

      if (!ht_pointer->next[tmp]) {
        ht_pointer->next[tmp] = (hufftree_t *)calloc(1, sizeof(hufftree_t));
        memset(ht_pointer->next[tmp], 0, sizeof(hufftree_t));
        ht_pointer = ht_pointer->next[tmp];
        continue;
      } else {
        ht_pointer = ht_pointer->next[tmp];
        continue;
      }

    }

  }

  fclose(Ftree);

  return hufftree;
}

/** 
 * writes a huffman tree 'hufftree' to a file 'path'. 'huffstat' is 
 * containing the huffman statistics the tree was created from.
 * if file can not be opened the return value is -1 otherwise its 0 
**/
int
hufftree_dump(const char *path, const hufftree_t *hufftree,
              const huffstat_t *huffstat)
{
  FILE               *Ftree = NULL;
  int                 i;
  char                code[1024], buffer[1024];

  code[0] = '\0';

  Ftree = fopen(path, "w");
  if (!Ftree) {
    return -1;
  }

  for (i = huffstat->size - 1; i >= 0; i--) {
    if (val_code(huffstat->elt[i].val, 0, hufftree, code) > 0) {
      snprintf(buffer, sizeof(buffer), "%i\t%s\n", huffstat->elt[i].val,
               code);
      fwrite(buffer, sizeof(char), strlen(buffer), Ftree);
    } else
      continue;
  }

  fclose(Ftree);

  return 0;
}

/** 
 * writes huffman statistics 'huffstat' to a file 'path'. if file can not be 
 * opened the return value is -1 otherwise its 0 
**/
int
huffstat_dump(const char *path, huffstat_t *huffstat)
{
  FILE               *Fstat = NULL;

  Fstat = fopen(path, "w");
  if (!Fstat)
    return -1;

  /* first write the array size into file */
  fwrite(&huffstat->size, sizeof(unsigned int), 1, Fstat);

  /* now write all the other data */
  fwrite(huffstat->elt, sizeof(struct elt_s), huffstat->size, Fstat);
  /* 
     fwrite(huffstat->freq, sizeof(*(huffstat->freq)), huffstat->size,
     Fstat); */

  fclose(Fstat);

  return 0;
}

/** 
 * reads huffman statistics from a file 'path'. If return value is NULL an 
 * error occured 
**/
huffstat_t         *
huffstat_read(const char *path)
{
  huffstat_t         *huffstat = NULL;
  FILE               *Fstat = NULL;

  huffstat = (huffstat_t *)calloc(1, sizeof(huffstat_t));

  Fstat = fopen(path, "r");
  if (!Fstat)
    return NULL;

  /* read the size first */
  fread(&huffstat->size, sizeof(unsigned int), 1, Fstat);

  /* allocate memory */

  huffstat->elt =
      (struct elt_s *)calloc(huffstat->size, sizeof(struct elt_s));
  /* 
     huffstat->val = (int *)calloc(huffstat->size, sizeof(int));
     huffstat->freq = (uint64 *)calloc(huffstat->size,
     sizeof(*(huffstat->freq))); */
  /* read array contents */
  fread(huffstat->elt, sizeof(struct elt_s), huffstat->size, Fstat);
  /* 
     fread(huffstat->val, sizeof(int), huffstat->size, Fstat);
     fread(huffstat->freq, sizeof(*(huffstat->freq)), huffstat->size, Fstat); */

  fclose(Fstat);

  return huffstat;
}

/** 
 * creates new blank huffman statistics with all values from 0 
 * to 'range_end'. All frequency information is set to 1 
**/
void
make_allstat(huffstat_t *huffstat, const int range_begin, const int range_end)
{
  int                 val;
  unsigned int        i;

  if (range_end < range_begin)
    return;

  /* range_begin is ignored for better performance */
  val = 0;
  huffstat->size = range_end + 1;

  huffstat->elt = (elt_t *)calloc(huffstat->size, sizeof(elt_t));
  /* 
     huffstat->val = (int *)calloc(huffstat->size, sizeof(int));
     huffstat->freq = (uint64 *)calloc(huffstat->size,
     sizeof(*(huffstat->freq))); */

  for (i = 0; (i < huffstat->size) && (val <= range_end); i++, val++) {
    huffstat->elt[i].val = val;
    huffstat->elt[i].freq = 1;
  }

}

/** 
 * inserts a new value 'value' into huffman statistics structure 'huffstat'. 
 * If 'value' is already present the freq entry is just increased. In this
 * version all values MUST be present!(use 'make_allstat()'!!!!)
**/
void
insert_value(huffstat_t *huffstat, int value)
{
  int                 val;

  /* wraparound for 12 bit data types */
  val = value & 4095;

  huffstat->elt[val].freq++;
  return;
}

/** 
 * sorts a huffman statistics structure beginning from the lowest frequented
 * value. this algorithm is just for testing and should be replaced by a
 * better one (qsort) 
**/
#if 0
void
sort_huffstat(huffstat_t *huffstat)
{
  unsigned int        i;

  for (i = 0; i < huffstat->size; i++) {
    unsigned int        min_pos = i, c;
    uint64              tmp_v, tmp_f;

    /* find minimum */
    for (c = i; c < huffstat->size; c++) {
      if (huffstat->freq[c] < huffstat->freq[min_pos])
        min_pos = c;
    }

    /* now switch the minimum with the i-th Element in the arrays */
    if (min_pos == i)
      continue;
    tmp_f = huffstat->freq[i];
    huffstat->freq[i] = huffstat->freq[min_pos];
    huffstat->freq[min_pos] = tmp_f;

    tmp_v = huffstat->val[i];
    huffstat->val[i] = huffstat->val[min_pos];
    huffstat->val[min_pos] = tmp_v;

  }
}

#endif

/** 
 * sorts a huffman statistics structure beginning from the lowest frequented
 * value. this algorithm works like quicksort 
**/
#if 0
void
quick_sort_huffstat(huffstat_t *huffstat, unsigned int al, unsigned int ar)
{
  unsigned int        left = al, right = ar;
  uint64              pivo = huffstat->freq[(al + ar) / 2];
  uint64              utmp;
  int                 itmp;


  do {
    while (huffstat->freq[left] < pivo)
      left++;
    while (huffstat->freq[right] > pivo)
      right--;

    if (left <= right) {
      utmp = huffstat->freq[left];
      huffstat->freq[left] = huffstat->freq[right];
      huffstat->freq[right] = utmp;

      itmp = huffstat->val[left];
      huffstat->val[left] = huffstat->val[right];
      huffstat->val[right] = itmp;

      left++;
      right--;
    }


  } while (left < right);

  if (al < right)
    quick_sort_huffstat(huffstat, al, right);
  if (left < ar)
    quick_sort_huffstat(huffstat, left, ar);
}
#endif


/** 
 * creates and returns a huffman tree from a huffman statistics structure 
 * 'huffstat'. 'huffstat' MUST be sorted with the value of lowest frequency 
 * at the beginning! If it's impossible to create a tree NULL is returned 
**/
hufftree_t         *
make_hufftree(huffstat_t *huffstat)
{
  hufftree_t        **the_hufftree = NULL;
  hufftree_t        **tmp_tree = NULL;
  hufftree_t         *ret;

  /* size of huffstat too low */
  if (huffstat->size < 2) {
    return NULL;
  }

  /* allocate enough tree members for all elements */
  {
    unsigned int        i;

    tmp_tree = (hufftree_t **)calloc(huffstat->size, sizeof(hufftree_t *));

    for (i = 0; i < huffstat->size; i++) {
      tmp_tree[i] = (hufftree_t *)calloc(1, sizeof(hufftree_t));
    }
  }

  /* fill all values from the arrays into the tree members */
  {
    unsigned int        i;

    for (i = 0; i < huffstat->size; i++) {
      (*tmp_tree[i]).val = (unsigned int)huffstat->elt[i].val;
      (*tmp_tree[i]).freq = huffstat->elt[i].freq;

      (*tmp_tree[i]).next[0] = NULL;
      (*tmp_tree[i]).next[1] = NULL;
    }
  }

  /* now allocate all additional tree members */
  {
    unsigned int        i;

    the_hufftree = (hufftree_t **)calloc(huffstat->size - 1,
                                         sizeof(hufftree_t *));

    for (i = 0; i < huffstat->size - 1; i++) {
      the_hufftree[i] = (hufftree_t *)calloc(1, sizeof(hufftree_t));
    }
  }
  /* initialize all new members */
  {
    unsigned int        i;

    for (i = 0; i < huffstat->size - 1; i++) {
      (*the_hufftree[i]).val = 0;
      (*the_hufftree[i]).freq = 0;

      (*the_hufftree[i]).next[0] = NULL;
      (*the_hufftree[i]).next[1] = NULL;
    }
  }

  /* make the tree */
  {
    int                 a;
    unsigned int        c = 0, e = huffstat->size - 2;
    hufftree_t         *min_1 = NULL, *min_2 = NULL;


    for (a = huffstat->size - 2; a >= 0; a--) {
      unsigned int        d;

      if (c < huffstat->size) {
        min_1 = tmp_tree[c];

        if (c < huffstat->size - 1)
          min_2 = tmp_tree[c + 1];
        else {
          min_2 = the_hufftree[e];
          e--;
        }
      } else {
        if (e > 1) {
          min_1 = the_hufftree[e];
          min_2 = the_hufftree[e - 1];
          e -= 2;
        }

      }

      /* get the 2 elemtents with the lowest frequency */
      /* 
         for (d = c + 1; d < huffstat->size; d++) {

         if (tmp_tree[d]->freq < min_1->freq) { min_2 = min_1; min_1 =
         tmp_tree[d]; continue; }

         if (tmp_tree[d]->freq < min_2->freq) { min_2 = tmp_tree[d];
         continue; } } */

      for (d = e; (d > (unsigned int)a) && (d <= huffstat->size - 2)
           && (c < huffstat->size); d--) {
        if (the_hufftree[d]->freq < min_1->freq) {
          if (min_2->freq >= min_1->freq) {
            min_2 = min_1;
            min_1 = the_hufftree[d];
            c--;
          } else {
            min_1 = min_2;
            min_2 = the_hufftree[d];
            c -= 2;
          }
          e = d - 1;
          continue;
        }

        if ((the_hufftree[d]->freq < min_2->freq)
            && (the_hufftree[d] != min_1)) {
          min_2 = the_hufftree[d];
          e = d - 1;
          c--;
          continue;
        }
      }


      (*the_hufftree[a]).next[0] = min_1;
      (*the_hufftree[a]).next[1] = min_2;

      (*the_hufftree[a]).freq = min_1->freq + min_2->freq;

      c += 2;

    }

  }

  ret = the_hufftree[0];


  free(the_hufftree);
  free(tmp_tree);


  return ret;
}

/** 
 * calculates the bitlength per value of a given tree with frequency 
 * information. This function should not be used anymore!
 * Use val_code() and the huffstat structure instead, this will be 
 * more exactly 
**/
void
calc_bit_per_val(hufftree_t *hufftree, const unsigned int level, float *bits,
                 uint64 *freq)
{

  uint64              l_freq = 0;
  float               l_bit = 0;
  uint64              r_freq = 0;
  float               r_bit = 0;

  if (hufftree->next[0]) {
    calc_bit_per_val(hufftree->next[0], level + 1, &l_bit, &l_freq);
  } else {
    l_bit = level;
    l_freq = hufftree->freq;
  }

  if (hufftree->next[1]) {
    calc_bit_per_val(hufftree->next[1], level + 1, &r_bit, &r_freq);
  } else {
    r_bit = level;
    r_freq = hufftree->freq;
  }

  *freq = l_freq + r_freq;
  *bits = (l_bit * (float)l_freq + r_bit * (float)r_freq) / (float)(*freq);

}

/** 
 * destroys a tree and frees all memory 
**/
void
free_tree(hufftree_t *hufftree)
{
  if ((!hufftree->next[0]) && (!hufftree->next[1]))
  {
    free_not_null(hufftree);
    return;
  }

  if ((!hufftree->next[0]->next[0]) && (!hufftree->next[0]->next[1]))
    free_not_null(hufftree->next[0]);
  else
    free_tree(hufftree->next[0]);

  if ((!hufftree->next[1]->next[0]) && (!hufftree->next[1]->next[1]))
    free_not_null(hufftree->next[1]);
  else
    free_tree(hufftree->next[1]);
  
  if ((!hufftree->next[0]) && (!hufftree->next[1]))
  {
    free_not_null(hufftree);
  }
  
}

/**
 * compares 2 huffstat frequency entries and returns an integer less than, equal to
 * or greater than zero if the 'freq' of first argument is greater than, equal or
 * less than the second.
**/
int
compar_huffstat(const void *e1, const void *e2)
{
  elt_t              *elt1 = (elt_t *)e1, *elt2 = (elt_t *)e2;

  if (elt1->freq > elt2->freq)
    return 1;
  if (elt1->freq < elt2->freq)
    return -1;

  return 0;
}

} //namespace CS
