/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 * This is an implementation of the EPSM algorithm in S. Faro and O. M. Kulekci.
 * It includes corrections to the original implementations gently provided by Jorma Tarhio and Jan Holub
 */

#include "include/define.h"
#include "include/main.h"
#include <memory.h>
#include <smmintrin.h>
#include <inttypes.h>

#define HASHSIZE 16

typedef union
{
    __m128i *symbol16;
    unsigned char *symbol;
} TEXT;

typedef union
{
    __m128i v;
    unsigned int ui[4];
    unsigned short int us[8];
    unsigned char uc[16];
} VectorUnion;

typedef struct list
{
    struct list *next;
    int pos;
} LIST;

int search16(unsigned char *pattern, int patlen, unsigned char *x, int textlen)
{
    LIST *flist[2048]; // 11-bit hash gives the best result according to our tests
    LIST *t;

    unsigned int i, filter, shift = (patlen / 4) - 1;
    if (shift < 1) shift = 1; // Ensure shift is at least 1
    unsigned int crc, seed = 123456789, mask;
    unsigned int *ptr32;
    unsigned int *lastchunk;
    unsigned char *charPtr;
    int count = 0, tmppatlen = (patlen / 4) * 4;
    mask = 2047; // 11-bit mask

    BEGIN_PREPROCESSING
    memset(flist, 0, sizeof(LIST *) * 2048);

    for (i = 1; i < tmppatlen - 3; i++) // Process 4-byte chunks
    {
        ptr32 = (unsigned int *)(&pattern[i]);
        crc = _mm_crc32_u32(seed, *ptr32);
        filter = (unsigned int)(crc & mask);

        if (flist[filter] == 0)
        {
            flist[filter] = (LIST *)malloc(sizeof(LIST));
            if (flist[filter])
            {
                flist[filter]->next = 0;
                flist[filter]->pos = i;
            }
        }
        else
        {
            t = flist[filter];
            while (t->next != 0)
                t = t->next;
            t->next = (LIST *)malloc(sizeof(LIST));
            if (t->next)
            {
                t = t->next;
                t->next = 0;
                t->pos = i;
            }
        }
    }
    END_PREPROCESSING

    BEGIN_SEARCHING
    lastchunk = (unsigned int *)&x[((textlen - tmppatlen) / 4) * 4 + 1];
    ptr32 = (unsigned int *)&x[(shift - 1) * 4];

    crc = _mm_crc32_u32(seed, *ptr32);
    filter = (unsigned int)(crc & mask);
    if (flist[filter])
    {
        charPtr = (unsigned char *)ptr32;
        t = flist[filter];
        while (t)
        {
            if (t->pos <= 4 * (shift - 1))
            {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
            }
            t = t->next;
        }
    }
    ptr32 += shift;

    crc = _mm_crc32_u32(seed, *ptr32);
    filter = (unsigned int)(crc & mask);
    if (flist[filter])
    {
        charPtr = (unsigned char *)ptr32;
        t = flist[filter];
        while (t)
        {
            if (t->pos <= 4 * (2 * shift - 1))
            {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
            }
            t = t->next;
        }
    }
    ptr32 += shift;

    while (ptr32 < lastchunk)
    {
        crc = _mm_crc32_u32(seed, *ptr32);
        filter = (unsigned int)(crc & mask);

        if (flist[filter])
        {
            charPtr = (unsigned char *)ptr32;
            t = flist[filter];
            while (t)
            {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
                t = t->next;
            }
        }
        ptr32 += shift;
    }

    ptr32 -= shift;
    charPtr = (unsigned char *)ptr32;
    charPtr += tmppatlen - 1; // First unchecked position where P may end

    while (charPtr < &x[textlen - 1])
    {
        if (0 == memcmp(pattern, charPtr - tmppatlen + 1, patlen))
            count++;
        charPtr++;
    }
    END_SEARCHING

    return count;
}

int search(unsigned char *pattern, int patlen, unsigned char *x, int textlen)
{
    if (patlen <= 7)
        return 0;
    else
        return search16(pattern, patlen, x, textlen);
}
