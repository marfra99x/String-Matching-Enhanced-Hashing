#include "include/define.h"
#include "include/main.h" 
#include <smmintrin.h>
#include <inttypes.h>

#define TYPE 64
#define CON 16
#define MASK ((1<<CON)-1)

#define Q 8
#define cc uint16_t
#define HASH(j) (y[j]<<10) + (y[j-1]<<8) + (y[j-2]<<6) + (y[j-3]<<4) + (y[j-4]<<2) + y[j-5]
#define READ(XX,II) *((uint64_t*)(XX+II-Q+1))


#define LIM 256*256

int search(unsigned char *x, int m, unsigned char *y, int n) {
    int j, count=0;
    int f[LIM];
    uint16_t h;

    BEGIN_PREPROCESSING
    for (j=0;j<LIM;j++) f[j]=-1;
    END_PREPROCESSING

    BEGIN_SEARCHING
    j=Q-1;
    while (j < n) {
        h = HASH(j);
        if (f[h]>-1 && READ(y,j) != READ(y,f[h])) count++;
        if (f[h]<0) f[h]=j;
        j++;
    }
    END_SEARCHING
    return count;
}

