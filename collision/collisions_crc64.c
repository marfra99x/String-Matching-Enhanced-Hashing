#include "include/define.h"
#include "include/main.h" 
#include <smmintrin.h>
#include <inttypes.h>

#define TYPE 64
#define CON 16
#define MASK ((1<<CON)-1)

#define Q 8
#define HASH(XX,II) _mm_crc32_u64(123456789,*((uint64_t*)(XX+II-Q+1)))&MASK
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
        h = HASH(y,j);
        if (f[h]>-1 && READ(y,j) != READ(y,f[h])) count++;
        if (f[h]<0) f[h]=j;
        j++;
    }
    END_SEARCHING
    if(count==0)return 1;
    return count;
}

