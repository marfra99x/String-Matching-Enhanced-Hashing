#include "include/define.h"
#include "include/main.h" // defines the search interface for time and statistics.
#include <smmintrin.h>
#include <inttypes.h>

#define Q 8
#define c8 uint64_t
#define HASH(XX,II) _mm_crc32_u64(seed,*((uint64_t*)(XX+II-7)))

/*
 * Pre-process q-gram factors of the pattern.
 */
int preprocessingQ(unsigned char *x, int m, char *F) {
    int i,j;
    unsigned short h;
    uint64_t seed = 123456789;
    for(i=0; i<256*256; i++) F[i] = FALSE;
    for(i=Q-1; i<m; i++) {
        //int stop = (i-fact+1)>=Q-1?(i-fact+1):Q-1;
        h = HASH(x,i);
        F[h] = TRUE;
    }
}


int search(unsigned char *x, int m, unsigned char *y, int n) {
    int i, j, p, k, count, test;
    char F[256*256];
    unsigned short h;
    uint64_t seed = 123456789;

    if(m<Q) return -1;

    BEGIN_PREPROCESSING
    /* Preprocessing */
    int plen = m;
    if(m%Q!=0) m = m-(m%Q);
    preprocessingQ(x,m,F);
    for(i=0; i<m; i++) y[n+i]=x[i];
    END_PREPROCESSING

    BEGIN_SEARCHING
    /* Searching */
    count = 0;
    if( !memcmp(x,y,plen) ) count++;
    j=m;
    while (j < n) {
        h = HASH(y,j);
        i = j-m+1;
        while((test=F[h]) && j>i+Q-1) {
            j-=Q;
            h = HASH(y,j);
        }
        if(j==i+Q-1 && test) {
            if(memcmp(y+i,x,plen) == 0 && i<=n-plen)
                count++;
        }
        j+= m-Q+1;
    }
    END_SEARCHING
    return count;
}
