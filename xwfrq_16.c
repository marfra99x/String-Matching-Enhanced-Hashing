#include <stdint.h>
#include <smmintrin.h>
#include "include/define.h"
#include "include/main.h" // defines the search interface for time and statistics.
#define Q 16
#define mi __m128i

/*
 * Pre-process q-gram factors of the pattern.
 */
int preprocessingQ(unsigned char *x, int m, char *F) {
    int i,j;
    unsigned short h;
    mi y_ptr;
    for(i=0; i<256*256; i++) F[i] = FALSE;
    for(i=Q-1; i<m; i++) {
        y_ptr = _mm_loadu_si128(( mi *)(x+i-15));
        h=_mm_movemask_epi8(y_ptr<<6);
        F[h] = TRUE;
    }
}


int search(unsigned char *x, int m, unsigned char *y, int n) {
    int i, j, p, k, count, test;
    char F[256*256];
    unsigned short h;
    if(m<Q) return -1;

    BEGIN_PREPROCESSING
    /* Preprocessing */
    int plen = m;
    if(m%Q!=0) m = m-(m%Q);
    int mq = m-Q+1;
    preprocessingQ(x,m,F);
    for(i=0; i<m; i++) y[n+i]=x[i];
    END_PREPROCESSING

    BEGIN_SEARCHING
    /* Searching */
    mi y_ptr;
    count = 0;
    if( !memcmp(x,y,plen) ) count++;
    j=m;
    while (j < n) {
        y_ptr = _mm_loadu_si128(( mi *)(y+j-15));
        h=_mm_movemask_epi8(y_ptr<<6);
        i = j-m+1;
        while((test=F[h]) && j>i+Q-1) {
            j-=Q;
            y_ptr = _mm_loadu_si128(( mi *)(y+j-15));
            h=_mm_movemask_epi8(y_ptr<<6);
        }
        if(j==i+Q-1 && test) {
            if(memcmp(y+i,x,plen) == 0 && i<=n-plen)
                count++;
        }
        j+=mq;
    }
    END_SEARCHING
    return count;
}
