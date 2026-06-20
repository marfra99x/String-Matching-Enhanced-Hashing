#include "include/define.h"
#include "include/main.h" // defines the search interface for time and statistics.
//#include "include/GRAPH.h"
#define Q 8
#define cc uint16_t
#define HASH(j) (*((cc*)(y+j-1))<<6) + (*((cc*)(y+j-3))<<4) + (*((cc*)(y+j-5))<<2) + (*((cc*)(y+j-7)))

/*
 * Pre-process q-gram factors of the pattern.
 */
int preprocessingQ(unsigned char *x, int m, char *F) {
    int i,j;
    unsigned short h;
    //int fact = m < 8 ? m : 8;
    for(i=0; i<256*256; i++) F[i] = FALSE;
    for(i=Q-1; i<m; i++) {
        //int stop = (i-fact+1)>=Q-1?(i-fact+1):Q-1;
        h = (*((cc*)(x+i-1))<<6)+ (*((cc*)(x+i-3))<<4)
          + (*((cc*)(x+i-5))<<2) + (*((cc*)(x+i-7)));
        F[h] = TRUE;
        /* DEAD CODE 
        for(j=i-Q; j>=stop; j-=Q) {
            h = (h<<16) + (x[j]<<14) + (x[j-1]<<12) + (x[j-2]<<10) 
                 + (x[j-3]<<8) + (x[j-4]<<6) + (x[j-5]<<4) 
                 + (x[j-6]<<2) + x[j-7];
            F[h] = TRUE;
        }
        */
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
    //int mq = m-Q+1;
    preprocessingQ(x,m,F);
    for(i=0; i<m; i++) y[n+i]=x[i];
    END_PREPROCESSING

    BEGIN_SEARCHING
    /* Searching */
    count = 0;
    if( !memcmp(x,y,plen) ) count++;
    j=m;
    while (j < n) {
        h = HASH(j);
        i = j-m+1;
        while((test=F[h]) && j>i+Q-1) {
            j-=Q;
            h = HASH(j);
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
