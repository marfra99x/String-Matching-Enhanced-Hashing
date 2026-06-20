#include "include/define.h"
#include "include/main.h" 
#define Q 6
#define cc uint16_t
#define HASH(j) (*((cc*)(y+j-1))<<4) + (*((cc*)(y+j-3))<<2) + (*((cc*)(y+j-5)))

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
        h = (*((cc*)(x+i-1))<<4)+ (*((cc*)(x+i-3))<<2) + (*((cc*)(x+i-5)));
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
