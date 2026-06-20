#include "include/define.h"
#include "include/main.h" // defines the search interface for time and statistics.
//#include "include/GRAPH.h"
#define Q 4
#define cc uint16_t
#define HASH(j) (*((cc*)(y+j))<<2) + (*((cc*)(y+j-2)))

/*
 * Pre-process q-gram factors of the pattern.
 */
int preprocessingQ(unsigned char *x, int m, char *F) {
    int i,j;
    unsigned short h;
    //int fact = m < 8 ? m : 8;
    for(i=0; i<256*256; i++) F[i] = FALSE;
    for(i=Q-2; i<m-1; i++) {
        h = (*((cc*)(x+i))<<2)+ (*((cc*)(x+i-2)));
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
    //for(i=0; i<m; i++) y[n+i]=x[i];
    END_PREPROCESSING

    BEGIN_SEARCHING
    /* Searching */
    count = 0;
    if( !memcmp(x,y,plen) ) count++;
    j=m-1;
    while (j < n) {
        h = HASH(j);
        i = j-m+2;
        while((test=F[h]) && j>i+Q-2) {
            j-=Q;
            h = HASH(j);
        }
        if(j==i+Q-2 && test) {
            if(memcmp(y+i,x,plen) == 0 && i<=n-plen)
                count++;
        }
        j+= m-Q+1;
    }
    END_SEARCHING
    return count;
}
