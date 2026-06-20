#include "include/define.h"
#include "include/main.h"
#include <stdint.h>

#define CHARTYPE unsigned char
#define MAXPAT 51000
#define LIM 8192
#define W 64
#define Q 13

typedef uint64_t bv;

int search(unsigned char* pat, int m, unsigned char *t, int n) {

bv pre, B[LIM];
int patlen, a, b;
uint16_t mask=8191;

   if(m<Q) return -1;

   BEGIN_PREPROCESSING

    int i, j, r, x, y;
    uint16_t ch=0;
    bv d;
    int k, count = 0;
 
    pre = *((uint64_t*) (pat));

    for (i=0; i<LIM; i++) B[i]=0;

    x=0; if (m<=4096) x=1;
    a=(m-Q+1-x)/W+x; 
    r=(m-Q+1)/a; 
    if (r>W) r=W;
    y=W-r;
    b=r*a;

    i=m-1; 
    for (i=m-1; i>m-Q; i--)
       ch = ((ch << 1) + (unsigned)pat[i])&mask;
    i=m-Q;
    while (i>m-b-Q) { 
       for (j=0; j<a; j++) {
          ch = ((ch << 1) + (unsigned)pat[i])&mask;
          B[ch] |= (bv)1 << y;
          i--;
       }
       y++;
    }


   END_PREPROCESSING

   BEGIN_SEARCHING

#define F(z) ( (((( (((( (((( (((( (((( (((( \
(unsigned)*(z))<< 1) + \
(unsigned)*((z)-1)) << 1) + (unsigned)*((z)-2)) << 1) + \
(unsigned)*((z)-3)) << 1) + (unsigned)*((z)-4)) << 1) + \
(unsigned)*((z)-5)) << 1) + (unsigned)*((z)-6)) << 1) + \
(unsigned)*((z)-7)) << 1) + (unsigned)*((z)-8)) << 1) + \
(unsigned)*((z)-9)) << 1) + (unsigned)*((z)-10))<< 1) + \
(unsigned)*((z)-11))<< 1) + (unsigned)*((z)-12))&mask

    i = m-1;
    t[n] = 0;

   while (i<n) {
      d = B[F(t+i)];
      if (d==0) i+=b;
      else {
         j = i;
         do {i-=a; 
            d = (d<<1)&B[F(t+i)];}
         while (d);
         i += b; 
         if (i==j) {
            for (k=0; k<a; k++) 
               if(pre == *((uint64_t*) (t+i-m+1+k)))
                  if(memcmp(t+i-m+1+k,pat,m) == 0) 
                     count++;
            i+=a;
         }
      }
    }

  END_SEARCHING
  return count;
}

