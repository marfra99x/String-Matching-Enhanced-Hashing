#include <smmintrin.h>
#include "include/define.h"
#include "include/main.h"

#define LIM 65536
#define W 64
#define Q 8
#define QM 7 //=Q-1
#define H(XX,II) _mm_crc32_u64(54319,*((uint64_t*)(XX+II-QM)))

typedef uint64_t bv;

int search(unsigned char* pat, int m, unsigned char *t, int n) {

   bv d, pre, B[LIM];
   int i, j, k, r, a, b, x, y, count = 0;
   uint16_t h;

   if(m<Q) return -1;
   
   BEGIN_PREPROCESSING

   pre = *((uint64_t*) (pat));

   for (i=0; i<LIM; i++) B[i]=0;

   x=0; if (m<=4096) x=1;
   a=(m-Q+1-x)/W+x; 
   r=(m-Q+1)/a; 
   if (r>W) r=W;
   b=r*a;

   y=W-r;
   i=m-Q;
   while (i>m-b-Q) {
      for (j=1; j<=a; j++) {
         B[h=H(pat,i+QM)] |= (bv)1 << y;
         i--;  
      }
      y++;
   }
   END_PREPROCESSING

   BEGIN_SEARCHING

   i = m-1;
   t[n] = 0;

   while (i<n) {
      d = B[h=H(t,i)];
      if (d==0) i+=b;
      else {
         j = i;
         do {i-=a; 
            d = (d<<1)&B[h=H(t,i)];}
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

