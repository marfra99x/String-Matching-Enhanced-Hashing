// RSSB + QF + CRC64

#include <smmintrin.h>
#include "include/define.h"
#include "include/main.h"

#define LIM 65536
#define W 64
#define Q 8
#define QM 7 //=Q-1
#define H(XX,II) _mm_crc32_u64(54319,*((uint64_t*)(XX+II-QM)))

int search(unsigned char* pat, int m, unsigned char *t, int n) {

   uint16_t h, d, B[LIM];
   int i, j, a, b, qq, x, count = 0;

   if(m<15) return -1;
   
   BEGIN_PREPROCESSING

   for (i=0; i<LIM; i++) B[i]=0;

   x=0; if (m<=4096) x=1;
   a=(m-Q+1-x)/W+x; 
   qq = Q*2;
   if (a<Q*2) {a=Q; qq=Q;}
   if (a%qq>0) a-=(a%qq);
   b=((m-Q+1)/a)*a; 

   for (i=m-b; i<m; i++) 
       B[h=H(pat,i)] |= (1<<((m-i) % qq));

   END_PREPROCESSING

   BEGIN_SEARCHING

   i = m-1;

   while (i<n) {
      d = B[h=H(t,i)];
      if (d==0) i+=b;
      else {
         j = i;
         do {i-=a;
            d = d & B[h=H(t,i)];}
         while (d && (i>j-b));
         i += b;
         if (i==j) {
            if(memcmp(t+i-m+1,pat,m) == 0) count++;
            i=j+1;
         }
      }
    }

   END_SEARCHING
   return count;
}
