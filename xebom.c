#include "include/define.h"
#include "include/main.h"
#define XSIZE       4200			//maximal length of the pattern

int search(unsigned char *x, int m, unsigned char *y, int n) {
   int S[XSIZE], FT[SIGMA*SIGMA];
   int *trans[XSIZE];
   int i, j, p, q;
   int iMinus1, mMinus1, count;
   unsigned char c;
   count = 0;

   BEGIN_PREPROCESSING
   for(i=0; i<=m+1; i++) trans[i] = (int *)malloc (sizeof(int)*(SIGMA));
   for(i=0; i<=m+1; i++) for(j=0; j<SIGMA; j++) trans[i][j]=UNDEFINED;
   S[m] = m + 1;
   for (i = m; i > 0; --i) {
      iMinus1 = i - 1;
      c = x[iMinus1];
      trans[i][c] = iMinus1;
      p = S[i];
      while (p <= m && (q = trans[p][c]) ==  UNDEFINED) {
         trans[p][c] = iMinus1;
         p = S[p];
      }
      S[iMinus1] = (p == m + 1 ? m : q);
   }

   /* Construct the FirstTransition table */
   for(i=0; i<SIGMA; i++) {
      q = trans[m][i];
      for(j=0; j<SIGMA; j++)
         if(q>=0) FT[(i<<8)+j] = trans[q][j];
         else FT[(i<<8)+j] = UNDEFINED;
   }
   END_PREPROCESSING

   BEGIN_SEARCHING
   for(i=0; i<m; i++) y[n+i]=x[i];
   if( !memcmp(x,y,m) ) count++;
   j=m;
   mMinus1 = m-1;
   while (j<n) {
      while ( (FT[*((uint16_t*)(y+j-1))]) == UNDEFINED )j+=mMinus1;
      i = j-2;
      p = FT[*((uint16_t*)(y+j-1))];
      while ( (p = trans[p][y[i]]) != UNDEFINED ) i--;
      if (i < j-mMinus1 && j<n) {
         count++;
         i++;
      }
      j = i + m;
   }
   END_SEARCHING

   for(i=0; i<=m+1; i++) free(trans[i]);
   return count;
}


