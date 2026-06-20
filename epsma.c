/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 * This is an implementation of the EPSMA algorithm in
 *    M. A. Aydogmus and M. O. Kulekci: Optimizing packed string 
 *    matching on AVX2 platform. In: Revised Selected Papers of 
 *    VECPAR 2018, Lecture Notes in Computer Science, Vol. 11333,
 *    Springer, 2018, pp. 45-61.
 */

#include "include/define.h"
#include "include/main.h"
#include <memory.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>
//#include <wmmintrin.h> 
#include <inttypes.h>

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>


typedef union{
   __m256i* symbol32;
   unsigned char* symbol;
} aTEXT; 

typedef union{
              __m256i  v;  /// 256 bit
        unsigned  int  ui[8];
   unsigned short int  us[16];
        unsigned char  uc[32]; 
}aVectorUnion; 

typedef struct list
{
    struct list *next;
    int pos;
} LIST;


void print256_num(__m256i var)
{
   uint8_t *v = (uint8_t*) &var; 
   printf("Num: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i.\n", 
   v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15], 
   v[16],v[17],v[18],v[19],v[20],v[21],v[22],v[23],v[24],v[25],v[26],v[27],v[28],v[29],v[30],v[31] ); 
} 
void print256_char(__m256i var)
{
   uint8_t *v = (uint8_t*) &var; 
   printf("Char: %c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c.\n", 
   v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15], 
   v[16],v[17],v[18],v[19],v[20],v[21],v[22],v[23],v[24],v[25],v[26],v[27],v[28],v[29],v[30],v[31] ); 
}  

// ==================================================================================================
int search1(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{  // we exactly know patlen=1	
    
    __m256i* text = (__m256i*)x; 
    __m256i* endTx  = (__m256i*)(x+32*(textlen/32));
	
    int count=0; 
    __m256i t0, a ;
    aVectorUnion template0;
    unsigned int j, k, d_k, indx;

   BEGIN_PREPROCESSING
    for(j=0; j<32; j++)
    {
        template0.uc[j]=pattern[0];  
    }
    t0 = template0.v; 
   END_PREPROCESSING

   BEGIN_SEARCHING	
   while(text<endTx){
        a  =  _mm256_cmpeq_epi8( t0,*text );  
        j  =  _mm256_movemask_epi8(a);   

        count  += _mm_popcnt_u32( j ); 
        text++; // next chunk
     } 
	
   // now we are at the beginning of the last 32-byte block, perform naive check
    for(j=32*(textlen/32); j<textlen; j++){ 
		count += (x[j]==pattern[0]);
    }
   END_SEARCHING
  
    return count;
} 

// ==================================================================================================
int search2(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{  // we exactly know patlen=2

    __m256i* text = (__m256i*)x; 
    __m256i* endTx  = (__m256i*)(x+32*(textlen/32));
	
    __m256i t0,t1,a,b ;
    aVectorUnion template0,template1;
    unsigned int j,k,carry=0 ; 
    int count=0;
    unsigned char firstch = pattern[0], lastch = pattern[1];

   BEGIN_PREPROCESSING  	
    for(j=0; j<32; j++) 
    {
        template0.uc[j]=firstch;        
        template1.uc[j]=lastch;        
    }
    t0 = template0.v;
    t1 = template1.v;
   END_PREPROCESSING   
		 
   BEGIN_SEARCHING  
    while( text < endTx ){
        a  =  _mm256_cmpeq_epi8( t0,*text );   
        j  =  _mm256_movemask_epi8(a);  
  
        b  = _mm256_cmpeq_epi8(t1,*text);
        k  = _mm256_movemask_epi8(b);

        count  += _mm_popcnt_u32(  ( (j<<1)|(carry>>31) )   & k );
        carry = j & 0x80000000;  

        text++; 
    }  
   
    // now we are at the beginning of the last 32-byte block, perform naive check
    for(j=32*(textlen/32);j<textlen;j++) 
         count += ((x[j-1]==firstch) && (x[j]==lastch)); 
	 
  END_SEARCHING

    return count;
} 

// ==================================================================================================
int search3(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{  // we exactly know patlen=3

   BEGIN_PREPROCESSING 
    __m256i* text = (__m256i*)x; 
    __m256i* endTx  = (__m256i*)(x+32*(textlen/32));

    __m256i t0,t1,t2,a,b,c;
    aVectorUnion template0, template1, template2;
    unsigned int j,k,l,dj,dk,dl, carry0=0, carry1=0, dc0=0,dc1=1;
    int count=0;

    for(j=0; j<32; j++) 
    {
        template0.uc[j]=pattern[0];       
        template1.uc[j]=pattern[1];         
        template2.uc[j]=pattern[2];     
    }
    t0 = template0.v; 
    t1 = template1.v;
    t2 = template2.v;
  END_PREPROCESSING  


  BEGIN_SEARCHING 
    while(text<endTx){
        a     = _mm256_cmpeq_epi8(t0,*text);
        j     = _mm256_movemask_epi8(a);

        b     = _mm256_cmpeq_epi8(t1,*text);
        k     = _mm256_movemask_epi8(b);

        c     = _mm256_cmpeq_epi8(t2,*text);
        l     = _mm256_movemask_epi8(c);   

	    count  += _mm_popcnt_u32(  ((j<<2)|(carry0>>30)) & ((k<<1)|(carry1>>31)) & l );

	    carry0 = j & 0xC0000000; 
       carry1 = k & 0x80000000;
	
       text++;
    }

    // now we are at the beginning of the last 32-byte block, perform naive check
    for(j = 32*(textlen/32); j<textlen; j++) 
        count += ((x[j-2]==pattern[0]) && (x[j-1]==pattern[1]) && (x[j]==pattern[2]));
   END_SEARCHING 

    return count;
} 

// ==================================================================================================
int search4_M1(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{//we exactly know patlen=4

    __m256i* text = (__m256i*)x; 
    __m256i* endTx  = (__m256i*)(x+32*(textlen/32));

   BEGIN_PREPROCESSING
    __m256i t0,t1,t2,t3, a,b,c,d;
    aVectorUnion template0,template1,template2,template3;
    unsigned int j,k,l,m, carry0=0,carry1=0,carry2=0;
    int count=0;

    for(j=0; j<32; j++) 
    {
        template0.uc[j]=pattern[0];       
        template1.uc[j]=pattern[1];         
        template2.uc[j]=pattern[2];  
	    template3.uc[j]=pattern[3];    
    }
    t0 = template0.v; 
    t1 = template1.v;
    t2 = template2.v;
    t3 = template3.v;
  END_PREPROCESSING

  BEGIN_SEARCHING
    while(text<endTx){
        a     = _mm256_cmpeq_epi8(t0,*text);
        j     = _mm256_movemask_epi8(a);

        b     = _mm256_cmpeq_epi8(t1,*text);
        k     = _mm256_movemask_epi8(b);

        c     = _mm256_cmpeq_epi8(t2,*text);
        l     = _mm256_movemask_epi8(c); 

        d     = _mm256_cmpeq_epi8(t3,*text);
        m     = _mm256_movemask_epi8(d);   

	   count  += _mm_popcnt_u32(  ( (j<<3)|(carry0>>29)) & ((k<<2)|(carry1>>30)) 
							   	& ( (l<<1)|(carry2>>31)) & m );
	    carry0 =  j & 0xE0000000;
       carry1 =  k & 0xC0000000; 
       carry2 =  l & 0x80000000;
	
       text++;
	}

    // now we are at the beginning of the last 32-byte block, perform naive check
    for(j=32*(textlen/32);j<textlen;j++) 
        count += ((x[j-3]==pattern[0]) && (x[j-2]==pattern[1]) && 
        		  (x[j-1]==pattern[2]) && (x[j]==pattern[3]));
  END_SEARCHING

    return count;
}

// ==================================================================================================
// Type 2: Alternative Method 
int search4_M2(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    unsigned char* y0, ytemp;
    int i,j,k,count=1, lpix=0, ix;

     const __m256i permute = _mm256_setr_epi32(0, 1, 2, 0,  
     										   2, 3, 4, 0);
     __m256i zeros  = _mm256_setzero_si256();
     __m256i res, a, b, z, p256, arrText2, data, tempData, cmp;
     aVectorUnion P;
	
     int tempTxtLen = 32 * (textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);
     endTx--;

     int residualTxt = textlen - tempTxtLen; 
     uint32_t mask, result ; 

    BEGIN_PREPROCESSING  
     P.uc[0] = pattern[0];  P.uc[1] = pattern[1];  P.uc[2] = pattern[2];  P.uc[3] = pattern[3]; 
     P.uc[4] = pattern[4];  P.uc[5] = pattern[5];  P.uc[6] = pattern[6];  P.uc[7] = pattern[7];

     P.uc[16] = pattern[0];  P.uc[17] = pattern[1];  P.uc[18] = pattern[2];  P.uc[19] = pattern[3]; 
     P.uc[20] = pattern[4];  P.uc[21] = pattern[5];  P.uc[22] = pattern[6];  P.uc[23] = pattern[7];

     p256 = P.v;
    END_PREPROCESSING

  
    BEGIN_SEARCHING
     while(text < endTx)
     { 
	     /// First 16 char part 
         data   = _mm256_permutevar8x32_epi32(*text, permute);
         res    = _mm256_mpsadbw_epu8(data, p256,0);
         cmp    = _mm256_cmpeq_epi16(res, zeros);
         mask   = _mm256_movemask_epi8(cmp) ;
         count += _mm_popcnt_u32(mask);

	     /// Second 16 char part 
         tempData = _mm256_permute2f128_si256  (*text, *(text+1), 33); 
         data   = _mm256_permutevar8x32_epi32(tempData, permute);
	     res    = _mm256_mpsadbw_epu8(data, p256,0);
         cmp    = _mm256_cmpeq_epi16(res, zeros);
         mask =  _mm256_movemask_epi8(cmp) ;
         count += _mm_popcnt_u32(mask);
	
         text++;
     } 
     count = count / 2;

    // now we are at the beginning of the last 32-byte block, perform naive check
    for(j = (tempTxtLen-32); j < textlen; j++) {
        count += ( (x[j-3]==pattern[0]) && (x[j-2]==pattern[1]) 
        	   &&  (x[j-1]==pattern[2]) && (x[j]==pattern[3])  );  
    }
	
    END_SEARCHING 

     return count;
} 

// ==================================================================================================
int search5(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    BEGIN_PREPROCESSING 
     
     unsigned char* y0, ytemp;
     int i,j,k,count=0, lpix=0, ix;
     aVectorUnion P,zero;
     __m256i res,a,b,z,p;

     int tempTxtLen = 32*(textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);
     int residualTxt= textlen - tempTxtLen; 

     __m256i t0,t1,t2,t3,t4, sa,sb,sc,sd,se ; 
     aVectorUnion template0,template1,template2,template3,template4; 
     uint32_t sj,sk,sl,sm, sn,so,sp,sr,st, carry0=0,carry1=0,carry2=0,carry3=0;

     for(j=0; j<32; j++)  
     {
        template0.uc[j]= pattern[0];       
        template1.uc[j]= pattern[1];         
        template2.uc[j]= pattern[2];  
        template3.uc[j]= pattern[3];
	    template4.uc[j]= pattern[4];
     }
     t0 = template0.v; 
     t1 = template1.v;
     t2 = template2.v;
     t3 = template3.v;
     t4 = template4.v;
    END_PREPROCESSING

    BEGIN_SEARCHING
     while(text < endTx)
     {
        sa     = _mm256_cmpeq_epi8(t0,*text);
        sj     = _mm256_movemask_epi8(sa);

        sb     = _mm256_cmpeq_epi8(t1,*text);
        sk     = _mm256_movemask_epi8(sb); 

        sc     = _mm256_cmpeq_epi8(t2,*text);
        sl     = _mm256_movemask_epi8(sc); 

        sd     = _mm256_cmpeq_epi8(t3,*text);
        sm     = _mm256_movemask_epi8(sd); 

        se     = _mm256_cmpeq_epi8(t4,*text);
        sn     = _mm256_movemask_epi8(se);

	    count  += _mm_popcnt_u32( ((sj<<4)|(carry0>>28)) & ( (sk<<3)|(carry1>>29) ) 
				         & ((sl<<2)|(carry2>>30)) & ( (sm<<1)|(carry3>>31) ) & sn );

        carry0 =  sj & 0xF0000000;   carry1 =  sk & 0xE0000000;  
        carry2 =  sl & 0xC0000000;   carry3 =  sm & 0x80000000;  
	
        text++;
     } 

  /// now we are at the beginning of the last 32-byte block, perform naive check
    for(j=(tempTxtLen-32); j<textlen; j++) {
        count += (  (x[j-4]==pattern[0]) && (x[j-3]==pattern[1]) && (x[j-2]==pattern[2]) 
	         	 && (x[j-1]==pattern[3]) && (x[j]==pattern[4]) );  
    }
    
	END_SEARCHING

     return count;
} 

// ==================================================================================================
int search6(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
   
     unsigned char* y0, ytemp;
     int i,j,k,count=0, lpix=0, ix;
     aVectorUnion P,zero;
     __m256i res,a,b,z,p;

     int tempTxtLen = 32*(textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);
     int residualTxt= textlen - tempTxtLen; 

     __m256i t0,t1,t2,t3, t4,t5 , sa, sb, sc, sd, se, sf ;  
     aVectorUnion template0,template1,template2,template3,template4,template5 ;
     uint32_t sj,sk,sl,sm, sn,so,sp,sr,st, carry0=0,carry1=0,carry2=0,carry3=0,carry4=0;

    BEGIN_PREPROCESSING  
 	
     for(j=0; j<32; j++)  
     {
        template0.uc[j]= pattern[0];       
        template1.uc[j]= pattern[1];         
        template2.uc[j]= pattern[2];  
        template3.uc[j]= pattern[3];
	     template4.uc[j]= pattern[4];
	     template5.uc[j]= pattern[5];
     }
     t0 = template0.v; 
     t1 = template1.v;
     t2 = template2.v;
     t3 = template3.v;
     t4 = template4.v;
     t5 = template5.v;
    END_PREPROCESSING

 
    BEGIN_SEARCHING

     while(text < endTx)
     {
        sa     = _mm256_cmpeq_epi8(t0,*text);
        sj     = _mm256_movemask_epi8(sa);

        sb     = _mm256_cmpeq_epi8(t1,*text);
        sk     = _mm256_movemask_epi8(sb); 

        sc     = _mm256_cmpeq_epi8(t2,*text);
        sl     = _mm256_movemask_epi8(sc); 

        sd     = _mm256_cmpeq_epi8(t3,*text);
        sm     = _mm256_movemask_epi8(sd); 

        se     = _mm256_cmpeq_epi8(t4,*text);
        sn     = _mm256_movemask_epi8(se);

	    sf     = _mm256_cmpeq_epi8(t5,*text);
        so     = _mm256_movemask_epi8(sf);

        count  += _mm_popcnt_u32(  ((sj<<5)|(carry0>>27)) & ( (sk<<4)|(carry1>>28) ) & ((sl<<3)|(carry2>>29)) & ((sm<<2)|(carry3>>30)) & ( (sn<<1)|(carry4>>31) ) &  so );

        carry0 =  sj & 0xF8000000;  carry1 =  sk & 0xF0000000;  
        carry2 =  sl & 0xE0000000;  carry3 =  sm & 0xC0000000;  
        carry4 =  sn & 0x80000000;  
	
        text++;
     } 

    /// now we are at the beginning of the last 32-byte block, perform naive check
   for(j=(tempTxtLen-32); j<textlen; j++) {
        count += ( (x[j-5]==pattern[0]) && (x[j-4]==pattern[1]) && (x[j-3]==pattern[2]) 
	         	&& (x[j-2]==pattern[3]) && (x[j-1]==pattern[4]) && (x[j]==pattern[5]) );  
    }
		
	END_SEARCHING 
	
	return count;
}


// ==================================================================================================
int search8(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
     unsigned char* y0, ytemp;
     int i,j,k,count=0, lpix=0, ix;
     aVectorUnion P,zero;
     __m256i res,a,b,z,p;

     int tempTxtLen = 32*(textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);
     int residualTxt= textlen - tempTxtLen; 

 	 __m256i t0,t1,t2,t3, t4,t5, t6, t7, sa, sb, sc, sd, se, sf, sg, sh ;  
 	 aVectorUnion template0,template1,template2,template3,template4,template5,template6,template7;
	 uint32_t sj,sk,sl,sm,  sn,so,sp,sr,st ;
	 uint32_t carry0=0,carry1=0,carry2=0,carry3=0,carry4=0,carry5=0,carry6=0;

    BEGIN_PREPROCESSING  

     for(j=0; j<32; j++)  
     {
        template0.uc[j]= pattern[0];       
        template1.uc[j]= pattern[1];         
        template2.uc[j]= pattern[2];  
        template3.uc[j]= pattern[3];
        template4.uc[j]= pattern[4];
        template5.uc[j]= pattern[5];
        template6.uc[j]= pattern[6];
        template7.uc[j]= pattern[7];
     }
     t0 = template0.v; 
     t1 = template1.v;
     t2 = template2.v;
     t3 = template3.v;
     t4 = template4.v;
     t5 = template5.v;
     t6 = template6.v;
     t7 = template7.v;
     
     END_PREPROCESSING

 
     BEGIN_SEARCHING

     while(text < endTx)
     {
        sa     = _mm256_cmpeq_epi8(t0,*text);
        sj     = _mm256_movemask_epi8(sa);

        sb     = _mm256_cmpeq_epi8(t1,*text);
        sk     = _mm256_movemask_epi8(sb); 

        sc     = _mm256_cmpeq_epi8(t2,*text);
        sl     = _mm256_movemask_epi8(sc); 

        sd     = _mm256_cmpeq_epi8(t3,*text);
        sm     = _mm256_movemask_epi8(sd); 

        se     = _mm256_cmpeq_epi8(t4,*text);
        sn     = _mm256_movemask_epi8(se);

	    sf     = _mm256_cmpeq_epi8(t5,*text);
        so     = _mm256_movemask_epi8(sf);
		
		sg     = _mm256_cmpeq_epi8(t6,*text);
        sp     = _mm256_movemask_epi8(sg);
		
		sh     = _mm256_cmpeq_epi8(t7,*text);
		sr     = _mm256_movemask_epi8(sh);
								  
	    count += _mm_popcnt_u32( ((sj<<7)|(carry0>>25)) & ((sk<<6)|(carry1>>26)) 
			      	 & ((sl<<5)|(carry2>>27)) & ((sm<<4)|(carry3>>28)) & ((sn<<3)|(carry4>>29))  
			         & ((so<<2)|(carry5>>30)) & ((sp<<1)|(carry6>>31)) &   sr );

        carry0=sj&0xFE000000; carry1=sk&0xFC000000; carry2=sl&0xF8000000; carry3=sm&0xF0000000; 
        carry4=sn&0xE0000000; carry5=so&0xC0000000; carry6=sp&0x80000000;
		
        text++;
     } 
 
	 /// now we are at the beginning of the last 32-byte block, perform naive check
     for(j=(tempTxtLen-32); j<textlen; j++) {
        count += (  (x[j-7]==pattern[0]) && (x[j-6]==pattern[1]) && (x[j-5]==pattern[2]) 
	         	&& (x[j-4]==pattern[3]) && (x[j-3]==pattern[4]) && (x[j-2]==pattern[5])
	         	&& (x[j-1]==pattern[6]) && (x[j]==pattern[7]) );  
     }
    
	 END_SEARCHING 

     return count;
}

// ==================================================================================================
int  search8_15(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{	  
     unsigned char* y0 ;
     int i,j,k, count=0, ptlnVal = patlen-8;
     uint32_t mj, nj; 

     const __m256i permute  = _mm256_setr_epi32(0, 1, 2, 0,   2, 3, 4, 0);
     const __m256i permute2 = _mm256_setr_epi32(1, 2, 3, 0,   3, 4, 5, 0);
     __m256i zeros  = _mm256_setzero_si256();

     aVectorUnion P, P2;
     __m256i res, a, b, z, p256, p256_2, data, tempData, cmp;

     int tempTxtLen = 32*(textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);

     int residualTxt= textlen - tempTxtLen; 
	 
    BEGIN_PREPROCESSING  
     // 4 + 4 char 
     P.uc[0] = pattern[0];   P.uc[1] = pattern[1];   P.uc[2] = pattern[2];   P.uc[3]  = pattern[3]; 
     P.uc[16] = pattern[0];  P.uc[17] = pattern[1];  P.uc[18] = pattern[2];  P.uc[19] = pattern[3]; 
     p256 = P.v;

     P2.uc[0] = pattern[4];   P2.uc[1] = pattern[5];   P2.uc[2] = pattern[6];   P2.uc[3] = pattern[7]; 
     P2.uc[16] = pattern[4];  P2.uc[17] = pattern[5];  P2.uc[18] = pattern[6];  P2.uc[19] = pattern[7]; 
     p256_2 = P2.v;
     
    END_PREPROCESSING
  
  
    BEGIN_SEARCHING
    
    while(text < endTx)
     { 
	    ///--- First 16 byte chunk --- 
        data  =  _mm256_permutevar8x32_epi32(*text, permute);
        res   =  _mm256_mpsadbw_epu8(data, p256,0);
        cmp   =  _mm256_cmpeq_epi16(res, zeros);
        nj    =  _mm256_movemask_epi8(cmp) ;
    
	    if (nj>0) 
       {      
          y0 = (unsigned char*)(text); 

	      if ( (nj&3221225472)==3221225472 && !memcmp(pattern,y0+15,patlen) ) {count++;}
	       if ( (nj&805306368)==805306368 && !memcmp(pattern,y0+14,patlen) ) {count++;}
	       if ( (nj&201326592)==201326592 && !memcmp(pattern,y0+13,patlen) ) {count++;}
	       if ( (nj&50331648)==50331648 && !memcmp(pattern,y0+12,patlen) )  {count++;}
	       if ( (nj&12582912)==12582912 && !memcmp(pattern,y0+11,patlen) ) {count++;}
	       if ( (nj&3145728)==3145728 && !memcmp(pattern,y0+10,patlen) ) {count++;}
	       if ( (nj&786432)==786432 && !memcmp(pattern,y0+9,patlen) ) {count++;}
	       if ( (nj&196608)==196608 && !memcmp(pattern,y0+8,patlen) ) {count++;}
	       if ( (nj&49152)==49152 && !memcmp(pattern,y0+7,patlen) ) {count++;}
	       if ( (nj&12288)==12288 && !memcmp(pattern,y0+6,patlen) ) {count++;}
	       if ( (nj&3072)==3072 && !memcmp(pattern,y0+5,patlen) ) {count++;}
	       if ( (nj&768)==768 && !memcmp(pattern,y0+4,patlen) ) {count++;}
	       if ( (nj&192)==192 && !memcmp(pattern,y0+3,patlen) ) {count++;}
	       if ( (nj&48)==48 && !memcmp(pattern,y0+2,patlen) ) {count++;}
	       if ( (nj&12)==12 && !memcmp(pattern,y0+1,patlen) ) {count++;}
	       if ( (nj&3)==3 && !memcmp(pattern,y0,patlen) ) {count++;} 
       } 
 
        ///--- Second 16 byte chunk --- 
        tempData = _mm256_permute2f128_si256  (*text, *(text+1), 33); 
        data  = _mm256_permutevar8x32_epi32(tempData, permute);
	    res   = _mm256_mpsadbw_epu8(data, p256, 0);
        cmp   = _mm256_cmpeq_epi16(res, zeros);
        nj = _mm256_movemask_epi8(cmp) ;

	    if (nj>0)
       {  
          y0 = (unsigned char*)(text); 

	       if ( (nj&3221225472)==3221225472 && !memcmp(pattern,y0+31,patlen) ) {count++;}
	       if ( (nj&805306368)==805306368 && !memcmp(pattern,y0+30,patlen) ) {count++;}
	       if ( (nj&201326592)==201326592 && !memcmp(pattern,y0+29,patlen) ) {count++;}
	       if ( (nj&50331648)==50331648 && !memcmp(pattern,y0+28,patlen) )  {count++;}
	       if ( (nj&12582912)==12582912 && !memcmp(pattern,y0+27,patlen) ) {count++;}
	       if ( (nj&3145728)==3145728 && !memcmp(pattern,y0+26,patlen) ) {count++;}
	       if ( (nj&786432)==786432 && !memcmp(pattern,y0+25,patlen) ) {count++;}
	       if ( (nj&196608)==196608 && !memcmp(pattern,y0+24,patlen) ) {count++;}
	       if ( (nj&49152)==49152 && !memcmp(pattern,y0+23,patlen) ) {count++;}
	       if ( (nj&12288)==12288 && !memcmp(pattern,y0+22,patlen) ) {count++;}
	       if ( (nj&3072)==3072 && !memcmp(pattern,y0+21,patlen) ) {count++;}
	       if ( (nj&768)==768 && !memcmp(pattern,y0+20,patlen) ) {count++;}
	       if ( (nj&192)==192 && !memcmp(pattern,y0+19,patlen) ) {count++;}
	       if ( (nj&48)==48 && !memcmp(pattern,y0+18,patlen) ) {count++;}
	       if ( (nj&12)==12 && !memcmp(pattern,y0+17,patlen) ) {count++;}
	       if ( (nj&3)==3 && !memcmp(pattern,y0+16,patlen) ) {count++;} 
       }  

       text++;
    }

  END_SEARCHING 

    return count;
} 


int search16_32(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
     unsigned char* y0, ytemp;
     int i,j,k,count=0, ptlnVal = patlen-8, ci=0; // ptlnVal: reduced plen
     uint32_t mj,nj ; 

     const __m256i permute  = _mm256_setr_epi32(0, 1, 2, 0,   2, 3, 4, 0);
     const __m256i permute2 = _mm256_setr_epi32(1, 2, 3, 0,   3, 4, 5, 0);
     __m256i zeros  = _mm256_setzero_si256();

     aVectorUnion P,P2;
     __m256i res, a, b, z, p256, p256_2, data, tempData, cmp;

     int tempTxtLen = 32*(textlen/32); 
     __m256i* text = (__m256i*)x; 
     __m256i* endTx  = (__m256i*)(x+tempTxtLen);
     endTx--;
     int residualTxt= textlen - tempTxtLen; 
	 
	 // int counter1=0, counter2=0, counter3=0, counter4=0; 
	 
     
	BEGIN_PREPROCESSING  
		
     P.uc[0] = pattern[0];   P.uc[1] = pattern[1];   P.uc[2] = pattern[2];   P.uc[3] = pattern[3]; 
     P.uc[16] = pattern[0];  P.uc[17] = pattern[1];  P.uc[18] = pattern[2];  P.uc[19] = pattern[3]; 
     p256 = P.v;

     P2.uc[0] = pattern[4];   P2.uc[1] = pattern[5];   P2.uc[2] = pattern[6];   P2.uc[3] = pattern[7]; 
     P2.uc[16] = pattern[4];  P2.uc[17] = pattern[5];  P2.uc[18] = pattern[6];  P2.uc[19] = pattern[7]; 
     p256_2 = P2.v;
     
    END_PREPROCESSING


    BEGIN_SEARCHING
     while(text < endTx )
     { 
       ///--- First 16 byte chunk --- 
       data   = _mm256_permutevar8x32_epi32(*text, permute);
       res    = _mm256_mpsadbw_epu8(data, p256,0);
       cmp    = _mm256_cmpeq_epi16(res, zeros);
       mj = _mm256_movemask_epi8(cmp) ;
    
	   if (mj>0) 
       {    
          /// --- Search for next 4 char
          data  = _mm256_permutevar8x32_epi32(*text, permute2);  
          res   = _mm256_mpsadbw_epu8(data, p256_2, 0);
          cmp   = _mm256_cmpeq_epi16(res, zeros);
	      nj = mj & ( _mm256_movemask_epi8(cmp) );

          if (nj>0)
          { 
             y0 = (unsigned char*)(text); 

		     if ( (nj&3)==3 && !memcmp(&pattern[8],y0+8,ptlnVal) ) {count++;}  	
			 if ( (nj&12)==12 && !memcmp(&pattern[8],y0+9,ptlnVal) ) {count++;}
		     if ( (nj&48)==48 && !memcmp(&pattern[8],y0+10,ptlnVal) ) {count++;}
			 if ( (nj&192)==192 && !memcmp(&pattern[8],y0+11,ptlnVal) ) {count++;}
			 if ( (nj&768)==768 && !memcmp(&pattern[8],y0+12,ptlnVal) ) {count++;}
			 if ( (nj&3072)==3072 && !memcmp(&pattern[8],y0+13,ptlnVal) ) {count++;}
			 if ( (nj&12288)==12288 && !memcmp(&pattern[8],y0+14,ptlnVal) ) {count++;}
			 if ( (nj&49152)==49152 && !memcmp(&pattern[8],y0+15,ptlnVal) ) {count++;}
			 if ( (nj&196608)==196608 && !memcmp(&pattern[8],y0+16,ptlnVal) ) {count++;}
			 if ( (nj&786432)==786432 && !memcmp(&pattern[8],y0+17,ptlnVal) ) {count++;}
			 if ( (nj&3145728)==3145728 && !memcmp(&pattern[8],y0+18,ptlnVal) ) {count++;}
			 if ( (nj&12582912)==12582912 && !memcmp(&pattern[8],y0+19,ptlnVal) ) {count++;}
			 if ( (nj&50331648)==50331648 && !memcmp(&pattern[8],y0+20,ptlnVal) )  {count++;}
			 if ( (nj&201326592)==201326592 && !memcmp(&pattern[8],y0+21,ptlnVal) ) {count++;}
			 if ( (nj&805306368)==805306368 && !memcmp(&pattern[8],y0+22,ptlnVal) ) {count++;}
	         if ( (nj&3221225472)==3221225472 && !memcmp(&pattern[8],y0+23,ptlnVal) ) {count++;} 	         
           }
        }  

	    // ____________________________
        ///--- Second 16 byte chunk --- 
        tempData = _mm256_permute2f128_si256(*text, *(text+1), 33); 
        data  = _mm256_permutevar8x32_epi32(tempData, permute);
	    res   = _mm256_mpsadbw_epu8(data, p256, 0);
        cmp   = _mm256_cmpeq_epi16(res, zeros);
        mj = _mm256_movemask_epi8(cmp) ;

	    if (mj>0)
        {
          /// --- Search for next 4 char 
          data   = _mm256_permutevar8x32_epi32(tempData, permute2);  
          res  = _mm256_mpsadbw_epu8(data, p256_2, 0);
          cmp    = _mm256_cmpeq_epi16(res, zeros);
	      nj = mj & ( _mm256_movemask_epi8(cmp) );
 
   	      if (nj>0)
          {      
             y0 = (unsigned char*)(text); 
			   
	         if ( (nj&3)==3 && !memcmp(&pattern[8],y0+24,ptlnVal) ) {count++;} 
			 if ( (nj&12)==12 && !memcmp(&pattern[8],y0+25,ptlnVal) ) {count++;}
			 if ( (nj&48)==48 && !memcmp(&pattern[8],y0+26,ptlnVal) ) {count++;}
			 if ( (nj&192)==192 && !memcmp(&pattern[8],y0+27,ptlnVal) ) {count++;} 
			 if ( (nj&768)==768 && !memcmp(&pattern[8],y0+28,ptlnVal) ) {count++;}
			 if ( (nj&3072)==3072 && !memcmp(&pattern[8],y0+29,ptlnVal) ) {count++;}
			 if ( (nj&12288)==12288 && !memcmp(&pattern[8],y0+30,ptlnVal) ) {count++;}
		     if ( (nj&49152)==49152 && !memcmp(&pattern[8],y0+31,ptlnVal) ) {count++;}
			 if ( (nj&196608)==196608 && !memcmp(&pattern[8],y0+32,ptlnVal) ) {count++;}
			 if ( (nj&786432)==786432 && !memcmp(&pattern[8],y0+33,ptlnVal) ) {count++;}
			 if ( (nj&3145728)==3145728 && !memcmp(&pattern[8],y0+34,ptlnVal) ) {count++;}
		   	 if ( (nj&12582912)==12582912 && !memcmp(&pattern[8],y0+35,ptlnVal) ) {count++;}
			 if ( (nj&50331648)==50331648 && !memcmp(&pattern[8],y0+36,ptlnVal) ) {count++;}
			 if ( (nj&201326592)==201326592 && !memcmp(&pattern[8],y0+37,ptlnVal) ) {count++;}
			 if ( (nj&805306368)==805306368 && !memcmp(&pattern[8],y0+38,ptlnVal) ) {count++;}
	         if ( (nj&3221225472)==3221225472 && !memcmp(&pattern[8],y0+39,ptlnVal) ) {count++;}
	        }
        }  

      text++;
    }
  
  END_SEARCHING
  
  return count;
}  


// ==================================================================================================
int search(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    int count=0 ;

    if (patlen<2)        	    return search1(pattern, patlen, x, textlen); 
    if (patlen==2)       	    return search2(pattern, patlen, x, textlen); 
    if (patlen==2)       	    return search3(pattern, patlen, x, textlen); 
    if (patlen==4)       	    return search4_M1(pattern, patlen, x, textlen); 
    if (patlen==5)       	    return search5(pattern, patlen, x, textlen);  
    if (patlen==6)              return search6(pattern, patlen, x, textlen); 
   // if (patlen==8)              return search8(pattern, patlen, x, textlen); 	
   
    if (7<patlen && patlen<16)   return search8_15(pattern, patlen, x, textlen);
    if (15<patlen && patlen<33)  return search16_32(pattern, patlen, x, textlen); 

    return count;
}
