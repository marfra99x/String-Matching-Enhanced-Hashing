#include "include/define.h"
#include "include/main.h" // defines the search interface for time and statistics.
#include <smmintrin.h>
#include <inttypes.h>
#include <immintrin.h>

#define Q 8
#define c8 uint64_t
#define HASH(XX, II) crc_using_pclmulqdq_optimized(((const uint8_t *)((unsigned char *)(XX) + (II) - 7)), 16)


uint64_t crc_using_pclmulqdq(const uint8_t *data, size_t len) {
    __m128i crc = _mm_setzero_si128(); // Initial CRC value

    __m128i polynomial = _mm_set_epi64x(0, 0x1DB71064); // CRC32 polynomial
    
    while (len >= 8) {
        // Load next 8 bytes of data
        __m128i block = _mm_set_epi64x(0, *(uint64_t *)data);
        crc = _mm_xor_si128(crc, block);

        // Carry-less multiply CRC with polynomial
        crc = _mm_clmulepi64_si128(crc, polynomial, 0x00);

        data += 8;
        len -= 8;
    }

    // Final reduction (simplified, may need more for correctness)
    uint64_t result[2];
    _mm_storeu_si128((__m128i *)result, crc);
    return result[0] ^ result[1];
}

#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

uint64_t crc_using_pclmulqdq_optimized(const uint8_t *data, size_t len) {
    if (len == 0) return 0;

    __m128i crc = _mm_setzero_si128(); // Initial CRC value
    __m128i polynomial = _mm_set_epi64x(0, 0x1DB71064); // CRC32 polynomial

    // Process 16 bytes at a time if possible
    while (len >= 16) {
        __m128i block1 = _mm_loadu_si128((const __m128i *)data); // Load first 16 bytes
        crc = _mm_xor_si128(crc, block1);                       // XOR with CRC
        crc = _mm_clmulepi64_si128(crc, polynomial, 0x00);      // Carry-less multiply

        data += 16;
        len -= 16;
    }

    // Process remaining 8 bytes
    while (len >= 8) {
        uint64_t block = *(const uint64_t *)data;
        __m128i block_vec = _mm_set_epi64x(0, block);
        crc = _mm_xor_si128(crc, block_vec);
        crc = _mm_clmulepi64_si128(crc, polynomial, 0x00);

        data += 8;
        len -= 8;
    }

    // Handle remaining bytes (less than 8)
    if (len > 0) {
        uint64_t tail = 0;
        memcpy(&tail, data, len); // Safely copy remaining bytes
        __m128i tail_vec = _mm_set_epi64x(0, tail);
        crc = _mm_xor_si128(crc, tail_vec);
        crc = _mm_clmulepi64_si128(crc, polynomial, 0x00);
    }

    // Final reduction
    uint64_t result[2];
    _mm_storeu_si128((__m128i *)result, crc);
    return result[0] ^ result[1];
}


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
        _mm_prefetch((const char *)(y + j + 64), _MM_HINT_T0); // Prefetch
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
