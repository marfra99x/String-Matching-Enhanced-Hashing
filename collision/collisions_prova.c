#include "include/define.h"
#include "include/main.h" 
#include <smmintrin.h>
#include <inttypes.h>
#include <immintrin.h>


#define TYPE 64
#define CON 16
#define MASK ((1<<CON)-1)

#define Q 8
#define cc uint16_t
#define HASH(XX, II) crc_using_pclmulqdq_optimized(((const uint8_t *)((unsigned char *)(XX) + (II) - 5)), 16)
#define READ(XX,II) *((uint64_t*)(XX+II-Q+1))


#define LIM 256*256


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


int search(unsigned char *x, int m, unsigned char *y, int n) {
    int j, count=0;
    int f[LIM];
    uint16_t h;

    BEGIN_PREPROCESSING
    for (j=0;j<LIM;j++) f[j]=-1;
    END_PREPROCESSING

    BEGIN_SEARCHING
    j=Q-1;
    while (j < n) {
        h = HASH(y,j);
        if (f[h]>-1 && READ(y,j) != READ(y,f[h])) count++;
        if (f[h]<0) f[h]=j;
        j++;
    }
    END_SEARCHING
    return count;
}

