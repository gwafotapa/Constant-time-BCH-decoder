/**
 * @file gf_pclmul.c
 * Galois field implementation with multiplication using the pclmulqdq instruction
 */

#include "gf.h"

#include "optimizations.h"

#include <emmintrin.h>
#include <wmmintrin.h>

#ifdef GF_LUT_1024
#include "gf_lut_1024.h"
#else
/**
 * Powers of the chosen primitive element of GF(2^GF_M).
 * The last two elements are needed when multiplication is done by lookup tables
 * (for example if both elements to multiply are zero).
 */
static uint16_t exp[(1 << GF_M) + 2];

/**
 * Logarithm of elements of GF(2^GF_M).
 * The logarithm of 0 is set to 2^GF_M by convention.
 */
static uint16_t log[1 << GF_M];

/**
 * Generates exp and log lookup tables of GF(2^GF_M).
 */
void
gf_generate_M()
{
        gf_generate(exp, log, GF_M);
}
#endif

/**
 * Returns a primitive polynomial of GF(2^m) in hexadecimal form.
 * m must be between 3 and 15 (included).
 * @returns the primitive polynomial of GF(2^m)
 * @param[in] m Parameter of Galois field GF(2^m)
 */
uint16_t
gf_primitive_poly(unsigned int m)
{
        switch(m) {
        case 3:  return 0xB;
        case 4:  return 0x11;
        case 5:  return 0x25;
        case 6:  return 0x43;
        case 7:  return 0x83;
        case 8:  return 0x11D;
        case 9:  return 0x211;
        case 10: return 0x409;
        case 11: return 0x805;
        case 12: return 0x1053;
        case 13: return 0x201B;
        case 14: return 0x4143;
        case 15: return 0x8003;
        default:
                return 0;
        }
}

/**
 * Generates exp and log lookup tables of GF(2^m).
 * The logarithm of 0 is defined as 2^GF_M by convention. <br>
 * The last two elements of the exp table are needed when multiplication is done by lookup tables
 * (for example if both elements to multiply are zero).
 * @param[out] exp Array of size 2^GF_M + 2 receiving the powers of the primitive element
 * @param[out] log Array of size 2^GF_M receiving the logarithms of the elements of GF(2^m)
 * @param[in] m Parameter of Galois field GF(2^m)
 */
void
gf_generate(uint16_t            *exp,
            uint16_t            *log,
            unsigned int         m)
{
        uint16_t elt = 1;
        uint16_t alpha = 2; // primitive element of GF(2^GF_M)
        uint16_t gf_poly = gf_primitive_poly(m);
        
        for (size_t i = 0; i < (1 << m) - 1; ++i){
                exp[i] = elt;
                log[elt] = i;

                elt *= alpha;
                if (elt >= 1 << m)
                        elt ^= gf_poly;
        }

        exp[(1 << m) - 1] = 1;
        exp[1 << m] = 2;
        exp[(1 << m) + 1] = 4;
        log[0] = 1 << m; // by convention
}

/**
 * Returns the requested power of the primitive element of GF(2^GF_M).
 * @returns a^i
 */
uint16_t
gf_exp(uint16_t i)
{
        return exp[i];
}

/**
 * Returns the integer i such that elt = a^i
 * where a is the primitive element of GF(2^GF_M).
 *@returns the logarithm of the given element
 */
uint16_t
gf_log(uint16_t elt)
{
        return log[elt];
}

/**
 * Reduces polynomial x modulo primitive polynomial GF_POLY.
 * @returns x mod GF_POLY
 * @param[in] x Polynomial of degree less than 64
 * @param[in] deg_x The degree of polynomial x
 */
uint16_t
gf_reduce(uint64_t      x,
          size_t        deg_x)
{  
        // Compute the distance between the primitive polynomial first two set bits
        size_t lz1 = __builtin_clz(GF_POLY);
        size_t lz2 = __builtin_clz(GF_POLY ^ 1<<GF_M);
        size_t dist = lz2 - lz1;

        // Deduce the number of steps of reduction
        size_t steps = CEIL_DIV(deg_x - (GF_M-1), dist);

        // Reduce
        for (size_t i = 0; i < steps; ++i) {
                uint64_t mod = x >> GF_M;
                x &= (1<<GF_M) - 1;
                x ^= mod;

                size_t tz1 = 0;
                uint16_t rmdr = GF_POLY ^ 1;
                for (size_t j = __builtin_popcount(GF_POLY)-2; j; --j) {      
                        size_t tz2 = __builtin_ctz(rmdr);
                        size_t shift = tz2 - tz1;
                        mod <<= shift;
                        x ^= mod;
                        rmdr ^= 1 << tz2;
                        tz1 = tz2;
                }
        }

        return x;
}

/**
 * Multiplies two elements of GF(2^GF_M).
 * @returns the product a*b
 * @param[in] a Element of GF(2^GF_M)
 * @param[in] b Element of GF(2^GF_M)
 */
uint16_t
gf_mul(uint16_t a,
       uint16_t b)
{
        __m128i va = _mm_cvtsi32_si128(a);
        __m128i vb = _mm_cvtsi32_si128(b);
        __m128i vab = _mm_clmulepi64_si128(va, vb, 0);
        uint32_t ab = _mm_cvtsi128_si32(vab);

        return gf_reduce(ab, 2*(GF_M-1));
}

/**
 * Squares an element of GF(2^GF_M).
 * @returns a^2
 * @param[in] a Element of GF(2^GF_M)
 */
uint16_t
gf_square(uint16_t a)
{
        uint32_t b = a;
        uint32_t s = b & 1;     
        for (size_t i = 1; i < GF_M; ++i) {
                b <<= 1;
                s ^= b & (1 << 2*i);
        }

        return gf_reduce(s, 2*(GF_M-1));
}

#ifdef GF_ARITHMETIC_1024 // GF inverse implementation specific to GF(2^10)

/**
 * Computes the 4th power of an element of GF(2^GF_M).
 * @returns a^4
 * @param[in] a Element of GF(2^GF_M)
 */
uint16_t
gf_quad(uint64_t a)
{
        uint64_t q = a & 1;
        for (size_t i = 1; i < GF_M; ++i) {
                a <<= 3;
                q ^= a & (1ull << 4*i);
        }

        return gf_reduce(q, 4*(GF_M-1));
}

/**
 * Computes the inverse of an element of GF(2^10),
 * using a shorter chain of squares and multiplications than fast exponentiation.
 * @returns the inverse of a
 * @param[in] a Element of GF(2^10)
 */
uint16_t
gf_inverse(uint16_t a)
{  
        uint16_t p;
        uint16_t a2;

        a2 = gf_square(a);  // a^2
        a = gf_mul(a2, a);  // a^2.a
        p = gf_quad(a);     // a^8.a^4  
        a = gf_mul(p, a);   // a^8.a^4.a^2.a
        p = gf_quad(a);     // a^32.a^16.a^8.a^4
        p = gf_quad(p);     // a^128.a^64.a^32.a^16
        a = gf_mul(p, a);   // a^128.a^64.a^32.a^16.a^8.a^4.a^2.a
        p = gf_quad(a);     // a^512.a^256.a^128.a^64.a^32.a^16.a^8.a^4
        p = gf_mul(p, a2);  // a^-1

        return p;
}

#else // !defined GF_ARITHMETIC_1024

/**
 * Computes the inverse of an element of GF(2^GF_M) by fast exponentiation.
 * @returns the inverse of a
 * @param[in] a Element of GF(2^GF_M)
 */
uint16_t
gf_inverse(uint16_t a)
{
        size_t pow = (1 << GF_M) - 2;
        uint16_t inv = 1;

        do {
                if (pow & 1)
                        inv = gf_mul(inv, a);
                a = gf_square(a);
                pow >>= 1;
        } while (pow);

        return inv;
}

#endif // !defined GF_ARITHMETIC_1024

/**
 * Returns i modulo 2^GF_M-1.
 * i must be less than 2*(2^GF_M-1).
 * Therefore, the return value is either i or i-2^GF_M+1.
 * @returns i mod (2^GF_M-1)
 * @param[in] i The integer whose modulo is taken
 */
uint16_t
gf_mod(uint16_t i)
{
        uint16_t tmp = i - GF_MUL_ORDER;
  
        // mask = 0xffff if (i < GF_MUL_ORDER)
        int16_t mask = -(tmp >> 15);
  
        return tmp + (mask & GF_MUL_ORDER);
}
