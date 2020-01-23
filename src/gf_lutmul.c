/**
 * @file gf_lutmul.c
 * Galois field implementation with multiplication using lookup tables
 */

#include "gf.h"

#include "optimizations.h"

#ifdef GF_LUT_1024
#include "gf_lut_1024.h"
#else
/**
 * Powers of the chosen primitive element of GF(2^GF_M).
 * The last two elements are needed by the gf_mul function from gf_lutmul.c
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
 * The last two elements of the exp table are needed by the gf_mul function from gf_lutmul.c
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
                if(elt >= 1 << m)
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
 * @returns the logarithm of the given element
 */
uint16_t
gf_log(uint16_t elt)
{
        return log[elt];
}

/**
 * Multiplies nonzero element 'a' by element 'b'.
 * @returns the product a*b
 * @param[in] a First element of GF(2^GF_M) to multiply (cannot be zero)
 * @param[in] b Second element of GF(2^GF_M) to multiply (cannot be zero)
 */
uint16_t
gf_mul(uint16_t a,
       uint16_t b)
{
        // mask = 0xffff if neither a nor b is zero. Otherwise mask is 0.
        int16_t mask = ((log[a]|log[b]) >> GF_M) - 1;
        return mask & exp[gf_mod(log[a] + log[b])];
}

/**
 * Squares an element of GF(2^GF_M).
 * @returns a^2
 * @param[in] a Element of GF(2^GF_M)
 */
uint16_t
gf_square(uint16_t a)
{
        int16_t mask = (log[a] >> GF_M) - 1;
        return mask & exp[gf_mod(2*log[a])];
}

/**
 * Computes the inverse of an element of GF(2^GF_M).
 * @returns the inverse of a
 * @param[in] a Element of GF(2^GF_M)
 */
uint16_t
gf_inverse(uint16_t a)
{
        return exp[GF_MUL_ORDER - log[a]];
}

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
