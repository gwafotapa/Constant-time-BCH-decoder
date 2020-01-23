/**
 * @file bch_8191_7580_47.h
 * Parameters of the BCH code (8191, 7580, 47)
 */

#ifndef BCH_8191_7580_47_H
#define BCH_8191_7580_47_H

/** Parameter of Galois field GF(2^GF_M) */
#define GF_M             13

/** Order of the Galois field multiplicative group */
#define GF_MUL_ORDER     ((1 << GF_M) - 1)

/** Primitive polynomial the Galois field was built with */
#define GF_POLY          0X201B

/** Length of the code */
#define BCH_N            8191

/** Dimension of the code */
#define BCH_K            7580

/** Correction capacity of the code */
#define BCH_DELTA        47

/** Number of bytes needed to store BCH_N bits */
#define BCH_N_BYTES      CEIL_DIV(BCH_N, 8)

/** Number of bytes needed to store BCH_K bits */
#define BCH_K_BYTES      CEIL_DIV(BCH_K, 8)

/** Number of coefficients of the generator polynomial of the BCH code */
#define BCH_POLY_SIZE    (BCH_N - BCH_K + 1)

/** Logarithm of the smallest power of 2 greater than or equal to BCH_DELTA */
#define FFT_PARAM        6

/** Logarithm of the smallest power of 2 greater than or equal to the number of syndromes */
#define FFT_T_PARAM      (FFT_PARAM + 1)

/** Computes the ceiling of the division of a by b */
#define CEIL_DIV(a, b)   ((a)/(b) + ((a)%(b) == 0 ? 0 : 1))

#endif // BCH_8191_7580_47_H
