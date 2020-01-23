/**
 * @file optimizations.h
 * Optimization flags
 */

#ifndef OPTIMIZATIONS_H
#define OPTIMIZATIONS_H

/**
 * Precompute lookup tables of FFT constants.
 * Comment to disable.
 * If FFT_LUT_1024_64 is defined, this has no effect.
 */
#define COMPUTE_FFT_LUT

/**
 * Use hardcoded lookup tables of GF(1024).
 * Comment to disable.
 * If GF_M != 10, this has no effect.
 */
#define GF_LUT_1024

/**
 * Use an optimal basis for FFT betas.
 * Comment to use the basis a^{GF_M-1}, a^{GF_M-2}, ..., 1.
 * If GF_M != 10 or FFT_LUT_1024_64 is defined, this has no effect.
 */
#define FFT_BETAS_1024

/**
 * Use hardcoded lookup tables of FFT on GF(1024) with FFT_PARAM = 6.
 * Comment to disable.
 * If defined, overrides COMPUTES_FFT_LUT and FFT_BETAS_1024.
 * If GF_M != 10 or FFT_PARAM != 6, this has no effect.
 */
#define FFT_LUT_1024_64

/**
 * Use unrolled FFT on GF(1024) with FFT_PARAM = 6.
 * Comment to disable.
 * Only takes effect if GF_M = 10, FFT_PARAM = 6 and
 * either FFT_LUT_1024_64 or COMPUTE_FFT_LUT is defined
 */
#define FFT_UNROLL_1024_64

/**
 * Use unrolled FFT transpose on GF(1024) with FFT_PARAM = 6.
 * Comment to disable.
 * Only takes effect if GF_M = 10, FFT_PARAM = 6 and
 * either FFT_LUT_1024_64 or COMPUTE_FFT_LUT is defined
 */
#define FFT_T_UNROLL_1024_64

/**
 * Use optimized inversion of GF(1024).
 * Commment to disable.
 * Only takes effect if gf_pclmul.c is used and GF_M = 10.
 */
#define GF_ARITHMETIC_1024

#if GF_M != 10
#undef GF_LUT_1024
#undef FFT_BETAS_1024
#undef FFT_LUT_1024_64
#undef FFT_UNROLL_1024_64
#undef FFT_T_UNROLL_1024_64
#undef GF_ARITHMETIC_1024
#endif

#if FFT_PARAM != 6
#undef FFT_LUT_1024_64
#undef FFT_UNROLL_1024_64
#undef FFT_T_UNROLL_1024_64
#endif

#ifdef FFT_LUT_1024_64
#undef COMPUTE_FFT_LUT
#endif

#if (!defined FFT_LUT_1024_64) && (!defined COMPUTE_FFT_LUT)
#undef FFT_UNROLL_1024_64
#undef FFT_T_UNROLL_1024_64
#endif

#endif // OPTIMIZATIONS_H
