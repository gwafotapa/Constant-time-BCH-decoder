/**
 * @file fft.h
 * Header file of fft.c
 */

#ifndef FFT_H
#define FFT_H

#include "optimizations.h"

#include <stddef.h>
#include <stdint.h>

void
compute_fft_lut();

void
compute_fft_betas(uint16_t* betas);

void
compute_subset_sums(uint16_t    *subset_sums,
                    uint16_t    *set,
                    size_t       set_size);

#ifdef FFT_T_UNROLL_1024_64

void
radix_t_1024_4(uint16_t         *f,
               const uint16_t   *f0,
               const uint16_t   *f1);

void
radix_t_1024_8(uint16_t         *f,
               const uint16_t   *f0,
               const uint16_t   *f1);

void
radix_t_1024_16(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1);

void
radix_t_1024_32(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1);

void
radix_t_1024_64(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1);

void
radix_t_1024_128(uint16_t       *f,
                 const uint16_t *f0,
                 const uint16_t *f1);

void
fft_t_1024_2(uint16_t           *f,
             const uint16_t     *w);

void
fft_t_1024_4(uint16_t           *f,
             const uint16_t     *w,
             size_t              f_coeffs);

void
fft_t_1024_8(uint16_t           *f,
             const uint16_t     *w,
             size_t              f_coeffs);

void
fft_t_1024_16(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs);

void
fft_t_1024_32(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs);

void
fft_t_1024_64(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs);

void
fft_t_1024_128(uint16_t         *f,
               const uint16_t   *w,
               size_t            f_coeffs);

#else // !defined FFT_T_UNROLL_1024_64

void
radix_t(uint16_t        *f,
        const uint16_t  *f0,
        const uint16_t  *f1,
        unsigned int     m_f);

#if (defined COMPUTE_FFT_LUT) || (defined FFT_LUT_1024_64)

void
fft_t_rec(uint16_t              *f,
          const uint16_t        *w,
          size_t                 f_coeffs,
          unsigned int           m,
          unsigned int           m_f,
          const uint16_t        *gammas_sums,
          const uint16_t        *betas_pows,
          const uint16_t        *fft_t_final_betas_sums);

void
fft_t(uint16_t          *f,
      const uint16_t    *w,
      size_t             f_coeffs);

#else // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)

void
fft_t_rec(uint16_t              *f,
          const uint16_t        *w,
          size_t                 f_coeffs,
          unsigned int           m,
          unsigned int           m_f,
          const uint16_t        *betas);

void
fft_t(uint16_t          *f,
      const uint16_t    *w,
      size_t             f_coeffs);

#endif // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
#endif // !defined FFT_T_UNROLL_1024_64

#ifdef FFT_UNROLL_1024_64

void
radix_1024_4(uint16_t           *f0,
             uint16_t           *f1,
             const uint16_t     *f);

void
radix_1024_8(uint16_t           *f0,
             uint16_t           *f1,
             const uint16_t     *f);

void
radix_1024_16(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f);

void
radix_1024_32(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f);

void
radix_1024_64(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f);

void
fft_1024_2 (uint16_t            *w,
            const uint16_t      *f);

void
fft_1024_4(uint16_t     *w,
           uint16_t     *f,
           size_t        f_coeffs);

void
fft_1024_8(uint16_t     *w,
           uint16_t     *f,
           size_t        f_coeffs);

void
fft_1024_16(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs);

void
fft_1024_32(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs);

void
fft_1024_64(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs);

#else // !defined FFT_UNROLL_1024_64

void
radix(uint16_t          *f0,
      uint16_t          *f1,
      const uint16_t    *f,
      unsigned int       m_f);

#if (defined COMPUTE_FFT_LUT) || (defined FFT_LUT_1024_64)

void
fft_rec(uint16_t        *w,
        uint16_t        *f,
        size_t           f_coeffs,
        unsigned int     m,
        unsigned int     m_f,
        const uint16_t  *gammas_sums,
        const uint16_t  *betas_pows,
        const uint16_t  *fft_final_betas);

void
fft(uint16_t    *w,
    uint16_t    *f,
    size_t       f_coeffs);

#else // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)

void
fft_rec(uint16_t        *w,
        uint16_t        *f,
        size_t           f_coeffs,
        unsigned int     m,
        unsigned int     m_f,
        const uint16_t  *betas);

void
fft(uint16_t    *w,
    uint16_t    *f,
    size_t       f_coeffs);

#endif // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
#endif // !defined FFT_UNROLL_1024_64

void
fft_t_preprocess_bch_codeword(uint16_t          *w,
                              const uint8_t     *rcv);

void
fft_retrieve_bch_error_poly(uint8_t             *error_poly,
                            const uint16_t      *w);

#endif // FFT_H
