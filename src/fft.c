/**
 * @file fft.c
 * Implementation of the additive FFT and its transpose.
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.clemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 */

#include "fft.h"

#include "gf.h"

#include <string.h>

#ifdef FFT_LUT_1024_64 // look up tables for FFT constants on the field GF(2^10)
#include "fft_lut_1024_64.h"
#elif defined COMPUTE_FFT_LUT
/**
 * Constants for the additive FFT and its transpose.
 */
static uint16_t gammas_sums[(1 << GF_M-(FFT_T_PARAM-1)) * ((1 << FFT_T_PARAM-1)-1)];

/**
 * Constants for the additive FFT and its transpose.
 */
#if FFT_PARAM == 1
static const uint16_t *betas_pows = NULL;
#else
static uint16_t betas_pows[4 * (1 << FFT_T_PARAM-2) - 4];
#endif

/**
 * Constants for the additive FFT transpose.
 */
static uint16_t fft_t_final_betas_sums[1 << GF_M-(FFT_T_PARAM-1)]; // used in fft_t last step only

/**
 * Constants for the additive FFT.
 */
static uint16_t fft_final_betas[GF_M - (FFT_PARAM-1)]; // used in fft last step only

/**
 * Computes all constants involved in the additive FFT and its transpose.
 */
void
compute_fft_lut()
{
        uint16_t gammas[GF_M-1];
        compute_fft_betas(gammas);

        uint16_t *gs = gammas_sums;
        uint16_t *fb = fft_final_betas;
        uint16_t *fbs = fft_t_final_betas_sums;
        
#if FFT_PARAM == 1
        memcpy(fb, gammas, 2*(GF_M-1));
        fb[GF_M-1] = 1;

        compute_subset_sums(gs, gammas, GF_M-1);
        
        fbs[0] = 0;
        for (size_t k = 0; k < GF_M-FFT_PARAM; ++k) {
                gammas[k] = gf_square(gammas[k]) ^ gammas[k];
                for (size_t l = 0; l < (1 << k); ++l)
                        fbs[(1 << k) + l] = fbs[l] ^ gammas[k];
        }
#else
        uint16_t *bp = betas_pows;
  
        // Compute the first gammas subset sums.
        // Since the GF_M-th beta is 1, there's no twisting on the first step,
        // i.e. gammas are betas.
        compute_subset_sums(gs, gammas, GF_M-1);

        gs += 1 << (GF_M-1);
  
        // First call is done. Compute constants for the recursive calls
        for (size_t j = 1; j < FFT_PARAM; ++j) {
    
                // Compute deltas (which are the next betas)
                for (size_t k = 0; k < GF_M-j; ++k)
                        gammas[k] = gf_square(gammas[k]) ^ gammas[k];

                // Compute last beta powers
                uint16_t beta_m = gammas[GF_M-j-1];
                bp[0] = 1;
                for (size_t i = 1; i < (1 << FFT_T_PARAM-j); ++i)
                        bp[i] = gf_mul(bp[i-1], beta_m);

                if (j == FFT_PARAM-1) // Compute betas from the last fft call
                        memcpy(fb, gammas, 2 * (GF_M - (FFT_PARAM-1)));
    
                // Compute gammas
                for (size_t k = 0; k < GF_M-j-1; ++k)
                        gammas[k] = gf_mul(gammas[k], gf_inverse(beta_m));

                // Compute gammas subset sums
                compute_subset_sums(gs, gammas, GF_M-j-1);

                bp += 1 << (FFT_T_PARAM-j);
                gs += 1 << (GF_M-1-j);
        }

        // Compute betas sums from the last fft transpose call
        fbs[0] = 0;
        for (size_t k = 0; k < GF_M-FFT_PARAM; ++k) {
                gammas[k] = gf_square(gammas[k]) ^ gammas[k];
                for (size_t l = 0; l < (1 << k); ++l)
                        fbs[(1 << k) + l] = fbs[l] ^ gammas[k];
        }
#endif  
}

#endif // COMPUTE_FFT_LUT

/**
 * Computes the basis of betas (omitting 1) used in the additive FFT and its transpose.
 * @param[out] betas Array of size GF_M-1
 */
void
compute_fft_betas(uint16_t *betas)
{
        for (size_t i = 0; i < GF_M-1; ++i)
                betas[i] = 1 << GF_M-1-i;

        // Optimisation suggested by Bernstein et al. to avoid twisting.
        // Implemented for GF(2^10) only.
#ifdef FFT_BETAS_1024
        betas[7] = 2;
        betas[8] = 237;
#endif  
}

/**
 * Computes the subset sums of the given set.
 * The array subset_sums is such that its ith element is
 * the subset sum of the set elements given by the binary form of i.
 * @param[out] subset_sums Array of size 2^set_size receiving the subset sums
 * @param[in] set Array of set_size elements
 * @param[in] set_size Size of the array set
 */
void
compute_subset_sums(uint16_t    *subset_sums,
                    uint16_t    *set,
                    size_t       set_size)
{
        subset_sums[0] = 0;
        for (size_t i = 0; i < set_size; ++i)
                for (size_t j = 0; j < 1<<i; ++j)
                        subset_sums[(1 << i) + j] = set[i] ^ subset_sums[j];
}

#ifdef FFT_T_UNROLL_1024_64

/**
 * Last recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 4 elements
 * @param[in] f0 Array of 2 elements
 * @param[in] f1 Array of 2 elements
 */
void
radix_t_1024_4(uint16_t         *f,
               const uint16_t   *f0,
               const uint16_t   *f1)
{
        f[0] = f0[0];
        f[1] = f1[0];
        f[2] = f0[1] ^ f1[0];
        f[3] = f[2] ^ f1[1];
}

/**
 * Fourth recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 8 elements
 * @param[in] f0 Array of 4 elements
 * @param[in] f1 Array of 4 elements
 */
void
radix_t_1024_8(uint16_t         *f,
               const uint16_t   *f0,
               const uint16_t   *f1)
{
        f[0] = f0[0];
        f[1] = f1[0];
        f[2] = f0[1] ^ f1[0];
        f[3] = f[2] ^ f1[1];
        f[4] = f[2] ^ f0[2];
        f[5] = f[3] ^ f1[2];
        f[6] = f[4] ^ f0[3] ^ f1[2];
        f[7] = f[3] ^ f0[3] ^ f1[3];
}

/**
 * Third recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 16 elements
 * @param[in] f0 Array of 8 elements
 * @param[in] f1 Array of 8 elements
 */
void
radix_t_1024_16(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1)
{
        f[0] = f0[0];
        f[1] = f1[0];
        f[2] = f0[1] ^ f1[0];
        f[3] = f[2] ^ f1[1];
        f[4] = f[2] ^ f0[2];
        f[5] = f[3] ^ f1[2];
        f[6] = f[4] ^ f0[3] ^ f1[2];
        f[7] = f[3] ^ f0[3] ^ f1[3];
        f[8] = f[4] ^ f0[4];
        f[9] = f[5] ^ f1[4];
        f[10] = f[6] ^ f0[5] ^ f1[4];
        f[11] = f[7] ^ f0[5] ^ f1[4] ^ f1[5];
        f[12] = f[8] ^ f0[5] ^ f0[6] ^ f1[4];
        f[13] = f[7] ^ f[9] ^ f[11] ^ f1[6];
        f[14] = f[6] ^ f0[6] ^ f0[7] ^ f1[6];
        f[15] = f[7] ^ f0[7] ^ f1[7];
}

/**
 * Second recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 32 elements
 * @param[in] f0 Array of 16 elements
 * @param[in] f1 Array of 16 elements
 */
void
radix_t_1024_32(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1)
{
        uint16_t Q0[8];
        uint16_t Q1[8];
        uint16_t R0[8];
        uint16_t R1[8];
        
        memcpy(Q0, f0 + 8, 2*8);
        memcpy(Q1, f1 + 8, 2*8);
        memcpy(R0, f0, 2*8);
        memcpy(R1, f1, 2*8);
    
        uint16_t Q[2*8];
        uint16_t R[2*8];
        
        radix_t_1024_16(Q, Q0, Q1);
        radix_t_1024_16(R, R0, R1);

        memcpy(f, R, 4*8);
        memcpy(f + 3*8, Q + 8, 2*8);

        for (size_t i = 0; i < 8; ++i) {
                f[2*8 + i] = Q[i] ^ R[8 + i];
                f[3*8 + i] ^= f[2*8 + i];
        }
}

/**
 * First recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 64 elements
 * @param[in] f0 Array of 32 elements
 * @param[in] f1 Array of 32 elements
 */
void
radix_t_1024_64(uint16_t        *f,
                const uint16_t  *f0,
                const uint16_t  *f1)
{
        uint16_t Q0[16];
        uint16_t Q1[16];
        uint16_t R0[16];
        uint16_t R1[16];
        
        memcpy(Q0, f0 + 16, 2*16);
        memcpy(Q1, f1 + 16, 2*16);
        memcpy(R0, f0, 2*16);
        memcpy(R1, f1, 2*16);
    
        uint16_t Q[2*16];
        uint16_t R[2*16];
        
        radix_t_1024_32(Q, Q0, Q1);
        radix_t_1024_32(R, R0, R1);

        memcpy(f, R, 4*16);
        memcpy(f + 3*16, Q + 16, 2*16);

        for (size_t i = 0; i < 16; ++i) {
                f[2*16 + i] = Q[i] ^ R[16 + i];
                f[3*16 + i] ^= f[2*16 + i];
        }
}

/**
 * Transpose of the radix conversion on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 128 elements
 * @param[in] f0 Array of 64 elements
 * @param[in] f1 Array of 64 elements
 */
void
radix_t_1024_128(uint16_t       *f,
                 const uint16_t *f0,
                 const uint16_t *f1)
{
        uint16_t Q0[32];
        uint16_t Q1[32];
        uint16_t R0[32];
        uint16_t R1[32];
        
        memcpy(Q0, f0 + 32, 2*32);
        memcpy(Q1, f1 + 32, 2*32);
        memcpy(R0, f0, 2*32);
        memcpy(R1, f1, 2*32);
    
        uint16_t Q[2*32];
        uint16_t R[2*32];
        
        radix_t_1024_64(Q, Q0, Q1);
        radix_t_1024_64(R, R0, R1);

        memcpy(f, R, 4*32);
        memcpy(f + 3*32, Q + 32, 2*32);

        for (size_t i = 0; i < 32; ++i) {
                f[2*32 + i] = Q[i] ^ R[32 + i];
                f[3*32 + i] ^= f[2*32 + i];
        }
}

/**
 * Last recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 2 elements
 * @param[in] w Array of 16 elements
 */
void
fft_t_1024_2(uint16_t           *f,
             const uint16_t     *w)
{
        f[0] = 0;
        for (size_t i = 0; i < 16; ++i)
                f[0] ^= w[i];
        f[1] = 0;
        for (size_t j = 0; j < 4; ++j) {
                for (size_t k = 0; k < (1 << j); ++k) {
                        size_t index = (1 << j) + k;
                        f[1] ^= gf_mul(fft_t_final_betas_sums[index], w[index]);
                }
        }
}

/**
 * Fifth recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 4 elements
 * @param[in] w Array of 32 elements
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_4(uint16_t           *f,
             const uint16_t     *w,
             size_t              f_coeffs)
{
        uint16_t u[16];  
        uint16_t f0[2];
        uint16_t f1[2];

        if (f_coeffs <= 3) {
                f1[1] = 0;
                u[0] = w[0] ^ w[16];
                f1[0] = w[16];
                for (size_t i = 1; i < 16; ++i) {
                        u[i] = w[i] ^ w[16+i];
                        f1[0] ^= gf_mul(gammas_sums[512+256+128+64+32+i], u[i]) ^ w[16+i];
                }
                fft_t_1024_2(f0, u);
        }
        else {
                uint16_t v[16];
  
                u[0] = w[0] ^ w[16];
                v[0] = w[16];

                for (size_t i = 1; i < 16; ++i) {
                        u[i] = w[i] ^ w[16+i];
                        v[i] = gf_mul(gammas_sums[512+256+128+64+32+i], u[i]) ^ w[16+i];
                }

                fft_t_1024_2(f0, u);
                fft_t_1024_2(f1, v);
        }

        radix_t_1024_4(f, f0, f1);

        for (size_t i = 1; i < 4; ++i)
                f[i] = gf_mul(betas_pows[64+32+16+8+i], f[i]);
}

/**
 * Fourth recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 8 elements
 * @param[in] w Array of 64 elements
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_8(uint16_t           *f,
             const uint16_t     *w,
             size_t              f_coeffs)
{
        uint16_t u[32];
        uint16_t v[32];
  
        u[0] = w[0] ^ w[32];
        v[0] = w[32];

        for (size_t i = 1; i < 32; ++i) {
                u[i] = w[i] ^ w[32+i];
                v[i] = gf_mul(gammas_sums[512+256+128+64+i], u[i]) ^ w[32+i];
        }

        uint16_t f0[4];
        uint16_t f1[4];

        fft_t_1024_4(f0, u, (f_coeffs+1)/2);
        fft_t_1024_4(f1, v, f_coeffs/2);

        radix_t_1024_8(f, f0, f1);

        for (size_t i = 1; i < 8; ++i)
                f[i] = gf_mul(betas_pows[64+32+16+i], f[i]);
}

/**
 * Third recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 16 elements
 * @param[in] w Array of 128 elements
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_16(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs)
{
        uint16_t u[64];
        uint16_t v[64];
  
        u[0] = w[0] ^ w[64];
        v[0] = w[64];

        for (size_t i = 1; i < 64; ++i) {
                u[i] = w[i] ^ w[64+i];
                v[i] = gf_mul(gammas_sums[512+256+128+i], u[i]) ^ w[64+i];
        }

        uint16_t f0[8];
        uint16_t f1[8];

        fft_t_1024_8(f0, u, (f_coeffs+1)/2);
        fft_t_1024_8(f1, v, f_coeffs/2);

        radix_t_1024_16(f, f0, f1);

        for (size_t i = 1; i < 16; ++i)
                f[i] = gf_mul(betas_pows[64+32+i], f[i]);
}

/**
 * Second recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 32 elements
 * @param[in] w Array of 256 elements
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_32(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs)
{
        uint16_t u[128];
        uint16_t v[128];
  
        u[0] = w[0] ^ w[128];
        v[0] = w[128];

        for (size_t i = 1; i < 128; ++i) {
                u[i] = w[i] ^ w[128+i];
                v[i] = gf_mul(gammas_sums[512+256+i], u[i]) ^ w[128+i];
        }

        uint16_t f0[16];
        uint16_t f1[16];

        fft_t_1024_16(f0, u, (f_coeffs+1)/2);
        fft_t_1024_16(f1, v, f_coeffs/2);

        radix_t_1024_32(f, f0, f1);

        for (size_t i = 1; i < 32; ++i)
                f[i] = gf_mul(betas_pows[64+i], f[i]);
}

/**
 * First recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 64 elements
 * @param[in] w Array of 512 elements
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_64(uint16_t          *f,
              const uint16_t    *w,
              size_t             f_coeffs)
{
        uint16_t u[256];
        uint16_t v[256];
  
        u[0] = w[0] ^ w[256];
        v[0] = w[256];

        for (size_t i = 1; i < 256; ++i) {
                u[i] = w[i] ^ w[256+i];
                v[i] = gf_mul(gammas_sums[512+i], u[i]) ^ w[256+i];
        }

        uint16_t f0[32];
        uint16_t f1[32];
        
        fft_t_1024_32(f0, u, (f_coeffs+1)/2);
        fft_t_1024_32(f1, v, f_coeffs/2);

        radix_t_1024_64(f, f0, f1);

#ifndef FFT_BETAS_1024
        for (size_t i = 1; i < 64; ++i)
                f[i] = gf_mul(betas_pows[i], f[i]);
#endif
}

/**
 * Additive FFT Transpose for a family of 1024 elements of GF(2^10).
 * Computes the first f_coeffs syndromes of family w.
 * @param[out] f Array of 128 elements receiving the syndromes
 * @param[in] w Array of 1024 elements of GF(2^10) storing the family
 * @param[in] f_coeffs Length of f
 */
void
fft_t_1024_128(uint16_t         *f,
               const uint16_t   *w,
               size_t            f_coeffs)
{ 
        uint16_t u[512];
        uint16_t v[512];

        u[0] = w[0] ^ w[512];
        v[0] = w[512];

        for (size_t i = 1; i < 512; ++i) {
                u[i] = w[i] ^ w[512+i];
                v[i] = gf_mul(gammas_sums[i], u[i]) ^ w[512+i];
        }

        uint16_t f0[64];
        uint16_t f1[64];
        
        fft_t_1024_64(f0, u, (f_coeffs+1)/ 2);
        fft_t_1024_64(f1, v, f_coeffs/2);

        radix_t_1024_128(f, f0, f1);
}

#else // !defined FFT_T_UNROLL_1024_64

/**
 * Transpose of the linear radix conversion.
 * This is a direct transposition of the radix(...) function
 * implemented following the process of transposing a linear function as exposed by Bernstein, Chou and Schwabe here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 * @param[out] f Array of size a power of 2
 * @param[in] f0 Array half the size of f
 * @param[in] f1 Array half the size of f
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to the number of coefficients of f
 */
void
radix_t(uint16_t        *f,
        const uint16_t  *f0,
        const uint16_t  *f1,
        unsigned int     m_f)
{
        switch (m_f) {
        case 4:
                f[0] = f0[0];
                f[1] = f1[0];
                f[2] = f0[1] ^ f1[0];
                f[3] = f[2] ^ f1[1];
                f[4] = f[2] ^ f0[2];
                f[5] = f[3] ^ f1[2];
                f[6] = f[4] ^ f0[3] ^ f1[2];
                f[7] = f[3] ^ f0[3] ^ f1[3];
                f[8] = f[4] ^ f0[4];
                f[9] = f[5] ^ f1[4];
                f[10] = f[6] ^ f0[5] ^ f1[4];
                f[11] = f[7] ^ f0[5] ^ f1[4] ^ f1[5];
                f[12] = f[8] ^ f0[5] ^ f0[6] ^ f1[4];
                f[13] = f[7] ^ f[9] ^ f[11] ^ f1[6];
                f[14] = f[6] ^ f0[6] ^ f0[7] ^ f1[6];
                f[15] = f[7] ^ f0[7] ^ f1[7];
                return;

        case 3:
                f[0] = f0[0];
                f[1] = f1[0];
                f[2] = f0[1] ^ f1[0];
                f[3] = f[2] ^ f1[1];
                f[4] = f[2] ^ f0[2];
                f[5] = f[3] ^ f1[2];
                f[6] = f[4] ^ f0[3] ^ f1[2];
                f[7] = f[3] ^ f0[3] ^ f1[3];
                return;
  
        case 2:
                f[0] = f0[0];
                f[1] = f1[0];
                f[2] = f0[1] ^ f1[0];
                f[3] = f[2] ^ f1[1];
                return;

        case 1:
                f[0] = f0[0];
                f[1] = f1[0];
                return;
    
        default:;
                size_t n = 1 << (m_f-2);

                uint16_t Q0[n];
                uint16_t Q1[n];
                uint16_t R0[n];
                uint16_t R1[n];         
                
                memcpy(Q0, f0 + n, 2*n);
                memcpy(Q1, f1 + n, 2*n);
                memcpy(R0, f0, 2*n);
                memcpy(R1, f1, 2*n);
    
                uint16_t Q[2*n];
                uint16_t R[2*n];
                
                radix_t (Q, Q0, Q1, m_f-1);
                radix_t (R, R0, R1, m_f-1);

                memcpy(f, R, 4*n);
                memcpy(f + 2*n, R + n, 2*n);
                memcpy(f + 3*n, Q + n, 2*n);

                for (size_t i = 0; i < n; ++i) {
                        f[2*n + i] ^= Q[i];
                        f[3*n + i] ^= f[2*n + i];
                }
        }
}

#if (defined COMPUTE_FFT_LUT) || (defined FFT_LUT_1024_64)

/**
 * Recursively computes syndromes of family w.
 * This function is a subroutine of the function fft_t(...).
 * @param[out] f Array receiving the syndromes
 * @param[in] w Array storing the family
 * @param[in] f_coeffs Length of syndromes vector
 * @param[in] m 2^m is the smallest power of 2 greater or equal to the length of family w
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to the length of f
 * @param[in] gammas_sums Logarithms of subset sums of gammas needed at each step
 * @param[in] betas_pows Logarithms of powers of beta^m needed at each recursive step
 * @param[in] fft_t_final_betas_sums Logarithms of betas subset sums at the last step
 */
void
fft_t_rec(uint16_t              *f,
          const uint16_t        *w,
          size_t                 f_coeffs,
          unsigned int           m,
          unsigned int           m_f,
          const uint16_t        *gammas_sums,
          const uint16_t        *betas_pows,
          const uint16_t        *fft_t_final_betas_sums)
{
        // Step 1
        if (m_f == 1) {
                f[0] = 0;
                for (size_t i = 0; i < (1 << m); ++i)
                        f[0] ^= w[i];
                f[1] = 0;
                for (size_t j = 0; j < m; ++j) {
                        for (size_t k = 0; k < (1 << j); ++k) {
                                size_t index = (1 << j) + k;
                                f[1] ^= gf_mul(fft_t_final_betas_sums[index], w[index]);
                        }
                }

                return;
        }

        /* Step 6: Compute u and v from w (aka w)
         *
         * We had:
         * w[i] = u[i] + G[i].v[i]
         * w[k+i] = w[i] + v[i] = u[i] + (G[i]+1).v[i]
         * Transpose:
         * u[i] = w[i] + w[k+i]
         * v[i] = G[i].w[i] + (G[i]+1).w[k+i] = G[i].u[i] + w[k+i] */
        size_t k = 1 << m-1;    
        uint16_t u[k];
        uint16_t f0[1 << m_f-1];
        uint16_t f1[1 << m_f-1];

        if (f_coeffs <= 3) { // 3-coefficient polynomial f case
                // Step 5: Compute f0 from u and f1 from v
                f1[1] = 0;
                u[0] = w[0] ^ w[k];
                f1[0] = w[k];
                for (size_t i = 1; i < k; ++i) {
                        u[i] = w[i] ^ w[k+i];
                        f1[0] ^= gf_mul(gammas_sums[i], u[i]) ^ w[k+i];
                }
                fft_t_rec(f0, u, (f_coeffs+1)/2, m-1, m_f-1, gammas_sums + k, betas_pows + (1<<m_f), fft_t_final_betas_sums);
        }
        else {
                uint16_t v[k];
    
                u[0] = w[0] ^ w[k];
                v[0] = w[k];

                for (size_t i = 1; i < k; ++i) {
                        u[i] = w[i] ^ w[k+i];
                        v[i] = gf_mul(gammas_sums[i], u[i]) ^ w[k+i];
                }

                // Step 5: Compute f0 from u and f1 from v
                fft_t_rec(f0, u, (f_coeffs+1)/2, m-1, m_f-1, gammas_sums + k, betas_pows + (1<<m_f), fft_t_final_betas_sums);
                fft_t_rec(f1, v, f_coeffs/2, m-1, m_f-1, gammas_sums + k, betas_pows + (1<<m_f), fft_t_final_betas_sums);
        }

        // Step 4: Constants gammas and deltas have been precomputed
  
        // Step 3: Compute g from g0 and g1
        radix_t(f, f0, f1, m_f);
  
        // Step 2: compute f from g
        if (betas_pows[1] != 1) {
                for (size_t i = 1; i < (1 << m_f); ++i)
                        f[i] = gf_mul(betas_pows[i], f[i]);
        }
}

/**
 * Computes the syndromes f of the family w.
 * Since the syndromes linear map is the transpose of multipoint evaluation,
 * it uses exactly the same constants, either hardcoded or precomputed by compute_fft_lut(...). <br>
 * This follows directives from Bernstein, Chou and Schwabe given here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 * @param[out] f Array of 2*(FFT_T_PARAM) elements
 * @param[in] w Array of GF_MUL_ORDER+1 elements
 * @param[in] f_coeffs Length of syndromes vector f
 */
void
fft_t(uint16_t          *f,
      const uint16_t    *w,
      size_t             f_coeffs)
{
        // Transposed from Gao and Mateer algorithm
  
        /* Step 6: Compute u and v from w (aka w)
         *
         * We had:
         * w[i] = u[i] + G[i].v[i]
         * w[k+i] = w[i] + v[i] = u[i] + (G[i]+1).v[i]
         * Transpose:
         * u[i] = w[i] + w[k+i]
         * v[i] = G[i].w[i] + (G[i]+1).w[k+i] = G[i].u[i] + w[k+i] */
        size_t k = 1 << GF_M-1;  
        uint16_t u[k];
        uint16_t v[k];

        u[0] = w[0] ^ w[k];
        v[0] = w[k];

        for (size_t i = 1; i < k; ++i) {
                u[i] = w[i] ^ w[k+i];
                v[i] = gf_mul(gammas_sums[i], u[i]) ^ w[k+i];
        }

        // Step 5: Compute f0 from u and f1 from v
        uint16_t f0[1 << FFT_T_PARAM-1];
        uint16_t f1[1 << FFT_T_PARAM-1];
        
        fft_t_rec(f0, u, (f_coeffs+1)/2, GF_M-1, FFT_T_PARAM-1, gammas_sums + k, betas_pows, fft_t_final_betas_sums);
        fft_t_rec(f1, v, f_coeffs/2, GF_M-1, FFT_T_PARAM-1, gammas_sums + k, betas_pows, fft_t_final_betas_sums);

        // Step 4: Constants gammas and deltas have been precomputed
  
        // Step 3: Compute g from g0 and g1
        radix_t(f, f0, f1, FFT_T_PARAM);

        // Step 2: beta_m = 1 so f = g
}

#else // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)

/**
 * Recursively computes syndromes of family w.
 * This function is a subroutine of the function fft_t(...).
 * @param[out] f Array receiving the syndromes
 * @param[in] w Array storing the family
 * @param[in] f_coeffs Length of syndromes vector
 * @param[in] m 2^m is the smallest power of 2 greater or equal to the length of family w
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to the length of f
 * @param[in] betas FFT constants
 */
void
fft_t_rec(uint16_t              *f,
          const uint16_t        *w,
          size_t                 f_coeffs,
          unsigned int           m,
          unsigned int           m_f,
          const uint16_t        *betas)
{
        // Step 1
        if (m_f == 1) {
                f[0] = 0;
                for (size_t i = 0; i < (1 << m); ++i)
                        f[0] ^= w[i];
                f[1] = 0;
                uint16_t betas_sums[1 << m];
                betas_sums[0] = 0;
                for (size_t j = 0; j < m; ++j) {
                        for (size_t k = 0; k < (1 << j); ++k) {
                                size_t index = (1 << j) + k;
                                betas_sums[index] = betas_sums[k] ^ betas[j];
                                f[1] ^= gf_mul(betas_sums[index], w[index]);
                        }
                }
                return;
        }

        size_t k = 1 << m-1;

        // Compute gammas and deltas
        uint16_t gammas[m-1];
        uint16_t deltas[m-1];

        for (size_t i = 0; i < m-1; ++i) {
                gammas[i] = gf_mul(betas[i], gf_inverse(betas[m-1]));
                deltas[i] = gf_square(gammas[i]) ^ gammas[i];
        }

        // Compute gammas subset sums
        uint16_t gammas_sums[k];
        compute_subset_sums(gammas_sums, gammas, m-1);
        
        /* Step 6: Compute u and v from w (aka w)
         *
         * We had:
         * w[i] = u[i] + G[i].v[i]
         * w[k+i] = w[i] + v[i] = u[i] + (G[i]+1).v[i]
         * Transpose:
         * u[i] = w[i] + w[k+i]
         * v[i] = G[i].w[i] + (G[i]+1).w[k+i] = G[i].u[i] + w[k+i] */
        uint16_t u[k];
        uint16_t f0[1 << m_f-1];
        uint16_t f1[1 << m_f-1];

        if (f_coeffs <= 3) { // 3-coefficient polynomial f case
                // Step 5: Compute f0 from u and f1 from v
                f1[1] = 0;
                u[0] = w[0] ^ w[k];
                f1[0] = w[k];
                for (size_t i = 1; i < k; ++i) {
                        u[i] = w[i] ^ w[k+i];
                        f1[0] ^= gf_mul(gammas_sums[i], u[i]) ^ w[k+i];
                }
                fft_t_rec(f0, u, (f_coeffs+1)/2, m-1, m_f-1, deltas);
        }
        else {
                uint16_t v[k];
    
                u[0] = w[0] ^ w[k];
                v[0] = w[k];

                for (size_t i = 1; i < k; ++i) {
                        u[i] = w[i] ^ w[k+i];
                        v[i] = gf_mul(gammas_sums[i], u[i]) ^ w[k+i];
                }
    
                // Step 5: Compute f0 from u and f1 from v
                fft_t_rec(f0, u, (f_coeffs+1)/2, m-1, m_f-1, deltas);
                fft_t_rec(f1, v, f_coeffs/2, m-1, m_f-1, deltas);
        }
                     
        // Step 3: Compute g from g0 and g1
        radix_t(f, f0, f1, m_f);
  
        // Step 2: compute f from g
        if (betas[m-1] != 1) {
                uint16_t beta_m_pow = 1;
                for (size_t i = 1; i < (1 << m_f); ++i) {
                        beta_m_pow = gf_mul(beta_m_pow, betas[m-1]);
                        f[i] = gf_mul(beta_m_pow, f[i]);
                }
        }
}

/**
 * Computes the syndromes f of the family w.
 * Since the syndromes linear map is the transpose of multipoint evaluation,
 * it uses exactly the same constants, either hardcoded or precomputed by compute_fft_lut(...). <br>
 * This follows directives from Bernstein, Chou and Schwabe given here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 * @param[out] f Array of size 2*(FFT_T_PARAM) elements receiving the syndromes
 * @param[in] w Array of GF_MUL_ORDER+1 elements
 * @param[in] f_coeffs Length of syndromes vector f
 */
void
fft_t(uint16_t          *f,
      const uint16_t    *w,
      size_t             f_coeffs)
{
        // Transposed from Gao and Mateer algorithm

        uint16_t betas[GF_M-1];
        uint16_t betas_sums[1 << GF_M-1];
        
        compute_fft_betas(betas);
        compute_subset_sums(betas_sums, betas, GF_M-1);
        
        /* Step 6: Compute u and v from w (aka w)
         *
         * We had:
         * w[i] = u[i] + G[i].v[i]
         * w[k+i] = w[i] + v[i] = u[i] + (G[i]+1).v[i]
         * Transpose:
         * u[i] = w[i] + w[k+i]
         * v[i] = G[i].w[i] + (G[i]+1).w[k+i] = G[i].u[i] + w[k+i] */
        size_t k = 1 << GF_M-1;
        uint16_t u[k];
        uint16_t v[k];
        
        u[0] = w[0] ^ w[k];
        v[0] = w[k];
        for (size_t i = 1; i < k; ++i) {
                u[i] = w[i] ^ w[k+i];
                v[i] = gf_mul(betas_sums[i], u[i]) ^ w[k+i];
        }

        // Compute deltas
        uint16_t deltas[GF_M-1];
        for (size_t i = 0; i < GF_M-1; ++i)
                deltas[i] = gf_square(betas[i]) ^ betas[i];
  
        // Step 5: Compute f0 from u and f1 from v
        uint16_t f0[1 << FFT_T_PARAM-1];
        uint16_t f1[1 << FFT_T_PARAM-1];
        
        fft_t_rec(f0, u, (f_coeffs+1)/2, GF_M-1, FFT_T_PARAM-1, deltas);
        fft_t_rec(f1, v, f_coeffs/2, GF_M-1, FFT_T_PARAM-1, deltas);

        // Step 3: Compute g from g0 and g1
        radix_t(f, f0, f1, FFT_T_PARAM);

        // Step 2: beta_m = 1 so f = g
}

#endif // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
#endif // !defined FFT_T_UNROLL_1024_64

#ifdef FFT_UNROLL_1024_64

/**
 * Last recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 2 elements
 * @param[out] f1 Array of 2 elements
 * @param[in] f Array of 4 elements
 */
void
radix_1024_4(uint16_t           *f0,
             uint16_t           *f1,
             const uint16_t     *f)
{
        f0[0] = f[0];
        f0[1] = f[2]^f[3];
        f1[0] = f[1]^f0[1];
        f1[1] = f[3];
}

/**
 * Third recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 4 elements
 * @param[out] f1 Array of 4 elements
 * @param[in] f Array of 8 elements
 */
void
radix_1024_8(uint16_t           *f0,
             uint16_t           *f1,
             const uint16_t     *f)
{
        f0[0] = f[0];
        f0[2] = f[4]^f[6];
        f0[3] = f[6]^f[7];
        f1[1] = f[3]^f[5]^f[7];
        f1[2] = f[5]^f[6];
        f1[3] = f[7];
        f0[1] = f[2]^f0[2]^f1[1];
        f1[0] = f[1]^f0[1];
}

/**
 * Second recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 8 elements
 * @param[out] f1 Array of 8 elements
 * @param[in] f Array of 16 elements
 */
void
radix_1024_16(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f)
{
        f0[4] = f[8]^f[12];
        f0[6] = f[12]^f[14];
        f0[7] = f[14]^f[15];
        f1[5] = f[11]^f[13];
        f1[6] = f[13]^f[14];
        f1[7] = f[15];
        f0[5] = f[10]^f[12]^f1[5];
        f1[4] = f[9]^f[13]^f0[5];

        f0[0] = f[0];
        f1[3] = f[7]^f[11]^f[15];
        f0[3] = f[6]^f[10]^f[14]^f1[3];
        f0[2] = f[4]^f0[4]^f0[3]^f1[3];
        f1[1] = f[3]^f[5]^f[9]^f[13]^f1[3];
        f1[2] = f[3]^f1[1]^f0[3];
        f0[1] = f[2]^f0[2]^f1[1];
        f1[0] = f[1]^f0[1];
}

/**
 * First recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 16 elements
 * @param[out] f1 Array of 16 elements
 * @param[in] f Array of 32 elements
 */
void
radix_1024_32(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f)
{
        uint16_t Q[2*8];
        uint16_t R[2*8];

        memcpy(Q + 8, f + 3*8, 2*8);
        memcpy(R, f, 4*8);

        for (size_t i = 0; i < 8; ++i) {
                Q[i] = f[2*8 + i] ^ f[3*8 + i];
                R[8 + i] ^= Q[i];
        }

        uint16_t Q0[8];
        uint16_t Q1[8];
        uint16_t R0[8];
        uint16_t R1[8];

        radix_1024_16 (Q0, Q1, Q);
        radix_1024_16 (R0, R1, R);

        memcpy(f0, R0, 2*8);
        memcpy(f0 + 8, Q0, 2*8);
        memcpy(f1, R1, 2*8);
        memcpy(f1 + 8, Q1, 2*8);
}

/**
 * Computes the radix conversion of f,
 * that is f0 and f1 such that f(x) = f0(x^2-x) + x.f1(x^2-x)
 * @param[out] f0 Array of 32 elements
 * @param[out] f1 Array of 32 elements
 * @param[in] f Array of 64 elements
 */
void
radix_1024_64(uint16_t          *f0,
              uint16_t          *f1,
              const uint16_t    *f)
{
        uint16_t Q[2*16];
        uint16_t R[2*16];

        memcpy(Q + 16, f + 3*16, 2*16);
        memcpy(R, f, 4*16);

        for (size_t i = 0; i < 16; ++i) {
                Q[i] = f[2*16 + i] ^ f[3*16 + i];
                R[16 + i] ^= Q[i];
        }

        uint16_t Q0[16];
        uint16_t Q1[16];
        uint16_t R0[16];
        uint16_t R1[16];

        radix_1024_32 (Q0, Q1, Q);
        radix_1024_32 (R0, R1, R);
     
        memcpy(f0, R0, 2*16);
        memcpy(f0 + 16, Q0, 2*16);
        memcpy(f1, R1, 2*16);
        memcpy(f1 + 16, Q1, 2*16);
}

/**
 * Fifth and last recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 32 elements
 * @param[in] f Array of 2 elements
 */
void
fft_1024_2 (uint16_t            *w,
            const uint16_t      *f)
{
        uint16_t tmp[5];
        for (size_t i = 0; i < 5; ++i)
                tmp[i] = gf_mul(fft_final_betas[i], f[1]);
        w[0] = f[0];
        for (size_t j = 0; j < 5; ++j) {
                for (size_t k = 0; k < 1 << j; ++k) {
                        w[(1 << j) + k] = w[k] ^ tmp[j];
                }
        }
}

/**
 * Fourth recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 64 elements
 * @param[in] f Array of 4 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 4)
 */
void
fft_1024_4(uint16_t     *w,
           uint16_t     *f,
           size_t        f_coeffs)
{
        for (size_t i = 1; i < (1 << 2); ++i)
                f[i] = gf_mul(betas_pows[64+32+16+i], f[i]);
  
        uint16_t f0[1 << 2-1];
        uint16_t f1[1 << 2-1];

        radix_1024_4 (f0, f1, f);

        uint16_t u[1 << 6-1];
        uint16_t v[1 << 6-1];
        
        fft_1024_2 (u, f0);
        
        if (f_coeffs <= 3) {
                size_t k = 1 << 6-1;
                w[0] = u[0];
                w[k] = u[0] ^ f1[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+64+i], f1[0]);
                        w[k+i] = w[i] ^ f1[0];
                }
        }
        else {
                fft_1024_2 (v, f1);

                size_t k = 1 << 6-1;
                memcpy(w + k, v, 2*k);
                w[0] = u[0];
                w[k] ^= u[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+64+i], v[i]);
                        w[k+i] ^= w[i];
                }
        }
}

/**
 * Third recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 128 elements
 * @param[in] f Array of 8 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 8)
 */
void
fft_1024_8(uint16_t     *w,
           uint16_t     *f,
           size_t        f_coeffs)
{
        for (size_t i = 1; i < (1 << 3); ++i)
                f[i] = gf_mul(betas_pows[64+32+i], f[i]);

        uint16_t f0[1 << 3-1];
        uint16_t f1[1 << 3-1];

        radix_1024_8 (f0, f1, f);

        uint16_t u[1 << 7-1];
        uint16_t v[1 << 7-1];

        fft_1024_4 (u, f0, (f_coeffs+1)/2);
        fft_1024_4 (v, f1, f_coeffs/2);

        size_t k = 1 << 7-1;
        memcpy(w + k, v, 2*k);
        w[0] = u[0];
        w[k] ^= u[0];
        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+i], v[i]);
                w[k+i] ^= w[i];
        }
}

/**
 * Second recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 256 elements
 * @param[in] f Array of 16 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 16)
 */
void
fft_1024_16(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs)
{
        for (size_t i = 1; i < (1 << 4); ++i)
                f[i] = gf_mul(betas_pows[64+i], f[i]);

        uint16_t f0[1 << 4-1];
        uint16_t f1[1 << 4-1];

        radix_1024_16 (f0, f1, f);

        uint16_t u[1 << 8-1];
        uint16_t v[1 << 8-1];

        fft_1024_8 (u, f0, (f_coeffs+1)/2);
        fft_1024_8 (v, f1, f_coeffs/2);

        size_t k = 1 << 8-1;
        memcpy(w + k, v, 2*k);
        w[0] = u[0];
        w[k] ^= u[0];
        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(gammas_sums[512+256+i], v[i]);
                w[k+i] ^= w[i];
        }
}

/**
 * First recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 512 elements
 * @param[in] f Array of 32 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 32)
 */
void
fft_1024_32(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs)
{
#ifndef FFT_BETAS_1024
        for (size_t i = 1; i < (1 << 5); ++i)
                f[i] = gf_mul(betas_pows[i], f[i]);
#endif
  
        uint16_t f0[1 << 5-1];
        uint16_t f1[1 << 5-1];

        radix_1024_32 (f0, f1, f);

        uint16_t u[1 << 9-1];
        uint16_t v[1 << 9-1];

        fft_1024_16 (u, f0, (f_coeffs+1)/2);
        fft_1024_16 (v, f1, f_coeffs/2);

        size_t k = 1 << 9-1;
        memcpy(w + k, v, 2*k);
        w[0] = u[0];
        w[k] ^= u[0];
        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(gammas_sums[512+i], v[i]);
                w[k+i] ^= w[i];
        }
}

/**
 * Additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * Evaluates polynomial f of degree f_coeffs-1 (less than 64) on all elements of GF(2^10).
 * @param[out] w Array of 1024 elements
 * @param[in] f Array of 64 elements
 * @param[in] f_coeffs Number of coefficients of f (less than 64)
 */
void
fft_1024_64(uint16_t    *w,
            uint16_t    *f,
            size_t       f_coeffs)
{
        uint16_t f0[32];
        uint16_t f1[32];

        radix_1024_64 (f0, f1, f);

        uint16_t u[512];
        uint16_t v[512];

        fft_1024_32 (u, f0, (f_coeffs+1)/2);
        fft_1024_32 (v, f1, f_coeffs/2);

        size_t k = 1 << 10-1;
        memcpy(w + k, v, 2*k);
        w[0] = u[0];
        w[k] ^= u[0];

        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);    
                w[k+i] ^= w[i];
        }
}

#else // !defined FFT_UNROLL_1024_64

/**
 * Computes the radix conversion of a polynomial f in GF(2^m)[x].
 * Computes f0 and f1 such that f(x) = f0(x^2-x) + x.f1(x^2-x)
 * as proposed by Bernstein, Chou and Schwabe:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 * @param[out] f0 Array half the size of f
 * @param[out] f1 Array half the size of f
 * @param[in] f Array of size a power of 2
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to the number of coefficients of f
 */
void
radix(uint16_t          *f0,
      uint16_t          *f1,
      const uint16_t    *f,
      unsigned int       m_f)
{
        switch (m_f) {
        case 4:
                f0[4] = f[8]^f[12];
                f0[6] = f[12]^f[14];
                f0[7] = f[14]^f[15];
                f1[5] = f[11]^f[13];
                f1[6] = f[13]^f[14];
                f1[7] = f[15];
                f0[5] = f[10]^f[12]^f1[5];
                f1[4] = f[9]^f[13]^f0[5];

                f0[0] = f[0];
                f1[3] = f[7]^f[11]^f[15];
                f0[3] = f[6]^f[10]^f[14]^f1[3];
                f0[2] = f[4]^f0[4]^f0[3]^f1[3];
                f1[1] = f[3]^f[5]^f[9]^f[13]^f1[3];
                f1[2] = f[3]^f1[1]^f0[3];
                f0[1] = f[2]^f0[2]^f1[1];
                f1[0] = f[1]^f0[1];
                return;

        case 3:
                f0[0] = f[0];
                f0[2] = f[4]^f[6];
                f0[3] = f[6]^f[7];
                f1[1] = f[3]^f[5]^f[7];
                f1[2] = f[5]^f[6];
                f1[3] = f[7];
                f0[1] = f[2]^f0[2]^f1[1];
                f1[0] = f[1]^f0[1];
                return;
  
        case 2:
                f0[0] = f[0];
                f0[1] = f[2]^f[3];
                f1[0] = f[1]^f0[1];
                f1[1] = f[3];
                return;

        case 1:
                f0[0] = f[0];
                f1[0] = f[1];
                return;
  
        default:;
                size_t n = 1 << m_f-2;

                uint16_t Q[2*n];
                uint16_t R[2*n];

                memcpy(Q, f + 3*n, 2*n);
                memcpy(Q + n, f + 3*n, 2*n);
                memcpy(R, f, 4*n);

                for (size_t i = 0; i < n; ++i) {
                        Q[i] ^= f[2*n + i];
                        R[n + i] ^= Q[i];
                }

                uint16_t Q0[n];
                uint16_t Q1[n];
                uint16_t R0[n];
                uint16_t R1[n];

                radix(Q0, Q1, Q, m_f-1);
                radix(R0, R1, R, m_f-1);

                memcpy(f0, R0, 2*n);
                memcpy(f0 + n, Q0, 2*n);
                memcpy(f1, R1, 2*n);
                memcpy(f1 + n, Q1, 2*n);
        }
}

#if (defined COMPUTE_FFT_LUT) || (defined FFT_LUT_1024_64)

/**
 * Evaluates f at all subset sums of a given set.
 * This function is a subroutine of function fft(...).
 * @param[out] w Array receiving the evaluations
 * @param[in] f Array storing the coefficients of the polynomial to evaluate
 * @param[in] f_coeffs Number of coefficients of f
 * @param[in] m 2^m is the smallest power of 2 greater or equal to the length of family w
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to f_coeffs
 * @param[in] gammas_sums Logarithms of subset sums of gammas needed at each recursive step.
 * @param[in] betas_pows Logarithms of powers of beta^{m} needed at each recursive step
 * @param[in] fft_final_betas Logarithms of betas at the last step
 */
void
fft_rec(uint16_t        *w,
        uint16_t        *f,
        size_t           f_coeffs,
        unsigned int     m,
        unsigned int     m_f,
        const uint16_t  *gammas_sums,
        const uint16_t  *betas_pows,
        const uint16_t  *fft_final_betas)
{
        // Step 1
        if (m_f == 1) {
                uint16_t tmp[GF_M - (FFT_PARAM-1)];
                for (size_t i = 0; i < m; ++i)
                        tmp[i] = gf_mul(fft_final_betas[i], f[1]);
                w[0] = f[0];
                for (size_t j = 0; j < m; ++j) {
                        for (size_t k = 0; k < 1 << j; ++k) {
                                w[(1 << j) + k] = w[k] ^ tmp[j];
                        }
                }
                
                return;
        }
  
        // Step 2: compute g
        if (betas_pows[1] != 1) {
                for (size_t i = 1; i < (1 << m_f); ++i)
                        f[i] = gf_mul(betas_pows[i], f[i]);
        }

        // Step 3
        uint16_t f0[1 << m_f-1];
        uint16_t f1[1 << m_f-1];

        radix(f0, f1, f, m_f);

        // Step 4: constants gammas and deltas have been precomputed

        // Step 5
        size_t k = 1 << m-1;
        uint16_t u[1 << m-1];
        uint16_t v[1 << m-1];
        
        fft_rec(u, f0, (f_coeffs+1)/2, m-1, m_f-1, gammas_sums + k, betas_pows + (1 << m_f+1), fft_final_betas);
        
        if (f_coeffs <= 3) {  // 3-coefficient polynomial f case: f1 is constant
                w[0] = u[0];
                w[k] = u[0] ^ f1[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[i], f1[0]);
                        w[k+i] = w[i] ^ f1[0];
                }
        }
        else {
                fft_rec(v, f1, f_coeffs/2, m-1, m_f-1, gammas_sums + k, betas_pows + (1 << m_f+1), fft_final_betas);

                // Step 6
                memcpy(w + k, v, 2*k);
                w[0] = u[0];
                w[k] ^= u[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);
                        w[k+i] ^= w[i];
                }
        }
}

/**
 * Evaluates f on all fields elements using an additive FFT algorithm.
 * f_coeffs is the number of coefficients of f (one less than its degree). <br>
 * The FFT proceeds recursively to evaluate f at all subset sums of a basis B. <br>
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.clemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf <br>
 * Constants betas, gammas and deltas involved in the algorithm are either hardcoded or precomputed
 * by the subroutine compute_fft_lut(...). <br>
 * Note that on this first call (as opposed to the recursive calls to fft_rec), gammas are equal to betas,
 * meaning the first gammas subset sums are actually the subset sums of betas (except 1). <br>
 * Also note that f is altered during computation (twisted at each level).
 * @param[out] w Array receiving the images of the subset sums by f
 * @param[in] f Array of 2^FFT_PARAM elements
 * @param[in] f_coeffs Number coefficients of f (i.e. deg(f)+1)
 */
void
fft(uint16_t    *w,
    uint16_t    *f,
    size_t       f_coeffs)
{
        // Follows Gao and Mateer algorithm

        // Step 1
#if FFT_PARAM == 1
        uint16_t tmp[GF_M-1];
        for (size_t i = 0; i < GF_M-1; ++i)
                tmp[i] = gf_mul(f[1], fft_final_betas[i]);
        
        w[0] = f[0];
        size_t l = 1 << GF_M-1;
        w[l] = w[0] ^ f[1];

        for (size_t j = 0; j < GF_M-1; ++j) {
                for (size_t k = 0; k < 1 << j; ++k) {
                        size_t i = (1 << j) + k;      
                        w[i] = w[k] ^ tmp[j];
                        w[l+i] = w[i] ^ f[1];
                }
        }
    
        return;
#else        
        // Step 2: beta_m = 1, nothing to do
  
        // Step 3
        uint16_t f0[1 << FFT_PARAM-1];
        uint16_t f1[1 << FFT_PARAM-1];
        
        radix(f0, f1, f, FFT_PARAM);

        // Step 4: beta_m = 1, nothing to do
  
        // Step 5
        size_t k = 1 << GF_M-1;
        uint16_t u[k];
        uint16_t v[k];

        fft_rec(u, f0, (f_coeffs+1)/2, GF_M-1, FFT_PARAM-1, gammas_sums + k, betas_pows, fft_final_betas);
        fft_rec(v, f1, f_coeffs/2, GF_M-1, FFT_PARAM-1, gammas_sums + k, betas_pows, fft_final_betas);
  
        // Step 6, 7

        memcpy(w + k, v, 2*k);

        // Check if 0 is root
        w[0] = u[0];

        // Check if 1 is root
        w[k] ^= u[0];

        // Find other roots
        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);
                w[k+i] ^= w[i];
        }
#endif
}

#else // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)

/**
 * Evaluates f at all subset sums of a given set.
 * This function is a subroutine of the function fft(...).
 * @param[out] w Array
 * @param[in] f Array
 * @param[in] f_coeffs Number of coefficients of f
 * @param[in] m Number of betas
 * @param[in] m_f Number of coefficients of f (one more than its degree)
 * @param[in] betas FFT constants
 */
void
fft_rec(uint16_t        *w,
        uint16_t        *f,
        size_t           f_coeffs,
        unsigned int     m,
        unsigned int     m_f,
        const uint16_t  *betas)
{
        // Step 1
        if (m_f == 1) {
                uint16_t tmp[GF_M - (FFT_PARAM-1)];
                for (size_t i = 0; i < m; ++i)
                        tmp[i] = gf_mul(betas[i], f[1]);

                w[0] = f[0];
                for (size_t j = 0; j < m; ++j) {
                        for (size_t k = 0; k < 1 << j; ++k) {
                                w[(1 << j) + k] = w[k] ^ tmp[j];
                        }
                }
    
                return;
        }
  
        // Step 2: compute g
        if (betas[m-1] != 1) {
                uint16_t beta_m_pow = 1;
                for (size_t i = 1; i < (1 << m_f); ++i) {
                        beta_m_pow = gf_mul(beta_m_pow, betas[m-1]);
                        f[i] = gf_mul(beta_m_pow, f[i]);
                }
        }

        // Step 3
        uint16_t f0[1 << m_f-1];
        uint16_t f1[1 << m_f-1];

        radix(f0, f1, f, m_f);

        // Step 4: compute gammas and deltas
        uint16_t gammas[m-1];
        uint16_t deltas[m-1];

        for (size_t i = 0; i < m-1; ++i) {
                gammas[i] = gf_mul(betas[i], gf_inverse(betas[m-1]));
                deltas[i] = gf_square(gammas[i]) ^ gammas[i];
        }

        size_t k = 1 << m-1;
  
        // Compute gammas sums
        uint16_t gammas_sums[k];
        compute_subset_sums(gammas_sums, gammas, m-1);
        
        // Step 5
        uint16_t u[1 << m-1];
        uint16_t v[1 << m-1];

        fft_rec(u, f0, (f_coeffs+1)/2, m-1, m_f-1, deltas);
        
        if (f_coeffs <= 3) {  // 3-coefficient polynomial f case: f1 is constant
                w[0] = u[0];
                w[k] = u[0] ^ f1[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[i], f1[0]);
                        w[k+i] = w[i] ^ f1[0];
                }
        }
        else {
                fft_rec(v, f1, f_coeffs/2, m-1, m_f-1, deltas);

                // Step 6
                memcpy(w + k, v, 2*k);
                w[0] = u[0];
                w[k] ^= u[0];
                for (size_t i = 1; i < k; ++i) {
                        w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);
                        w[k+i] ^= w[i];
                }
        }
}

/**
 * Evaluates f on all fields elements using an additive FFT algorithm.
 * f_coeffs is the number of coefficients of f (one less than its degree). <br>
 * The FFT proceeds recursively to evaluate f at all subset sums of a basis B. <br>
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.clemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf <br>
 * Constants betas, gammas and deltas involved in the algorithm are either hardcoded or precomputed
 * by the subroutine compute_fft_lut(...). <br>
 * Note that on this first call (as opposed to the recursive calls to fft_rec), gammas are equal to betas,
 * meaning the first gammas subset sums are actually the subset sums of betas (except 1). <br>
 * Also note that f is altered during computation (twisted at each level).
 * @param[out] w Array
 * @param[in] f Array of 2^FFT_PARAM elements
 * @param[in] f_coeffs Number coefficients of f (i.e. deg(f)+1)
 */
void
fft(uint16_t    *w,
    uint16_t    *f,
    size_t       f_coeffs)
{
        // Follows Gao and Mateer algorithm
        uint16_t betas[GF_M-1];
        compute_fft_betas(betas);
        
        // Step 1
#if FFT_PARAM == 1
        uint16_t tmp[GF_M-1];
        for (size_t i = 0; i < GF_M-1; ++i)
                tmp[i] = gf_mul(f[1], betas[i]);
        
        w[0] = f[0];
        size_t l = 1 << GF_M-1;
        w[l] = w[0] ^ f[1];
        
        for (size_t j = 0; j < GF_M-1; ++j) {
                for (size_t k = 0; k < 1 << j; ++k) {
                        size_t i = (1 << j) + k;
                        w[i] = w[k] ^ tmp[j];    
                        w[i+l] = w[i] ^ f[1];
                }
        }
    
        return;
#else
        // Compute gammas sums
        uint16_t betas_sums[1 << GF_M-1];
        compute_subset_sums(betas_sums, betas, GF_M-1);
        
        // Step 2: beta_m = 1, nothing to do
  
        // Step 3
        uint16_t f0[1 << FFT_PARAM-1];
        uint16_t f1[1 << FFT_PARAM-1];

        radix(f0, f1, f, FFT_PARAM);
  
        // Step 4: Compute deltas
        uint16_t deltas[GF_M-1];
        for (size_t i = 0; i < GF_M-1; ++i)
                deltas[i] = gf_square(betas[i]) ^ betas[i];
  
        // Step 5
        size_t k = 1 << GF_M-1;
        uint16_t u[k];
        uint16_t v[k];

        fft_rec(u, f0, (f_coeffs+1)/2, GF_M-1, FFT_PARAM-1, deltas);
        fft_rec(v, f1, f_coeffs/2, GF_M-1, FFT_PARAM-1, deltas);
  
        // Step 6, 7 and error polynomial computation

        memcpy(w + k, v, 2*k);

        // Check if 0 is root
        w[0] = u[0];

        // Check if 1 is root
        w[k] ^= u[0];

        // Find other roots
        for (size_t i = 1; i < k; ++i) {
                w[i] = u[i] ^ gf_mul(betas_sums[i], v[i]);
                w[k+i] ^= w[i];
        }
#endif
}

#endif // (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
#endif // !defined FFT_UNROLL_1024_64

/**
 * Arranges the received word rcv in a form w such that applying the additive FFT transpose to w
 * yields the BCH syndromes of the received word rcv.
 * Since the received word rcv gives coefficients of the primitive element alpha, we twist accordingly. <br>
 * Furthermore, the additive FFT transpose needs elements indexed by their decomposition on the chosen basis,
 * so we apply the adequate permutation.
 * @param[out] w Array of size 2^GF_M
 * @param[in] rcv Array of size BCH_N_BYTES
 */
void
fft_t_preprocess_bch_codeword(uint16_t          *w,
                              const uint8_t     *rcv)
{
        // Unpack the received word rcv into array r
        uint16_t r[1 << GF_M];
        size_t i;
        for (i = 0; i < BCH_N_BYTES-(BCH_N%8!=0); ++i)
                for (size_t j = 0; j < 8; ++j)
                        r[8*i+j] = (rcv[i] >> j) & 1;

        for (size_t j = 0; j < BCH_N%8; ++j) // last byte
                r[8*i+j] = (rcv[i] >> j) & 1;

        memset(r+BCH_N, 0, 2*((1<<GF_M)-BCH_N)); // complete r with zeros

#if (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
        uint16_t gammas[GF_M-1];
        uint16_t gammas_sums[1 << GF_M-1];

        compute_fft_betas(gammas);
        compute_subset_sums(gammas_sums, gammas, GF_M-1);
#endif

        // Twist and permute r adequately to obtain w
        size_t k = 1 << GF_M-1;
        w[0] = 0;
        w[k] = -r[0] & 1;
        for (size_t i = 1; i < k; ++i) {
                w[i] = -r[gf_log(gammas_sums[i])] & gammas_sums[i];
                w[k+i] = -r[gf_log(gammas_sums[i]^1)] & (gammas_sums[i]^1);
        }
}

/**
 * Retrieves the error polynomial err from the evaluations w of the ELP
 * (Error Locator Polynomial) on all field elements.
 * @param[out] err Array of size BCH_N_BYTES
 * @param[in] w Array of size 2^GF_M
 */
void
fft_retrieve_bch_error_poly(uint8_t             *err,
                            const uint16_t      *w)
{
#if (!defined COMPUTE_FFT_LUT) && (!defined FFT_LUT_1024_64)
        uint16_t gammas[GF_M-1];
        uint16_t gammas_sums[1 << GF_M-1];

        compute_fft_betas(gammas);
        compute_subset_sums(gammas_sums, gammas, GF_M-1);
#endif

        size_t k = 1 << GF_M-1;
        err[0] ^= 1 ^ ((uint16_t)-w[0] >> 15);
        size_t index = GF_MUL_ORDER;
        uint8_t bit = 1 ^ ((uint16_t)-w[k] >> 15);
        err[index/8] ^= bit << (index%8);
  
        for (size_t i = 1; i < k; ++i) {
                index = GF_MUL_ORDER - gf_log(gammas_sums[i]);
                bit = 1 ^ ((uint16_t)-w[i] >> 15);
                err[index/8] ^= bit << (index%8);

                index = GF_MUL_ORDER - gf_log(gammas_sums[i] ^ 1);
                bit = 1 ^ ((uint16_t)-w[k+i] >> 15);
                err[index/8] ^= bit << (index%8);
        }
}
