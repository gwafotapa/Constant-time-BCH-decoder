/**
 * @file bch.c
 * Constant time implementation of BCH codes
 */

#include "bch.h"

#include "fft.h"
#include "gf.h"
#include "optimizations.h"

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/**
 * Returns i modulo the given modulus.
 * i must be less than 2*modulus.
 * Therefore, the return value is either i or i-modulus.
 * @returns i mod (modulus)
 * @param[in] i The integer whose modulo is taken
 * @param[in] modulus The modulus
 */
uint16_t
mod(uint16_t    i,
    uint16_t    modulus)
{
        uint16_t tmp = i - modulus;
  
        // mask = 0xffff if (i < GF_MUL_ORDER)
        int16_t mask = -(tmp >> 15);
  
        return tmp + (mask & modulus);  
}

/**
 * Computes the odd binary cyclotomic cosets modulo 2^m-1 for integers less than upper_bound.
 * The array cosets of size 2^m-1 is filled by placing at index i the coset representative of i.
 * @param[out] cosets Array receiving the coset representatives
 * @param[in] upper_bound The upper bound
 * @param[in] m Parameter of Galois field GF(2^m)
 */
void
compute_cyclotomic_cosets(uint16_t      *cosets,
                          uint16_t       upper_bound,
                          size_t         m)
{
        // Compute the odd cyclotomic classes 
        for (uint16_t i = 1; i < upper_bound; i += 2) {
                if (cosets[i] == 0) { // If i does not already belong to a class
                        cosets[i] = i;
                        uint16_t tmp = i;
                        size_t gf_mul_order = (1 << m) - 1;
                        size_t j = m;
                        while (--j) { // Complete i's class
                                tmp = mod(2*tmp, gf_mul_order);
                                cosets[tmp] = i;
                        }
                }
        }
}

/**
 * Computes the generator polynomial of the primitive BCH code with given parameters.
 * Code length is 2^m-1. <br>
 * Parameter t is the targeted correction capacity of the code
 * and receives the real correction capacity (which is at least equal to the target). <br>
 * exp and log are arrays giving antilog and log of GF(2^m) elements.
 * @returns the degree of the generator polynomial
 * @param[out] bch_poly Array of size (m*t + 1) receiving the coefficients of the generator polynomial
 * @param[in,out] t Targeted correction capacity; receives the real correction capacity
 * @param[in] m Parameter of Galois field GF(2^m)
 * @param[in] exp Antilog table of GF(2^m)
 * @param[in] log Log table of GF(2^m)
 */
size_t
compute_bch_poly(uint16_t       *bch_poly,
                 size_t         *t,
                 size_t          m,
                 const uint16_t *exp,
                 const uint16_t *log)
{
        size_t gf_mul_order = (1 << m) - 1;
        uint16_t cosets[gf_mul_order];
        memset(cosets, 0, 2 * gf_mul_order);
        compute_cyclotomic_cosets(cosets, 2 * *t, m);
  
        // Start with bch_poly(X) = 1
        bch_poly[0] = 1;
        size_t deg_bch_poly = 0;

        for (uint16_t i = 1; i < gf_mul_order; ++i) {
                if (cosets[i] == 0)
                        continue;
    
                // Multiply bch_poly(X) by X-a^i
                for (size_t j = deg_bch_poly; j; --j) {
                        int16_t mask = -((uint16_t)-bch_poly[j] >> 15);
                        bch_poly[j] = (mask & exp[mod(log[bch_poly[j]] + i, gf_mul_order)]) ^ bch_poly[j-1];
                }
                bch_poly[0] = exp[mod(log[bch_poly[0]] + i, gf_mul_order)];
                bch_poly[++deg_bch_poly] = 1;
        }

        // Determine the real correction capacity
        while (cosets[2 * *t + 1] != 0)
                ++*t;

        return deg_bch_poly;
}

/**
 * Unpacks the message msg to the array msg_unpacked where each byte stores a bit of the message.
 * @param[out] msg_unpacked Array of BCH_K_BYTES bytes receiving the unpacked message
 * @param[in] msg Array of BCH_K bytes storing the packed message
 */
void
unpack_message(uint8_t          *msg_unpacked,
               const uint8_t    *msg)
{ 
        for (size_t i = 0; i < (BCH_K_BYTES - (BCH_K%8 != 0)); ++i) {
                for (size_t j = 0; j < 8; ++j) {
                        msg_unpacked[j + 8*i] = (msg[i] >> j) & 0x01;
                }
        }

        for (size_t j = 0; j < BCH_K%8; ++j) {
                msg_unpacked[j + 8*(BCH_K_BYTES-1)] = (msg[BCH_K_BYTES-1] >> j) & 0x01;
        }
}

/**
 * Encodes the message msg to a codeword cdw using the generator polynomial bch_poly of the code.
 * @param[out] cdw Array of BCH_N bytes receiving the codeword
 * @param[in] bch_poly Array of at least BCH_POLY_SIZE bytes storing the generator polynomial
 * @param[in] msg Array of BCH_K bytes storing the message to encode
 */
void
lfsr_encode(uint8_t             *cdw,
            const uint8_t       *bch_poly,
            const uint8_t       *msg)
{
        uint8_t gate_value = 0;
  
        // Compute the parity-check digits
        for (int i = BCH_K-1; i >= 0; --i) {
                gate_value = msg[i] ^ cdw[BCH_N - BCH_K - 1];

                if (gate_value) {
                        for (size_t j = BCH_N - BCH_K - 1; j; --j)
                                cdw[j] = cdw[j-1] ^ bch_poly[j];
                } else {
                        for (size_t j = BCH_N - BCH_K - 1; j; --j)
                                cdw[j] = cdw[j-1];
                }

                cdw[0] = gate_value; 
        }
  
        // Add the message 
        memcpy(cdw + BCH_N - BCH_K, msg, BCH_K);
}

/**
 * Packs the codeword from an array cdw_unpacked where each byte stores a bit to a compact array cdw.
 * @param[out] cdw Array of BCH_N_BYTES bytes receiving the packed codeword
 * @param[in] cdw_unpacked Array of BCH_N bytes storing the unpacked codeword
 */
void
pack_codeword(uint8_t           *cdw,
              const uint8_t     *cdw_unpacked)
{
        for (size_t i = 0; i < (BCH_N_BYTES - (BCH_N%8 != 0)); ++i) {
                for (size_t j = 0; j < 8; ++j) {
                        cdw[i] |= cdw_unpacked[j + 8*i] << j;
                }
        }

        for (size_t j = 0; j < BCH_N%8; ++j) {
                cdw[BCH_N_BYTES - 1] |= cdw_unpacked[j + 8*(BCH_N_BYTES - 1)] << j;
        }
}

/**
 * Encodes a message msg of BCH_K bits to a BCH codeword cdw of BCH_N bits.
 * Following @cite lin1983error (Chapter 4 - Cyclic Codes),
 * we perform a systematic encoding using a linear (BCH_N - BCH_K)-stage shift register
 * with feedback connections based on the generator polynomial bch_poly of the BCH code.
 * @param[out] cdw Array of size BCH_N_BYTES receiving the encoded message
 * @param[in] bch_poly Array of size at least BCH_POLY_SIZE storing the generator polynomial
 * @param[in] msg Array of size BCH_K_BYTES storing the message
 */
void
bch_encode(uint8_t              *cdw,
           const uint8_t        *bch_poly,
           const uint8_t        *msg)
{
        uint8_t msg_unpacked[BCH_K];
        uint8_t cdw_unpacked[BCH_N] = {0};
        
        unpack_message(msg_unpacked, msg);
        lfsr_encode(cdw_unpacked, bch_poly, msg_unpacked);
        pack_codeword(cdw, cdw_unpacked);
}

/**
 * Computes the error locator polynomial (ELP) sigma.
 * This is a constant time implementation of Berlekamp's simplified algorithm (see @cite joiner1995decoding). <br>
 * We use the letter p for rho which is initialized at -1/2. <br>
 * The array X_sigma_p represents the polynomial X^(2(mu-rho))*sigma_p(X). <br>
 * Instead of maintaining a list of sigmas, we update in place both sigma and X_sigma_p. <br>
 * sigma_copy serves as a temporary save of sigma in case X_sigma_p needs to be updated. <br>
 * We can properly correct only if the degree of sigma does not exceed BCH_DELTA.
 * This means only the first BCH_DELTA + 1 coefficients of sigma are of value
 * and we only need to save its first BCH_DELTA - 1 coefficients.
 * @returns the degree of the ELP sigma
 * @param[out] sigma Array of size (at least) BCH_DELTA receiving the ELP
 * @param[in] syndromes Array of size (at least) 2*BCH_DELTA storing the syndromes
 */
size_t
compute_elp(uint16_t            *sigma,
            const uint16_t      *syndromes)
{
#if BCH_DELTA == 1
        sigma[0] = 1;
        sigma[1] = syndromes[0];
        size_t deg_sigma = (uint16_t)-syndromes[0] >> 15;
        
        return deg_sigma;
#else
        sigma[0] = 1;
        size_t deg_sigma = 0;
        size_t deg_sigma_p = 0;
        uint16_t sigma_copy[BCH_DELTA - 1] = {0};
        size_t deg_sigma_copy = 0;
        uint16_t X_sigma_p[BCH_DELTA + 1] = {0, 1};
        int pp = -1; // 2*rho
        uint16_t d_p = 1;
        uint16_t d = syndromes[0];
        
        for (size_t mu = 0; mu < BCH_DELTA; ++mu) {
                // Save sigma in case we need it to update X_sigma_p
                memcpy(sigma_copy, sigma, 2*(BCH_DELTA-1));
                deg_sigma_copy = deg_sigma;

                uint16_t dd = gf_mul(d, gf_inverse(d_p)); // 0 if (d == 0)
                for (size_t i = 1; (i <= 2*mu+1) && (i <= BCH_DELTA); ++i)
                        sigma[i] ^= gf_mul(dd, X_sigma_p[i]);
    
                size_t deg_X = 2*mu-pp; // 2*(mu-rho)
                size_t deg_X_sigma_p = deg_X + deg_sigma_p;
                
                // mask1 = 0xffff if (d != 0) and 0 otherwise
                int16_t mask1 = -((uint16_t)-d >> 15);

                // mask2 = 0xffff if (deg_X_sigma_p > deg_sigma) and 0 otherwise
                int16_t mask2 = -((uint16_t)(deg_sigma - deg_X_sigma_p) >> 15);

                // mask12 = 0xffff if the deg_sigma increased and 0 otherwise
                int16_t mask12 = mask1 & mask2;
                deg_sigma = (mask12 & deg_X_sigma_p) ^ (~mask12 & deg_sigma);

                if (mu == BCH_DELTA-1)
                        break;

                // Update pp, d_p and X_sigma_p if needed
                pp = (mask12 & (2*mu)) ^ (~mask12 & pp);
                d_p = (mask12 & d) ^ (~mask12 & d_p);
                for (size_t i = BCH_DELTA-1; i; --i)
                        X_sigma_p[i+1] = (mask12 & sigma_copy[i-1]) ^ (~mask12 & X_sigma_p[i-1]);
                X_sigma_p[1] = 0;
                X_sigma_p[0] = 0;
                deg_sigma_p = (mask12 & deg_sigma_copy) ^ (~mask12 & deg_sigma_p);

                // Compute the next discrepancy
                d = syndromes[2*mu+2];
                for (size_t i = 1; (i <= 2*mu+1) && (i <= BCH_DELTA); ++i)
                        d ^= gf_mul(sigma[i], syndromes[2*mu+2-i]);
        }
  
        return deg_sigma;
#endif
}

/**
 * Retrieves the message msg from the codeword cdw.
 * Since we performed a systematic encoding, the message is the last BCH_K bits of the codeword.
 * @param[out] msg Array of size BCH_K_BYTES receiving the message
 * @param[in] cdw Array of size BCH_N_BYTES storing the codeword
 */
void
message_from_codeword(uint8_t           *msg,
                      const uint8_t     *cdw)
{
        size_t val = BCH_N - BCH_K;

        uint8_t mask1 = 0xff << val%8;
        uint8_t mask2 = 0xff >> 8-val%8;
        size_t index = val/8;

        for (size_t i = 0; i < BCH_K_BYTES - 1; ++i) {
                uint8_t msg1 = (cdw[index] & mask1) >> val%8;
                uint8_t msg2 = (cdw[++index] & mask2) << 8-val%8;
                msg[i] = msg1 | msg2;
        }

        // Last byte
        if ((BCH_K%8 == 0) || (8-val%8 < BCH_K%8)) { // 8-val%8 is the number of bits given by m1
                uint8_t msg1 = (cdw[index] & mask1) >> val%8;
                uint8_t msg2 = (cdw[++index] & mask2) << 8-val%8;
                msg[BCH_K_BYTES-1] = msg1 | msg2;
        }
        else {
                uint8_t msg1 = (cdw[index] & mask1) >> val%8;
                msg[BCH_K_BYTES-1] = msg1;
        }
}

/**
 * Computes the 2^BCH_DELTA syndromes from the received vector rcv,
 * that is the sum of powers of alpha weighted by rcv's coefficients.
 * To do so, we use the additive FFT transpose, which takes as input a family w of GF(2^GF_M) elements
 * and outputs the weighted power sums of these w. <br>
 * Therefore, this requires twisting and applying a permutation before feeding rcv to the fft transpose. <br>
 * For more details see Berstein, Chou and Schawbe's explanations:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 * @param[out] syndromes Array of size 2^(FFT_T_PARAM) receiving the 2*BCH_DELTA syndromes
 * @param[in] rcv Array of size BCH_N_BYTES storing the received word
 */
void
compute_syndromes(uint16_t      *syndromes,
                  const uint8_t *rcv)
{
        uint16_t w[1 << GF_M];
        fft_t_preprocess_bch_codeword(w, rcv);
  
#ifdef FFT_T_UNROLL_1024_64
        fft_t_1024_128(syndromes, w, 2*BCH_DELTA);
#else
        fft_t(syndromes, w, 2*BCH_DELTA);
#endif
}

/**
 * Computes the error polynomial err from the error locator polynomial sigma.
 * See function fft(...) for more details.
 * @param[out] err Array of BCH_N_BYTES elements receiving the error polynomial
 * @param[in] sigma Array of 2^FFT_PARAM elements storing the error locator polynomial
 * @param[in] deg_sigma The degree of sigma
 */
void
compute_roots(uint8_t   *err,
              uint16_t  *sigma,
              size_t     deg_sigma)
{
        // w will receive the evaluation of sigma in all field elements
        uint16_t w[1 << GF_M] = {0};
  
#ifdef FFT_UNROLL_1024_64
        fft_1024_64(w, sigma, BCH_DELTA+1);
#else
        fft(w, sigma, BCH_DELTA+1);
#endif
  
        fft_retrieve_bch_error_poly(err, w);
}

/**
 * Adds two vectors of the given size (XORing coordinate-wise).
 * @param[out] sum Array receiving the sum
 * @param[in] v1 Array storing the first vector to add
 * @param[in] v2 Array storing the second vector to add
 * @param[in] size Size of the vectors to be added
 */
void
vect_add(uint8_t        *sum,
         const uint8_t  *v1,
         const uint8_t  *v2,
         size_t          size)
{
        for (size_t i = 0; i < size; ++i)
                sum[i] = v1[i] ^ v2[i];
}

/**
 * Decodes the received word in four steps:
 *    <ol>
 *    <li> The first step, done by additive FFT transpose, is the computation of the 2 * BCH_DELTA syndromes.
 *    <li> The second step is the computation of the error-locator polynomial sigma.
 *    <li> The third step, done by additive FFT, is finding the error-locator numbers by calculating the roots of the polynomial sigma and takings their inverses.
 *    <li> The fourth step is the correction of the errors in the received polynomial.
 *    </ol>
 * For a more complete picture on BCH decoding, see Shu. Lin and Daniel J. Costello in Error Control Coding: Fundamentals and Applications @cite lin1983error
 * @param[out] msg Array of size BCH_K_BYTES receiving the decoded message
 * @param[in] rcv Array of size BCH_N_BYTES storing the received word
 */
void
bch_decode(uint8_t      *msg,
           uint8_t      *rcv)
{
#ifndef GF_LUT_1024
        gf_generate_M();
#endif

#ifdef COMPUTE_FFT_LUT
        compute_fft_lut();
#endif

        // Calculate the 2 * BCH_DELTA syndromes
        uint16_t syndromes[1 << FFT_T_PARAM];
        compute_syndromes(syndromes, rcv);

        // Compute the error locator polynomial sigma.
        // sigma's degree is at most BCH_DELTA but the FFT requires the extra room.
        uint16_t sigma[1 << FFT_PARAM] = {0};
        size_t deg_sigma = compute_elp(sigma, syndromes);

        // Compute the error polynomial err
        uint8_t err[(1 << GF_M) / 8] = {0};
        compute_roots(err, sigma, deg_sigma);
  
        // Add the error polynomial to the received polynomial 
        vect_add(rcv, rcv, err, BCH_N_BYTES);
  
        // Retrieve the message from the decoded codeword
        message_from_codeword(msg, rcv);
}
