/**
 * @file time_once.c
 * Measures time of one execution of a BCH decoding
 */

#include "bch.h"
#include "fft.h"
#include "gf.h"
#include "optimizations.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Reads the time from the time stamp counter.
 * @returns time read in CPU cycles
 */
uint64_t
cpu_cycles()
{
        uint64_t tsc;
        __asm__ __volatile__("rdtsc           \n\t"
                             "shlq $32, %%rdx \n\t"
                             "orq %%rdx, %%rax"
                             : "=a" (tsc)
                             :
                             : "rdx");
  
        return tsc;
}

/**
 * Converts an hexadecimal string into a byte string.
 * Inspired by https://gist.github.com/xsleonard/7341172
 * @returns -1 if the input string has an incorrect length and 0 otherwise
 * @param[out] uint8str Array of size BCH_N_BYTES receiving the byte string
 * @param[in@ hexstr Array of size 2*BCH_N_BYTES+1 storing the input hexadecimal string
 */
int
hexstr_to_uint8str(uint8_t      *uint8str,
                   char const   *hexstr)
{
        size_t len = strlen(hexstr);
        if (len != 2 * BCH_N_BYTES)
                return -1;
  
        for (size_t i = 0, j = 0; j < BCH_N_BYTES; i += 2, ++j)
                uint8str[j] = (hexstr[i] % 32 + 9) % 25 * 16 + (hexstr[i+1] % 32 + 9) % 25;

        return 0;
}

int
main(int         argc,
     char       *argv[])
{
        if (argc != 2) {
                fprintf(stderr, "error: argument expected.\n");
                exit(1);
        }

        uint8_t rcv[BCH_N_BYTES]; // received vector
        if (hexstr_to_uint8str(rcv, argv[1]) == -1) {
                fprintf(stderr, "error: wrong received vector length.\n");
                exit(1);
        }
        uint8_t msg[BCH_K_BYTES] = {0}; // decoded message
  
        // Time measures (in CPU cycles) of the five different parts of the function bch_decode.
        uint64_t cycles[5] = {0};

  
        /*** Beginning of bch_decode function ***/

        // Time precomputation step
        uint64_t clock = cpu_cycles();
  
#ifndef GF_LUT_1024
        gf_generate_M();
#endif
  
#ifdef COMPUTE_FFT_LUT
        compute_fft_lut();
#endif
  
        cycles[0] = cpu_cycles() - clock;

        // Time syndromes generation
        clock = cpu_cycles();
        uint16_t syndromes[1 << FFT_T_PARAM] = {0};
        compute_syndromes(syndromes, rcv);
        cycles[1] = cpu_cycles() - clock;

        // Time error locator polynomial computation
        clock = cpu_cycles();
        uint16_t sigma[1 << FFT_PARAM] = {0};
        size_t deg_sigma = compute_elp(sigma, syndromes);
        cycles[2] = cpu_cycles() - clock;
    
        // Time roots search
        clock = cpu_cycles();
        uint8_t err[(1 << GF_M) / 8] = {0};
        compute_roots(err, sigma, deg_sigma);
        cycles[3] = cpu_cycles() - clock;
   
        // Time final calculations and clean up
        clock = cpu_cycles();
        vect_add(rcv, rcv, err, BCH_N_BYTES);
        message_from_codeword(msg, rcv);
        cycles[4] = cpu_cycles() - clock;

        /*** End of bch_decode function ***/

  
        // Print results
        printf("\t\t\tCPU cycles\n\n");
        printf("Precomputation step:\t%9" PRIu64 "\n", cycles[0]);
        printf("Syndrome computation:\t%9" PRIu64 "\n", cycles[1]);
        printf("elp computation:\t%9" PRIu64 "\n", cycles[2]);
        printf("Roots search:\t\t%9" PRIu64 "\n", cycles[3]);
        printf("Final calculations:\t%9" PRIu64 "\n\n", cycles[4]);
        printf("Total:\t\t\t%9" PRIu64 "\n",
               cycles[0] + cycles[1] + cycles[2] + cycles[3] + cycles[4]);
}
