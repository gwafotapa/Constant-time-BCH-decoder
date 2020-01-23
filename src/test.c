/**
 * @file test.c
 * Simple test of BCH encoding and decoding
 */

#include "bch.h"
#include "fft.h"
#include "gf.h"
#include "optimizations.h"

#include <fcntl.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BLUE    "\x1b[1;34m"
#define CYAN    "\x1b[1;36m"
#define GREEN   "\x1b[1;32m"
#define RED     "\x1b[1;31m"
#define RESET   "\x1b[0m"

#define EXIT(message, ...)                                              \
        do {                                                            \
                fprintf(stderr, "error: " message ".\n", __VA_ARGS__);	\
                exit(1);                                                \
        } while (0)                                                     \

static int urandom_fd; /**< File descriptor of '/dev/urandom' */

/**
 * Generates a random non-negative integer less than upper_bound.
 * @returns the randomly generated integer
 * @param upper_bound Upper bound on the generated integer
 */
uint16_t
random_int_less_than(uint16_t upper_bound)
{
        uint16_t weight = 0;
        do {
                int res = read(urandom_fd, &weight, 2);
                if (res == -1)
                        EXIT("%s", "failed to read from urandom");
        } while (weight >= (1UL<<16) - (1UL<<16)%upper_bound);
        weight %= upper_bound;

        return weight;
}

/**
 * Generates a random vector of the given weight.
 * @param[out] vec Array receiving the generated vector
 * @param[in] size Size of the generated vector
 * @param[in] weight Weight of the generated vector
 */
void
vec_weighted(uint8_t	*vec,
	     size_t	 size,
	     size_t	 weight)
{
        size_t rand_bytes_size = 3 * weight;
        unsigned char rand_bytes[rand_bytes_size];
        memset(rand_bytes, 0, rand_bytes_size);
        uint32_t support[weight];

        int res = read(urandom_fd, rand_bytes, rand_bytes_size);
        if (res == -1)
                EXIT("%s", "failed to read from urandom");
  
        unsigned long j = 0;
        uint32_t rand_index;
  
        for (size_t i = 0; i < weight; ++i) {
                do {
                        if (j == rand_bytes_size) {
                                res = read(urandom_fd, rand_bytes, rand_bytes_size);
                                if (res == -1)
                                        EXIT("%s", "failed to read from urandom");
                                j = 0;
                        }

                        rand_index  = (uint32_t)rand_bytes[j++] << 16;
                        rand_index |= (uint32_t)rand_bytes[j++] << 8;
                        rand_index |= rand_bytes[j++];

                } while (rand_index >= (1<<24) - (1<<24)%size);

                rand_index %= size;

                bool new_index = true;
                for (size_t k = 0; k < i; ++k) {
                        if (support[k] == rand_index) {
                                new_index = false;
                                --i;
                                break;
                        }
                }

                if (new_index)
                        support[i] = rand_index;
        }

        for (size_t i = 0; i < weight; ++i) {
                size_t byte = support[i] / 8;
                size_t pos = support[i] % 8;
                vec[byte] |= 1 << pos;
        }
}

/** Array reversing the bits of the given byte */
static const uint8_t reverse_bits[256] = {
	0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
	0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
	0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
	0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
	0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
	0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
	0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
	0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
	0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
	0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
	0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
	0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
	0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
	0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
	0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
	0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

/**
 * Prints the given vector to the given stream.
 * @param[out] stream Stream the vector is output to
 * @param[in] vec Vector to print
 * @param[in] size Size of the vector vec
 */
void
vec_print(FILE		*stream,
	  const uint8_t	*vec,
	  size_t	 size)
{
	for (size_t i = 0; i < size; ++i)
		fprintf(stream, "%02X", reverse_bits[vec[i]]);
	printf("\n");
}

int
main(int        argc,
     char      *argv[])
{
        int decodings = 1;
        if (argc == 2)
                decodings = strtol(argv[1], NULL, 10);

        urandom_fd = open("/dev/urandom", O_RDONLY);
        if (urandom_fd == -1)
                EXIT("%s", "cannot open /dev/urandom");

	printf(CYAN "Encoding and decoding of BCH code (%u,%u,%u)\n" RESET,
	       BCH_N, BCH_K, BCH_DELTA);
	
        for (size_t i = 0; i < decodings; ++i) {
		if (decodings > 1)
			printf(CYAN "\n\n### Test %zu ###\n" RESET, i);
    
                // Generate a random plaintext
                uint8_t msg1[BCH_K_BYTES] = {0};
                int res = read(urandom_fd, msg1, BCH_K_BYTES);
                if (res == -1)
                        EXIT("%s", "failed to read from urandom");
		if (BCH_K % 8 != 0)
			msg1[BCH_K_BYTES-1] &= (1 << BCH_K%8) - 1;
		printf(BLUE "\nRandom message: " RESET);
		vec_print(stdout, msg1, BCH_K_BYTES);

		// Encode it
                uint8_t rcv[BCH_N_BYTES] = {0};
                bch_encode(rcv, bch_poly, msg1);
		printf(BLUE "\nEncoded message: " RESET);
		vec_print(stdout, rcv, BCH_N_BYTES);
		
                // Add a random correctible error
                uint16_t weight = random_int_less_than(BCH_DELTA + 1);
		uint8_t err1[BCH_N_BYTES] = {0};
                vec_weighted(err1, BCH_N, weight);
                vect_add(rcv, rcv, err1, BCH_N_BYTES);
		printf(BLUE "\nRandom error: " RESET);
		vec_print(stdout, err1, BCH_N_BYTES);
		printf(BLUE "\nReceived word: " RESET);
		vec_print(stdout, rcv, BCH_N_BYTES);

#ifndef GF_LUT_1024
		gf_generate_M();
#endif

#ifdef COMPUTE_FFT_LUT
		compute_fft_lut();
#endif
	
		// Calculate the 2 * BCH_DELTA syndromes
		uint16_t syndromes[1 << FFT_T_PARAM];
		compute_syndromes(syndromes, rcv);
		printf(BLUE "\nSyndromes: " RESET);
		for (size_t i = 0; i < 2*BCH_DELTA; ++i)
			printf("%u ", syndromes[i]);
		
		// Compute the error-locator polynomial (ELP)
		uint16_t sigma[1 << FFT_PARAM] = {0};
		size_t deg_sigma = compute_elp(sigma, syndromes);
		printf(BLUE "\n\nError-locator polynomial: " RESET);
		bool first_coeff = true;
		if (sigma[0]) {
			printf ("%u", sigma[0]);
			first_coeff = false;
		}
		for (size_t i = 1; i <= deg_sigma; ++i) {
			if (sigma[i] == 0)
				continue;
			if (!first_coeff)
				printf(" + ");
			first_coeff = false;
			if (sigma[i] != 1)
				printf("%u", sigma[i]);
			if (i == 1)
				printf("x");
			else
				printf("x^%zu", i);
		}
		if (first_coeff)
			printf("0");
		
		// Compute the error polynomial
		uint8_t err2[(1 << GF_M) / 8] = {0};
		compute_roots(err2, sigma, deg_sigma);
		printf(BLUE "\n\nError-locator numbers: " RESET);
		for (size_t i = 0; i < BCH_N; ++i)
			if (err2[i/8] & (1 << (i%8)))
				printf("%zu ", i);
		printf(BLUE "\n\nDecoded error: " RESET);
		vec_print(stdout, err2, BCH_N_BYTES);
		
		// Add the error polynomial to the received polynomial 
		vect_add(rcv, rcv, err2, BCH_N_BYTES);
		printf(BLUE "\nDecoded codeword: " RESET);
		vec_print(stdout, rcv, BCH_N_BYTES);
		
		// Retrieve the message from the decoded codeword and check the result
		uint8_t msg2[BCH_K_BYTES] = {0};
		message_from_codeword(msg2, rcv);
		printf(BLUE "\nDecoded message: " RESET);
                if (memcmp(msg1, msg2, BCH_K_BYTES-(BCH_K%8!=0)) ||
                    ((msg1[BCH_K_BYTES-1] ^ msg2[BCH_K_BYTES-1]) & ((uint8_t)0xff >> 8-BCH_K%8))) {
			printf(RED);
			vec_print(stdout, msg2, BCH_K_BYTES);
			printf(RESET);
			exit(1);
		} else {
			printf(GREEN);
			vec_print(stdout, msg2, BCH_K_BYTES);
			printf(RESET);
		}
        }

        close(urandom_fd);
}
