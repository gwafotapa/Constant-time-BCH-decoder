/**
 * @file time.c
 * Measures execution time of BCH decoding
 */

#include "bch.h"

#include <fcntl.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>

#define PIECES          6       /**< pieces of codes to measure */
#define IGNORE          10      /**< default number of decodings ignored as warmup */
#define DECODINGS       10      /**< default number of decodings */
#define REPEAT          10      /**< default repeat amount for each decoding */
#define REPEAT_MAX      1000    /**< maximum value of REPEAT constant */
#define TOLERANCE       0.01    /**< default distance to mean (in percent) above which times are printed in red */
#define WEIGHTS         14      /**< default weight distribution size */
#define WEIGHTS_MAX     20      /**< maximum value of WEIGHTS constant */

#define BLUE    "\x1b[1;34m"
#define CYAN    "\x1b[1;36m"
#define YELLOW  "\x1b[1;33m"
#define GREEN   "\x1b[1;32m"
#define RED     "\x1b[1;31m"
#define RESET   "\x1b[0m"

#define BINARY  (strrchr(*argv, '/') + 1)

#define EXIT(message, ...)                                              \
        do {                                                            \
                fprintf(stderr, "error: " message ".\n", __VA_ARGS__);  \
                exit(1);                                                \
        } while (0)                                                     \

static int urandom_fd; /**< File descriptor of '/dev/urandom' */

/**
 * Generates a BCH_N-bit random vector of the supplied weight.
 * @param[out] vec Array of size BCH_N_BYTES
 * @param[in] weight Weight of the output vector vec
 */
void
vec_weighted(uint8_t    *vec,
             uint16_t    weight)
{
        unsigned long random_bytes_size = 3 * weight;
        unsigned char rand_bytes[random_bytes_size];
        memset(rand_bytes, 0, random_bytes_size);
        uint32_t support[weight];

        int res = read(urandom_fd, rand_bytes, random_bytes_size);
        if (res == -1)
                EXIT("%s", "failed to read from urandom");
  
        unsigned long j = 0;
        uint32_t random_index;
  
        for (size_t i = 0; i < weight; ++i) {
                do {
                        if (j == random_bytes_size) {
                                res = read(urandom_fd, rand_bytes, random_bytes_size);
                                if (res == -1)
                                        EXIT("%s", "failed to read from urandom");
                                j = 0;
                        }

                        random_index  = (uint32_t)rand_bytes[j++] << 16;
                        random_index |= (uint32_t)rand_bytes[j++] << 8;
                        random_index |= rand_bytes[j++];

                } while (random_index >= 16777216 - 16777216 % BCH_N);

                random_index %= BCH_N;

                bool new_index = true;
                for (size_t k = 0; k < i; ++k)
                        if (support[k] == random_index) {
                                new_index = false;
                                --i;
                                break;
                        }

                if (new_index)
                        support[i] = random_index;
        }

        for (size_t i = 0; i < weight; ++i) {
                size_t byte = support[i] / 8;
                size_t pos = support[i] % 8;
                vec[byte] |= 1 << pos;
        }
}

/**
 * Converts a byte string into an hexadecimal one.
 * @param[out] hexstr Array of size 2*BCH_N_BYTES+1 receiving the conversion
 * @param[in] uint8str Array of size BCH_N_BYTES to be converted
 */
void
uint8str_to_hexstr(char                 *hexstr,
                   uint8_t const        *uint8str)
{
        for (size_t i = 0; i < BCH_N_BYTES; ++i)
                sprintf(hexstr + 2*i, "%02x", uint8str[i]);
        hexstr[2 * BCH_N_BYTES] = '\0';
}

/**
 * Displays the measured execution times.
 * @param[in] str Name of the measured piece of code
 * @param[in] min Array of minimum execution times
 * @param[in] mean Array of mean execution times
 * @param[in] max Array of maximum execution times
 * @param[in] index Line number holding the times in the min, mean and max arrays
 * @param[in] tolerance Times whose distance to mean time (in percent) is greater than tolerance are printed in red
 * @param[in] weights Number of columns in arrays min, mean and max
 */
void
display(const char      *str,
        const uint64_t  *minima,
        const uint64_t  *means,
        const uint64_t  *maxima,
        size_t           index,
        float            tolerance,
        size_t           weights)
{
        printf("%s min:\t", str);
        for (size_t i = 0; i < weights+1; ++i) {
                if (i == weights) {
                        if (minima[(WEIGHTS_MAX+1) * index + i] < (1 - tolerance) * means[(WEIGHTS_MAX+1) * index + i])
                                printf(RED);
                        else 
                                printf(GREEN);
                }
                printf(" %9" PRIu64 RESET, minima[(WEIGHTS_MAX+1) * index + i]);
        }
  
        printf("\n%s mean:\t", str);
        for (size_t i = 0; i < weights; ++i)
                printf(" %9" PRIu64, means[(WEIGHTS_MAX+1) * index + i]);
        printf(BLUE " %9" PRIu64 "\n" RESET, means[(WEIGHTS_MAX+1) * index + weights]);
  
        printf("%s max:\t", str);
        for (size_t i = 0; i < weights+1; ++i) {
                if (i == weights) {
                        if (maxima[(WEIGHTS_MAX+1) * index + i] > (1 + tolerance) * means[(WEIGHTS_MAX+1) * index + i])
                                printf(RED);
                        else
                                printf(GREEN);
                }
                printf(" %9" PRIu64 RESET, maxima[(WEIGHTS_MAX+1) * index + i]);
        }
        printf("\n\n");
}

/**
 * Generates a random codeword with errors.
 * @param[out] cdw Array of size BCH_N_BYTES receiving the erroneous codeword
 * @param[in] weight Number of errors introduced in the codeword
 */
void
erroneous_codeword(uint8_t      *cdw,
                   size_t        weight)
{
        memset(cdw, 0, BCH_N_BYTES);
  
        // Generate a random plaintext and encode it
        uint8_t m[BCH_K_BYTES] = {0};
        int res = read(urandom_fd, m, BCH_K_BYTES);
        if (res == -1)
                EXIT("%s", "failed to read from urandom");
        bch_encode(cdw, bch_poly, m);

        uint8_t err[BCH_N_BYTES] = {0};
        vec_weighted(err, weight);
        vect_add(cdw, cdw, err, BCH_N_BYTES);
}

/**
 * Decodes the word rcv repeats times recording minimum execution times for each piece of code in array cycles.
 * @param[out] cycles Array receiving the minimum execution times
 * @param[in] rcv Word to decode
 * @param[in] repeats Amount of decodings to perform
 * @param[in] binary Name of this binary
 */
void
repeat_decoding(uint64_t        *cycles,
                const uint8_t   *rcv,
                size_t           repeats,
                const char      *binary)
{
        // Build command (message is passed to binary as hex string)
        char rcv_hex[2 * BCH_N_BYTES + 1];
        uint8str_to_hexstr(rcv_hex, rcv);
        char cmd[100 + 2 * BCH_N_BYTES + 1];
        sprintf(cmd, "%s_once %s", binary, rcv_hex);
      
        for (size_t repeat = 0; repeat < repeats; ++repeat) {
                char line[50];
                FILE *output = popen(cmd, "r");
                size_t j = 0;
                while (fgets(line, 50, output)) {
                        char *ptr = strchr(line, ':');
                        if (ptr != NULL) {
                                uint64_t measure = strtoull(++ptr, NULL, 10);
                                if ((cycles[j] == 0) || (measure < cycles[j]))
                                        cycles[j] = measure;
                                ++j;
                        }
                }
                pclose(output);
        }
}

/**
 * Updates columns col of 2D-arrays cycles_minima, cycles_means and cycles_maxima with values from array cycles.
 * @param[in,out] cycles_minima Array holding minimum execution times
 * @param[in,out] cycles_means Array holding mean execution times
 * @param[in,out] cycles_maxima Array holding maximum execution times
 * @param[in] cycles Array holding decoding times of a word
 * @param[in] repeats Number of execution times held in array cycles
 * @param[in] col Decoded error weight of decodings whose times is held in array cycles
 */
void
update_column(uint64_t          *cycles_minima,
              uint64_t          *cycles_means,
              uint64_t          *cycles_maxima,
              const uint64_t    *cycles,
              size_t             repeats,
              size_t             col)
{
        for (size_t j = 0; j < PIECES; ++j) {
                size_t index = j*(WEIGHTS_MAX+1) + col;
                cycles_means[index] += cycles[j];
                if ((cycles_minima[index] == 0) || (cycles[j] < cycles_minima[index]))
                        cycles_minima[index] = cycles[j];
                if (cycles[j] > cycles_maxima[index])
                        cycles_maxima[index] = cycles[j];
        }
}

/**
 * Compares the decoding time of word with minimum and maximum and updates them as needed.
 * @param[in,out] word_minimum Word for which cycles_minimum was obtained
 * @param[in,out] word_maximum Word for which cycles_maximum was obtained
 * @param[in,out] cycles_minimum Global minimum decoding time
 * @param[in,out] cycles_maximum Global maximum decoding time
 * @prarm[in] word Word for which cycles_word was obtained
 * @param[in] cycles_word Decoding time of word
 */
void
update_extremums(uint8_t        *word_minimum,
                 uint8_t        *word_maximum,
                 uint64_t       *cycles_minimum,
                 uint64_t       *cycles_maximum,
                 const uint8_t  *word,
                 uint64_t        cycles_word)
{
        if ((*cycles_minimum == 0) || (cycles_word < *cycles_minimum)) {
                *cycles_minimum = cycles_word;
                memcpy(word_minimum, word, BCH_N_BYTES);
        }
        if (cycles_word > *cycles_maximum) {
                *cycles_maximum = cycles_word;
                memcpy(word_maximum, word, BCH_N_BYTES);
        }
}

/**
 * Returns a random non negative integer less than the supplied bound.
 * @returns an integer less than upper_bound
 * @param[in] upper_bound Upper bound
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

int
main(int          argc,
     char       **argv)
{
        struct option long_options[] = {
                { "decodings", required_argument, NULL, 'd' },
                { "help", no_argument, NULL, 'h' },
                { "ignore", required_argument, NULL, 'i' },
                { "random", no_argument, NULL, 'R' },
                { "repeat", required_argument, NULL, 'r' },
                { "tolerance", required_argument, NULL, 't' },
                { "totals", no_argument, NULL, 'T' },
                { "weights", required_argument, NULL, 'w' },
                { NULL, 0, NULL, 0 }
        };

        bool total_times_only = false;
        bool random = false;
        size_t ignore = IGNORE;
        size_t repeats = REPEAT;
        size_t decodings = DECODINGS;
        float tolerance = TOLERANCE;
        size_t weights_count = WEIGHTS;
        char opt;
  
        while ((opt = getopt_long(argc, argv, "d:hi:r:Rt:Tw:", long_options, NULL)) != -1) {
                switch (opt) {
                case 'd':
                        decodings = strtoul(optarg, NULL, 10);
                        if (decodings == 0)
                                EXIT("%s", "decodings count cannot be 0");
                        break;
        
                case 'h':
                        printf("Usage: %s[OPTION]...\n"
                               "Measures running time of the BCH code decoder for code (%u,%u,%u).\n"
                               "First does and ignores IGNORE decodings as warmup.\n"
                               "Then does and measures DECODINGS decodings for each different weight.\n"
                               "Each decoding is done REPEAT times and the minimum execution time is taken as estimate.\n"
                               "\n"
                               "Mandatory arguments to long options are mandatory for short options too.\n"
                               "  -d, --decodings=DECODINGS\tnumber of decodings to run (default: %u)\n"
                               "  -h, --help\t\t\tdisplay this help\n"
                               "  -i, --ignore=IGNORE\t\tignore first IGNORE measures (default: %u)\n"
                               "  -R, --random\t\t\tdecode random words\n"
                               "  -r, --repeat=REPEAT\t\trepeat count for minimum computation (default: %u)\n"
                               "  -t, --tolerance=TOLERANCE\tredflag value for min and max (default: %.2f)\n"
                               "  -T, --totals\t\t\tdisplay total times only\n"
                               "  -w, --weights=WEIGHTS\t\tnumber of different weights (default: %u)\n",
                               BINARY, BCH_N, BCH_K, BCH_DELTA, DECODINGS, IGNORE, REPEAT, TOLERANCE, WEIGHTS);
                        exit(0);

                case 'i':
                        ignore = strtoul(optarg, NULL, 10);
                        break;

                case 'r':
                        repeats = strtoul(optarg, NULL, 10);
                        if (repeats == 0)
                                EXIT("%s", "repeat count cannot be 0");
                        if (repeats > REPEAT_MAX)
                                EXIT("repeat count cannot exceed %u", REPEAT_MAX);
                        break;

                case 'R':
                        random = true;
                        break;
        
                case 't':
                        tolerance = strtof(optarg, NULL);
                        break;
        
                case 'T':
                        total_times_only = true;
                        break;

                case 'w':
                        weights_count = strtoul(optarg, NULL, 10);
                        if (weights_count < 3)
                                EXIT("%s", "weights cannot be less than 3");
                        if (weights_count > WEIGHTS_MAX)
                                EXIT("weights cannot exceed %u", WEIGHTS_MAX);
                        break;
        
                default:
                        exit(1);
                }
        }

        // Build the array of weights.
        // Spread the weights from 0 to BCH_DELTA and add a last one greater than BCH_DELTA.
        size_t weights[WEIGHTS_MAX];
        weights[0] = 0;
        weights[1] = 1;
        for (size_t i = 2; i < weights_count-1; ++i)
                weights[i] = weights[i-1] + (BCH_DELTA-weights[i-1]) / (weights_count-2-i+1);
        weights[weights_count-1] = BCH_DELTA + (weights[weights_count-2] - weights[weights_count-3]);

        // Minimums, means and maximums running times of the different pieces of code for each weight
        uint64_t cycles_minima[PIECES][WEIGHTS_MAX + 1] = {0};  // WEIGHTS_MAX + 1 because of the total column
        uint64_t cycles_means[PIECES][WEIGHTS_MAX + 1] = {0};
        uint64_t cycles_maxima[PIECES][WEIGHTS_MAX + 1] = {0};

        // Global minimum and maximum execution times and the words giving them
        uint64_t cycles_minimum = 0;
        uint64_t cycles_maximum = 0;
        uint8_t word_minimum[BCH_N_BYTES];
        uint8_t word_maximum[BCH_N_BYTES];
  
        urandom_fd = open("/dev/urandom", O_RDONLY);
        if (urandom_fd == -1)
                EXIT("%s", "cannot open /dev/urandom");

        // Display test intro message
        printf(BLUE "Decodings of BCH code (%u,%u,%u)\n" RESET, BCH_N, BCH_K, BCH_DELTA);
        printf("Warmup decodings:    %7zu\n", ignore);
        if (random)
                printf("Measured decodings:  %7zu\n", decodings);
        else
                printf("Decodings per error: %7zu\n", decodings);
        printf("Repeat per decoding: %7zu\n", repeats);
        printf(YELLOW "Running test"), fflush(stdout);

        size_t cnt = 0; // number of decodings done
  
        for (size_t try = 0; try < ignore; ++try, cnt+=repeats) {
                // Generate a random error
                uint16_t weight = random_int_less_than(2 * BCH_DELTA);

                // Generate a random codeword with weight errors
                uint8_t word[BCH_N_BYTES];
                erroneous_codeword(word, weight);

                // Decode it ignore times
                uint64_t cycles_word[PIECES] = {0};
                repeat_decoding(cycles_word, word, repeats, *argv);
        }

        if (random) {
                weights_count = 0;
    
                for (size_t try = 0; try < decodings; ++try, cnt+=repeats) {      
                        // Server we run on will kick us if we're idle too long
                        if (cnt > 10000) {
                                printf("."), fflush(stdout);
                                cnt -= 10000;
                        }

                        // Generate a random error
                        uint16_t weight = random_int_less_than(2*BCH_DELTA);
      
                        // Generate a random codeword with weight errors
                        uint8_t word[BCH_N_BYTES];
                        erroneous_codeword(word, weight);
      
                        // Decode it repeat times
                        uint64_t cycles_word[PIECES] = {0};
                        repeat_decoding(cycles_word, word, repeats, *argv);

                        // Update min, mean and max
                        update_column((uint64_t *)cycles_minima,
                                      (uint64_t *)cycles_means,
                                      (uint64_t *)cycles_maxima,
                                      cycles_word,
                                      repeats,
                                      0);

                        update_extremums(word_minimum,
                                         word_maximum,
                                         &cycles_minimum,
                                         &cycles_maximum,
                                         word,
                                         cycles_word[PIECES-1]);
                }
    
                for (size_t j = 0; j < PIECES; ++j)
                        cycles_means[j][0] /= decodings;

                printf(RESET "\n");
        } else {
                for (size_t i = 0; i < weights_count; ++i) {
                        for (size_t try = 0; try < decodings; ++try, cnt+=repeats) {
                                // Server we run on will kick us if we're idle too long
                                if (cnt > 10000) {
                                        printf("."), fflush(stdout);
                                        cnt = 0;
                                }
      
                                // Generate a random codeword with weight errors
                                uint8_t word[BCH_N_BYTES];
                                erroneous_codeword(word, weights[i]);

                                // Decode it repeat times
                                uint64_t cycles_word[PIECES] = {0};
                                repeat_decoding(cycles_word, word, repeats, *argv);
                                
                                // Update min, mean and max from column i
                                update_column((uint64_t *)cycles_minima,
                                              (uint64_t *)cycles_means,
                                              (uint64_t *)cycles_maxima,
                                              cycles_word,
                                              repeats,
                                              i);
      
                                update_extremums(word_minimum,
                                                 word_maximum,
                                                 &cycles_minimum,
                                                 &cycles_maximum,
                                                 word,
                                                 cycles_word[PIECES-1]);
                        }
    
                        for (size_t j = 0; j < PIECES; ++j)
                                cycles_means[j][i] /= decodings;

                        // Update min, mean and max from last column 'All'
                        for (size_t j = 0; j < PIECES; ++j) {
                                cycles_means[j][weights_count] += cycles_means[j][i];
                                bool new_global_minimum = (cycles_minima[j][weights_count] == 0) ||
                                        (cycles_minima[j][i] < cycles_minima[j][weights_count]);
                                if (new_global_minimum)
                                        cycles_minima[j][weights_count] = cycles_minima[j][i];
                                bool new_global_maximum = cycles_maxima[j][i] > cycles_maxima[j][weights_count];
                                if (new_global_maximum)
                                        cycles_maxima[j][weights_count] = cycles_maxima[j][i];
                        }
                }

                for (size_t j = 0; j < PIECES; ++j)
                        cycles_means[j][weights_count] /= weights_count;
  
                printf(RESET "\n");
  
                printf("Errors:\t\t");
                for (size_t i = 0; i < weights_count; ++i)
                        printf(" %9zu", weights[i]);
                printf(CYAN " %9s\n\n" RESET, "All");
        }
  
        close(urandom_fd);

        if (!total_times_only) {
                display("LUT",
                        (uint64_t *)cycles_minima,
                        (uint64_t *)cycles_means,
                        (uint64_t *)cycles_maxima,
                        0, tolerance, weights_count);
                display("Syndromes",
                        (uint64_t *)cycles_minima,
                        (uint64_t *)cycles_means,
                        (uint64_t *)cycles_maxima,
                        1, tolerance, weights_count);
                display("ELP",
                        (uint64_t *)cycles_minima,
                        (uint64_t *)cycles_means,
                        (uint64_t *)cycles_maxima,
                        2, tolerance, weights_count);
                display("Roots",
                        (uint64_t *)cycles_minima,
                        (uint64_t *)cycles_means,
                        (uint64_t *)cycles_maxima,
                        3, tolerance, weights_count);
                display("Cleanup",
                        (uint64_t *)cycles_minima,
                        (uint64_t *)cycles_means,
                        (uint64_t *)cycles_maxima,
                        4, tolerance, weights_count);
                printf("\n");
        }
        display("Total",
                (uint64_t *)cycles_minima,
                (uint64_t *)cycles_means,
                (uint64_t *)cycles_maxima,
                5, tolerance, weights_count);

        // Rerun test for min and max

        uint64_t cycles[PIECES] = {0};
        repeat_decoding(cycles, word_minimum, repeats, *argv);
        printf("Rerun %zu decodings for total minimum: %" PRIu64 "\n", repeats, cycles[PIECES-1]);

        memset(cycles, 0, PIECES * sizeof(uint64_t));
        repeat_decoding(cycles, word_maximum, repeats, *argv);
        printf("Rerun %zu decodings for total maximum: %" PRIu64 "\n", repeats, cycles[PIECES-1]);
}
