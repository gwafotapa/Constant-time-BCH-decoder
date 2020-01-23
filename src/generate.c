/**
 * @file generate.c
 * Generates parameter files of BCH codes
 */

#include "bch.h"
#include "gf.h"

#include <getopt.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define BINARY (strrchr(*argv, '/') + 1)

#define EXIT(message, ...)                                              \
        do {                                                            \
                fprintf(stderr, "%s: " message ".\n", BINARY, __VA_ARGS__); \
                exit(1);                                                \
        } while (0)                                                     \
                
int
main (int       argc,
      char     *argv[])
{
        struct option long_options[] = {
                { "help", no_argument, NULL, 'h' },
                { "no-files", no_argument, NULL, 'n'},
                { "shorten", required_argument, NULL, 's' },
                { NULL, 0, NULL, 0 }
        };
                
        bool generate_header = true;
        size_t cut = 0;
        char opt;
  
        while ((opt = getopt_long(argc, argv, "hns:", long_options, NULL)) != -1) {
                switch (opt)
                        {
                        case 'h':
                                printf("Usage: %s M T\n"
                                       "Compute the largest dimension K of the primitive BCH code (2^M-1, K, DELTA >= T) "
                                       "and generate corresponding files.\n"
                                       "2 < M < 16 supported.\n"
                                       "\n"
                                       "Mandatory arguments to long options are mandatory for short options too.\n"
                                       "  -h, --help\t\t\tdisplay this help\n"
                                       "  -n, --no-files\t\tdisplay results but do not generate files\n"
                                       "  -s, --shorten=CUT\t\tshorten the code by discarding the first CUT coordinates\n",
                                       BINARY);
                                exit(0);

                        case 'n':
                                generate_header = false;
                                break;

                        case 's':
                                cut = strtoul(optarg, NULL, 10);
                                break;
                        
                        default:
                                exit(1);
                        }
        }
        
        if (argc - optind != 2)
                EXIT("%s", "two arguments expected.\n"
                     "Use '--help' for more info");
        
        size_t m = strtoul(argv[optind], NULL, 10);
        if ((m < 3) || (15 < m))
                EXIT("M = %s not supported. M must verify 2 < M < 16", argv[optind]);

        size_t bch_n = (1 << m) - 1; // length
        size_t bch_t = strtoul(argv[++optind], NULL, 10); // correction capacity
        if ((1 << m-1) <= bch_t)
                EXIT("%s", "T must be less than 2^(M-1)");

        size_t bch_poly_size = m * bch_t + 1; // upper bound on redundancy
        uint16_t *bch_poly = calloc(bch_poly_size, sizeof *bch_poly);
        uint16_t *exp = malloc(((1 << m) + 2) * sizeof *exp);
        uint16_t *log = malloc((1 << m) * sizeof *log);

        gf_generate(exp, log, m);
        size_t deg_bch_poly = compute_bch_poly(bch_poly, &bch_t, m, exp, log);

        size_t bch_k = bch_n - deg_bch_poly; // code dimension
        if (cut > bch_k-1)
                EXIT("cannot shorten that much: K = %zu", bch_k);
        bch_n -= cut;
        bch_k -= cut;

        printf("Code length:              %5zu\n", bch_n);
        printf("Code dimension:           %5zu\n", bch_k);
        printf("Code correction capacity: %5zu\n", bch_t);

        if (!generate_header)
                goto Cleanup;

        uint16_t gf_poly = gf_primitive_poly(m);
        size_t fft_param = 0; // log of the smallest power of 2 greater than correction capacity
        size_t tmp = bch_t;
        do ++fft_param; while (tmp >>=1);

        // Create filenames
        char *output_dirname = malloc(30 * sizeof *output_dirname);
        if (output_dirname == NULL)
                EXIT("%s", "cannot allocate memory for output directory name");
        sprintf(output_dirname, "bch_%zu_%zu_%zu", bch_n, bch_k, bch_t);
                
        char *bch_parameters_filename = malloc(30 * sizeof *bch_parameters_filename);
        if (bch_parameters_filename == NULL)
                EXIT("%s", "cannot allocate memory for the BCH parameters filename");
        sprintf(bch_parameters_filename, "bch_%zu_%zu_%zu.h", bch_n, bch_k, bch_t);
                
        char *bch_poly_filename = malloc(30 * sizeof *bch_poly_filename);
        if (bch_poly_filename == NULL)
                EXIT("%s", "cannot allocate memory for the generator polynomial filename");
        sprintf(bch_poly_filename, "bch_%zu_%zu_%zu_poly.h", bch_n, bch_k, bch_t);

        // Create and enter directory
        int error_status = mkdir(output_dirname, 0700);
        if (error_status == -1)
                EXIT("cannot create directory %s", output_dirname);
        error_status = chdir(output_dirname);
        if (error_status == -1)
                EXIT("cannot enter directory %s", output_dirname);

        // Create bch parameters file
        FILE *bch_parameters_file = fopen(bch_parameters_filename, "w");
        if (bch_parameters_file == NULL)
                EXIT("cannot open file %s", bch_parameters_filename);
        fprintf(bch_parameters_file, "/**\n");
        fprintf(bch_parameters_file, " * @file bch_%zu_%zu_%zu.h\n", bch_n, bch_k, bch_t);
        fprintf(bch_parameters_file, " * Parameters of the BCH code (%zu, %zu, %zu)\n", bch_n, bch_k, bch_t);
        fprintf(bch_parameters_file, " */\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "#ifndef BCH_%zu_%zu_%zu_H\n", bch_n, bch_k, bch_t);
        fprintf(bch_parameters_file, "#define BCH_%zu_%zu_%zu_H\n", bch_n, bch_k, bch_t);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Parameter of Galois field GF(2^GF_M) */\n");
        fprintf(bch_parameters_file, "#define GF_M             %zu\n", m);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Order of the Galois field multiplicative group */\n");
        fprintf(bch_parameters_file, "#define GF_MUL_ORDER     ((1 << GF_M) - 1)\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Primitive polynomial the Galois field was built with */\n");
        fprintf(bch_parameters_file, "#define GF_POLY          %#X\n", gf_poly);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Length of the code */\n");
        fprintf(bch_parameters_file, "#define BCH_N            %zu\n", bch_n);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Dimension of the code */\n");
        fprintf(bch_parameters_file, "#define BCH_K            %zu\n", bch_k);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Correction capacity of the code */\n");
        fprintf(bch_parameters_file, "#define BCH_DELTA        %zu\n", bch_t);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Number of bytes needed to store BCH_N bits */\n");
        fprintf(bch_parameters_file, "#define BCH_N_BYTES      CEIL_DIV(BCH_N, 8)\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Number of bytes needed to store BCH_K bits */\n");
        fprintf(bch_parameters_file, "#define BCH_K_BYTES      CEIL_DIV(BCH_K, 8)\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Number of coefficients of the generator polynomial of the BCH code */\n");
        fprintf(bch_parameters_file, "#define BCH_POLY_SIZE    (BCH_N - BCH_K + 1)\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Logarithm of the smallest power of 2 greater than or equal to BCH_DELTA */\n");
        fprintf(bch_parameters_file, "#define FFT_PARAM        %zu\n", fft_param);
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Logarithm of the smallest power of 2 greater than or equal to the number of syndromes */\n");
        fprintf(bch_parameters_file, "#define FFT_T_PARAM      (FFT_PARAM + 1)\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "/** Computes the ceiling of the division of a by b */\n");
        fprintf(bch_parameters_file, "#define CEIL_DIV(a, b)   ((a)/(b) + ((a)%%(b) == 0 ? 0 : 1))\n");
        fprintf(bch_parameters_file, "\n");
        fprintf(bch_parameters_file, "#endif // BCH_%zu_%zu_%zu_H\n", bch_n, bch_k, bch_t);
        fclose(bch_parameters_file);

        // Create generator polynomial file
        FILE *bch_poly_file = fopen(bch_poly_filename, "w");
        if (bch_poly_file == NULL)
                EXIT("cannot open file %s", bch_poly_filename);
        fprintf(bch_poly_file, "/**\n");
        fprintf(bch_poly_file, " * @file bch_%zu_%zu_%zu_poly.h\n", bch_n, bch_k, bch_t);
        fprintf(bch_poly_file, " * Generator polynomial of the BCH code (%zu, %zu, %zu)\n", bch_n, bch_k, bch_t);
        fprintf(bch_poly_file, " */\n");
        fprintf(bch_poly_file, "\n");
        fprintf(bch_poly_file, "#ifndef BCH_%zu_%zu_%zu_POLY_H\n", bch_n, bch_k, bch_t);
        fprintf(bch_poly_file, "#define BCH_%zu_%zu_%zu_POLY_H\n", bch_n, bch_k, bch_t);
        fprintf(bch_poly_file, "\n");
        fprintf(bch_poly_file, "#include <stdint.h>\n");
        fprintf(bch_poly_file, "\n");
        fprintf(bch_poly_file, "/** Coefficients of the generator polynomial of the code in ascending order */\n");
        fprintf(bch_poly_file, "static const uint8_t bch_poly[%zu] = { ", deg_bch_poly + 1);
        for (size_t i = 0; i < deg_bch_poly; ++i)
                fprintf(bch_poly_file, "%u,", bch_poly[i]);
        fprintf(bch_poly_file, "%u };\n", bch_poly[deg_bch_poly]);      
        fprintf(bch_poly_file, "\n");
        fprintf(bch_poly_file, "#endif // BCH_%zu_%zu_%zu_POLY_H\n", bch_n, bch_k, bch_t);
        fclose(bch_poly_file);

        // Display closing message and clean up
        printf("BCH files generated in directory %s.\n", output_dirname);
        printf("Run 'make %s BCHDIR=%s' to use it.\n",
               (m <= 12) ? "lutmul" : "pclmul",
               output_dirname);
        free(output_dirname);
        free(bch_parameters_filename);
        free(bch_poly_filename);
 Cleanup:
        free(bch_poly);
        free(exp);
        free(log);
}
