/**
 * @file gf.h
 * Header file of gf_lutmul.c and gf_pclmul.c
 */

#ifndef GF_H
#define GF_H

#include <stddef.h>
#include <stdint.h>

uint16_t
gf_primitive_poly(unsigned int m);

void
gf_generate(uint16_t            *exp,
            uint16_t            *log,
            unsigned int         m);

void
gf_generate_M();

uint16_t
gf_exp(uint16_t i);

uint16_t
gf_log(uint16_t elt);

uint16_t
gf_reduce(uint64_t      x,
          size_t        deg_x);

uint16_t
gf_mul(uint16_t a,
       uint16_t b);

uint16_t
gf_square(uint16_t a);

uint16_t
gf_quad(uint64_t a);

uint16_t
gf_inverse(uint16_t a);

uint16_t
gf_mod(uint16_t i);

#endif // GF_H
