/**
 * @file bch.h
 * Header file of bch.c
 */

#ifndef BCH_H
#define BCH_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

uint16_t
mod(uint16_t    i,
    uint16_t    modulus);

void
compute_cyclotomic_cosets(uint16_t      *cosets,
                          uint16_t       upper_bound,
                          size_t         m);

size_t
compute_bch_poly(uint16_t       *bch_poly,
                 size_t         *t,
                 size_t          m,
                 const uint16_t *exp,
                 const uint16_t *log);

void
unpack_message(uint8_t          *msg_unpacked,
               const uint8_t    *msg);

void
lfsr_encode(uint8_t             *cdw,
            const uint8_t       *bch_poly,
            const uint8_t       *msg);

void
pack_codeword(uint8_t           *cdw,
              const uint8_t     *cdw_unpacked);

void
bch_encode(uint8_t              *cdw,
           const uint8_t        *bch_poly,
           const uint8_t        *msg);

size_t
compute_elp(uint16_t            *sigma,
            const uint16_t      *syndromes);

void
message_from_codeword(uint8_t           *msg,
                      const uint8_t     *cdw);

void
compute_syndromes(uint16_t      *syndromes,
                  const uint8_t *rcv);

void
compute_roots(uint8_t   *err,
              uint16_t  *sigma,
              size_t     deg_sigma);

void
vect_add(uint8_t        *sum,
         const uint8_t  *v1,
         const uint8_t  *v2,
         size_t          size);

void
bch_decode(uint8_t      *msg,
           uint8_t      *rcv);

#endif // BCH_H
