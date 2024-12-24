/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2015 Marco Costalba, Joona Kiiski, Tord Romstad
  Copyright (C) 2015-2016 Marco Costalba, Joona Kiiski, Gary Linscott, Tord Romstad

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MISC_H
#define MISC_H


#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>

#include "types.h"


// prefetch() preloads the given address in L1/L2 cache. This is
// a non-blocking function that doesn't stall the CPU waiting for data
// to be loaded from memory, which can be quite slow.

SMALL static void prefetch(void* addr) { __builtin_prefetch(addr); }

typedef int64_t TimePoint;  // A value in milliseconds

static TimePoint now(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return 1000 * (uint64_t) tv.tv_sec + (uint64_t) tv.tv_usec / 1000;
}

struct PRNG {
    uint64_t s;
};

typedef struct PRNG PRNG;

void     prng_init(PRNG* rng, uint64_t seed);
uint64_t prng_rand(PRNG* rng);
uint64_t prng_sparse_rand(PRNG* rng);

static uint64_t mul_hi64(uint64_t a, uint64_t b) {
#if defined(__GNUC__)
    __extension__ typedef unsigned __int128 uint128;
    return ((uint128) a * (uint128) b) >> 64;
#else
    uint64_t aL = (uint32_t) a, aH = a >> 32;
    uint64_t bL = (uint32_t) b, bH = b >> 32;
    uint64_t c1 = (aL * bL) >> 32;
    uint64_t c2 = aH * bL + c1;
    uint64_t c3 = aL * bH + (uint32_t) c2;
    return aH * bH + (c2 >> 32) + (c3 >> 32);
#endif
}

#endif
