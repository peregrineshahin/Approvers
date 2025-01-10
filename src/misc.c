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

#include "misc.h"

#ifndef KAGGLE
    #include <stdio.h>
    #include <stdint.h>
    #include <stdatomic.h>
    #include <math.h>
#endif

// xorshift64star Pseudo-Random Number Generator
// This class is based on original code written and dedicated
// to the public domain by Sebastiano Vigna (2014).
// It has the following characteristics:
//
//  -  Outputs 64-bit numbers
//  -  Passes Dieharder and SmallCrush test batteries
//  -  Does not require warm-up, no zeroland to escape
//  -  Internal state is a single 64-bit integer
//  -  Period is 2^64 - 1
//  -  Speed: 1.60 ns/call (Core i7 @3.40GHz)
//
// For further analysis see
//   <http://vigna.di.unimi.it/ftp/papers/xorshift.pdf>

void prng_init(PRNG* rng, uint64_t seed) { rng->s = seed; }

uint64_t prng_rand(PRNG* rng) {
    uint64_t s = rng->s;

    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    rng->s = s;

    return s * 2685821657736338717LL;
}

uint64_t prng_sparse_rand(PRNG* rng) {
    uint64_t r1 = prng_rand(rng);
    uint64_t r2 = prng_rand(rng);
    uint64_t r3 = prng_rand(rng);
    return r1 & r2 & r3;
}

#ifndef KAGGLE

typedef struct {
    atomic_int_fast64_t data[6];
} DebugInfo;

    #define MAX_DEBUG_SLOTS 32

DebugInfo hit[MAX_DEBUG_SLOTS];
DebugInfo mean[MAX_DEBUG_SLOTS];
DebugInfo stdev[MAX_DEBUG_SLOTS];
DebugInfo correl[MAX_DEBUG_SLOTS];

void dbg_hit_on(int cond, int slot) {
    atomic_fetch_add(&hit[slot].data[0], 1);
    if (cond)
    {
        atomic_fetch_add(&hit[slot].data[1], 1);
    }
}

void dbg_mean_of(int64_t value, int slot) {
    atomic_fetch_add(&mean[slot].data[0], 1);
    atomic_fetch_add(&mean[slot].data[1], value);
}

void dbg_stdev_of(int64_t value, int slot) {
    atomic_fetch_add(&stdev[slot].data[0], 1);
    atomic_fetch_add(&stdev[slot].data[1], value);
    atomic_fetch_add(&stdev[slot].data[2], value * value);
}

void dbg_correl_of(int64_t value1, int64_t value2, int slot) {
    atomic_fetch_add(&correl[slot].data[0], 1);
    atomic_fetch_add(&correl[slot].data[1], value1);
    atomic_fetch_add(&correl[slot].data[2], value1 * value1);
    atomic_fetch_add(&correl[slot].data[3], value2);
    atomic_fetch_add(&correl[slot].data[4], value2 * value2);
    atomic_fetch_add(&correl[slot].data[5], value1 * value2);
}

void dbg_print() {
    int64_t n;
    for (int i = 0; i < MAX_DEBUG_SLOTS; ++i)
    {
        if ((n = atomic_load(&hit[i].data[0])))
        {
            double hitRate = 100.0 * (double) atomic_load(&hit[i].data[1]) / n;
            fprintf(stderr, "Hit #%d: Total %lld Hits %lld Hit Rate (%%) %.2f\n", i, n,
                    atomic_load(&hit[i].data[1]), hitRate);
        }
    }

    for (int i = 0; i < MAX_DEBUG_SLOTS; ++i)
    {
        if ((n = atomic_load(&mean[i].data[0])))
        {
            double meanValue = (double) atomic_load(&mean[i].data[1]) / n;
            fprintf(stderr, "Mean #%d: Total %lld Mean %.2f\n", i, n, meanValue);
        }
    }

    for (int i = 0; i < MAX_DEBUG_SLOTS; ++i)
    {
        if ((n = atomic_load(&stdev[i].data[0])))
        {
            double meanValue = (double) atomic_load(&stdev[i].data[1]) / n;
            double variance =
              ((double) atomic_load(&stdev[i].data[2]) / n) - (meanValue * meanValue);
            double stddev = sqrt(variance);
            fprintf(stderr, "Stdev #%d: Total %lld Stdev %.2f\n", i, n, stddev);
        }
    }

    for (int i = 0; i < MAX_DEBUG_SLOTS; ++i)
    {
        if ((n = atomic_load(&correl[i].data[0])))
        {
            double meanX       = (double) atomic_load(&correl[i].data[1]) / n;
            double meanXX      = (double) atomic_load(&correl[i].data[2]) / n;
            double meanY       = (double) atomic_load(&correl[i].data[3]) / n;
            double meanYY      = (double) atomic_load(&correl[i].data[4]) / n;
            double meanXY      = (double) atomic_load(&correl[i].data[5]) / n;
            double numerator   = meanXY - meanX * meanY;
            double denominator = sqrt((meanXX - meanX * meanX) * (meanYY - meanY * meanY));
            double correlation = numerator / denominator;
            fprintf(stderr, "Correl. #%d: Total %lld Coefficient %.2f\n", i, n, correlation);
        }
    }
}
#endif
