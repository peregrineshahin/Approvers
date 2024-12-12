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

#include <fcntl.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
    #include <windows.h>
#else
    #include <sys/mman.h>
#endif

#include "misc.h"
#include "thread.h"

// Version number. If Version is left empty, then compile date in the format
// DD-MM-YY and show in engine_info.
char Version[] = "";

#ifndef _WIN32
pthread_mutex_t ioMutex = PTHREAD_MUTEX_INITIALIZER;
#else
HANDLE ioMutex;
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

ssize_t getline(char** lineptr, size_t* n, FILE* stream) {
    if (*n == 0)
        *lineptr = malloc(*n = 100);

    int    c = 0;
    size_t i = 0;
    while ((c = getc(stream)) != EOF)
    {
        (*lineptr)[i++] = c;
        if (i == *n)
            *lineptr = realloc(*lineptr, *n += 100);
        if (c == '\n')
            break;
    }
    (*lineptr)[i] = 0;
    return i;
}

#ifdef _WIN32
typedef SIZE_T(WINAPI* GLPM)(void);
// The following two functions were taken from mingw_lock.c

void __cdecl _lock(int locknum);
void __cdecl _unlock(int locknum);
    #define _STREAM_LOCKS 16
    #define _IOLOCKED 0x8000
typedef struct {
    FILE             f;
    CRITICAL_SECTION lock;
} _FILEX;
#endif
