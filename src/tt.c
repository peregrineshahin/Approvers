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

#ifndef DKAGGLE
    #include <inttypes.h>
#endif

#include <stdio.h>
#include <string.h>  // For memset
#ifndef _WIN32
    #include <sys/mman.h>
#endif

#include "bitboard.h"
#include "settings.h"
#include "thread.h"
#include "tt.h"
#include "types.h"
#include "uci.h"

TranspositionTable TT;  // Our global transposition table

// tt_free() frees the allocated transposition table memory.

SMALL void tt_free(void) {
#ifdef _WIN32
    if (TT.mem)
        VirtualFree(TT.mem, 0, MEM_RELEASE);
#else
    if (TT.mem)
        munmap(TT.mem, TT.allocSize);
#endif
    TT.mem = NULL;
}


// tt_allocate() allocates the transposition table, measured in megabytes.

SMALL void tt_allocate(size_t mbSize) {
    TT.clusterCount = mbSize * 1024 * 1024 / sizeof(Cluster);
    size_t size     = TT.clusterCount * sizeof(Cluster);

#ifdef _WIN32

    TT.mem = NULL;

    if (!TT.mem)
        TT.mem = VirtualAlloc(NULL, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    #ifndef DKAGGLE
    if (!TT.mem)
        goto failed;
    #endif

    TT.table = (Cluster*) TT.mem;

#else /* Unix */

    size_t alignment = 1;
    size_t allocSize = size + alignment - 1;

    #if defined(__APPLE__) && defined(VM_FLAGS_SUPERPAGE_SIZE_2MB)

    if (!TT.mem)
        TT.mem = mmap(NULL, allocSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    #else

    TT.mem = mmap(NULL, allocSize, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    #endif

    TT.allocSize = allocSize;
    TT.table     = (Cluster*) ((((uintptr_t) TT.mem) + alignment - 1) & ~(alignment - 1));
    #ifndef DKAGGLE
    if (!TT.mem)
        goto failed;
    #endif

    #if defined(__linux__) && defined(MADV_HUGEPAGE)

    #endif
#endif

    // Clear the TT table to page in the memory immediately. This avoids
    // an initial slow down during the first second or minutes of the search.
    tt_clear();
    return;

#ifndef DKAGGLE
failed:
    fprintf(stderr,
            "Failed to allocate %" PRIu64 "MB for "
            "transposition table.\n",
            (uint64_t) mbSize);
    exit(EXIT_FAILURE);
#endif
}


// tt_clear() initialises the entire transposition table to zero.

SMALL void tt_clear(void) {
    if (TT.table)
    {
        memset((uint8_t*) TT.table, 0, TT.clusterCount * sizeof(Cluster));
    }
}


// tt_probe() looks up the current position in the transposition table.
// It returns true and a pointer to the TTEntry if the position is found.
// Otherwise, it returns false and a pointer to an empty or least valuable
// TTEntry to be replaced later. The replace value of an entry is
// calculated as its depth minus 8 times its relative age. TTEntry t1 is
// considered more valuable than TTEntry t2 if its replace value is greater
// than that of t2.

SMALL TTEntry* tt_probe(Key key, bool* found) {
    TTEntry* tte   = tt_first_entry(key);
    uint16_t key16 = key;  // Use the low 16 bits as key inside the cluster

    for (int i = 0; i < ClusterSize; i++)
        if (tte[i].key16 == key16 || !tte[i].depth8)
        {
            //      if ((tte[i].genBound8 & 0xF8) != TT.generation8 && tte[i].key16)
            tte[i].genBound8 = TT.generation8 | (tte[i].genBound8 & 0x7);  // Refresh
            *found           = tte[i].depth8;
            return &tte[i];
        }

    // Find an entry to be replaced according to the replacement strategy
    TTEntry* replace = tte;
    for (int i = 1; i < ClusterSize; i++)
        // Due to our packed storage format for generation and its cyclic
        // nature we add 263 (256 is the modulus plus 7 to keep the unrelated
        // lowest three bits from affecting the result) to calculate the entry
        // age correctly even after generation8 overflows into the next cycle.
        if (replace->depth8 - ((263 + TT.generation8 - replace->genBound8) & 0xF8)
            > tte[i].depth8 - ((263 + TT.generation8 - tte[i].genBound8) & 0xF8))
            replace = &tte[i];

    *found = false;
    return replace;
}