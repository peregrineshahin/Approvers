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

#include <stdio.h>
#include <string.h>  // For memset
#ifndef _WIN32
    #include <sys/mman.h>
#endif

#include "tt.h"
#include "types.h"

TranspositionTable TT;  // Our global transposition table

// Frees the allocated transposition table memory.
void tt_free(void) {
#ifdef _WIN32
    if (TT.mem)
        VirtualFree(TT.mem, 0, MEM_RELEASE);
#else
    if (TT.mem)
        munmap(TT.mem, TT.allocSize);
#endif
    TT.mem = NULL;
}


// Allocates the transposition table, measured in megabytes.
void tt_allocate(size_t mbSize) {
    TT.count    = mbSize * 1024 * 1024 / sizeof(TTEntry);
    size_t size = TT.count * sizeof(TTEntry);

#ifdef _WIN32
    TT.mem   = VirtualAlloc(NULL, size, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    TT.table = (TTEntry*) TT.mem;
#else /* Unix */
    TT.mem = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    TT.allocSize = size;
    TT.table     = (TTEntry*) ((uintptr_t) TT.mem & ~0);
#endif

    // Clear the TT table to page in the memory immediately. This avoids
    // an initial slow down during the first second or minutes of the search.
    tt_clear();
}


// Initialises the entire transposition table to zero.
void tt_clear(void) {
    if (TT.table)
        memset((uint8_t*) TT.table, 0, TT.count * sizeof(TTEntry));
}


// Looks up the current position in the transposition table.
// It returns true and a pointer to the TTEntry if the position is found.
// Otherwise, it returns false and a pointer to an empty or least valuable
// TTEntry to be replaced later. The replace value of an entry is
// calculated as its depth minus 8 times its relative age. TTEntry t1 is
// considered more valuable than TTEntry t2 if its replace value is greater
// than that of t2.
TTEntry* tt_probe(Key key, bool* found) {
    TTEntry* tte = tt_entry(key);
    *found       = tte->key16 == (uint16_t) key;
    return tte;
}

void tte_save(TTEntry* tte, Key k, Value v, bool pv, int b, Depth d, Move m, Value ev) {
    // Don't overwrite more valuable entries
    if (!(b == BOUND_EXACT || tte->key16 != (uint16_t) k
          || d - DEPTH_OFFSET + 4 + 2 * pv > tte->depth8))
        return;

    // Preserve any existing move for the same position
    if (m || (uint16_t) k != tte->key16)
        tte->move16 = (uint16_t) m;

    tte->key16   = (uint16_t) k;
    tte->depth8  = (uint8_t) (d - DEPTH_OFFSET);
    tte->bound8  = (uint8_t) (((uint8_t) pv << 2) | b);
    tte->value16 = (int16_t) v;
    tte->eval16  = (int16_t) ev;
}
