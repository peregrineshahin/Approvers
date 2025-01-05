/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2015 Marco Costalba, Joona Kiiski, Tord Romstad
  Copyright (C) 2015-2017 Marco Costalba, Joona Kiiski, Gary Linscott, Tord Romstad

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

#include "movepick.h"
#include "search.h"
#include "thread.h"

ThreadStruct Thread = {0};

void thread_init() {
    Thread.testPonder = 0;

    Position* pos        = calloc(sizeof(Position), 1);
    pos->mainHistory     = calloc(sizeof(ButterflyHistory), 1);
    pos->captureHistory  = calloc(sizeof(CapturePieceToHistory), 1);
    pos->corrHists       = calloc(sizeof(CorrectionHistory), 1);
    pos->stackAllocation = calloc(63 + (MAX_PLY + 110) * sizeof(Stack), 1);
    pos->moveList        = calloc(10000 * sizeof(ExtMove), 1);

    pos->stack    = (Stack*) (((uintptr_t) pos->stackAllocation + 0x3f) & ~0x3f);
    pos->contHist = calloc(sizeof(ContinuationHistoryStat), 1);
#pragma clang loop unroll(disable)
    for (int pc = 0; pc < 15; pc++)
#pragma clang loop unroll(disable)
        for (int sq = 0; sq < 64; sq++)
            (*pos->contHist)[0][0][pc][sq] = -1;

    Thread.pos = pos;

    search_init();
}

void thread_exit() {
    Position* pos = Thread.pos;

    free(pos->mainHistory);
    free(pos->captureHistory);
    free(pos->stackAllocation);
    free(pos->moveList);
    free(pos->corrHists);
    free(pos->contHist);
    free(pos);
}
