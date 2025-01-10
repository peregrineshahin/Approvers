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

#ifndef SEARCH_H
#define SEARCH_H

#include "misc.h"
#include "position.h"
#include "types.h"

// LimitsType struct stores information sent by the caller about the analysis required.
struct LimitsType {
    int       time[2];
    int       inc[2];
    int       depth;
    uint64_t  nodes;
    TimePoint startTime;
};

typedef struct LimitsType LimitsType;

extern LimitsType Limits;

#ifndef KAGGLE
static int use_time_management(void) { return Limits.time[WHITE] || Limits.time[BLACK]; }
#endif

void  search_init(void);
void  search_clear(void);
void  start_thinking(Position* root);
void  prepare_for_search(Position* root);
Value qsearch(Position* pos, Stack* ss, Value alpha, Value beta, Depth depth);
Value search(Position* pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode, int NT);

#endif
