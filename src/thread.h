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

#ifndef THREAD_H
#define THREAD_H

#include "types.h"

void thread_search(Position* pos);
void mainthread_search(void);

typedef struct ThreadStruct ThreadStruct;

struct ThreadStruct {
    Position* pos;
    double    previousTimeReduction;
    Value     previousScore;
    Value     iterValue[4];
    bool      ponder, stop, increaseDepth;
    // Flag for testing pondering outside the Kaggle environment
    bool testPonder;
};

extern ThreadStruct Thread;

void thread_init();
void thread_exit(void);

#endif
