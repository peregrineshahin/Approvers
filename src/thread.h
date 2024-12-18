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

#define MAX_THREADS 1

void thread_search(Position* pos);

// MainThread struct seems to exist mostly for easy move.

struct MainThread {
    double previousTimeReduction;
    Value  previousScore;
    Value  iterValue[4];
};

typedef struct MainThread MainThread;

extern MainThread mainThread;

void mainthread_search(void);


// ThreadPool struct handles all the threads-related stuff like init,
// starting, parking and, most importantly, launching a thread. All the
// access to threads data is done through this class.

struct ThreadPool {
    Position* pos[MAX_THREADS];
    int       numThreads;
    bool      ponder, stop, increaseDepth;
};

typedef struct ThreadPool ThreadPool;

void threads_init(void);
void threads_exit(void);
void threads_set_number(int num);

extern ThreadPool Threads;

static Position* threads_main(void) { return Threads.pos[0]; }

extern CounterMoveHistoryStat** cmhTables;
extern int                      numCmhTables;

#endif
