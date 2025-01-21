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

#include <float.h>
#include "search.h"
#include "timeman.h"

int MoveOverhead = 10;

extern int tm_v13;
extern int tm_v14;
extern int tm_v15;
extern int tm_v16;
extern int tm_v17;
extern int tm_v18;
extern int tm_v19;
extern int tm_v20;

struct TimeManagement Time;  // Our global time management struct

SMALL double my_sqrt(double x) {
    if (x <= 0)
        return 0;
    double guess = x / 2.0;
#pragma clang loop unroll(disable)
    for (int i = 0; i < 10; ++i)
        guess = (guess + x / guess) / 2.0;
    return guess;
}

// Called at the beginning of the search and calculates
// the time bounds allowed for the current game ply.
void time_init(Color us, int ply) {
    // opt_scale is a percentage of available time to use for the current move.
    // max_scale is a multiplier applied to optimumTime.
    double opt_scale, max_scale;

    Time.startTime = Limits.startTime;

    int mtg = 50;
    if (Limits.time[us] < 1000 && (double) mtg / Limits.time[us] > 0.05)
        mtg = (int) (Limits.time[us] * 0.05);

    // Make sure that timeLeft > 0 since we may use it as a divisor
    TimePoint timeLeft = max(1, Limits.time[us] - MoveOverhead * (2 + mtg));

    // If there is a healthy increment, timeLeft can exceed actual available
    // game time for the current move, so also cap to 20% of available game time.
    opt_scale = min(tm_v13 / 10000.0 + my_sqrt(ply + tm_v14 / 100.0) / (double) tm_v15,
                    tm_v16 / 100.0 * Limits.time[us] / (double) timeLeft);
    max_scale = min(tm_v17 / 100.0, tm_v18 / 100.0 + ply / (tm_v19 / 100.0));

    // Never use more than 80% of the available time for this move
    Time.optimumTime = (Limits.inc[us] / 2) + opt_scale * timeLeft;
    Time.maximumTime = (Limits.inc[us] / 2) +
      min(tm_v20 / 1000.0 * Limits.time[us] - MoveOverhead, max_scale * Time.optimumTime);
}
