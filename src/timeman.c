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
#include "uci.h"

struct TimeManagement Time;  // Our global time management struct

double my_sqrt(double x) {
    if (x <= 0)
        return 0;
    double guess = x / 2.0;
    for (int i = 0; i < 10; ++i)
        guess = (guess + x / guess) / 2.0;
    return guess;
}

// tm_init() is called at the beginning of the search and calculates the
// time bounds allowed for the current game ply. We currently support:
// 1) x basetime (+z increment)
// 2) x moves in y seconds (+z increment)

void time_init(Color us, int ply) {
    int moveOverhead = option_value(OPT_MOVE_OVERHEAD);

    // opt_scale is a percentage of available time to use for the current move.
    // max_scale is a multiplier applied to optimumTime.
    double opt_scale, max_scale;

    Time.startTime = Limits.startTime;

    // Maximum move horizon of 50 moves
    int mtg = 50;

    // Make sure that timeLeft > 0 since we may use it as a divisor
    TimePoint timeLeft =
      max(1, Limits.time[us] + Limits.inc[us] * (mtg - 1) - moveOverhead * (2 + mtg));

    timeLeft = 100 * timeLeft / 100;

    // x basetime (+z increment)
    // If there is a healthy increment, timeLeft can exceed actual available
    // game time for the current move, so also cap to 20% of available game time.
    opt_scale = min(0.008 + my_sqrt(ply + 3.0) / 250.0, 0.2 * Limits.time[us] / (double) timeLeft);
    max_scale = min(7.0, 4.0 + ply / 12.0);

    // Never use more than 80% of the available time for this move
    Time.optimumTime = opt_scale * timeLeft;
    Time.maximumTime = min(0.8 * Limits.time[us] - moveOverhead, max_scale * Time.optimumTime);

    if (option_value(OPT_PONDER))
        Time.optimumTime += Time.optimumTime / 4;
}
