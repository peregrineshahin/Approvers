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

double my_exp(double x) {
    double result = 1.0;
    double term = 1.0;
    for (int n = 1; n < 20; ++n) {
        term *= x / n;
        result += term;
    }
    return result;
}

double my_ln(double x) {
    double result = 0.0;
    double term = (x - 1) / (x + 1);
    double term_squared = term * term;
    double numerator = term;
    for (int n = 1; n < 20; ++n) {
        result += numerator / (2 * n - 1);
        numerator *= term_squared;
    }
    return 2 * result;
}

double my_log10(double x) {
    return my_ln(x) / my_ln(10.0);
}

double my_pow(double base, double exp) {
    return my_exp(exp * my_ln(base));
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

    // If less than one second, gradually reduce mtg
    if (Limits.time[us] < 1000 && (double)mtg / Limits.inc[us] >= 1000)
        mtg = Limits.time[us] * 0.05;

    // Make sure timeLeft is > 0 since we may use it as a divisor
    TimePoint timeLeft =
      max(1, Limits.time[us] + Limits.inc[us] * (mtg - 1) - moveOverhead * (2 + mtg));

    // Calculate time constants based on current time left.
    double logTimeInSec = my_log10(Limits.time[us] / 1000.0);
    double optConstant  = min(0.00308 + 0.000319 * logTimeInSec, 0.00506);
    double maxConstant  = max(3.39 + 3.01 * logTimeInSec, 2.93);

    opt_scale = min(0.0122 + my_pow(ply + 2.95, 0.462) * optConstant, 0.213 * Limits.time[us] / timeLeft);
    max_scale = min(6.64, maxConstant + ply / 12.0);

    // Never use more than 80% of the available time for this move
    Time.optimumTime = opt_scale * timeLeft;
    Time.maximumTime = min(0.825 * Limits.time[us] - moveOverhead, max_scale * Time.optimumTime) - 10;

    if (option_value(OPT_PONDER))
        Time.optimumTime += Time.optimumTime / 4;
}
