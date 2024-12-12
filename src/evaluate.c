/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2015 Marco Costalba, Joona Kiiski, Tord Romstad
  Copyright (C) 2015-2018 Marco Costalba, Joona Kiiski, Gary Linscott, Tord Romstad

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

#include "bitboard.h"
#include "evaluate.h"
#include "nnue.h"
#include "position.h"

#define MAX_PHASE 62

Value evaluate(Position* pos) {
    Value v = nnue_evaluate(pos);

    const int phase = 9 * popcount(pieces_p(QUEEN)) + 5 * popcount(pieces_p(ROOK))
                    + 3 * popcount(pieces_p(BISHOP) | popcount(pieces_p(KNIGHT)));


    v -= v * (MAX_PHASE - phase) / (5 * MAX_PHASE);

    // Damp down the evaluation linearly when shuffling
    v = v * (100 - rule50_count()) / 100;

    // v = (v / 16) * 16;
    // v = (stm() == WHITE ? v : -v) + Tempo

    return clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);
}
