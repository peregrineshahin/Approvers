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

extern int eval_scale;

Value evaluate(Position* pos) {

    int scale = 1700 + 16 * (piece_count(WHITE, PAWN) + piece_count(BLACK, PAWN));

    Value nnue = nnue_evaluate(pos) * scale / 2048;

    // Damp down the evaluation linearly when shuffling
    Value v = nnue * (100 - rule50_count()) / 100;

    v = eval_scale * v / 100;
    // v = (v / 16) * 16;
    // v = (stm() == WHITE ? v : -v) + Tempo

    return clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);
}
