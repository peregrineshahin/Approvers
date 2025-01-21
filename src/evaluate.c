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

#include "evaluate.h"
#include "nnue.h"
#include "position.h"

extern int eval_scale;
extern int mat_scale;
extern int mat_n;
extern int mat_b;
extern int mat_r;
extern int mat_q;

Value evaluate(Position* pos, bool isRootStm) {
    Value v = nnue_evaluate(pos);

    v = eval_scale * v / 100;

    if (isRootStm)
        v += 40;
    else
        v -= 40;

    int non_pawn_material = mat_n * popcount(pieces_p(KNIGHT)) + mat_b * popcount(pieces_p(BISHOP))
                          + mat_r * popcount(pieces_p(ROOK)) + mat_q * popcount(pieces_p(QUEEN));

    v = v * (mat_scale + non_pawn_material) / 32768;

    // Damp down the evaluation linearly when shuffling
    v = v * (100 - rule50_count()) / 100;

    return clamp(v, VALUE_MATED_IN_MAX_PLY + 1, VALUE_MATE_IN_MAX_PLY - 1);
}
