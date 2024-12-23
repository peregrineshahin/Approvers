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

#include "types.h"

Value PieceValue[16] = {0, PawnValue, KnightValue, BishopValue, RookValue, QueenValue};

uint32_t NonPawnPieceValue[16];

// init() initializes piece-square tables: the white halves of the tables
// are copied from Bonus[] adding the piece value, then the black  halves
// of the tables are initialized by flipping and changing the sign of the
// white scores.

SMALL void psqt_init(void) {
    for (int pt = PAWN; pt <= KING; pt++)
    {
        PieceValue[make_piece(BLACK, pt)] = PieceValue[pt];
    }
    union {
        uint16_t val[2];
        uint32_t combi;
    } tmp;
    NonPawnPieceValue[W_PAWN] = NonPawnPieceValue[B_PAWN] = 0;
    for (int pt = KNIGHT; pt < KING; pt++)
    {
        tmp.val[0]                = PieceValue[pt];
        tmp.val[1]                = 0;
        NonPawnPieceValue[pt]     = tmp.combi;
        tmp.val[0]                = 0;
        tmp.val[1]                = PieceValue[pt];
        NonPawnPieceValue[pt + 8] = tmp.combi;
    }
}
