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

#ifndef BITBOARD_H
#define BITBOARD_H

#include "types.h"
void bitboards_init(void);

#define FileABB 0x0101010101010101ULL
#define FileBBB (FileABB << 1)
#define FileCBB (FileABB << 2)
#define FileDBB (FileABB << 3)
#define FileEBB (FileABB << 4)
#define FileFBB (FileABB << 5)
#define FileGBB (FileABB << 6)
#define FileHBB (FileABB << 7)

#define Rank1BB 0xFFULL
#define Rank2BB (Rank1BB << (8 * 1))
#define Rank3BB (Rank1BB << (8 * 2))
#define Rank4BB (Rank1BB << (8 * 3))
#define Rank5BB (Rank1BB << (8 * 4))
#define Rank6BB (Rank1BB << (8 * 5))
#define Rank7BB (Rank1BB << (8 * 6))
#define Rank8BB (Rank1BB << (8 * 7))

extern Bitboard BetweenBB[64][64];
extern Bitboard LineBB[64][64];
extern Bitboard PseudoAttacks[8][64];
extern Bitboard PawnAttacks[2][64];


static __attribute__((pure)) Bitboard sq_bb(Square s) { return 1LL << s; }

static uint64_t more_than_one(Bitboard b) { return b & (b - 1); }

static Bitboard rank_bb(Rank r) { return Rank1BB << (r << 3); }

static Bitboard file_bb(File f) { return FileABB << f; }


// Moves a bitboard one or two steps as specified by the direction.
static Bitboard shift_bb(int direction, Bitboard b) {
    return direction == NORTH         ? b << 8
         : direction == SOUTH         ? b >> 8
         : direction == NORTH + NORTH ? b << 16
         : direction == SOUTH + SOUTH ? b >> 16
         : direction == EAST          ? (b & ~FileHBB) << 1
         : direction == WEST          ? (b & ~FileABB) >> 1
         : direction == NORTH_EAST    ? (b & ~FileHBB) << 9
         : direction == SOUTH_EAST    ? (b & ~FileHBB) >> 7
         : direction == NORTH_WEST    ? (b & ~FileABB) << 7
         : direction == SOUTH_WEST    ? (b & ~FileABB) >> 9
                                      : 0;
}


// Returns a bitboard representing the squares in the semi-open
// segment between the squares s1 and s2 (excluding s1 but including s2).
static Bitboard between_bb(Square s1, Square s2) { return BetweenBB[s1][s2]; }


// Returns true if the squares s1, s2 and s3 are aligned either
// on a straight or on a diagonal line.
static uint64_t aligned(Move m, Square s) { return ((Bitboard*) LineBB)[m & 4095] & sq_bb(s); }

static unsigned distance_f(Square x, Square y) {
    unsigned f1 = file_of(x), f2 = file_of(y);
    return f1 < f2 ? f2 - f1 : f1 - f2;
}

static unsigned distance_r(Square x, Square y) {
    unsigned r1 = rank_of(x), r2 = rank_of(y);
    return r1 < r2 ? r2 - r1 : r1 - r2;
}

// Returns the distance between x and y, defined as the number of steps
// for a king in x to reach y. Works with squares, ranks, files.
static int distance(Square x, Square y) { return max(distance_f(x, y), distance_r(x, y)); }

#define attacks_bb_queen(s, occupied) \
    (attacks_bb_bishop((s), (occupied)) | attacks_bb_rook((s), (occupied)))

#if defined(MAGIC_PLAIN)
    #include "magic-plain.h"
#elif defined(AVX2_BITBOARD)
    #include "avx2-bitboard.h"
#endif

static Bitboard attacks_bb(int pt, Square s, Bitboard occupied) {

    switch (pt)
    {
    case BISHOP :
        return attacks_bb_bishop(s, occupied);
    case ROOK :
        return attacks_bb_rook(s, occupied);
    case QUEEN :
        return attacks_bb_queen(s, occupied);
    default :
        return PseudoAttacks[pt][s];
    }
}


// Counts the number of non-zero bits in a bitboard.
static int popcount(Bitboard b) { return __builtin_popcountll(b); }


// Returns the least significant bit in a non-zero bitboard.
static int lsb(Bitboard b) { return __builtin_ctzll(b); }

// Finds and clears the least significant bit in a non-zero bitboard.
static Square pop_lsb(Bitboard* b) {
    const Square s = lsb(*b);
    *b &= *b - 1;
    return s;
}

#endif
