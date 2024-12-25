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

extern uint8_t SquareDistance[64][64];

extern Bitboard FileBB[8];
extern Bitboard RankBB[8];
extern Bitboard BetweenBB[64][64];
extern Bitboard LineBB[64][64];
extern Bitboard PseudoAttacks[8][64];
extern Bitboard PawnAttacks[2][64];


static __attribute__((pure)) Bitboard sq_bb(Square s) { return 1LL << s; }

static uint64_t more_than_one(Bitboard b) { return b & (b - 1); }


// rank_bb() and file_bb() return a bitboard representing all the squares on
// the given file or rank.

static Bitboard rank_bb(Rank r) { return RankBB[r]; }

static Bitboard rank_bb_s(Square s) { return RankBB[rank_of(s)]; }

static Bitboard file_bb(File f) { return FileBB[f]; }

static Bitboard file_bb_s(Square s) { return FileBB[file_of(s)]; }


// shift_bb() moves a bitboard one step along direction Direction.
static Bitboard shift_bb(int Direction, Bitboard b) {
    return Direction == NORTH         ? b << 8
         : Direction == SOUTH         ? b >> 8
         : Direction == NORTH + NORTH ? b << 16
         : Direction == SOUTH + SOUTH ? b >> 16
         : Direction == EAST          ? (b & ~FileHBB) << 1
         : Direction == WEST          ? (b & ~FileABB) >> 1
         : Direction == NORTH_EAST    ? (b & ~FileHBB) << 9
         : Direction == SOUTH_EAST    ? (b & ~FileHBB) >> 7
         : Direction == NORTH_WEST    ? (b & ~FileABB) << 7
         : Direction == SOUTH_WEST    ? (b & ~FileABB) >> 9
                                      : 0;
}


// between_bb() returns a bitboard representing all the squares between
// the two given ones. For instance, between_bb(SQ_C4, SQ_F7) returns a
// bitboard with the bits for square d5 and e6 set. If s1 and s2 are not
// on the same rank, file or diagonal, 0 is returned.

static Bitboard between_bb(Square s1, Square s2) { return BetweenBB[s1][s2]; }


// aligned() returns true if square s is on the line determined by move m.

static uint64_t aligned(Move m, Square s) { return ((Bitboard*) LineBB)[m & 4095] & sq_bb(s); }


// distance() functions return the distance between x and y, defined as
// the number of steps for a king in x to reach y. Works with squares,
// ranks, files.

static int distance(Square x, Square y) { return SquareDistance[x][y]; }

static unsigned distance_f(Square x, Square y) {
    unsigned f1 = file_of(x), f2 = file_of(y);
    return f1 < f2 ? f2 - f1 : f1 - f2;
}

static unsigned distance_r(Square x, Square y) {
    unsigned r1 = rank_of(x), r2 = rank_of(y);
    return r1 < r2 ? r2 - r1 : r1 - r2;
}

#define attacks_bb_queen(s, occupied) \
    (attacks_bb_bishop((s), (occupied)) | attacks_bb_rook((s), (occupied)))

#if defined(MAGIC_FANCY)
    #include "magic-fancy.h"
#elif defined(MAGIC_PLAIN)
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


// popcount() counts the number of non-zero bits in a bitboard.

static int popcount(Bitboard b) { return __builtin_popcountll(b); }


// lsb() return the least significant bit in a non-zero
// bitboard.

static int lsb(Bitboard b) { return __builtin_ctzll(b); }

// pop_lsb() finds and clears the least significant bit in a non-zero
// bitboard.

static Square pop_lsb(Bitboard* b) {
    const Square s = lsb(*b);
    *b &= *b - 1;
    return s;
}

#endif
