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

#ifndef POSITION_H
#define POSITION_H

#include <stddef.h>  // For offsetof()
#include <string.h>

#include "bitboard.h"
#include "nnue.h"
#include "types.h"


extern const char PieceToChar[];
extern Key        matKey[16];

struct Zob {
    Key psq[16][64];
    Key enpassant[8];
    Key castling[16];
    Key side, noPawns;
};

extern struct Zob zob;

void zob_init(void);

// Stack struct stores information needed to restore a Position struct to
// its previous state when we retract a move.

struct Stack {
    // Copied when making a move
    Key pawnKey;
    Key materialKey;
    Key nonPawnKey[2];
    union {
        struct {
            uint8_t pliesFromNull;
            uint8_t rule50;
        };
        uint16_t plyCounters;
    };
    uint8_t castlingRights;

    // Not copied when making a move
    uint8_t    capturedPiece;
    uint8_t    epSquare;
    Key        key;
    Bitboard   checkersBB;
    PVariation pv;

    // Original search stack data
    PieceToHistory* continuationHistory;
    Move            currentMove;
    Move            excludedMove;
    Move            killers[2];
    Value           staticEval;
    Value           statScore;
    int             cutoffCnt;
    int             moveCount;
    bool            ttPv;
    uint8_t         ply;

    // MovePicker data
    uint8_t  stage;
    uint8_t  recaptureSquare;
    Depth    depth;
    Move     ttMove;
    Value    threshold;
    Move     mpKillers[2];
    ExtMove *cur, *endMoves, *endBadCaptures;

    // CheckInfo data
    Bitboard blockersForKing[2];
    union {
        struct {
            Bitboard pinnersForKing[2];
        };
        struct {
            Bitboard dummy;            // pinnersForKing[WHITE]
            Bitboard checkSquares[7];  // element 0 is pinnersForKing[BLACK]
        };
    };
    Square ksq;

    Accumulator accumulator;
};

typedef struct Stack Stack;

#define StateCopySize offsetof(Stack, capturedPiece)
#define StateSize offsetof(Stack, continuationHistory)
#define SStackBegin(st) (&st.continuationHistory)
#define SStackSize (offsetof(Stack, stage) - offsetof(Stack, continuationHistory))


// Position struct stores information regarding the board representation as
// pieces, side to move, hash keys, castling info, etc. The search uses
// the functions do_move() and undo_move() on a Position struct to traverse
// the search tree.

struct Position {
    Stack* st;
    // Board / game representation.
    Bitboard byTypeBB[7];  // no reason to allocate 8 here
    Bitboard byColorBB[2];
    Color    sideToMove;
    uint8_t  board[64];
    uint8_t  castlingRightsMask[64];
    uint16_t gamePly;

    ExtMove* moveList;

    // Relevant mainly to the search of the root position.
    Stack*   stack;
    uint64_t nodes;
    Depth    rootDepth;
    Depth    completedDepth;

    // Pointers to thread-specific tables.
    ButterflyHistory*        mainHistory;
    CapturePieceToHistory*   captureHistory;
    CorrectionHistory*       matCorrHist;
    CorrectionHistory*       pawnCorrHist;
    CorrectionHistory*       prevMoveCorrHist;
    CorrectionHistory*       wNonPawnCorrHist;
    CorrectionHistory*       bNonPawnCorrHist;
    ContinuationHistoryStat* contHist;

    // Thread-control data.
    uint64_t bestMoveChanges;

    void* stackAllocation;
};

// FEN string input/output
void pos_set(Position* pos, char* fen);

//PURE Bitboard attackers_to_occ(const Position *pos, Square s, Bitboard occupied);
PURE Bitboard slider_blockers(const Position* pos, Bitboard sliders, Square s, Bitboard* pinners);

PURE bool is_legal(const Position* pos, Move m);
PURE bool is_pseudo_legal(const Position* pos, Move m);
PURE bool gives_check_special(const Position* pos, Stack* st, Move m);

// Doing and undoing moves
void        do_move(Position* pos, Move m, int givesCheck);
void        undo_move(Position* pos, Move m);
void        do_null_move(Position* pos);
static void undo_null_move(Position* pos);

// Static exchange evaluation
PURE bool see_test(const Position* pos, Move m, int value);

PURE Key  key_after(const Position* pos, Move m);
PURE bool is_draw(const Position* pos);

// Position representation
#define pieces() (pos->byTypeBB[0])
#define pieces_p(p) (pos->byTypeBB[p])
#define pieces_pp(p1, p2) (pos->byTypeBB[p1] | pos->byTypeBB[p2])
#define pieces_c(c) (pos->byColorBB[c])
#define pieces_cp(c, p) (pieces_p(p) & pieces_c(c))
#define pieces_cpp(c, p1, p2) (pieces_pp(p1, p2) & pieces_c(c))
#define piece_on(s) (pos->board[s])
#define ep_square() (pos->st->epSquare)
#define is_empty(s) (!piece_on(s))
#define square_of(c, p) (lsb(pieces_cp(c, p)))

// Castling
#define can_castle_cr(cr) (pos->st->castlingRights & (cr))
#define can_castle_c(c) can_castle_cr((WHITE_OO | WHITE_OOO) << (2 * (c)))

// Checking
#define checkers() (pos->st->checkersBB)

// Attacks to/from a given square
#define attackers_to(s) attackers_to_occ(pos, s, pieces())
#define attacks_from_pawn(s, c) (PawnAttacks[c][s])
#define attacks_from_knight(s) (PseudoAttacks[KNIGHT][s])
#define attacks_from_bishop(s) attacks_bb_bishop(s, pieces())
#define attacks_from_rook(s) attacks_bb_rook(s, pieces())
#define attacks_from_queen(s) (attacks_from_bishop(s) | attacks_from_rook(s))
#define attacks_from_king(s) (PseudoAttacks[KING][s])
#define attacks_from(pc, s) attacks_bb(pc, s, pieces())

// Properties of moves
#define moved_piece(m) (piece_on(from_sq(m)))
#define captured_piece() (pos->st->capturedPiece)

// Accessing hash keys
#define key() (pos->st->key)
#define material_key() (pos->st->materialKey)
#define pawn_key() (pos->st->pawnKey)
#define w_nonpawn_key() (pos->st->nonPawnKey[WHITE])
#define b_nonpawn_key() (pos->st->nonPawnKey[BLACK])

// Other properties of the position
#define stm() (pos->sideToMove)
#define game_ply() (pos->gamePly)
#define nodes_searched() (pos->nodes)
#define rule50_count() (pos->st->rule50)

static Bitboard blockers_for_king(const Position* pos, Color c) {
    return pos->st->blockersForKing[c];
}

static bool is_capture_or_promotion(const Position* pos, Move m) {
    return type_of_m(m) != NORMAL ? type_of_m(m) != CASTLING : !is_empty(to_sq(m));
}

static bool is_capture(const Position* pos, Move m) {
    // Castling is encoded as "king captures the rook"
    return (!is_empty(to_sq(m)) && type_of_m(m) != CASTLING) || type_of_m(m) == ENPASSANT;
}

static bool gives_check(const Position* pos, Stack* st, Move m) {
    return type_of_m(m) == NORMAL && !(blockers_for_king(pos, !stm()) & pieces_c(stm()))
           ? (bool) (st->checkSquares[type_of_p(moved_piece(m))] & sq_bb(to_sq(m)))
           : gives_check_special(pos, st, m);
}

static bool non_pawn_material(const Position* pos) {
    return pieces_cpp(stm(), PAWN, KING) != pieces_c(stm());
}

static bool low_material(const Position* pos) {
    const Bitboard rooks = pieces_p(ROOK);
    return pieces_p(QUEEN) == 0 && !more_than_one(rooks)
        && popcount(pieces_pp(KNIGHT, BISHOP)) + 2 * (rooks != 0) <= 3;
}

void pos_set_check_info(Position* pos);

// undo_null_move is used to undo a null move.

static void undo_null_move(Position* pos) {

    pos->st--;
    pos->sideToMove = !pos->sideToMove;
}

// attackers_to() computes a bitboard of all pieces which attack a given
// square. Slider attacks use the occupied bitboard to indicate occupancy.

static Bitboard attackers_to_occ(const Position* pos, Square s, Bitboard occupied) {
    return (attacks_from_pawn(s, BLACK) & pieces_cp(WHITE, PAWN))
         | (attacks_from_pawn(s, WHITE) & pieces_cp(BLACK, PAWN))
         | (attacks_from_knight(s) & pieces_p(KNIGHT))
         | (attacks_bb_rook(s, occupied) & pieces_pp(ROOK, QUEEN))
         | (attacks_bb_bishop(s, occupied) & pieces_pp(BISHOP, QUEEN))
         | (attacks_from_king(s) & pieces_p(KING));
}

#endif
