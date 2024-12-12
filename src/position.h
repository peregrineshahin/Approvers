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


#ifndef _WIN32
    #include <pthread.h>
#endif
#include <stdatomic.h>
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

void psqt_init(void);
void zob_init(void);

// Stack struct stores information needed to restore a Position struct to
// its previous state when we retract a move.

struct Stack {
    // Copied when making a move
    Key pawnKey;
    Key materialKey;
    union {
        struct {
            Score psq;
            union {
                uint16_t nonPawnMaterial[2];
                uint32_t nonPawn;
            };
        };
        uint64_t psqnpm;
    };
    union {
        struct {
            uint8_t pliesFromNull;
            uint8_t rule50;
        };
        uint16_t plyCounters;
    };
    uint8_t castlingRights;

    // Not copied when making a move
    uint8_t  capturedPiece;
    uint8_t  epSquare;
    Key      key;
    Bitboard checkersBB;

    // Original search stack data
    Move*           pv;
    PieceToHistory* history;
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
    uint8_t  mp_ply;
    Move     countermove;
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
#define StateSize offsetof(Stack, pv)
#define SStackBegin(st) (&st.pv)
#define SStackSize (offsetof(Stack, countermove) - offsetof(Stack, pv))


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
    uint8_t  pieceCount[16];
    uint8_t  pieceList[256];
    uint8_t  index[64];
    uint8_t  castlingRightsMask[64];
    uint8_t  castlingRookSquare[16];
    Bitboard castlingPath[16];
    Key      rootKeyFlip;
    uint16_t gamePly;
    bool     hasRepeated;

    ExtMove* moveList;

    // Relevant mainly to the search of the root position.
    RootMoves* rootMoves;
    Stack*     stack;
    uint64_t   nodes;
    uint64_t   tbHits;
    uint64_t   ttHitAverage;
    int        pvIdx, pvLast;
    int        nmpMinPly;
    Color      nmpColor;
    Depth      rootDepth;
    Depth      completedDepth;
    Score      contempt;

    // Pointers to thread-specific tables.
    CounterMoveStat*        counterMoves;
    ButterflyHistory*       history;
    LowPlyHistory*          lowPlyHistory;
    CapturePieceToHistory*  captureHistory;
    correction_history_t    corrHistory;
    PawnEntry*              pawnTable;
    CounterMoveHistoryStat* counterMoveHistory;

    // Thread-control data.
    uint64_t    bestMoveChanges;
    atomic_bool resetCalls;
    int         callsCnt;
    int         action;
    int         threadIdx;
#ifndef _WIN32
    pthread_t       nativeThread;
    pthread_mutex_t mutex;
    pthread_cond_t  sleepCondition;
#else
    HANDLE nativeThread;
    HANDLE startEvent, stopEvent;
#endif
    void* stackAllocation;
};

// FEN string input/output
void pos_set(Position* pos, char* fen);
void pos_fen(const Position* pos, char* fen);

//PURE Bitboard attackers_to_occ(const Position *pos, Square s, Bitboard occupied);
PURE Bitboard slider_blockers(const Position* pos, Bitboard sliders, Square s, Bitboard* pinners);

PURE bool is_legal(const Position* pos, Move m);
PURE bool is_pseudo_legal(const Position* pos, Move m);
PURE bool gives_check_special(const Position* pos, Stack* st, Move m);

// Doing and undoing moves
void        do_move(Position* pos, Move m, int givesCheck);
void        undo_move(Position* pos, Move m);
void        do_null_move(Position* pos);
INLINE void undo_null_move(Position* pos);

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
#define piece_count(c, p) (pos->pieceCount[8 * (c) + (p)] - (8 * (c) + (p)) * 16)
#define piece_list(c, p) (&pos->pieceList[16 * (8 * (c) + (p))])
#define square_of(c, p) (pos->pieceList[16 * (8 * (c) + (p))])
#define loop_through_pieces(c, p, s) \
    const uint8_t* pl = piece_list(c, p); \
    while ((s = *pl++) != SQ_NONE)
#define piece_count_mk(c, p) (((material_key()) >> (20 * (c) + 4 * (p) + 4)) & 15)

// Castling
#define can_castle_cr(cr) (pos->st->castlingRights & (cr))
#define can_castle_c(c) can_castle_cr((WHITE_OO | WHITE_OOO) << (2 * (c)))
#define can_castle_any() (pos->st->castlingRights)
#define castling_impeded(cr) (pieces() & pos->castlingPath[cr])
#define castling_rook_square(cr) (pos->castlingRookSquare[cr])

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

// Other properties of the position
#define stm() (pos->sideToMove)
#define game_ply() (pos->gamePly)
#define nodes_searched() (pos->nodes)
#define rule50_count() (pos->st->rule50)
#define psq_score() (pos->st->psq)
#define non_pawn_material_c(c) (pos->st->nonPawnMaterial[c])
#define non_pawn_material() (non_pawn_material_c(WHITE) + non_pawn_material_c(BLACK))
#define pawns_only() (!pos->st->nonPawn)

INLINE Bitboard blockers_for_king(const Position* pos, Color c) {
    return pos->st->blockersForKing[c];
}

INLINE bool is_discovery_check_on_king(const Position* pos, Color c, Move m) {
    return pos->st->blockersForKing[c] & sq_bb(from_sq(m));
}

INLINE bool pawn_passed(const Position* pos, Color c, Square s) {
    return !(pieces_cp(!c, PAWN) & passed_pawn_span(c, s));
}

INLINE bool advanced_pawn_push(const Position* pos, Move m) {
    return type_of_p(moved_piece(m)) == PAWN && relative_rank_s(stm(), from_sq(m)) > RANK_4;
}

INLINE bool opposite_bishops(const Position* pos) {
    return piece_count(WHITE, BISHOP) == 1 && piece_count(BLACK, BISHOP) == 1
        && (pieces_p(BISHOP) & DarkSquares) && (pieces_p(BISHOP) & ~DarkSquares);
}

INLINE bool is_capture_or_promotion(const Position* pos, Move m) {
    return type_of_m(m) != NORMAL ? type_of_m(m) != CASTLING : !is_empty(to_sq(m));
}

INLINE bool is_capture(const Position* pos, Move m) {
    // Castling is encoded as "king captures the rook"
    return (!is_empty(to_sq(m)) && type_of_m(m) != CASTLING) || type_of_m(m) == ENPASSANT;
}

INLINE bool gives_check(const Position* pos, Stack* st, Move m) {
    return type_of_m(m) == NORMAL && !(blockers_for_king(pos, !stm()) & pieces_c(stm()))
           ? (bool) (st->checkSquares[type_of_p(moved_piece(m))] & sq_bb(to_sq(m)))
           : gives_check_special(pos, st, m);
}

void pos_set_check_info(Position* pos);

// undo_null_move is used to undo a null move.

INLINE void undo_null_move(Position* pos) {

    pos->st--;
    pos->sideToMove = !pos->sideToMove;
}

// attackers_to() computes a bitboard of all pieces which attack a given
// square. Slider attacks use the occupied bitboard to indicate occupancy.

INLINE Bitboard attackers_to_occ(const Position* pos, Square s, Bitboard occupied) {
    return (attacks_from_pawn(s, BLACK) & pieces_cp(WHITE, PAWN))
         | (attacks_from_pawn(s, WHITE) & pieces_cp(BLACK, PAWN))
         | (attacks_from_knight(s) & pieces_p(KNIGHT))
         | (attacks_bb_rook(s, occupied) & pieces_pp(ROOK, QUEEN))
         | (attacks_bb_bishop(s, occupied) & pieces_pp(BISHOP, QUEEN))
         | (attacks_from_king(s) & pieces_p(KING));
}

#endif
