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


#include "movegen.h"
#include "position.h"
#include "types.h"


INLINE ExtMove* make_promotions(ExtMove* list, Square to, Square ksq, const int Type, const int D) {
    if (Type == CAPTURES || Type == EVASIONS || Type == NON_EVASIONS)
    {
        (list++)->move = make_promotion(to - D, to, QUEEN);
        if (attacks_from_knight(to) & sq_bb(ksq))
            (list++)->move = make_promotion(to - D, to, KNIGHT);
    }

    if (Type == QUIETS || Type == EVASIONS || Type == NON_EVASIONS)
    {
        (list++)->move = make_promotion(to - D, to, ROOK);
        (list++)->move = make_promotion(to - D, to, BISHOP);
        if (!(attacks_from_knight(to) & sq_bb(ksq)))
            (list++)->move = make_promotion(to - D, to, KNIGHT);
    }

    return list;
}


INLINE ExtMove* generate_pawn_moves(
  const Position* pos, ExtMove* list, Bitboard target, const Color Us, const int Type) {
    // Compute our parametrized parameters at compile time, named according to
    // the point of view of white side.
    const Color    Them     = (Us == WHITE ? BLACK : WHITE);
    const Bitboard TRank8BB = (Us == WHITE ? Rank8BB : Rank1BB);
    const Bitboard TRank7BB = (Us == WHITE ? Rank7BB : Rank2BB);
    const Bitboard TRank3BB = (Us == WHITE ? Rank3BB : Rank6BB);
    const int      Up       = (Us == WHITE ? NORTH : SOUTH);
    const int      Right    = (Us == WHITE ? NORTH_EAST : SOUTH_WEST);
    const int      Left     = (Us == WHITE ? NORTH_WEST : SOUTH_EAST);

    Bitboard emptySquares;

    Bitboard pawnsOn7    = pieces_cp(Us, PAWN) & TRank7BB;
    Bitboard pawnsNotOn7 = pieces_cp(Us, PAWN) & ~TRank7BB;

    Bitboard enemies = (Type == EVASIONS   ? pieces_c(Them) & target
                        : Type == CAPTURES ? target
                                           : pieces_c(Them));

    // Single and double pawn pushes, no promotions
    if (Type != CAPTURES)
    {
        emptySquares = (Type == QUIETS || Type == QUIET_CHECKS ? target : ~pieces());

        Bitboard b1 = shift_bb(Up, pawnsNotOn7) & emptySquares;
        Bitboard b2 = shift_bb(Up, b1 & TRank3BB) & emptySquares;

        if (Type == EVASIONS)
        {  // Consider only blocking squares
            b1 &= target;
            b2 &= target;
        }

        if (Type == QUIET_CHECKS)
        {
            Stack* st = pos->st;
            b1 &= attacks_from_pawn(st->ksq, Them);
            b2 &= attacks_from_pawn(st->ksq, Them);

            // Add pawn pushes which give discovered check. This is possible only
            // if the pawn is not on the same file as the enemy king, because we
            // don't generate captures. Note that a possible discovery check
            // promotion has been already generated amongst the captures.
            Bitboard dcCandidatesQuiets = blockers_for_king(pos, Them) & pawnsNotOn7;
            if (dcCandidatesQuiets)
            {
                Bitboard dc1 =
                  shift_bb(Up, dcCandidatesQuiets) & emptySquares & ~file_bb_s(st->ksq);
                Bitboard dc2 = shift_bb(Up, dc1 & TRank3BB) & emptySquares;

                b1 |= dc1;
                b2 |= dc2;
            }
        }

        while (b1)
        {
            Square to      = pop_lsb(&b1);
            (list++)->move = make_move(to - Up, to);
        }

        while (b2)
        {
            Square to      = pop_lsb(&b2);
            (list++)->move = make_move(to - Up - Up, to);
        }
    }

    // Promotions and underpromotions
    if (pawnsOn7 && (Type != EVASIONS || (target & TRank8BB)))
    {
        if (Type == CAPTURES)
            emptySquares = ~pieces();

        if (Type == EVASIONS)
            emptySquares &= target;

        Bitboard b1 = shift_bb(Right, pawnsOn7) & enemies;
        Bitboard b2 = shift_bb(Left, pawnsOn7) & enemies;
        Bitboard b3 = shift_bb(Up, pawnsOn7) & emptySquares;

        while (b1)
            list = make_promotions(list, pop_lsb(&b1), pos->st->ksq, Type, Right);

        while (b2)
            list = make_promotions(list, pop_lsb(&b2), pos->st->ksq, Type, Left);

        while (b3)
            list = make_promotions(list, pop_lsb(&b3), pos->st->ksq, Type, Up);
    }

    // Standard and en-passant captures
    if (Type == CAPTURES || Type == EVASIONS || Type == NON_EVASIONS)
    {
        Bitboard b1 = shift_bb(Right, pawnsNotOn7) & enemies;
        Bitboard b2 = shift_bb(Left, pawnsNotOn7) & enemies;

        while (b1)
        {
            Square to      = pop_lsb(&b1);
            (list++)->move = make_move(to - Right, to);
        }

        while (b2)
        {
            Square to      = pop_lsb(&b2);
            (list++)->move = make_move(to - Left, to);
        }

        if (ep_square() != 0)
        {

            // An en passant capture can be an evasion only if the checking piece
            // is the double pushed pawn and so is in the target. Otherwise this
            // is a discovery check and we are forced to do otherwise.
            if (Type == EVASIONS && !(target & sq_bb(ep_square() - Up)))
                return list;

            b1 = pawnsNotOn7 & attacks_from_pawn(ep_square(), Them);


            while (b1)
                (list++)->move = make_enpassant(pop_lsb(&b1), ep_square());
        }
    }

    return list;
}


INLINE ExtMove* generate_moves(const Position* pos,
                               ExtMove*        list,
                               Bitboard        target,
                               const Color     Us,
                               const int       Pt,
                               const bool      Checks) {

    Square from;

    loop_through_pieces(Us, Pt, from) {
        if (Checks)
        {
            if ((Pt == BISHOP || Pt == ROOK || Pt == QUEEN)
                && !(PseudoAttacks[Pt][from] & target & pos->st->checkSquares[Pt]))
                continue;

            if (blockers_for_king(pos, !Us) & sq_bb(from))
                continue;
        }

        Bitboard b = attacks_from(Pt, from) & target;

        if (Checks)
            b &= pos->st->checkSquares[Pt];

        while (b)
            (list++)->move = make_move(from, pop_lsb(&b));
    }

    return list;
}


INLINE ExtMove*
generate_all(const Position* pos, ExtMove* list, Bitboard target, const Color Us, const int Type) {
    const int  OO     = make_castling_right(Us, KING_SIDE);
    const int  OOO    = make_castling_right(Us, QUEEN_SIDE);
    const bool Checks = Type == QUIET_CHECKS;

    list = generate_pawn_moves(pos, list, target, Us, Type);
    list = generate_moves(pos, list, target, Us, KNIGHT, Checks);
    list = generate_moves(pos, list, target, Us, BISHOP, Checks);
    list = generate_moves(pos, list, target, Us, ROOK, Checks);
    list = generate_moves(pos, list, target, Us, QUEEN, Checks);

    if (Type != QUIET_CHECKS && Type != EVASIONS)
    {
        Square   ksq = square_of(Us, KING);
        Bitboard b   = attacks_from_king(ksq) & target;
        while (b)
            (list++)->move = make_move(ksq, pop_lsb(&b));

        if (Type != CAPTURES && can_castle_c(Us))
        {
            if (!castling_impeded(OO) && can_castle_cr(OO))
                (list++)->move = make_castling(ksq, castling_rook_square(OO));

            if (!castling_impeded(OOO) && can_castle_cr(OOO))
                (list++)->move = make_castling(ksq, castling_rook_square(OOO));
        }
    }

    return list;
}


// generate_captures() generates all pseudo-legal captures plus queen and
// checking knight promotions.
//
// generate_quiets() generates all pseudo-legal non-captures and
// underpromotions (except checking knight promotions).
//
// generate_non_evasions() generates all pseudo-legal captures and
// non-captures.

ExtMove* generate(const Position* pos, ExtMove* list, const int Type) {

    Color us = stm();

    Bitboard target = Type == CAPTURES     ? pieces_c(!us)
                    : Type == QUIETS       ? ~pieces()
                    : Type == NON_EVASIONS ? ~pieces_c(us)
                                           : 0;

    return generate_all(pos, list, target, us, Type);
}


// generate_quiet_checks() generates all pseudo-legal non-captures and
// knight underpromotions that give check.
NOINLINE ExtMove* generate_quiet_checks(const Position* pos, ExtMove* list) {

    Color    us = stm();
    Bitboard dc = blockers_for_king(pos, !us) & pieces_c(us);

    while (dc)
    {
        Square from = pop_lsb(&dc);
        int    pt   = type_of_p(piece_on(from));

        if (pt == PAWN)
            continue;  // Will be generated together with direct checks

        Bitboard b = attacks_from(pt, from) & ~pieces();

        if (pt == KING)
            b &= ~PseudoAttacks[QUEEN][pos->st->ksq];

        while (b)
            (list++)->move = make_move(from, pop_lsb(&b));
    }

    return us == WHITE ? generate_all(pos, list, ~pieces(), WHITE, QUIET_CHECKS)
                       : generate_all(pos, list, ~pieces(), BLACK, QUIET_CHECKS);
}


// generate_evasions() generates all pseudo-legal check evasions when the
// side to move is in check.
NOINLINE ExtMove* generate_evasions(const Position* pos, ExtMove* list) {

    Color    us            = stm();
    Square   ksq           = square_of(us, KING);
    Bitboard sliderAttacks = 0;
    Bitboard sliders       = checkers() & ~pieces_pp(KNIGHT, PAWN);

    // Find all the squares attacked by slider checkers. We will remove them
    // from the king evasions in order to skip known illegal moves, which
    // avoids any useless legality checks later on.
    while (sliders)
    {
        Square checksq = pop_lsb(&sliders);
        sliderAttacks |= LineBB[ksq][checksq] ^ sq_bb(checksq);
    }

    // Generate evasions for king, capture and non capture moves
    Bitboard b = attacks_from_king(ksq) & ~pieces_c(us) & ~sliderAttacks;
    while (b)
        (list++)->move = make_move(ksq, pop_lsb(&b));

    if (more_than_one(checkers()))
        return list;  // Double check, only a king move can save the day

    // Generate blocking evasions or captures of the checking piece
    Square   checksq = lsb(checkers());
    Bitboard target  = between_bb(ksq, checksq) | sq_bb(checksq);

    return us == WHITE ? generate_all(pos, list, target, WHITE, EVASIONS)
                       : generate_all(pos, list, target, BLACK, EVASIONS);
}


// generate_legal() generates all the legal moves in the given position
NOINLINE ExtMove* generate_legal(const Position* pos, ExtMove* list) {
    Color    us     = stm();
    Bitboard pinned = blockers_for_king(pos, us) & pieces_c(us);
    Square   ksq    = square_of(us, KING);
    ExtMove* cur    = list;

    list = checkers() ? generate_evasions(pos, list) : generate(pos, list, NON_EVASIONS);
    while (cur != list)
        if ((pinned || from_sq(cur->move) == ksq || type_of_m(cur->move) == ENPASSANT)
            && !is_legal(pos, cur->move))
            cur->move = (--list)->move;
        else
            ++cur;

    return list;
}
