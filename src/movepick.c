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


#include "movepick.h"

extern int mp_v1;
extern int mp_v2;
extern int mp_v3;
extern int mp_v4;
extern int mp_v5;
extern int mp_v6;
extern int mp_v7;
extern int mp_v8;
extern int mp_v9;
extern int mp_v10;
extern int mp_v11;

// An insertion sort which sorts moves in descending order up to and
// including a given limit. The order of moves smaller than the limit is
// left unspecified.

static void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {
    for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; p++)
        if (p->value >= limit)
        {
            ExtMove tmp = *p, *q;
            *p          = *++sortedEnd;
            for (q = sortedEnd; q != begin && (q - 1)->value < tmp.value; q--)
                *q = *(q - 1);
            *q = tmp;
        }
}


// pick_best() finds the best move in the range (begin, end).

static Move pick_best(ExtMove* begin, ExtMove* end) {
    ExtMove *p, *q;

    for (p = begin, q = begin + 1; q < end; q++)
        if (q->value > p->value)
            p = q;
    Move m       = p->move;
    int  v       = p->value;
    *p           = *begin;
    begin->value = v;

    return m;
}


// score() assigns a numerical value to each move in a move list. The moves with
// highest values will be picked first.

static void score_captures(const Position* pos) {
    Stack*                 st      = pos->st;
    CapturePieceToHistory* history = pos->captureHistory;

    // Winning and equal captures in the main search are ordered by MVV,
    // preferring captures near our with a good history.

    for (ExtMove* m = st->cur; m < st->endMoves; m++)
        m->value =
          PieceValue[piece_on(to_sq(m->move))] * 6
          + (*history)[moved_piece(m->move)][to_sq(m->move)][type_of_p(piece_on(to_sq(m->move)))];
}

SMALL
static void score_quiets(const Position* pos) {
    Stack*            st      = pos->st;
    ButterflyHistory* history = pos->mainHistory;

    PieceToHistory* contHist0 = (st - 1)->continuationHistory;
    PieceToHistory* contHist1 = (st - 2)->continuationHistory;
    PieceToHistory* contHist2 = (st - 4)->continuationHistory;
    PieceToHistory* contHist3 = (st - 6)->continuationHistory;

    Color c = stm();

    for (ExtMove* m = st->cur; m < st->endMoves; m++)
    {
        uint32_t move = m->move & 4095;
        Square   to   = move & 63;
        Square   from = move >> 6;
        m->value =
          (mp_v4 * (*history)[c][move] + mp_v5 * (*contHist0)[piece_on(from)][to]
           + mp_v6 * (*contHist1)[piece_on(from)][to] + mp_v7 * (*contHist2)[piece_on(from)][to]
           + mp_v8 * (*contHist3)[piece_on(from)][to])
          / 100;
    }
}

static void score_evasions(const Position* pos) {
    Stack* st = pos->st;
    // Try captures ordered by MVV/LVA, then non-captures ordered by
    // stats heuristics.

    ButterflyHistory* history   = pos->mainHistory;
    PieceToHistory*   contHist0 = (st - 1)->continuationHistory;
    Color             c         = stm();

    for (ExtMove* m = st->cur; m < st->endMoves; m++)
        if (is_capture(pos, m->move))
            m->value = PieceValue[piece_on(to_sq(m->move))] - type_of_p(moved_piece(m->move));
        else
            m->value = (*history)[c][from_to(m->move)]
                     + (*contHist0)[moved_piece(m->move)][to_sq(m->move)] - (1 << 28);
}


// next_move() returns the next pseudo-legal move to be searched.

Move next_move(const Position* pos, bool skipQuiets) {
    Stack* st = pos->st;
    Move   move;

    switch (st->stage)
    {

    case ST_MAIN_SEARCH :
    case ST_EVASION :
    case ST_QSEARCH :
    case ST_PROBCUT :
        st->endMoves = (st - 1)->endMoves;
        st->stage++;
        return st->ttMove;

    case ST_CAPTURES_INIT :
        st->endBadCaptures = st->cur = (st - 1)->endMoves;
        st->endMoves                 = generate(pos, st->cur, CAPTURES);
        score_captures(pos);
        st->stage++;
        /* fallthrough */

    case ST_GOOD_CAPTURES :
        while (st->cur < st->endMoves)
        {
            move = pick_best(st->cur++, st->endMoves);
            if (move != st->ttMove)
            {
                if (see_test(pos, move, -mp_v1 * (st->cur - 1)->value / mp_v2))
                    return move;

                // Losing capture, move it to the beginning of the array.
                (st->endBadCaptures++)->move = move;
            }
        }
        st->stage++;

        // First killer move.
        move = st->mpKillers[0];
        if (move && move != st->ttMove && is_pseudo_legal(pos, move) && !is_capture(pos, move))
            return move;
        /* fallthrough */

    case ST_KILLERS :
        st->stage++;
        move = st->mpKillers[1];  // Second killer move.
        if (move && move != st->ttMove && is_pseudo_legal(pos, move) && !is_capture(pos, move))
            return move;
        /* fallthrough */

    case ST_KILLERS_2 :
        st->stage++;
        move = st->countermove;
        if (move && move != st->ttMove && move != st->mpKillers[0] && move != st->mpKillers[1]
            && is_pseudo_legal(pos, move) && !is_capture(pos, move))
            return move;
        /* fallthrough */

    case ST_QUIET_INIT :
        if (!skipQuiets)
        {
            st->cur      = st->endBadCaptures;
            st->endMoves = generate(pos, st->cur, QUIETS);
            score_quiets(pos);
            partial_insertion_sort(st->cur, st->endMoves, -mp_v3 * st->depth);
        }
        st->stage++;
        /* fallthrough */

    case ST_QUIET :
        if (!skipQuiets)
            while (st->cur < st->endMoves)
            {
                move = (st->cur++)->move;
                if (move != st->ttMove && move != st->mpKillers[0] && move != st->mpKillers[1]
                    && move != st->countermove)
                    return move;
            }
        st->stage++;
        st->cur = (st - 1)->endMoves;  // Return to bad captures.
                                       /* fallthrough */

    case ST_BAD_CAPTURES :
        if (st->cur < st->endBadCaptures)
            return (st->cur++)->move;
        break;

    case ST_EVASIONS_INIT :
        st->cur      = (st - 1)->endMoves;
        st->endMoves = generate(pos, st->cur, EVASIONS);
        score_evasions(pos);
        st->stage++;

    case ST_ALL_EVASIONS :
        while (st->cur < st->endMoves)
        {
            move = pick_best(st->cur++, st->endMoves);
            if (move != st->ttMove)
                return move;
        }
        break;

    case ST_QCAPTURES_INIT :
        st->cur      = (st - 1)->endMoves;
        st->endMoves = generate(pos, st->cur, CAPTURES);
        score_captures(pos);
        st->stage++;

    case ST_QCAPTURES :
        while (st->cur < st->endMoves)
        {
            move = pick_best(st->cur++, st->endMoves);
            if (move != st->ttMove
                && (st->depth > DEPTH_QS_RECAPTURES || to_sq(move) == st->recaptureSquare))
                return move;
        }
        if (st->depth <= DEPTH_QS_NO_CHECKS)
            break;
        st->cur      = (st - 1)->endMoves;
        st->endMoves = generate(pos, st->cur, QUIET_CHECKS);
        st->stage++;
        /* fallthrough */

    case ST_QCHECKS :
        while (st->cur < st->endMoves)
        {
            move = (st->cur++)->move;
            if (move != st->ttMove)
                return move;
        }
        break;

    case ST_PROBCUT_INIT :
        st->cur      = (st - 1)->endMoves;
        st->endMoves = generate(pos, st->cur, CAPTURES);
        score_captures(pos);
        st->stage++;
        /* fallthrough */

    case ST_PROBCUT_2 :
        while (st->cur < st->endMoves)
        {
            move = pick_best(st->cur++, st->endMoves);
            if (move != st->ttMove && see_test(pos, move, st->threshold))
                return move;
        }
        break;

    default :;
    }

    return 0;
}
