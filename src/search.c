/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2020 The Stockfish developers

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


#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#include "evaluate.h"
#include "misc.h"
#include "movegen.h"
#include "movepick.h"
#include "search.h"
#include "settings.h"
#include "timeman.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"

#ifndef KAGGLE
Parameter parameters[255];
int       parameters_count = 0;

    #define STRINGIFY(x) #x
    #define PARAM(Name, Value) \
        int                                               Name = Value; \
        __attribute__((constructor(__LINE__ + 100))) void ctor_##Name(void) { \
            Parameter* param = &parameters[parameters_count++]; \
            strncpy(param->name, STRINGIFY(Name), sizeof(param->name)); \
            param->value = &Name; \
        }
#else
    #define PARAM(Name, Value) int Name = Value;
#endif


PARAM(nmp_v1, 764)
PARAM(nmp_v2, 56)
PARAM(nmp_v3, 165)
PARAM(nmp_v4, 182)
PARAM(nmp_v5, 25312)
PARAM(nmp_v6, 27)
PARAM(nmp_v7, 23)
PARAM(nmp_v8, 92)
PARAM(nmp_v9, 194)
PARAM(qmo_v1, 357)
PARAM(qmo_v2, 1118)
PARAM(qmo_v3, 1113)
PARAM(rz_v1, 565)
PARAM(rz_v2, 115)
PARAM(ft_v1, 210)
PARAM(rd_v1, 578)
PARAM(rd_v2, 1251)
PARAM(rd_v3, 854)
PARAM(rd_init_v1, 2857)
PARAM(d_v1, 18)
PARAM(iir_v1, 670)
PARAM(iir_v2, 238)
PARAM(cbp_v1, 353)
PARAM(cbp_v2, 6)
PARAM(fpp_v1, 706)
PARAM(fpp_v2, 217)
PARAM(fpp_v3, 184)
PARAM(sqsee_v1, 2856)
PARAM(sfpc_v1, 630)
PARAM(sfpc_v2, 197)
PARAM(sfpc_v3, 183)
PARAM(sfpc_v4, 230)
PARAM(scsee_v1, 199)
PARAM(se_v1, 575)
PARAM(se_v2, 177)
PARAM(se_v5, 3755)
PARAM(se_v6, 112)
PARAM(se_v7, 240)
PARAM(se_v8, 257)
PARAM(prb_v1, 134)
PARAM(prb_v2, 45)
PARAM(rfp_v1, 933)
PARAM(lmr_v3, 3792)
PARAM(lmr_v4, 99)
PARAM(lmr_v5, 88)
PARAM(lmr_v6, 96)
PARAM(lmr_v7, 106)
PARAM(lmr_v8, 14626)
PARAM(fmc_v1, 305)
PARAM(fmc_v2, 298)
PARAM(fmc_v3, 238)
PARAM(hb_v1, 690)
PARAM(hb_v2, 198)
PARAM(hb_v3, 155)
PARAM(hb_v4, 2362)
PARAM(hm_v1, 699)
PARAM(hm_v2, 187)
PARAM(hm_v3, 167)
PARAM(hm_v4, 1581)
PARAM(asd_v1, 493)
PARAM(ses_v1, 396)
PARAM(qsf_v1, 177)
PARAM(ch_v1, 1849)
PARAM(ch_v2, 180)
PARAM(ch_v3, 262)
PARAM(ch_v4, 86)
PARAM(ch_v5, 101)
PARAM(tempo, 47)
PARAM(mp_v1, 70)
PARAM(mp_v2, 1064)
PARAM(mp_v3, 2839)
PARAM(mp_v4, 126)
PARAM(mp_v5, 242)
PARAM(mp_v6, 236)
PARAM(mp_v7, 200)
PARAM(mp_v8, 99)
PARAM(mp_v9, 107)
PARAM(eval_scale, 90)
PARAM(mat_scale, 26211)
PARAM(pcmb_v1, 113)
PARAM(pcmb_v2, 484)
PARAM(pcmb_v3, 38)
PARAM(pcmb_v4, 165)
PARAM(pcmb_v5, 786)
PARAM(pcmb_v6, 117)
PARAM(pcmb_v7, 101)
PARAM(pcmb_v8, 131)
PARAM(pcmb_v9, 232)
PARAM(pcmb_v10, 110)
PARAM(pcmb_v11, 146)
PARAM(r_v1, 1039)
PARAM(r_v2, 2278)
PARAM(r_v3, 861)
PARAM(r_v4, 1033)
PARAM(r_v5, 1094)
PARAM(r_v6, 1158)
PARAM(r_v7, 1056)
PARAM(r_v8, 2101)
PARAM(r_v9, 1852)
PARAM(r_v10, 872)
PARAM(r_v11, 924)
PARAM(r_v12, 1145)
PARAM(r_v13, 1016)
PARAM(lce_v1, 2280)
PARAM(qb_v1, 199)
PARAM(cms_v1, 29304)
PARAM(hu_v1, 10217)
PARAM(cpth_v1, 11598)
PARAM(kb_v1, 29)

PARAM(tm_v1, 329)
PARAM(tm_v2, 555)
PARAM(tm_v3, 657)
PARAM(tm_v4, 809)
PARAM(tm_v5, 52)
PARAM(tm_v6, 155)
PARAM(tm_v7, 869)
PARAM(tm_v8, 187)
PARAM(tm_v9, 101)
PARAM(tm_v10, 147)
PARAM(tm_v11, 216)
PARAM(tm_v12, 57)
PARAM(tm_v13, 89)
PARAM(tm_v14, 296)
PARAM(tm_v15, 248)
PARAM(tm_v16, 21)
PARAM(tm_v17, 674)
PARAM(tm_v18, 416)
PARAM(tm_v19, 1186)
PARAM(tm_v20, 801)
PARAM(tm_v21, 1069)
PARAM(tm_v22, 225)
PARAM(tm_v23, 990)
PARAM(mtg, 50)

LimitsType Limits;

enum {
    NonPV,
    PV
};

static int futility_margin(Depth d, bool improving) { return ft_v1 * (d - improving); }

// Reductions lookup tables, initialized at startup
static int Reductions[MAX_MOVES];  // [depth or moveNumber]

static Depth reduction(int i, Depth d, int mn) {
    int r = Reductions[d] * Reductions[mn];
    return (r + rd_v1) / rd_v2 + (!i && r > rd_v3);
}

static int futility_move_count(bool improving, Depth depth) {
    //  return (3 + depth * depth) / (2 - improving);
    return improving ? fmc_v1 / 100 + depth * depth
                     : (fmc_v2 / 100 + depth * depth) / (fmc_v3 / 100);
}

// History and stats update bonus, based on depth
static Value stat_bonus(Depth d) { return min((hb_v1 / 100 * d + hb_v2) * d - hb_v3, hb_v4); }

// History and stats update malus, based on depth
static Value stat_malus(Depth d) { return min((hm_v1 / 100 * d + hm_v2) * d - hm_v3, hm_v4); }

// Add a small random component to draw evaluations to keep search dynamic
// and to avoid three-fold blindness. (Yucks, ugly hack)
static Value value_draw(Position* pos) { return VALUE_DRAW + 2 * (pos->nodes & 1) - 1; }

static Value value_to_tt(Value v, int ply);
static Value value_from_tt(Value v, int ply, int r50c);
static void  update_pv(Move* pv, Move move, Move* childPv);
static void  update_continuation_histories(Stack* ss, Piece pc, Square s, int bonus);
Value        to_corrected(Position* pos, Value rawEval);
static void
add_correction_history(CorrectionHistory hist, Color side, Key key, Depth depth, int32_t diff);
static void update_quiet_stats(const Position* pos, Stack* ss, Move move, int bonus);
static void
update_capture_stats(const Position* pos, Move move, Move* captures, int captureCnt, int bonus);
static void check_time(void);
static void stable_sort(RootMove* rm, int num);
static void uci_print_pv(Position* pos, Depth depth, Value alpha, Value beta);
static int  extract_ponder_from_tt(RootMove* rm, Position* pos);

SMALL double my_log(double x) {
    double result = 0.0;
    while (x >= 2.0)
    {
        x /= 2.718281828459045;  // Divide by e
        result += 1.0;
    }
    x -= 1.0;
    double term = x, sum = 0.0;
    for (int n = 1; n < 10; n++)
    {
        sum += (n % 2 ? 1 : -1) * term / n;  // Alternating series
        term *= x;
    }
    return result + sum;
}

// search_init() is called during startup to initialize various lookup tables

SMALL void search_init(void) {
    for (int i = 1; i < MAX_MOVES; i++)
        Reductions[i] = rd_init_v1 / 100.0 * my_log(i);
}


// search_clear() resets search state to zero, to obtain reproducible results

SMALL void search_clear(void) {
    if (!settings.ttSize)
    {
        delayedSettings.clear = true;
        return;
    }

    Time.availableNodes = 0;

    tt_clear();

    Position* pos = Thread.pos;
    stats_clear(pos->mainHistory);
    stats_clear(pos->captureHistory);
    stats_clear(pos->contHist);
    stats_clear(pos->matCorrHist);
    stats_clear(pos->pawnCorrHist);

    for (int c = 0; c < 2; c++)
        for (int j = 0; j < 16; j++)
            for (int k = 0; k < 64; k++)
                (*pos->contHist)[c][0][j][k] = CounterMovePruneThreshold - 1;

    Thread.previousScore         = VALUE_INFINITE;
    Thread.previousTimeReduction = 1;
}

// mainthread_search() is called by the main thread when the program
// receives the UCI 'go' command. It searches from the root position and
// outputs the "bestmove".

void mainthread_search(void) {
    Position* pos = Thread.pos;
    Color     us  = stm();
    time_init(us, game_ply());
    tt_new_search();
    char buf[16];

    Thread.pos->bestMoveChanges = 0;
    thread_search(pos);
    Thread.previousScore = pos->rootMoves->move[0].score;

    printf("bestmove %s\n", uci_move(buf, pos->rootMoves->move[0].pv[0]));
    fflush(stdout);

    if (!IsKaggle && !Thread.testPonder)
        return;

    // Start pondering right after the best move has been printed if we can
    if (pos->rootMoves->move[0].pvSize >= 2
        || extract_ponder_from_tt(&pos->rootMoves->move[0], pos))
    {
        Thread.ponder = true;
        Thread.stop   = false;

        const Move bestMove = pos->rootMoves->move[0].pv[0];
        const Move ponder   = pos->rootMoves->move[0].pv[1];

        do_move(pos, bestMove, gives_check(pos, pos->st, bestMove));
        do_move(pos, ponder, gives_check(pos, pos->st, ponder));

        pos->completedDepth = 0;
        pos->rootDepth      = 0;
        pos->pvLast         = 0;

        prepare_for_search(pos);
        thread_search(pos);

        Thread.ponder = false;
        Thread.stop   = true;
    }
}


// thread_search() is the main iterative deepening loop. It calls search()
// repeatedly with increasing depth until the allocated thinking time has
// been consumed, the user stops the search, or the maximum search depth is
// reached.

void thread_search(Position* pos) {
    Value  bestValue, alpha, beta, delta;
    Move   pv[3];
    Move   lastBestMove      = 0;
    Depth  lastBestMoveDepth = 0;
    double timeReduction = 1.0, totBestMoveChanges = 0;
    int    iterIdx = 0;

    Stack* ss = pos->st;  // At least the seventh element of the allocated array.
    for (int i = -7; i < 3; i++)
    {
        memset(SStackBegin(ss[i]), 0, SStackSize);
    }
    (ss - 1)->endMoves = pos->moveList;

    for (int i = -7; i < 0; i++)
    {
        ss[i].continuationHistory = &(*pos->contHist)[0][0];  // Use as sentinel
        ss[i].staticEval          = VALUE_NONE;
        ss[i].checkersBB          = 0;
    }

    for (int i = 0; i <= MAX_PLY; i++)
        ss[i].ply = i;
    ss->pv = pv;

    ss->accumulator.needs_refresh = 1;

    bestValue = delta = alpha = -VALUE_INFINITE;
    beta                      = VALUE_INFINITE;
    pos->completedDepth       = 0;

    int value = Thread.previousScore == VALUE_INFINITE ? VALUE_ZERO : Thread.previousScore;
    for (int i = 0; i < 4; i++)
        Thread.iterValue[i] = value;

    RootMoves* rm = pos->rootMoves;

    // Iterative deepening loop until requested to stop or the target depth
    // is reached.
    while ((pos->rootDepth += 2) < MAX_PLY && !Thread.stop
           && !(Limits.depth && pos->rootDepth > Limits.depth))
    {
        // Age out PV variability metric
        totBestMoveChanges /= 2;

        // Save the last iteration's scores before first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
        for (int idx = 0; idx < rm->size; idx++)
            rm->move[idx].previousScore = rm->move[idx].score;

        pos->pvLast = rm->size;

        // Reset aspiration window starting size
        if (pos->rootDepth >= 4)
        {
            Value previousScore = rm->move[0].previousScore;
            delta               = d_v1;
            alpha               = max(previousScore - delta, -VALUE_INFINITE);
            beta                = min(previousScore + delta, VALUE_INFINITE);
        }

        // Start with a small aspiration window and, in the case of a fail
        // high/low, re-search with a bigger window until we're not failing
        // high/low anymore.
        while (true)
        {
            Depth adjustedDepth = max(1, pos->rootDepth);
            bestValue           = search(pos, ss, alpha, beta, adjustedDepth, false, true);

            // Bring the best move to the front. It is critical that sorting
            // is done with a stable algorithm because all the values but the
            // first and eventually the new best one are set to -VALUE_INFINITE
            // and we want to keep the same order for all the moves except the
            // new PV that goes to the front. Note that in case of MultiPV
            // search the already searched PV lines are preserved.
            stable_sort(&rm->move[0], rm->size);

            // If search has been stopped, we break immediately. Sorting and
            // writing PV back to TT is safe because RootMoves is still
            // valid, although it refers to the previous iteration.
            if (Thread.stop)
                break;

            // In case of failing low/high increase aspiration window and
            // re-search, otherwise exit the loop.
            if (bestValue <= alpha)
            {
                beta  = (alpha + beta) / 2;
                alpha = max(bestValue - delta, -VALUE_INFINITE);
            }
            else
                break;

            delta += delta / 4 + asd_v1 / 100;
        }

#ifndef KAGGLE
        if (!Thread.ponder)
            uci_print_pv(pos, pos->rootDepth, alpha, beta);
#endif

        if (!Thread.stop)
            pos->completedDepth = pos->rootDepth;

        if (rm->move[0].pv[0] != lastBestMove)
        {
            lastBestMove      = rm->move[0].pv[0];
            lastBestMoveDepth = pos->rootDepth;
        }

// Do we have time for the next iteration? Can we stop searching now?
#ifndef KAGGLE
        if (use_time_management() && !Thread.stop)
#else
        if (!Thread.stop)
#endif
        {
            double fallingEval = (tm_v1 + tm_v2 / 100.0 * (Thread.previousScore - bestValue)
                                  + tm_v3 / 100.0 * (Thread.iterValue[iterIdx] - bestValue))
                               / (double) tm_v4;
            fallingEval = clamp(fallingEval, tm_v5 / 100.0, tm_v6 / 100.0);

            // If the best move is stable over several iterations, reduce time
            // accordingly
            timeReduction =
              lastBestMoveDepth + tm_v7 / 100 < pos->completedDepth ? tm_v8 / 100.0 : tm_v9 / 100.0;
            double reduction =
              (tm_v10 / 100.0 + Thread.previousTimeReduction) / (tm_v11 / 100.0 * timeReduction);

            // Use part of the gained time from a previous stable move for this move
            totBestMoveChanges += Thread.pos->bestMoveChanges;
            Thread.pos->bestMoveChanges = 0;

            double bestMoveInstability =
              tm_v21 / 1000.0
              + max(1.0, tm_v22 / 100.0 - tm_v23 / 100.0 / (pos->rootDepth)) * totBestMoveChanges;

            double totalTime =
              rm->size == 1 ? 0 : time_optimum() * fallingEval * reduction * bestMoveInstability;

            // Stop the search if we have exceeded the totalTime (at least 1ms)
            if (time_elapsed() > totalTime)
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "stop".
                if (!Thread.ponder)
                    Thread.stop = true;
            }
            else
                Thread.increaseDepth = !(Thread.increaseDepth && !Thread.ponder
                                         && time_elapsed() > totalTime * tm_v12 / 100.0);
        }

        Thread.iterValue[iterIdx] = bestValue;
        iterIdx                   = (iterIdx + 1) & 3;
    }

    Thread.previousTimeReduction = timeReduction;
}

// search() is the main search function template for both PV
// and non-PV nodes
Value search(
  Position* pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode, const int NT) {
    const bool PvNode   = NT == PV;
    const bool rootNode = PvNode && ss->ply == 0;

    // Dive into quiescense search when the depth reaches zero
    if (depth <= 0)
        return qsearch(pos, ss, alpha, beta, 0, PvNode, checkers());

    Move     pv[3], capturesSearched[32], quietsSearched[32];
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, excludedMove, bestMove;
    Depth    extension, newDepth;
    Value    bestValue, value, ttValue, eval, rawEval, probCutBeta;
    bool     ttHit, givesCheck, improving;
    bool     captureOrPromotion, inCheck, moveCountPruning;
    bool     ttCapture;
    Piece    movedPiece;
    int      moveCount, captureCount, quietCount;

    // Step 1. Initialize node
    inCheck   = checkers();
    moveCount = captureCount = quietCount = ss->moveCount = 0;
    bestValue                                             = -VALUE_INFINITE;

    // Check for the available remaining time
    if (--pos->callsCnt <= 0)
    {
        pos->callsCnt = 1024;
        check_time();
    }

    if (!rootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (Thread.stop || is_draw(pos) || ss->ply >= MAX_PLY)
            return ss->ply >= MAX_PLY && !inCheck ? evaluate(pos) : value_draw(pos);
    }

    (ss + 1)->ttPv         = false;
    (ss + 1)->excludedMove = bestMove = 0;
    (ss + 2)->killers[0] = (ss + 2)->killers[1] = 0;
    (ss + 2)->cutoffCnt                         = 0;
    Square prevSq = move_is_ok((ss - 1)->currentMove) ? to_sq((ss - 1)->currentMove) : SQ_NONE;

    // Initialize statScore to zero for the grandchildren of the current
    // position. So the statScore is shared between all grandchildren and only
    // the first grandchild starts with startScore = 0. Later grandchildren
    // start with the last calculated statScore of the previous grandchild.
    // This influences the reduction rules in LMR which are based on the
    // statScore of the parent position.
    if (rootNode)
        (ss + 4)->statScore = 0;
    else
        (ss + 2)->statScore = 0;

    // Step 4. Transposition table lookup. We don't want the score of a
    // partial search to overwrite a previous full search TT value, so we
    // use a different position key in case of an excluded move.
    excludedMove = ss->excludedMove;
    posKey       = !excludedMove ? key() : key() ^ make_key(excludedMove);
    tte          = tt_probe(posKey, &ttHit);
    ttValue      = ttHit ? value_from_tt(tte_value(tte), ss->ply, rule50_count()) : VALUE_NONE;
    ttMove       = rootNode ? pos->rootMoves->move[0].pv[0] : ttHit ? tte_move(tte) : 0;
    if (!excludedMove)
        ss->ttPv = PvNode || (ttHit && tte_is_pv(tte));

    // At non-PV nodes we check for an early TT cutoff.
    if (!PvNode && ttHit && tte_depth(tte) >= depth
        && ttValue != VALUE_NONE  // Possible in case of TT access race.
        && (ttValue >= beta ? (tte_bound(tte) & BOUND_LOWER) : (tte_bound(tte) & BOUND_UPPER)))
    {
        // If ttMove is quiet, update move sorting heuristics on TT hit.
        if (ttMove)
        {
            if (ttValue >= beta)
            {
                if (!is_capture_or_promotion(pos, ttMove))
                    update_quiet_stats(pos, ss, ttMove, stat_bonus(depth));

                // Extra penalty for early quiet moves of the previous ply
                if ((ss - 1)->moveCount <= 2 && !captured_piece() && prevSq != SQ_NONE)
                    update_continuation_histories(ss - 1, piece_on(prevSq), prevSq,
                                                  -stat_malus(depth + 1));
            }
            // Penalty for a quiet ttMove that fails low
            else if (!is_capture_or_promotion(pos, ttMove))
            {
                int penalty = -stat_malus(depth);
                history_update(*pos->mainHistory, stm(), ttMove, penalty);
                update_continuation_histories(ss, moved_piece(ttMove), to_sq(ttMove), penalty);
            }
        }
        if (rule50_count() < 90)
            return ttValue;
    }

    // Step 6. Static evaluation of the position
    if (inCheck)
    {
        // Skip early pruning when in check
        ss->staticEval = eval = rawEval = VALUE_NONE;
        improving                       = false;
        goto moves_loop;
    }
    else if (ttHit)
    {
        // Never assume anything about values stored in TT
        if ((rawEval = tte_eval(tte)) == VALUE_NONE)
            rawEval = evaluate(pos);

        eval = ss->staticEval = to_corrected(pos, rawEval);

        if (eval == VALUE_DRAW)
            eval = value_draw(pos);

        // Can ttValue be used as a better position evaluation?
        if (ttValue != VALUE_NONE
            && (tte_bound(tte) & (ttValue > eval ? BOUND_LOWER : BOUND_UPPER)))
            eval = ttValue;
    }
    else
    {
        if ((ss - 1)->currentMove != MOVE_NULL)
            rawEval = evaluate(pos);
        else
            rawEval = -(ss - 1)->staticEval + tempo;

        eval = ss->staticEval = to_corrected(pos, rawEval);

        tte_save(tte, posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_NONE, 0, rawEval);
    }

    // Step 7. Razoring
    if (!rootNode && depth <= rz_v2 / 100 && eval <= alpha - rz_v1)
        return qsearch(pos, ss, alpha, beta, 0, PvNode, false);


    improving = (ss - 2)->staticEval == VALUE_NONE
                ? (ss->staticEval > (ss - 4)->staticEval || (ss - 4)->staticEval == VALUE_NONE)
                : ss->staticEval > (ss - 2)->staticEval;

    if (prevSq != SQ_NONE && !(ss - 1)->checkersBB && !captured_piece())
    {
        int bonus = clamp(-depth * qmo_v1 / 100 * ((ss - 1)->staticEval + ss->staticEval - tempo),
                          -qmo_v2, qmo_v3);
        history_update(*pos->mainHistory, !stm(), (ss - 1)->currentMove, bonus);
    }

    // Step 8. Futility pruning: child node
    if (!PvNode && depth < rfp_v1 / 100 && eval - futility_margin(depth, improving) >= beta
        && eval < VALUE_KNOWN_WIN)  // Do not return unproven wins
        return eval;                // - futility_margin(depth); (do not do the right thing)

    // Step 9. Null move search
    if (!PvNode && (ss - 1)->currentMove != MOVE_NULL && (ss - 1)->statScore < nmp_v5
        && eval >= beta && eval >= ss->staticEval
        && ss->staticEval >= beta - nmp_v6 * depth - nmp_v7 * improving + nmp_v8 * ss->ttPv + nmp_v9
        && !excludedMove && non_pawn_material_c(stm()))
    {
        // Null move dynamic reduction based on depth and value
        Depth R = (nmp_v1 + nmp_v2 * depth) / nmp_v3 + min((eval - beta) / nmp_v4, 3);

        ss->currentMove         = MOVE_NULL;
        ss->continuationHistory = &(*pos->contHist)[0][0];

        do_null_move(pos);
        ss->endMoves    = (ss - 1)->endMoves;
        Value nullValue = -search(pos, ss + 1, -beta, -beta + 1, depth - R, !cutNode, false);
        undo_null_move(pos);

        if (nullValue >= beta)
            return nullValue > VALUE_MATE_IN_MAX_PLY ? beta : nullValue;
    }

    probCutBeta = beta + prb_v1 - prb_v2 * improving;

    // Step 10. ProbCut
    // If we have a good enough capture and a reduced search returns a value
    // much above beta, we can (almost) safely prune the previous move.
    if (!PvNode && depth > 4 && abs(beta) < VALUE_MATE_IN_MAX_PLY
        && !(ttHit && tte_depth(tte) >= depth - 3 && ttValue != VALUE_NONE
             && ttValue < probCutBeta))
    {
        if (ttHit && tte_depth(tte) >= depth - 3 && ttValue != VALUE_NONE && ttValue >= probCutBeta
            && ttMove && is_capture_or_promotion(pos, ttMove))
            return probCutBeta;

        mp_init_pc(pos, ttMove, probCutBeta - ss->staticEval);
        int  probCutCount = 2 + 2 * cutNode;
        bool ttPv         = ss->ttPv;
        ss->ttPv          = false;

        while ((move = next_move(pos, 0)) && probCutCount)
            if (move != excludedMove && is_legal(pos, move))
            {
                captureOrPromotion = true;
                probCutCount--;

                ss->currentMove         = move;
                ss->continuationHistory = &(*pos->contHist)[moved_piece(move)][to_sq(move)];
                givesCheck              = gives_check(pos, ss, move);
                do_move(pos, move, givesCheck);

                // Perform a preliminary qsearch to verify that the move holds
                value = -qsearch(pos, ss + 1, -probCutBeta, -probCutBeta + 1, 0, false, givesCheck);

                // If the qsearch held, perform the regular search
                if (value >= probCutBeta)
                    value = -search(pos, ss + 1, -probCutBeta, -probCutBeta + 1, depth - 4,
                                    !cutNode, false);
                undo_move(pos, move);
                if (value >= probCutBeta)
                {
                    if (!(ttHit && tte_depth(tte) >= depth - 3 && ttValue != VALUE_NONE))
                        tte_save(tte, posKey, value_to_tt(value, ss->ply), ttPv, BOUND_LOWER,
                                 depth - 3, move, rawEval);
                    return value;
                }
            }
        ss->ttPv = ttPv;
    }

    // Step 11. If the position is not in TT, decrease depth by 2
    if (PvNode && depth >= max(iir_v1 / 100, 3) && !ttMove)
        depth -= max(iir_v2 / 100, 2);

moves_loop:  // When in check search starts from here.
  ;          // Avoid a compiler warning. A label must be followed by a statement.
    PieceToHistory* contHist0 = (ss - 1)->continuationHistory;
    PieceToHistory* contHist1 = (ss - 2)->continuationHistory;
    PieceToHistory* contHist2 = (ss - 4)->continuationHistory;

    mp_init(pos, ttMove, depth, ss->ply);

    value            = bestValue;
    moveCountPruning = false;
    ttCapture        = ttMove && is_capture_or_promotion(pos, ttMove);

    // Step 12. Loop through moves
    // Loop through all pseudo-legal moves until no moves remain or a beta
    // cutoff occurs
    while ((move = next_move(pos, moveCountPruning)))
    {
        if (move == excludedMove)
            continue;

        // Check for legality just before making the move
        if (!is_legal(pos, move))
            continue;

        ss->moveCount = ++moveCount;

        if (PvNode && ss->ply <= 2)
            (ss + 1)->pv = NULL;

        extension          = 0;
        captureOrPromotion = is_capture_or_promotion(pos, move);
        movedPiece         = moved_piece(move);

        givesCheck = gives_check(pos, ss, move);

        // Calculate new depth for this move
        newDepth = depth - 1;

        Depth r = reduction(improving, depth, moveCount);

        // Step 13. Pruning at shallow depth
        if (!rootNode && non_pawn_material_c(stm()) && bestValue > VALUE_MATED_IN_MAX_PLY)
        {
            // Skip quiet moves if movecount exceeds our FutilityMoveCount threshold
            moveCountPruning = moveCount >= futility_move_count(improving, depth);

            // Reduced depth of the next LMR search
            int lmrDepth = max(newDepth - r, 0);

            if (!captureOrPromotion && !givesCheck)
            {
                // Countermoves based pruning
                if (lmrDepth < cbp_v1 / 100
                                 + ((ss - 1)->statScore > cbp_v2 / 100 || (ss - 1)->moveCount == 1)
                    && (*contHist0)[movedPiece][to_sq(move)] < CounterMovePruneThreshold
                    && (*contHist1)[movedPiece][to_sq(move)] < CounterMovePruneThreshold)
                    continue;

                // Futility pruning: parent node
                if (lmrDepth < fpp_v1 / 100 && !inCheck
                    && ss->staticEval + fpp_v2 + fpp_v3 * lmrDepth <= alpha)
                    continue;

                // Prune moves with negative SEE at low depths and below a decreasing
                // threshold at higher depths.
                if (!see_test(pos, move, -(sqsee_v1 / 100 * lmrDepth * lmrDepth)))
                    continue;
            }
            else
            {
                // Futility pruning for captures
                if (!givesCheck && lmrDepth < sfpc_v1 / 100
                    && !(PvNode && abs(bestValue) < sfpc_v2 / 100) && !inCheck
                    && ss->staticEval + sfpc_v3 + sfpc_v4 * lmrDepth
                           + PieceValue[type_of_p(piece_on(to_sq(move)))]
                         <= alpha)
                    continue;

                // See based pruning
                if (!see_test(pos, move, -scsee_v1 * depth))
                    continue;
            }
        }

        // Step 14. Extensions

        // Singular extension search. If all moves but one fail low on a search
        // of (alpha-s, beta-s), and just one fails high on (alpha, beta), then
        // that move is singular and should be extended. To verify this we do a
        // reduced search on all the other moves but the ttMove and if the
        // result is lower than ttValue minus a margin, then we extend the ttMove.
        if (depth >= se_v1 / 100 && move == ttMove && !rootNode
            && !excludedMove  // No recursive singular search
                              /* &&  ttValue != VALUE_NONE implicit in the next condition */
            && abs(ttValue) < VALUE_KNOWN_WIN && (tte_bound(tte) & BOUND_LOWER)
            && tte_depth(tte) >= depth - ses_v1 / 100)
        {
            Value singularBeta  = ttValue - se_v2 / 100 * depth;
            Depth singularDepth = (depth - 1) / 2;
            ss->excludedMove    = move;
            Move k1 = ss->mpKillers[0], k2 = ss->mpKillers[1];
            value = search(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode, false);
            ss->excludedMove = 0;

            if (value < singularBeta)
            {
                extension = 1;
                if (!PvNode && value < singularBeta - se_v5 / 100)
                    extension = 2;
            }

            // Multi-cut pruning. Our ttMove is assumed to fail high, and now we
            // failed high also on a reduced search without the ttMove. So we
            // assume that this expected cut-node is not singular, i.e. multiple
            // moves fail high. We therefore prune the whole subtree by returning
            // a soft bound.
            else if (singularBeta >= beta)
            {
                return singularBeta;
            }

            // If the eval of ttMove is greater than beta we also check whether
            // there is another move that pushes it over beta. If so, we prune.
            else if (ttValue >= beta)
            {
                // Fix up our move picker data
                mp_init(pos, ttMove, depth, ss->ply);
                ss->stage++;
                ss->mpKillers[0] = k1;
                ss->mpKillers[1] = k2;

                ss->excludedMove = move;
                value = search(pos, ss, beta - 1, beta, (depth + se_v7 / 100) / (se_v8 / 100),
                               cutNode, false);
                ss->excludedMove = 0;

                if (value >= beta)
                {
                    return beta;
                }
            }
            else if (cutNode)
            {
                extension -= se_v6 / 100;
            }

            // The call to search_NonPV with the same value of ss messed up our
            // move picker data. So we fix it.
            mp_init(pos, ttMove, depth, ss->ply);
            ss->stage++;
            ss->mpKillers[0] = k1;
            ss->mpKillers[1] = k2;
        }

        // Last capture extension
        else if (PieceValue[captured_piece()] > PawnValue && non_pawn_material() <= lce_v1)
            extension = 1;

        // Late irreversible move extension
        if (move == ttMove && rule50_count() > 80
            && (captureOrPromotion || type_of_p(movedPiece) == PAWN))
            extension = 2;

        // Add extension to new depth
        newDepth += extension;

        // Speculative prefetch as early as possible
        prefetch(tt_first_entry(key_after(pos, move)));

        // Update the current move (this must be done after singular extension
        // search)
        ss->currentMove         = move;
        ss->continuationHistory = &(*pos->contHist)[movedPiece][to_sq(move)];

        // Step 15. Make the move.
        do_move(pos, move, givesCheck);

        r = r * r_v1;

        // Step 16. Reduced depth search (LMR). If the move fails high it will be
        // re-searched at full depth.
        if (depth >= 2 && moveCount > 1 + 2 * rootNode
            && (!captureOrPromotion || cutNode || !ss->ttPv))
        {

            // Decrease reduction if position is or has been on the PV
            if (ss->ttPv)
                r -= r_v2;

            if (!captureOrPromotion)
            {
                // Increase reduction if ttMove is a capture
                if (ttCapture)
                    r += r_v6;

                if ((ss + 1)->cutoffCnt > 3)
                    r += r_v7;

                // Increase reduction for cut nodes
                if (cutNode)
                    r += r_v8;

                // Decrease reduction for moves that escape a capture. Filter out
                // castling moves, because they are coded as "king captures rook" and
                // hence break make_move().
                else if (type_of_m(move) == NORMAL && !see_test(pos, reverse_move(move), 0))
                    r -= r_v9 + r_v10 * (ss->ttPv - (type_of_p(movedPiece) == PAWN));

                ss->statScore = (*contHist0)[movedPiece][to_sq(move)]
                              + (*contHist1)[movedPiece][to_sq(move)]
                              + (*contHist2)[movedPiece][to_sq(move)]
                              + (*pos->mainHistory)[!stm()][from_to(move)] - lmr_v3;

                // Decrease/increase reduction by comparing with opponent's stat score.
                if (ss->statScore >= -lmr_v4 && (ss - 1)->statScore < -lmr_v5)
                    r -= r_v11;

                else if ((ss - 1)->statScore >= -lmr_v6 && ss->statScore < -lmr_v7)
                    r += r_v12;

                // Decrease/increase reduction for moves with a good/bad history.
                r -= ss->statScore / lmr_v8 * r_v13;
            }

            Depth d = clamp(newDepth - r / 1000, 1, newDepth);
            value   = -search(pos, ss + 1, -(alpha + 1), -alpha, d, true, false);

            if (value > alpha && d < newDepth)
            {
                value = -search(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode, false);
                if (!captureOrPromotion)
                {
                    int bonus = value > alpha ? stat_bonus(newDepth) : -stat_malus(newDepth);

                    if (move == ss->killers[0])
                        bonus += bonus * kb_v1 / 100;

                    update_continuation_histories(ss, movedPiece, to_sq(move), bonus);
                }
            }
        }
        // Step 17. Full depth search when LMR is skipped or fails high.
        else if (!PvNode || moveCount > 1)
        {
            value = -search(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode, false);
        }

        // For PV nodes only, do a full PV search on the first move or after a fail
        // high (in the latter case search only if value < beta), otherwise let the
        // parent node fail low with value <= alpha and try another move.
        if (PvNode && (moveCount == 1 || value > alpha))
        {
            if (ss->ply <= 2)
            {
                (ss + 1)->pv    = pv;
                (ss + 1)->pv[0] = 0;
            }

            // Extend move from transposition table if we are about to dive into qsearch.
            if (move == ttMove && ss->ply <= pos->rootDepth * 2)
                newDepth = max(newDepth, 1);

            value = -search(pos, ss + 1, -beta, -alpha, newDepth, false, true);
        }

        // Step 18. Undo move
        undo_move(pos, move);

        // Step 19. Check for a new best move
        // Finished searching the move. If a stop occurred, the return value of
        // the search cannot be trusted, and we return immediately without
        // updating best move, PV and TT.
        if (Thread.stop)
        {
            return 0;
        }

        if (rootNode)
        {
            RootMove* rm = NULL;
            for (int idx = 0; idx < pos->rootMoves->size; idx++)
                if (pos->rootMoves->move[idx].pv[0] == move)
                {
                    rm = &pos->rootMoves->move[idx];
                    break;
                }

            // PV move or new best move ?
            if (moveCount == 1 || value > alpha)
            {
                rm->score  = value;
                rm->pvSize = 1;

                for (Move* m = (ss + 1)->pv; *m; ++m)
                    rm->pv[rm->pvSize++] = *m;

                // We record how often the best move has been changed in each
                // iteration. This information is used for time management: When
                // the best move changes frequently, we allocate some more time.
                if (moveCount > 1)
                    pos->bestMoveChanges++;
            }
            else
                // All other moves but the PV are set to the lowest value: this is
                // not a problem when sorting because the sort is stable and the
                // move position in the list is preserved - just the PV is pushed up.
                rm->score = -VALUE_INFINITE;
        }

        if (value > bestValue)
        {
            bestValue = value;

            if (value > alpha)
            {
                bestMove = move;

                if (PvNode && !rootNode && ss->ply <= 2)  // Update pv even in fail-high case
                    update_pv(ss->pv, move, (ss + 1)->pv);

                if (value >= beta)
                {
                    ss->cutoffCnt += !ttMove + (extension < 2);
                    ss->statScore = 0;
                    break;
                }

                alpha = value;
            }
        }

        if (move != bestMove && moveCount < 32)
        {
            if (captureOrPromotion)
                capturesSearched[captureCount++] = move;
            else
                quietsSearched[quietCount++] = move;
        }
    }


    // The following condition would detect a stop only after move loop has
    // been completed. But in this case bestValue is valid because we have
    // fully searched our subtree, and we can anyhow save the result in TT.
    /*
  if (Threads.stop)
    return VALUE_DRAW;
  */

    if (!PvNode && bestValue >= beta && abs(bestValue) < VALUE_MATE_IN_MAX_PLY
        && abs(beta) < VALUE_MATE_IN_MAX_PLY && abs(alpha) < VALUE_MATE_IN_MAX_PLY)
        bestValue = (bestValue * depth + beta) / (depth + 1);

    // Step 20. Check for mate and stalemate
    // All legal moves have been searched and if there are no legal moves,
    // it must be a mate or a stalemate. If we are in a singular extension
    // search then return a fail low score.
    if (!moveCount)
        bestValue = excludedMove ? alpha : inCheck ? mated_in(ss->ply) : VALUE_DRAW;
    else if (bestMove)
    {
        // Quiet best move: update move sorting heuristics
        if (!is_capture_or_promotion(pos, bestMove))
        {
            int bonus = stat_bonus(depth + (bestValue > beta + qb_v1));
            update_quiet_stats(pos, ss, bestMove, bonus);

            // Decrease all the other played quiet moves
            for (int i = 0; i < quietCount; i++)
            {
                history_update(*pos->mainHistory, stm(), quietsSearched[i], -bonus);
                update_continuation_histories(ss, moved_piece(quietsSearched[i]),
                                              to_sq(quietsSearched[i]), -bonus);
            }
        }

        update_capture_stats(pos, bestMove, capturesSearched, captureCount, stat_bonus(depth + 1));

        // Extra penalty for a quiet TT or main killer move in previous ply
        // when it gets refuted
        if ((prevSq != SQ_NONE && (ss - 1)->moveCount == 1
             || (ss - 1)->currentMove == (ss - 1)->killers[0])
            && !captured_piece())
            update_continuation_histories(ss - 1, piece_on(prevSq), prevSq, -stat_malus(depth + 1));
    }
    // Bonus for prior countermove that caused the fail low
    else if (!captured_piece() && prevSq != SQ_NONE)
    {
        int bonus = pcmb_v1 * (depth > pcmb_v2 / 100) + pcmb_v3 * !(PvNode || cutNode)
                  + pcmb_v4 * ((ss - 1)->moveCount > pcmb_v5 / 100)
                  + pcmb_v6 * (!inCheck && bestValue <= ss->staticEval - pcmb_v7);

        // Proportional to "how much damage we have to undo"
        bonus += min(-(ss - 1)->statScore / pcmb_v8, pcmb_v9);

        bonus = max(bonus, 0);
        update_continuation_histories(ss - 1, piece_on(prevSq), prevSq,
                                      stat_bonus(depth) * bonus / pcmb_v10);

        history_update(*pos->mainHistory, !stm(), (ss - 1)->currentMove,
                       stat_bonus(depth) * bonus / pcmb_v11);
    }

    // If no good move is found and the previous position was ttPv, then the
    // previous opponent move is probably good and the new position is added
    // to the search tree
    if (bestValue <= alpha)
        ss->ttPv = ss->ttPv || ((ss - 1)->ttPv && depth > 3);

    if (!excludedMove)
        tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                 bestValue >= beta    ? BOUND_LOWER
                 : PvNode && bestMove ? BOUND_EXACT
                                      : BOUND_UPPER,
                 depth, bestMove, rawEval);

    // Adjust correction history
    if (!inCheck && (!bestMove || !is_capture_or_promotion(pos, bestMove))
        && !(bestValue >= beta && bestValue <= ss->staticEval)
        && !(!bestMove && bestValue >= ss->staticEval))
    {
        add_correction_history(*pos->matCorrHist, stm(), material_key(), depth,
                               bestValue - ss->staticEval);
        add_correction_history(*pos->pawnCorrHist, stm(), pawn_key(), depth,
                               bestValue - ss->staticEval);
    }

    return bestValue;
}

// qsearch() is the quiescence search function template, which is
// called by the main search function with zero depth, or recursively with
// further decreasing depth per call.
Value qsearch(Position*  pos,
              Stack*     ss,
              Value      alpha,
              Value      beta,
              Depth      depth,
              const int  NT,
              const bool InCheck) {
    const bool PvNode = NT == PV;

    Move     pv[3];
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, bestMove;
    Value    bestValue, value, rawEval, ttValue, futilityValue, futilityBase;
    bool     ttHit, pvHit, givesCheck;
    Depth    ttDepth;
    int      moveCount;

    if (PvNode && ss->ply <= 2)
    {
        (ss + 1)->pv = pv;
        ss->pv[0]    = 0;
    }

    bestMove  = 0;
    moveCount = 0;

    // Check for an instant draw or if the maximum ply has been reached
    if (is_draw(pos) || ss->ply >= MAX_PLY)
        return ss->ply >= MAX_PLY && !InCheck ? evaluate(pos) : VALUE_DRAW;

    // Decide whether or not to include checks: this fixes also the type of
    // TT entry depth that we are going to use. Note that in qsearch we use
    // only two types of depth in TT: DEPTH_QS_CHECKS or DEPTH_QS_NO_CHECKS.
    ttDepth = InCheck || depth >= DEPTH_QS_CHECKS ? DEPTH_QS_CHECKS : DEPTH_QS_NO_CHECKS;

    // Transposition table lookup
    posKey  = key();
    tte     = tt_probe(posKey, &ttHit);
    ttValue = ttHit ? value_from_tt(tte_value(tte), ss->ply, rule50_count()) : VALUE_NONE;
    ttMove  = ttHit ? tte_move(tte) : 0;
    pvHit   = ttHit && tte_is_pv(tte);

    if (!PvNode && ttHit && tte_depth(tte) >= ttDepth
        && ttValue != VALUE_NONE  // Only in case of TT access race
        && (ttValue >= beta ? (tte_bound(tte) & BOUND_LOWER) : (tte_bound(tte) & BOUND_UPPER)))
        return ttValue;

    // Evaluate the position statically
    if (InCheck)
    {
        rawEval        = VALUE_NONE;
        ss->staticEval = VALUE_NONE;
        bestValue = futilityBase = -VALUE_INFINITE;
    }
    else
    {
        if (ttHit)
        {
            // Never assume anything about values stored in TT
            if ((rawEval = tte_eval(tte)) == VALUE_NONE)
                rawEval = evaluate(pos);

            ss->staticEval = bestValue = to_corrected(pos, rawEval);


            // Can ttValue be used as a better position evaluation?
            if (ttValue != VALUE_NONE
                && (tte_bound(tte) & (ttValue > bestValue ? BOUND_LOWER : BOUND_UPPER)))
                bestValue = ttValue;
        }
        else
        {
            rawEval =
              (ss - 1)->currentMove != MOVE_NULL ? evaluate(pos) : -(ss - 1)->staticEval + tempo;

            ss->staticEval = bestValue = to_corrected(pos, rawEval);
        }

        // Stand pat. Return immediately if static value is at least beta
        if (bestValue >= beta)
        {
            if (!ttHit)
                tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), false, BOUND_LOWER,
                         DEPTH_NONE, 0, rawEval);

            return bestValue;
        }

        if (PvNode && bestValue > alpha)
            alpha = bestValue;

        futilityBase = bestValue + qsf_v1;
    }

    ss->continuationHistory = &(*pos->contHist)[0][0];

    // Initialize move picker data for the current position, and prepare
    // to search the moves. Because the depth is <= 0 here, only captures,
    // queen promotions and checks (only if depth >= DEPTH_QS_CHECKS) will
    // be generated.

    Square prevSq = move_is_ok((ss - 1)->currentMove) ? to_sq((ss - 1)->currentMove) : SQ_NONE;

    mp_init_q(pos, ttMove, depth, prevSq);

    // Loop through the moves until no moves remain or a beta cutoff occurs
    while ((move = next_move(pos, 0)))
    {
        // Check for legality just before making the move
        if (!is_legal(pos, move))
            continue;
        moveCount++;

        givesCheck = gives_check(pos, ss, move);

        // Futility pruning
        if (bestValue > VALUE_MATED_IN_MAX_PLY && !givesCheck && futilityBase > -VALUE_KNOWN_WIN
            && type_of_m(move) != PROMOTION)
        {
            if (moveCount > 2)
                continue;

            futilityValue = futilityBase + PieceValue[piece_on(to_sq(move))];

            if (futilityValue <= alpha)
            {
                bestValue = max(bestValue, futilityValue);
                continue;
            }

            if (futilityBase <= alpha && !see_test(pos, move, 1))
            {
                bestValue = max(bestValue, futilityBase);
                continue;
            }
        }

        // Do not search moves with negative SEE values
        if (bestValue > VALUE_MATED_IN_MAX_PLY && !see_test(pos, move, 0))
            continue;

        // Speculative prefetch as early as possible
        prefetch(tt_first_entry(key_after(pos, move)));

        ss->currentMove         = move;
        bool captureOrPromotion = is_capture_or_promotion(pos, move);
        ss->continuationHistory = &(*pos->contHist)[moved_piece(move)][to_sq(move)];

        if (!captureOrPromotion && bestValue > VALUE_MATED_IN_MAX_PLY
            && (*(ss - 1)->continuationHistory)[moved_piece(move)][to_sq(move)]
                 < CounterMovePruneThreshold
            && (*(ss - 2)->continuationHistory)[moved_piece(move)][to_sq(move)]
                 < CounterMovePruneThreshold)
            continue;

        // Make and search the move
        do_move(pos, move, givesCheck);

        value = -qsearch(pos, ss + 1, -beta, -alpha, depth - 1, PvNode, givesCheck);
        undo_move(pos, move);

        // Check for a new best move
        if (value > bestValue)
        {
            bestValue = value;

            if (value > alpha)
            {
                bestMove = move;

                if (value >= beta)
                    break;

                alpha = value;
            }
        }
    }

    // All legal moves have been searched. A special case: If we're in check
    // and no legal moves were found, it is checkmate.
    if (InCheck && bestValue == -VALUE_INFINITE)
        return mated_in(ss->ply);  // Plies to mate from the root

    if (abs(bestValue) < VALUE_MATE_IN_MAX_PLY && bestValue >= beta)
        bestValue = (3 * bestValue + beta) / 4;

    tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), pvHit,
             bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, ttDepth, bestMove, rawEval);

    return bestValue;
}

#define rm_lt(m1, m2) \
    ((m1).score != (m2).score ? (m1).score < (m2).score : (m1).previousScore < (m2).previousScore)

// stable_sort() sorts RootMoves from highest-scoring move to lowest-scoring
// move while preserving order of equal elements.
static void stable_sort(RootMove* rm, int num) {
    int i, j;

    for (i = 1; i < num; i++)
        if (rm_lt(rm[i - 1], rm[i]))
        {
            RootMove tmp = rm[i];
            rm[i]        = rm[i - 1];
            for (j = i - 1; j > 0 && rm_lt(rm[j - 1], tmp); j--)
                rm[j] = rm[j - 1];
            rm[j] = tmp;
        }
}

// value_to_tt() adjusts a mate score from "plies to mate from the root" to
// "plies to mate from the current position". Non-mate scores are unchanged.
// The function is called before storing a value in the transposition table.

static Value value_to_tt(Value v, int ply) {
    return v >= VALUE_MATE_IN_MAX_PLY ? v + ply : v <= VALUE_MATED_IN_MAX_PLY ? v - ply : v;
}


// value_from_tt() is the inverse of value_to_tt(): It adjusts a mate score
// from the transposition table (which refers to the plies to mate/be mated
// from current position) to "plies to mate/be mated from the root".

static Value value_from_tt(Value v, int ply, int r50c) {
    if (v == VALUE_NONE)
        return VALUE_NONE;

    if (v >= VALUE_MATE_IN_MAX_PLY)
        return (VALUE_MATE - v > 99 - r50c) ? VALUE_MATE_IN_MAX_PLY - 1 : v - ply;

    if (v <= VALUE_MATED_IN_MAX_PLY)
        return (VALUE_MATE + v > 99 - r50c) ? VALUE_MATED_IN_MAX_PLY + 1 : v + ply;

    return v;
}

// update_pv() adds current move and appends child pv[]

static void update_pv(Move* pv, Move move, Move* childPv) {
    for (*pv++ = move; childPv && *childPv;)
        *pv++ = *childPv++;
    *pv = 0;
}

// differential.
static void
add_correction_history(CorrectionHistory hist, Color side, Key key, Depth depth, int32_t diff) {
    int16_t* entry      = &hist[side][key % CORRECTION_HISTORY_ENTRY_NB];
    int32_t  newWeight  = min(ch_v1 / 100, 1 + depth);
    int32_t  scaledDiff = diff * ch_v2;
    int32_t  update     = *entry * (ch_v3 - newWeight) + scaledDiff * newWeight;
    // Clamp entry in-bounds.
    *entry = max(-CORRECTION_HISTORY_MAX, min(CORRECTION_HISTORY_MAX, update / ch_v3));
}

Value to_corrected(Position* pos, Value rawEval) {
    int32_t mch = ch_v4 * (*pos->matCorrHist)[stm()][material_key() % CORRECTION_HISTORY_ENTRY_NB];
    int32_t pch = ch_v5 * (*pos->pawnCorrHist)[stm()][pawn_key() % CORRECTION_HISTORY_ENTRY_NB];
    Value   v   = rawEval + (pch + mch) / 100 / ch_v2;
    v           = clamp(v, -VALUE_MATE_IN_MAX_PLY, VALUE_MATE_IN_MAX_PLY);
    return v;
}

// update_continuation_histories() updates countermove and follow-up move history.

static void update_continuation_histories(Stack* ss, Piece pc, Square s, int bonus) {
    if (move_is_ok((ss - 1)->currentMove))
        update_contHist(*(ss - 1)->continuationHistory, pc, s, bonus);

    if (move_is_ok((ss - 2)->currentMove))
        update_contHist(*(ss - 2)->continuationHistory, pc, s, bonus);

    if (ss->checkersBB)
        return;

    if (move_is_ok((ss - 4)->currentMove))
        update_contHist(*(ss - 4)->continuationHistory, pc, s, bonus);

    if (move_is_ok((ss - 6)->currentMove))
        update_contHist(*(ss - 6)->continuationHistory, pc, s, bonus);
}

// update_capture_stats() updates move sorting heuristics when a new capture
// best move is found

static void
update_capture_stats(const Position* pos, Move move, Move* captures, int captureCnt, int bonus) {
    Piece moved_piece = moved_piece(move);
    int   captured    = type_of_p(piece_on(to_sq(move)));

    if (is_capture_or_promotion(pos, move))
        cpth_update(*pos->captureHistory, moved_piece, to_sq(move), captured, bonus);

    // Decrease all the other played capture moves
    for (int i = 0; i < captureCnt; i++)
    {
        moved_piece = moved_piece(captures[i]);
        captured    = type_of_p(piece_on(to_sq(captures[i])));
        cpth_update(*pos->captureHistory, moved_piece, to_sq(captures[i]), captured, -bonus);
    }
}

// update_quiet_stats() updates killers, history, countermove and countermove
// plus follow-up move history when a new quiet best move is found.

static void update_quiet_stats(const Position* pos, Stack* ss, Move move, int bonus) {
    if (ss->killers[0] != move)
    {
        ss->killers[1] = ss->killers[0];
        ss->killers[0] = move;
    }

    Color c = stm();
    history_update(*pos->mainHistory, c, move, bonus);
    update_continuation_histories(ss, moved_piece(move), to_sq(move), bonus);

    if (type_of_p(moved_piece(move)) != PAWN)
        history_update(*pos->mainHistory, c, reverse_move(move), -bonus);
}

static int peak_stdin() {
#ifndef WIN32
    fd_set         rf = {0};
    struct timeval tv = {0, 0};
    FD_SET(fileno(stdin), &rf);
    select(fileno(stdin) + 1, &rf, NULL, NULL, &tv);
    return FD_ISSET(fileno(stdin), &rf);
#else
    static HANDLE hIn;
    static int    init = 0, pipe = 0;
    DWORD         dw;

    if (!init++)
    {
        hIn  = GetStdHandle(STD_INPUT_HANDLE);
        pipe = !GetConsoleMode(hIn, &dw);
        if (!pipe)
        {
            SetConsoleMode(hIn, dw & ~(ENABLE_MOUSE_INPUT | ENABLE_WINDOW_INPUT));
            FlushConsoleInputBuffer(hIn);
        }
    }

    if (pipe)
        return PeekNamedPipe(hIn, NULL, 0, NULL, &dw, NULL) ? dw : 1;

    GetNumberOfConsoleInputEvents(hIn, &dw);
    return dw > 1 ? dw : 0;
#endif
}

// check_time() is used to print debug info and, more importantly, to detect
// when we are out of available time and thus stop the search.
static void check_time(void) {
    if (Thread.ponder)
    {
        if (peak_stdin())
            Thread.stop = 1;
        return;
    }

    TimePoint elapsed = time_elapsed();
#ifndef KAGGLE
    if ((use_time_management() && elapsed > time_maximum() - 10)
        || (Limits.nodes && Thread.pos->nodes >= Limits.nodes))
        Thread.stop = 1;
#else
    if (elapsed > time_maximum() - 10)
        Thread.stop = 1;
#endif
}

// uci_print_pv() prints PV information according to the UCI protocol.
// UCI requires that all (if any) unsearched PV lines are sent with a
// previous search score.

static void uci_print_pv(Position* pos, Depth depth, Value alpha, Value beta) {
    TimePoint  elapsed        = time_elapsed() + 1;
    RootMoves* rm             = pos->rootMoves;
    uint64_t   nodes_searched = Thread.pos->nodes;
    char       buf[16];

    int i = 0;

    bool updated = rm->move[i].score != -VALUE_INFINITE;

    Depth d = updated ? depth : max(1, depth - 1);
    Value v = updated ? rm->move[i].score : rm->move[i].previousScore;

    if (v == -VALUE_INFINITE)
        v = VALUE_ZERO;

    printf("info depth %d score %s nodes %" PRIu64 " nps %" PRIu64 " time %" PRIi64 " pv", d,
           uci_value(buf, v), nodes_searched, nodes_searched * 1000 / elapsed, elapsed);

    for (int idx = 0; idx < rm->move[i].pvSize; idx++)
        printf(" %s", uci_move(buf, rm->move[i].pv[idx]));
    printf("\n");

    fflush(stdout);
}

static int extract_ponder_from_tt(RootMove* rm, Position* pos) {
    if (!rm->pv[0])
        return 0;

    do_move(pos, rm->pv[0], gives_check(pos, pos->st, rm->pv[0]));

    bool     ttHit;
    TTEntry* tte = tt_probe(key(), &ttHit);
    if (ttHit && is_pseudo_legal(pos, tte_move(tte)))
    {
        rm->pv[1]  = tte_move(tte);
        rm->pvSize = 2;
    }

    undo_move(pos, rm->pv[0]);
    return rm->pvSize > 1;
}

// start_thinking() wakes up the main thread to start a new search,
// then returns immediately.

void start_thinking(Position* root) {
    prepare_for_search(root);
    mainthread_search();
}

SMALL void prepare_for_search(Position* root) {
    Thread.stop          = false;
    Thread.increaseDepth = true;

    // Generate all legal moves.
    ExtMove  list[MAX_MOVES];
    ExtMove* end = generate_pseudo_legal(root, list);

    Position* pos  = Thread.pos;
    pos->rootDepth = 0;
    pos->nodes     = 0;

    RootMoves* rm = pos->rootMoves;

    rm->size = end - list;
    for (int i = 0; i < rm->size; i++)
    {
        rm->move[i].pvSize        = 1;
        rm->move[i].pv[0]         = list[i].move;
        rm->move[i].score         = -VALUE_INFINITE;
        rm->move[i].previousScore = -VALUE_INFINITE;
    }
    memcpy(pos, root, offsetof(Position, moveList));

    // Copy enough of the root State buffer.
    int n = max(7, root->st->pliesFromNull);
    for (int i = 0; i <= n; i++)
        memcpy(&pos->stack[i], &root->st[i - n], StateSize);
    pos->st                 = pos->stack + n;
    (pos->st - 1)->endMoves = pos->moveList;

    pos_set_check_info(pos);
}
