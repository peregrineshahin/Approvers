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
#include "movepick.h"
#include "search.h"
#include "timeman.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"

#ifndef KAGGLE
Parameter parameters[255];
int       parameters_count = 0;

    #define STRINGIFY(x) #x
    #define PARAM(Name, Value, Step) \
        int                                               Name = Value; \
        __attribute__((constructor(__LINE__ + 100))) void ctor_##Name(void) { \
            Parameter* param = &parameters[parameters_count++]; \
            strncpy(param->name, STRINGIFY(Name), sizeof(param->name)); \
            param->value = &Name; \
            param->step  = Step; \
        }
#else
    #define PARAM(Name, Value, Step) int Name = Value;
#endif

// Disabled parameters
PARAM(nmp_v1, 764, 0)
PARAM(nmp_v2, 56, 0)
PARAM(nmp_v3, 165, 0)
PARAM(nmp_v4, 182, 0)
PARAM(mp_v1, 70, 0)
PARAM(mp_v2, 1064, 0)
PARAM(r_v1, 1056, 0)

// Search parameters
PARAM(ttct_v1, 138, 12.0)
PARAM(ttct_v2, 131, 12.0)
PARAM(qmo_v1, 430, 38.0)
PARAM(qmo_v2, 1274, 120.0)
PARAM(qmo_v3, 918, 120.0)
PARAM(rz_v1, 6, 0.5)
PARAM(rz_v2, 354, 16.0)
PARAM(rz_v3, 242, 10.0)
PARAM(ft_v1, 87, 12.0)
PARAM(ft_v2, 16, 0.5)
PARAM(nmp_v5, 25301, 150.0)
PARAM(nmp_v6, 27, 3.0)
PARAM(nmp_v8, 91, 6.0)
PARAM(nmp_v9, 204, 12.0)
PARAM(rd_v1, 578, 60.0)
PARAM(rd_v2, 1242, 48.0)
PARAM(rd_v3, 908, 72.0)
PARAM(rd_init_v1, 3667, 150.0)
PARAM(d_v1, 18, 1.2)
PARAM(cbp_v2, 4, 9.6)
PARAM(cbp_v3, -11, 9.6)
PARAM(cbp_v4, -8, 9.6)
PARAM(fpp_v1, 7, 0.5)
PARAM(fpp_v2, 219, 24.0)
PARAM(fpp_v3, 172, 21.6)
PARAM(fpp_v4, 67, 7.2)
PARAM(fpp_v5, 70, 7.2)
PARAM(sqsee_v1, 27, 2.4)
PARAM(scsee_v1, 207, 12.0)
PARAM(se_v2, 154, 15.0)
PARAM(se_v5, 27, 3.6)
PARAM(prb_v1, 120, 14.4)
PARAM(prb_v2, 47, 4.8)
PARAM(lmr_v3, 3689, 300.0)
PARAM(lmr_v8, 14516, 150.0)
PARAM(hb_v1, 925, 75.0)
PARAM(hb_v2, 199, 18.0)
PARAM(hb_v3, 145, 24.0)
PARAM(hb_v4, 2510, 180.0)
PARAM(hm_v1, 716, 75.0)
PARAM(hm_v2, 201, 18.0)
PARAM(hm_v3, 125, 24.0)
PARAM(hm_v4, 1780, 180.0)
PARAM(cnht_v1, 1007, 80.0)
PARAM(cnht_v2, 1038, 80.0)
PARAM(cnht_v3, 987, 80.0)
PARAM(cnht_v4, 975, 80.0)
PARAM(asd_v1, 4, 0.5)
PARAM(qsf_v1, 193, 18.0)
PARAM(ch_v1, 22, 3.6)
PARAM(ch_v2, 175, 18.0)
PARAM(ch_v3, 248, 24.0)
PARAM(ch_v4, 96, 12.0)
PARAM(ch_v5, 137, 12.0)
PARAM(ch_v6, 131, 12.0)
PARAM(ch_v7, 119, 12.0)
PARAM(ch_v8, 121, 12.0)
PARAM(tempo, 45, 4.8)
PARAM(mp_v3, 2623, 180.0)
PARAM(mp_v4, 188, 12.0)
PARAM(mp_v5, 328, 12.0)
PARAM(mp_v6, 293, 12.0)
PARAM(mp_v7, 248, 12.0)
PARAM(mp_v8, 134, 12.0)
PARAM(eval_scale, 92, 3.0)
PARAM(mat_scale, 25957, 420.0)
PARAM(mat_n, 670, 60.0)
PARAM(mat_b, 715, 78.0)
PARAM(mat_r, 1477, 120.0)
PARAM(mat_q, 2567, 240.0)
PARAM(pcmb_v1, 79, 12.0)
PARAM(pcmb_v4, 154, 12.0)
PARAM(pcmb_v5, 7, 0.5)
PARAM(pcmb_v6, 110, 12.0)
PARAM(pcmb_v7, 95, 12.0)
PARAM(pcmb_v8, 130, 8.4)
PARAM(pcmb_v9, 252, 30.0)
PARAM(pcmb_v10, 125, 8.4)
PARAM(pcmb_v11, 143, 8.4)
PARAM(r_v2, 2454, 250.0)
PARAM(r_v6, 1244, 150.0)
PARAM(r_v7, 1087, 150.0)
PARAM(r_v8, 1325, 150.0)
PARAM(r_v9, 920, 150.0)
PARAM(r_v13, 998, 50.0)
PARAM(ded_v1, 63, 7.2)
PARAM(qb_v1, 182, 18.0)
PARAM(qb_v2, 185, 18.0)
PARAM(cms_v1, 29166, 300.0)
PARAM(hu_v1, 10294, 300.0)
PARAM(cpth_v1, 11627, 300.0)

// Time management parameters
PARAM(tm_v1, 329, 28)
PARAM(tm_v2, 555, 63)
PARAM(tm_v3, 657, 63)
PARAM(tm_v4, 809, 40)
PARAM(tm_v5, 52, 7)
PARAM(tm_v6, 155, 7)
PARAM(tm_v7, 869, 75)
PARAM(tm_v8, 187, 17)
PARAM(tm_v9, 101, 11)
PARAM(tm_v10, 147, 7)
PARAM(tm_v11, 216, 16)
PARAM(tm_v13, 89, 4)
PARAM(tm_v14, 296, 14)
PARAM(tm_v15, 248, 14)
PARAM(tm_v16, 21, 2)
PARAM(tm_v17, 674, 29)
PARAM(tm_v18, 416, 14)
PARAM(tm_v19, 1186, 14)
PARAM(tm_v20, 801, 14)
PARAM(tm_v21, 100, 4)
PARAM(tm_v22, 170, 7)

LimitsType Limits;

enum {
    NonPV,
    PV
};

static int futility_margin(Depth d, bool improving) {
    return ft_v1 * (d - improving) + ft_v2 * d * d;
}

// Reductions lookup tables, initialized at startup
static int Reductions[MAX_MOVES];  // [depth or moveNumber]

static Depth reduction(int i, Depth d, int mn) {
    int r = Reductions[d] * Reductions[mn];
    return (r + rd_v1) / rd_v2 + (!i && r > rd_v3);
}

static int futility_move_count(bool improving, Depth depth) {
    return improving ? 3 + depth * depth : (2 + depth * depth) / 2;
}

// History and stats update bonus, based on depth
static Value stat_bonus(Depth d) { return min((hb_v1 * d / 128 + hb_v2) * d - hb_v3, hb_v4); }

// History and stats update malus, based on depth
static Value stat_malus(Depth d) { return min((hm_v1 * d / 128 + hm_v2) * d - hm_v3, hm_v4); }

static Value value_to_tt(Value v, int ply);
static Value value_from_tt(Value v, int ply, int r50c);
static void  update_continuation_histories(Stack* ss, Piece pc, Square s, int bonus);
Value        to_corrected(Position* pos, Value unadjustedStaticEval);
static void  update_correction_histories(const Position* pos, Depth depth, int32_t diff);
static void  update_quiet_stats(const Position* pos, Stack* ss, Move move, int bonus);
static void
update_capture_stats(const Position* pos, Move move, Move* captures, int captureCnt, Depth depth);
static void check_time(void);
static void uci_print_pv(Position* pos, Depth depth);

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

// Called during startup to initialize various lookup tables
SMALL void search_init(void) {
    for (int i = 1; i < MAX_MOVES; i++)
        Reductions[i] = (int) (rd_init_v1 / 128.0 * my_log(i));
}


// Resets search state to zero, to obtain reproducible results
SMALL void search_clear(void) {
    Time.availableNodes = 0;

    tt_clear();

    Position* pos = Thread.pos;
    stats_clear(pos->mainHistory);
    stats_clear(pos->captureHistory);
    stats_clear(pos->contHist);

#pragma clang loop unroll(disable)
    for (int pc = 0; pc < 7; pc++)
#pragma clang loop unroll(disable)
        for (int sq = 0; sq < 64; sq++)
            (*pos->contHist)[0][0][0][pc][sq] = -1;

    Thread.previousScore = VALUE_INFINITE;
}

// Called by the main thread when the program receives the UCI 'go' command.
// It searches from the root position and outputs the "bestmove".
void mainthread_search(void) {
    Position* pos = Thread.pos;
    Color     us  = stm();
    time_init(us, game_ply());
    tt_new_search();
    char buf[16];

    Thread.pos->bestMoveChanges = 0;
    thread_search(pos);
    Thread.previousScore = pos->st->pv.score;

    printf("bestmove %s\n", uci_move(buf, pos->st->pv.line[0]));
    fflush(stdout);

    if (!IsKaggle && !Thread.testPonder)
        return;

    // Start pondering right after the best move has been printed if we can
    if (pos->st->pv.length >= 1)
    {
        Thread.ponder = true;
        Thread.stop   = false;

        const Move bestMove = pos->st->pv.line[0];
        do_move(pos, bestMove, gives_check(pos, pos->st, bestMove));

        pos->completedDepth = 0;
        pos->rootDepth      = 0;

        prepare_for_search(pos);
        thread_search(pos);

        Thread.ponder = false;
        Thread.stop   = true;
    }
}


// Main iterative deepening loop. It calls search() repeatedly with increasing
// depth until the allocated thinking time has been consumed, the user stops
// the search, or the maximum search depth is reached.
void thread_search(Position* pos) {
    Value  bestValue, alpha, beta, delta;
    Move   lastMove           = MOVE_NONE;
    Depth  pvStability        = 0;
    double totBestMoveChanges = 0;
    int    iterIdx            = 0;

    Stack* ss = pos->st;  // At least the seventh element of the allocated array.
#pragma clang loop unroll(disable)
    for (int i = -7; i < 3; i++)
    {
        memset(SStackBegin(ss[i]), 0, SStackSize);
    }
    (ss - 1)->endMoves = pos->moveList;

#pragma clang loop unroll(disable)
    for (int i = -7; i < 0; i++)
    {
        ss[i].continuationHistory = &(*pos->contHist)[0][0][0];  // Use as sentinel
        ss[i].staticEval          = VALUE_NONE;
        ss[i].checkersBB          = 0;
    }

#pragma clang loop unroll(disable)
    for (int i = 0; i <= MAX_PLY; i++)
        ss[i].ply = i;

    ss->accumulator.needs_refresh = 1;

    bestValue = delta = alpha = -VALUE_INFINITE;
    beta                      = VALUE_INFINITE;
    pos->completedDepth       = 0;

    int value = Thread.previousScore == VALUE_INFINITE ? VALUE_ZERO : Thread.previousScore;
#pragma clang loop unroll(disable)
    for (int i = 0; i < 4; i++)
        Thread.iterValue[i] = value;

    PVariation* pv = &pos->st->pv;

    // Iterative deepening loop until requested to stop or the target depth
    // is reached.
    while ((pos->rootDepth += 2) < MAX_PLY && !Thread.stop
           && !(Limits.depth && pos->rootDepth > Limits.depth))
    {
        // Age out PV variability metric
        totBestMoveChanges /= 2;

        // Reset aspiration window starting size
        if (pos->rootDepth >= 4)
        {
            Value previousScore = pv->score;
            delta               = d_v1;
            alpha               = max(previousScore - delta, -VALUE_INFINITE);
            beta                = min(previousScore + delta, VALUE_INFINITE);
        }

        while (true)
        {
            bestValue = search(pos, ss, alpha, beta, pos->rootDepth, false, true);

            // If search has been stopped, we break immediately
            if (Thread.stop)
                break;

            if (bestValue > alpha)
                break;

            beta  = (alpha + beta) / 2;
            alpha = max(bestValue - delta, -VALUE_INFINITE);
            delta += delta / 4 + asd_v1;
        }

#ifndef KAGGLE
        if (!Thread.ponder)
            uci_print_pv(pos, pos->rootDepth);
#endif

        if (!Thread.stop)
            pos->completedDepth = pos->rootDepth;

        pvStability = pv->line[0] == lastMove ? min(pvStability + 1, 8) : 0;
        lastMove    = pv->line[0];

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

            double pvFactor = 1.2 - 0.05 * pvStability;

            // Use part of the gained time from a previous stable move for this move
            totBestMoveChanges += Thread.pos->bestMoveChanges;
            Thread.pos->bestMoveChanges = 0;

            double bestMoveInstability = (tm_v21 / 100.0) + (tm_v22 / 100.0) * totBestMoveChanges;

            double totalTime = time_optimum() * fallingEval * pvFactor * bestMoveInstability;

            // Stop the search if we have exceeded the totalTime (at least 1ms)
            if (time_elapsed() > totalTime)
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "stop".
                if (!Thread.ponder)
                    Thread.stop = true;
            }
        }

        Thread.iterValue[iterIdx] = bestValue;
        iterIdx                   = (iterIdx + 1) & 3;
    }
}

// Main search function for both PV and non-PV nodes
Value search(
  Position* pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode, const int NT) {
    const bool PvNode   = NT == PV;
    const bool rootNode = PvNode && ss->ply == 0;

    ss->pv.length = 0;

    // Check if we have an upcoming move which draws by repetition, or if the
    // opponent had an alternative move earlier to this position.
    if (pos->st->pliesFromNull >= 3 && alpha < VALUE_DRAW && !rootNode
        && has_game_cycle(pos, ss->ply))
    {
        alpha = VALUE_DRAW;
        if (alpha >= beta)
            return alpha;
    }

    // Dive into quiescence search when the depth reaches zero
    if (depth <= 0)
        return qsearch(pos, ss, alpha, beta, 0);

    Move     capturesSearched[32], quietsSearched[32];
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, excludedMove, bestMove;
    Depth    extension, newDepth;
    Value    bestValue, value, ttValue, eval, unadjustedStaticEval, probCutBeta;
    bool     ttHit, givesCheck, improving;
    bool     capture, moveCountPruning;
    bool     ttCapture;
    int      moveCount, captureCount, quietCount;

    // Step 1. Initialize node
    moveCount = captureCount = quietCount = ss->moveCount = 0;
    bestValue                                             = -VALUE_INFINITE;

    // Check for the available remaining time
    if (pos->completedDepth >= 1 && (pos->nodes & 1023) == 0)
        check_time();

    // Step 2. Check for aborted search and immediate draw
    if (!rootNode && (Thread.stop || is_draw(pos) || ss->ply >= MAX_PLY))
        return ss->ply >= MAX_PLY && !ss->checkersBB ? evaluate(pos) : VALUE_DRAW;

    (ss + 1)->ttPv         = false;
    (ss + 1)->excludedMove = bestMove = 0;
    (ss + 2)->killers[0] = (ss + 2)->killers[1] = 0;
    (ss + 2)->cutoffCnt                         = 0;
    Square prevSq = move_is_ok((ss - 1)->currentMove) ? to_sq((ss - 1)->currentMove) : SQ_NONE;
    ss->statScore = 0;


    // Step 3. Transposition table lookup
    excludedMove = ss->excludedMove;
    posKey       = key();
    tte          = tt_probe(posKey, &ttHit);
    ttValue      = ttHit ? value_from_tt(tte_value(tte), ss->ply, rule50_count()) : VALUE_NONE;
    ttMove       = ttHit ? tte_move(tte) : 0;
    ttCapture    = ttMove && capture_stage(pos, ttMove);

    if (!excludedMove)
        ss->ttPv = PvNode || (ttHit && tte_is_pv(tte));

    // At non-PV nodes we check for an early TT cutoff
    if (!PvNode && ttHit && tte_depth(tte) >= depth && !excludedMove
        && tte_bound(tte) & (ttValue >= beta ? BOUND_LOWER : BOUND_UPPER))
    {
        // If ttMove is quiet, update move sorting heuristics on TT hit
        if (ttMove && ttValue >= beta)
        {
            if (!capture_stage(pos, ttMove))
                update_quiet_stats(pos, ss, ttMove, ttct_v1 * stat_bonus(depth) / 128);

            // Extra penalty for early quiet moves of the previous ply
            if ((ss - 1)->moveCount <= 2 && !captured_piece() && prevSq != SQ_NONE)
                update_continuation_histories(ss - 1, piece_on(prevSq), prevSq,
                                              ttct_v2 * -stat_malus(depth + 1) / 128);
        }

        // Partial workaround for the graph history interaction problem
        // For high rule50 counts don't produce transposition table cutoffs.
        if (rule50_count() < 90)
            return ttValue;
    }

    // Step 4. Static evaluation of the position
    if (ss->checkersBB)
    {
        // Skip early pruning when in check
        unadjustedStaticEval = eval = ss->staticEval = VALUE_NONE;
        improving                                    = false;
        goto moves_loop;
    }
    else if (excludedMove)
    {
        // Providing the hint that this node's accumulator will be used often
        evaluate(pos);
        unadjustedStaticEval = eval = ss->staticEval;
    }
    else if (ttHit)
    {
        // Never assume anything about values stored in TT
        if ((unadjustedStaticEval = tte_eval(tte)) == VALUE_NONE)
            unadjustedStaticEval = evaluate(pos);

        eval = ss->staticEval = to_corrected(pos, unadjustedStaticEval);

        // ttValue can be used as a better position evaluation
        if (ttValue != VALUE_NONE
            && (tte_bound(tte) & (ttValue > eval ? BOUND_LOWER : BOUND_UPPER)))
            eval = ttValue;
    }
    else
    {
        if ((ss - 1)->currentMove != MOVE_NULL)
            unadjustedStaticEval = evaluate(pos);
        else
            unadjustedStaticEval = -(ss - 1)->staticEval + tempo;

        eval = ss->staticEval = to_corrected(pos, unadjustedStaticEval);

        tte_save(tte, posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_NONE, 0,
                 unadjustedStaticEval);
    }

    improving = (ss - 2)->staticEval == VALUE_NONE
                ? (ss->staticEval > (ss - 4)->staticEval || (ss - 4)->staticEval == VALUE_NONE)
                : ss->staticEval > (ss - 2)->staticEval;

    if (prevSq != SQ_NONE && !(ss - 1)->checkersBB && !captured_piece())
    {
        int bonus = clamp(-depth * qmo_v1 / 128 * ((ss - 1)->staticEval + ss->staticEval - tempo),
                          -qmo_v2, qmo_v3);
        history_update(*pos->mainHistory, !stm(), (ss - 1)->currentMove, bonus);
    }

    // Step 5. Razoring
    if (!PvNode && !improving && depth < rz_v1 && eval < alpha - rz_v2 - rz_v3 * depth * depth)
    {
        value = qsearch(pos, ss, alpha - 1, alpha, 0);
        if (value < alpha)
            return value;
    }

    // Step 6. Futility pruning: child node
    if (!PvNode && eval - futility_margin(depth, improving) >= beta && eval < VALUE_MATE_IN_MAX_PLY
        && beta > -VALUE_MATE_IN_MAX_PLY)
        return eval;

    // Step 7. Null move search
    if (!PvNode && (ss - 1)->currentMove != MOVE_NULL && (ss - 1)->statScore < nmp_v5
        && eval >= beta && eval >= ss->staticEval
        && ss->staticEval >= beta - nmp_v6 * depth + nmp_v8 * ss->ttPv + nmp_v9 && !excludedMove
        && non_pawn_material(pos))
    {
        // Null move dynamic reduction based on depth and value
        Depth R = (nmp_v1 + nmp_v2 * depth) / nmp_v3 + min((eval - beta) / nmp_v4, 3) + ttCapture;

        ss->currentMove         = MOVE_NULL;
        ss->continuationHistory = &(*pos->contHist)[0][0][0];

        do_null_move(pos);
        ss->endMoves    = (ss - 1)->endMoves;
        Value nullValue = -search(pos, ss + 1, -beta, -beta + 1, depth - R, !cutNode, false);
        undo_null_move(pos);

        if (nullValue >= beta)
            return nullValue > VALUE_MATE_IN_MAX_PLY ? beta : nullValue;
    }

    probCutBeta = beta + prb_v1 - prb_v2 * improving;

    // Step 8. ProbCut
    // If we have a good enough capture and a reduced search returns a value
    // much above beta, we can (almost) safely prune the previous move
    if (!PvNode && depth > 4 && abs(beta) < VALUE_MATE_IN_MAX_PLY
        && !(ttHit && tte_depth(tte) >= depth - 3 && ttValue != VALUE_NONE
             && ttValue < probCutBeta))
    {

        if (ttHit && tte_depth(tte) >= depth - 3 && ttValue != VALUE_NONE && ttValue >= probCutBeta
            && ttMove && capture_stage(pos, ttMove))
            return probCutBeta;

        mp_init_pc(pos, ttMove, probCutBeta - ss->staticEval);
        int  probCutCount = 2 + 2 * cutNode;
        bool ttPv         = ss->ttPv;
        ss->ttPv          = false;

        while ((move = next_move(pos, 0)) && probCutCount)
            if (move != excludedMove && is_legal(pos, move))
            {
                capture = true;
                probCutCount--;

                ss->currentMove = move;
                ss->continuationHistory =
                  &(*pos->contHist)[stm()][type_of_p(moved_piece(move))][to_sq(move)];
                givesCheck = gives_check(pos, ss, move);
                do_move(pos, move, givesCheck);

                // Perform a preliminary qsearch to verify that the move holds
                value = -qsearch(pos, ss + 1, -probCutBeta, -probCutBeta + 1, 0);

                // If the qsearch held, perform the regular search
                if (value >= probCutBeta)
                    value = -search(pos, ss + 1, -probCutBeta, -probCutBeta + 1, depth - 4,
                                    !cutNode, false);
                undo_move(pos, move);
                if (value >= probCutBeta)
                {
                    tte_save(tte, posKey, value_to_tt(value, ss->ply), ttPv, BOUND_LOWER, depth - 3,
                             move, unadjustedStaticEval);
                    return value - (probCutBeta - beta);
                }
            }
        ss->ttPv = ttPv;
    }

    // Step 9. Internal iterative reductions
    if ((PvNode || cutNode) && depth >= (1 + cutNode) * 6 && !ttMove)
        depth -= 2;

moves_loop:  // When in check search starts from here.
  ;          // Avoid a compiler warning. A label must be followed by a statement.
    PieceToHistory* contHist0 = (ss - 1)->continuationHistory;
    PieceToHistory* contHist1 = (ss - 2)->continuationHistory;
    PieceToHistory* contHist2 = (ss - 4)->continuationHistory;

    mp_init(pos, ttMove, depth);

    value            = bestValue;
    moveCountPruning = false;

    // Step 10. Loop through all pseudo-legal moves until no moves remain or a beta cutoff occurs
    while ((move = next_move(pos, moveCountPruning)))
    {
        if (move == excludedMove)
            continue;

        // Check for legality
        if (!is_legal(pos, move))
            continue;

        ss->moveCount = ++moveCount;

        PieceType movedPiece = type_of_p(moved_piece(move));

        extension  = 0;
        capture    = capture_stage(pos, move);
        givesCheck = gives_check(pos, ss, move);

        // Calculate new depth for this move
        newDepth = depth - 1;

        Depth r = reduction(improving, depth, moveCount);

        // Step 11. Pruning at shallow depth
        if (!rootNode && non_pawn_material(pos) && bestValue > VALUE_MATED_IN_MAX_PLY)
        {
            // Skip quiet moves if movecount exceeds our FutilityMoveCount threshold
            moveCountPruning = moveCount >= futility_move_count(improving, depth);

            // Reduced depth of the next LMR search
            int lmrDepth = max(newDepth - r, 0);

            if (capture || givesCheck)
            {
                // SEE based pruning
                if (!see_test(pos, move, -scsee_v1 * depth))
                    continue;
            }
            else
            {
                // Countermoves based pruning
                if (lmrDepth < 3 + ((ss - 1)->statScore > cbp_v2 || (ss - 1)->moveCount == 1)
                    && (*contHist0)[movedPiece][to_sq(move)] < cbp_v3
                    && (*contHist1)[movedPiece][to_sq(move)] < cbp_v4)
                    continue;

                // Futility pruning: parent node
                if (lmrDepth < fpp_v1 && !ss->checkersBB
                    && ss->staticEval + (bestValue < ss->staticEval - fpp_v4 ? fpp_v2 : fpp_v5)
                           + fpp_v3 * lmrDepth
                         <= alpha)
                    continue;

                // Prune moves with negative SEE at low depths and below a decreasing
                // threshold at higher depths
                if (!see_test(pos, move, -(sqsee_v1 * lmrDepth * lmrDepth)))
                    continue;
            }
        }

        // Step 12. Extensions

        // Singular extension search. If all moves but one fail low on a search
        // of (alpha-s, beta-s), and just one fails high on (alpha, beta), then
        // that move is singular and should be extended. To verify this we do a
        // reduced search on all the other moves but the ttMove and if the
        // result is lower than ttValue minus a margin, then we extend the ttMove.
        if (depth >= 5 && move == ttMove && !rootNode
            && !excludedMove  // No recursive singular search
                              /* &&  ttValue != VALUE_NONE implicit in the next condition */
            && abs(ttValue) < VALUE_MATE_IN_MAX_PLY && (tte_bound(tte) & BOUND_LOWER)
            && tte_depth(tte) >= depth - 3)
        {
            Value singularBeta  = ttValue - se_v2 * depth / 128;
            Depth singularDepth = newDepth / 2;
            ss->excludedMove    = move;
            Move k1 = ss->mpKillers[0], k2 = ss->mpKillers[1];
            value = search(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode, false);
            ss->excludedMove = 0;

            if (value < singularBeta)
            {
                extension = 1;
                if (!PvNode && value < singularBeta - se_v5)
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
            else if (cutNode || ttValue >= beta)
                extension--;

            // The call to search_NonPV with the same value of ss messed up our
            // move picker data. So we fix it.
            mp_init(pos, ttMove, depth);
            ss->stage++;
            ss->mpKillers[0] = k1;
            ss->mpKillers[1] = k2;
        }

        // Last capture extension
        else if (PieceValue[captured_piece()] > PawnValue && low_material(pos))
            extension = 1;

        // Add extension to new depth
        newDepth += extension;

        // Step 13. Make the move.
        do_move(pos, move, givesCheck);

        // Speculative prefetch as early as possible
        prefetch(tt_first_entry(key()));

        // Update the current move (this must be done after singular extension search)
        ss->currentMove         = move;
        ss->continuationHistory = &(*pos->contHist)[stm()][movedPiece][to_sq(move)];

        r = r * r_v1;

        // Step 14. Late move reductions (LMR)
        if (depth >= 2 && moveCount > 1 + 2 * rootNode && (!capture || cutNode || !ss->ttPv))
        {
            // Decrease reduction if position is or has been on the PV
            if (ss->ttPv)
                r -= r_v2;

            if (cutNode)
                r += r_v8 + r_v9 * !capture;

            if (capture)
                ss->statScore = 0;
            else
            {
                // Increase reduction if ttMove is a capture
                if (ttCapture)
                    r += r_v6;

                if ((ss + 1)->cutoffCnt > 3)
                    r += r_v7;

                ss->statScore = (*contHist0)[movedPiece][to_sq(move)]
                              + (*contHist1)[movedPiece][to_sq(move)]
                              + (*contHist2)[movedPiece][to_sq(move)]
                              + (*pos->mainHistory)[!stm()][from_to(move)] - lmr_v3;
            }

            // Decrease/increase reduction for moves with a good/bad history.
            r -= ss->statScore / lmr_v8 * r_v13;

            Depth d = clamp(newDepth - r / 1000, 1, newDepth);
            value   = -search(pos, ss + 1, -(alpha + 1), -alpha, d, true, false);

            if (value > alpha && d < newDepth)
            {
                // Adjust full-depth search based on LMR results - if the result was
                // good enough search deeper, if it was bad enough search shallower.
                const bool doDeeperSearch    = value > bestValue + ded_v1;
                const bool doShallowerSearch = value < bestValue + newDepth;

                newDepth += doDeeperSearch - doShallowerSearch;

                if (newDepth > d)
                    value = -search(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode, false);

                if (!capture)
                {
                    int bonus = value > alpha ? stat_bonus(newDepth) : -stat_malus(newDepth);
                    update_continuation_histories(ss, make_piece(0, movedPiece), to_sq(move),
                                                  bonus);
                }
            }
        }
        // Step 15. Full-depth search when LMR is skipped or fails high.
        else if (!PvNode || moveCount > 1)
        {
            value = -search(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode, false);
        }

        // For PV nodes only, do a full PV search on the first move or after a fail
        // high (in the latter case search only if value < beta), otherwise let the
        // parent node fail low with value <= alpha and try another move.
        if (PvNode && (moveCount == 1 || value > alpha))
        {
            // Extend move from transposition table if we are about to dive into qsearch.
            if (move == ttMove && ss->ply <= pos->rootDepth * 2)
                newDepth = max(newDepth, 1);

            value = -search(pos, ss + 1, -beta, -alpha, newDepth, false, true);
        }

        // Step 16. Undo move
        undo_move(pos, move);

        // Step 17. Check for a new best move
        // Finished searching the move. If a stop occurred, the return value of
        // the search cannot be trusted, and we return immediately without
        // updating best move, PV and TT.
        if (Thread.stop)
            return 0;

        if (rootNode)
        {
            if (moveCount == 1 || value > alpha)
            {
                ss->pv.score = value;

                // We record how often the best move has been changed in each
                // iteration. This information is used for time management.
                if (moveCount > 1)
                    pos->bestMoveChanges++;
            }
        }

        if (value > bestValue)
        {
            bestValue = value;

            if (value > alpha)
            {
                bestMove = move;

                if (PvNode && ss->ply < 5)
                {
                    ss->pv.line[0] = move;
                    ss->pv.length  = (ss + 1)->pv.length + 1;
                    memcpy(ss->pv.line + 1, (ss + 1)->pv.line, sizeof(Move) * (ss + 1)->pv.length);
                }

                if (value >= beta)
                {
                    ss->cutoffCnt += !ttMove + (extension < 2);
                    break;
                }

                alpha = value;
            }
        }

        if (move != bestMove && moveCount < 32)
        {
            if (capture)
                capturesSearched[captureCount++] = move;
            else
                quietsSearched[quietCount++] = move;
        }
    }

    // Step 18. Check for mate and stalemate
    // All legal moves have been searched and if there are no legal moves,
    // it must be a mate or a stalemate. If we are in a singular extension
    // search then return a fail low score.
    if (!moveCount)
        bestValue = excludedMove ? alpha : ss->checkersBB ? mated_in(ss->ply) : VALUE_DRAW;
    else if (bestMove)
    {
        // Quiet best move: update move sorting heuristics
        if (!capture_stage(pos, bestMove))
        {
            int bonus = stat_bonus(depth + (bestValue > beta + qb_v1));
            int malus = stat_malus(depth + (bestValue > beta + qb_v2));

            update_quiet_stats(pos, ss, bestMove, bonus);

#pragma clang loop unroll(disable)
            // Decrease all the other played quiet moves
            for (int i = 0; i < quietCount; i++)
            {
                history_update(*pos->mainHistory, stm(), quietsSearched[i], -malus);
                update_continuation_histories(ss, moved_piece(quietsSearched[i]),
                                              to_sq(quietsSearched[i]), -malus);
            }
        }

        update_capture_stats(pos, bestMove, capturesSearched, captureCount, depth + 1);

        // Extra penalty for a quiet TT or main killer move in previous ply when it gets refuted
        if ((prevSq != SQ_NONE && (ss - 1)->moveCount == 1
             || (ss - 1)->currentMove == (ss - 1)->killers[0])
            && !captured_piece())
            update_continuation_histories(ss - 1, piece_on(prevSq), prevSq, -stat_malus(depth + 1));
    }
    // Bonus for prior countermove that caused the fail low
    else if (!captured_piece() && prevSq != SQ_NONE)
    {
        int bonus = pcmb_v1 * (depth > 4) + pcmb_v4 * ((ss - 1)->moveCount > pcmb_v5)
                  + pcmb_v6 * (!ss->checkersBB && bestValue <= ss->staticEval - pcmb_v7)
                  + 119 * (!(ss - 1)->checkersBB && bestValue <= -(ss - 1)->staticEval - 83);

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

    // Write gathered information in transposition table. Note that the
    // static evaluation is saved as it was before correction history.
    if (!excludedMove)
        tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                 bestValue >= beta    ? BOUND_LOWER
                 : PvNode && bestMove ? BOUND_EXACT
                                      : BOUND_UPPER,
                 depth, bestMove, unadjustedStaticEval);

    // Adjust correction history
    if (!ss->checkersBB && (!bestMove || !capture_stage(pos, bestMove))
        && !(bestValue >= beta && bestValue <= ss->staticEval)
        && !(!bestMove && bestValue >= ss->staticEval))
    {
        update_correction_histories(pos, depth, bestValue - ss->staticEval);
    }

    return bestValue;
}

// Quiescence search function, which is called by the main search function
// with zero depth, or recursively with further decreasing depth per call.
Value qsearch(Position* pos, Stack* ss, Value alpha, Value beta, Depth depth) {
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, bestMove;
    Value    bestValue, value, unadjustedStaticEval, ttValue, futilityValue, futilityBase;
    bool     ttHit, pvHit, givesCheck;
    Depth    ttDepth;
    int      moveCount;

    // Step 1. Initialize node
    bestMove  = 0;
    moveCount = 0;

    // Step 2. Check for an immediate draw or maximum ply reached
    if (is_draw(pos) || ss->ply >= MAX_PLY)
        return ss->ply >= MAX_PLY && !ss->checkersBB ? evaluate(pos) : VALUE_DRAW;

    // Decide whether to include checks: this fixes also the type of
    // TT entry depth that we are going to use. Note that in qsearch we use
    // only two types of depth in TT: DEPTH_QS_CHECKS or DEPTH_QS_NO_CHECKS.
    ttDepth = ss->checkersBB || depth >= DEPTH_QS_CHECKS ? DEPTH_QS_CHECKS : DEPTH_QS_NO_CHECKS;

    // Step 3. Transposition table lookup
    posKey  = key();
    tte     = tt_probe(posKey, &ttHit);
    ttValue = ttHit ? value_from_tt(tte_value(tte), ss->ply, rule50_count()) : VALUE_NONE;
    ttMove  = ttHit ? tte_move(tte) : 0;
    pvHit   = ttHit && tte_is_pv(tte);

    if (ttHit && tte_depth(tte) >= ttDepth
        && tte_bound(tte) & (ttValue >= beta ? BOUND_LOWER : BOUND_UPPER))
        return ttValue;

    // Step 4. Static evaluation of the position
    if (ss->checkersBB)
    {
        unadjustedStaticEval = VALUE_NONE;
        ss->staticEval       = VALUE_NONE;
        bestValue = futilityBase = -VALUE_INFINITE;
    }
    else
    {
        if (ttHit)
        {
            // Never assume anything about values stored in TT
            if ((unadjustedStaticEval = tte_eval(tte)) == VALUE_NONE)
                unadjustedStaticEval = evaluate(pos);

            ss->staticEval = bestValue = to_corrected(pos, unadjustedStaticEval);

            // ttValue can be used as a better position evaluation
            if (ttValue != VALUE_NONE
                && tte_bound(tte) & (ttValue > bestValue ? BOUND_LOWER : BOUND_UPPER))
                bestValue = ttValue;
        }
        else
        {
            unadjustedStaticEval =
              (ss - 1)->currentMove != MOVE_NULL ? evaluate(pos) : -(ss - 1)->staticEval + tempo;

            ss->staticEval = bestValue = to_corrected(pos, unadjustedStaticEval);
        }

        // Stand pat. Return immediately if static value is at least beta
        if (bestValue >= beta)
        {
            if (!ttHit)
                tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), false, BOUND_LOWER,
                         DEPTH_NONE, 0, unadjustedStaticEval);

            return bestValue;
        }

        if (bestValue > alpha)
            alpha = bestValue;

        futilityBase = ss->staticEval + qsf_v1;
    }

    ss->continuationHistory = &(*pos->contHist)[0][0][0];

    Square prevSq = move_is_ok((ss - 1)->currentMove) ? to_sq((ss - 1)->currentMove) : SQ_NONE;

    // Initialize move picker data for the current position, and prepare to search
    // the moves. Because the depth is <= 0 here, only captures, queen promotions
    // and checks (only if depth >= DEPTH_QS_CHECKS) will be generated.
    mp_init_q(pos, ttMove, depth, prevSq);

    // Step 5. Loop through the moves until no moves remain or a beta cutoff occurs
    while ((move = next_move(pos, 0)))
    {
        // Check for legality just before making the move
        if (!is_legal(pos, move))
            continue;

        givesCheck = gives_check(pos, ss, move);
        moveCount++;

        // Step 6. Pruning
        if (bestValue > VALUE_MATED_IN_MAX_PLY)
        {
            // Futility pruning
            if (!givesCheck && futilityBase > -VALUE_MATE_IN_MAX_PLY
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
            if (!see_test(pos, move, 0))
                continue;
        }

        PieceType movedPiece = type_of_p(moved_piece(move));

        // Step 7. Make and search the move
        do_move(pos, move, givesCheck);

        // Speculative prefetch as early as possible
        prefetch(tt_first_entry(key()));

        ss->currentMove         = move;
        ss->continuationHistory = &(*pos->contHist)[stm()][movedPiece][to_sq(move)];

        value = -qsearch(pos, ss + 1, -beta, -alpha, depth - 1);
        undo_move(pos, move);

        // Step 8. Check for a new best move
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

    // Step 9. Check for mate
    // All legal moves have been searched. A special case: if we are
    // in check and no legal moves were found, it is checkmate.
    if (ss->checkersBB && bestValue == -VALUE_INFINITE)
        return mated_in(ss->ply);  // Plies to mate from the root

    if (abs(bestValue) < VALUE_MATE_IN_MAX_PLY && bestValue >= beta)
        bestValue = (3 * bestValue + beta) / 4;

    // Save gathered info in transposition table. The static evaluation
    // is saved as it was before adjustment by correction history.
    tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), pvHit,
             bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, ttDepth, bestMove,
             unadjustedStaticEval);

    return bestValue;
}


// Adjusts a mate score from "plies to mate from the root" to
// "plies to mate from the current position". Non-mate scores are unchanged.
// The function is called before storing a value in the transposition table.
static Value value_to_tt(Value v, int ply) {
    return v >= VALUE_MATE_IN_MAX_PLY ? v + ply : v <= VALUE_MATED_IN_MAX_PLY ? v - ply : v;
}


// Inverse of value_to_tt(): It adjusts a mate score from the transposition
// table (which refers to the plies to mate/be mated from current position)
// to "plies to mate/be mated from the root".
static Value value_from_tt(Value v, int ply, int r50c) {
    if (v >= VALUE_MATE_IN_MAX_PLY)
        return (VALUE_MATE - v > 99 - r50c) ? VALUE_MATE_IN_MAX_PLY - 1 : v - ply;

    if (v <= VALUE_MATED_IN_MAX_PLY)
        return (VALUE_MATE + v > 99 - r50c) ? VALUE_MATED_IN_MAX_PLY + 1 : v + ply;

    return v;
}


static void update_correction_histories(const Position* pos, Depth depth, int32_t diff) {
    Key keys[] = {material_key(), pawn_key(), prev_move_key(), w_nonpawn_key(), b_nonpawn_key()};

#pragma clang loop unroll(disable)
    for (size_t i = 0; i < CORRECTION_HISTORY_NB; i++)
    {
        int16_t* entry = &(*pos->corrHists)[stm()][i][keys[i] & CORRECTION_HISTORY_MASK];

        int32_t newWeight  = min(ch_v1, 1 + depth);
        int32_t scaledDiff = diff * ch_v2;
        int32_t update     = *entry * (ch_v3 - newWeight) + scaledDiff * newWeight;

        *entry = max(-CORRECTION_HISTORY_MAX, min(CORRECTION_HISTORY_MAX, update / ch_v3));
    }
}

Value to_corrected(Position* pos, Value unadjustedStaticEval) {
    Key keys[]    = {material_key(), pawn_key(), prev_move_key(), w_nonpawn_key(), b_nonpawn_key()};
    int weights[] = {ch_v4, ch_v5, ch_v6, ch_v7, ch_v8};

    int32_t correction = 0;
    for (size_t i = 0; i < CORRECTION_HISTORY_NB; i++)
        correction += weights[i] * (*pos->corrHists)[stm()][i][keys[i] & CORRECTION_HISTORY_MASK];

    Value v = unadjustedStaticEval + correction / 128 / ch_v2;

    // Damp down the evaluation linearly when shuffling
    v = v * (100 - rule50_count()) / 100;

    return clamp(v, -VALUE_MATE_IN_MAX_PLY, VALUE_MATE_IN_MAX_PLY);
}

// Updates histories of the move pairs formed by moves
// at ply -1, -2, -4, and -6 with current move.
static void update_continuation_histories(Stack* ss, Piece pc, Square s, int bonus) {
    PieceType pt = type_of_p(pc);

    if (move_is_ok((ss - 1)->currentMove))
        update_contHist(*(ss - 1)->continuationHistory, pt, s, cnht_v1 * bonus / 1024);

    if (move_is_ok((ss - 2)->currentMove))
        update_contHist(*(ss - 2)->continuationHistory, pt, s, cnht_v2 * bonus / 1024);

    if (ss->checkersBB)
        return;

    if (move_is_ok((ss - 4)->currentMove))
        update_contHist(*(ss - 4)->continuationHistory, pt, s, cnht_v3 * bonus / 1024);

    if (move_is_ok((ss - 6)->currentMove))
        update_contHist(*(ss - 6)->continuationHistory, pt, s, cnht_v4 * bonus / 1024);
}

// Updates move sorting heuristics when a new capture best move is found
static void
update_capture_stats(const Position* pos, Move move, Move* captures, int captureCnt, Depth depth) {
    Piece moved_piece = moved_piece(move);
    int   captured    = type_of_p(piece_on(to_sq(move)));

    if (capture_stage(pos, move))
        cpth_update(*pos->captureHistory, moved_piece, to_sq(move), captured, stat_bonus(depth));

    Value malus = -stat_malus(depth);

#pragma clang loop unroll(disable)
    for (int i = 0; i < captureCnt; i++)
    {
        moved_piece = moved_piece(captures[i]);
        captured    = type_of_p(piece_on(to_sq(captures[i])));
        cpth_update(*pos->captureHistory, moved_piece, to_sq(captures[i]), captured, malus);
    }
}

// Updates move sorting heuristics when a new quiet best move is found
static void update_quiet_stats(const Position* pos, Stack* ss, Move move, int bonus) {
    if (ss->killers[0] != move)
    {
        ss->killers[1] = ss->killers[0];
        ss->killers[0] = move;
    }

    Color c = stm();
    history_update(*pos->mainHistory, c, move, bonus);
    update_continuation_histories(ss, moved_piece(move), to_sq(move), bonus);
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

// Used to detect when we are out of available time and thus stop the search.
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

// Prints PV information according to the UCI protocol.
static void uci_print_pv(Position* pos, Depth depth) {
    TimePoint   elapsed        = time_elapsed() + 1;
    PVariation* pv             = &pos->st->pv;
    uint64_t    nodes_searched = Thread.pos->nodes;
    char        buf[16];

    printf("info depth %d score %s nodes %" PRIu64 " nps %" PRIu64 " time %" PRIi64 " pv", depth,
           uci_value(buf, pv->score), nodes_searched, nodes_searched * 1000 / elapsed, elapsed);
#pragma clang loop unroll(disable)
    for (int i = 0; i < pv->length; i++)
        printf(" %s", uci_move(buf, pv->line[i]));
    printf("\n");

    fflush(stdout);
}


// Wakes up the main thread to start a new search, then returns immediately.
void start_thinking(Position* root) {
    prepare_for_search(root);
    mainthread_search();
}

SMALL void prepare_for_search(Position* root) {
    Thread.stop = false;

    Position* pos  = Thread.pos;
    pos->rootDepth = 0;
    pos->nodes     = 0;

    root->st->pv.length = 0;
    memcpy(pos, root, offsetof(Position, moveList));

    // Copy enough of the root State buffer.
    int n = max(7, root->st->pliesFromNull);
#pragma clang loop unroll(disable)
    for (int i = 0; i <= n; i++)
        memcpy(&pos->stack[i], &root->st[i - n], StateSize);
    pos->st                 = pos->stack + n;
    (pos->st - 1)->endMoves = pos->moveList;

    pos_set_check_info(pos);
}
