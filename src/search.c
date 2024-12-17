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

#define load_rlx(x) atomic_load_explicit(&(x), memory_order_relaxed)
#define store_rlx(x, y) atomic_store_explicit(&(x), y, memory_order_relaxed)

int nmp_v1     = 773;
int nmp_v2     = 72;
int nmp_v3     = 193;
int nmp_v4     = 206;
int nmp_v5     = 23275;
int nmp_v6     = 29;
int nmp_v7     = 26;
int nmp_v8     = 88;
int nmp_v9     = 164;
int lph_v1     = 1149;
int lph_v2     = 562;
int qmo_v1     = 380;
int qmo_v2     = 1081;
int qmo_v3     = 1040;
int rz_v1      = 517;
int ft_v1      = 217;
int rd_v1      = 538;
int rd_v2      = 1151;
int rd_v3      = 888;
int sb_v1      = 13;
int sb_v2      = 30;
int sb_v3      = 17;
int sb_v4      = 136;
int rd_init_v1 = 2319;
int d_v1       = 21;
int iir_v1     = 614;
int iir_v2     = 202;
int cbp_v1     = 428;
int cbp_v2     = 6;
int fpp_v1     = 690;
int fpp_v2     = 298;
int fpp_v3     = 168;
int fpp_v4     = 27912;
int sqsee_v1   = 2937;
int sqsee_v2   = 1790;
int sch_v1     = 99;
int sch_v2     = -1;
int sfpc_v1    = 617;
int sfpc_v2    = 179;
int sfpc_v3    = 159;
int sfpc_v4    = 238;
int scsee_v1   = 219;
int se_v1      = 654;
int se_v2      = 372;
int se_v3      = 229;
int se_v4      = 283;
int se_v5      = 7696;
int se_v6      = 101;
int prb_v1     = 173;
int prb_v2     = 51;
int rfp_v1     = 840;
int lmr_v1     = 1019;
int lmr_v2     = 1287;
int lmr_v3     = 4593;
int lmr_v4     = 98;
int lmr_v5     = 98;
int lmr_v6     = 127;
int lmr_v7     = 138;
int lmr_v8     = 15862;
int fmc_v1     = 293;
int fmc_v2     = 285;
int fmc_v3     = 203;
int asd_v1     = 464;
int ses_v1     = 323;
int qsf_v1     = 151;
int ch_v1      = 1760;
int ch_v2      = 261;
int ch_v3      = 292;
int tempo      = 59;
int mp_v1      = 67;
int mp_v2      = 1099;
int mp_v3      = 3041;
int mp_v4      = 110;
int mp_v5      = 209;
int mp_v6      = 204;
int mp_v7      = 195;
int mp_v8      = 107;
int mp_v9      = 98;
int mp_v10     = 377;
int mp_v11     = 251;
int mg_pawn    = 128;
int eg_pawn    = 212;
int mg_knight  = 626;
int eg_knight  = 911;
int mg_bishop  = 797;
int eg_bishop  = 988;
int mg_rook    = 1263;
int eg_rook    = 1487;
int mg_queen   = 2564;
int eg_queen   = 2740;
int eval_scale = 95;

LimitsType Limits;

extern char lastFen[256];

static int base_ct;

// Different node types, used as template parameter
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
static Value stat_bonus(Depth depth) { return min((6 * depth + 229) * depth - 215, 2000); }

// Add a small random component to draw evaluations to keep search dynamic
// and to avoid three-fold blindness. (Yucks, ugly hack)
static Value value_draw(Position* pos) { return VALUE_DRAW + 2 * (pos->nodes & 1) - 1; }

static Value   value_to_tt(Value v, int ply);
static Value   value_from_tt(Value v, int ply, int r50c);
static void    update_pv(Move* pv, Move move, Move* childPv);
static void    update_cm_stats(Stack* ss, Piece pc, Square s, int bonus);
static int32_t get_correction(correction_history_t hist, Color side, Key materialKey);
static void    add_correction_history(
     correction_history_t hist, Color side, Key materialKey, Depth depth, int32_t diff);
static void update_quiet_stats(const Position* pos, Stack* ss, Move move, int bonus);
static void
update_capture_stats(const Position* pos, Move move, Move* captures, int captureCnt, int bonus);
static void check_time(void);
static void stable_sort(RootMove* rm, int num);
static void uci_print_pv(Position* pos, Depth depth, Value alpha, Value beta);

// search_init() is called during startup to initialize various lookup tables

double my_log(double x) {
    if (x <= 0)
        return -1e9;  // Handle log(0) or negative input
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

void search_init(void) {
    for (int i = 1; i < MAX_MOVES; i++)
        Reductions[i] = rd_init_v1 / 100.0 * my_log(i);
}


// search_clear() resets search state to zero, to obtain reproducible results

void search_clear(void) {
    if (!settings.ttSize)
    {
        delayedSettings.clear = true;
        return;
    }

    Time.availableNodes = 0;

    tt_clear();
    for (int i = 0; i < numCmhTables; i++)
        if (cmhTables[i])
        {
            stats_clear(cmhTables[i]);
            for (int c = 0; c < 2; c++)
                for (int j = 0; j < 16; j++)
                    for (int k = 0; k < 64; k++)
                        (*cmhTables[i])[c][0][j][k] = CounterMovePruneThreshold - 1;
        }

    Position* pos = Threads.pos[0];
    stats_clear(pos->counterMoves);
    stats_clear(pos->history);
    stats_clear(pos->captureHistory);
    stats_clear(pos->corrHistory);

    mainThread.previousScore         = VALUE_INFINITE;
    mainThread.previousTimeReduction = 1;
}

// mainthread_search() is called by the main thread when the program
// receives the UCI 'go' command. It searches from the root position and
// outputs the "bestmove".

void mainthread_search(void) {
    Position* pos = Threads.pos[0];
    Color     us  = stm();
    time_init(us, game_ply());
    tt_new_search();
    char buf[16];

    if (pos->rootMoves->size > 0)
    {
        Move bookMove = 0;


        for (int i = 0; i < pos->rootMoves->size; i++)
            if (pos->rootMoves->move[i].pv[0] == bookMove)
            {
                RootMove tmp            = pos->rootMoves->move[0];
                pos->rootMoves->move[0] = pos->rootMoves->move[i];
                pos->rootMoves->move[i] = tmp;
                break;
            }

        Threads.pos[0]->bestMoveChanges = 0;
        thread_search(pos);  // Let's start searching!
    }

    if (pos->rootMoves->size == 0)
    {
        pos->rootMoves->move[0].pv[0]  = 0;
        pos->rootMoves->move[0].pvSize = 1;
        pos->rootMoves->size++;
        printf("info depth 0 score %s\n", uci_value(buf, checkers() ? -VALUE_MATE : VALUE_DRAW));
        fflush(stdout);
    }

    mainThread.previousScore = pos->rootMoves->move[0].score;

    printf("bestmove %s", uci_move(buf, pos->rootMoves->move[0].pv[0]));
    printf("\n");
    fflush(stdout);

#ifdef KAGGLE
    // Start pondering right after the best move has been printed if we can
    if (pos->rootMoves->move[0].pvSize >= 2)
    {
        Threads.ponder = true;
        Threads.stop   = false;

        const Move bestMove = pos->rootMoves->move[0].pv[0];
        const Move ponder   = pos->rootMoves->move[0].pv[1];

        char command[2048];
        snprintf(command, 2048, "%s moves %s %s", lastFen, uci_move(buf, bestMove),
                 uci_move(buf, ponder));

        position(pos, command);

        prepare_for_search(pos, true);
        thread_search(pos);

        Threads.ponder = false;
        Threads.stop   = true;
    }
#endif
}


// thread_search() is the main iterative deepening loop. It calls search()
// repeatedly with increasing depth until the allocated thinking time has
// been consumed, the user stops the search, or the maximum search depth is
// reached.

void thread_search(Position* pos) {
    Value  bestValue, alpha, beta, delta;
    Move   pv[MAX_PLY + 1];
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
        ss[i].history    = &(*pos->counterMoveHistory)[0][0];  // Use as sentinel
        ss[i].staticEval = VALUE_NONE;
    }

    for (int i = 0; i <= MAX_PLY; i++)
        ss[i].ply = i;
    ss->pv = pv;

    ss->accumulator.needs_refresh = 1;

    bestValue = delta = alpha = -VALUE_INFINITE;
    beta                      = VALUE_INFINITE;
    pos->completedDepth       = 0;

    if (pos->threadIdx == 0)
    {
        if (mainThread.previousScore == VALUE_INFINITE)
            for (int i = 0; i < 4; i++)
                mainThread.iterValue[i] = VALUE_ZERO;
        else
            for (int i = 0; i < 4; i++)
                mainThread.iterValue[i] = mainThread.previousScore;
    }

    RootMoves* rm                 = pos->rootMoves;
    int        searchAgainCounter = 0;

    // Iterative deepening loop until requested to stop or the target depth
    // is reached.
    while ((pos->rootDepth += 2) < MAX_PLY && !Threads.stop
           && !(Limits.depth && pos->threadIdx == 0 && pos->rootDepth > Limits.depth))
    {
        // Age out PV variability metric
        if (pos->threadIdx == 0)
            totBestMoveChanges /= 2;

        // Save the last iteration's scores before first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
        for (int idx = 0; idx < rm->size; idx++)
            rm->move[idx].previousScore = rm->move[idx].score;

        int pvFirst = 0, pvLast = 0;

        if (!Threads.increaseDepth)
            searchAgainCounter++;

        int pvIdx = 0;

        pos->pvIdx = pvIdx;
        if (pvIdx == pvLast)
        {
            pvFirst = pvLast;
            for (pvLast++; pvLast < rm->size; pvLast++)
                if (rm->move[pvLast].tbRank != rm->move[pvFirst].tbRank)
                    break;
            pos->pvLast = pvLast;
        }

        // Reset aspiration window starting size
        if (pos->rootDepth >= 4)
        {
            Value previousScore = rm->move[pvIdx].previousScore;
            delta               = d_v1;
            alpha               = max(previousScore - delta, -VALUE_INFINITE);
            beta                = min(previousScore + delta, VALUE_INFINITE);
        }

        // Start with a small aspiration window and, in the case of a fail
        // high/low, re-search with a bigger window until we're not failing
        // high/low anymore.
        int failedHighCnt = 0;
        while (true)
        {
            Depth adjustedDepth = max(1, pos->rootDepth - failedHighCnt - searchAgainCounter);
            bestValue           = search(pos, ss, alpha, beta, adjustedDepth, false, true);

            // Bring the best move to the front. It is critical that sorting
            // is done with a stable algorithm because all the values but the
            // first and eventually the new best one are set to -VALUE_INFINITE
            // and we want to keep the same order for all the moves except the
            // new PV that goes to the front. Note that in case of MultiPV
            // search the already searched PV lines are preserved.
            stable_sort(&rm->move[pvIdx], pvLast - pvIdx);

            // If search has been stopped, we break immediately. Sorting and
            // writing PV back to TT is safe because RootMoves is still
            // valid, although it refers to the previous iteration.
            if (Threads.stop)
                break;

            // In case of failing low/high increase aspiration window and
            // re-search, otherwise exit the loop.
            if (bestValue <= alpha)
            {
                beta  = (alpha + beta) / 2;
                alpha = max(bestValue - delta, -VALUE_INFINITE);

                failedHighCnt = 0;
            }
            else
            {
                rm->move[pvIdx].bestMoveCount++;
                break;
            }

            delta += delta / 4 + asd_v1 / 100;
        }

        // Sort the PV lines searched so far and update the GUI
        stable_sort(&rm->move[pvFirst], pvIdx - pvFirst + 1);

#ifndef KAGGLE
        uci_print_pv(pos, pos->rootDepth, alpha, beta);
#endif

        if (!Threads.stop)
            pos->completedDepth = pos->rootDepth;

        if (rm->move[0].pv[0] != lastBestMove)
        {
            lastBestMove      = rm->move[0].pv[0];
            lastBestMoveDepth = pos->rootDepth;
        }

        if (pos->threadIdx != 0)
            continue;

        // Do we have time for the next iteration? Can we stop searching now?
        if (use_time_management() && !Threads.stop)
        {
            double fallingEval = (318 + 6 * (mainThread.previousScore - bestValue)
                                  + 6 * (mainThread.iterValue[iterIdx] - bestValue))
                               / 825.0;
            fallingEval = clamp(fallingEval, 0.5, 1.5);

            // If the best move is stable over several iterations, reduce time
            // accordingly
            timeReduction    = lastBestMoveDepth + 9 < pos->completedDepth ? 1.92 : 0.95;
            double reduction = (1.47 + mainThread.previousTimeReduction) / (2.32 * timeReduction);

            // Use part of the gained time from a previous stable move for this move
            totBestMoveChanges += Threads.pos[0]->bestMoveChanges;
            Threads.pos[0]->bestMoveChanges = 0;

            double bestMoveInstability = 1 + totBestMoveChanges;

            double totalTime =
              rm->size == 1 ? 0 : time_optimum() * fallingEval * reduction * bestMoveInstability;

            // Stop the search if we have exceeded the totalTime (at least 1ms)
            if (time_elapsed() > totalTime)
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "stop".
                if (Threads.ponder)
                {}
                else
                    Threads.stop = true;
            }
            else
                Threads.increaseDepth =
                  !(Threads.increaseDepth && !Threads.ponder && time_elapsed() > totalTime * 0.58);
        }

        mainThread.iterValue[iterIdx] = bestValue;
        iterIdx                       = (iterIdx + 1) & 3;
    }

    if (pos->threadIdx != 0)
        return;

    mainThread.previousTimeReduction = timeReduction;
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

    Move     pv[MAX_PLY + 1], capturesSearched[32], quietsSearched[64];
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, excludedMove, bestMove;
    Depth    extension, newDepth;
    Value    bestValue, value, ttValue, eval, rawEval, probCutBeta;
    bool     ttHit, formerPv, givesCheck, improving;
    bool     captureOrPromotion, inCheck, moveCountPruning;
    bool     ttCapture, singularQuietLMR;
    Piece    movedPiece;
    int      moveCount, captureCount, quietCount;

    // Step 1. Initialize node
    inCheck   = checkers();
    moveCount = captureCount = quietCount = ss->moveCount = 0;
    bestValue                                             = -VALUE_INFINITE;

    // Check for the available remaining time
    if (load_rlx(pos->resetCalls))
    {
        store_rlx(pos->resetCalls, false);
        pos->callsCnt = Limits.nodes ? min(1024, Limits.nodes / 1024) : 1024;
    }
    if (--pos->callsCnt <= 0)
    {
        store_rlx(Threads.pos[0]->resetCalls, true);
        check_time();
    }

    if (!rootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (load_rlx(Threads.stop) || is_draw(pos) || ss->ply >= MAX_PLY)
            return ss->ply >= MAX_PLY && !inCheck ? evaluate(pos) : value_draw(pos);
    }


    (ss + 1)->ttPv         = false;
    (ss + 1)->excludedMove = bestMove = 0;
    (ss + 2)->killers[0] = (ss + 2)->killers[1] = 0;
    (ss + 2)->cutoffCnt                         = 0;
    Square prevSq                               = to_sq((ss - 1)->currentMove);

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
    ttMove       = rootNode ? pos->rootMoves->move[pos->pvIdx].pv[0] : ttHit ? tte_move(tte) : 0;
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
                if ((ss - 1)->moveCount <= 2 && !captured_piece())
                    update_cm_stats(ss - 1, piece_on(prevSq), prevSq, -stat_bonus(depth + 1));
            }
            // Penalty for a quiet ttMove that fails low
            else if (!is_capture_or_promotion(pos, ttMove))
            {
                int penalty = -stat_bonus(depth);
                history_update(*pos->history, stm(), ttMove, penalty);
                update_cm_stats(ss, moved_piece(ttMove), to_sq(ttMove), penalty);
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

        eval = ss->staticEval = rawEval + get_correction(pos->corrHistory, stm(), material_key());

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

        eval = ss->staticEval = rawEval + get_correction(pos->corrHistory, stm(), material_key());

        tte_save(tte, posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_NONE, 0, rawEval);
    }

    // Step 7. Razoring
    if (!rootNode && depth == 1 && eval <= alpha - rz_v1)
        return qsearch(pos, ss, alpha, beta, 0, PvNode, false);


    improving = (ss - 2)->staticEval == VALUE_NONE
                ? (ss->staticEval > (ss - 4)->staticEval || (ss - 4)->staticEval == VALUE_NONE)
                : ss->staticEval > (ss - 2)->staticEval;

    if (move_is_ok((ss - 1)->currentMove) && !(ss - 1)->checkersBB && !captured_piece())
    {
        int bonus = clamp(-depth * qmo_v1 / 100 * ((ss - 1)->staticEval + ss->staticEval - tempo),
                          -qmo_v2, qmo_v3);
        history_update(*pos->history, !stm(), (ss - 1)->currentMove, bonus);
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

        ss->currentMove = MOVE_NULL;
        ss->history     = &(*pos->counterMoveHistory)[0][0];

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
    if (!PvNode && depth > 4 && abs(beta) < VALUE_TB_WIN_IN_MAX_PLY
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

                ss->currentMove = move;
                ss->history     = &(*pos->counterMoveHistory)[moved_piece(move)][to_sq(move)];
                givesCheck      = gives_check(pos, ss, move);
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
    PieceToHistory* cmh  = (ss - 1)->history;
    PieceToHistory* fmh  = (ss - 2)->history;
    PieceToHistory* fmh2 = (ss - 4)->history;
    PieceToHistory* fmh3 = (ss - 6)->history;

    mp_init(pos, ttMove, depth, ss->ply);

    value            = bestValue;
    singularQuietLMR = moveCountPruning = false;
    ttCapture                           = ttMove && is_capture_or_promotion(pos, ttMove);
    formerPv                            = ss->ttPv && !PvNode;

    // Step 12. Loop through moves
    // Loop through all pseudo-legal moves until no moves remain or a beta
    // cutoff occurs
    while ((move = next_move(pos, moveCountPruning)))
    {

        if (move == excludedMove)
            continue;

        // At root obey the "searchmoves" option and skip moves not listed
        // inRoot Move List. As a consequence any illegal move is also skipped.
        // In MultiPV mode we also skip PV moves which have been already
        // searched.
        if (rootNode)
        {
            int idx;
            for (idx = pos->pvIdx; idx < pos->pvLast; idx++)
                if (pos->rootMoves->move[idx].pv[0] == move)
                    break;
            if (idx == pos->pvLast)
                continue;
        }

        // Check for legality just before making the move
        if (!rootNode && !is_legal(pos, move))
            continue;

        ss->moveCount = ++moveCount;

        if (PvNode)
            (ss + 1)->pv = NULL;

        extension          = 0;
        captureOrPromotion = is_capture_or_promotion(pos, move);
        movedPiece         = moved_piece(move);

        givesCheck = gives_check(pos, ss, move);

        // Calculate new depth for this move
        newDepth = depth - 1;

        // Step 13. Pruning at shallow depth
        if (!rootNode && non_pawn_material_c(stm()) && bestValue > VALUE_TB_LOSS_IN_MAX_PLY)
        {
            // Skip quiet moves if movecount exceeds our FutilityMoveCount threshold
            moveCountPruning = moveCount >= futility_move_count(improving, depth);

            // Reduced depth of the next LMR search
            int lmrDepth = max(newDepth - reduction(improving, depth, moveCount), 0);

            if (!captureOrPromotion && !givesCheck)
            {
                // Countermoves based pruning
                if (lmrDepth < cbp_v1 / 100
                                 + ((ss - 1)->statScore > cbp_v2 / 100 || (ss - 1)->moveCount == 1)
                    && (*cmh)[movedPiece][to_sq(move)] < CounterMovePruneThreshold
                    && (*fmh)[movedPiece][to_sq(move)] < CounterMovePruneThreshold)
                    continue;

                // Futility pruning: parent node
                if (lmrDepth < fpp_v1 / 100 && !inCheck
                    && ss->staticEval + fpp_v2 + fpp_v3 * lmrDepth <= alpha
                    && (*cmh)[movedPiece][to_sq(move)] + (*fmh)[movedPiece][to_sq(move)]
                           + (*fmh2)[movedPiece][to_sq(move)] + (*fmh3)[movedPiece][to_sq(move)] / 2
                         < fpp_v4)
                    continue;

                // Prune moves with negative SEE at low depths and below a decreasing
                // threshold at higher depths.
                if (!see_test(pos, move,
                              -(sqsee_v1 / 100 - min(lmrDepth, sqsee_v2 / 100)) * lmrDepth
                                * lmrDepth))
                    continue;
            }
            else
            {
                // Capture history based pruning when the move doesn't give check
                if (!givesCheck && lmrDepth < sch_v1 / 100
                    && (*pos->captureHistory)[movedPiece][to_sq(move)]
                                             [type_of_p(piece_on(to_sq(move)))]
                         < sch_v2)
                    continue;

                // Futility pruning for captures
                if (!givesCheck && lmrDepth < sfpc_v1 / 100
                    && !(PvNode && abs(bestValue) < sfpc_v2 / 100)
                    && *PieceValue[MG][type_of_p(movedPiece)]
                         >= *PieceValue[MG][type_of_p(piece_on(to_sq(move)))]
                    && !inCheck
                    && ss->staticEval + sfpc_v3 + sfpc_v4 * lmrDepth
                           + *PieceValue[MG][type_of_p(piece_on(to_sq(move)))]
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
            Value singularBeta  = ttValue - ((formerPv + se_v2 / 100) * depth) / (se_v3 / 100);
            Depth singularDepth = (depth - 1 + se_v4 / 100 * formerPv) / 2;
            ss->excludedMove    = move;
            Move cm             = ss->countermove;
            Move k1 = ss->mpKillers[0], k2 = ss->mpKillers[1];
            value = search(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode, false);
            ss->excludedMove = 0;

            if (value < singularBeta)
            {
                extension = 1;
                if (!PvNode && value < singularBeta - se_v5 / 100)
                    extension = 2;

                singularQuietLMR = !ttCapture;
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
                ss->countermove  = cm;  // pedantic
                ss->mpKillers[0] = k1;
                ss->mpKillers[1] = k2;

                ss->excludedMove = move;
                value            = search(pos, ss, beta - 1, beta, (depth + 3) / 2, cutNode, false);
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
            ss->countermove  = cm;  // pedantic
            ss->mpKillers[0] = k1;
            ss->mpKillers[1] = k2;
        }

        // Last capture extension
        else if (*PieceValue[EG][captured_piece()] > eg_pawn && non_pawn_material() <= 2 * mg_rook)
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
        ss->currentMove = move;
        ss->history     = &(*pos->counterMoveHistory)[movedPiece][to_sq(move)];

        // Step 15. Make the move.
        do_move(pos, move, givesCheck);
        // HACK: Fix bench after introduction of 2-fold MultiPV bug
        if (rootNode)
            pos->st[-1].key ^= pos->rootKeyFlip;

        // Step 16. Reduced depth search (LMR). If the move fails high it will be
        // re-searched at full depth.
        if (depth >= 3 && moveCount > 1 + 2 * rootNode
            && (!captureOrPromotion || cutNode || (!PvNode && !formerPv)))
        {
            Depth r = reduction(improving, depth, moveCount);

            // Decrease reduction at non-check cut nodes for second move at low
            // depths
            if (cutNode && depth <= lmr_v1 / 100 && moveCount <= 2 && !inCheck)
                r--;

            // Decrease reduction if position is or has been on the PV
            if (ss->ttPv)
                r -= 2;

            if (moveCountPruning && !formerPv)
                r++;

            // Decrease reduction if opponent's move count is high
            if ((ss - 1)->moveCount > lmr_v2 / 100)
                r--;

            // Decrease reduction if ttMove has been singularly extended
            if (singularQuietLMR)
                r -= 1 + formerPv;

            if (!captureOrPromotion)
            {
                // Increase reduction if ttMove is a capture
                if (ttCapture)
                    r++;

                if ((ss + 1)->cutoffCnt > 3)
                    r++;

                // Increase reduction for cut nodes
                if (cutNode)
                    r += 2;

                ss->statScore = (*cmh)[movedPiece][to_sq(move)] + (*fmh)[movedPiece][to_sq(move)]
                              + (*fmh2)[movedPiece][to_sq(move)]
                              + (*pos->history)[!stm()][from_to(move)] - lmr_v3;

                // Decrease/increase reduction for moves with a good/bad history.
                r -= ss->statScore / lmr_v8;
            }

            Depth d = clamp(newDepth - r, 1, newDepth);
            value   = -search(pos, ss + 1, -(alpha + 1), -alpha, d, true, false);

            if (value > alpha && d != newDepth)
            {

                value = -search(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode, false);
                if (!captureOrPromotion)
                {
                    int bonus = value > alpha ? stat_bonus(newDepth) : -stat_bonus(newDepth);

                    if (move == ss->killers[0])
                        bonus += bonus / 4;

                    update_cm_stats(ss, movedPiece, to_sq(move), bonus);
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
        if (PvNode && (moveCount == 1 || (value > alpha && (rootNode || value < beta))))
        {
            (ss + 1)->pv    = pv;
            (ss + 1)->pv[0] = 0;

            // Extend move from transposition table if we are about to dive into qsearch.
            if (move == ttMove && ss->ply <= pos->rootDepth * 2)
                newDepth = max(newDepth, 1);

            value = -search(pos, ss + 1, -beta, -alpha, newDepth, false, true);
        }

        // Step 18. Undo move
        // HACK: Fix bench after introduction of 2-fold MultiPV bug
        if (rootNode)
            pos->st[-1].key ^= pos->rootKeyFlip;
        undo_move(pos, move);


        // Step 19. Check for a new best move
        // Finished searching the move. If a stop occurred, the return value of
        // the search cannot be trusted, and we return immediately without
        // updating best move, PV and TT.
        if (load_rlx(Threads.stop))
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

                ss->cutoffCnt += !ttMove + (extension < 2);

                if (PvNode && !rootNode)  // Update pv even in fail-high case
                    update_pv(ss->pv, move, (ss + 1)->pv);

                if (PvNode && value < beta)  // Update alpha! Always alpha < beta
                    alpha = value;
                else
                {
                    ss->statScore = 0;
                    break;
                }
            }
        }

        if (move != bestMove)
        {
            if (captureOrPromotion && captureCount < 32)
                capturesSearched[captureCount++] = move;

            else if (!captureOrPromotion && quietCount < 64)
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

    if (!PvNode && bestValue >= beta && abs(bestValue) < VALUE_TB_WIN_IN_MAX_PLY
        && abs(beta) < VALUE_TB_WIN_IN_MAX_PLY && abs(alpha) < VALUE_TB_WIN_IN_MAX_PLY)
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
            int bonus = stat_bonus(depth + (bestValue > beta + mg_pawn));
            update_quiet_stats(pos, ss, bestMove, bonus);

            // Decrease all the other played quiet moves
            for (int i = 0; i < quietCount; i++)
            {
                history_update(*pos->history, stm(), quietsSearched[i], -bonus);
                update_cm_stats(ss, moved_piece(quietsSearched[i]), to_sq(quietsSearched[i]),
                                -bonus);
            }
        }

        update_capture_stats(pos, bestMove, capturesSearched, captureCount, stat_bonus(depth + 1));

        // Extra penalty for a quiet TT or main killer move in previous ply
        // when it gets refuted
        if (((ss - 1)->moveCount == 1 || (ss - 1)->currentMove == (ss - 1)->killers[0])
            && !captured_piece())
            update_cm_stats(ss - 1, piece_on(prevSq), prevSq, -stat_bonus(depth + 1));
    }
    // Bonus for prior countermove that caused the fail low
    else if ((depth >= 3 || PvNode) && !captured_piece())
        update_cm_stats(ss - 1, piece_on(prevSq), prevSq, stat_bonus(depth));

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
        add_correction_history(pos->corrHistory, stm(), material_key(), depth,
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

    Move     pv[MAX_PLY + 1];
    TTEntry* tte;
    Key      posKey;
    Move     ttMove, move, bestMove;
    Value    bestValue, value, rawEval, ttValue, futilityValue, futilityBase, oldAlpha;
    bool     ttHit, pvHit, givesCheck;
    Depth    ttDepth;
    int      moveCount;

    if (PvNode)
    {
        oldAlpha     = alpha;  // To flag BOUND_EXACT when eval above alpha and no available moves
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

            ss->staticEval = bestValue =
              rawEval + get_correction(pos->corrHistory, stm(), material_key());


            // Can ttValue be used as a better position evaluation?
            if (ttValue != VALUE_NONE
                && (tte_bound(tte) & (ttValue > bestValue ? BOUND_LOWER : BOUND_UPPER)))
                bestValue = ttValue;
        }
        else
        {
            rawEval =
              (ss - 1)->currentMove != MOVE_NULL ? evaluate(pos) : -(ss - 1)->staticEval + tempo;

            ss->staticEval = bestValue =
              rawEval + get_correction(pos->corrHistory, stm(), material_key());
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

    ss->history = &(*pos->counterMoveHistory)[0][0];

    // Initialize move picker data for the current position, and prepare
    // to search the moves. Because the depth is <= 0 here, only captures,
    // queen promotions and checks (only if depth >= DEPTH_QS_CHECKS) will
    // be generated.
    mp_init_q(pos, ttMove, depth, to_sq((ss - 1)->currentMove));

    // Loop through the moves until no moves remain or a beta cutoff occurs
    while ((move = next_move(pos, 0)))
    {

        givesCheck = gives_check(pos, ss, move);

        moveCount++;

        // Futility pruning
        if (!InCheck && !givesCheck && futilityBase > -VALUE_KNOWN_WIN
            && !advanced_pawn_push(pos, move))
        {

            if (moveCount > 2)
                continue;

            futilityValue = futilityBase + *PieceValue[EG][piece_on(to_sq(move))];

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
        if (!InCheck && !see_test(pos, move, 0))
            continue;

        // Speculative prefetch as early as possible
        prefetch(tt_first_entry(key_after(pos, move)));

        // Check for legality just before making the move
        if (!is_legal(pos, move))
        {
            moveCount--;
            continue;
        }

        ss->currentMove         = move;
        bool captureOrPromotion = is_capture_or_promotion(pos, move);
        ss->history             = &(*pos->counterMoveHistory)[moved_piece(move)][to_sq(move)];

        if (!captureOrPromotion && moveCount
            && (*(ss - 1)->history)[moved_piece(move)][to_sq(move)] < CounterMovePruneThreshold
            && (*(ss - 2)->history)[moved_piece(move)][to_sq(move)] < CounterMovePruneThreshold)
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

                if (PvNode)  // Update pv even in fail-high case
                    update_pv(ss->pv, move, (ss + 1)->pv);

                if (PvNode && value < beta)  // Update alpha here!
                    alpha = value;
                else
                    break;  // Fail high
            }
        }
    }

    // All legal moves have been searched. A special case: If we're in check
    // and no legal moves were found, it is checkmate.
    if (InCheck && bestValue == -VALUE_INFINITE)
        return mated_in(ss->ply);  // Plies to mate from the root

    if (abs(bestValue) < VALUE_TB_WIN_IN_MAX_PLY && bestValue >= beta)
        bestValue = (3 * bestValue + beta) / 4;

    tte_save(tte, posKey, value_to_tt(bestValue, ss->ply), pvHit,
             bestValue >= beta                ? BOUND_LOWER
             : PvNode && bestValue > oldAlpha ? BOUND_EXACT
                                              : BOUND_UPPER,
             ttDepth, bestMove, rawEval);

    return bestValue;
}

#define rm_lt(m1, m2) \
    ((m1).tbRank != (m2).tbRank ? (m1).tbRank < (m2).tbRank \
     : (m1).score != (m2).score ? (m1).score < (m2).score \
                                : (m1).previousScore < (m2).previousScore)

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

    return v >= VALUE_TB_WIN_IN_MAX_PLY ? v + ply : v <= VALUE_TB_LOSS_IN_MAX_PLY ? v - ply : v;
}


// value_from_tt() is the inverse of value_to_tt(): It adjusts a mate score
// from the transposition table (which refers to the plies to mate/be mated
// from current position) to "plies to mate/be mated from the root".

static Value value_from_tt(Value v, int ply, int r50c) {
    if (v == VALUE_NONE)
        return VALUE_NONE;

    if (v >= VALUE_TB_WIN_IN_MAX_PLY)
    {
        if (v >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - v > 99 - r50c)
            return VALUE_MATE_IN_MAX_PLY - 1;
        return v - ply;
    }

    if (v <= VALUE_TB_LOSS_IN_MAX_PLY)
    {
        if (v <= VALUE_MATED_IN_MAX_PLY && VALUE_MATE + v > 99 - r50c)
            return VALUE_MATED_IN_MAX_PLY + 1;
        return v + ply;
    }

    return v;
}

// update_pv() adds current move and appends child pv[]

static void update_pv(Move* pv, Move move, Move* childPv) {
    for (*pv++ = move; childPv && *childPv;)
        *pv++ = *childPv++;
    *pv = 0;
}

// differential.
static void add_correction_history(
  correction_history_t hist, Color side, Key materialKey, Depth depth, int32_t diff) {
    int32_t* entry      = &hist[side][materialKey % CORRECTION_HISTORY_ENTRY_NB];
    int32_t  newWeight  = min(ch_v1 / 100, 1 + depth);
    int32_t  scaledDiff = diff * ch_v2;
    int32_t  update     = *entry * (ch_v3 - newWeight) + scaledDiff * newWeight;
    // Clamp entry in-bounds.
    *entry = max(-CORRECTION_HISTORY_MAX, min(CORRECTION_HISTORY_MAX, update / ch_v3));
}

// Get the correction history differential for the given side and materialKey.
static int32_t get_correction(correction_history_t hist, Color side, Key materialKey) {
    return hist[side][materialKey % CORRECTION_HISTORY_ENTRY_NB] / ch_v2;
}

// update_cm_stats() updates countermove and follow-up move history.

static void update_cm_stats(Stack* ss, Piece pc, Square s, int bonus) {
    if (move_is_ok((ss - 1)->currentMove))
        cms_update(*(ss - 1)->history, pc, s, bonus);

    if (move_is_ok((ss - 2)->currentMove))
        cms_update(*(ss - 2)->history, pc, s, bonus);

    if (ss->checkersBB)
        return;

    if (move_is_ok((ss - 4)->currentMove))
        cms_update(*(ss - 4)->history, pc, s, bonus);

    if (move_is_ok((ss - 6)->currentMove))
        cms_update(*(ss - 6)->history, pc, s, bonus);
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
    history_update(*pos->history, c, move, bonus);
    update_cm_stats(ss, moved_piece(move), to_sq(move), bonus);

    if (type_of_p(moved_piece(move)) != PAWN)
        history_update(*pos->history, c, reverse_move(move), -bonus);

    if (move_is_ok((ss - 1)->currentMove))
    {
        Square prevSq                                  = to_sq((ss - 1)->currentMove);
        (*pos->counterMoves)[piece_on(prevSq)][prevSq] = move;
    }
}


// check_time() is used to print debug info and, more importantly, to detect
// when we are out of available time and thus stop the search.

static void check_time(void) {
    TimePoint elapsed = time_elapsed();

    // An engine may not stop pondering until told so by the GUI
    if (Threads.ponder)
        return;

    if ((use_time_management() && elapsed > time_maximum() - 10)
        || (Limits.movetime && elapsed >= Limits.movetime)
        || (Limits.nodes && Threads.pos[0]->nodes >= Limits.nodes))
        Threads.stop = 1;
}

// uci_print_pv() prints PV information according to the UCI protocol.
// UCI requires that all (if any) unsearched PV lines are sent with a
// previous search score.

static void uci_print_pv(Position* pos, Depth depth, Value alpha, Value beta) {
    TimePoint  elapsed        = time_elapsed() + 1;
    RootMoves* rm             = pos->rootMoves;
    int        pvIdx          = pos->pvIdx;
    uint64_t   nodes_searched = Threads.pos[0]->nodes;
    char       buf[16];

    int i = 0;

    bool updated = rm->move[i].score != -VALUE_INFINITE;

    Depth d = updated ? depth : max(1, depth - 1);
    Value v = updated ? rm->move[i].score : rm->move[i].previousScore;

    if (v == -VALUE_INFINITE)
        v = VALUE_ZERO;

    printf("info depth %d score %s", d, uci_value(buf, v));

    if (i == pvIdx)
        printf("%s", v >= beta ? " lowerbound" : v <= alpha ? " upperbound" : "");

    printf(" nodes %" PRIu64 " nps %" PRIu64, nodes_searched, nodes_searched * 1000 / elapsed);

    printf(" time %" PRIi64 " pv", elapsed);

    for (int idx = 0; idx < rm->move[i].pvSize; idx++)
        printf(" %s", uci_move(buf, rm->move[i].pv[idx]));
    printf("\n");

    fflush(stdout);
}

// start_thinking() wakes up the main thread to start a new search,
// then returns immediately.

void start_thinking(Position* root, bool ponderMode) {
    if (Threads.searching)
        thread_wait_until_sleeping(threads_main());

    prepare_for_search(root, ponderMode);
    thread_wake_up(threads_main(), THREAD_SEARCH);
}

void prepare_for_search(Position* root, bool ponderMode) {
    Threads.stop          = false;
    Threads.increaseDepth = true;
    Threads.ponder        = ponderMode;

    // Generate all legal moves.
    ExtMove    list[MAX_MOVES];
    ExtMove*   end   = generate_legal(root, list);
    RootMoves* moves = Threads.pos[0]->rootMoves;
    moves->size      = end - list;
    for (int i = 0; i < moves->size; i++)
        moves->move[i].pv[0] = list[i].move;

    Position* pos  = Threads.pos[0];
    pos->rootDepth = 0;
    pos->nodes     = 0;
    RootMoves* rm  = pos->rootMoves;
    rm->size       = end - list;
    for (int i = 0; i < rm->size; i++)
    {
        rm->move[i].pvSize        = 1;
        rm->move[i].pv[0]         = moves->move[i].pv[0];
        rm->move[i].score         = -VALUE_INFINITE;
        rm->move[i].previousScore = -VALUE_INFINITE;
        rm->move[i].bestMoveCount = 0;
        rm->move[i].tbRank        = moves->move[i].tbRank;
    }
    memcpy(pos, root, offsetof(Position, moveList));
    // Copy enough of the root State buffer.
    int n = max(7, root->st->pliesFromNull);
    for (int i = 0; i <= n; i++)
        memcpy(&pos->stack[i], &root->st[i - n], StateSize);
    pos->st                 = pos->stack + n;
    (pos->st - 1)->endMoves = pos->moveList;
    pos_set_check_info(pos);

    Threads.searching = true;
}
