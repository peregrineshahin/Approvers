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

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "search.h"
#include "settings.h"
#include "thread.h"
#include "uci.h"

#ifndef KAGGLE
extern void benchmark();

extern Parameter parameters[255];
extern int       parameters_count;
#endif

// FEN string of the initial position, normal chess
char StartFEN[] = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";

// position() is called when the engine receives the "position" UCI
// command. The function sets up the position described in the given FEN
// string ("fen") or the starting position ("startpos") and then makes
// the moves given in the following move list ("moves").

SMALL void position(Position* pos, char* str) {
    char  fen[128];
    char* moves;

    moves = strstr(str, "moves");
    if (moves)
    {
        if (moves > str)
            moves[-1] = 0;
        moves += 5;
    }

    if (strncmp(str, "fen", 3) == 0)
    {
        strncpy(fen, str + 4, 127);
        fen[127] = 0;
    }
#ifndef DKAGGLE
    else if (strncmp(str, "startpos", 8) == 0)
        strcpy(fen, StartFEN);
#endif
    else
        return;

    pos->st = pos->stack + 100;  // Start of circular buffer of 100 slots.
    pos_set(pos, fen);

    // Parse move list (if any).
    if (moves)
    {
        int ply = 0;

        for (moves = strtok(moves, " \t"); moves; moves = strtok(NULL, " \t"))
        {
            Move m = uci_to_move(pos, moves);
            if (!m)
                break;
            do_move(pos, m, gives_check(pos, pos->st, m));
            pos->gamePly++;
            // Roll over if we reach 100 plies.
            if (++ply == 100)
            {
                memcpy(pos->st - 100, pos->st, StateSize);
                pos->st -= 100;
                pos_set_check_info(pos);
                ply -= 100;
            }
        }

        // Make sure that is_draw() never tries to look back more than 99 ply.
        // This is enough, since 100 ply history means draw by 50-move rule.
        if (pos->st->pliesFromNull > 99)
            pos->st->pliesFromNull = 99;

        // Now move some of the game history at the end of the circular buffer
        // in front of that buffer.
        int k = (pos->st - (pos->stack + 100)) - max(7, pos->st->pliesFromNull);
        for (; k < 0; k++)
            memcpy(pos->stack + 100 + k, pos->stack + 200 + k, StateSize);
    }

    pos->rootKeyFlip        = pos->st->key;
    (pos->st - 1)->endMoves = pos->moveList;

    // Clear history position keys that have not yet repeated. This ensures
    // that is_draw() does not flag as a draw the first repetition of a
    // position coming before the root position. In addition, we set
    // pos->hasRepeated to indicate whether a position has repeated since
    // the last irreversible move.
    for (int k = 0; k <= pos->st->pliesFromNull; k++)
    {
        int l;
        for (l = k + 4; l <= pos->st->pliesFromNull; l += 2)
            if ((pos->st - k)->key == (pos->st - l)->key)
                break;
        if (l <= pos->st->pliesFromNull)
            pos->hasRepeated = true;
        else
            (pos->st - k)->key = 0;
    }
    pos->rootKeyFlip ^= pos->st->key;
    pos->st->key ^= pos->rootKeyFlip;
}


// setoption() is called when the engine receives the "setoption" UCI
// command. The function updates the UCI option ("name") to the given
// value ("value").

void setoption(char* str) {
    char* name  = strstr(str, "name") + 5;
    char* value = strstr(name, "value");
    if (value)
    {
        char* p = value - 1;
        while (isblank(*p))
            p--;
        p[1] = 0;
        value += 5;
        while (isblank(*value))
            value++;
    }

#ifndef KAGGLE
    for (int i = 0; i < parameters_count; i++)
    {
        if (strcmp(parameters[i].name, name) == 0)
        {
            *parameters[i].value = atoi(value);
            return;
        }
    }
#endif

    if (strcmp("Hash", name) == 0)
    {
        delayedSettings.ttSize = atoi(value);
        return;
    }

    if (strcmp("Pondering", name) == 0)
    {
        Thread.testPonder = strcmp("true", value) == 0;
        return;
    }

    fprintf(stderr, "No such option: %s\n", name);
}


// go() is called when engine receives the "go" UCI command. The function sets
// the thinking time and other parameters from the input string, then starts
// the search.

static void go(Position* pos, char* str) {
    bool ponderMode = false;

    process_delayed_settings();

    Limits           = (struct LimitsType) {0};
    Limits.startTime = now();  // As early as possible!

    for (char* token = strtok(str, " \t"); token; token = strtok(NULL, " \t"))
    {
        if (strcmp(token, "wtime") == 0)
            Limits.time[WHITE] = atoi(strtok(NULL, " \t"));
        else if (strcmp(token, "btime") == 0)
            Limits.time[BLACK] = atoi(strtok(NULL, " \t"));
        else if (strcmp(token, "winc") == 0)
            Limits.inc[WHITE] = atoi(strtok(NULL, " \t"));
        else if (strcmp(token, "binc") == 0)
            Limits.inc[BLACK] = atoi(strtok(NULL, " \t"));
#ifndef DKAGGLE
        else if (strcmp(token, "depth") == 0)
            Limits.depth = atoi(strtok(NULL, " \t"));
        else if (strcmp(token, "nodes") == 0)
            Limits.nodes = atoi(strtok(NULL, " \t"));
#endif
    }

    start_thinking(pos, ponderMode);
}


// uci_loop() waits for a command from stdin, parses it and calls the
// appropriate function. Also intercepts EOF from stdin to ensure
// gracefully exiting if the GUI dies unexpectedly. When called with some
// command line arguments, e.g. to run 'bench', once the command is
// executed the function returns immediately. In addition to the UCI ones,
// also some additional debug commands are supported.

SMALL void uci_loop(int argc, char** argv) {
#ifndef KAGGLE
    if (argc == 2 && strcmp(argv[1], "bench") == 0)
    {
        benchmark();
        return;
    }
#endif

    Position pos;
    char     fen[strlen(StartFEN) + 1];

    // Allocate 215 Stack slots.
    // Slots 100-200 form a circular buffer to be filled with game moves.
    // Slots 0-99 make room for prepending the part of game history relevant
    // for repetition detection.
    // Slots 201-214 may be used by TB root probing.
    pos.stackAllocation = malloc(63 + 215 * sizeof(Stack));
    pos.stack           = (Stack*) (((uintptr_t) pos.stackAllocation + 0x3f) & ~0x3f);
    pos.moveList        = malloc(1000 * sizeof(ExtMove));
    pos.st              = pos.stack + 100;
    pos.st[-1].endMoves = pos.moveList;

    char cmd[4096] = {0};

    delayedSettings.ttSize = 1;

    strcpy(fen, StartFEN);
    pos_set(&pos, fen);
    pos.rootKeyFlip = pos.st->key;

    while (get_input(cmd))
    {
        char* token = cmd;
        while (isblank(*token))
            token++;

        char* str = token;
        while (*str && !isblank(*str))
            str++;

        if (*str)
        {
            *str++ = 0;
            while (isblank(*str))
                str++;
        }

        // At this point we have received a command, so stop pondering
        Thread.ponder = false;
        Thread.stop   = true;

        if (strcmp(token, "quit") == 0)
            break;

        if (strcmp(token, "uci") == 0)
        {
#ifndef KAGGLE
            printf("id name\n");
            printf("option name Threads type spin default 1 min 1 max 2\n");
            printf("option name Hash type spin default 1 min 1 max 16\n");
            printf("option name Pondering type string\n");

            for (int i = 0; i < parameters_count; i++)
                printf("option name %s type string\n", parameters[i].name);

            printf("uciok\n");
            fflush(stdout);
#endif
        }
        else if (strcmp(token, "ucinewgame") == 0)
        {
            process_delayed_settings();
            search_clear();
        }
        else if (strcmp(token, "isready") == 0)
        {
            process_delayed_settings();
            printf("readyok\n");
            fflush(stdout);
        }
        else if (strcmp(token, "go") == 0)
            go(&pos, str);
        else if (strcmp(token, "position") == 0)
            position(&pos, str);
        else if (strcmp(token, "setoption") == 0)
            setoption(str);
    }

    free(pos.stackAllocation);
    free(pos.moveList);
}


// uci_value() converts a Value to a string suitable for use with the UCI
// protocol specification:
//
// cp <x>    The score from the engine's point of view in centipawns.
// mate <y>  Mate in y moves, not plies. If the engine is getting mated
//           use negative values for y.

char* uci_value(char* str, Value v) {
    if (abs(v) < VALUE_MATE_IN_MAX_PLY)
        sprintf(str, "cp %d", v);
    else
        sprintf(str, "mate %d", (v > 0 ? VALUE_MATE - v + 1 : -VALUE_MATE - v) / 2);

    return str;
}


// uci_square() converts a Square to a string in algebraic notation
// (g1, a7, etc.)

char* uci_square(char* str, Square s) {
    str[0] = 'a' + file_of(s);
    str[1] = '1' + rank_of(s);
    str[2] = 0;

    return str;
}


// uci_move() converts a Move to a string in coordinate notation (g1f3, a7a8q).
// Internally all castling moves are always encoded as 'king captures rook'.

char* uci_move(char* str, Move m) {
    char   buf1[8], buf2[8];
    Square from = from_sq(m);
    Square to   = to_sq(m);

    if (type_of_m(m) == CASTLING)
        to = make_square(to > from ? FILE_G : FILE_C, rank_of(from));

    strcat(strcpy(str, uci_square(buf1, from)), uci_square(buf2, to));

    if (type_of_m(m) == PROMOTION)
    {
        str[strlen(str) + 1] = 0;
        str[strlen(str)]     = " pnbrqk"[promotion_type(m)];
    }

    return str;
}


// uci_to_move() converts a string representing a move in coordinate
// notation (g1f3, a7a8q) to the corresponding legal Move, if any.

Move uci_to_move(const Position* pos, char* str) {
    if (strlen(str) == 5)
        str[4] = tolower(str[4]);

    ExtMove  list[MAX_MOVES];
    ExtMove* last = generate_legal(pos, list);

    char buf[16];

    for (ExtMove* m = list; m < last; m++)
        if (strcmp(str, uci_move(buf, m->move)) == 0)
            return m->move;

    return 0;
}

int get_input(char* str) {
    if (fgets(str, 4096, stdin) == NULL)
        return 0;

    char* ptr = strchr(str, '\n');
    if (ptr != NULL)
        *ptr = '\0';

    ptr = strchr(str, '\r');
    if (ptr != NULL)
        *ptr = '\0';

    return 1;
}
