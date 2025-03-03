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
#include "thread.h"
#include "tt.h"
#include "uci.h"

#ifndef KAGGLE
extern void benchmark();

extern Parameter parameters[255];
extern int       parameters_count;

extern alignas(64) int16_t in_biases[L1SIZE];
#endif

// FEN string of the initial position, normal chess
char StartFEN[] = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";

char KagglePosition[4096] = {0};

// position() is called when the engine receives the "position" UCI
// command. The function sets up the position described in the given FEN
// string ("fen") or the starting position ("startpos") and then makes
// the moves given in the following move list ("moves").
SMALL void position(Position* pos, char* str) {
    // Start of circular buffer of 100 slots.
    pos->st = pos->stack + 100;

    if (strncmp(str, "fen", 3) == 0)
        pos_set(pos, str + 4);
#ifndef KAGGLE
    else if (strncmp(str, "startpos", 8) == 0)
        pos_set(pos, StartFEN);
#endif

    // Parse move list (if any).
    char* moves = strstr(str, "moves");
    if (moves)
    {
        if (moves > str)
            moves[-1] = 0;
        moves += 5;

        int ply = 0;
#pragma clang loop unroll(disable)
        for (moves = strtok(moves, " \t"); moves; moves = strtok(NULL, " \t"))
        {
            Move move = uci_to_move(pos, moves);
            do_move(pos, move, gives_check(pos, pos->st, move));
            pos->accumulator--;
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
#pragma clang loop unroll(disable)
        for (; k < 0; k++)
            memcpy(pos->stack + 100 + k, pos->stack + 200 + k, StateSize);
    }
}


// Called when the engine receives the "setoption" UCI. The function
// updates the UCI option ("name") to the given value ("value").
void setoption(char* str) {
#ifndef KAGGLE
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

    for (int i = 0; i < parameters_count; i++)
    {
        if (strcmp(parameters[i].name, name) == 0)
        {
            *parameters[i].value = atoi(value);
            return;
        }
    }

    if (strcmp("Pondering", name) == 0)
    {
        Thread.testPonder = strcmp("true", value) == 0;
        return;
    }

    if (strstr(name, "inb_v"))
    {
        int i        = atoi(name + 5);
        in_biases[i] = atoi(value);
        return;
    }

    fprintf(stderr, "No such option: %s\n", name);
#endif
}


// Called when engine receives the "go" UCI command. The function sets
// the thinking time and other parameters from the input string, then
// starts the search.
static void go(char* str) {
    Limits           = (struct LimitsType) {0};
    Limits.startTime = now();  // As early as possible!

#pragma clang loop unroll(disable)
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
#ifndef KAGGLE
        else if (strcmp(token, "depth") == 0)
            Limits.depth = atoi(strtok(NULL, " \t"));
        else if (strcmp(token, "nodes") == 0)
            Limits.nodes = atoi(strtok(NULL, " \t"));
#endif
    }

    start_thinking();
}


// Waits for a command from stdin, parses it and calls the
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

    char cmd[4096] = {0};
    char fen[strlen(StartFEN) + 1];

    Position* pos = Thread.pos;

    strcpy(fen, StartFEN);
    pos_set(pos, fen);

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

        if (strcmp(token, "go") == 0)
            go(str);
        if (strcmp(token, "kpos") == 0)
        {
            char* lastMove = strstr(str, "moves") + 6;
            if (!*lastMove || KagglePosition[0] == 0)
            {
                lastMove[-1] = 0;
                strcpy(KagglePosition, str);
            }
            else
            {
                strcat(KagglePosition, " ");
                strcat(KagglePosition, lastMove);
            }

            char* buffer = strdup(KagglePosition);
            position(pos, buffer);
            free(buffer);
        }
#ifndef KAGGLE
        else if (strcmp(token, "position") == 0)
            position(pos, str);
        else if (strcmp(token, "uci") == 0)
        {
            printf("id name\n");
            printf("option name Threads type spin default 1 min 1 max 2\n");
            printf("option name Hash type spin default 1 min 1 max 16\n");
            printf("option name Pondering type string\n");

            for (int i = 0; i < parameters_count; i++)
                printf("option name %s type string\n", parameters[i].name);

            for (int i = 0; i < L1SIZE; i++)
                printf("option name inb_v%d type string\n", i);

            printf("uciok\n");
            fflush(stdout);
        }
        else if (strcmp(token, "params") == 0)
        {
            for (int i = 0; i < parameters_count; i++)
            {
                const int v   = *parameters[i].value;
                const int min = (v > 0) ? 1 : v * 2;
                const int max = (v > 0) ? v * 2 : -1;

                printf("%s, int, %d, %d, %d, %.2f, 0.002\n", parameters[i].name, v, min, max,
                       (double) (max - min) / 20.0);
            }
        }
        else if (strcmp(token, "nnparams") == 0)
        {
            for (int i = 0; i < L1SIZE; i++)
                printf("inb_v%d, int, %d, -128, 127, %.3f, 0.002\n", i, in_biases[i],
                       abs(in_biases[i]) / 20.0);
        }
        else if (strcmp(token, "ucinewgame") == 0)
        {
            search_clear();
        }
        else if (strcmp(token, "isready") == 0)
        {
            printf("readyok\n");
            fflush(stdout);
        }
        else if (strcmp(token, "setoption") == 0)
            setoption(str);
        else if (strcmp(token, "quit") == 0)
            break;
#endif
    }
}


// Converts a Value to a string suitable for use with the UCI protocol
// specification:
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


// Converts a Square to a string in algebraic notation (g1, a7, etc.)
char* uci_square(char* str, Square s) {
    str[0] = 'a' + file_of(s);
    str[1] = '1' + rank_of(s);
    str[2] = 0;

    return str;
}


// Converts a Move to a string in coordinate notation (g1f3, a7a8q).
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


// Converts a string representing a move in coordinate notation (g1f3, a7a8q)
// to the corresponding legal Move, if any.
Move uci_to_move(const Position* pos, char* str) {
    ExtMove  list[MAX_MOVES];
    ExtMove* last = generate_pseudo_legal(pos, list);

    char buf[16];

#pragma clang loop unroll(disable)
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
