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

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "evaluate.h"
#include "misc.h"
#include "position.h"
#include "search.h"
#include "settings.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"

static char *Defaults[] = {
  "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1"
};

// benchmark() runs a simple benchmark by letting Stockfish analyze a set
// of positions for a given limit each. There are six optional parameters:
// - Transposition table size. Default is 16 MB.
// - Number of search threads to use. Default is 1 thread.
// - Limit value for each search. Default is (depth) 13.
// - File name with the positions to search in FEN format. The default
//   positions are listed above.
// - Type of the limit value: depth (default), time (in msecs), nodes.
// - Evaluation: classical, nnue (hybrid), pure (NNUE only), mixed (default).

void benchmark(Position *current, char *str)
{
  char *token;
  char **fens;
  int numFens;

  Limits = (struct LimitsType){ 0 };

  int ttSize      = (token = strtok(str , " ")) ? atoi(token)  : 16;
  int threads     = (token = strtok(NULL, " ")) ? atoi(token)  : 1;
  int64_t limit   = (token = strtok(NULL, " ")) ? atoll(token) : 13;
  char *fenFile   = (token = strtok(NULL, " ")) ? token        : "default";
  char *limitType = (token = strtok(NULL, " ")) ? token        : "depth";

  delayedSettings.ttSize = ttSize;
  delayedSettings.numThreads = threads;
  process_delayed_settings();
  search_clear();

  if (strcmp(limitType, "time") == 0)
    Limits.movetime = limit; // movetime is in millisecs
  else if (strcmp(limitType, "nodes") == 0)
    Limits.nodes = limit;
  else if (strcmp(limitType, "mate") == 0)
    Limits.mate = limit;
  else
    Limits.depth = limit;

  if (strcasecmp(fenFile, "default") == 0) {
    fens = Defaults;
    numFens = sizeof(Defaults) / sizeof(char *);
  }
  else if (strcasecmp(fenFile, "current") == 0) {
    fens = malloc(sizeof(*fens));
    fens[0] = malloc(128);
    pos_fen(current, fens[0]);
    numFens = 1;
  }
  else {
    int maxFens = 100;
    numFens = 0;
    FILE *F = fopen(fenFile, "r");
    if (!F) {
      return;
    }
    fens = malloc(maxFens * sizeof(*fens));
    fens[0] = NULL;
    size_t length = 0;
    while (getline(&fens[numFens], &length, F) > 0) {
      numFens++;
      if (numFens == maxFens) {
        maxFens += 100;
        fens = realloc(fens, maxFens * sizeof(*fens));
      }
      fens[numFens] = NULL;
      length = 0;
    }
    fclose(F);
  }

  uint64_t nodes = 0;
  Position pos;
  memset(&pos, 0, sizeof(pos));
  pos.stackAllocation = malloc(63 + 217 * sizeof(*pos.stack));
  pos.stack = (Stack *)(((uintptr_t)pos.stackAllocation + 0x3f) & ~0x3f);
  pos.st = pos.stack + 7;
  pos.moveList = malloc(10000 * sizeof(*pos.moveList));
  TimePoint elapsed = now();

  int numOpts = 0;
  for (int i = 0; i < numFens; i++)
    if (strncmp(fens[i], "setoption ", 9) == 0)
      numOpts++;

  for (int i = 0, j = 0; i < numFens; i++) {
    char buf[128];

    if (strncmp(fens[i], "setoption ", 9) == 0) {
      strncpy(buf, fens[i] + 10, 127 - 10);
      buf[127] = 0;
      setoption(buf);
      continue;
    }

    strcpy(buf, "fen ");
    strncat(buf, fens[i], 127 - 4);
    buf[127] = 0;

    position(&pos, buf);

    fprintf(stderr, "\nPosition: %d/%d\n", ++j, numFens - numOpts);

      Limits.startTime = now();
      start_thinking(&pos, false);
      thread_wait_until_sleeping(threads_main());
      nodes += threads_nodes_searched();
  }

  elapsed = now() - elapsed + 1; // Ensure positivity to avoid a 'divide by zero'

  fprintf(stderr, "\n==========================="
                  "\nTotal time (ms) : %" PRIu64
                  "\nNodes searched  : %" PRIu64
                  "\nNodes/second    : %" PRIu64 "\n",
                  elapsed, nodes, 1000 * nodes / elapsed);

  if (fens != Defaults) {
    for (int i = 0; i < numFens; i++)
      free(fens[i]);
    free(fens);
  }
  free(pos.stackAllocation);
  free(pos.moveList);
}
