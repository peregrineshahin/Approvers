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


#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#ifndef _WIN32
#include <sys/mman.h>
#endif

#include "evaluate.h"
#include "misc.h"
#include "search.h"
#include "settings.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"

// 'On change' actions, triggered by an option's value change
static void on_clear_hash(Option *opt)
{
  (void)opt;

  if (settings.ttSize)
    search_clear();
}

static void on_hash_size(Option *opt)
{
  delayedSettings.ttSize = opt->value;
}

static void on_threads(Option *opt)
{
  delayedSettings.numThreads = opt->value;
}

static void on_large_pages(Option *opt)
{
  delayedSettings.largePages = opt->value;
}

#ifdef IS_64BIT
#define MAXHASHMB 33554432
#else
#define MAXHASHMB 2048
#endif

static Option optionsMap[] = {
  { "Contempt", OPT_TYPE_SPIN, 24, -100, 100, NULL, NULL, 0, NULL },
  { "Threads", OPT_TYPE_SPIN, 1, 1, MAX_THREADS, NULL, on_threads, 0, NULL },
  { "Hash", OPT_TYPE_SPIN, 1, 1, MAXHASHMB, NULL, on_hash_size, 0, NULL },
  { "Clear Hash", OPT_TYPE_BUTTON, 0, 0, 0, NULL, on_clear_hash, 0, NULL },
  { "Ponder", OPT_TYPE_CHECK, 0, 0, 0, NULL, NULL, 0, NULL },
  { "Move Overhead", OPT_TYPE_SPIN, 10, 0, 5000, NULL, NULL, 0, NULL },
  { "LargePages", OPT_TYPE_CHECK, 1, 0, 0, NULL, on_large_pages, 0, NULL },
  { 0 }
};

// options_init() initializes the UCI options to their hard-coded default
// values.

void options_init()
{
#ifdef _WIN32
  // Disable the LargePages option if the machine does not support it.
  if (!large_pages_supported())
    optionsMap[OPT_LARGE_PAGES].type = OPT_TYPE_DISABLED;
#endif
#if defined(__linux__) && !defined(MADV_HUGEPAGE)
  optionsMap[OPT_LARGE_PAGES].type = OPT_TYPE_DISABLED;
#endif
  for (Option *opt = optionsMap; opt->name != NULL; opt++) {
    if (opt->type == OPT_TYPE_DISABLED)
      continue;
    switch (opt->type) {
    case OPT_TYPE_CHECK:
    case OPT_TYPE_SPIN:
      opt->value = opt->def;
    case OPT_TYPE_BUTTON:
      break;
    }
    if (opt->onChange)
      opt->onChange(opt);
  }
}

static const char *optTypeStr[] = {
  "check", "spin", "button"
};

// print_options() prints all options in the format required by the
// UCI protocol.

void print_options(void) {
  for (Option *opt = optionsMap; opt->name != NULL; opt++) {
    if (opt->type == OPT_TYPE_DISABLED)
      continue;
    printf("option name %s type %s", opt->name, optTypeStr[opt->type]);
    switch (opt->type) {
    case OPT_TYPE_CHECK:
      printf(" default %s", opt->def ? "true" : "false");
      break;
    case OPT_TYPE_SPIN:
      printf(" default %d min %d max %d", opt->def, opt->minVal, opt->maxVal);
    case OPT_TYPE_BUTTON:
      break;
    }
    printf("\n");
  }
  fflush(stdout);
}

int option_value(int optIdx)
{
  return optionsMap[optIdx].value;
}

const char *option_string_value(int optIdx)
{
  return optionsMap[optIdx].valString;
}

const char *option_default_string_value(int optIdx)
{
  return optionsMap[optIdx].defString;
}

void option_set_value(int optIdx, int value)
{
  Option *opt = &optionsMap[optIdx];

  opt->value = value;
  if (opt->onChange)
    opt->onChange(opt);
}

bool option_set_by_name(char *name, char *value)
{
  for (Option *opt = optionsMap; opt->name != NULL; opt++) {
    if (opt->type == OPT_TYPE_DISABLED)
      continue;
    if (strcasecmp(opt->name, name) == 0) {
      int val;
      switch (opt->type) {
      case OPT_TYPE_CHECK:
        if (strcmp(value, "true") == 0)
          opt->value = 1;
        else if (strcmp(value, "false") == 0)
          opt->value = 0;
        else
          return true;
        break;
      case OPT_TYPE_SPIN:
        val = atoi(value);
        if (val < opt->minVal || val > opt->maxVal)
          return true;
        opt->value = val;
      case OPT_TYPE_BUTTON:
        break;
      }
      if (opt->onChange)
        opt->onChange(opt);
      return true;
    }
  }

  return false;
}
