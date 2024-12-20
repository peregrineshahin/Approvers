# Stockfish, a UCI chess playing engine derived from Glaurung 2.1
# Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
# Copyright (C) 2008-2015 Marco Costalba, Joona Kiiski, Tord Romstad
# Copyright (C) 2015-2016 Marco Costalba, Joona Kiiski, Gary Linscott, Tord Romstad
#
# Stockfish is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stockfish is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
CC      = clang
CFLAGS  = -fcommon -O3 -ffunction-sections -fdata-sections -mavx2

ifeq ($(comp),$(filter $(comp),gcc clang icc))
    ifeq ($(KERNEL),Darwin)
        CFLAGS += -mdynamic-no-pic
    endif
endif
LDFLAGS = -Wl,--gc-sections
OBJ     = uci.o settings.o nnue.o benchmark.o bitboard.o evaluate.o    \
        main.o misc.o movegen.o movepick.o position.o psqt.o search.o    \
        thread.o timeman.o tt.o
EXE     = cfish

all: build

build: clean $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS)

kaggle: CFLAGS += -DKAGGLE
kaggle: clean $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXE) $(LDFLAGS)

bench: build
	./$(EXE) bench

clean:
	@rm -f $(EXE) $(OBJ)

format:
	find . -name "*.c" -o -name "*.h" | xargs $(CLANG-FORMAT) -i -style=file