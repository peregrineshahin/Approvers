cmake_minimum_required(VERSION 3.30)
project(Peregrine C)

set(CMAKE_C_STANDARD 11)

# Enable this to use the Kaggle environment
# add_definitions(-DKAGGLE)

include_directories(src)

add_executable(Peregrine
        src/benchmark.c
        src/bitbase.c
        src/bitboard.c
        src/bitboard.h
        src/config.h
        src/evaluate.c
        src/evaluate.h
        src/main.c
        src/misc.c
        src/misc.h
        src/movegen.c
        src/movegen.h
        src/movepick.c
        src/movepick.h
        src/pawns.c
        src/pawns.h
        src/position.c
        src/position.h
        src/psqt.c
        src/search.c
        src/search.h
        src/settings.c
        src/settings.h
        src/thread.c
        src/thread.h
        src/timeman.c
        src/timeman.h
        src/tt.c
        src/tt.h
        src/types.h
        src/uci.c
        src/uci.h
        src/ucioption.c
        src/incbin.h
        src/nnue.c
        src/nnue.h)
