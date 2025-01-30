cmake_minimum_required(VERSION 3.30)
project(Approvers C)

set(CMAKE_C_STANDARD 11)

add_compile_options(-mavx2)

include_directories(src)

add_executable(Approvers
        src/benchmark.c
        src/bitboard.c
        src/bitboard.h
        src/evaluate.c
        src/evaluate.h
        src/main.c
        src/misc.c
        src/misc.h
        src/movegen.c
        src/movegen.h
        src/movepick.c
        src/movepick.h
        src/position.c
        src/position.h
        src/search.c
        src/search.h
        src/thread.c
        src/thread.h
        src/timeman.c
        src/timeman.h
        src/tt.c
        src/tt.h
        src/types.h
        src/uci.c
        src/uci.h
        src/incbin.h
        src/nnue.c
        src/nnue.h)
