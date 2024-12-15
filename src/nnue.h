#ifndef NNUE_H
#define NNUE_H

#include <stdalign.h>

#include "types.h"

#define INSIZE 768
#define L1SIZE 64
#define OUTSIZE 1

#define QA 192
#define QB 64
#define SCALE 400

enum {
    ACC_EMPTY,
    ACC_COMPUTED,
    ACC_INIT
};

typedef struct Accumulator Accumulator;

struct Accumulator {
    int state[2];
    alignas(64) int16_t values[2][L1SIZE];
};

typedef struct DirtyPiece DirtyPiece;

struct DirtyPiece {
    int    len;
    Piece  piece[3];
    Square from[3];
    Square to[3];
};

typedef struct IndexList IndexList;

struct IndexList {
    unsigned size;
    unsigned values[32];
};

void nnue_init();

Value nnue_evaluate(Position* pos);

#endif  //NNUE_H
