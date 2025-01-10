#ifndef NNUE_H
#define NNUE_H

#include <stdalign.h>

#include "types.h"

#define INSIZE 768
#define L1SIZE 64

#define QA 192
#define QB 64
#define SCALE 400

typedef struct Accumulator Accumulator;

enum {
    ACC_EMPTY,
    ACC_FORCED,
    ACC_COMPUTED,
};

struct Accumulator {
    uint8_t state;
    alignas(64) int16_t values[2][L1SIZE];
};

typedef struct DirtyPiece DirtyPiece;

struct DirtyPiece {
    int    len;
    Piece  piece[3];
    Square from[3];
    Square to[3];
};

void nnue_init();

Value nnue_evaluate(Position* pos);

#endif  //NNUE_H
