#ifndef NNUE_H
#define NNUE_H

#include <stdalign.h>

#include "types.h"

#define INSIZE 768
#define L1SIZE 64

#define QA 101
#define QB 160
#define SCALE 340

typedef struct Accumulator Accumulator;

struct Accumulator {
    bool needs_refresh;
    alignas(64) int16_t values[2][L1SIZE];
};

void nnue_init();

void  nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq);
void  nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq);
Value nnue_evaluate(Position* pos);

#endif  //NNUE_H
