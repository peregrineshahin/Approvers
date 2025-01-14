#ifndef NNUE_H
#define NNUE_H

#include <stdalign.h>

#include "types.h"

#define FT_SIZE 768
#define L1_SIZE 64
#define L2_SIZE 8
#define L3_SIZE 16

#define SCALE 400

typedef struct Accumulator Accumulator;

struct Accumulator {
    bool needs_refresh;
    alignas(64) float values[2][L1_SIZE];
};

void nnue_init();

void  nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq);
void  nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq);
Value nnue_evaluate(Position* pos);

#endif  //NNUE_H
