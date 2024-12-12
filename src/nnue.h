#ifndef NNUE_H
#define NNUE_H

#include <stdalign.h>

#include "types.h"

#define INSIZE  768
#define L1SIZE  64
#define OUTSIZE 1

#define QA    192
#define QB     64
#define SCALE 400

typedef struct  Accumulator Accumulator;

struct Accumulator {
    alignas(64) int16_t values[2][L1SIZE];
};

void nnue_init();
Value nnue_evaluate(const Position *pos);

#endif //NNUE_H
