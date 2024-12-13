#ifndef EVALUATE_H
#define EVALUATE_H

#include "types.h"

enum {
    Tempo = 28
};

Value evaluate(Position* pos);
void hint_nnue(Position* pos);
#endif // EVALUATE_H
