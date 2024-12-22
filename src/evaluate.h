#ifndef EVALUATE_H
#define EVALUATE_H

#include "types.h"

enum {
    Tempo = 28
};

Value evaluate(Position* pos);
int simple_eval(Position* pos, Color c);
#endif  // EVALUATE_H
