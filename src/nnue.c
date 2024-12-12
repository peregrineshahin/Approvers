#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../default.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t l1_weights[L1SIZE * OUTSIZE * 2];

alignas(64) int16_t in_biases[L1SIZE];
alignas(64) int16_t l1_biases[OUTSIZE];

void nnue_init() {
    int16_t *data = (int16_t *) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
        in_weights[i] = *(data++);

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data++);

    for (int i = 0; i < L1SIZE * OUTSIZE * 2; i++)
        l1_weights[i] = *(data++);

    for (int i = 0; i < OUTSIZE; i++)
        l1_biases[i] = *(data++);
}

static int make_index(PieceType pt, Color c, Square sq, Color side) {
    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static Value activate(const int x) {
    const int v = max(0, min(QA, x));
    return v * v;
}

static Value output_transform(const Accumulator *acc, const Position *pos) {
    const int16_t *stm = acc->values[pos->sideToMove];
    const int16_t *nstm = acc->values[!pos->sideToMove];

    int output = 0;
    for (int i = 0; i < L1SIZE; i++) {
        output += activate(stm[i]) * (int) l1_weights[i];
        output += activate(nstm[i]) * (int) l1_weights[i + L1SIZE];
    }
    return (output / QA + l1_biases[0]) * SCALE / (QA * QB);
}

static void build_accumulator(Accumulator *acc, const Position *pos, Color side) {
    memcpy(acc->values[side], in_biases, sizeof(acc->values[side]));

    for (int c = WHITE; c <= BLACK; c++) {
        for (int pt = PAWN; pt <= KING; pt++) {
            Bitboard pieces = pieces_cp(c, pt);
            while (pieces) {
                const int idx = make_index(pt, c, pop_lsb(&pieces), side);
                for (int i = 0; i < L1SIZE; i++)
                    acc->values[side][i] += in_weights[idx * L1SIZE + i];
            }
        }
    }
}

Value nnue_evaluate(const Position *pos) {
    Accumulator acc = {0};
    build_accumulator(&acc, pos, WHITE);
    build_accumulator(&acc, pos, BLACK);

    return output_transform(&acc, pos);
}
