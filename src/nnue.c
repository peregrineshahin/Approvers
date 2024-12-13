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
    int8_t *data8 = (int8_t *) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
        in_weights[i] = *(data8++);

    int16_t *data16 = (int16_t *) data8;

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data16++);

    for (int i = 0; i < L1SIZE * OUTSIZE * 2; i++)
        l1_weights[i] = *(data16++);

    for (int i = 0; i < OUTSIZE; i++)
        l1_biases[i] = *(data16++);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (file_of(ksq) > 3)
        sq ^= 7;

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

    Square ksq = square_of(side, KING);
    for (int c = WHITE; c <= BLACK; c++) {
        for (int pt = PAWN; pt <= KING; pt++) {
            Bitboard pieces = pieces_cp(c, pt);
            while (pieces) {
                const int idx = make_index(pt, c, pop_lsb(&pieces), ksq, side);
                for (int i = 0; i < L1SIZE; i++)
                    acc->values[side][i] += in_weights[idx * L1SIZE + i];
            }
        }
    }
}

void nnue_add_piece(Accumulator *acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++) {
        acc->values[WHITE][i] += in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] += in_weights[black * L1SIZE + i];
    }
}

void nnue_remove_piece(Accumulator *acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++) {
        acc->values[WHITE][i] -= in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] -= in_weights[black * L1SIZE + i];
    }
}

Value nnue_evaluate(Position *pos) {
    Accumulator *acc = &pos->st->accumulator;

    if (acc->needs_refresh) {
        build_accumulator(acc, pos, WHITE);
        build_accumulator(acc, pos, BLACK);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
