#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../default.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t l1_weights[L1SIZE * 2];

alignas(64) int16_t in_biases[L1SIZE];
alignas(64) int16_t l1_bias;

SMALL void nnue_init() {
    int8_t* data = (int8_t*) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
    {
        int x = i / L1SIZE;
        if (!(x < 8 || (56 <= x && x < 64) || (384 <= x && x < 392) || (440 <= x && x < 448)
              || (320 <= x && x < 384 && (x - 320) % 8 > 3)))
            in_weights[i] = *(data++);
    }

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data++);

    for (int i = 0; i < L1SIZE * 2; i++)
        l1_weights[i] = *(data++);

    l1_bias = *((int16_t*) data);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static Value screlu(const int x) {
    const int v = max(0, min(QA, x));
    return v * v;
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    const int16_t* stm  = acc->values[pos->sideToMove];
    const int16_t* nstm = acc->values[!pos->sideToMove];

    int output = 0;
    for (int i = 0; i < L1SIZE; i++)
    {
        output += screlu(stm[i]) * (int) l1_weights[i];
        output += screlu(nstm[i]) * (int) l1_weights[i + L1SIZE];
    }

    return (output / QA + l1_bias) * SCALE / (QA * QB);
}

static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
    memcpy(acc->values[side], in_biases, L1SIZE * sizeof(int16_t));

    const Square ksq = square_of(side, KING);
    for (Bitboard pieces = pieces(); pieces;)
    {
        const Square sq = pop_lsb(&pieces);
        const Piece  pc = piece_on(sq);

        const int index = make_index(type_of_p(pc), color_of(pc), sq, ksq, side) * L1SIZE;

        for (int i = 0; i < L1SIZE; i++)
            acc->values[side][i] += in_weights[index + i];
    }
}

void nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1SIZE;

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] += in_weights[white + i];
        acc->values[BLACK][i] += in_weights[black + i];
    }
}

void nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1SIZE;

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] -= in_weights[white + i];
        acc->values[BLACK][i] -= in_weights[black + i];
    }
}

Value nnue_evaluate(Position* pos) {
    Accumulator* acc = &pos->st->accumulator;

    if (acc->needs_refresh)
    {
        build_accumulator(acc, pos, WHITE);
        build_accumulator(acc, pos, BLACK);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
