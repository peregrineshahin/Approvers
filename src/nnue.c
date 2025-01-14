#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../raw.nnue");

alignas(64) float in_weights[INSIZE * L1SIZE];
alignas(64) float l1_weights[L1SIZE * 2];

alignas(64) float in_biases[L1SIZE];
alignas(64) float l1_bias;

SMALL void nnue_init() {
    float* data = (float*) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
        in_weights[i] = *(data++);

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data++);

    for (int i = 0; i < L1SIZE * 2; i++)
        l1_weights[i] = *(data++);

    l1_bias = *(data);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static float screlu(const float x) {
    const float v = max(0.0, min(1.0, x));
    return v * v;
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    const float* stm  = acc->values[pos->sideToMove];
    const float* nstm = acc->values[!pos->sideToMove];

    float output = 0;
    for (int i = 0; i < L1SIZE; i++)
    {
        output += screlu(stm[i]) * l1_weights[i];
        output += screlu(nstm[i]) * l1_weights[i + L1SIZE];
    }

    return (Value) ((output + l1_bias) * SCALE);
}

static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
    memcpy(acc->values[side], in_biases, L1SIZE * sizeof(float));

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
