#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../morelayers.nnue");

alignas(64) float ft_weights[FT_SIZE * L1_SIZE];
alignas(64) float ft_biases[L1_SIZE];

alignas(64) float l1_weights[L1_SIZE][L2_SIZE];
alignas(64) float l1_biases[L2_SIZE];

alignas(64) float l2_weights[L2_SIZE][L3_SIZE];
alignas(64) float l2_biases[L3_SIZE];

alignas(64) float l3_weights[L3_SIZE];
alignas(64) float l3_biases[1];

SMALL void nnue_init() {
    float* data = (float*) gNetworkData;

    for (int i = 0; i < FT_SIZE * L1_SIZE; i++)
        ft_weights[i] = *(data++);

    for (int i = 0; i < L1_SIZE; i++)
        ft_biases[i] = *(data++);

    for (int i = 0; i < L1_SIZE; i++)
        for (int j = 0; j < L2_SIZE; j++)
            l1_weights[i][j] = *(data++);

    for (int i = 0; i < L2_SIZE; i++)
        l1_biases[i] = *(data++);

    for (int i = 0; i < L2_SIZE; i++)
        for (int j = 0; j < L3_SIZE; j++)
            l2_weights[i][j] = *(data++);

    for (int i = 0; i < L3_SIZE; i++)
        l2_biases[i] = *(data++);

    for (int i = 0; i < L3_SIZE; i++)
        l3_weights[i] = *(data++);

    l3_biases[0] = *data;
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static float crelu(const float x) { return clamp(x, 0.0, 1.0); }

static Value output_transform(const Accumulator* acc, const Position* pos) {
    float ftOutput[L1_SIZE] = {0};
    float l1Output[L2_SIZE] = {0};
    float l2Output[L3_SIZE] = {0};

    // accumulators -> l1 (pairwise CReLU)
    for (int flip = 0; flip <= 1; flip++)
    {
        const float* input = acc->values[pos->sideToMove ^ flip];
        for (int i = 0; i < L1_SIZE / 2; ++i)
        {
            const float left  = crelu(input[i]);
            const float right = crelu(input[i + L1_SIZE / 2]);

            ftOutput[flip * L1_SIZE / 2 + i] = left * right;
        }
    }

    // l1 -> l2 (CReLU)
    for (int i = 0; i < L2_SIZE; i++)
    {
        float sum = l1_biases[i];
        for (int j = 0; j < L1_SIZE; j++)
            sum += ftOutput[j] * l1_weights[j][i];

        l1Output[i] = crelu(sum);
    }

    // l2 -> l3 (CReLU)
    for (int i = 0; i < L3_SIZE; i++)
    {
        float sum = l2_biases[i];
        for (int j = 0; j < L2_SIZE; j++)
            sum += l1Output[j] * l2_weights[j][i];

        l2Output[i] = crelu(sum);
    }

    // l3 -> output (linear)
    float sum = 0;
    for (int i = 0; i < L3_SIZE; ++i)
        sum += l2Output[i] * l3_weights[i];

    const float output = sum + l3_biases[0];
    return (Value) (output * SCALE);
}

static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
    memcpy(acc->values[side], ft_biases, L1_SIZE * sizeof(float));

    const Square ksq = square_of(side, KING);
    for (Bitboard pieces = pieces(); pieces;)
    {
        const Square sq = pop_lsb(&pieces);
        const Piece  pc = piece_on(sq);

        const int index = make_index(type_of_p(pc), color_of(pc), sq, ksq, side) * L1_SIZE;

        for (int i = 0; i < L1_SIZE; i++)
            acc->values[side][i] += ft_weights[index + i];
    }
}

void nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1_SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1_SIZE;

    for (int i = 0; i < L1_SIZE; i++)
    {
        acc->values[WHITE][i] += ft_weights[white + i];
        acc->values[BLACK][i] += ft_weights[black + i];
    }
}

void nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1_SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1_SIZE;

    for (int i = 0; i < L1_SIZE; i++)
    {
        acc->values[WHITE][i] -= ft_weights[white + i];
        acc->values[BLACK][i] -= ft_weights[black + i];
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
