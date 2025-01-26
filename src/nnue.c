#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../pairwise-05.bin");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t in_biases[L1SIZE];

alignas(64) int16_t l1_weights[L1SIZE];
alignas(64) int16_t l1_bias;

SMALL void nnue_init() {
    int16_t* data = (int16_t*) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
        in_weights[i] = *(data++);

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data++);

    for (int i = 0; i < L1SIZE; i++)
        l1_weights[i] = *(data++);

    l1_bias = *((int16_t*) data);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    int halfWidth = L1SIZE / 2;
    int output    = 0;

    for (int flip = 0; flip <= 1; flip++)
    {
        int16_t* input   = acc->values[pos->sideToMove ^ flip];
        int16_t* weights = &l1_weights[flip * halfWidth];

        for (int i = 0; i < halfWidth; i++)
        {
            int32_t v1 = clamp(input[i], 0, QA);
            int32_t v2 = clamp(input[i + halfWidth], 0, QA);

            output += v1 * v2 * weights[i];
        }
    }

    return (output / QA + l1_bias) * SCALE / (QA * QB);
}

static void refresh_accumulator(Accumulator* acc, const Position* pos, Color side) {
    const __m256i* biases                 = (__m256i*) in_biases;
    __m256i        registers[L1SIZE / 16] = {
      biases[0],
      biases[1],
      biases[2],
      biases[3],
    };

    const Square ksq = square_of(side, KING);
    for (Bitboard pieces = pieces(); pieces;)
    {
        const Square sq = pop_lsb(&pieces);
        const Piece  pc = piece_on(sq);

        const int      index   = make_index(type_of_p(pc), color_of(pc), sq, ksq, side);
        const __m256i* weights = (__m256i*) &in_weights[index * L1SIZE];

        registers[0] = _mm256_add_epi16(registers[0], weights[0]);
        registers[1] = _mm256_add_epi16(registers[1], weights[1]);
        registers[2] = _mm256_add_epi16(registers[2], weights[2]);
        registers[3] = _mm256_add_epi16(registers[3], weights[3]);
    }

    __m256i* values = (__m256i*) &acc->values[side];
    values[0]       = registers[0];
    values[1]       = registers[1];
    values[2]       = registers[2];
    values[3]       = registers[3];
}

void nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1SIZE;

    for (int i = 0; i < L1SIZE; i += 16)
    {
        __m256i w_acc = _mm256_load_si256((__m256i*) &acc->values[WHITE][i]);
        __m256i b_acc = _mm256_load_si256((__m256i*) &acc->values[BLACK][i]);

        w_acc = _mm256_adds_epi16(w_acc, _mm256_load_si256((__m256i*) &in_weights[white + i]));
        b_acc = _mm256_adds_epi16(b_acc, _mm256_load_si256((__m256i*) &in_weights[black + i]));

        _mm256_store_si256((__m256i*) &acc->values[WHITE][i], w_acc);
        _mm256_store_si256((__m256i*) &acc->values[BLACK][i], b_acc);
    }
}

void nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE) * L1SIZE;
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK) * L1SIZE;

    for (int i = 0; i < L1SIZE; i += 16)
    {
        __m256i w_acc = _mm256_load_si256((__m256i*) &acc->values[WHITE][i]);
        __m256i b_acc = _mm256_load_si256((__m256i*) &acc->values[BLACK][i]);

        w_acc = _mm256_subs_epi16(w_acc, _mm256_load_si256((__m256i*) &in_weights[white + i]));
        b_acc = _mm256_subs_epi16(b_acc, _mm256_load_si256((__m256i*) &in_weights[black + i]));

        _mm256_store_si256((__m256i*) &acc->values[WHITE][i], w_acc);
        _mm256_store_si256((__m256i*) &acc->values[BLACK][i], b_acc);
    }
}

Value nnue_evaluate(Position* pos) {
    Accumulator* acc = &pos->st->accumulator;

    if (acc->needs_refresh)
    {
        refresh_accumulator(acc, pos, WHITE);
        refresh_accumulator(acc, pos, BLACK);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
