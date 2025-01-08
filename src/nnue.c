#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../expt-07.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t l1_weights[L1SIZE * 2] = {
  -41, -41, 39,  -44, -41, -42, -36, 36,  42,  -33, 31,  43,  41,  44,  -43, -42, 45,  -42, -41,
  -39, 43,  -42, -38, 39,  42,  -39, -43, -44, 43,  -41, 41,  -42, -44, 43,  41,  42,  42,  -38,
  -41, -42, 41,  42,  -37, -32, 41,  42,  -41, 38,  -43, 39,  40,  -38, -41, 41,  -42, 38,  39,
  -42, 37,  40,  41,  40,  -42, -42, 43,  42,  -41, 41,  41,  41,  41,  -37, -31, 34,  -31, -44,
  -44, -44, 43,  41,  -42, 42,  42,  42,  -43, 42,  40,  -44, -44, 41,  42,  42,  -41, 43,  -42,
  43,  43,  -33, -40, -43, -42, 40,  41,  41,  -39, -42, 42,  35,  -42, -42, 41,  -42, 44,  -40,
  -42, 43,  41,  -40, 40,  -41, -42, 40,  -34, -42, -41, -42, 41,  42};

alignas(64) int16_t in_biases[L1SIZE] = {
  124, 75,  96,  -92, -48, 58, 11,  63,  -50, 16, 36, 100, 61,  72,  3,  52,
  68,  38,  0,   8,   51,  1,  16,  43,  -4,  57, 33, 36,  29,  21,  78, 48,
  98,  0,   34,  22,  100, 14, 50,  30,  -5,  57, 4,  -19, -46, 67,  48, 108,
  -3,  -46, -19, 8,   -12, 11, -43, -62, -61, 23, 20, 83,  60,  122, 61, 75};
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
        *(data++);

    for (int i = 0; i < L1SIZE * 2; i++)
        *(data++);

    l1_bias = *((int16_t*) data);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (file_of(ksq) > 3)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static int reduce_add(__m256i vector) {
    __m256i v1 = _mm256_hadd_epi32(vector, vector);
    __m256i v2 = _mm256_hadd_epi32(v1, v1);
    return _mm256_extract_epi32(v2, 0) + _mm256_extract_epi32(v2, 4);
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    const __m256i min    = _mm256_setzero_si256();
    const __m256i max    = _mm256_set1_epi16(QA);
    __m256i       vector = _mm256_setzero_si256();

    for (int flip = 0; flip <= 1; flip++)
    {
        __m256i* input   = (__m256i*) acc->values[pos->sideToMove ^ flip];
        __m256i* weights = (__m256i*) &l1_weights[flip * L1SIZE];

        for (int i = 0; i < L1SIZE / 16; ++i)
        {
            __m256i v = _mm256_min_epi16(_mm256_max_epi16(input[i], min), max);
            __m256i w = _mm256_mullo_epi16(v, weights[i]);
            vector    = _mm256_add_epi32(vector, _mm256_madd_epi16(w, v));
        }
    }

    const Value output = reduce_add(vector);
    return (output / QA + l1_bias) * SCALE / (QA * QB);
}

static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
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
        build_accumulator(acc, pos, WHITE);
        build_accumulator(acc, pos, BLACK);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
