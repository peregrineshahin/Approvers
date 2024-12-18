#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../default.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t l1_weights[L1SIZE * OUTSIZE * 2] = {
  -43, -43, 30,  -45, -44, -41, -41, 35,  43,  -43, 19,  39,  36,  41,  -40, -44, 41,  -39, -39,
  -40, 41,  -42, -35, 41,  41,  -38, -44, -36, 38,  -41, 43,  -40, -41, 41,  42,  40,  39,  -41,
  -41, -41, 37,  43,  -42, -27, 42,  42,  -42, 37,  -40, 41,  39,  -34, -42, 38,  -41, 39,  39,
  -41, 33,  43,  38,  39,  -36, -41, 41,  43,  -41, 40,  42,  40,  42,  -33, -28, 42,  -23, -42,
  -44, -41, 42,  38,  -44, 45,  43,  41,  -44, 41,  43,  -41, -41, 44,  45,  44,  -38, 41,  -41,
  40,  42,  -37, -41, -41, -39, 44,  42,  44,  -39, -41, 42,  30,  -42, -41, 41,  -43, 42,  -42,
  -42, 42,  40,  -42, 40,  -42, -40, 43,  -31, -44, -43, -42, 44,  44,
};

alignas(64) int16_t in_biases[L1SIZE] = {
  123, 49, 93, -91, -55, 41,  17,  47, -65, 1,  28,  93,  77,  68, 6,  -23, 67, 35,  14, 74, 34, -2,
  7,   42, -7, 64,  79,  30,  24,  17, 89,  62, 92,  0,   19,  20, 81, 11,  14, 14,  -7, 48, -6, 4,
  -37, 65, 79, 37,  0,   -49, -29, 27, 96,  14, -49, -60, -53, 41, 18, 81,  39, 127, 35, 61,
};
alignas(64) int16_t l1_biases[OUTSIZE] = {533};

void nnue_init() {
    int8_t* data8 = (int8_t*) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
    {
        int x = i / L1SIZE;
        if (!(x < 8 || (56 <= x && x < 64) || (384 <= x && x < 392) || (440 <= x && x < 448)
              || (320 <= x && x < 384 && (x - 320) % 8 > 3)))
            in_weights[i] = *(data8++);
    }

    // int16_t* data16 = (int16_t*) data8;
    //
    // for (int i = 0; i < L1SIZE; i++)
    //     in_biases[i] = *(data16++);
    //
    // for (int i = 0; i < L1SIZE * OUTSIZE * 2; i++)
    //     l1_weights[i] = *(data16++);
    //
    // for (int i = 0; i < OUTSIZE; i++)
    //     l1_biases[i] = *(data16++);
}

static Value forward(const int16_t* acc, const int16_t* weights) {
    const __m256i min    = _mm256_setzero_si256();
    const __m256i max    = _mm256_set1_epi16(QA);
    __m256i       vector = _mm256_setzero_si256();

    for (int i = 0; i < L1SIZE; i += 16)
    {
        __m256i v = _mm256_load_si256((__m256i*) (acc + i));
        v         = _mm256_min_epi16(_mm256_max_epi16(v, min), max);

        const __m256i w       = _mm256_load_si256((__m256i*) (weights + i));
        const __m256i product = _mm256_madd_epi16(_mm256_mullo_epi16(v, w), v);

        vector = _mm256_add_epi32(vector, product);
    }

    const __m128i upper_half = _mm256_extracti128_si256(vector, 1);
    const __m128i lower_half = _mm256_castsi256_si128(vector);

    const __m128i sum_128 = _mm_add_epi32(upper_half, lower_half);
    const __m128i sum_64  = _mm_add_epi32(_mm_unpackhi_epi64(sum_128, sum_128), sum_128);

    const __m128i shuffled = _mm_shuffle_epi32(sum_64, 1);
    const __m128i sum      = _mm_add_epi32(shuffled, sum_64);

    return _mm_cvtsi128_si32(sum);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (file_of(ksq) > 3)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    const int16_t* stm  = acc->values[pos->sideToMove];
    const int16_t* nstm = acc->values[!pos->sideToMove];

    Value output = forward(stm, l1_weights) + forward(nstm, l1_weights + L1SIZE);
    return (output / QA + l1_biases[0]) * SCALE / (QA * QB);
}

static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
    memcpy(acc->values[side], in_biases, sizeof(acc->values[side]));

    Square ksq = square_of(side, KING);
    for (int c = WHITE; c <= BLACK; c++)
    {
        for (int pt = PAWN; pt <= KING; pt++)
        {
            Bitboard pieces = pieces_cp(c, pt);
            while (pieces)
            {
                const int idx = make_index(pt, c, pop_lsb(&pieces), ksq, side);
                for (int i = 0; i < L1SIZE; i++)
                    acc->values[side][i] += in_weights[idx * L1SIZE + i];
            }
        }
    }
}

void nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] += in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] += in_weights[black * L1SIZE + i];
    }
}

void nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] -= in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] -= in_weights[black * L1SIZE + i];
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
