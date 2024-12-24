#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../small.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t l1_weights[L1SIZE * OUTSIZE * 2];

alignas(64) int16_t in_biases[L1SIZE];
alignas(64) int16_t l1_biases[OUTSIZE];

SMALL void nnue_init() {
    int8_t* data8 = (int8_t*) gNetworkData;

    for (int i = 0; i < INSIZE * L1SIZE; i++)
    {
        int x = i / L1SIZE;
        if (!(x < 8 || (56 <= x && x < 64) || (384 <= x && x < 392) || (440 <= x && x < 448)
              || (320 <= x && x < 384 && (x - 320) % 8 > 3)))
            in_weights[i] = *(data8++);
    }

    int16_t* data16 = (int16_t*) data8;

    for (int i = 0; i < L1SIZE; i++)
        in_biases[i] = *(data16++);

    for (int i = 0; i < L1SIZE * OUTSIZE * 2; i++)
        l1_weights[i] = *(data16++);

    for (int i = 0; i < OUTSIZE; i++)
        l1_biases[i] = *(data16++);
}

SMALL static Value forward(const int16_t* acc, const int16_t* weights) {
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

SMALL static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (file_of(ksq) > 3)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

SMALL static Value output_transform(const Accumulator* acc, const Position* pos) {
    const int16_t* stm  = acc->values[pos->sideToMove];
    const int16_t* nstm = acc->values[!pos->sideToMove];

    Value output = forward(stm, l1_weights) + forward(nstm, l1_weights + L1SIZE);
    return (output / QA + l1_biases[0]) * SCALE / (QA * QB);
}

SMALL static void build_accumulator(Accumulator* acc, const Position* pos, Color side) {
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

SMALL void nnue_add_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] += in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] += in_weights[black * L1SIZE + i];
    }
}

SMALL void nnue_remove_piece(Accumulator* acc, Piece pc, Square sq, Square wksq, Square bksq) {
    const int white = make_index(type_of_p(pc), color_of(pc), sq, wksq, WHITE);
    const int black = make_index(type_of_p(pc), color_of(pc), sq, bksq, BLACK);

    for (int i = 0; i < L1SIZE; i++)
    {
        acc->values[WHITE][i] -= in_weights[white * L1SIZE + i];
        acc->values[BLACK][i] -= in_weights[black * L1SIZE + i];
    }
}

SMALL Value nnue_evaluate(Position* pos) {
    Accumulator* acc = &pos->st->accumulator;

    if (acc->needs_refresh)
    {
        build_accumulator(acc, pos, WHITE);
        build_accumulator(acc, pos, BLACK);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
