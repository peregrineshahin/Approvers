#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

INCBIN(Network, "../default.nnue");

alignas(64) int16_t in_weights[INSIZE * L1SIZE];
alignas(64) int16_t in_biases[L1SIZE];

alignas(64) int16_t l1_weights[BUCKETS][L1SIZE * 2];
alignas(64) int16_t l1_biases[BUCKETS];

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

    for (int i = 0; i < BUCKETS; i++)
        for (int j = 0; j < L1SIZE * 2; j++)
            l1_weights[i][j] = *(data++);

    int16_t* data16 = (int16_t*) data;

    for (int i = 0; i < BUCKETS; i++)
        l1_biases[i] = *(data16++);
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
        sq ^= 7;

    return 384 * (c != side) + 64 * (pt - 1) + (side == WHITE ? sq : sq ^ 56);
}

static Value output_transform(const Accumulator* acc, const Position* pos) {
    const int bucket = (popcount(pieces()) - 2) / 4;

    const __m256i min    = _mm256_setzero_si256();
    const __m256i max    = _mm256_set1_epi16(QA);
    __m256i       vector = _mm256_setzero_si256();

    for (int flip = 0; flip <= 1; flip++)
    {
        __m256i* input   = (__m256i*) acc->values[pos->sideToMove ^ flip];
        __m256i* weights = (__m256i*) &l1_weights[bucket][flip * L1SIZE];

        for (int i = 0; i < L1SIZE / 16; ++i)
        {
            __m256i v = _mm256_min_epi16(_mm256_max_epi16(input[i], min), max);
            __m256i w = _mm256_mullo_epi16(v, weights[i]);
            vector    = _mm256_add_epi32(vector, _mm256_madd_epi16(w, v));
        }
    }

    __m256i v1     = _mm256_hadd_epi32(vector, vector);
    __m256i v2     = _mm256_hadd_epi32(v1, v1);
    Value   output = _mm256_extract_epi32(v2, 0) + _mm256_extract_epi32(v2, 4);

    return (output / QA + l1_biases[bucket]) * SCALE / (QA * QB);
}

static void refresh_accumulators(Accumulator* acc, const Position* pos) {
    memcpy(acc->values[WHITE], in_biases, sizeof(in_biases));
    memcpy(acc->values[BLACK], in_biases, sizeof(in_biases));

    const Square wksq = square_of(WHITE, KING);
    const Square bksq = square_of(BLACK, KING);

    for (Bitboard pieces = pieces(); pieces;)
    {
        const Square sq = pop_lsb(&pieces);
        const Piece  pc = piece_on(sq);
        nnue_add_piece(acc, pc, sq, wksq, bksq);
    }
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
    Accumulator* acc = pos->accumulator;

    if (acc->needs_refresh)
    {
        refresh_accumulators(acc, pos);
        acc->needs_refresh = false;
    }

    return output_transform(acc, pos);
}
