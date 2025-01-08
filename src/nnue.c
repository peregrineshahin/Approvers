#include <immintrin.h>
#include <stdalign.h>
#include <stdint.h>

#include "incbin.h"
#include "nnue.h"

#include "bitboard.h"
#include "position.h"

static FinnyEntry FinnyTable[2];

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

    for (int i = 0; i < 2; i++)
    {
        memcpy(&FinnyTable[i].accumulator.values[0], in_biases, sizeof(in_biases));
        memcpy(&FinnyTable[i].accumulator.values[1], in_biases, sizeof(in_biases));
    }
}

static int make_index(PieceType pt, Color c, Square sq, Square ksq, Color side) {
    if (ksq & 4)
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
    Square ksq = square_of(side, KING);

    // 64 bucket
    // FinnyEntry* entry = &FinnyTable[ksq];

    // 16 buckets
    // FinnyEntry* entry = &FinnyTable[rank_of(ksq) + 8 * (file_of(ksq) > 3)];

    // 2 buckets
    FinnyEntry* entry = &FinnyTable[file_of(ksq) > 3];

    __m256i* values                 = (__m256i*) entry->accumulator.values[side];
    __m256i  registers[L1SIZE / 16] = {values[0], values[1], values[2], values[3]};

    for (int c = WHITE; c <= BLACK; c++)
    {
        for (int pt = PAWN; pt <= KING; pt++)
        {
            Bitboard pieces = pieces_cp(c, pt);
            Bitboard adds   = pieces & ~entry->occupancy[side][c][pt - 1];
            Bitboard subs   = ~pieces & entry->occupancy[side][c][pt - 1];

            while (adds)
            {
                const int      idx     = make_index(pt, c, pop_lsb(&adds), ksq, side) * L1SIZE;
                const __m256i* weights = (__m256i*) &in_weights[idx];

                registers[0] = _mm256_adds_epi16(registers[0], weights[0]);
                registers[1] = _mm256_adds_epi16(registers[1], weights[1]);
                registers[2] = _mm256_adds_epi16(registers[2], weights[2]);
                registers[3] = _mm256_adds_epi16(registers[3], weights[3]);
            }

            while (subs)
            {
                const int      idx     = make_index(pt, c, pop_lsb(&subs), ksq, side) * L1SIZE;
                const __m256i* weights = (__m256i*) &in_weights[idx];

                registers[0] = _mm256_subs_epi16(registers[0], weights[0]);
                registers[1] = _mm256_subs_epi16(registers[1], weights[1]);
                registers[2] = _mm256_subs_epi16(registers[2], weights[2]);
                registers[3] = _mm256_subs_epi16(registers[3], weights[3]);
            }

            entry->occupancy[side][c][pt - 1] = pieces;
        }
    }

    values[0] = registers[0];
    values[1] = registers[1];
    values[2] = registers[2];
    values[3] = registers[3];

    memcpy(acc->values[side], entry->accumulator.values[side], sizeof(acc->values[side]));
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
