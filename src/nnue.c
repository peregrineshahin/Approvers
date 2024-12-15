#include <immintrin.h>
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

static Value forward(const int16_t* acc, const int16_t* weights) {
    const __m256i min = _mm256_setzero_si256();
    const __m256i max = _mm256_set1_epi16(QA);
    __m256i vector = _mm256_setzero_si256();

    for (int i = 0; i < L1SIZE; i += 16) {
        __m256i v = _mm256_load_si256((__m256i *)(acc + i));
        v = _mm256_min_epi16(_mm256_max_epi16(v, min), max);

        const __m256i w = _mm256_load_si256((__m256i *)(weights + i));
        const __m256i product = _mm256_madd_epi16(_mm256_mullo_epi16(v, w), v);

        vector = _mm256_add_epi32(vector, product);
    }

    const __m128i upper_half = _mm256_extracti128_si256(vector, 1);
    const __m128i lower_half = _mm256_castsi256_si128(vector);

    const __m128i sum_128 = _mm_add_epi32(upper_half, lower_half);
    const __m128i sum_64 = _mm_add_epi32(_mm_unpackhi_epi64(sum_128, sum_128), sum_128);

    const __m128i shuffled = _mm_shuffle_epi32(sum_64, 1);
    const __m128i sum = _mm_add_epi32(shuffled, sum_64);

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

static void append_changed_indices(const Position* pos, const Color c, const DirtyPiece* dp, IndexList* removed, IndexList* added) {
    Square ksq = square_of(c, KING);
    for (int i = 0; i < dp->len; i++)
    {
        Piece pc = dp->piece[i];
        if (dp->from[i] != SQ_NONE)
            removed->values[removed->size++] = make_index(type_of_p(pc), color_of(pc), dp->from[i], ksq, c);

        if (dp->to[i] != SQ_NONE)
            added->values[added->size++] = make_index(type_of_p(pc), color_of(pc), dp->to[i], ksq, c);
    }
}

static void append_active_indices(const Position* pos, const Color c, IndexList* active) {
    Square   ksq = square_of(c, KING);
    Bitboard bb  = pieces();

    while (bb)
    {
        Square s  = pop_lsb(&bb);
        Piece  p  = piece_on(s);
        active->values[active->size++] = make_index(type_of_p(p), color_of(p), s, ksq, c);
    }
}

static void update_accumulator(const Position* pos, const Color c) {
    Stack* st   = pos->st;
    int    gain = popcount(pieces()) - 2;
    while (st->accumulator.state[c] == ACC_EMPTY)
    {
        DirtyPiece* dp = &st->dirtyPiece;
        if ((gain -= dp->len + 1) < 0)
            break;

        if (dp->piece[0] == make_piece(c, KING) && file_of(dp->from[0]) != file_of(dp->to[0]))
            break;

        st--;
    }

    if (st->accumulator.state[c] == ACC_COMPUTED)
    {
        if (st == pos->st)
            return;

        IndexList added[2], removed[2];
        added[0].size = added[1].size = removed[0].size = removed[1].size = 0;

        append_changed_indices(pos, c, &(st + 1)->dirtyPiece, &removed[0], &added[0]);
        for (Stack* st2 = st + 2; st2 <= pos->st; st2++)
            append_changed_indices(pos, c, &st2->dirtyPiece, &removed[1], &added[1]);

        (st + 1)->accumulator.state[c] = ACC_COMPUTED;
        pos->st->accumulator.state[c]  = ACC_COMPUTED;

        Stack* stack[3] = {st + 1, st + 1 == pos->st ? NULL : pos->st, NULL};

        for (unsigned l = 0; stack[l]; l++)
        {
            memcpy(&stack[l]->accumulator.values[c], &st->accumulator.values[c], L1SIZE * sizeof(int16_t));
            st = stack[l];

            for (unsigned k = 0; k < removed[l].size; k++)
            {
                unsigned       index  = removed[l].values[k];
                const unsigned offset = L1SIZE * index;

                for (unsigned j = 0; j < L1SIZE; j++)
                    st->accumulator.values[c][j] -= in_weights[offset + j];
            }

            for (unsigned k = 0; k < added[l].size; k++)
            {
                unsigned       index  = added[l].values[k];
                const unsigned offset = L1SIZE * index;

                for (unsigned j = 0; j < L1SIZE; j++)
                    st->accumulator.values[c][j] += in_weights[offset + j];
            }
        }
    }
    else
    {
        Accumulator* accumulator = &pos->st->accumulator;
        accumulator->state[c]    = ACC_COMPUTED;
        IndexList active;
        active.size = 0;
        append_active_indices(pos, c, &active);

        memcpy(accumulator->values[c], in_biases, L1SIZE * sizeof(int16_t));

        for (unsigned k = 0; k < active.size; k++)
        {
            unsigned index  = active.values[k];
            unsigned offset = L1SIZE * index;

            for (unsigned j = 0; j < L1SIZE; j++)
                accumulator->values[c][j] += in_weights[offset + j];
        }
    }
}

Value nnue_evaluate(Position* pos) {
    Accumulator* acc = &pos->st->accumulator;

    update_accumulator(pos, WHITE);
    update_accumulator(pos, BLACK);

    return output_transform(acc, pos);
}
