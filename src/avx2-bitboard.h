#include "types.h"
#include <immintrin.h>

extern __m256i queen_mask_v4[64][2];
extern __m256i bishop_mask_v4[64];
extern __m128i rook_mask_NS[64];
extern uint8_t rook_attacks_EW[64 * 8];

#undef attacks_bb_queen

// attacks_bb() returns a bitboard representing all the squares attacked
// by a piece of type Pt (bishop or rook) placed on 's'. The helper
// magic_index() looks up the index using the 'magic bitboards' approach.

// avx2/sse2 versions of BLSMSK (https://www.chessprogramming.org/BMI1#BLSMSK)
static __m256i blsmsk64x4(__m256i y) {
    return _mm256_xor_si256(_mm256_add_epi64(y, _mm256_set1_epi64x(-1)), y);
}

static __m128i blsmsk64x2(__m128i x) {
    return _mm_xor_si128(_mm_add_epi64(x, _mm_set1_epi64x(-1)), x);
}

static Bitboard attacks_bb_queen(Square s, Bitboard occupied) {
    const __m256i occupied4 = _mm256_set1_epi64x(occupied);
    const __m256i lmask     = queen_mask_v4[s][0];
    const __m256i rmask     = queen_mask_v4[s][1];
    __m256i       slide4, rslide;
    __m128i       slide2;

    // Left bits: set mask bits lower than occupied LS1B
    slide4 = _mm256_and_si256(occupied4, lmask);
    slide4 = _mm256_and_si256(blsmsk64x4(slide4), lmask);
    // Right bits: set shadow bits lower than occupied MS1B (6 bits max)
    rslide = _mm256_and_si256(occupied4, rmask);
    rslide =
      _mm256_or_si256(_mm256_srlv_epi64(rslide, _mm256_set_epi64x(14, 18, 16, 2)),  // PP Fill
                      _mm256_srlv_epi64(rslide, _mm256_set_epi64x(7, 9, 8, 1)));
    rslide = _mm256_or_si256(
      _mm256_srlv_epi64(rslide, _mm256_set_epi64x(28, 36, 32, 4)),
      _mm256_or_si256(rslide, _mm256_srlv_epi64(rslide, _mm256_set_epi64x(14, 18, 16, 2))));
    // add mask bits higher than blocker
    slide4 = _mm256_or_si256(slide4, _mm256_andnot_si256(rslide, rmask));

    // OR 4 vectors
    slide2 = _mm_or_si128(_mm256_castsi256_si128(slide4), _mm256_extracti128_si256(slide4, 1));
    return _mm_cvtsi128_si64(_mm_or_si128(slide2, _mm_unpackhi_epi64(slide2, slide2)));
}

static Bitboard attacks_bb_rook(Square s, Bitboard occupied) {
    // flip vertical to simulate MS1B by LS1B
    const __m128i swapl2h   = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i       occupied2 = _mm_shuffle_epi8(_mm_cvtsi64_si128(occupied), swapl2h);
    const __m128i mask      = rook_mask_NS[s];
    // set mask bits lower than occupied LS1B
    __m128i       slide2 = _mm_and_si128(blsmsk64x2(_mm_and_si128(occupied2, mask)), mask);
    const __m128i swaph2l =
      _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 8, 9, 10, 11, 12, 13, 14, 15);
    Bitboard slides = _mm_cvtsi128_si64(_mm_or_si128(slide2, _mm_shuffle_epi8(slide2, swaph2l)));

    // East-West: from precomputed table
    int r8 = rank_of(s) * 8;
    slides |= (Bitboard) (rook_attacks_EW[((occupied >> r8) & 0x7e) * 4 + file_of(s)]) << r8;
    return slides;
}

static Bitboard attacks_bb_bishop(Square s, Bitboard occupied) {
    // flip vertical to simulate MS1B by LS1B
    const __m128i swapl2h   = _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 3, 2, 1, 0);
    __m128i       occupied2 = _mm_shuffle_epi8(_mm_cvtsi64_si128(occupied), swapl2h);

    __m128i mask_lo = _mm_load_si128((__m128i*) &bishop_mask_v4[s]);
    __m128i mask_hi = _mm_load_si128((__m128i*) &bishop_mask_v4[s] + 1);

    __m128i slide_lo = _mm_and_si128(blsmsk64x2(_mm_and_si128(occupied2, mask_lo)), mask_lo);
    __m128i slide_hi = _mm_and_si128(blsmsk64x2(_mm_and_si128(occupied2, mask_hi)), mask_hi);

    __m128i slide2 = _mm_or_si128(slide_lo, slide_hi);

    __m128i swaph2l = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 8, 9, 10, 11, 12, 13, 14, 15);
    return _mm_cvtsi128_si64(_mm_or_si128(slide2, _mm_shuffle_epi8(slide2, swaph2l)));
}
