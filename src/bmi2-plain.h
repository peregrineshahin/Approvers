#include <immintrin.h>

extern Bitboard  RookMasks[64];
extern Bitboard  BishopMasks[64];
extern Bitboard* RookAttacks[64];
extern Bitboard* BishopAttacks[64];

static unsigned bmi2_index_bishop(Square s, Bitboard occupied) {
    return (unsigned) _pext_u64(occupied, BishopMasks[s]);
}

static unsigned bmi2_index_rook(Square s, Bitboard occupied) {
    return (unsigned) _pext_u64(occupied, RookMasks[s]);
}

static Bitboard attacks_bb_bishop(Square s, Bitboard occupied) {
    return BishopAttacks[s][bmi2_index_bishop(s, occupied)];
}

static Bitboard attacks_bb_rook(Square s, Bitboard occupied) {
    return RookAttacks[s][bmi2_index_rook(s, occupied)];
}
