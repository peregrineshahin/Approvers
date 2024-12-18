extern Bitboard  RookMasks[64];
extern Bitboard  RookMagics[64];
extern uint8_t   RookShifts[64];
extern Bitboard  BishopMasks[64];
extern Bitboard  BishopMagics[64];
extern uint8_t   BishopShifts[64];
extern Bitboard* RookAttacks[64];
extern Bitboard* BishopAttacks[64];

// attacks_bb() returns a bitboard representing all the squares attacked
// by a // piece of type Pt (bishop or rook) placed on 's'. The helper
// magic_index() looks up the index using the 'magic bitboards' approach.

static unsigned magic_index_bishop(Square s, Bitboard occupied) {
    return (unsigned) (((occupied & BishopMasks[s]) * BishopMagics[s]) >> BishopShifts[s]);
}

static unsigned magic_index_rook(Square s, Bitboard occupied) {
    return (unsigned) (((occupied & RookMasks[s]) * RookMagics[s]) >> RookShifts[s]);
}

static Bitboard attacks_bb_bishop(Square s, Bitboard occupied) {
    return BishopAttacks[s][magic_index_bishop(s, occupied)];
}

static Bitboard attacks_bb_rook(Square s, Bitboard occupied) {
    return RookAttacks[s][magic_index_rook(s, occupied)];
}
