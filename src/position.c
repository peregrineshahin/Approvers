
#include <ctype.h>
#include <inttypes.h>
#include <string.h>

#include "bitboard.h"
#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "tt.h"

static void set_castling_right(Position* pos, Color c, Square rfrom);
static void set_state(Position* pos, Stack* st);

struct Zob zob;

Key matKey[16] = {0ULL,
                  0x5ced000000000101ULL,
                  0xe173000000001001ULL,
                  0xd64d000000010001ULL,
                  0xab88000000100001ULL,
                  0x680b000001000001ULL,
                  0x0000000000000001ULL,
                  0ULL,
                  0ULL,
                  0xf219000010000001ULL,
                  0xbb14000100000001ULL,
                  0x58df001000000001ULL,
                  0xa15f010000000001ULL,
                  0x7c94100000000001ULL,
                  0x0000000000000001ULL,
                  0ULL};

const char PieceToChar[] = " PNBRQK  pnbrqk";

static void put_piece(Position* pos, Color c, Piece piece, Square s) {
    pos->board[s] = piece;
    pos->byTypeBB[0] |= sq_bb(s);
    pos->byTypeBB[type_of_p(piece)] |= sq_bb(s);
    pos->byColorBB[c] |= sq_bb(s);
    pos->index[s]                 = pos->pieceCount[piece]++;
    pos->pieceList[pos->index[s]] = s;
}

static void remove_piece(Position* pos, Color c, Piece piece, Square s) {
    // WARNING: This is not a reversible operation.
    pos->byTypeBB[0] ^= sq_bb(s);
    pos->byTypeBB[type_of_p(piece)] ^= sq_bb(s);
    pos->byColorBB[c] ^= sq_bb(s);
    /* board[s] = 0;  Not needed, overwritten by the capturing one */
    Square lastSquare                      = pos->pieceList[--pos->pieceCount[piece]];
    pos->index[lastSquare]                 = pos->index[s];
    pos->pieceList[pos->index[lastSquare]] = lastSquare;
    pos->pieceList[pos->pieceCount[piece]] = SQ_NONE;
}

static void move_piece(Position* pos, Color c, Piece piece, Square from, Square to) {
    // index[from] is not updated and becomes stale. This works as long as
    // index[] is accessed just by known occupied squares.
    Bitboard fromToBB = sq_bb(from) ^ sq_bb(to);
    pos->byTypeBB[0] ^= fromToBB;
    pos->byTypeBB[type_of_p(piece)] ^= fromToBB;
    pos->byColorBB[c] ^= fromToBB;
    pos->board[from]               = 0;
    pos->board[to]                 = piece;
    pos->index[to]                 = pos->index[from];
    pos->pieceList[pos->index[to]] = to;
}


// Calculate CheckInfo data.

static void set_check_info(Position* pos) {
    Stack* st = pos->st;

    st->blockersForKing[WHITE] =
      slider_blockers(pos, pieces_c(BLACK), square_of(WHITE, KING), &st->pinnersForKing[WHITE]);
    st->blockersForKing[BLACK] =
      slider_blockers(pos, pieces_c(WHITE), square_of(BLACK, KING), &st->pinnersForKing[BLACK]);

    Color them = !stm();
    st->ksq    = square_of(them, KING);

    st->checkSquares[PAWN]   = attacks_from_pawn(st->ksq, them);
    st->checkSquares[KNIGHT] = attacks_from_knight(st->ksq);
    st->checkSquares[BISHOP] = attacks_from_bishop(st->ksq);
    st->checkSquares[ROOK]   = attacks_from_rook(st->ksq);
    st->checkSquares[QUEEN]  = st->checkSquares[BISHOP] | st->checkSquares[ROOK];
    st->checkSquares[KING]   = 0;
}

static Key H1(Key h) { return h & 0x1fff; }

static Key H2(Key h) { return (h >> 16) & 0x1fff; }

// zob_init() initializes at startup the various arrays used to compute
// hash keys.

void zob_init(void) {

    PRNG rng;
    prng_init(&rng, 1070372);

    for (int c = 0; c < 2; c++)
        for (int pt = PAWN; pt <= KING; pt++)
            for (Square s = 0; s < 64; s++)
                zob.psq[make_piece(c, pt)][s] = prng_rand(&rng);

    for (int f = 0; f < 8; f++)
        zob.enpassant[f] = prng_rand(&rng);

    for (int cr = 0; cr < 16; cr++)
        zob.castling[cr] = prng_rand(&rng);

    zob.side    = prng_rand(&rng);
    zob.noPawns = prng_rand(&rng);
}


// pos_set() initializes the position object with the given FEN string.
// This function is not very robust - make sure that input FENs are correct,
// this is assumed to be the responsibility of the GUI.

void pos_set(Position* pos, char* fen) {
    unsigned char col, row, token;
    Square        sq = SQ_A8;

    Stack* st = pos->st;
    memset(pos, 0, offsetof(Position, moveList));
    pos->st = st;
    memset(st, 0, StateSize);
    for (int i = 0; i < 256; i++)
        pos->pieceList[i] = SQ_NONE;
    for (int i = 0; i < 16; i++)
        pos->pieceCount[i] = 16 * i;

    // Piece placement
    while ((token = *fen++) && token != ' ')
    {
        if (token >= '0' && token <= '9')
            sq += token - '0';  // Advance the given number of files
        else if (token == '/')
            sq -= 16;
        else
        {
            for (int piece = 0; piece < 16; piece++)
                if (PieceToChar[piece] == token)
                {
                    put_piece(pos, color_of(piece), piece, sq++);
                    break;
                }
        }
    }

    // Active color
    token           = *fen++;
    pos->sideToMove = token == 'w' ? WHITE : BLACK;
    token           = *fen++;

    // Castling availability. Compatible with 3 standards: Normal FEN
    // standard, Shredder-FEN that uses the letters of the columns on which
    // the rooks began the game instead of KQkq and also X-FEN standard
    // that, in case of Chess960, // if an inner rook is associated with
    // the castling right, the castling tag is replaced by the file letter
    // of the involved rook, as for the Shredder-FEN.
    while ((token = *fen++) && !isspace(token))
    {
        Square rsq;
        int    c    = islower(token) ? BLACK : WHITE;
        Piece  rook = make_piece(c, ROOK);

        token = toupper(token);

        if (token == 'K')
            for (rsq = relative_square(c, SQ_H1); piece_on(rsq) != rook; --rsq)
                ;
        else if (token == 'Q')
            for (rsq = relative_square(c, SQ_A1); piece_on(rsq) != rook; ++rsq)
                ;
        else
            continue;

        set_castling_right(pos, c, rsq);
    }

    // En passant square. Ignore if no pawn capture is possible.
    if (((col = *fen++) && (col >= 'a' && col <= 'h'))
        && ((row = *fen++) && (row == (stm() == WHITE ? '3' : '6'))))
    {
        st->epSquare = make_square(col - 'a', row - '1');

        if (!(attackers_to(st->epSquare) & pieces_cp(stm(), PAWN)))
            st->epSquare = 0;
    }
    else
        st->epSquare = 0;

    // Halfmove clock and fullmove number
    st->rule50   = strtol(fen, &fen, 10);
    pos->gamePly = strtol(fen, NULL, 10);

    // Convert from fullmove starting from 1 to ply starting from 0,
    // handle also common incorrect FEN with fullmove = 0.
    pos->gamePly = max(2 * (pos->gamePly - 1), 0) + (stm() == BLACK);

    set_state(pos, st);
}


// set_castling_right() is a helper function used to set castling rights
// given the corresponding color and the rook starting square.

static void set_castling_right(Position* pos, Color c, Square rfrom) {
    Square kfrom = square_of(c, KING);
    int    cs    = kfrom < rfrom ? KING_SIDE : QUEEN_SIDE;
    int    cr    = (WHITE_OO << ((cs == QUEEN_SIDE) + 2 * c));

    Square kto = relative_square(c, cs == KING_SIDE ? SQ_G1 : SQ_C1);
    Square rto = relative_square(c, cs == KING_SIDE ? SQ_F1 : SQ_D1);

    pos->st->castlingRights |= cr;

    pos->castlingRightsMask[kfrom] |= cr;
    pos->castlingRightsMask[rfrom] |= cr;

    for (Square s = min(rfrom, rto); s <= max(rfrom, rto); s++)
        if (s != kfrom && s != rfrom)
            pos->castlingPath[cr] |= sq_bb(s);

    for (Square s = min(kfrom, kto); s <= max(kfrom, kto); s++)
        if (s != kfrom && s != rfrom)
            pos->castlingPath[cr] |= sq_bb(s);
}


// set_state() computes the hash keys of the position, and other data
// that once computed is updated incrementally as moves are made. The
// function is only used when a new position is set up, and to verify
// the correctness of the Stack data when running in debug mode.

static void set_state(Position* pos, Stack* st) {
    st->key = st->materialKey = 0;
    st->pawnKey               = zob.noPawns;
    st->nonPawn               = 0;
    st->psq                   = 0;

    st->checkersBB = attackers_to(square_of(stm(), KING)) & pieces_c(!stm());

    set_check_info(pos);

    for (Bitboard b = pieces(); b;)
    {
        Square s  = pop_lsb(&b);
        Piece  pc = piece_on(s);
        st->key ^= zob.psq[pc][s];
    }

    if (st->epSquare != 0)
        st->key ^= zob.enpassant[file_of(st->epSquare)];

    if (stm() == BLACK)
        st->key ^= zob.side;

    st->key ^= zob.castling[st->castlingRights];

    for (Bitboard b = pieces_p(PAWN); b;)
    {
        Square s = pop_lsb(&b);
        st->pawnKey ^= zob.psq[piece_on(s)][s];
    }

    for (PieceType pt = PAWN; pt <= KING; pt++)
    {
        st->materialKey += piece_count(WHITE, pt) * matKey[8 * WHITE + pt];
        st->materialKey += piece_count(BLACK, pt) * matKey[8 * BLACK + pt];
    }

    for (PieceType pt = KNIGHT; pt <= QUEEN; pt++)
        for (int c = 0; c < 2; c++)
            st->nonPawn += piece_count(c, pt) * NonPawnPieceValue[make_piece(c, pt)];
}


// Turning slider_blockers() into an static function was slower, even
// though it should only add a single slightly optimised copy to evaluate().
#if 1
// slider_blockers() returns a bitboard of all pieces that are blocking
// attacks on the square 's' from 'sliders'. A piece blocks a slider if
// removing that piece from the board would result in a position where
// square 's' is attacked. Both pinned pieces and discovered check
// candidates are slider blockers and are calculated by calling this
// function.

Bitboard slider_blockers(const Position* pos, Bitboard sliders, Square s, Bitboard* pinners) {
    Bitboard result = 0, snipers;
    *pinners        = 0;

    // Snipers are sliders that attack square 's'when a piece removed.
    snipers = ((PseudoAttacks[ROOK][s] & pieces_pp(QUEEN, ROOK))
               | (PseudoAttacks[BISHOP][s] & pieces_pp(QUEEN, BISHOP)))
            & sliders;
    Bitboard occupancy = pieces() ^ snipers;

    while (snipers)
    {
        Square   sniperSq = pop_lsb(&snipers);
        Bitboard b        = between_bb(s, sniperSq) & occupancy;

        if (!more_than_one(b))
        {
            result |= b;
            if (b & pieces_c(color_of(piece_on(s))))
                *pinners |= sq_bb(sniperSq);
        }
    }
    return result;
}
#endif


#if 0
// attackers_to() computes a bitboard of all pieces which attack a given
// square. Slider attacks use the occupied bitboard to indicate occupancy.

Bitboard attackers_to_occ(const Position *pos, Square s, Bitboard occupied)
{
  return  (attacks_from_pawn(s, BLACK)    & pieces_cp(WHITE, PAWN))
        | (attacks_from_pawn(s, WHITE)    & pieces_cp(BLACK, PAWN))
        | (attacks_from_knight(s)         & pieces_p(KNIGHT))
        | (attacks_bb_rook(s, occupied)   & pieces_pp(ROOK,   QUEEN))
        | (attacks_bb_bishop(s, occupied) & pieces_pp(BISHOP, QUEEN))
        | (attacks_from_king(s)           & pieces_p(KING));
}
#endif


// is_legal() tests whether a pseudo-legal move is legal

bool is_legal(const Position* pos, Move m) {

    Color  us   = stm();
    Square from = from_sq(m);
    Square to   = to_sq(m);

    // En passant captures are a tricky special case. Because they are rather
    // uncommon, we do it simply by testing whether the king is attacked after
    // the move is made.
    if (unlikely(type_of_m(m) == ENPASSANT))
    {
        Square   ksq      = square_of(us, KING);
        Square   capsq    = to ^ 8;
        Bitboard occupied = pieces() ^ sq_bb(from) ^ sq_bb(capsq) ^ sq_bb(to);

        return !(attacks_bb_rook(ksq, occupied) & pieces_cpp(!us, QUEEN, ROOK))
            && !(attacks_bb_bishop(ksq, occupied) & pieces_cpp(!us, QUEEN, BISHOP));
    }

    // Check legality of castling moves.
    if (unlikely(type_of_m(m) == CASTLING))
    {
        // to > from works both for standard chess and for Chess960.
        to       = relative_square(us, to > from ? SQ_G1 : SQ_C1);
        int step = to > from ? WEST : EAST;

        for (Square s = to; s != from; s += step)
            if (attackers_to(s) & pieces_c(!us))
                return false;

        return true;
    }

    // If the moving piece is a king, check whether the destination
    // square is attacked by the opponent. Castling moves are checked
    // for legality during move generation.
    if (pieces_p(KING) & sq_bb(from))
        return !(attackers_to(to) & pieces_c(!us));

    // A non-king move is legal if and only if it is not pinned or it
    // is moving along the ray towards or away from the king.
    return !(blockers_for_king(pos, us) & sq_bb(from)) || aligned(m, square_of(us, KING));
}


// is_pseudo_legal() takes a random move and tests whether the move is
// pseudo legal. It is used to validate moves from TT that can be corrupted
// due to SMP concurrent access or hash position key aliasing.

#if 0
int is_pseudo_legal_old(Position *pos, Move m)
{
  int us = stm();
  Square from = from_sq(m);
  Square to = to_sq(m);
  Piece pc = moved_piece(m);

  // Use a slower but simpler function for uncommon cases
  if (type_of_m(m) != NORMAL) {
    ExtMove list[MAX_MOVES];
    ExtMove *last = generate_legal(pos, list);
    for (ExtMove *p = list; p < last; p++)
      if (p->move == m)
        return 1;
    return 0;
  }

  // Is not a promotion, so promotion piece must be empty
  if (promotion_type(m) - KNIGHT != 0)
    return 0;

  // If the 'from' square is not occupied by a piece belonging to the side to
  // move, the move is obviously not legal.
  if (pc == 0 || color_of(pc) != us)
    return 0;

  // The destination square cannot be occupied by a friendly piece
  if (pieces_c(us) & sq_bb(to))
    return 0;

  // Handle the special case of a pawn move
  if (type_of_p(pc) == PAWN) {
    // We have already handled promotion moves, so destination
    // cannot be on the 8th/1st rank.
    if (!((to + 0x08) & 0x30))
      return 0;

    if (   !(attacks_from_pawn(from, us) & pieces_c(!us) & sq_bb(to)) // Not a capture
        && !((from + pawn_push(us) == to) && is_empty(to))       // Not a single push
        && !( (from + 2 * pawn_push(us) == to)              // Not a double push
           && (rank_of(from) == relative_rank(us, RANK_2))
           && is_empty(to)
           && is_empty(to - pawn_push(us))))
      return 0;
  }
  else if (!(attacks_from(pc, from) & sq_bb(to)))
    return 0;

  // Evasions generator already takes care to avoid some kind of illegal moves
  // and legal() relies on this. We therefore have to take care that the same
  // kind of moves are filtered out here.
  if (checkers()) {
    if (type_of_p(pc) != KING) {
      // Double check? In this case a king move is required
      if (more_than_one(checkers()))
        return 0;

      // Our move must be a blocking evasion or a capture of the checking piece
      if (!((between_bb(lsb(checkers()), square_of(us, KING)) | checkers()) & sq_bb(to)))
        return 0;
    }
    // In case of king moves under check we have to remove king so as to catch
    // invalid moves like b1a1 when opposite queen is on c1.
    else if (attackers_to_occ(pos, to, pieces() ^ sq_bb(from)) & pieces_c(!us))
      return 0;
  }

  return 1;
}
#endif

bool is_pseudo_legal(const Position* pos, Move m) {
    Color  us   = stm();
    Square from = from_sq(m);

    if (!(pieces_c(us) & sq_bb(from)))
        return false;

    if (unlikely(type_of_m(m) == CASTLING))
    {
        if (checkers())
            return false;
        ExtMove  list[MAX_MOVES];
        ExtMove* end = generate(pos, list, QUIETS);
        for (ExtMove* p = list; p < end; p++)
            if (p->move == m)
                return is_legal(pos, m);
        return false;
    }

    Square to = to_sq(m);
    if (pieces_c(us) & sq_bb(to))
        return false;

    PieceType pt = type_of_p(piece_on(from));
    if (pt != PAWN)
    {
        if (type_of_m(m) != NORMAL)
            return false;
        switch (pt)
        {
        case KNIGHT :
            if (!(attacks_from_knight(from) & sq_bb(to)))
                return false;
            break;
        case BISHOP :
            if (!(attacks_from_bishop(from) & sq_bb(to)))
                return false;
            break;
        case ROOK :
            if (!(attacks_from_rook(from) & sq_bb(to)))
                return false;
            break;
        case QUEEN :
            if (!(attacks_from_queen(from) & sq_bb(to)))
                return false;
            break;
        case KING :
            if (!(attacks_from_king(from) & sq_bb(to)))
                return 0;
            // is_legal() does not remove the "from" square from the "occupied"
            // bitboard when checking that the king is not in check on the "to"
            // square. So we need to be careful here.
            if (checkers() && (attackers_to_occ(pos, to, pieces() ^ sq_bb(from)) & pieces_c(!us)))
                return false;
            return true;
        default :
            break;
        }
    }
    else
    {
        if (likely(type_of_m(m) == NORMAL))
        {
            if (!((to + 0x08) & 0x30))
                return false;
            if (!(attacks_from_pawn(from, us) & pieces_c(!us) & sq_bb(to))
                && !((from + pawn_push(us) == to) && is_empty(to))
                && !((from + 2 * pawn_push(us) == to)
                     && (rank_of(from) == relative_rank(us, RANK_2)) && is_empty(to)
                     && is_empty(to - pawn_push(us))))
                return false;
        }
        else if (likely(type_of_m(m) == PROMOTION))
        {
            // No need to test for pawn to 8th rank.
            if (!(attacks_from_pawn(from, us) & pieces_c(!us) & sq_bb(to))
                && !((from + pawn_push(us) == to) && is_empty(to)))
                return false;
        }
        else
            return to == ep_square() && (attacks_from_pawn(from, us) & sq_bb(to));
    }
    if (checkers())
    {
        // Again we need to be a bit careful.
        if (more_than_one(checkers()))
            return false;
        if (!((between_bb(lsb(checkers()), square_of(us, KING)) | checkers()) & sq_bb(to)))
            return false;
    }
    return true;
}

#if 0
int is_pseudo_legal(Position *pos, Move m)
{
  int r1 = is_pseudo_legal_old(pos, m);
  int r2 = is_pseudo_legal_new(pos, m);
  if (r1 != r2) {
    printf("old: %d, new: %d\n", r1, r2);
    printf("old: %d\n", is_pseudo_legal_old(pos, m));
    printf("new: %d\n", is_pseudo_legal_new(pos, m));
exit(1);
  }
  return r1;
}
#endif


// gives_check_special() is invoked by gives_check() if there are
// discovered check candidates or the move is of a special type

bool gives_check_special(const Position* pos, Stack* st, Move m) {

    Square from = from_sq(m);
    Square to   = to_sq(m);

    if ((blockers_for_king(pos, !stm()) & sq_bb(from)) && !aligned(m, st->ksq))
        return true;

    switch (type_of_m(m))
    {
    case NORMAL :
        return st->checkSquares[type_of_p(piece_on(from))] & sq_bb(to);

    case PROMOTION :
        return attacks_bb(promotion_type(m), to, pieces() ^ sq_bb(from)) & sq_bb(st->ksq);

    case ENPASSANT : {
        if (st->checkSquares[PAWN] & sq_bb(to))
            return true;
        Square capsq = make_square(file_of(to), rank_of(from));
        //    Bitboard b = pieces() ^ sq_bb(from) ^ sq_bb(capsq) ^ sq_bb(to);
        Bitboard b = inv_sq(inv_sq(inv_sq(pieces(), from), to), capsq);
        return (attacks_bb_rook(st->ksq, b) & pieces_cpp(stm(), QUEEN, ROOK))
            || (attacks_bb_bishop(st->ksq, b) & pieces_cpp(stm(), QUEEN, BISHOP));
    }
    case CASTLING : {
        // Castling is encoded as 'King captures the rook'
        Square rto = relative_square(stm(), to > from ? SQ_F1 : SQ_D1);
        return (PseudoAttacks[ROOK][rto] & sq_bb(st->ksq))
            && (attacks_bb_rook(rto, pieces() ^ sq_bb(from)) & sq_bb(st->ksq));
    }
    default :
        return false;
    }
}


// do_move() makes a move. The move is assumed to be legal.
#include "evaluate.h"
void do_move(Position* pos, Move m, int givesCheck) {
    //if(pos->nodes==7381)
    //printf("hey\n");
    //printf("n = %lu, key = %08lx, m = %d, eval = %d\n", pos->nodes, key(), (int)m, !checkers() ? evaluate(pos) : 0);

    Key key = key() ^ zob.side;

    // Copy some fields of the old state to our new Stack object except the
    // ones which are going to be recalculated from scratch anyway and then
    // switch our state pointer to point to the new (ready to be updated)
    // state.
    Stack* st = ++pos->st;
    memcpy(st, st - 1, (StateCopySize + 7) & ~7);

    Accumulator* acc = &st->accumulator;
    memcpy(acc, &(st - 1)->accumulator, sizeof(st->accumulator));

    // Increment ply counters. Note that rule50 will be reset to zero later
    // on in case of a capture or a pawn move.
    st->plyCounters += 0x101;  // Increment both rule50 and pliesFromNull

    Color  us       = stm();
    Color  them     = !us;
    Square from     = from_sq(m);
    Square to       = to_sq(m);
    Piece  piece    = piece_on(from);
    Piece  captured = type_of_m(m) == ENPASSANT ? make_piece(them, PAWN) : piece_on(to);
    Square wksq     = square_of(WHITE, KING);
    Square bksq     = square_of(BLACK, KING);

    if (type_of_p(piece) == KING && file_of(from) > 3 != file_of(to) > 3)
        acc->needs_refresh = 1;

    if (unlikely(type_of_m(m) == CASTLING))
    {


        Square rfrom, rto;

        int kingSide = to > from;
        rfrom        = to;  // Castling is encoded as "king captures friendly rook"
        rto          = relative_square(us, kingSide ? SQ_F1 : SQ_D1);
        to           = relative_square(us, kingSide ? SQ_G1 : SQ_C1);

        // Remove both pieces first since squares could overlap in Chess960
        remove_piece(pos, us, piece, from);
        remove_piece(pos, us, captured, rfrom);
        pos->board[from] = pos->board[rfrom] = 0;
        put_piece(pos, us, piece, to);
        put_piece(pos, us, captured, rto);

        nnue_remove_piece(acc, piece, from, wksq, bksq);
        nnue_remove_piece(acc, captured, rfrom, wksq, bksq);

        nnue_add_piece(acc, piece, to, wksq, bksq);
        nnue_add_piece(acc, captured, rto, wksq, bksq);

        key ^= zob.psq[captured][rfrom] ^ zob.psq[captured][rto];
        captured = 0;
    }

    else if (captured)
    {
        Square capsq = to;

        // If the captured piece is a pawn, update pawn hash key. Otherwise,
        // update non-pawn material.
        if (type_of_p(captured) == PAWN)
        {
            if (unlikely(type_of_m(m) == ENPASSANT))
            {
                capsq ^= 8;


                pos->board[capsq] = 0;  // Not done by remove_piece()
            }

            st->pawnKey ^= zob.psq[captured][capsq];
        }
        else
            st->nonPawn -= NonPawnPieceValue[captured];

        nnue_remove_piece(acc, captured, capsq, wksq, bksq);

        // Update board and piece lists
        remove_piece(pos, them, captured, capsq);

        // Update material hash key
        key ^= zob.psq[captured][capsq];
        st->materialKey -= matKey[captured];

        // Reset ply counters
        st->plyCounters = 0;
    }

    // Set captured piece
    st->capturedPiece = captured;

    // Update hash key
    key ^= zob.psq[piece][from] ^ zob.psq[piece][to];

    // Reset en passant square
    if (unlikely((st - 1)->epSquare != 0))
        key ^= zob.enpassant[file_of((st - 1)->epSquare)];
    st->epSquare = 0;

    // Update castling rights if needed
    if (st->castlingRights && (pos->castlingRightsMask[from] | pos->castlingRightsMask[to]))
    {
        //    uint32_t cr = pos->castlingRightsMask[from] | pos->castlingRightsMask[to];
        //    key ^= zob.castling[st->castlingRights & cr];
        //    st->castlingRights &= ~cr;
        key ^= zob.castling[st->castlingRights];
        st->castlingRights &= ~(pos->castlingRightsMask[from] | pos->castlingRightsMask[to]);
        key ^= zob.castling[st->castlingRights];
    }

    // Move the piece. The tricky Chess960 castling is handled earlier.
    if (likely(type_of_m(m) != CASTLING))
    {
        move_piece(pos, us, piece, from, to);

        nnue_remove_piece(acc, piece, from, wksq, bksq);
        nnue_add_piece(acc, piece, to, wksq, bksq);
    }

    // If the moving piece is a pawn do some special extra work
    if (type_of_p(piece) == PAWN)
    {
        // Set en-passant square if the moved pawn can be captured
        if ((to ^ from) == 16 && (attacks_from_pawn(to ^ 8, us) & pieces_cp(them, PAWN)))
        {
            st->epSquare = to ^ 8;
            key ^= zob.enpassant[file_of(st->epSquare)];
        }
        else if (type_of_m(m) == PROMOTION)
        {
            Piece promotion = make_piece(us, promotion_type(m));

            remove_piece(pos, us, piece, to);
            put_piece(pos, us, promotion, to);

            nnue_remove_piece(acc, piece, to, wksq, bksq);
            nnue_add_piece(acc, promotion, to, wksq, bksq);

            // Update hash keys
            key ^= zob.psq[piece][to] ^ zob.psq[promotion][to];
            st->pawnKey ^= zob.psq[piece][to];
            st->materialKey += matKey[promotion] - matKey[piece];

            // Update material
            st->nonPawn += NonPawnPieceValue[promotion];
        }

        // Update pawn hash key
        st->pawnKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];

        // Reset ply counters.
        st->plyCounters = 0;
    }

    // Update the key with the final value
    st->key = key;

    // Calculate checkers bitboard (if move gives check)
#if 1
    st->checkersBB = givesCheck ? attackers_to(square_of(them, KING)) & pieces_c(us) : 0;
#else
    st->checkersBB = 0;
    if (givesCheck)
    {
        if (type_of_m(m) != NORMAL || ((st - 1)->blockersForKing[them] & sq_bb(from)))
            st->checkersBB = attackers_to(square_of(them, KING)) & pieces_c(us);
        else
            st->checkersBB = (st - 1)->checkSquares[piece & 7] & sq_bb(to);
    }
#endif

    pos->sideToMove = !pos->sideToMove;
    pos->nodes++;

    set_check_info(pos);
}


// undo_move() unmakes a move. When it returns, the position should
// be restored to exactly the same state as before the move was made.

void undo_move(Position* pos, Move m) {

    pos->sideToMove = !pos->sideToMove;

    Color  us   = stm();
    Square from = from_sq(m);
    Square to   = to_sq(m);
    Piece  pc   = piece_on(to);

    if (unlikely(type_of_m(m) == PROMOTION))
    {

        remove_piece(pos, us, pc, to);
        pc = make_piece(us, PAWN);
        put_piece(pos, us, pc, to);
    }

    if (unlikely(type_of_m(m) == CASTLING))
    {
        Square rfrom, rto;
        int    kingSide = to > from;
        rfrom           = to;  // Castling is encoded as "king captures friendly rook"
        rto             = relative_square(us, kingSide ? SQ_F1 : SQ_D1);
        to              = relative_square(us, kingSide ? SQ_G1 : SQ_C1);
        Piece king      = make_piece(us, KING);
        Piece rook      = make_piece(us, ROOK);

        // Remove both pieces first since squares could overlap in Chess960
        remove_piece(pos, us, king, to);
        remove_piece(pos, us, rook, rto);
        pos->board[to] = pos->board[rto] = 0;
        put_piece(pos, us, king, from);
        put_piece(pos, us, rook, rfrom);
    }
    else
    {
        move_piece(pos, us, pc, to, from);  // Put the piece back at the source square

        if (pos->st->capturedPiece)
        {
            Square capsq = to;

            if (unlikely(type_of_m(m) == ENPASSANT))
            {
                capsq ^= 8;
            }

            put_piece(pos, !us, pos->st->capturedPiece, capsq);  // Restore the captured piece
        }
    }

    // Finally, point our state pointer back to the previous state
    pos->st--;
}


// do_null_move() is used to do a null move

void do_null_move(Position* pos) {

    Stack* st = ++pos->st;
    memcpy(st, st - 1, (StateSize + 7) & ~7);

    memcpy(&st->accumulator, &(st - 1)->accumulator, sizeof(st->accumulator));

    if (unlikely(st->epSquare))
    {
        st->key ^= zob.enpassant[file_of(st->epSquare)];
        st->epSquare = 0;
    }

    st->key ^= zob.side;
    prefetch(tt_first_entry(st->key));

    st->rule50++;
    st->pliesFromNull = 0;

    pos->sideToMove = !pos->sideToMove;

    set_check_info(pos);
}

// See position.h for undo_null_move()


// key_after() computes the new hash key after the given move. Needed
// for speculative prefetch. It does not recognize special moves like
// castling, en-passant and promotions.

Key key_after(const Position* pos, Move m) {
    Square from     = from_sq(m);
    Square to       = to_sq(m);
    Piece  pc       = piece_on(from);
    Piece  captured = piece_on(to);
    Key    k        = key() ^ zob.side;

    if (captured)
        k ^= zob.psq[captured][to];

    return k ^ zob.psq[pc][to] ^ zob.psq[pc][from];
}


// Test whether SEE >= value.
bool see_test(const Position* pos, Move m, int value) {
    if (unlikely(type_of_m(m) != NORMAL))
        return 0 >= value;

    Square   from = from_sq(m), to = to_sq(m);
    Bitboard occ;

    int swap = PieceValue[piece_on(to)] - value;
    if (swap < 0)
        return false;

    swap = PieceValue[piece_on(from)] - swap;
    if (swap <= 0)
        return true;

    occ                = pieces() ^ sq_bb(from) ^ sq_bb(to);
    Color    stm       = color_of(piece_on(from));
    Bitboard attackers = attackers_to_occ(pos, to, occ), stmAttackers;
    bool     res       = true;

    while (true)
    {
        stm = !stm;
        attackers &= occ;
        if (!(stmAttackers = attackers & pieces_c(stm)))
            break;
        if ((stmAttackers & blockers_for_king(pos, stm)) && (pos->st->pinnersForKing[stm] & occ))
            stmAttackers &= ~blockers_for_king(pos, stm);
        if (!stmAttackers)
            break;
        res = !res;
        Bitboard bb;
        if ((bb = stmAttackers & pieces_p(PAWN)))
        {
            if ((swap = PawnValue - swap) < res)
                break;
            occ ^= bb & -bb;
            attackers |= attacks_bb_bishop(to, occ) & pieces_pp(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces_p(KNIGHT)))
        {
            if ((swap = KnightValue - swap) < res)
                break;
            occ ^= bb & -bb;
        }
        else if ((bb = stmAttackers & pieces_p(BISHOP)))
        {
            if ((swap = BishopValue - swap) < res)
                break;
            occ ^= bb & -bb;
            attackers |= attacks_bb_bishop(to, occ) & pieces_pp(BISHOP, QUEEN);
        }
        else if ((bb = stmAttackers & pieces_p(ROOK)))
        {
            if ((swap = RookValue - swap) < res)
                break;
            occ ^= bb & -bb;
            attackers |= attacks_bb_rook(to, occ) & pieces_pp(ROOK, QUEEN);
        }
        else if ((bb = stmAttackers & pieces_p(QUEEN)))
        {
            if ((swap = QueenValue - swap) < res)
                break;
            occ ^= bb & -bb;
            attackers |= (attacks_bb_bishop(to, occ) & pieces_pp(BISHOP, QUEEN))
                       | (attacks_bb_rook(to, occ) & pieces_pp(ROOK, QUEEN));
        }
        else  // KING
            return (attackers & ~pieces_c(stm)) ? !res : res;
    }

    return res;
}


// is_draw() tests whether the position is drawn by 50-move rule or by
// repetition. It does not detect stalemates.

SMALL
bool is_draw(const Position* pos) {
    Stack* st = pos->st;

    if (unlikely(st->rule50 > 99))
    {
        if (!checkers())
            return true;
        return generate_legal(pos, (st - 1)->endMoves) != (st - 1)->endMoves;
    }

    // st->pliesFromNull is reset both on null moves and on zeroing moves.
    int e = st->pliesFromNull - 4;
    if (e >= 0)
    {
        Stack* stp = st - 2;
        for (int i = 0; i <= e; i += 2)
        {
            stp -= 2;
            if (stp->key == st->key)
                return true;  // Draw at first repetition
        }
    }

    return false;
}


void pos_set_check_info(Position* pos) { set_check_info(pos); }
