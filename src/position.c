
#include <ctype.h>
#include <inttypes.h>
#include <string.h>

#include "bitboard.h"
#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "tt.h"

static void set_state(Position* pos, Stack* st);

struct Zob zob;

const char PieceToChar[] = " PNBRQK  pnbrqk";

static void put_piece(Position* pos, Color c, Piece piece, Square s) {
    pos->board[s] = piece;
    pos->byTypeBB[0] |= sq_bb(s);
    pos->byTypeBB[type_of_p(piece)] |= sq_bb(s);
    pos->byColorBB[c] |= sq_bb(s);
}

static void remove_piece(Position* pos, Color c, Piece piece, Square s) {
    // WARNING: This is not a reversible operation.
    pos->board[s] = 0;
    pos->byTypeBB[0] ^= sq_bb(s);
    pos->byTypeBB[type_of_p(piece)] ^= sq_bb(s);
    pos->byColorBB[c] ^= sq_bb(s);
}

static void move_piece(Position* pos, Color c, Piece piece, Square from, Square to) {
    Bitboard fromToBB = sq_bb(from) ^ sq_bb(to);
    pos->byTypeBB[0] ^= fromToBB;
    pos->byTypeBB[type_of_p(piece)] ^= fromToBB;
    pos->byColorBB[c] ^= fromToBB;
    pos->board[from] = 0;
    pos->board[to]   = piece;
}


// Sets king attacks to detect if a move gives check.
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

Key H1(Key h) { return h & 0x1fff; }

Key H2(Key h) { return (h >> 16) & 0x1fff; }

static Key      cuckoo[8192];
static uint16_t cuckooMove[8192];

// Initializes at startup the various arrays used to compute hash keys.
SMALL void zob_init(void) {
    PRNG rng;
    prng_init(&rng, 1070372);

#pragma clang loop unroll(disable)
    for (int c = 0; c < 2; c++)
#pragma clang loop unroll(disable)
        for (int pt = PAWN; pt <= KING; pt++)
#pragma clang loop unroll(disable)
            for (Square s = 0; s < 64; s++)
                zob.psq[make_piece(c, pt)][s] = prng_rand(&rng);
#pragma clang loop unroll(disable)
    for (int f = 0; f < 8; f++)
        zob.enpassant[f] = prng_rand(&rng);
#pragma clang loop unroll(disable)
    for (int cr = 0; cr < 16; cr++)
        zob.castling[cr] = prng_rand(&rng);

    zob.side    = prng_rand(&rng);
    zob.noPawns = prng_rand(&rng);

    // Prepare the cuckoo tables
    int count = 0;
#pragma clang loop unroll(disable)
    for (int c = 0; c < 2; c++)
#pragma clang loop unroll(disable)
        for (int pt = PAWN; pt <= KING; pt++)
        {
            int pc = make_piece(c, pt);
#pragma clang loop unroll(disable)
            for (Square s1 = 0; s1 < 64; s1++)
#pragma clang loop unroll(disable)
                for (Square s2 = s1 + 1; s2 < 64; s2++)
                    if (PseudoAttacks[pt][s1] & sq_bb(s2))
                    {
                        //            Move move = between_bb(s1, s2) ? make_move(s1, s2)
                        //                                           : make_move(SQ_C3, SQ_D5);
                        Move     move = make_move(s1, s2);
                        Key      key  = zob.psq[pc][s1] ^ zob.psq[pc][s2] ^ zob.side;
                        uint32_t i    = H1(key);
                        while (true)
                        {
                            Key tmpKey    = cuckoo[i];
                            cuckoo[i]     = key;
                            key           = tmpKey;
                            Move tmpMove  = cuckooMove[i];
                            cuckooMove[i] = move;
                            move          = tmpMove;
                            if (!move)
                                break;
                            i = (i == H1(key)) ? H2(key) : H1(key);
                        }
                        count++;
                    }
        }
}


// Initializes the position object with the given FEN string.
SMALL void pos_set(Position* pos, char* fen) {
    unsigned char col, row, token;
    Square        sq = SQ_A8;

    Stack* st = pos->st;
    memset(pos, 0, offsetof(Position, moveList));
    pos->st = st;
    memset(st, 0, StateSize);

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
    pos->sideToMove = *fen++ == 'b';
    token           = *fen++;

    // Castling availability.
    while (((token = *fen++)) && !isspace(token))
    {
        const Color c = islower(token);

        if ((token & ~32) > 64)
        {
            const int    rights = WHITE_OO << (!(token & 2) + 2 * c);
            const Square rfrom  = relative_square(c, token & 2 ? SQ_H1 : SQ_A1);
            const Square kfrom  = square_of(c, KING);

            pos->st->castlingRights |= rights;
            pos->castlingRightsMask[kfrom] |= rights;
            pos->castlingRightsMask[rfrom] |= rights;
        }
    }

    // En passant square. Ignore if no pawn capture is possible.
    if (((col = *fen++) && (col >= 'a' && col <= 'h'))
        && ((row = *fen++) && (row == (stm() == WHITE ? '6' : '3'))))
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


// Computes the hash keys of the position, and other data that once
// computed is updated incrementally as moves are made.
// The function is only used when a new position is set up.
SMALL static void set_state(Position* pos, Stack* st) {
    st->key               = 0;
    st->pawnKey           = zob.noPawns;
    st->nonPawnKey[WHITE] = st->nonPawnKey[BLACK] = 0;
    st->majorKey = st->minorKey = 0;

    st->checkersBB = attackers_to(square_of(stm(), KING)) & pieces_c(!stm());

    set_check_info(pos);

    for (Bitboard b = pieces(); b;)
    {
        Square    s  = pop_lsb(&b);
        Piece     pc = piece_on(s);
        PieceType pt = type_of_p(pc);

        st->key ^= zob.psq[pc][s];

        if (pt == PAWN)
            st->pawnKey ^= zob.psq[piece_on(s)][s];
        else
        {
            st->nonPawnKey[color_of(pc)] ^= zob.psq[pc][s];

            if (pt == ROOK || pt == QUEEN || pt == KING)
                st->majorKey ^= zob.psq[pc][s];
            if (pt == BISHOP || pt == KNIGHT || pt == KING)
                st->minorKey ^= zob.psq[pc][s];
        }
    }

    if (st->epSquare != 0)
        st->key ^= zob.enpassant[file_of(st->epSquare)];

    if (stm() == BLACK)
        st->key ^= zob.side;

    st->key ^= zob.castling[st->castlingRights];
}


// Turning slider_blockers() into an static function was slower, even
// though it should only add a single slightly optimised copy to evaluate().
#if 1
// Returns a bitboard of all pieces that are blocking attacks on
// the square 's' from 'sliders'. A piece blocks a slider if
// removing that piece from the board would result in a position where
// square 's' is attacked. Both pinned pieces and discovered check
// candidates are slider blockers and are calculated by calling this
// function.
Bitboard slider_blockers(const Position* pos, Bitboard sliders, Square s, Bitboard* pinners) {
    Bitboard blockers = 0, snipers;
    *pinners          = 0;

    // Snipers are sliders that attack square 's'when a piece removed.
    snipers = ((PseudoAttacks[ROOK][s] & pieces_pp(QUEEN, ROOK))
               | (PseudoAttacks[BISHOP][s] & pieces_pp(QUEEN, BISHOP)))
            & sliders;
    Bitboard occupancy = pieces() ^ snipers;

    while (snipers)
    {
        Square   sniperSq = pop_lsb(&snipers);
        Bitboard b        = between_bb(s, sniperSq) & occupancy;

        if (b && !more_than_one(b))
        {
            blockers |= b;
            if (b & pieces_c(color_of(piece_on(s))))
                *pinners |= sq_bb(sniperSq);
        }
    }
    return blockers;
}
#endif


// Tests whether a pseudo-legal move is legal
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
        return !(attackers_to_occ(pos, to, pieces() ^ sq_bb(from)) & pieces_c(!us));

    // A non-king move is legal if and only if it is not pinned or it
    // is moving along the ray towards or away from the king.
    return !(blockers_for_king(pos, us) & sq_bb(from)) || aligned(m, square_of(us, KING));
}


// Takes a random move and tests whether the move is pseudo legal.
// It is used to validate moves from TT that can be corrupted
// due to hash position key aliasing.
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
        if (!(between_bb(square_of(us, KING), lsb(checkers())) & sq_bb(to)))
            return false;
    }
    return true;
}


// Invoked by gives_check() if there are discovered
// check candidates or the move is of a special type
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
        Square   capsq = make_square(file_of(to), rank_of(from));
        Bitboard b     = pieces() ^ sq_bb(from) ^ sq_bb(capsq) ^ sq_bb(to);
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


// Makes a move, and saves all information necessary
// to a StateInfo object. The move is assumed to be legal. Pseudo-legal
// moves should be filtered out before this function is called.
void do_move(Position* pos, Move m, int givesCheck) {
    Key key = key() ^ zob.side;

    // Copy some fields of the old state to our new Stack object except the
    // ones which are going to be recalculated from scratch anyway and then
    // switch our state pointer to point to the new (ready to be updated)
    // state.
    Stack* st = ++pos->st;
    memcpy(st, st - 1, (StateCopySize + 7) & ~7);

    Accumulator* acc = &st->accumulator;
    DirtyPiece*  dp  = &st->dirtyPiece;
    dp->len          = 1;

    // Increment ply counters. Note that rule50 will be reset to zero later
    // on in case of a capture or a pawn move.
    st->plyCounters += 0x101;  // Increment both rule50 and pliesFromNull

    Color  us       = stm();
    Color  them     = !us;
    Square from     = from_sq(m);
    Square to       = to_sq(m);
    Piece  piece    = piece_on(from);
    Piece  captured = type_of_m(m) == ENPASSANT ? make_piece(them, PAWN) : piece_on(to);

    if (unlikely(type_of_m(m) == CASTLING))
    {
        int    kingSide = to > from;
        Square rfrom    = to;  // Castling is encoded as "king captures friendly rook"
        Square rto      = relative_square(us, kingSide ? SQ_F1 : SQ_D1);
        to              = relative_square(us, kingSide ? SQ_G1 : SQ_C1);

        // Remove both pieces first since squares could overlap in Chess960
        remove_piece(pos, us, piece, from);
        remove_piece(pos, us, captured, rfrom);
        put_piece(pos, us, piece, to);
        put_piece(pos, us, captured, rto);

        dp->len      = 2;
        dp->piece[1] = captured;
        dp->from[1]  = rfrom;
        dp->to[1]    = rto;

        key ^= zob.psq[captured][rfrom] ^ zob.psq[captured][rto];
        st->nonPawnKey[us] ^= zob.psq[captured][rfrom] ^ zob.psq[captured][rto];
        st->majorKey ^= zob.psq[captured][rfrom] ^ zob.psq[captured][rto];
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
                capsq ^= 8;

            st->pawnKey ^= zob.psq[captured][capsq];
        }
        else
        {
            st->nonPawnKey[them] ^= zob.psq[captured][capsq];

            if (type_of_p(captured) >= ROOK)
                st->majorKey ^= zob.psq[captured][capsq];
            else
                st->minorKey ^= zob.psq[captured][capsq];
        }

        dp->len      = 2;
        dp->piece[1] = captured;
        dp->from[1]  = capsq;
        dp->to[1]    = SQ_NONE;

        // Update board and piece lists
        remove_piece(pos, them, captured, capsq);

        // Update hash key
        key ^= zob.psq[captured][capsq];

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
        key ^= zob.castling[st->castlingRights];
        st->castlingRights &= ~(pos->castlingRightsMask[from] | pos->castlingRightsMask[to]);
        key ^= zob.castling[st->castlingRights];
    }

    dp->piece[0] = piece;
    dp->from[0]  = from;
    dp->to[0]    = to;

    // Move the piece. The tricky Chess960 castling is handled earlier.
    if (likely(type_of_m(m) != CASTLING))
        move_piece(pos, us, piece, from, to);

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

            dp->to[0]          = SQ_NONE;
            dp->piece[dp->len] = promotion;
            dp->from[dp->len]  = SQ_NONE;
            dp->to[dp->len]    = to;
            dp->len++;

            // Update hash keys
            key ^= zob.psq[piece][to] ^ zob.psq[promotion][to];
            st->pawnKey ^= zob.psq[piece][to];

            if (type_of_p(promotion) >= ROOK)
                st->majorKey ^= zob.psq[promotion][to];
            else
                st->minorKey ^= zob.psq[promotion][to];
        }

        // Update pawn hash key
        st->pawnKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];

        // Reset ply counters.
        st->plyCounters = 0;
    }
    else
    {
        st->nonPawnKey[us] ^= zob.psq[piece][from] ^ zob.psq[piece][to];

        if (type_of_p(piece) == KING)
        {
            st->majorKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];
            st->minorKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];
        }
        else if (type_of_p(piece) >= ROOK)
            st->majorKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];
        else
            st->minorKey ^= zob.psq[piece][from] ^ zob.psq[piece][to];
    }

    // Update the key with the final value
    st->key = key;

    // Calculate checkers bitboard (if move gives check)
    st->checkersBB = givesCheck ? attackers_to(square_of(them, KING)) & pieces_c(us) : 0;

    pos->sideToMove = !pos->sideToMove;
    pos->nodes++;

    set_check_info(pos);

    if (type_of_p(piece) == KING && (from & 4) != (to & 4))
        acc->needs_refresh = true;
    else
        update_accumulator(acc, pos);
}


// Unmakes a move. When it returns, the position should
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


// Used to do a "null move": it flips
// the side to move without executing any move on the board.
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


// Must be used to undo a "null move"
void undo_null_move(Position* pos) {
    pos->st--;
    pos->sideToMove = !pos->sideToMove;
}


// Tests if the SEE (Static Exchange Evaluation) value of move
// is greater or equal to the given threshold. We'll use an
// algorithm similar to alpha-beta pruning with a null window.
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


// Tests whether the position is drawn by 50-move rule or by repetition.
// It does not detect stalemates.
SMALL bool is_draw(const Position* pos) {
    Stack* st = pos->st;

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


// Tests whether there has been at least one repetition
// of positions since the last capture or pawn move.
bool has_game_cycle(const Position* pos, int ply) {
    unsigned int j;

    int end = pos->st->pliesFromNull;

    Key    originalKey = pos->st->key;
    Stack* stp         = pos->st - 1;

    for (int i = 3; i <= end; i += 2)
    {
        stp -= 2;

        Key moveKey = originalKey ^ stp->key;
        if ((j = H1(moveKey), cuckoo[j] == moveKey) || (j = H2(moveKey), cuckoo[j] == moveKey))
        {
            Move m = cuckooMove[j];
            if (!((((Bitboard*) BetweenBB)[m] ^ sq_bb(to_sq(m))) & pieces()))
            {
                if (ply > i
                    || color_of(piece_on(is_empty(from_sq(m)) ? to_sq(m) : from_sq(m))) == stm())
                    return true;
            }
        }
    }
    return false;
}

void pos_set_check_info(Position* pos) { set_check_info(pos); }
