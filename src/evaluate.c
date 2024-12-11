/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2008 Tord Romstad (Glaurung author)
  Copyright (C) 2008-2015 Marco Costalba, Joona Kiiski, Tord Romstad
  Copyright (C) 2015-2018 Marco Costalba, Joona Kiiski, Gary Linscott, Tord Romstad

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "bitboard.h"
#include "evaluate.h"
#include "material.h"
#include "pawns.h"

// Struct EvalInfo contains various information computed and collected
// by the evaluation functions.
struct EvalInfo {
  MaterialEntry *me;
  PawnEntry *pe;
  Bitboard mobilityArea[2];

  // attackedBy[color][piece type] is a bitboard representing all squares
  // attacked by a given color and piece type. A special "piece type" which
  // is also calculated is ALL_PIECES.
  Bitboard attackedBy[2][8];

  // attackedBy2[color] are the squares attacked by 2 pieces of a given
  // color, possibly via x-ray or by one pawn and one piece. Diagonal
  // x-ray through pawn or squares attacked by 2 pawns are not explicitly
  // added.
  Bitboard attackedBy2[2];

  // kingRing[color] are the squares adjacent to the king plus some other
  // very near squares, depending on king position.
  Bitboard kingRing[2];

  // kingAttackersCount[color] is the number of pieces of the given color
  // which attack a square in the kingRing of the enemy king.
  int kingAttackersCount[2];

  // kingAttackersWeight[color] is the sum of the "weights" of the pieces
  // of the given color which attack a square in the kingRing of the enemy
  // king. The weights of the individual piece types are given by the
  // elements in the KingAttackWeights array.
  int kingAttackersWeight[2];

  // kingAttacksCount[color] is the number of attacks by the given color to
  // squares directly adjacent to the enemy king. Pieces which attack more
  // than one square are counted multiple times. For instance, if there is
  // a white knight on g5 and black's king is on g8, this white knight adds
  // 2 to kingAttacksCount[WHITE].
  int kingAttacksCount[2];
};

typedef struct EvalInfo EvalInfo;

static const Bitboard KingFlank[8] = {
  QueenSide ^ FileDBB, QueenSide, QueenSide, CenterFiles,
  CenterFiles, KingSide, KingSide, KingSide ^ FileEBB
};

// Thresholds for lazy and space evaluation
enum {
  LazyThreshold1 =  1400,
  LazyThreshold2 =  1300,
  SpaceThreshold = 12222,
};

// KingAttackWeights[PieceType] contains king attack weights by piece type
static const int KingAttackWeights[8] = { 0, 0, 81, 52, 44, 10 };

// Penalties for enemy's safe checks
enum {
  QueenSafeCheck  = 772,
  RookSafeCheck   = 1084,
  BishopSafeCheck = 645,
  KnightSafeCheck = 792
};

#define V(v) (Value)(v)
#define S(mg,eg) make_score(mg,eg)

// MobilityBonus[PieceType-2][attacked] contains bonuses for middle and
// end game, indexed by piece type and number of attacked squares in the
// mobility area.
static const Score MobilityBonus[4][32] = {
  // Knight
  { S(-62,-81), S(-53,-56), S(-12,-31), S( -4,-16), S(  3,  5), S( 13, 11),
    S( 22, 17), S( 28, 20), S( 33, 25) },
  // Bishop
  { S(-48,-59), S(-20,-23), S( 16, -3), S( 26, 13), S( 38, 24), S( 51, 42),
    S( 55, 54), S( 63, 57), S( 63, 65), S( 68, 73), S( 81, 78), S( 81, 86),
    S( 91, 88), S( 98, 97) },
  // Rook
  { S(-60,-78), S(-20,-17), S(  2, 23), S(  3, 39), S(  3, 70), S( 11, 99), // Rook
    S( 22,103), S( 31,121), S( 40,134), S( 40,139), S( 41,158), S( 48,164),
    S( 57,168), S( 57,169), S( 62,172) },
  // Queen
  { S(-30,-48), S(-12,-30), S( -8, -7), S( -9, 19), S( 20, 40), S( 23, 55), // Queen
    S( 23, 59), S( 35, 75), S( 38, 78), S( 53, 96), S( 64, 96), S( 65,100),
    S( 65,121), S( 66,127), S( 67,131), S( 67,133), S( 72,136), S( 72,141),
    S( 77,147), S( 79,150), S( 93,151), S(108,168), S(108,168), S(108,171),
    S(110,182), S(114,182), S(114,192), S(116,219) }
};

// RookOnFile[semiopen/open] contains bonuses for each rook when there is
// no friendly pawn on the rook file.
static const Score RookOnFile[2] = { S(19, 7), S(48, 27) };

// ThreatByMinor/ByRook[attacked PieceType] contains bonuses according to
// which piece type attacks which one. Attacks on lesser pieces which are
// pawn defended are not considered.
static const Score ThreatByMinor[8] = {
  S(0, 0), S(5, 32), S(55, 41), S(77, 56), S(89,119), S(79,162)
};

static const Score ThreatByRook[8] = {
  S(0, 0), S(3, 44), S(37, 68), S(42, 60), S( 0, 39), S(58, 43)
};

// PassedRank[mg/eg][Rank] contains midgame and endgame bonuses for passed
// pawns. We don't use a Score because we process the two components
// independently.
static const Value PassedRank[][8] = {
  { V(0), V( 9), V(15), V(17), V(64), V(171), V(277) },
  { V(0), V(28), V(31), V(39), V(70), V(177), V(260) }
};

// PassedFile[File] contains a bonus according to the file of a passed pawn
static const Score PassedFile[8] = {
  S( 0,  0), S(11,  8), S(22, 16), S(33, 24),
  S(33, 24), S(22, 16), S(11,  8), S( 0,  0)
};

// Assorted bonuses and penalties used by evaluation
static const Score BadOutpost          = S( -7, 36);
static const Score BishopKingProtector = S(  6,  9);
static const Score BishopOnKingRing    = S( 24,  0);
static const Score BishopOutpost       = S( 31, 23);
static const Score BishopPawns         = S(  3,  7);
static const Score BishopXRayPawns     = S(  4,  5);
static const Score CorneredBishop      = S( 50, 50);
static const Score FlankAttacks        = S(  8,  0);
static const Score Hanging             = S( 69, 36);
static const Score KnightKingProtector = S(  8,  9);
static const Score KnightOnQueen       = S( 16, 11);
static const Score KnightOutpost       = S( 56, 34);
static const Score LongDiagonalBishop  = S( 45,  0);
static const Score MinorBehindPawn     = S( 18,  3);
static const Score PawnlessFlank       = S( 17, 95);
static const Score ReachableOutpost    = S( 31, 22);
static const Score RestrictedPiece     = S(  7,  7);
static const Score RookOnKingRing      = S( 16,  0);
static const Score RookOnQueenFile     = S(  6, 11);
static const Score SliderOnQueen       = S( 60, 18);
static const Score ThreatByKing        = S( 24, 89);
static const Score ThreatByPawnPush    = S( 48, 39);
static const Score ThreatBySafePawn    = S(173, 94);
static const Score TrappedRook         = S( 55, 13);
static const Score WeakQueen           = S( 56, 15);
static const Score WeakQueenProtection = S( 14,  0);

#undef S
#undef V

// eval_init() initializes king and attack bitboards for a given color
// adding pawn attacks. To be done at the beginning of the evaluation.

INLINE void evalinfo_init(const Position *pos, EvalInfo *ei, const Color Us)
{
  const Color Them = (Us == WHITE ? BLACK : WHITE);
  const int   Down = (Us == WHITE ? SOUTH : NORTH);
  const Bitboard LowRanks = (Us == WHITE ? Rank2BB | Rank3BB
                                         : Rank7BB | Rank6BB);

  const Square ksq = square_of(Us, KING);

  Bitboard dblAttackByPawn = pawn_double_attacks_bb(pieces_cp(Us, PAWN), Us);

  // Find our pawns on the first two ranks, and those which are blocked
  Bitboard b = pieces_cp(Us, PAWN) & (shift_bb(Down, pieces()) | LowRanks);

  // Squares occupied by those pawns, by our king or queen, by blockers to
  // attacks on our king or controlled by enemy pawns are excluded from the
  // mobility area
  ei->mobilityArea[Us] = ~(b | pieces_cpp(Us, KING, QUEEN) | blockers_for_king(pos, Us) | ei->pe->pawnAttacks[Them]);

  // Initialise attackedBy[] for kings and pawns
  b = ei->attackedBy[Us][KING] = attacks_from_king(square_of(Us, KING));
  ei->attackedBy[Us][PAWN] = ei->pe->pawnAttacks[Us];
  ei->attackedBy[Us][0] = b | ei->attackedBy[Us][PAWN];
  ei->attackedBy2[Us]   = (b & ei->attackedBy[Us][PAWN]) | dblAttackByPawn;

  // Init our king safety tables only if we are going to use them
  Square s = make_square(clamp(file_of(ksq), FILE_B, FILE_G),
                         clamp(rank_of(ksq), RANK_2, RANK_7));
  ei->kingRing[Us] = PseudoAttacks[KING][s] | sq_bb(s);

  ei->kingAttackersCount[Them] = popcount(ei->kingRing[Us] & ei->pe->pawnAttacks[Them]);
  ei->kingAttacksCount[Them] = ei->kingAttackersWeight[Them] = 0;

  // Remove from kingRing[] the squares defended by two pawns
  ei->kingRing[Us] &= ~dblAttackByPawn;
}


// evaluate_piece() assigns bonuses and penalties to the pieces of a given
// color and type.

INLINE Score evaluate_pieces(const Position *pos, EvalInfo *ei, Score *mobility,
    const Color Us, const int Pt)
{
  const Color Them  = (Us == WHITE ? BLACK      : WHITE);
  const int   Down  = (Us == WHITE ? SOUTH      : NORTH);
  const Bitboard OutpostRanks = (Us == WHITE ? Rank4BB | Rank5BB | Rank6BB
                                             : Rank5BB | Rank4BB | Rank3BB);

  Bitboard b, bb;
  Square s;
  Score score = SCORE_ZERO;

  ei->attackedBy[Us][Pt] = 0;

  loop_through_pieces(Us, Pt, s) {
    // Find attacked squares, including x-ray attacks for bishops and rooks
    b = Pt == BISHOP ? attacks_bb_bishop(s, pieces() ^ pieces_p(QUEEN))
      : Pt == ROOK ? attacks_bb_rook(s,
                              pieces() ^ pieces_p(QUEEN) ^ pieces_cp(Us, ROOK))
                   : attacks_from(Pt, s);

    if (blockers_for_king(pos, Us) & sq_bb(s))
      b &= LineBB[square_of(Us, KING)][s];

    ei->attackedBy2[Us] |= ei->attackedBy[Us][0] & b;
    ei->attackedBy[Us][Pt] |= b;
    ei->attackedBy[Us][0] |= b;

    if (b & ei->kingRing[Them]) {
      ei->kingAttackersCount[Us]++;
      ei->kingAttackersWeight[Us] += KingAttackWeights[Pt];
      ei->kingAttacksCount[Us] += popcount(b & ei->attackedBy[Them][KING]);
    }
    else if (Pt == ROOK && (file_bb_s(s) & ei->kingRing[Them]))
      score += RookOnKingRing;
    else if (Pt == BISHOP && (attacks_bb_bishop(s, pieces_p(PAWN)) & ei->kingRing[Them]))
      score += BishopOnKingRing;

    int mob = popcount(b & ei->mobilityArea[Us]);

    mobility[Us] += MobilityBonus[Pt - 2][mob];

    if (Pt == BISHOP || Pt == KNIGHT) {
      // Bonus if the piece is on an outpost square or can reach one.
      // Reduced bonus for knights (BadOutpost) if it has few relevant targets.
      bb = OutpostRanks & ei->attackedBy[Us][PAWN] & ~ei->pe->pawnAttacksSpan[Them];
      bb = OutpostRanks & (ei->attackedBy[Us][PAWN] | shift_bb(Down, pieces_p(PAWN)))
                        & ~ei->pe->pawnAttacksSpan[Them];
      Bitboard targets = pieces_c(Them) & ~pieces_p(PAWN);
      if (   Pt == KNIGHT
          && (bb & sq_bb(s) & ~CenterFiles) // on a side outpost
          && !(b & targets)                 // no relevant attacks
          && (!more_than_one(targets & (sq_bb(s) & QueenSide ? QueenSide : KingSide))))
        score += BadOutpost;
      else if (bb & sq_bb(s))
        score += Pt == KNIGHT ? KnightOutpost : BishopOutpost;

      else if (Pt == KNIGHT && bb & b & ~pieces_c(Us))
        score += ReachableOutpost;

      // Knight and Bishop bonus for being right behind a pawn
      if (shift_bb(Down, pieces_p(PAWN)) & sq_bb(s))
        score += MinorBehindPawn;

      // Penalty if the minor is far from the king
      score -= (Pt == KNIGHT ? KnightKingProtector : BishopKingProtector) *
                  distance(s, square_of(Us, KING));

      if (Pt == BISHOP) {
        // Penalty according to number of pawns on the same color square as
        // the bishop, bigger when the center files are blocked with pawns
        // and smaller when the bishop is outside the pawn chain
        Bitboard blocked = pieces_cp(Us, PAWN) & shift_bb(Down, pieces());

        score -= BishopPawns * pawns_on_same_color_squares(ei->pe, Us, s)
                             * (!(ei->attackedBy[Us][PAWN] & sq_bb(s)) + popcount(blocked & CenterFiles));

        // Penalty for all enemy pawns x-rayed
        score -= BishopXRayPawns * popcount(PseudoAttacks[BISHOP][s] & pieces_cp(Them, PAWN));

        // Bonus for bishop on a long diagonal which can "see" both center
        // squares
        if (more_than_one(attacks_bb_bishop(s, pieces_p(PAWN)) & Center))
          score += LongDiagonalBishop;

        // An important Chess960 pattern: a cornered bishop blocked by a
        // friendly pawn diagonally in front of it is a very serious problem,
        // especially when that pawn is also blocked.
        if (   is_chess960()
            && (s == relative_square(Us, SQ_A1) || s == relative_square(Us, SQ_H1)))
        {
          Square d = pawn_push(Us) + (file_of(s) == FILE_A ? EAST : WEST);
          if (piece_on(s + d) == make_piece(Us, PAWN))
            score -=  piece_on(s + d + pawn_push(Us))             ? CorneredBishop * 4
                    : piece_on(s + d + d) == make_piece(Us, PAWN) ? CorneredBishop * 2
                                                                  : CorneredBishop;
        }
      }
    }

    if (Pt == ROOK) {
      // Bonus for rook on the same file as a queen
      if (file_bb_s(s) & pieces_p(QUEEN))
        score += RookOnQueenFile;

      // Bonus for rook on an open or semi-open file
      if (semiopen_file(ei->pe, Us, file_of(s)))
        score += RookOnFile[!!semiopen_file(ei->pe, Them, file_of(s))];

      // Penalty when trapped by the king, even more if the king cannot castle
      else if (mob <= 3) {
        File kf = file_of(square_of(Us, KING));

        if ((kf < FILE_E) == (file_of(s) < kf))
          score -= TrappedRook * (1 + !can_castle_c(Us));
      }
    }

    if (Pt == QUEEN) {
      // Penalty if any relative pin or discovered attack against the queen
      Bitboard pinners;
      if (slider_blockers(pos, pieces_cpp(Them, ROOK, BISHOP), s, &pinners))
          score -= WeakQueen;
    }
  }

  return score;
}

// evaluate_king() assigns bonuses and penalties to a king of a given color.

INLINE Score evaluate_king(const Position *pos, EvalInfo *ei, Score *mobility,
    const Color Us)
{
  const Color Them = Us == WHITE ? BLACK : WHITE;
  const Bitboard Camp =  Us == WHITE
                       ? AllSquares ^ Rank6BB ^ Rank7BB ^ Rank8BB
                       : AllSquares ^ Rank1BB ^ Rank2BB ^ Rank3BB;

  const Square ksq = square_of(Us, KING);
  Bitboard weak, b1, b2, b3, safe, unsafeChecks = 0;
  Bitboard rookChecks, queenChecks, bishopChecks, knightChecks;
  int kingDanger = 0;

  // King shelter and enemy pawns storm
  Score score = Us == WHITE ? king_safety_white(ei->pe, pos, ksq)
                            : king_safety_black(ei->pe, pos, ksq);

  // Attacked squares defended at most once by our queen or king
  weak =  ei->attackedBy[Them][0]
        & ~ei->attackedBy2[Us]
        & ( ~ei->attackedBy[Us][0]
           | ei->attackedBy[Us][KING] | ei->attackedBy[Us][QUEEN]);

  // Analyse the safe enemy's checks which are possible on next move
  safe  = ~pieces_c(Them);
  safe &= ~ei->attackedBy[Us][0] | (weak & ei->attackedBy2[Them]);

  b1 = attacks_bb_rook(ksq, pieces() ^ pieces_cp(Us, QUEEN));
  b2 = attacks_bb_bishop(ksq, pieces() ^ pieces_cp(Us, QUEEN));

  // Enemy rooks checks
  rookChecks = b1 & ei->attackedBy[Them][ROOK] & safe;
  if (rookChecks)
    kingDanger += more_than_one(rookChecks) ? RookSafeCheck * 175/100
                                            : RookSafeCheck;
  else
    unsafeChecks |= b1 & ei->attackedBy[Them][ROOK];

  queenChecks =  (b1 | b2) & ei->attackedBy[Them][QUEEN] & safe
               & ~(ei->attackedBy[Us][QUEEN] | rookChecks);
  if (queenChecks)
    kingDanger += more_than_one(queenChecks) ? QueenSafeCheck * 145/100
                                             : QueenSafeCheck;

  // Enemy bishops checks: we count them only if they are from squares from
  // which we can't give a queen check, because queen checks are more valuable.
  bishopChecks =  b2 & ei->attackedBy[Them][BISHOP] & safe
                & ~queenChecks;
  if (bishopChecks)
    kingDanger += more_than_one(bishopChecks) ? BishopSafeCheck * 3/2
                                              : BishopSafeCheck;
  else
    unsafeChecks |= b2 & ei->attackedBy[Them][BISHOP];

  // Enemy knights checks
  knightChecks = attacks_from_knight(ksq) & ei->attackedBy[Them][KNIGHT];
  if (knightChecks & safe)
    kingDanger += more_than_one(knightChecks & safe) ? KnightSafeCheck * 162/100 
                                                     : KnightSafeCheck;
  else
    unsafeChecks |= knightChecks;

  // Find the squares that opponent attacks in our king flank, the squares
  // which they attack twice in that flank, and the squares that we defend.
  b1 = ei->attackedBy[Them][0] & KingFlank[file_of(ksq)] & Camp;
  b2 = b1 & ei->attackedBy2[Them];
  b3 = ei->attackedBy[Us][0] & KingFlank[file_of(ksq)] & Camp;

  int kingFlankAttack = popcount(b1) + popcount(b2);
  int kingFlankDefense = popcount(b3);

  kingDanger +=  ei->kingAttackersCount[Them] * ei->kingAttackersWeight[Them]
               + 185 * popcount(ei->kingRing[Us] & weak)
               + 148 * popcount(unsafeChecks)
               +  98 * popcount(blockers_for_king(pos, Us))
               +  69 * ei->kingAttacksCount[Them]
               +   3 * kingFlankAttack * kingFlankAttack / 8
               +       mg_value(mobility[Them] - mobility[Us])
               - 873 * !pieces_cp(Them, QUEEN)
               - 100 * !!(ei->attackedBy[Us][KNIGHT] & ei->attackedBy[Us][KING])
               -   6 * mg_value(score) / 8
               -   4 * kingFlankDefense
               +  37;

  // Transform the kingDanger units into a Score, and subtract it from
  // the evaluation
  if (kingDanger > 100)
    score -= make_score(kingDanger * kingDanger / 4096, kingDanger / 16);

  // Penalty when our king is on a pawnless flank
  if (!(pieces_p(PAWN) & KingFlank[file_of(ksq)]))
    score -= PawnlessFlank;

  // Penalty if king flank is under attack, potentially moving toward the king
  score -= FlankAttacks * kingFlankAttack;

  return score;
}


// evaluate_threats() assigns bonuses according to the types of the
// attacking and the attacked pieces.

INLINE Score evaluate_threats(const Position *pos, EvalInfo *ei, const Color Us)
{
  const Color Them = (Us == WHITE ? BLACK : WHITE);
  const int   Up   = (Us == WHITE ? NORTH : SOUTH);
  const Bitboard TRank3BB = (Us == WHITE ? Rank3BB : Rank6BB);

  enum { Minor, Rook };

  Bitboard b, weak, defended, nonPawnEnemies, stronglyProtected, safe;
  Score score = SCORE_ZERO;

  // Non-pawn enemies
  nonPawnEnemies = pieces_c(Them) & ~pieces_p(PAWN);

  // Squares strongly protected by the opponent, either because they attack the
  // square with a pawn or because they attack the square twice and we don't.
  stronglyProtected =  ei->attackedBy[Them][PAWN]
                     | (ei->attackedBy2[Them] & ~ei->attackedBy2[Us]);

  // Non-pawn enemies, strongly protected
  defended = nonPawnEnemies & stronglyProtected;

  // Enemies not strongly protected and under our attack
  weak = pieces_c(Them) & ~stronglyProtected & ei->attackedBy[Us][0];

  // Add a bonus according to the kind of attacking pieces
  if (defended | weak) {
    b = (defended | weak) & (ei->attackedBy[Us][KNIGHT] | ei->attackedBy[Us][BISHOP]);
    while (b)
      score += ThreatByMinor[piece_on(pop_lsb(&b)) - 8 * Them];

    b = weak & ei->attackedBy[Us][ROOK];
    while (b)
      score += ThreatByRook[piece_on(pop_lsb(&b)) - 8 * Them];

    if (weak & ei->attackedBy[Us][KING])
      score += ThreatByKing;

    b =  ~ei->attackedBy[Them][0]
       | (nonPawnEnemies & ei->attackedBy2[Us]);
    score += Hanging * popcount(weak & b);

    // Additional bonus if weak piece is only protected by a queen
    score += WeakQueenProtection * popcount(weak & ei->attackedBy[Them][QUEEN]);
  }

  // Bonus for restricting their piece moves
  b =   ei->attackedBy[Them][0]
     & ~stronglyProtected
     &  ei->attackedBy[Us][0];
  score += RestrictedPiece * popcount(b);

  // Protected or unattacked squares
  safe = ~ei->attackedBy[Them][0] | ei->attackedBy[Us][0];

  // Bonus for attacking enemy pieces with our relatively safe pawns
  b = pieces_cp(Us, PAWN) & safe;
  b = pawn_attacks_bb(b, Us) & nonPawnEnemies;
  score += ThreatBySafePawn * popcount(b);

  // Find the squares reachable by a single pawn push
  b  = shift_bb(Up, pieces_cp(Us, PAWN)) & ~pieces();
  b |= shift_bb(Up, b & TRank3BB) & ~pieces();

  // Keep only those squares which are relatively safe
  b &= ~ei->attackedBy[Them][PAWN] & safe;

  // Bonus for safe pawn threats on the next move
  b = pawn_attacks_bb(b, Us) & nonPawnEnemies;
  score += ThreatByPawnPush * popcount(b);

  // Bonus for impending threats against enemy queen
  if (piece_count(Them, QUEEN) == 1) {
    bool queenImbalance = !pieces_cp(Us, QUEEN);

    Square s = square_of(Them, QUEEN);
    safe =   ei->mobilityArea[Us]
          & ~pieces_cp(Us, PAWN)
          & ~stronglyProtected;

    b = ei->attackedBy[Us][KNIGHT] & attacks_from_knight(s);

    score += KnightOnQueen * popcount(b & safe) * (1 + queenImbalance);

    b =  (ei->attackedBy[Us][BISHOP] & attacks_from_bishop(s))
       | (ei->attackedBy[Us][ROOK  ] & attacks_from_rook(s));

    score += SliderOnQueen * popcount(b & safe & ei->attackedBy2[Us]) * (1 + queenImbalance);
  }

  return score;
}


// Helper function
INLINE int capped_distance(Square s1, Square s2)
{
  return min(distance(s1, s2), 5);
}

// evaluate_passed() evaluates the passed pawns and candidate passed
// pawns of the given color.

INLINE Score evaluate_passed(const Position *pos, EvalInfo *ei, const Color Us)
{
  const Color Them = (Us == WHITE ? BLACK : WHITE);
  const int   Up   = (Us == WHITE ? NORTH : SOUTH);
  const int   Down = (Us == WHITE ? SOUTH : NORTH);

  Bitboard b, bb, squaresToQueen, unsafeSquares, blockedPassers, helpers;
  Score score = SCORE_ZERO;

  b = ei->pe->passedPawns[Us];

  blockedPassers = b & shift_bb(Down, pieces_cp(Them, PAWN));
  if (blockedPassers) {
    helpers =  shift_bb(Up, pieces_cp(Us, PAWN))
             & ~pieces_c(Them)
             & (~ei->attackedBy2[Them] | ei->attackedBy[Us][0]);

    // Remove blocked candidate passers that don't have help to pass
    b &=  ~blockedPassers
        | shift_bb(WEST, helpers)
        | shift_bb(EAST, helpers);
  }

  while (b) {
    Square s = pop_lsb(&b);

    int r = relative_rank_s(Us, s);

    Value mbonus = PassedRank[MG][r], ebonus = PassedRank[EG][r];

    if (r > RANK_3) {
      int w = 5 * r - 13;
      Square blockSq = s + Up;

      // Adjust bonus based on the king's proximity
      ebonus += ( (capped_distance(square_of(Them, KING), blockSq) * 19) / 4
                 - capped_distance(square_of(Us, KING), blockSq) * 2 ) * w;

      // If blockSq is not the queening square then consider also a second push
      if (r != RANK_7)
        ebonus -= capped_distance(square_of(Us, KING), blockSq + Up) * w;

      // If the pawn is free to advance, then increase the bonus
      if (is_empty(blockSq)) {
        // If there is a rook or queen attacking/defending the pawn from behind,
        // consider all the squaresToQueen. Otherwise consider only the squares
        // in the pawn's path attacked or occupied by the enemy.
        squaresToQueen = forward_file_bb(Us, s);
        unsafeSquares = passed_pawn_span(Us, s);

        bb = forward_file_bb(Them, s) & pieces_pp(ROOK, QUEEN);

        if (!(pieces_c(Them) & bb))
          unsafeSquares &= ei->attackedBy[Them][0];

        // If there are no enemy attacks on passed pawn span, assign a big
        // bonus. Otherwise, assign a smaller bonus if the path to queen is
        // not attacked and an even smaller bonus if it is attacked but
        // block square is not.
        int k =  !unsafeSquares                    ? 35
               : !(unsafeSquares & squaresToQueen) ? 20
               : !(unsafeSquares & sq_bb(blockSq)) ? 9 : 0;

        // Assign a larger bonus if the block square is defended
        if ((pieces_c(Us) & bb) || (ei->attackedBy[Us][0] & sq_bb(blockSq)))
          k += 5;

        mbonus += k * w, ebonus += k * w;
      }
    } // r > RANK_3

    score += make_score(mbonus, ebonus) - PassedFile[file_of(s)];
  }

  return score;
}


// evaluate_space() computes the space evaluation for a given side. The
// space evaluation is a simple bonus based on the number of safe squares
// available for minor pieces on the central four files on ranks 2--4. Safe
// squares one, two or three squares behind a friendly pawn are counted
// twice. Finally, the space bonus is multiplied by a weight. The aim is to
// improve play on game opening.

INLINE Score evaluate_space(const Position *pos, EvalInfo *ei, const Color Us)
{
  // Early exit if, for example, bot queens or 6 minor pieces have been
  // exchanged
  if (non_pawn_material() < SpaceThreshold)
    return SCORE_ZERO;

  const Color Them = (Us == WHITE ? BLACK : WHITE);
  const int   Down = (Us == WHITE ? SOUTH : NORTH);
  const Bitboard SpaceMask =
    Us == WHITE ? (FileCBB | FileDBB | FileEBB | FileFBB) & (Rank2BB | Rank3BB | Rank4BB)
                : (FileCBB | FileDBB | FileEBB | FileFBB) & (Rank7BB | Rank6BB | Rank5BB);

  // Find the available squares for our pieces inside the SpaceMask area
  Bitboard safe =   SpaceMask
                 & ~pieces_cp(Us, PAWN)
                 & ~ei->attackedBy[Them][PAWN];

  // Find all squares which are at most three squares behind some friendly pawn
  Bitboard behind = pieces_cp(Us, PAWN);
  behind |= shift_bb(Down, behind);
  behind |= shift_bb(Down + Down, behind);

  int bonus = popcount(safe) + popcount(behind & safe & ~ei->attackedBy[Them][0]);
  int weight = popcount(pieces_c(Us)) - 3 + min(ei->pe->blockedCount, 9);
  Score score = make_score(bonus * weight * weight / 16, 0);

  return score;
}


// evaluate_winnable() adusts the mg and eg score components based on the
// known attacking/defending status of the players.
// A single value is derived from the mg and eg values and returned.
INLINE Value evaluate_winnable(const Position *pos, EvalInfo *ei, Score score)
{
  int outflanking =  distance_f(square_of(WHITE, KING), square_of(BLACK, KING))
                   - distance_r(square_of(WHITE, KING), square_of(BLACK, KING));

  bool pawnsOnBothFlanks =   (pieces_p(PAWN) & QueenSide)
                          && (pieces_p(PAWN) & KingSide);

  bool almostUnwinnable =   outflanking < 0
                         && !pawnsOnBothFlanks;

  bool infiltration =   rank_of(square_of(WHITE, KING)) > RANK_4
                     || rank_of(square_of(BLACK, KING)) < RANK_5;

  // Compute the complexity bonus for the attacking side
  int complexity =   9 * ei->pe->passedCount
                  + 12 * popcount(pieces_p(PAWN))
                  +  9 * outflanking
                  + 21 * pawnsOnBothFlanks
                  + 24 * infiltration
                  + 51 * !non_pawn_material()
                  - 43 * almostUnwinnable
                  - 110;

  Value mg = mg_value(score);
  Value eg = eg_value(score);

  // Now apply the bonus: note that we find the attacking side by extracting
  // the sign of the midgame or endgame values, and that we carefully cap the
  // bonus so that the midgame and endgame scores will never change sign after
  // the bonus.
  int u = ((mg > 0) - (mg < 0)) * clamp(complexity + 50, -abs(mg), 0);
  int v = ((eg > 0) - (eg < 0)) * max(complexity, -abs(eg));

  mg += u;
  eg += v;

  // Compute the scale factor for the winning side
  Color strongSide = eg > VALUE_DRAW ? WHITE : BLACK;
  int sf = material_scale_factor(ei->me, pos, strongSide);

  // If scale is not already specific, scale down via general heuristics
  if (sf == SCALE_FACTOR_NORMAL) {
    if (opposite_bishops(pos)) {
      if (   non_pawn_material_c(WHITE) == BishopValueMg
          && non_pawn_material_c(BLACK) == BishopValueMg)
        sf = 18 + 4 * popcount(ei->pe->passedPawns[strongSide]);
      else
        sf = 22 + 3 * popcount(pieces_c(strongSide));
    } else if (   non_pawn_material_c(WHITE) == RookValueMg
               && non_pawn_material_c(BLACK) == RookValueMg
               && piece_count(strongSide, PAWN) - piece_count(!strongSide, PAWN) <= 1
               && !(KingSide & pieces_cp(strongSide, PAWN)) != !(QueenSide & pieces_cp(strongSide, PAWN))
               && (attacks_from_king(square_of(!strongSide, KING)) & pieces_cp(!strongSide, PAWN)))
      sf = 36;
    else if (popcount(pieces_p(QUEEN)) == 1)
      sf = 37 + 3 * (pieces_cp(WHITE, QUEEN) ? piece_count(BLACK, BISHOP) + piece_count(BLACK, KNIGHT)
                                             : piece_count(WHITE, BISHOP) + piece_count(WHITE, KNIGHT));
    else
      sf = min(sf, 36 + 7 * piece_count(strongSide, PAWN));
  }

  // Interpolate between the middlegame and the scaled endgame score
  v =  mg * ei->me->gamePhase
     + eg * (PHASE_MIDGAME - ei->me->gamePhase) * sf / SCALE_FACTOR_NORMAL;
  v /= PHASE_MIDGAME;

  return v;
}


// evaluate_classical() is the classical evaluation function. It returns
// a static evaluation of the position from the point of view of the side
// to move.

static Value evaluate_classical(const Position *pos)
{

  Score mobility[2] = { SCORE_ZERO, SCORE_ZERO };
  Value v;
  EvalInfo ei;

  // Probe the material hash table
  ei.me = material_probe(pos);

  // If we have a specialized evaluation function for the current material
  // configuration, call it and return.
  if (material_specialized_eval_exists(ei.me))
    return material_evaluate(ei.me, pos);
    
  // Initialize score by reading the incrementally updated scores included
  // in the position struct (material + piece square tables) and the
  // material imbalance. Score is computed internally from the white point
  // of view.
  Score score = psq_score() + material_imbalance(ei.me) + pos->contempt;

  // Probe the pawn hash table
  ei.pe = pawn_probe(pos);
  score += ei.pe->score;

  // Early exit if score is high
#define lazy_skip(v) (abs(mg_value(score) + eg_value(score)) / 2 > v + non_pawn_material() / 64)
  if (lazy_skip(LazyThreshold1))
    goto make_v;

  // Initialize attack and king safety bitboards.
  evalinfo_init(pos, &ei, WHITE);
  evalinfo_init(pos, &ei, BLACK);

  // Evaluate all pieces but king and pawns
  score +=  evaluate_pieces(pos, &ei, mobility, WHITE, KNIGHT)
          - evaluate_pieces(pos, &ei, mobility, BLACK, KNIGHT)
          + evaluate_pieces(pos, &ei, mobility, WHITE, BISHOP)
          - evaluate_pieces(pos, &ei, mobility, BLACK, BISHOP)
          + evaluate_pieces(pos, &ei, mobility, WHITE, ROOK)
          - evaluate_pieces(pos, &ei, mobility, BLACK, ROOK)
          + evaluate_pieces(pos, &ei, mobility, WHITE, QUEEN)
          - evaluate_pieces(pos, &ei, mobility, BLACK, QUEEN);

  score += mobility[WHITE] - mobility[BLACK];

  // Evaluate kings after all other pieces because we need full attack
  // information when computing the king safety evaluation.
  score +=  evaluate_king(pos, &ei, mobility, WHITE)
          - evaluate_king(pos, &ei, mobility, BLACK);

  // Evaluate passed pawns, we need full attack information including king
  score +=  evaluate_passed(pos, &ei, WHITE)
          - evaluate_passed(pos, &ei, BLACK);

  if (lazy_skip(LazyThreshold2))
    goto make_v;

  // Evaluate tactical threats, we need full attack information including king
  score +=  evaluate_threats(pos, &ei, WHITE)
          - evaluate_threats(pos, &ei, BLACK);

  // Evaluate space for both sides, only during opening
  score +=  evaluate_space(pos, &ei, WHITE)
          - evaluate_space(pos, &ei, BLACK);

make_v:
  // Derive single value from the mg and eg parts of the score
  v = evaluate_winnable(pos, &ei, score);

  // Evaluation grain
  v = (v / 16) * 16;

  // Side to move point of view
  v = (stm() == WHITE ? v : -v) + Tempo;

  return v;
}

Value evaluate(const Position *pos)
{
  Value v;

  v = evaluate_classical(pos);

  // Damp down the evalation linearly when shuffling
  v = v * (100 - rule50_count()) / 100;

  return clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);
}
