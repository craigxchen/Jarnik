# Working notes for round9/README.md (coordinator synthesis; not a report)

## Regime remark (README remark, unrefereed)
- (T') of L4 sec 1 (referee: survives): growth needs only, for every c_0, SOME A = A(c_0) with
  no member of F^all_{3M^A} with p_j <= M^{h_0}, h_0 = 2c_0 - 3 + o(1), having W <= c_0 M log M.
- Every residue-class cap of round7 README sec 2.1 rows 2, 3, 4, 4a uses p_j > X:
  - sharp.md Def 3.2 / Thm 5.1: p_j > X = 3M^A, r = b(M-1), b >= 33;
  - construct.md Thm 2: p_j >= P_* >= 160 M^2 L/C with L > Q_0 = X, r = 9M - 8;
  - lower.md (16+32A) at M^A; construct 9(2+max(K,1)) at M^K.
- p3lower Lemma D.0: all p_j > X => W >= r log X. So these caps have W >= 8(M-1) log X, constant
  c(A) >= 8A (construct) or 33A (sharp) -> infinity with A.
- Hence none of the residue caps is uniform in A, and none blocks an argument that chooses A after
  c_0 (the (T') scheme). Under (T') the barrier status of rows 2-4 in regime (1.4) is OPEN
  (not refuted: caps at each fixed A remain correct). Row 1 (A_0) cap is X-free, so stays.
- The old dictionary "capped => no argument using its inputs proves growth" is valid for arguments
  with A fixed independently of c_0.
- sharp.md sec 7 (X up to N^A) is row 8; the (T') regime is X = M^{A(c_0)}, polynomial but with
  exponent depending on c_0, and all p_j <= X (so level data at q = p_j enter, sharp.md sec 7).

## Elementary heuristic for regime (1.4)
- independent labels: E[log m_(v,s)] ~ sum_{q<=X} log q / m_q ~ log X = A log M + O(1), versus the
  budget w(v)/2 <= k h_0 log M / 2 for ||v||_1 = k. So for k < 2A/h_0 every single monomial fails
  in expectation (not only the tail). Genuine labels: m_(v,s) | D_s(A_v), an integer of size
  ~ e^{w(v)/2}, so demand <= budget exactly.
- Consistent with referee_L4 sec 6.3: Construction D in regime (1.4) violates 93-134 of 100-144
  norm-2 and 8-10 of 10-12 norm-1 monomials; genuine 0.
- Counting across monomials does not force (representation multiplicity binom(2k,k) cancels the
  4^k gain) -- consistent with p3Dlower Lemma 5.1.

## Next step candidates
1. Regime-(1.4) L4 at norm 1 and 2 (finite, exact): characterise parity-respecting systems with
   p_j <= M^{h_0}, X = 3M^A, A >> h_0, that satisfy the per-unit height condition for all
   ||v||_1 <= 2. Question: must such columns be F1-coherent with Gaussian integers of norm ~ p_j?
   Norm 1: v = e_j: collision set of column j must have log-mass <= (1/2) log p_j + O(1), out of
   log X. That is: ell_j(q) avoids the unit values at all q in a set of log-mass >= log X - ...
2. Uniform-in-A caps of A_X (point-level) in regime (1.4): construct or refute.

## Decisions after reading lower.md Lemma 4.1
- A-dependence of the A_X caps comes from lower.md Lemma 4.1: auxiliary primes ell_x in (X, Y], needed
  so that varpi_x is a unit mod every q <= X; then m_c^gen | Im(Omega_c^4)/4 <= Y^(2||c||_1), Y > X.
- Repair sketch (R2): take ell_x <= M^{O(1)}, ell_x > M+1, distinct, not among p_j; at the single bad
  modulus q = ell_x of point x choose rho_x(q) of norm N freely: distinct from all rho_y mod q (pairs: 0
  demand there), generic otherwise (non-pair cost O(||c||_1 log M) w.h.p. by a union bound).
  Level data at q = p_j <= X: same-level pairs only; choose injective. Then character demands are
  X-independent; expect A_X cap at constant independent of A (p_j ~ M^{3+o(1)} as in construct.md with
  L ~ M log M). NOT written out; flag as sketch.
- P3-char: the small coherent point classes log(varpi_x/conj varpi_x) + c_q need (ell_x N_x / q)
  constant in x for every q in Q_1 <= X; for M random-like signs this fails at almost every q (prob
  1 - 2^(1-M)); repair would be needed at almost all of Q_1, i.e. no coherence left. So the sketch
  does not transfer.
- In regime (1.4) coherent parity-respecting columns have genuine height (Thm A); with Linnik the
  rest of Norm zeta_j is a square up to 2, so ell_j = genuine + 2 log(xi_j/conj xi_j) (+F1 ops):
  coherent side = genuine-type labels = P3-G territory. Route 3 and route 4 meet in regime (1.4).
  (L4.md sec 1 already says "certificates of genuine height remain, exactly the genuine-type labels".)

## Prop H (hybrid systems), elementary
Claim: fix genuine primes p_1..p_r (=1 mod 4), K >= 1, Q_K := 8K(2r+1)^K + 1. Define ell_j(q) = genuine
log(pi_j/conj pi_j) for q <= Q_K (q not | N), and for q > Q_K choose ell_1(q),...,ell_r(q) greedily in
their forced parity classes, avoiding: for every v with v_j != 0, supp v in {1..j}, ||v||_1 <= K, every
unit s: <v, ell(q)> = s iota mod m_q (level 1).
Count: #v <= (2r+1)^K; for fixed v, s the equation v_j ell_j = c mod m has <= gcd(v_j, m) <= K solutions;
4 units -> forbidden <= 4K(2r+1)^K < m_q/2 = size of parity class (m_q >= q-1 > 8K(2r+1)^K).
Hence: no monomial of norm <= K collides at any level at any q > Q_K. Its per-unit demand is the
genuine demand restricted to q <= Q_K, <= w(v)/2 + O(1): beta-good up to norm K with beta of genuine
labels, per-unit height condition holds for ||v||_1 <= K. At each greedy step >= half of the parity
class is admissible -> labels above Q_K arbitrary within that freedom (e.g. random), so no system-level
coherence is forced above Q_K.
Consequence: an inverse theorem certifying coherence up to Y from norms <= K needs Q_K >= Y, i.e.
K >~ log Y / log(2r+1). Theorem A needs Y >= C_L (8 p_j Norm zeta)^L, p_j Norm zeta <= M^(h_0+h):
K >~ L (h_0 + h) log M / log(2r+1) ~ L (h_0 + h) for r ~ M.
Replaces (in regime (1.4)) the withdrawn "must use norm >= 3 / deep 2-adic moduli".
