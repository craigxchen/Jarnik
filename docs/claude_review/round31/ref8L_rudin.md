# Referee report, block "rudin": Part VIII.2 of Astra's answer

Referee: adversarial referee for block VIII.2 (two local classes, nonunit square roots, Rudin).
Material refereed: `astra_lemma_VIII/answer.md`, lines 1778-1804 (VIII.2), with VIII.1 (lines 1748-1776) read only to
check that its hypotheses are met. Scripts (my own; Astra's `viii_check.py` was read but not imported or run) are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref8L_rudin/`.

## 0. Results

| item | statement | verdict | label |
|---|---|---|---|
| (1) | reduction to Q = q/g; two classes with unit difference at every p^f ∥ Q; diameter bound | correct | Proved + Computed |
| (2) | square-block identity (VIII.3): columns, determinant, valuations, blocks, chain, cofactors | correct | Proved + Computed |
| (3) | cofactor bound log\|t_ij\| ≤ D_ij/2 − W/4 + (ε/2)W | correct (it is an identity plus the diameter bound) | Proved |
| (4) | nonunit square-root reduction | correct | Proved + Computed (exhaustive, q ≤ 3000) |
| (5) | Rudin deduction and citation | correct, with a citation correction | Proved; citation corrected |
| (6) | consistency with constructions | consistent; nothing found at exponent ≥ 2 for B ≥ 5 | Computed + Heuristic |

**Main conclusion.** VIII.2 is correct as a conditional implication, given VIII.1. Specifically, the restriction of (L_8) to
data whose assigned blocks are pairwise coprime perfect squares, with the general nonzero integer cofactors that (L_8)
allows, implies the following three statements:
* the two-class statement at length q^{1/2+ε};
* Cilleruelo–Granville's Conjecture 9, for all a, units and nonunits alike;
* Rudin's conjecture σ(k) = O(√k).

The square-block identity, which Astra asked to have refereed first, is exact. I checked it by hand in every case and by
622,664 exact pair checks. I have not refereed the proof of VIII.1 (that is another block). I checked only that every
hypothesis of VIII.1 holds for the data VIII.2 builds.

**One correction (citation).** The text I read is arXiv math/0608109v1. In it, the implication "Conjecture 9 ⇒
Conjecture 1 (Rudin)" is **Theorem 11**, not Theorem 12. Theorem 12 there is a different implication: Conjecture 12, a
statement about lattice points on circles with a congruence condition, implies Conjecture 1. I could not open the
"final text" (Conjavi.pdf) that Astra cites, so I cannot say whether its numbering differs. Astra's Rudin argument is
exactly the proof of Theorem 11 in v1.

## 1. Definitions used

* For q ≥ 1 and each p^e ∥ q, Ω(p^e) is a prescribed set of at most two residues mod p^e. Ω(q) is the CRT set.
* X is a finite set of integers in Ω(q) lying in a real interval of length q^{1/2+ε}, with 0 < ε < 1/2. Fix x_* ∈ X,
  put g = gcd(q, {x − x_* : x ∈ X}), Q = q/g, and Y = {(x − x_*)/g : x ∈ X}.
* For an ordered distinct 8-tuple (y_0, …, y_7) of Y, with slot 0 the anchor:
  * d_i = gcd(Q, y_i − y_0), so d_0 = Q, and h_ij = gcd(Q, y_j − y_i);
  * P_i = (Q/d_i, (y_i − y_0)/d_i)^T;
  * T_p = {i ∈ [7] : y_i ≢ y_0 mod p^f} for each p^f ∥ Q;
  * n_T = ∏_{p : T_p = T} p^{2f} for nonempty T, and n_T = 1 if no prime has T_p = T.
* Measure space: one atom of mass 2f log p for each p^f ∥ Q. Then W = 2 log Q. The sign of y at that atom is its class
  mod p^f. D_ij is the mass of the atoms where y_i and y_j have different signs.

## 2. Item (1): the reduction (Proved; Computed)

**Statement (correct as written).** Let p^e ∥ q and s = v_p(g).
* If s = e, then p ∤ Q.
* If s < e, then Y occupies **exactly two** classes mod p^{e−s}, and their difference is a unit mod p.

*Proof.* Let r_1, r_2 be the occupied classes of X mod p^e. Exactly one of the following holds:
* only one class is occupied (this includes r_1 = r_2, two prescribed residues that coincide mod p^e);
* two classes are occupied, with v_p(r_1 − r_2) = σ < e.

In the first case p^e divides every x − x_*, so s = e. In the second case s = σ. The prime-to-p part g' of g is a unit
at p. Hence (x − x_*)/g ≡ 0 or (r_2 − r_1)/(p^σ g') mod p^{e−σ}, and the second value is a unit.

The proof does not depend on p being odd: p = 2 is the same computation. Coincidence of the prescribed residues modulo a
proper divisor p^σ of p^e is exactly the case s = σ, and the reduced classes are then 0 and a unit mod p^{e−σ}. Q = 1
can occur; it falls under "bounded Q".

**Diameter.** diam Y ≤ q^{1/2+ε}/g = Q^{1/2+ε} g^{−1/2+ε} ≤ Q^{1/2+ε}, since ε < 1/2. The map x ↦ (x − x_*)/g is
injective, so |Y| = |X|. Two cases are immediate:
* bounded Q gives |Y| ≤ Q^{1/2+ε} + 1;
* X = ∅ needs nothing.

**Computed** (`check_block.py`). I tested 2463 random instances: q is a product of 1 to 5 prime powers from
{2, …, 17}, with 2^e up to 2^6 and other p^e up to p^4. Residue modes were single, random pair, coalescing mod p^s
(s = 1, …, e−1) and unit difference. X is a random subset (8 to 40 marks) of Ω(q) lifted by random multiples of q.
* At every p^f ∥ Q, Y has exactly two classes, and their difference is a unit.
* 1483 primes of q dropped out of Q (only one class occupied).
* 3890 layers were coalescing.
* 1582 instances had a 2-adic layer in Q.

## 3. Item (2): the square-block identity (Proved; Computed)

All of the following hold for every ordered distinct 8-tuple, with any element as anchor.

1. **Primitivity.** gcd(Q/d_i, (y_i − y_0)/d_i) = 1 because d_i = gcd(Q, y_i − y_0). P_0 = (Q/Q, 0) = (1, 0). The
   slope of P_i is (y_i − y_0)/Q, so Φ is the image of (y_0, …, y_7) under one affine map. This is VIII.1, hypothesis 1.
2. **Determinant.** det(P_i, P_j) = Q(y_j − y_i)/(d_i d_j) = ((y_j − y_i)/h_ij) · (Q h_ij/(d_i d_j)). This is an
   identity.
3. **Valuations.** Fix p^f ∥ Q. By §2, v_p(d_i) = f if y_i is in the anchor's class and 0 otherwise, and
   v_p(h_ij) = f if y_i, y_j are in the same class and 0 otherwise. So v_p(Q h_ij/(d_i d_j)) equals:

   | position of i, j | valuation |
   |---|---|
   | both in the anchor class | f + f − f − f = 0 |
   | both opposite the anchor | f + f − 0 − 0 = 2f |
   | split | f + 0 − f − 0 = 0 |
   | i = 0, j in the anchor class | f + f − f − f = 0 |
   | i = 0, j opposite | f + 0 − f − 0 = 0 |

   Primes not dividing Q do not occur. Hence Q h_ij/(d_i d_j) = ∏_{T ∋ i,j} n_T exactly. This includes the following
   cases:
   * tuples in which all seven non-anchor marks lie in the anchor class at p (T_p = ∅, no block, every valuation 0);
   * all seven lie opposite (T_p = [7], the top block, which (L_8) allows and does not balance);
   * singleton sides (T_p = {i}, which occurs in no minor).
4. **Cofactors.** t_ij = (y_j − y_i)/h_ij is a nonzero integer, because the marks are distinct and h_ij divides
   y_j − y_i. For the anchor pairs, h_0j = d_j, so t_0j = (y_j − y_0)/d_j = Y_j, and the empty product of blocks is 1.
   This is the convention t_0j = Y_j of the expert note, and H_Σ in VIII.1 equals the H_Σ of the note.
5. **Blocks.** Each prime of Q contributes to exactly one n_T. So the n_T are pairwise coprime squares. At every prime the
   sets T with p | n_T are one set T_p (one chain plus an empty chain), so the chain condition holds trivially.
6. **Masses.** log n_T is the mass of cut T, as VIII.1 hypothesis 2 requires, and the empty pattern gets no block. The
   separating mass is D_ij = Σ_{p separates i,j} 2f log p = 2 log(Q/h_ij). This is the integer identity
   ∏_{sep} p^{2f} · h_ij² = Q².
7. **Sign functions.** VIII.1 needs one fixed sign function per mark. Each p contributes the same two classes of Y for
   every tuple (§2), so the class of y mod p^f is a fixed sign, and each tuple's cut set is read off relative to its
   anchor.

**Computed** (`check_block.py`). 22,238 ordered 8-tuples and 622,664 pairs were checked in exact arithmetic. Each of the
following was asserted for every tuple or pair:
* primitivity and P_0 = (1, 0);
* slope = affine image;
* det = t·∏ n_T, and det = t·(Q h/(d d));
* the valuation table, prime by prime;
* t ≠ 0 and t_0j = Y_j;
* pairwise coprime square blocks and the chain condition;
* the D_ij identity;
* t² h² = (y_j − y_i)².

Tuples were random, plus targeted ones (all non-anchor marks in the anchor class, or all opposite, at a chosen prime).
Coverage: 4090 tuples had a top block, 2121 a singleton side, 3792 an empty side, and 4106 had two primes sharing a
side. No assertion failed.

## 4. Item (3): the cofactor bound (Proved)

log|t_ij| = log|y_j − y_i| − log h_ij. Also D_ij/2 − W/4 + (ε/2)W = log Q − log h_ij − ½ log Q + ε log Q =
(½ + ε) log Q − log h_ij. So (VIII.1) with λ = ε/2 is **equivalent** to |y_j − y_i| ≤ Q^{1/2+ε}, which is the diameter
bound of §2. The lower bound 0 ≤ log|t_ij| holds because t_ij is a nonzero integer.

VIII.1 then applies with ε = 2λ_* (and ε < 1/2), for W = 2 log Q ≥ W_*. Bounded W is covered by the diameter. All
quantifiers are absolute: λ_*, M_*, W_* come from VIII.1, and VIII.1 fixes them before q, the residues and the interval.

## 5. Item (4): nonunit square roots (Proved; Computed)

Every step is correct. For p^e ∥ q, put v = min(e, v_p(a)).
* **v = e.** x² ≡ 0 mod p^e iff v_p(x) ≥ ⌈e/2⌉. Remove p^{⌈e/2⌉}. Then p^e ≤ p^{2⌈e/2⌉}, with strict inequality
  for odd e.
* **v < e.** If v_p(x²) ≠ v, then v_p(x² − a) = min(v_p(x²), v) < e. So a root forces v = 2t, and every root has
  valuation exactly t. Remove p^t and keep p^{e−2t}.
* Let u be the product of the removed factors and Q_1 the product of the kept factors. Then q ≤ u²Q_1, u divides
  every root, and q/u ≡ 0 mod Q_1, so x ↦ x/u mod Q_1 is well defined on root classes mod q. At p^f ∥ Q_1, write
  u = p^t u' with u' a unit at p. Then (x/u)² ≡ (a/p^{2t}) u'^{−2}, a unit square congruence.
* **Odd p.** z_1² ≡ z_2² with z_1, z_2 units, and p cannot divide both z_1 ± z_2. So z_1 ≡ ±z_2 mod p^f: at most two
  classes, with unit difference.
* **p = 2, f ≥ 2** (Astra says f ≥ 3; the argument works for f ≥ 2). One of z_1 ± z_2 has 2-adic valuation exactly 1,
  so z_1 ≡ ±z_2 mod 2^{f−1}. For f ≥ 3 there are 4 classes mod 2^f and 2 mod 2^{f−1}. For f = 1, 2 there are at most
  2 unit classes.
* Q_2 = Q_1/gcd(Q_1, 2) has at most two classes at every prime. At 2 the two classes ±z differ by 2z, which is **not**
  a unit. This is harmless, because VIII.2 allows any two prescribed residues and performs its own reduction.
* An interval of length √q becomes one of length √q/u ≤ √Q_1 ≤ √(2Q_2) < 2√Q_2. That is two intervals of length
  √Q_2 ≤ Q_2^{1/2+ε}, so VIII.2 applies at modulus Q_2. The example a = 0, q = m² (roots = multiples of m) is correct.

**Computed** (`check_sqrt.py`). The check is exhaustive over every q ≤ 400 and every a with a root. That is 26,000
classes; for q ≤ 300 it is 14,979, which agrees with Astra's count. I extended it to q ≤ 3000 (1,244,208 classes).
Every assertion above holds. Four classes at 2^f ∥ Q_1 occur exactly in the f ≥ 3 layers (30,214 of them for
q ≤ 3000), and all of them collapse to two classes in Q_2. Side datum: the largest number of roots in a closed interval
of length √q is 3 for all q ≤ 3000, first attained at q = 105, a = 46.

**Side remark (Proved, minor).** Cilleruelo–Granville's first proof of their Theorem 10 (v1, lines 572-596) does the same
factor removal. Two of its statements are slightly inaccurate:
* "(a, q) must be a square" fails when p^e | a with e odd (q = p, a = 0);
* "|Ω(q')| ≤ 2" fails at 2^f with f ≥ 3 (q = 8, a = 1 has four roots).

Astra's ⌈e/2⌉ and the passage to Q_2 repair both points.

## 6. Item (5): the Rudin deduction and the literature (Proved; citation corrected)

**What I read.** I read the text extraction `ref6L_pell/dl/cg07.txt` (arXiv math/0608109v1, 3 Aug 2006) in place, lines
1-60 and 466-800. These cover:
* the definition of Ω(n) by CRT;
* Theorems 8, 9, 10, Conjecture 9, Theorem 11, Conjectures 10, 11, 12, Theorem 12 with its proof, and the flowchart.

WebFetch of dms.umontreal.ca (Conjavi.pdf) and of arxiv.org failed: neither host resolves in this session. One
WebSearch returned only snippets. Those snippets give the venue (CRM Proceedings and Lecture Notes 43, 2007,
pp. 241-262) and the abstract. They did not show the numbering of the final text.

**What the v1 text says.**
* **Conjecture 9.** "There exists a constant N such that there are no more than N solutions
  0 < x_1 < ⋯ < x_N < x_1 + b^{1/2} to the equation x_i² ≡ a (mod b), for any given a and b." It is stated for
  **all a**, nonunits included. Astra's "O(1) roots in any interval of length √q, all a" is the same statement, up to
  a closed versus half-open interval, which changes the count by at most one.
* **Theorem 11.** "Conjecture 9 implies Conjecture 1." Its proof uses
  (x_ℓ − x_1)² ≤ (x_ℓ + x_1)(x_ℓ − x_1) ≤ (k − 1)b, and is exactly Astra's last paragraph.
* **Conjecture 10.** At most d residues per prime power; then for every ε there is N(d, ε) such that at most N(d, ε)
  integers 0 ≤ n < b^{1−ε} lie in Ω(b). The authors say they "strongly believe" the version with exponent 1/d.
* **Theorem 9** (unconditional, pigeonhole). At most k elements in any interval of length n^{α_d(k)}, with
  α_2(k) = 1/2 − O(1/k).
* **Theorem 12.** "Conjecture 12 implies Conjecture 1." Conjecture 12 says that if a_i² + b_i² = n and
  a_i² ≡ a_1² (mod q) for i ≤ m, then q = O(n^{1−δ}).

**Assessment.**
1. The VIII.2 two-class statement (d = 2, length q^{1/2+ε}) lies strictly between the two forms of Conjecture 10:
   * it is implied by Conjecture 10 for d = 2, at any exponent 1 − ε' > 1/2 + ε;
   * it is strictly stronger than the "1/d" critical form at d = 2.

   Astra's description ("a fixed power beyond the critical version"; "their full conjecture asks for lengths up to
   q^{1−ε}") is accurate. Arbitrary intervals and the range [0, b^{1−ε}) are equivalent by translating the residues.
2. The Rudin deduction is correct. In a nonconstant k-term progression with step of absolute value q:
   * distinct squares have distinct nonnegative roots;
   * all the roots solve one congruence x² ≡ a mod q;
   * their span is at most √((k−1)q);
   * 1 + ⌈√(k−1)⌉ intervals of length √q cover them (⌈√(k−1)⌉ already suffice when k ≥ 2).

   So σ(k) = O(√k).
3. **Correction.** Astra's sentence "Conjecture 9 and Theorem 12 in the authors' final text", and the brief's phrase
   "by Cilleruelo-Granville's Theorem 12", should read "Conjecture 9 and Theorem 11 (arXiv v1 numbering)". In v1,
   Theorem 12 is the lattice-point implication. If the published text renumbers, the citation must name the
   statement, not the number. I could not check the published numbering.
4. CG's Theorem 12 is itself relevant to Astra's admission that she "could not deduce Rudin from the Gaussian short-arc
   conjecture alone". CG deduce Rudin from Conjecture 12, which is a congruence-restricted lattice-point statement on
   circles, not the uniform short-arc conjecture (L_8) is about. So that gap is consistent with the literature, as far
   as I read it.

**Corrected full statement of VIII.2 (status: Proved, conditional on VIII.1 and the restricted (L_8)).** Assume (L_8) for
data whose assigned blocks n_T are pairwise coprime perfect squares (one block per prime), with arbitrary nonzero
integer cofactors t_ij. Then there are absolute B and ε ∈ (0, 1/2) such that the following hold.
* (a) For every q ≥ 1 and every choice of at most two residues mod p^e at each p^e ∥ q, every real interval of length
  q^{1/2+ε} contains at most B integers lying in the prescribed sets at every prime.
* (b) For all integers a and q ≥ 1, every interval of length √q contains O(1) roots of x² ≡ a (mod q). This is
  Cilleruelo–Granville's Conjecture 9, for all a.
* (c) σ(k) = O(√k) (Rudin's conjecture), by Cilleruelo–Granville's Theorem 11 (arXiv v1 numbering).

## 7. Item (6): constructions and consistency (Computed; Heuristic)

**Setting.** For a finite set X of B ≥ 3 integers, let Q_max(X) = ∏_p p^{e_p}, where e_p is the largest e such that X
meets at most two classes mod p^e. Every q | Q_max(X) is admissible, and Q_max is the best. Put
E(X) = log Q_max / log diam X. VIII.2(a) says: for some B_0 and ε, E(X) < 2/(1+2ε) < 2 whenever |X| > B_0 and Q_max is
large. "Approaching 2 from the wrong side" would mean exponents ≥ 2 − o(1) for growing B.

**Bounds.**
* Pigeonhole, which is CG Theorem 9 rewritten: E ≤ C(B, 2)/m(B), where
  m(B) = C(⌊B/2⌋, 2) + C(⌈B/2⌉, 2). That is 2 + 2/(B−2) for even B and 2 + 2/(B−1) for odd B.
* **Heuristic count.** Assign moduli N_π to the bipartitions π of the B points, and set Q = ∏ N_π. For fixed Q, the
  configurations with y_0 = 0 form a lattice of index Q^{B−2} in Z^{B−1}. So about L^{B−1} Q^{−(B−2)} of them lie in
  [0, L]^{B−1}. Summing over Q ≥ L^E gives about L^{B−1−E(B−3)}. So generic arithmetic gives E up to
  (B−1)/(B−3), and a single fixed Q gives only the Minkowski exponent (B−1)/(B−2). This count is a heuristic and no
  proof uses it.
* **Linear families** (Proved, elementary; it bounds only the contributions forced by the linear structure). Suppose
  X = {c·u : c ∈ C} with C ⊂ Z^r fixed and u a parameter vector. A bipartition forces a divisor of size |u| only if
  both sides lie on parallel lines. Two such directions force |C| ≤ 4, since the four intersection points of two
  pairs of parallel lines are all there is. Hence:
  * exponent 3 is possible at B = 3 ({0, a, a+b}, Q = ab(a+b));
  * exponent 2 at B = 4 ({0, u, v, u+v}, Q = uv);
  * for B ≥ 5, linear structure forces at most one large factor, and any further gain must come from coincidences of
    values.

**Results.** Exh = exhaustive over all X ⊂ [0, L] with min 0 and max L, for every L ≤ L_max (`search_exh.py`). LLL =
random bipartition designs with random prime moduli of 8 to 22 bits, lattice built by CRT, my own exact LLL, small
combinations of the reduced basis, diam ≥ 1000 (`search_lll.py`). For the LLL runs, Q_max is a lower bound using all
design primes and all primes < 200. I spot-checked it by full factorisation and found agreement.

| B | pigeonhole | heuristic (B−1)/(B−3) | best exh (L_max) | best LLL (diam) |
|---|---|---|---|---|
| 3 | 3 | — (unbounded count; true max 3) | 2.818 at L = 2000, {0, 999, 2000} | 2.938 (1.07·10^10) |
| 4 | 3 | 3 | **2.137** at L = 276, {0, 21, 136, 276}, Q = 164220 (450) | 1.956 (7.7·10^6) |
| 5 | 2.5 | 2 | 1.648 at L = 110, {0, 33, 35, 68, 110} (140) | 1.645 (3.1·10^10) |
| 6 | 2.5 | 1.667 | 1.521 at L = 89 (90) | 1.658 (7.4·10^5), symmetric set {0, 127199, 279756, 464012, 616569, 743768} |
| 7 | 2.333 | 1.5 | 1.152 at L = 35 (48) | 1.281 (5.8·10^13) |
| 8 | 2.333 | 1.4 | 1.088 at L = 31 (38) | 1.212 (9.0·10^24) |
| 9 | 2.25 | 1.333 | 1.000 at L = 30 (33) | 1.153 (3 trials) |
| 10 | 2.25 | 1.286 | 0.784 at L = 29 (30) | 1.146 (2.6·10^39) |

**Findings.**
* Exponents above 2 occur only at B ≤ 4: B = 3 tends to 3, and B = 4 reaches 2.137 at small L. B = 4 is far from B_0,
  so this does not contradict VIII.2.
* For B ≥ 5 every exponent found is below 1.66, and for B ≥ 7 below 1.29. The LLL values sit at or just above the
  Minkowski exponent (B−1)/(B−2) and below the heuristic (B−1)/(B−3).
* Nothing approaches 2 from the wrong side.
* The known structured construction (CG Theorem 8(3): x(x − 1) ≡ 0 mod p_1⋯p_k, with all 2^k elements in |x| < εn) has
  exponent about 1 with B = 2^k unbounded. It is consistent with VIII.2(a), and it is the reason CG's Conjecture 10
  needs 1 − ε rather than 1.

**Scope.** These are finite computations. They cannot test the asymptotic statement, and B = 5 is heuristically
borderline (heuristic exponent exactly 2).

## 8. Files

* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref8L_rudin/check_block.py`
  checks items (1)-(3). Output: 2463 instances, 22,238 tuples, 622,664 pairs, 0 failures.
* `.../ref8L_rudin/check_sqrt.py` checks item (4). Run as `python3 check_sqrt.py 400` and `python3 check_sqrt.py 3000`.
* `.../ref8L_rudin/search_exh.py` and `search_lll.py` produce the item (6) tables. The outputs are `exh_a.txt`,
  `exh_b.txt`, `lll_a.txt`, `lll_b.txt`, `lll_9.txt`, `lll_10.txt`.
