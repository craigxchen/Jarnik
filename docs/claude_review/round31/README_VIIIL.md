# Judge report, round 31: Part VIII of the Question B thread (Astra, "Two arithmetic consequences of Question A and a moving-place obstruction")

Source: `astra_lemma_VIII/answer.md` from line 1742 (`## Part VIII.`) to the end (line 1929), with `final8.txt`,
`viii_check.py` and the question `reply8.txt`.
Referee reports in this directory, all seven present: `ref8L_extract.md` (VIII.1), `ref8L_rudin.md` (VIII.2),
`ref8L_lenstra.md` (VIII.3), `ref8L_subspace.md` (VIII.4), `ref8L_growth.md` (VIII.5-VIII.6), `ref8L_attack.md`
(attempted refutation of (L_8) through the consequences), `ref8L_meaning.md` (novelty, literature, consequences).
Record consulted: `expert_note/expert_note.md` §§1, 2, 5; `growth/ptolemy.md` §§3(iv), 6, 7; the previous judge report
`round31/README_VIIL.md` for format; `ref6L_pell/dl/cg07.txt` (Cilleruelo-Granville, arXiv math/0608109v1) for the
citation question.
Judge scripts (my own code; I did not import Astra's or any referee's script), in `README_VIIIL/`:
* `judge8L_a.py` (outputs `out_judge8L_a_s1.txt`, `_s2.txt`, `_s3.txt`): the square-block identity (VIII.3), the divisor
  columns (VIII.5), equality of the moduli point with the configuration of the marks (all 70 cross-ratios), the chain
  condition, and the record's estimate (3.2) of `growth/ptolemy.md` §3(iv) on the same data. `judge8.py` is an earlier copy
  with a slow divisor enumeration; it was superseded and its output was not kept.
* `judge8L_b.py` (output `out_judge8L_b.txt`): the Walsh expansion and moment identities of VIII.1, a Schwartz-Zippel
  grid test, the numerical margins, the attack referee's five-point exponent-2 family, and the 16-divisor family.
* `judge8L_c.py` (output `out_judge8L_c.txt`): a spot check of the attack referee's four-point bound.

Claims of my own are labelled **[J]** and are unrefereed.

**Notation used throughout.**
* **(L_8)** is the matrix statement of `expert_note.md` §2: there are c > 0, eta_0 > 0 and a proper Zariski-closed
  Z ⊂ Mbar_{0,8} such that every datum (primitive integer columns P_0 = (1,0)^T, P_1, ..., P_7; positive integers n_T for
  every nonempty T ⊆ {1..7}; det(P_i,P_j) = t_ij ∏_{T∋i,j} n_T with t_ij nonzero integers and t_0j = Y_j; the chain
  condition at every prime) that is eta-balanced at scale w on the 119 visible blocks, with eta <= eta_0 and moduli point
  Phi ∉ Z, satisfies H_Sigma = Σ_{0<=i<j<=7} log|t_ij| >= c w - O(1).
* **(L_8)^sq** is (L_8) asserted only for data in which each prime divides exactly one block n_T and that block is an even
  power of the prime (so the blocks are pairwise coprime perfect squares), with arbitrary nonzero integer cofactors
  (cofactors may be divisible by block primes, so a block need not equal the contact integer d_T(Phi)).
* **V_8** is Vojta's Main Conjecture on Mbar_{0,8} over Q with D = 0, A = K + B, at eps = 1/43 (`growth/ptolemy.md` §6).
  Theorem D of the record derives (L_8) from V_8 with c = 6/69, eta_0 = 1/381, Z = Z_{1/43}.
* **TR(theta)**: there is B such that for every q >= 1 and every choice of at most two residues modulo each exact prime
  power p^e || q, every real interval of length q^theta contains at most B integers that lie in a prescribed residue at
  every p^e || q.
* **SQ** (Cilleruelo-Granville Conjecture 9): uniformly in a and q, every interval of length sqrt(q) contains O(1) roots
  of x^2 ≡ a (mod q).
* **RUD** (Rudin's conjecture): every nonconstant k-term arithmetic progression of integers contains O(sqrt k) squares.
* **LEN(theta)**: there is K such that for all n, q >= n^theta and r with gcd(r,q) = 1, at most K positive divisors of n are
  ≡ r (mod q). Lenstra (1984) proved LEN(alpha) for each alpha > 1/4 with K depending on alpha.
* **QS(m, delta, f, H_0)** is round 8's hypothesis (`round8/weakest.md` lines 690-691): for every finite set P of split
  primes, the solutions of the system (5.1) of `weakest.md` with height H >= H_0 (H_0 an absolute constant) lie in at most
  f(|P|) proper linear subspaces. **QS^prod** is the record's relative-threshold variant (`lit18/README.md` §2.3), with
  cutoff H >= max(H_0, (∏_{p∈P} p)^gamma), 0 < gamma < 1/2.
* **U** is the uniform Cilleruelo-Granville bound M(C) < infinity for every C.

## 0. Verdict (five lines)

1. Nothing in Part VIII is false. No referee returned FALSE on a claim of Astra; the two FALSE verdicts (meaning) are on sentences **in the record**, and the GAP verdicts are on literature status (lenstra) and on the plausibility of a hypothesis Astra does not claim (growth). I uphold all of them.
2. **"(L_8) implies Rudin's conjecture" is a theorem**, for the matrix form (L_8) of the expert note with the linear lower bound c w, anchor 0 and unsymmetrised Z, and already for its restriction (L_8)^sq. I checked myself [J] the two steps Astra asked to have prioritised (the square-block identity (VIII.3) and the escape from the fixed exceptional set in VIII.1) and that the datum is admissible for the (L_8) that the record states and that Theorem D derives from V_8. Consequently V_8 implies Rudin's conjecture; we found no earlier derivation.
3. Full (L_8), hence V_8, implies LEN(1/4 - kappa) for a fixed kappa > 0 (kappa = 10^-17 under V_8). This conflicts with nothing known: no construction of unboundedly many divisors in an invertible class at any fixed exponent is known, and the best known gives 16. No referee refuted or weakened (L_8).
4. Round 8's QS with an absolute height cutoff is **false** in every dimension (VIII.4; also on equal norms in dimension 5), now proved and refereed with the sharper bound f(s) >= exp(s^{1-delta/(m-1)-o(1)}). The growth statements VIII.5-VIII.6 are correct implications from non-standard relative-height counting hypotheses that are open, untestable, and, in the strongest case, at least as strong as U.
5. For the integer problem nothing moves: the growth rate (2/3+eps) log R/loglog R and U are unchanged, and Question B is unchanged. Question A is now known to contain Rudin's conjecture and a divisor bound past the Lenstra/Coppersmith barrier; it should no longer be presented as a realistic standalone target. This is information about difficulty, not progress.

## 1. Per block: referee-corrected statements

### 1.1 VIII.1, the extraction lemma (referee: extract; re-derived by lenstra and meaning; checked by the judge)

**Proved (statement as corrected).**
> Assume (L_8) in the matrix form of `expert_note.md` §2, with constants c > 0, eta_0 > 0, C_0 >= 0 (the O(1)) and a
> proper Zariski-closed Z ⊂ Mbar_{0,8}, anchor in slot 0 and Z not symmetrised. Fix a nonzero polynomial F in 8 variables,
> of total degree d, vanishing on the closure in A^8 of the preimage of Z ∩ M_{0,8} under the ordered-configuration map
> (x_0, ..., x_7) ↦ [x_0 : ... : x_7]. Then there are lambda_* > 0, M_* and W_*, depending only on (c, eta_0, C_0, d), with
> the following property. Let X ⊂ Q be finite with |X| = M >= 8, let (Omega, mu) be a finite measure space of total mass
> W > 0, and let each x ∈ X carry a measurable sign function b_x : Omega -> {+1, -1}; write D_xy = mu{b_x ≠ b_y}. For an
> ordered tuple and nonempty T ⊆ {1..7}, the cut T is the set of atoms where b_{x_i} ≠ b_{x_0} exactly for i ∈ T. Let
> lambda <= lambda_*. Suppose every ordered tuple of 8 distinct elements of X comes with an (L_8) datum such that
> (1) the columns are primitive, P_0 = (1,0)^T, the chain condition holds, and the moduli point Phi equals the class of
> (x_0, ..., x_7) in M_{0,8}, that is, the columns represent the marks after **one common Mobius map** (which may depend on
> the tuple); (2) log n_T = mu(cut T) for every nonempty T; (3) 0 <= log|t_ij| <= D_{x_i x_j}/2 - W/4 + lambda W for all 28
> pairs, anchor pairs included. Then W >= W_* implies M <= M_*. The same holds for a restriction of (L_8) if all data
> produced lie in the restricted class.

Corrections, all wording (extract §6): the dependence of the constants; read lambda_+ in (VIII.2); the condition
"28a/4 < c/256" is redundant; M_* >= max(4(28+d), 2.67·10^15) under the Vojta numbers.

**Which version of (L_8) is used (extract §2.4; upheld).** The matrix (L_8) with the **linear** bound c w. It does not
follow from (L_8^f) with f = o(w), because the lemma's upper bound on H_Sigma is linear in W with a small coefficient.
(L_8^pc) and (L_8^c)^cat do not suffice: both concern canonical contact integers, and the lattice form also needs pairwise
coprime cofactors of size e^{w(1±eta)}. (L_8^c)^all suffices, since it implies (L_8).

**[J] Judge check of the escape from the fixed exceptional set.** The points that carry the weight:
1. **Preimage is proper.** U = A^8 minus the big diagonal; phi: U -> M_{0,8} is the quotient by PGL_2, a surjective
   morphism. Z ∩ M_{0,8} is a proper closed subset of the irreducible M_{0,8}, so phi^{-1}(Z ∩ M_{0,8}) is a proper
   closed subset of U, and its closure in the irreducible A^8 is proper. So some nonzero F of degree d = d(Z) vanishes on
   it. d depends on Z only. A distinct tuple has Phi ∈ M_{0,8} (not on the boundary), so Phi ∈ Z forces F(x) = 0.
2. **Grid bound.** Schwartz-Zippel holds on X^8 for any finite X in a field: P(F = 0) <= d/M, by induction on the
   variables through the top coefficient. [J] Computed: for the degree-2 polynomial vanishing on the preimage of the
   harmonic cross-ratio locus of slots 0-3, the exact probability on three 30-element grids (consecutive integers,
   squares, reciprocals) is 0.0044-0.0061 against d/M = 0.067 (`out_judge8L_b.txt`).
3. **Order of quantifiers.** The balance event is computed on X^8 with repeated rows allowed (product measure); the
   hypotheses are invoked only for the tuple finally chosen, which is distinct and off Z. P(balanced and distinct and
   F ≠ 0) >= 3/4 - (28 + d)/M > 0 for M > 4(28 + d). a, lambda_*, M_*, W_* are fixed from (c, eta_0, C_0, d) before X,
   the measure space, the modulus, the parent or the grid. No exceptional set depends on q, n, S or a template. The
   standing caution of the record is respected.
4. **Common Mobius map.** In VIII.2 the slopes are (y_i - y_0)/Q and in VIII.3 they are (d_i/a - 1)/q with a = d_0: in
   each case one affine map applied to all eight marks, depending on the tuple. [J] Computed: for every tuple tested in
   `judge8L_a.py`, all 70 cross-ratios of the columns (as points Y/X of P^1) equal the cross-ratios of the marks exactly.
   A coordinatewise transformation that is not one Mobius map would break the escape; none occurs.

**[J] Judge check of the Gram/moment step.** G_xy = 1 - 2D_xy/W <= 4 lambda (integrality |t| >= 1 in (3)); E G >= 0;
E G^2 <= E|G| <= 2 E G_+ <= 2/M + 8 lambda_+; for even |A| >= 2, E_x m_A^2 = E_{u,v} g^{|A|} <= E_{u,v} g^2 = E G^2 with
g(u,v) = M^{-1} Σ_x b_x(u) b_x(v); the anchored pattern masses are 2^{-7}(1 + Σ_{A even, nonempty} ± m_A) with 127
summands; Markov gives probability >= 3/4 that all 128 masses are within aW/128 of w = W/128; exactly 64 patterns separate
any two slots, so D_ij <= (1+a)W/2 and log|t_ij| <= (a/4 + lambda)W; 28(a/4 + lambda_*) < c/256 contradicts
H_Sigma >= cW/128 - C_0 once W >= 256 C_0/c. Computed [J]: the Walsh expansion and the identity E m_A^2 = E g^{|A|} hold
exactly (Fractions) on 102,400 pattern checks over random weighted systems, with the 64-pattern bound and (VIII.2) never
violated (`out_judge8L_b.txt`). Margins [J], exact rationals: with c = 6/69, eta_0 = 1/381, a = 2·10^-5,
lambda_* = 10^-16 one has a < eta_0, 28(a/4 + lambda_*) = 1.4000·10^-4 < c/256 = 3.3967·10^-4, 127^2·8·lambda_* < (a/4)^2,
and the 2/M term needs M > 2.667·10^15.

**Computed (referee extract).** Hadamard, Paley, random, biased and clustered sign systems: the one-sided Gram bound with
E G >= 0 is what forces balance; balance at a = 0.05 is reached only at M of order 10^8. Consistent with the lemma.

**Verdict.** Correct. Gap: none. Withdrawn: none.

### 1.2 VIII.2, two residues, square roots and Rudin (referee: rudin; cross-checked by meaning and the judge)

**Proved (statement as corrected), conditional on (L_8)^sq.**
> There are absolute B and eps ∈ (0, 1/2) such that:
> (a) TR(1/2 + eps) holds with constant B (two prescribed residues per exact prime power, which may coincide modulo a
> proper divisor of p^e);
> (b) SQ holds: for all integers a and q >= 1, every interval of length sqrt(q) contains O(1) roots of x^2 ≡ a (mod q),
> nonunits included;
> (c) RUD holds: sigma(k) = O(sqrt k), by Cilleruelo-Granville **Theorem 11 (arXiv math/0608109v1 numbering)**.
> Under V_8 one may take eps = 10^-16.

**Correction (citation).** Astra cites "Conjecture 9 and Theorem 12". [J] I read `cg07.txt` lines 640-700: in v1,
Theorem 11 is "Conjecture 9 implies Conjecture 1", and Theorem 12 is "Conjecture 12 implies Conjecture 1", where
Conjecture 12 is a congruence condition a_i^2 ≡ a_1^2 (mod q) on m representations n = a_i^2 + b_i^2. The published
numbering could not be checked by anyone (no network access to the host). The citation must name the statement
("their Conjecture 9 implies Rudin"), not the number.

**[J] Judge check of the square-block identity (VIII.3).** With Q = q/g, y = (x - x_*)/g, d_i = gcd(Q, y_i - y_0),
d_0 = Q, h_ij = gcd(Q, y_j - y_i), P_i = (Q/d_i, (y_i - y_0)/d_i)^T:
* P_i is primitive and P_0 = (1, 0)^T;
* det(P_i, P_j) = Q(y_j - y_i)/(d_i d_j) = ((y_j - y_i)/h_ij)·(Q h_ij/(d_i d_j));
* at p^f || Q the two occupied classes of y differ by a unit (reduction step), so v_p(d_i) = f exactly when y_i is in the
  anchor's class and v_p(h_ij) = f exactly when y_i, y_j are in one class; hence v_p(Q h_ij/(d_i d_j)) = 2f when both i,
  j lie opposite the anchor and 0 in the other four cases (both in the anchor class; split; i = 0 with j in either class);
* so with n_T = ∏_{p : T_p = T} p^{2f} (T_p the slots opposite the anchor at p, no block when T_p is empty) one has
  det(P_i, P_j) = t_ij ∏_{T∋i,j} n_T with t_ij = (y_j - y_i)/h_ij a nonzero integer and t_0j = Y_j;
* each prime divides one block, so the chain condition holds trivially and the blocks are pairwise coprime squares;
* log n_T is the mass of cut T for the atom system "one atom of mass 2f log p per p^f || Q, sign = class of y mod p^f";
  D_ij = 2 log(Q/h_ij), W = 2 log Q;
* hypothesis (3) of VIII.1 with lambda = eps/2 is **equivalent** to |y_j - y_i| <= Q^{1/2+eps}, which is the diameter
  bound Q^{1/2+eps} g^{-1/2+eps} <= Q^{1/2+eps} (uses eps < 1/2).
Computed [J] (`judge8L_a.py`, three seeds): 1,980 ordered eight-tuples and 55,440 pairs from random two-class
configurations with coincident and coalescing residues and with powers of 2 up to 2^4, every assertion above exact, 0
failures; all 70 cross-ratios per tuple agree with those of the marks.

**[J] Admissibility for the (L_8) of the record and of Theorem D.** The datum has every feature the expert note's (L_8)
asks for (blocks for all 127 T, singletons and top block included; primitive columns; nonzero cofactors; chain condition).
Balance is proved by VIII.1 on all 128 patterns, so both the 119-visible-block reading (expert note) and the all-blocks
reading (`ptolemy.md` §7) are met without the top-block inflation argument. Theorem D's proof applies to such data: its
only arithmetic input is (3.2), |Σ_p l_T(p) log p - log n_T| <= 2 H_Sigma, whose proof (`ptolemy.md` §3(iv)) uses
primitivity, the factorisation det = t ∏ n_T and the two-chain condition, and absorbs cofactors divisible by block primes
through |l_T(p) - v_p(n_T)| <= 2 max v_p(t_ij). Computed [J]: on the same square-block data I computed the p-adic tree
lengths l_T(p) from the Gromov products v_p(det) at every prime of the data and every one of the 119 visible T; the bound
|l_T(p) - v_p(n_T)| <= 2 max v_p(t_ij) held in all 2,294,201 instances, and l_T(p) ≠ v_p(n_T) in 40,291 of them. So the
blocks are genuinely non-canonical in a positive fraction of cases, and Theorem D still covers them. The record's (L_8)
admits non-canonical blocks and cofactors sharing primes (`expert_note.md` §2 imposes no coprimality on t_ij).

**Proved (nonunit reduction; rudin item 4).** Correct as written. Minor: the 2-adic sign argument works for f >= 2, not only
f >= 3; at 2 the two classes of Q_2 differ by a nonunit, which VIII.2 allows. This also repairs two small slips in CG's
first proof of their Theorem 10 (rudin §5, side remark). Computed (rudin): exhaustive for q <= 3000 (1,244,208 classes);
the largest number of roots in an interval of length sqrt(q) is 3 (meaning reproduced this).

**Computed (rudin item 6) and reconciliation with attack.** Rudin's finite searches found best exponents
log Q_max/log diam of 1.648 (B = 5) and 1.658 (B = 6); the attack referee proves polynomial families with B = 5 and B = 6
at exponent exactly 2 in the limit. These do not conflict: the families approach 2 slowly. [J] I evaluated the attack
referee's B = 5 family {0, α(α+β), β(α+2β), α^2+2αβ+2β^2, (α+β)(α+2β)} exactly: ratio 1.9075, 1.9535, 1.9721 at
α ≈ 10^6, 10^12, 10^20. Exponent 2 at fixed B <= 6 is consistent with TR(1/2+eps), which bounds only B > B_0.

**Verdict.** Correct, conditional on VIII.1 and (L_8)^sq. Citation corrected. Gap: none.

### 1.3 VIII.3, divisors in residue classes (referee: lenstra; checked by the judge)

**Proved (statement as corrected), conditional on (L_8).**
> There are absolute lambda > 0 and B < infinity such that for every n >= 2, q >= 1, r with gcd(r, q) = 1 and real
> 0 < A < B_1 with q sqrt(A B_1)/(B_1 - A) >= n^{1/4 - lambda}, at most B positive divisors of n lie in [A, B_1] and are
> ≡ r (mod q). Consequently, for every kappa with 0 < kappa <= 2 lambda/3, at most K = B·ceil(1/kappa) positive divisors
> of n are ≡ r (mod q) whenever q >= n^{1/4 - kappa}. This is a cardinality statement and gives no enumeration algorithm.
> Under V_8, kappa = 10^-17 is admissible.

Corrections: the condition kappa < 1/8 is superfluous; the minimal condition is kappa <= 2 lambda/3 (Astra's
kappa < lambda/2 meets it). The letter B is used both for the count and inside B_1.

**[J] Judge check of the columns (VIII.5).** With a = d_0, g_i = gcd(a, d_i), P_i = (a/g_i, (d_i - a)/(q g_i))^T:
integrality uses gcd(g_i, q) = 1 (from gcd(r, q) = 1); primitivity holds since the gcd divides gcd(a/g_i, d_i/g_i) = 1;
det(P_i, P_j) = t_ij a h_ij/(g_i g_j) with t_ij = (d_j - d_i)/(q h_ij) a nonzero integer; at p the block part has valuation
b_0 + min(b_i, b_j) - min(b_0, b_i) - min(b_0, b_j), the number of layers at which both i, j are opposite slot 0; the
cuts at layers s <= b_0 are {i : b_i < s} (increasing) and at s > b_0 are {i : b_i >= s} (decreasing), two chains inside
the disjoint sets {b_i < b_0} and {b_i > b_0}. Computed [J] (`judge8L_a.py`, three seeds):
2,708 ordered eight-tuples, 75,824 pairs, with n having 2 to 5 prime factors of exponent up to 4 and, in 30% of instances, a
prime dividing q (non-invertible divisors excluded): divisibility, primitivity, nonzero cofactors, t_0j = Y_j, the chain
condition by brute-force partition search, and the 70 cross-ratios, all exact, 0 failures; (3.2) held in all 3,414,943
instances, with l_T(p) ≠ v_p(n_T) in 55,945.

**Literature status (lenstra item 4; GAP upheld).** From search snippets only (all full-text fetches failed): Lenstra
1984 proves LEN(alpha) for alpha > 1/4 and at most 11 divisors for q > n^{1/3}; Coppersmith-Howgrave-Graham-Nagaraj 2008
give O((alpha - 1/4)^{-3/2}) constructively; unboundedness is known only for non-invertible classes (Lenstra's divisors of
n^2 that are ≡ n mod 2n). No source seen proves boundedness at any alpha <= 1/4 or constructs unboundedly many divisors in an
invertible class at any fixed alpha > 0. The statement "LEN(1/4 - kappa) is open" is therefore a **literature statement
from snippets, not confirmed against full texts**.

**Proved (lenstra, attack; reproduced [J]).** For primes or merely pairwise coprime factors 1 + iq with t! | q,
n = ∏_{i<=t}(1 + iq) has 2^t divisors ≡ 1 (mod q) with log q/log n -> 1/t from below. [J] At t = 4, q = 24·10^3 and
24·10^8: 16 divisors, alpha = 0.2317 and 0.2411. So K >= 16 at 1/4 - kappa, which is bounded. Computed (lenstra):
exhaustive n <= 3·10^6, records 20 at alpha >= 0.20, 11 at 0.25, 6 at 1/3; no growth.

**Verdict.** Correct, conditional on (L_8). Literature status open from snippets.

### 1.4 VIII.4, the absolute-threshold Subspace hypothesis (referee: subspace)

**Proved (statement as corrected).**
> Fix k >= 2, 0 < delta < k - 1, H_0 >= 1, 0 < theta < 1 - delta/(k-1) and 0 < c < 1 - theta - delta/(k-1). Let P_X be the
> primes p ≡ 1 (mod 4) in [X, 2X] and s = |P_X|. For all X >= X_0(k, delta, H_0, theta, c), any family of proper linear
> subspaces of Q(i)^k covering all solutions of round 8's system (5.1) with m = k, P = P_X and H >= H_0 has at least
> exp(c s^theta log s) members. Since s(X) moves by at most 1, QS(m, delta, f, H_0) forces f(s) >= exp(c s^theta log s)
> for all large s. So QS is false for every m >= 2, 0 < delta < m - 1, H_0, and every f = exp(o(s^{1-delta/(m-1)})).
> The same holds for thresholds H_0(s) <= exp(2b(1 - o(1)) s^theta log s), b = delta/(k-1).

> **Equal norms (k = 5).** For 0 < delta < 2 and 0 < theta < 1 - delta/2, the same lower bound holds for covers of the
> solutions whose five coordinates have exactly equal complex modulus.

> **Referee extension (subspace §5; Proved, not judge-checked).** The equal-norm obstruction holds throughout
> 0 < delta < (k-1)/2 for every 2 <= k <= 6; for k >= 7 the method covers only delta < 2(k-1)/5.

Corrections:
1. Round 8 proposed f = o(s), not a polynomial f. Refuting every sub-exponential f is stronger, so the conclusion stands;
   Astra's phrase "polynomial bound proposed in round 8" should read "round 8's QS (which asked for f = o(s)) and every
   sub-exponential weakening".
2. The refutation is not new: lit18 Prop A (unrefereed) already gave exp(s^{1/2-o(1)}). Astra's exponent
   s^{1-delta/(m-1)-o(1)} is at least as good everywhere and strictly better for delta < (m-1)/2; her proof is
   elementary, and it is now refereed. lit18 Prop B (thresholds H_0(s), linear count) stays unrefereed.
3. The record's statement that an escaping S-integer input "must use equality of norms, or a threshold tied to ∏p" is
   too weak: equal norms alone do not escape at a fixed threshold. What escapes the known constructions is full support or
   a threshold tied to ∏_{p∈P} p.

Scope (Proved): both constructions have log H = o(s log X) < log ∏_{p∈P} p, so QS^prod is untouched. Computed (subspace):
exact place-by-place checks at X = 1000, s = 67.

**Verdict.** Correct. Novelty: sharper exponent and equal-norm version; the bare refutation was already in the record.

### 1.5 VIII.5-VIII.6, conditional growth (referee: growth)

**Proved (incidence bound).**
> Let X ⊂ Z[i] be M distinct marks of norm N on an arc of angular width at most C N^{-1/4}, k >= 2, M >= 2k. A ratio
> u ≠ 1 is realised by at most 1 + C^2 ordered pairs of X; a proper linear subspace of Q(i)^k contains at most
> max(2, 1 + C^2) M^{k-2} ordered distinct k-tuples of X; so the least number of proper subspaces covering them lies
> between c_{k,C} M^2 and M(M - 1).

**Proved (conditional growth).**
> Fix k >= 2, 0 < delta < (k-1)/2, 0 < gamma < 1/2, H_0 > 0 and f. If for every finite set P of split primes the solutions
> of round 8's (5.1) with H >= H_0 (∏_{p∈P} p)^gamma lie in at most f(|P|) proper subspaces, then every critical arc has
> M <= max(2k, c_{k,C}^{-1/2} f(s)^{1/2}) once N' >= N_0(C, k, delta, gamma, H_0), with s = omega(N') <=
> (2 + o(1)) log R/loglog R. Hence f(s) = o(s^2) gives M = o(log R/loglog R), and f(s) << (1 + s)^{2-xi} gives
> M << (1 + log R/loglog R)^{1 - xi/2}.

**Proved (VIII.6).**
> Fix 0 < gamma < 1/2, an integer l >= 1 and 0 < delta < 1/(2l). Let F_delta(P) count u ∈ Q(i) with |u| = 1, u ≠ 1, a unit
> at every finite place outside P, H(u:1) >= H_0 (∏_{p∈P} p)^gamma and |u - 1|^2 <= H(u:1)^{-delta}. If
> F_delta(P) << (1 + |P|)^A, then critical arcs satisfy M << 1 + (log R/loglog R)^{A/(2l)}. If one A works for all small
> fixed delta, then M << (log R/loglog R)^tau for every tau > 0. Up to four unit choices, F_delta(P) counts b = x + iy with
> gcd(b, conj b) = 1, all primes of h = x^2 + y^2 in P, h >= H_0 (∏ p)^gamma and 0 < |eps b - conj b|^2 <= h^{1-delta}; the
> four left sides are 4y^2, 4x^2, 2(x+y)^2, 2(x-y)^2.

**Gap (growth item 2'; not a gap in Astra's proof).** The proof uses exactly gamma < 1/2; whether gamma >= 1/2 is
excluded is undetermined.

**Heuristic + Computed (growth item 5; GAP on plausibility upheld).** F_delta(P) << |P|^A is not refuted. A random-angle
model predicts F_delta(P) = O_{delta,gamma}(1), but only with constants of size about exp(c exp(C/(gamma delta))), and
exact counts for the first s split primes (s = 3..9) double with each prime, within 5% of the model. No finite computation
can test the hypothesis. It is not implied by abc (truncated Vojta on P^1 gives only rad h >= h^{delta/2 - eps}), by U, or by
any conjecture the referee found. With A = 0 it implies U through VIII.6, so it is not weaker than the target.
QS^prod(2, delta, f) implies F_delta(P) <= f(|P|).

**Not checked (growth).** Whether unions of unrelated clusters refute the subquadratic hypothesis of VIII.5: heuristically
no (cross-cluster products need r <= (k-1)/(2 gamma delta) clusters).

**Verdict.** All proofs correct. Novelty: the M^{k-2} incidence bound (QS^prod now needs f = o(s^2) instead of o(s)) and
the VIII.6 torus count with l-fold amplification. The hypotheses are non-standard and untestable.

### 1.6 Attempted refutation of (L_8) (referee: attack)

**NOT REFUTED.** No exact family meets the hypotheses of VIII.1, or the conclusions TR(1/2+eps) or LEN(1/4-kappa) in the
failing direction, with unboundedly many points.

**Proved (attack).**
* CG Conjecture 10 with d = 2 implies TR(1/2 + eps) for every eps < 1/2. So refuting (L_8) through TR first requires
  refuting that published conjecture. [J] Wording: attack's sentence that CG's "1/d" belief is "of the same strength as"
  TR(1/2+eps) is inaccurate; TR(1/2+eps) implies TR(1/2), and no converse is known (rudin and meaning say so correctly).
* With R(B) the supremum of log q_2/log L over B-point two-class configurations, (L_8) implies limsup_B R(B) < 2 (under V_8,
  R(B) < 2 - 4·10^-16 for large B). **Correction to the referee brief (FALSE upheld):** a rate R(B) = 2 - c/B would refute
  TR(1/2 + eps), hence (L_8); it is not "consistent with (L_8)".
* A weighted Prouhet-Tarry-Escott dictionary produces two-class families with exponent k/(k-1-r); proved families give
  R(B) >= 2 for B <= 6 and R(B) >= 1 + (2 + o(1))/log_2 B in general; the best exponents decrease towards 1.
* For four integers, q_2 <= sqrt(2) L^{5/2}. [J] Spot check: no violation on 297,639 four-point sets (exhaustive for
  L < 120), maximum ratio 1.877 (`out_judge8L_c.txt`).
* In characteristic 2, affine F_2-hyperplanes of F_{2^m} give 2^{m+1} points at exponent exactly 2 (Sylvester-Hadamard);
  the constant field grows. Not judge-checked.

**Computed (attack).** Record families have at most 10 marks, or Gram entries tending to 1, or no integer realisation.
Best integer trade-off found: M = C(k, k/2) at lambda = (k-4)/(4k), so M = 6 as lambda -> 0.

### 1.7 Novelty, literature, consequences (referee: meaning)

* **Proved (by reading the record).** "(L_8) implies RUD" and "(L_8) implies LEN(1/4 - kappa)" are new to the record.
* **FALSE (record sentences; upheld [J]).** `overview/overview.md` line 179-180 ("arcs => Rudin-type statements, not the
  reverse"), `lit18/direct.md` line 389 ("CG prove only the direction arcs => squares-in-AP and Rudin-type statements")
  and `local_smoothing/chains.md` row F ("arc / short-interval boundedness conjectures imply Rudin-type conjectures")
  misdescribe CG. [J] In `cg07.txt`, Rudin follows from Conjecture 9 (Theorem 11) and from Conjecture 12 (Theorem 12);
  both are congruence statements, and no arc conjecture is linked to Rudin.
* **Literature (bounded search; arxiv.org unreachable).** No earlier derivation of RUD, SQ or TR(theta >= 1/2) from Vojta,
  Bombieri-Lang or abc was found. Known conditional bound: Bombieri-Lang with Caporaso-Harris-Mazur gives sigma(k) << k^{5/7}
  (CG, read); unconditional best k^{3/5+o(1)} (Bombieri-Zannier). "V_8 implies Rudin" should be stated as "we found no
  earlier derivation".
* **Computed + Proved (meaning R6, not judge-checked).** The Liouville exponents of the two-residue problem (CG Theorem
  10) and of the arc problem (Cilleruelo-Cordoba) coincide for every k (exact for k <= 59).
* **Open.** U => RUD; RUD, SQ or TR(theta) => U; Question B => RUD; (L_8^pc) => RUD; TR(1/2+eps) => CG Conjecture 12.
* **Assessment (upheld).** Part VIII is no evidence against (L_8): every consequence is believed true. It shows that
  Question A is at least as hard as Rudin's conjecture.

## 2. What is new relative to the record and to the literature, and correct

New to the record and correct:
1. VIII.1, a probabilistic extraction lemma turning a sign system with one-sided Gram bound into a balanced (L_8) datum off
   any fixed Z (Proved). It is the shared engine of VIII.2 and VIII.3; its mechanism (Walsh expansion plus E G >= 0) is in
   substance item 532 of the record made abstract.
2. **(L_8)^sq implies TR(1/2 + eps), SQ for all a, and RUD** (Proved, conditional). With Theorem D: **V_8 implies Rudin's
   conjecture.** New to the record; no earlier derivation found in the literature (bounded search). The step
   TR(1/2 + eps) => SQ => RUD is Cilleruelo-Granville's (their §4 remark, the first proof of their Theorem 10, Theorem
   11); Astra's new step is (L_8)^sq => TR(1/2 + eps).
3. **(L_8) implies LEN(1/4 - kappa)** and the interval version (VIII.4), conditional (Proved). New to the record; V_8 thus
   implies a divisor bound past the Lenstra/Coppersmith-Howgrave-Graham-Nagaraj barrier. Not found in the literature
   (snippets only).
4. QS at a fixed threshold is false with f(s) >= exp(s^{1-delta/(m-1)-o(1)}), and in dimension 5 (referee: 2 <= k <= 6)
   even on equal norms (Proved). The bare refutation was lit18 Prop A (exponent 1/2, unrefereed); the sharper exponent and
   the equal-norm version are new, and the fixed-threshold refutation is now refereed.
5. The arc incidence bound M^{k-2} per hyperplane, which lowers QS^prod's requirement from o(s) to o(s^2) (Proved).
6. The VIII.6 reduction to a one-variable count of norm-one P-units near 1, with l-fold amplification (Proved).
7. Referee results: the Prouhet-Tarry-Escott dictionary and exponent-2 two-class families at B = 5, 6; q_2 <= sqrt2 L^{5/2}
   for B = 4; the characteristic-2 Hadamard family; the 2^t-divisor family (K >= 16 at 1/4 - kappa); the coincidence of
   Liouville exponents of the two-residue and arc problems; the correction of the record's "CG: arcs => Rudin".

Not new: TR => SQ => RUD (Cilleruelo-Granville); the bare refutation of fixed-threshold QS (lit18 Prop A); the relative
cutoff hypothesis itself (the record's QS^prod); the affine-independence lemma used in VIII.6 (parallel thread's
Theorem 1.2 mechanism).

## 3. Rulings

**(a) "(L_8) implies Rudin's conjecture" is a THEOREM**, for exactly this version:
> (L_8) in the matrix form of `expert_note.md` §2, with the linear lower bound H_Sigma >= c w - O(1), anchor 0 and
> unsymmetrised Z, asserted for data whose blocks n_T are pairwise coprime perfect squares (one block per prime, assigned
> blocks not required to equal the contact integers) and whose cofactors are arbitrary nonzero integers, implies that
> every nonconstant k-term integer arithmetic progression contains O(sqrt k) squares.

It is **not** shown from (L_8^f) (sublinear f), from (L_8^c)^cat, or from the lattice form (L_8^pc) of the Question A box;
it does follow from (L_8^c)^all and from V_8 (Theorem D applies to these data; checked [J], §1.2). Every constant is fixed
before q, the residues and the interval; the exceptional set is the fixed Z of (L_8). The implication is one-directional:
Question A is not equivalent to Rudin or to a two-residue problem. The constants are ineffective through Z and C_0.

**(b) The divisor statement below 1/4 is a THEOREM conditional on full (L_8)** (and on V_8 with kappa = 10^-17). It
conflicts with nothing known. Known: LEN(alpha) for every alpha > 1/4 (Lenstra; CHN); unboundedness only in non-invertible
classes; the best invertible-class construction near 1/4 gives 16 divisors; exhaustive data to 3·10^6 show no growth;
random models predict boundedness for every fixed alpha > 0. Whether LEN(1/4) or LEN(1/4 - kappa) holds is open as far as
the referees could establish, from snippets only (GAP on literature confirmation). The statement uses gcd(r, q) = 1
essentially, and it supplies no algorithm.

**(c) No referee refuted or weakened (L_8).** The attack referee found every known integer family consistent with it and
with CG Conjecture 10. The consequences (TR(1/2+eps), SQ, RUD, LEN(1/4-kappa)) are all believed true, so Part VIII adds no
evidence against (L_8) and none for it beyond consistency. It does give two new falsification tests: an integer family of
two-class configurations with log q_2/log L >= 2 - o(1) and unboundedly many points, or unboundedly many divisors in one
invertible class at q >= n^{1/4 - kappa} for every kappa, would refute (L_8) and V_8 (in the form Theorem D uses). Under the
Vojta margins (10^-16) neither test is practical.

**(d) Round 8's hypothesis QS(m, delta, f, H_0) is FALSE** for every m >= 2, 0 < delta < m - 1 and absolute H_0, with
f(s) >= exp(s^{1-delta/(m-1)-o(1)}) forced; round 8's Theorem 5.2 remains a correct implication from a false hypothesis.
Nothing unconditional depended on it. The relative-threshold QS^prod is untouched and stays open and non-standard; the
record's "do not spend proof effort on it" stands.

**(e) Conditional growth statements: correct implications; hypotheses non-standard, open and untested.**
* VIII.5: QS^prod with f(s) = o(s^2) gives M = o(log R/loglog R). The hypothesis is the record's QS^prod with a weaker
  required bound. Not refuted; no standard conjecture implies it.
* VIII.6: F_delta(P) << |P|^A with A independent of delta gives M << (log R)^tau for every tau. Not refuted; heuristically
  true even with A = 0, with constants beyond any computation; with A = 0 it implies U. It is not "weaker standard" input.
* None is false, none is plausible in a testable sense, none is derived from abc or Bombieri-Lang. The record's sentence
  "no intermediate growth bound is known from a weaker standard conjecture" **stands as written**.

**(f) Consequences for the integer problem. Blunt: nothing moves.**
* Growth rate (2/3 + eps) log R/loglog R (Theorem 1.1): **unchanged**.
* Uniform bound U: **unchanged**. U is not known to imply Rudin (meaning §5: the sum-of-two-squares route puts both
  coordinates in one interval, the unconditional case). The hardness is a property of (L_8), not of U.
* Question A: **status open, now known to be at least as hard as Rudin's conjecture** (open since 1960, best exponent 3/5)
  even in the restriction (L_8)^sq, and, in full, as hard as a divisor bound past the 1/4 barrier where lattice methods stop.
  It is not a realistic standalone target. Its value is as a precise Vojta-predicted statement and a falsification object.
  The Question A box currently invites a proof of "(L_8), or at least its lattice form (L_8^pc)"; the hardness
  result applies to (L_8)^sq, not to (L_8^pc), so (L_8^pc) is not shown to be Rudin-hard.
* Question B ((L_5°)): **unchanged**. Part VIII proves nothing about it; it is not known to imply Rudin, and that is no
  evidence that it is easier.

## 4. Recommended edits to the expert note

1. **§2, new paragraph after "Theorems around (L_8)".** Insert:
   > **Strength of $(L_8)$** (round 31, Part VIII of the Question B thread; **proved** conditional implications, refereed,
   > `round31/README_VIIIL.md`). $(L_8)$ implies Rudin's conjecture: every nonconstant $k$-term integer progression
   > contains $O(\sqrt k)$ squares. The proof uses $(L_8)$, with its linear bound $c\,w$, only on data whose blocks are
   > pairwise coprime perfect squares, one block per prime, with arbitrary nonzero cofactors; it does not use $(L_8^f)$,
   > $(L_8^c)^{\rm cat}$ or the lattice form $(L_8^{pc})$. It passes through three statements. (i) If at most two residues
   > are prescribed modulo each exact prime power $p^e\parallel q$, any interval of length $q^{1/2+\epsilon}$ contains
   > $O(1)$ integers in the prescribed classes, for one absolute $\epsilon>0$ (the case $d=2$ of Cilleruelo–Granville's
   > Conjecture 10 at a fixed exponent above their critical $1/2$). (ii) Hence any interval of length $\sqrt q$ contains
   > $O(1)$ roots of $x^2\equiv a\pmod q$, for every $a$ (their Conjecture 9). (iii) Hence Rudin's conjecture (their
   > implication "Conjecture 9 $\Rightarrow$ Conjecture 1", Theorem 11 in arXiv math/0608109v1). Full $(L_8)$ also implies
   > that an invertible residue class modulo $q\ge n^{1/4-\kappa}$ contains $O(1)$ divisors of $n$, for one absolute
   > $\kappa>0$: past the exponent $1/4$ of Lenstra and Coppersmith–Howgrave-Graham–Nagaraj (a cardinality bound, no
   > algorithm). Under $V_8$ one may take $\epsilon=10^{-16}$, $\kappa=10^{-17}$, so Vojta's conjecture for the single
   > variety $\overline M_{0,8}$ implies Rudin's conjecture; we found no earlier derivation (Bombieri–Lang gives
   > $\sigma(k)\ll k^{5/7}$, Cilleruelo–Granville). None of these consequences is known to fail; they measure the strength
   > of $(L_8)$, not its plausibility. The uniform bound itself is not known to imply Rudin's conjecture.
2. **Question A box (line 241).** Append: "Any proof of $(L_8)$, even restricted to pairwise coprime square blocks, proves
   Rudin's conjecture; full $(L_8)$ also proves a divisor bound below Lenstra's exponent $1/4$. The lattice form
   $(L_8^{pc})$ is not known to imply either."
3. **§1, line 72.** After "no intermediate growth bound is known from a weaker standard conjecture." insert: "In the other
   direction, the sharpest Vojta-type form of the problem, $(L_8)$, implies Rudin's conjecture on squares in arithmetic
   progressions (§2); the uniform bound itself is not known to imply any recognised open problem."
4. **§5 item 5 (line 722-723).** Replace "These give at best $M\le\omega(N')$ (lit18 Props A, B; **unrefereed**)." with:
   "With an absolute height threshold, polynomial or even sub-exponential exceptional-subspace counts are false in every
   dimension: covering needs $\ge\exp(s^{1-\delta/(m-1)-o(1)})$ subspaces, also on exactly equal norms in dimensions
   $2\le k\le6$ (round 31, Part VIII.4 of the Question B thread; **proved**, refereed; earlier lit18 Prop A,
   $\exp(s^{1/2-o(1)})$). With thresholds depending on $|S|$ only, lit18 Prop B (**unrefereed**) still caps them at
   $M\le\omega(N')$. A threshold relative to $\prod_{p\in P}p$ (QS$^{\rm prod}$) is open and non-standard; on one critical
   arc it would suffice with $f(s)=o(s^2)$ subspaces (Part VIII.5), and a polynomial count of norm-one $P$-units near $1$
   with exponent independent of $\delta$ would give $M\ll(\log R)^\tau$ for every $\tau>0$ (Part VIII.6); with exponent
   $0$ that count already implies the uniform bound."
5. **§5 item 12.** Add the bullet:
   > * *Hardness of Question A* (Part VIII of the Question B thread; refereed, `round31/README_VIIIL.md`). $(L_8)$
   >   restricted to pairwise coprime square blocks implies Cilleruelo–Granville's two-residue statement at length
   >   $q^{1/2+\epsilon}$, their Conjecture 9 for all $a$, and Rudin's conjecture; full $(L_8)$ implies $O(1)$ divisors of
   >   $n$ in an invertible class modulo $q\ge n^{1/4-\kappa}$. New to the record and not found in the literature
   >   (bounded search). No referee refuted $(L_8)$; the best integer two-class families reach exponent $2$ only for at most
   >   six points, and the best divisor family gives $16$ divisors near $1/4$. Round 8's QS with an absolute height
   >   threshold is false in every dimension. The growth rate, the uniform bound and Question B are unchanged.
   Also extend the judge list and scores of item 12 with `README_VIIIL.md` and 2/10.
6. **Elsewhere in the record (not the note).** In `overview/overview.md` line 179-180, `lit18/direct.md` line 389 and
   `local_smoothing/chains.md` row F, replace "arcs ⇒ Rudin-type statements" by "two congruence conjectures (their
   Conjectures 9 and 12) ⇒ Rudin (Theorems 11 and 12 of arXiv v1); their arc conjectures 13–15 are linked only to an
   $L^4$ conjecture for squares (Theorem 16)". In `round8/README.md` lines 331 and 392 and `round8/weakest.md` lines 37,
   725-726, 759-765 and 948, mark QS as false for every m at fixed H_0. In `overview/overview.md` line 1730,
   `lit18/README.md` line 135 and `lit18/check_tools.md` lines 40 and 321, replace "equality of norms, or a threshold tied
   to ∏p" by "full support, or a threshold tied to ∏p".

Sentences the note must not contain: "the uniform bound implies Rudin"; "Cilleruelo–Granville show that arcs imply Rudin";
"Vojta's conjecture is known to imply Rudin" (say "we found no earlier derivation"); "Question A is equivalent to Rudin";
"(L_8^{pc}) implies Rudin"; "Part VIII is evidence against Vojta"; "Part VIII gives growth from a weaker hypothesis".

## 5. Open sub-questions, ranked by honest probability of progress

1. **Bookkeeping** (P ≈ 0.95). Record §4 and correct the three "arcs ⇒ Rudin" sentences. *First step:* the edits above.
2. **Confirm the literature status of LEN(1/4) and of "V_8 ⇒ Rudin"** (P ≈ 0.7 of a clear answer, once full texts are
   reachable). *First step:* with network access to arxiv.org, ams.org and pub.math.leidenuniv.nl, read the closing remarks
   of Lenstra 1984, §§1 and 6 of Coppersmith-Howgrave-Graham-Nagaraj 2008, the open-questions section of Hales
   arXiv:2410.05030, and the published Cilleruelo-Granville text for its theorem numbering.
3. **Does TR(1/2 + eps) imply CG Conjecture 12?** (P ≈ 0.4; calibration only.) *First step:* write the coordinates
   a_i <= sqrt(n) with q >= n^{1-delta} as roots of a_i^2 ≡ a_1^2 (mod q) and handle nonunit a_1 by the VIII.2 reduction.
4. **Is (L_8^pc) Rudin-hard?** (P ≈ 0.25.) This decides whether the Question A box's fallback target is also hard.
   *First step:* test whether the square-block data of VIII.2 can be modified so that cofactors are pairwise coprime and of
   size e^{w(1±eta)} and blocks equal contact integers; the obstacle is that VIII.2's cofactors are small and share primes.
5. **Sharpness of gamma < 1/2 in VIII.5** (P ≈ 0.3; low value). *First step:* search for arcs of arc-gcd 1 containing a
   pair of height O(N'^{1/2}) as N' grows.
6. **The characteristic-0 exponent-2 question in the node model** (P ≈ 0.2): is the number of 0/1 patterns with k/2 - 1
   vanishing weighted moments bounded uniformly in k? A bounded answer is a structured analogue of TR(1/2); growth would
   threaten (L_8). *First step:* a full search for r = 3, 4 beyond two-block seeding.
7. **A growth bound from a standard conjecture** (P ≈ 0.03). Part VIII confirms that the relative-height counts needed are
   non-standard and at least as strong as U in the extreme case. *First step:* none worth spending effort on.
8. **Question A or Question B through this thread** (P ≈ 0). Proving Question A now requires proving Rudin's conjecture.

**Recommendation.** Stop the Question A hardness thread after items 1-2: its answer is in, and it is negative for
feasibility. Do not ask the colleague to attack (L_8) further. If the exchange continues, the one question worth asking
the colleague is the converse direction for the target that matters:
> Is Question B ((L_5°)) or the uniform bound itself known or provable to imply any recognised open problem (Rudin, CG
> Conjecture 12, a two-residue statement, or a uniform S-unit count)? Or, conversely, can (L_8^pc), the lattice form
> that every proof of (L_8) must cover, be shown Rudin-hard?
Otherwise, close the thread.

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| VIII.1 extraction lemma | 1.5 | correct and uniform, escape from Z checked [J]; a tool for hardness implications, no bound |
| VIII.2 two residues, square roots, Rudin | 2.5 | correct, square-block identity checked [J]; first recognised hard consequence of Question A; information about difficulty only |
| VIII.3 divisors below 1/4 | 2 | correct, columns checked [J]; past the Lenstra barrier, conditional; literature status from snippets |
| VIII.4 fixed-threshold QS false | 1 | correct, sharper than lit18 Prop A, equal-norm version new; closes a hypothesis nobody relied on |
| VIII.5-VIII.6 conditional growth | 1 | correct implications; hypotheses non-standard and untestable, at least as strong as U at A = 0 |
| attack (referee) | 1 | no refutation; useful corrected test criterion and families; consistent with everything |
| meaning (referee) | 1 | correct novelty map; corrects a false sentence in the record |
| **Part VIII overall** | **2** | sound, nothing false, one citation fix; the first hardness theorem for Question A; no movement on the growth rate, U or Question B |

Literature: only `cg07.txt` (Cilleruelo-Granville arXiv v1) was read in full text, by the rudin, attack and meaning referees and,
for the theorem numbering, by me. Everything said about Lenstra, Coppersmith-Howgrave-Graham-Nagaraj, Hales,
Shkredov-Solymosi and the published Cilleruelo-Granville text rests on search snippets. My own checks are the escape from Z
and the Gram/moment step (§1.1), the square-block identity and its admissibility for Theorem D (§1.2), the divisor columns
(§1.3), the reconciliation of the five-point exponent (§1.2), and the four-point spot check (§1.6).
