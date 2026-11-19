# Referee report, round 31, block "square"

Scope: `astra_lemma/answer.md` sections 3, 4 and 6. That covers the relations modulo squares (4), (5), (6), the
reduction to indices 1, 3, 5, the Zariski-dense parity family, the squarefree bounds (7), the ten-pair-cut angle
model, the fifteen-cut index-one profiles, identity (8) and the "exact loss" paragraph. Sections 2 and 5 belong to
other blocks. Where this block uses section 2, I say so and treat that part as conditional.

My scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_square/`:
* `gauss_util.py`: Gaussian arithmetic, profiles, lattice index (Smith form), rank mod 2 including the unit bit.
* `combin.py`: the combinatorial lemmas, exact LPs over pair distributions, the Smith form of the pair matrix, and
  the M(c) census.
* `angle_model.py`: an independent sampling of the ten-pair-cut angle model with high precision, checking both gap
  families.
* `identity8_family.py`: identity (8) on 200 actual tuples with multi-prime blocks; exact checks of the parity
  family; the synthetic index distribution.
* `windows.py`: actual circles. Consecutive 5-point windows, gcd-normalised. Checks (1) on every pair and (4) on
  every squareclass-dependent window.

I wrote all of them myself and did not reuse Astra's code. I also read Astra's `check_profiles.py`, and its
assertions agree with mine.

## 0. Summary verdict

Every mathematical statement in the block is correct. Two need a small repair of wording or hypothesis, and none
is false. Read against the record, nothing in the block goes past Liouville. (4), (5), (6) and (7) are
consequences of the monomial gaps of `round4/outside.md` Lemma 0.1: (5) is the gap applied to the half-character
`u`, where `s = t^2` and `t` is the monomial of `u`. So these bounds are fake-satisfied and fall under Lemma K1
(`round20/kill_gauss.md` §6.1). They improve the exponent only on **thin** classes:
* squareclass-dependent profiles (exponent 4, or 5 when balanced);
* squarefree norms with at most 9 or at most 6 prime factors (exponent 7/2 or 4);
* lattices of index at least 7 (through section 2).

The generic class, index 1, is untouched. Astra says this itself. The ten-pair-cut angle model is correct. It is a
routine instance of the record's fake constructions, placed on the Cilleruelo–Córdoba-extremal profile. **None of
this changes the growth rate, the uniform bound, or the status of Questions A or B or route 5.** The block is
correct and peripheral.

## 1. Relations modulo squares, (4) and (5)

**(1) and the structure of U(Q).** If `s ∈ Q(i)` has `|s| = 1`, the norm-one condition forces `k_p + l_p = 0` at
each split prime, and no 2-adic or inert part survives. So `U(Q) = μ4 × ⊕ Z g_p` with `g_p = π_p/π̄_p`, and this
group is free modulo μ4. Then `D(s) = |A|^2` and `|s-1| ≥ D^{-1/2}` because `εA - Ā` is a nonzero Gaussian integer.
Correct.

**(5).** Let `s = t^2 ≠ 1` with `t ∈ U(Q)`. Then `|t-1|^2 + |t+1|^2 = 4`, so `max(|t-1|, |t+1|) ≥ √2`. Each factor
is at least `D(t)^{-1/2}`, using (1) for `t` and for `-t`, with `D(-t) = D(t)`. Since `D(t^2) = D(t)^2`, this gives
`|s-1| ≥ √2 D(s)^{-1/4}`. Correct. In my 20000 random exact-to-60-digit trials the minimum ratio was 2.83. That only
confirms the sign, since random squares are far from 1.

**(4).**
* *Even support.* Adding the anchor to an odd support is legitimate because `y_5 = 1`. Any two odd-on-support
  integer lifts with sum zero differ by an even vector, so every such word is a square. Correct.
* *Two marks.* `D(y_i/y_j) ≤ H` gives `Δ ≥ √2 H^{-1/4}`, which is stronger than stated.
* *Four marks.*
  * All three matchings `(ij|kl)` are admissible lifts, and `|s-1| ≤ 2Δ` by telescoping.
  * The averaging step needs `Σ_matchings |a_i+a_j-a_k-a_l| ≤ 3e_p` on `[0,e_p]^4`. This holds because the sum is
    convex in each variable separately, so its maximum is at a vertex. The vertex values are `0,3,2,3,0`. I checked
    this by hand and exhaustively for `e_p ≤ 7` (`combin.py`).
  * So some matching has `log D(s) ≤ W`, and if `s ≠ 1` then (5) gives `Δ ≥ 2^{-1/2}H^{-1/4}`.
* *A word equal to 1.* This is an exact affine dependence. The common-mean lemma uses `|X-Y| ≤ X+Y-2XY` on
  `[0,1]^2`: when `X ≥ Y` this is `2Y(1-X) ≥ 0`. Taking expectations gives `2μ(1-μ) ≤ 1/2`, so some cross pair has
  `d_ij ≤ W/2` and `Δ ≥ H^{-1/4}`. The pair is distinct because the supports are disjoint. Correct. An adversarial LP
  search reaches exactly 1/2, so the constant is attained.
* *Units.* The argument never uses the units except through the squareclass hypothesis. "Arbitrary Gaussian units"
  is therefore correct, **provided** the torsion bit (whether `ε_j/ε_5` lies in `±i`) is part of the dependence, as
  Astra says.

**Mod-2 criterion in one unit class.** When all `ε_j` agree, `∏ g_p^{m_p}` is a square iff all `m_p` are even, by
freeness modulo μ4; note `-1 = i^2`. So the hypothesis is exactly that the `4×n` matrix of the `b_p` has rank below
4 over F_2. Correct.

**Actual circles** (`windows.py`, 464 gcd-normalised five-windows, `N ≤ 4·10^11`). Inequality (1) held on every
pair. (4) held on all 232 squareclass-dependent windows, with min `Δ H^{1/4} = 19.2`. Also confirmed:
* even index always gives mod-2 deficiency, and odd index never does;
* no window was in one unit class. The unit bit matters in practice: 80 windows are exponent-deficient but not
  squareclass-deficient once the units are counted.

These windows are not near-extremal, so the run is a consistency check, not a stress test.

**Verdict on (4)/(5): CORRECT.** Novelty: the record has no mod-squares bound. Its "parity lemma"
(`round23/split.md` §1.4, `round24/theory.md` Lemma 2.1, Prop 5.1) concerns `h_{E_A} = h_{E_{A^c}}`, which is a
different statement. The trick `|t^2-1| = |t-1||t+1|` is classical. More important, (5) **is** the record's
monomial gap at `u = v/2`. Fakes satisfy it, and it sits inside the Liouville class of Lemma K1. Exponent 4 is
exactly the critical scale `Δ ≍ H^{-1/4}`, i.e. arcs of length about `√R` in the reduced radius. It is not better
than critical.

## 2. Balanced version (6) and the H_Σ bound

* *Cuts.* Every level `t ∈ [1, e_p]` gives a nontrivial cut because min = 0 and max = e_p. Hence `Σ_A log n_A = W`,
  and each pair is separated by exactly `2^3 = 8` of the 15 cuts (checked).
* *No affine dependence.* `8(1-η) > 15(1+η)/2` iff `η < 1/31`. So every `d_ij > W/2`, and the common-mean lemma
  rules out all affine dependences, including rank below 4. Correct.
* *Layer decomposition.* `|Σ c_i a_ip| ≤ Σ_t |Σ_{i∈A_t} c_i|`, with either side allowed since `Σc = 0`. Checked on
  5000 random profiles.
* *Costs.* Every four-mark word costs 12, and I checked all 120 placements, including words through mark 5. The
  two-mark word costs 8.
* *Constants.* `log D(s) ≤ 12(1+η)w`, so `Δ ≥ |s-1|/2 ≥ 2^{-1/2}e^{-3(1+η)w}`. With `W ≥ 15(1-η)w` the exponent is
  `5(1-η)/(1+η)`. Correct.
* *H_Σ.* In one unit class, `|y_a - y_b| = 2|Im A| e^{-d_ab/2}`, so `|t_ab| = |Im A| ≥ 1`, consistent with the
  expert note's definition. The diameter pair then gives
  `log|t| ≥ -3(1+η)w - (1/2)log 2 + 4(1-η)w - log 2 = (1-7η)w - O(1)`. The other nine `|t|` are at least 1. So
  `H_Σ ≥ (1-7η)w - O(1)`. Correct.

**Verdict: CORRECT.** It is new in the record, and it is a Liouville-type bound on a congruence-thin class.

## 3. "Reduction to three indices 1, 3, 5"

**Arithmetic.**
* With `g = 94/25`, the bracket `4(1-η) - 15(1+η)·25/94 = (1-751η)/94`, which is at least `1/188` for
  `η ≤ 1/1502`.
* `B_7 = (4/5)·490^{1/4} > (4/5)·4.7 = 3.76`, since `4.7^4 = 487.97 < 490`. Also `94/25 > 124/33 > 15/4`.
* For an even index, `Z^4/L` has even order, so some surjection `Z^4 → Z/2` kills `L`, and the columns do not span
  F_2^4. For an odd index, `L + 2Z^4 = Z^4`, so parity never helps at index 1, 3 or 5.
* Balance excludes rank below 4 through the common-mean lemma.

The bookkeeping is right.

**What it gives.** Assume the section 2 theorem for `d_0 = 7`, which is another block. Then `(L_5°)`, in the
one-unit-class formulation with `η_0 = 1/1502` and `c = min(1/188, 1 - 7η_0) = 1/188`, holds outside `Z(7, 94/25)`
for every balanced primitive tuple with `I ∉ {1,3,5}`. The parity part (all even `I`) needs no exceptional set and
no section 2.

**What it does not give.**
1. Nothing at `I = 1, 3, 5`.
2. Indices 3 and 5 cannot be handled by the same cover argument. `B_3 ≈ 3.03` and `B_5 ≈ 3.46` are below 15/4, the
   threshold needed for balanced `(L_5°)`. 7 is the first degree at which the amplified McKinnon–Roth constant
   `0.8(70d)^{1/4}` passes 15/4.
3. No `A_5(δ)` statement for arbitrary profiles at small index.
4. No unit-uniform form. Astra correctly flags the torsion coordinate.

**Does index 1 contain essentially all growth-regime configurations? Yes.**
* Index greater than 1 needs a character `Z^4 → Z/q` that kills every column. For a balanced tuple, at least 15
  threshold layers spread over its primes, so that is a congruence coincidence of codimension growing with the
  number of primes.
* In my synthetic run (`identity8_family.py` (c), random profiles), the fraction with `I > 1` falls from 59% at 5
  primes to 3.5% at 8 primes, about 2% at 10 and 0% at 14–15 primes.
* Every explicit index-one family has `I = 1`: the 15-cut squarefree profiles, the ten-pair-cut profile, and
  balanced profiles from a single simple singleton cut at each mark.
* Exponent-two structure (for example squares of clusters) is thin, and it is not short-arc anyway.

So the reduction removes thin sets and leaves the generic case whole. "It suffices to treat three indices" is
literally true, but the remaining case is the whole problem.

**Verdict: CORRECT, conditional on §2 for `I ≥ 7`.** The wording "reduction" overstates the progress, and the
report should say that index 1 is generic.

## 4. The Zariski-dense, balanced, full-rank parity family

The construction checks out:
* There is one prime per nonempty `A ⊆ [4]`, with exponent `e_A` on `A` and 0 elsewhere. That makes the tuple
  primitive.
* The singleton columns give rank 4.
* The odd columns `1_{12}, 1_{13}, 1_{14}` span the parity hyperplane, so the mod-2 rank is 3. The unique relation
  has support `{1,2,3,4}`.
* In one unit class with unit 1, (6) applies.

I built two instances exactly, at `w = 60` and `w = 120`, with singleton exponent 6, pair exponents 1 or 2, primes up
to about `e^{120}`, and Cornacchia for the large primes. They are primitive, with rank 4, mod-2 rank 3 and
`η ≤ 4·10^{-4}`, and the layer bound behind (6) holds for all three words.

**Correction to the density step.** As written, `e_A log p_A ∈ [w, w+√w]` does not by itself let the four
singleton exponents vary over sets whose size tends to infinity. If `p_A ≈ e^{w/e_A}`, the interval holds only
`O(√w/log p_A)` admissible exponents. The argument needs the singleton primes small, say fixed or with
`log p_{A} = o(√w)`, and their even exponents of order `w/log p_A`. Then each factor has about `√w/(2 log p_A) → ∞`
elements. Also, `g_p^e` are distinct because `g_p` has infinite order, and a nonzero Laurent polynomial cannot
vanish on a grid whose factors exceed its degrees. With that choice, made explicit, the statement holds, with
`η ≍ w^{-1/2}`.

**Verdict: CORRECT WITH CORRECTIONS.** The singleton primes must be fixed, or satisfy `log p = o(√w)`, with growing
even exponents. Significance: it only shows that the parity class is not contained in a proper closed set, so (6)
is not vacuous. It says nothing about archimedean proximity, as Astra notes.

## 5. Squarefree bounds (7): at most nine and at most six primes

* *Certificate for (7).* For `E = {1,2}`, the seven pairs meeting `E` give degrees 4, 4, 2, 2, 2. The cut `E`
  crosses 6 of them and every other cut at most 4. I checked this by the degree formula and exhaustively for all 10
  two-mark cuts. Hence `min d ≤ (4W + 2w_E)/7`, and `Δ ≥ H^{-2/7}` (exponent 7/2) when some two-mark cut is absent.
  This applies to all profiles, not only squarefree ones.
* *At most nine primes.* Each prime of a squarefree norm gives one cut, so with at most 9 primes some two-mark cut is
  absent. Correct.
* *At most six primes.* There are two cases:
  * Two disjoint absent edges: their four connecting pairs form a 4-cycle. Only the two absent cuts cross all four,
    and other cuts cross 0 or 2, since the number of cycle edges a cut crosses is even.
  * Absent edges pairwise intersecting: with at least four of them, they form a star at some vertex `v`, since a
    triangle has only 3 edges. Any further absent edge would be disjoint from a star edge, which is the first case.
    So the six present cuts are exactly the six edges away from `v`, and each crosses 2 of the 4 star edges.

  In both cases `min d ≤ W/2`, so the exponent is 4. Correct.

**Exact LP check** (`combin.py`). This is the minimum over distributions on pairs of the maximum crossing fraction:
* all 15 cuts: 3/5 (Cilleruelo–Córdoba);
* the ten two-mark cuts: 3/5;
* one two-mark cut missing: exactly 4/7;
* every 6-cut support: at most 1/2;
* some 7-cut support: 4/7;
* every 9-cut support: at most 4/7.

So both constants are optimal for pair-Liouville averaging. "Six" cannot be raised to seven at exponent 4 by this
method.

**Verdict: CORRECT.** Novelty: a routine refinement of the record's `Σ d_ab = Σ|A||A^c| log n_A ≤ 6W'`
(`round20/gauss.md` §4.4). Significance: nil for growth. Near-extremal clusters have `Θ(M)` prime factors
(`overview/overview.md`, around line 536), and balanced five-tuples need all 15 cuts.

## 6. The ten-pair-cut angle model at 10/3

* **Lattice.** In `Z^5/Z1` the columns are `1_E`. Their differences give `e_i - e_j`. The quotient by the sum-zero
  image is `Z/5`, where `1_E` has class 2, a generator. So the columns span everything and `I = 1`. My Smith form of
  the anchored `4×10` matrix is `(1,1,1,1)`. The kernel torus is therefore connected, and the row lattice is
  saturated: an integral `v` in the real row space is in the integral row lattice. Correct.
* **Characters.** The model angle is `θ_x = Σ_E (1_{x∈E} - 1_{5∈E}) φ_E`, so the point angles are `2θ_x`, consistent
  with `round4/outside.md` (1.1). For `v = B_0^T c`, `<v,φ> ≡ Σ c_jθ_j (mod 2π)`, and `v_E = c_i + c_j` with
  `Σ_{1..5} c = 0`. Correct.
* **Classification.**
  * `2M(c) ≥ 3||c||_1`, from `Σ_{j≠i}(c_i+c_j) = 3c_i`.
  * Roots have `M = 6`. The three norm-4 patterns give 12, 12 and 8. Norm at least 6 gives `M ≥ 9`.
  * My census over `c ∈ [-4,4]^4` confirms that roots have `M = 6` and everything else `M ≥ 8`. The counts are 20
    roots and 30 vectors at `M = 8`, then 12, 14, 16.
* **Gaps — both families are required, and both are handled.** `round4/outside.md` Lemma 0.1, which `round8/lw.md`
  Thm 6.1 item 2 uses, has two families: the axial gap `dist(<v,φ>, (π/2)Z) ≥ arcsin e^{-w/2}` and the diagonal gap
  `dist(<v,φ>, π/4 + (π/2)Z) ≥ arcsin(e^{-w/2}/√2)`. The expert note §1 displays only the axial one. Astra imposes
  both. Its probability bounds cover the union, inside radius `arcsin e^{-w/2}` about `(π/4)Z`.
  * Roots: deterministic. The angles lie in `[9s, 41s]` (`[9.5s, 40.5s]` for roots through mark 5). Since
    `e^{-w/2} ≤ X^{-3} = s`, both gaps hold, with a factor of about 9 to spare.
  * Non-roots: one free coordinate on an interval of length `s`. The probability is `O((1+s^{-1})X^{-M/2}) = O(X^{3-M/2})`,
    and `M ≥ 8` gives `X^{-1}`. The count is polynomial in `M`. Correct.
  * Off the row lattice: Haar-uniform, total `O(Σ_{v≠0} X^{-||v||_1/2}) = O(X^{-1/2})`. Correct.
* **Avoiding Z.** `θ` has a density on an open 4-box, and `U(R)^4` is Zariski dense in `U^4`. So a proper closed `Z`
  meets the box in measure zero. Correct.
* **Numerics** (`angle_model.py`). With `X = 10^6`, ten primes `≡ 1 (4)` in `[X, 1.001X]`, 80-digit arithmetic and
  Haar sampling on the fibre, I ran 6 samples. Each checked 3960 vectors: all `v` with `||v||_1 ≤ 3` and all row
  vectors with `c ∈ [-3,3]^4`. There were 0 failures in either family. `Δ/s ≈ 80`, and `log Δ/log H = -0.268`,
  which tends to `-0.3` as `X → ∞` (the offset is `log 80`).

**Scope of "sharp".** The model attains the Liouville and Cilleruelo–Córdoba value `Δ ≍ H^{-3/10}`, i.e. reduced
diameter about `R^{2/5}`. The LP above shows that pair averaging on this profile gives exactly 3/5. So the
statement is "gap-only inputs at index one give exactly 10/3". It is **not** a statement about actual circles.
Whether five lattice points lie on arcs `O(R^{2/5})` infinitely often is open (`overview/overview.md` Cor 3.8).
Actual circles also satisfy residue constraints the model ignores: Theorem 3.6 gives `6 | L_0` and
`σ∏s_ij ≤ K^{10}/128`. Those affect constants only.

**Verdict: CORRECT.**
* *Novelty.* The ingredients are all in the record:
  * the fake construction (`round8/lw.md` Thm 6.1 random-choice proof; `round20/gauss.md` §5 item 1 and
    `kill_gauss.md` §5 item 1, balanced 5-row fakes at `R''^{7/15}`);
  * the 10/3 value (Cilleruelo–Córdoba, `round20/gauss.md` §4.4; `round22/lit.md`; `τ_e = 10/3`, `round23/width.md`
    Thm 2.5);
  * the ten 2|3 cuts profile as the expected uniform witness (`round23/width.md` §5 item 2, for sections,
    conditional on (W)).

  An explicit unbalanced index-one angle model at 10/3 with both gap families is not written down in the record. It
  is a routine consequence.
* *Significance.* Peripheral. It confirms that the index-1 case cannot be settled by gaps plus the index. That is
  Lemma K1 again.

## 7. Fifteen balanced simple cuts at index one

The four singleton columns are `e_1, …, e_4`, so `I = 1`. The weights lie in `[log X, log X + log 2]`, and the tuple
is primitive. Trivially **CORRECT**. These profiles also have mod-2 rank 4, so (6) is unavailable, as section 6
says.

## 8. Identity (8)

Write `z_i = G·P_i/P̄_i` with `G = ∏ β̄_E`. Then:
* `P̄_iP_j = n_ij Q̄_iQ_j`, so `t_ij = det(Q_i, Q_j) ∈ Z`;
* `D(y_i/y_j) = |P_i|^2|P_j|^2/n_ij^2`, since `Q_iQ̄_j` is conjugate-primitive by disjoint prime supports;
* `|y_i - y_j| = 2|det|/(|P_i||P_j|)`, which matches the record's `t_ab` (`|ζ_a-ζ_b| = 2|t|R e^{-d/2}`).

**The identity.** In a frame with all `Y_i ≠ 0`, write `det(P_i,P_j) = Y_iY_j(s_i - s_j)`. Then
`F(P_i)/∏_{j≠i}det(P_i,P_j) = F(s_i,1)/(∏_k Y_k ∏_{j≠i}(s_i-s_j))`. Summing over `i` gives the degree-4 Lagrange
coefficient of a cubic, which is 0. Take `F = |P|^2 P` and use `|P_i|^2 = ∏_{j≠i} n_ij`. Clearing denominators gives
integer coefficients of absolute value at most `∏|t| = e^{H_Σ}`. Correct.

**Check.** It verified exactly (`Fraction` arithmetic) on 200 actual primitive tuples, with blocks of one or two
Gaussian primes to powers 1 or 2 and random orientations. Integrality, the norm identity and
`|y_i - y_j| = 2|t_ij|e^{-d_ij/2}` (to 300 digits) were verified too.

**Verdict: CORRECT.**
* *Novelty.* Mildly new. The record's Ptolemy/Plücker system explicitly never sees the Gaussian norms `|P_i|^2`
  (`growth/ptolemy.md` §0 and Thm B(iii)), and (8) does use them.
* *Significance.* Nil as it stands. It is a polynomial identity with forced divisibility, the class that
  `round20/gauss.md` §4.4 places under the auxiliary-polynomial cap, and Astra itself says it gives no bound.

## 9. Section 6, "exact loss"

At `L = Z^4`, the lift has degree 1 and `L* = Z^4`. Both statements are trivially correct. "All monomial gaps and the
index permit 10/3" is the content of section 6 of this report. The section 5 remark is outside this block.
**CORRECT.**

## 10. Standing cautions

* *Fixed-S exceptional sets.* None are used in this block; sections 3–4 have no exceptional set at all. The section 2
  dependence in the reduction is flagged as conditional.
* *Gaussian gcd and denominator costs.* Kept: primitivity, one unit class, and the torsion bit.
* *Drifting phase.* The angle model keeps the common phase free (`θ_0`), and Astra's `Δ` is phase-free.
* *Template-dependent thresholds.* The thresholds `η < 1/31` and `η ≤ 1/1502` and the exponents are absolute.
