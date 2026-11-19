# Referee report, block "index" (round 31): astra_free_II/answer.md §II.4

**Scope.** I checked every claim in §II.4, lines 545–605 of the frozen `astra_free_II/answer.md`:
* the normalized counting theorem, with its pieces (real affine independence, the bounded minor, the two free rows, parity);
* the effect on the grid test;
* the Paley double core and its HNF certificate;
* the Baker–Wüstholz removal of the growth hypothesis;
* "(L_5°) on index one suffices";
* the calibration paragraph;
* the selection theorem (m/145 − K_Q);
* the finite Gini certificate.

I also checked the Part I §5 extraction this block relies on, the logical chain to uniformity, and how the
statement compares with the record's (L_5°).

**Scripts** are in `round31/ref2_index/`. They are my own code; I used no sympy HNF, except as an oracle in a
self-test.
* `lat.py` computes lattice indices by integer Euclidean elimination. `test_lat.py` validates it against the gcd of
  minors on 400 random cases, with 0 mismatches.
* `paley.py` checks the skew double core for q = 3, 7, 11, 19, 23, 31.
* `gini_cert.py` reproduces the Gini certificate. It uses the CDF formula, enumerates all of F_q^5, and checks
  that the mixture GMD is quadratic.
* `ineq.py` gives numerical checks of every inequality used.
* `clusters.py` and `clusters.out` cover the record's genuine clusters.

All runs used `prlimit --as=3.5e9`.

## Verdicts

| # | claim | verdict |
|---|---|---|
| 1 | Real affine independence of all rows from pair separation | CORRECT (known in substance) |
| 2 | e_p log p ≤ 4W/m; bounded minor (2D)^4 at W ≤ b m log m | CORRECT (first half is in the record) |
| 3 | Two free rows: O_q(m^3) rank losses mod a fixed odd q, with no hypothesis on W | CORRECT, new |
| 4 | Parity: no five-tuple loses rank mod 2 | CORRECT (mod-2 Sidon content already in the parallel thread) |
| 5 | Normalized counting theorem: all but O_b(m^3) tuples have index one | CORRECT |
| 6 | Effect on the grid test | CORRECT (immediate) |
| 7 | Paley double core: proper restrictions have index 1, the whole profile has index M/2 | CORRECT; reproduced for M = 4, 8, 12, 20, 24, 32 |
| 8 | Baker–Wüstholz gives an absolute index bound on genuine shortened clusters, hence O(m^3) at any growth rate | CORRECT (rests on lemma §2.2, refereed CORRECT in ref_baker.md) |
| 9 | Index-one (L_5°) suffices for the uniform bound; index-one (5) suffices for growth | CORRECT |
| 10 | Calibration: H ≤ N^{31/33}; 15/16 < 31/33 | CORRECT |
| 11 | Selection theorem: ≥ m/145 − K_Q points with every five-tuple of index one | CORRECT |
| 12 | Finite Gini certificate 315/256, 239/128, 633/256 | CORRECT (reproduced). Not used by the final proof |

No GAP and no FALSE item in this block.

## 1. Real affine independence (line 554)

**Pair gap.** The pair gap gives w(a_i − a_j) > W/2:
* The phase difference is within Δ of 2πZ, so ⟨a_i−a_j, φ⟩ lies within Δ/2 of πZ.
* The axial gap gives e^{−w/2} ≤ sin(Δ/2) ≤ ¼e^{−W/4}.

**Independence.** Suppose there were an affine dependence. Split it into positive and negative parts and normalize
both to probability laws P and N with equal barycentre and disjoint supports. At each prime, X, Y ∈ [0,1] have the
same mean μ. Then |x−y| ≤ x+y−2xy, because the difference of the two sides is 2·min(x,y)(1−max(x,y)) ≥ 0. So
E|X−Y| ≤ 2μ(1−μ) ≤ ½. Weighting by e_p log p gives a cross pair, of distinct rows, at distance ≤ W/2. That is a
contradiction.

**Verdict.** CORRECT.

**Novelty.** Low. The record already has affine independence of actual clusters (`round5/referee_construct.md:233`,
via `linear_allocation_affine_rigidity.md` for C ≤ √2; `round10/determinant.md` Lemma 3.2). The probabilistic proof
from pair separation alone is neat and also applies to gap-admissible angle models.

## 2. Bounded minor at the growth scale (lines 556–558)

**The e_p bound.**
* Choose s ≥ m/2 rows in one closed half of [0, e_p]. Their weighted range is W′ ≤ W − (e_p/2) log p.
* Each threshold separates ≤ s²/4 of their pairs, so the mean pair distance is ≤ sW′/(2(s−1)).
* Since this mean is > W/2, we get e_p log p < 2W/s ≤ 4W/m.

This is the record's `round8/stability.md` Thm 3.3(2), which has the same 4W/M.

**Small and large primes.**
* Primes p < √m contribute ≤ π(√m)·4W/m ≤ 4W/√m.
* For p ≥ √m, e_p ≤ 4W/(m · ½ log m) ≤ 8b.

**Bounded relation.** Suppose the rank on the retained primes were r < 4.
* Cramer's rule gives a zero-sum integer relation c with ‖c‖₁ ≤ O(r!(2D)^r).
* The full residual v ≠ 0 by §1. It is supported on small primes, so w(v) ≤ ‖c‖₁·4W/√m.
* Its phase is within ‖c‖₁Δ/2 of πZ, because the coefficients sum to zero and the common phase cancels exactly.
  It is not set to zero.
* The gap is violated once m is large. W → ∞ because W ≥ (m−1) log 5.

**Minor and index.** Hadamard on 4×4 minors with entries in [−D, D] gives ≤ (2D)^4. The index divides every
nonzero minor.

**Verdict.** CORRECT.

## 3. Two free rows, fixed odd q (line 560): CORRECT

I went through it line by line.

**Rank loss.** Rank loss mod q is equivalent to a nonzero zero-sum c ∈ F_q^5 with Σc_i a_{T_i} ≡ 0 coordinatewise.

**The branch argument.**
* Fix c and the rows in three slots. Use centered representatives h, k, with L = |h|+|k| ≤ q−1.
* The integer sums s = h a_x + k a_y lie in one class r mod q.
* exp(2i⟨(s−r)/q, φ⟩)^q = e^{i(hθ_x+kθ_y)}·const lies in an arc of width LΔ. The centre is (h+k)θ_0 plus fixed
  data: it is common to all (x, y) for this c and these fixed slots, so the drift is handled correctly.
* Inside one of the q inverse arcs, v = (s−s′)/q is integral and nonzero. Its angle is within LΔ/(2q) < Δ/2 of πZ.
  So w(s−s′) > qW/2.
* The sums lie in a box of weighted width LW, and the threshold count gives mean distance ≤ nLW/(2(n−1)).
  Therefore n < q/(q−L) ≤ q.

**Representations.** Each exact sum has ≤ 2 ordered representations. The signed measures hδ_x + kδ_y must agree,
by affine independence of all rows.

**Total.** ≤ 2q² completions per (c, fixed slots), so O_q(m^3) tuples. No hypothesis on W is used. The argument is
valid for angle models with the axial gap.

## 4. Parity (line 566): CORRECT

**Two marks.** For two marks, v = (a_x−a_y)/2 has w(v) ≤ W/2.

**Four marks.** Average over the three pairings, layer by layer:
* minority of size 0 gives 0;
* minority of size 1 gives 1 for every pairing;
* minority of size 2 gives 2, 0, 0, average 2/3.

So some pairing gives w(v) ≤ W/2.

**Phase.** ⟨v, φ⟩ lies within Δ/2 of (π/2)Z (for two marks, within Δ/4). This contradicts the gap with room to
spare.

**Novelty.** The four-mark case is the "clusters are Sidon mod 2" observation the coordinator attributes to the
other thread (reply2.txt item 3). Stating it as "no rank loss mod 2" is new but immediate.

## 5. Counting theorem (line 552): CORRECT

It combines §2 (primes dividing the index are ≤ (2D)^4), §3 (O_q(m^3) for each odd q) and §4 (q = 2 is empty).
Rank over Q is always 4, by §1. The index is intrinsic: it is invariant under conjugating π_p, changing the anchor,
or dividing by the tuple gcd, since constant columns contribute zero.

The earlier draft (`ii_index_count.md`) claimed O(m^4). The final O(m^3) is right.

## 6. Effect on the grid test (line 568)

O(m^3) < d·m^4 eventually. CORRECT and immediate.

## 7. Paley double core (lines 570–572): CORRECT, reproduced

**The column difference.** For each label a, the two binary columns differ by e_a − e_0. They have opposite signs
at the border and at a.

**Proper restrictions.** For a proper restriction S:
* If the border is not in S, restriction gives every e_a.
* If the border is in S, an omitted a gives −e_0, and then every e_b follows.

**Whole profile.** Every column B has ΣB = M/2. The column sums of K vanish: −1 + 1 + Σχ = 0. So
L + Z1 = {x : Σx ∈ (M/2)Z}, and the index is M/2.

**My computation (`paley.py`).**

| q | M | whole index | five-subsets with index ≠ 1 | proper subsets (size ≥ 2) with index ≠ 1 |
|---|---|---|---|---|
| 3 | 4 | 2 | — | 0/10 |
| 7 | 8 | 4 | 0/56 | 0/246 |
| 11 | 12 | 6 | 0/792 | 0/4082 |
| 19 | 20 | 10 | 0/15504 | not run |
| 23 | 24 | 12 | 0/42504 | not run |
| 31 | 32 | 16 | 0/201376 | not run |

* In every case the minimum pair Hamming distance is M > W/2 = q.
* K Kᵀ = M·I was checked.

**Assessment.** A nice illustration, but peripheral. The record already knew that the double core's characters
have (M/2)λ integral (`round8/lw.md` §5.1). The fact that this is invisible on every proper sub-tuple is new and
explains why five-point index arguments cannot see it.

## 8. Baker–Wüstholz removal (lines 574–581): CORRECT

**The deduction.**
* H | N.
* δ_T ≤ ½N^{−1/4} ≤ ½H^{−1/4}, so −log δ_T / W_T ≥ ¼.
* Then (2) of the lemma thread, with an absolute constant, bounds I absolutely for large W_T. For bounded W_T,
  I ≤ 24(W_T/log 5)^4.
* Hence every prime dividing an index is ≤ I_abs, and the counting theorem gives O(m^3) with an absolute constant
  at any growth rate.

**The completeness sketch.** It repeats lemma §2.2: Minkowski, n ≤ I, 2h(u) ≤ 4‖v‖_∞W_T, |u^{4n}−1| ≤ 4‖c‖₁δ_T,
‖c‖₁ ≤ 8I. Every step was verified in `ref_baker.md` §§1–6, and I re-read it.

**Source check.** The cited statement is the standard two-logarithm specialization of Baker–Wüstholz (Crelle 442,
1993, 19–62): log|Λ| > −c(n, d)·h′(α_1)···h′(α_n)·log B, with linear dependence on log B. I confirmed the paper and
this shape only through secondary sources (arXiv:1611.04664; the Cambridge monograph listing). I did not open the
paper or check the explicit constant. The deduction only needs an absolute constant, and BW is effective.

**Standing cautions.** No fixed-S exceptional set is used, the gcd is accounted for, and the phase drift is handled.

## 9. "(L_5°) on index one suffices" (line 583): CORRECT

**The form used.** It is exactly the record's (L_5°) (expert_note §3; gauss.md §3, Thm 3.4′), with its hypothesis
restricted to tuples of index one:
* fixed proper closed Z ⊂ T_5;
* absolute η, γ, w_0;
* one unit class, primitive;
* all fifteen cuts in [(1−η)w, (1+η)w];
* conclusion H_Σ = Σ log|t_ab| ≥ γw.

Index-one (L_5°) is formally weaker than (L_5°), so the implication is stronger than Thm 3.4′.

**Chain.**
1. Subdivide, keep one unit class, divide by the gcd.
2. Part I §5, which I re-checked and which ref_sieve.md also found CORRECT, gives ≫ m^5 ordered distinct tuples
   with all 31 Walsh coefficients O(m^{−1/2}).
   * E f_J² ≤ 2/m via ‖G‖ ≤ 2 for |J| = 1, and tr(G²)/m² for |J| ≥ 2.
   * Hence relative cut error O(m^{−1/2}) with w = W/16.
   * H_Σ = O(W/m), from 0 ≤ log t_ij ≤ L_ij/2 and Σ_{i<j} L_ij ≤ mW/4.
   * Here t_ij and d_ij are unchanged by tuple-gcd reduction.
3. Remove O(m^3) tuples of higher index (§8). The constant is absolute, which uniformity needs. For angle models
   only the growth-scale version is available, but the uniform statement concerns genuine points.
4. Schwartz–Zippel: zeros of the homogenized Laurent polynomial vanishing on Z number ≤ deg·m^4.
5. A surviving tuple has H_Σ = O(W/m) < γW/16 once m is large. This is a contradiction.

**The growth variant.** The growth variant (5) restricted to index one is likewise correct, using the O_b(m^3)
count at W ≤ b m log m.

## 10. Calibration (line 585): CORRECT

* A_5(1/31) requires δ_T ≤ H^{−(1+2/31)/4} = H^{−33/124}.
* The cluster gives ½N^{−1/4}. This suffices iff H ≤ N^{31/33}.
* The typical W_T is (15/16 + O(m^{−1/2}))W, since the two constant patterns carry 2/32.
* 15·33 = 495 < 496 = 31·16, so the slack is 1/528.

## 11. Selection theorem (lines 589–603): CORRECT

**Eligibility.** At η = 1/200, 8(1+η)/(15(1−η)) = 0.53869 < 13/24 = 0.54167. So G_ij ≥ −0.0774 ≥ −1/12, and
balanced tuples are eligible. W ≥ W_T is used in the right direction.

**Scaling.** The mean of the centered square of a uniform nonzero residue is q(q+1)/12. I checked this exactly for
q = 3, 5, 7, 11, 13. Averaging over scalars gives S_2 ≤ 5q(q+1)/12.

**Lower bound.** For fixed c, the arc centre (Σc_i)θ_0 is common. The branches have width LΔ/q with L/q ≤ 5/2, so
two distinct sums in one branch have w > qW/2. The mean is ≥ (1 − 1/k)(qW/2 − O(q)).

**Upper bound.**
* The layer decomposition and the triangle inequality give the first step.
* GMD(ν_ℓ) ≤ GMD(½ν_ℓ + ½ν̃_ℓ) ≤ GMD(pooled), by concavity of 2∫F(1−F) and invariance of GMD under reflection.
* GMD ≤ (2/√3)·sd follows from E|Z−Z′| = 2∫(2u−1)(H(u)−EZ) du and Cauchy–Schwarz. Numerically, the worst ratio
  over 20,000 random laws was 1.1402 against a bound of 1.1547.
* The flip makes every mean equal to Σc/2. Var = ¼cᵀG_Tc; I checked this identity on a random profile to 1e−15.
* cᵀGc ≤ S_2 + (1/12)(L² − S_2) ≤ (4/3)S_2, since L² ≤ 5S_2.
* So GMD² ≤ 4S_2/9 ≤ 5q(q+1)/27 ≤ 20q²/81 < q²/4. The middle inequality needs q ≥ 3.

**Conclusion of the counting step.** k ≤ 1/(1 − 2√(20/81)) ≈ 161.5 for large W, with a threshold depending only on
q ≤ Q. The active row set is determined by the exact sum, through the signed measures and full affine independence.
Deleting O(Q^5·161·5) marks removes every eligible tuple with rank loss at any q ≤ Q.

**Independent set.**
* ‖I − G‖ ≤ 1, since I − G is nonnegative with λ_max ≤ 1, and Perron–Frobenius.
* So Σ_{j≠i} G_ij² ≤ 1, and the degree in {|G_ij| > 1/12} is < 144.
* Greedy selection gives ≥ m/145 points.
* Random test: the maximum row sum was 0.99996.

**Inputs.** For genuine clusters Q is absolute, by §8. For angle models Q = Q_b in the growth regime.

## 12. Finite Gini certificate (line 605): CORRECT

My own code (CDF formula; full enumeration of F_q^5; quadraticity in r checked at 3 extra points) gives:

| q | classes | value | worst c |
|---|---|---|---|
| 3 | 4 | 315/256 | (1,1,1,1,2) |
| 5 | 8 | 239/128 | (1,1,2,2,4) |
| 7 | 13 | 633/256 | (1,1,3,4,5) |

All are below q/2.

This concerns the balanced-mixture variant from the draft. The final proof uses the eligibility/S_2 bound instead,
so the certificate is a side check only.

## Numerical test on the record's clusters (`clusters.out`)

**Result.** No explicit genuine cluster in `growth/families.md` satisfies the theorem's hypotheses, so the theorem
cannot be tested on genuine data:

| cluster | width over the required ½e^{−W/4} | pairs at sep ≤ W/2 | affine rank | 5-tuple index distribution |
|---|---|---|---|---|
| 8-point record (R² = 32988780417125) | 11× too long | 12/28 | 4 | {0: 24, 1: 32} |
| 7-point record (235370525) | 10.5× too long | 8/21 | — | {0: 3, 1: 7, 3: 5, 4: 6} |
| golden n = 4 | 11.8× too long | 12/28 | 4 | {0: 24, 1: 32} |

In the record's data, the shortened genuine clusters have ≤ 4 points (C ≈ 0.50 at 3 points, 0.88 at 4 points).
The cyclotomic five-row family (C → 0) is a single tuple. At block level its index is 1; I did not factor it at
10^3428.

**What this does show.** Without shortening, the theorem's conclusions genuinely fail:
* rank deficiency (index 0) occurs;
* index 4, i.e. rank loss mod 2, occurs;
* index 3 occurs.

So the hypotheses are doing real work. The only available test of the theorem is the fake double core (§7), which
is consistent with it.

## Novelty and significance

**New to the record.** grep finds no prior five-tuple index discussion. New items:
* the O(m^3) count;
* the two-slot argument;
* the selection theorem;
* the double-core observation that every proper restriction has index 1.

**Already known.** Affine independence, the e_p log p ≤ 4W/m bound and the mod-2 Sidon fact are known or
essentially known.

**Significance: a clean negative clarification for Question B; peripheral to the rates.**
* Shortened clusters live at index one up to O(m^3) tuples, and there is a positive-proportion subcluster with
  every five-tuple of index one. So the parallel thread's high-index cover theorem (A_5(1/31) for I ≥ 7) and its
  "indices 1, 3, 5" reduction contribute nothing to the uniform bound or to growth.
* The whole open content of (L_5°) for clusters is the index-one case. That case already contains every hard
  model: fakes, the ten-pair-cut angle model at 10/3, and balanced squarefree profiles.
* "Index-one (L_5°) suffices" is formally weaker than (L_5°) but not easier in substance.
* Nothing here changes the growth rate (2/3+ε) log R/log log R, the uniform bound, F1/F2, or route 5.
