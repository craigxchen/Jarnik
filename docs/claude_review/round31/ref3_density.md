# Referee report, round 31, Part III, block "density" (astra_free_III §III.9)

Scope: §III.9 of `astra_free_III/answer.md`, plus the parts of §III.1 (the exact test (III.1), the matching
factors A,B,C, the "each triple forces at most one subset" uniqueness) and §III.7 (the bound (III.5), the
exposure argument) that its proof uses. The parallel derivation in `astra_free_III/iii_graph_notes.md` lines 141-153
was read too. Scripts are in `round31/ref3_density/` (`thresholds.py`, `experiment.py`).

## Summary verdict

The proof of III.9 is **correct**. I checked it line by line. There is one presentational repair: the exposure
step needs the sample in uniformly random order, not in increasing order. There are two scope clarifications: which
circles are covered, and the fact that the exceptional set depends on C. All the arithmetic (F_M, (III.7) at
M=17/18 for D=4, 13/14 for large D, 37/38 for D=3/2) was reproduced with exact rationals.

**But the headline numbers are not the best ones available, and for D>4 the statement is weaker than an
elementary fact.** A one-line first-moment count over **pairs** proves the following. For D>4 almost every parent
in the family has **no two lattice points** on any arc of length C√R. For D=4 (indeed for 8/3<D≤4) it has
at most **2**. For every D>2 the bound is a constant. Astra's quartet argument beats this pair count only for
1<D≲2.08, and it is genuinely needed only for 1<D≤2 (e.g. 37 marks at D=3/2). The bounds 17 (D=4) and 13
(large D) are therefore superseded by 2 and 1. The result concerns a density-zero, specially shaped family of
parents. It says nothing about the worst case and has no effect on the growth rate, the uniform bound,
Questions A or B, or route 5. **Peripheral.**

## 1. Definitions used by the proof (III.1, III.5): CORRECT

* One-unit convention: lattice points of norm N = ∏p^{e_p} (split p) are z_i = u·Q̄·P_i/P̄_i with
  P_i = ∏π_p^{a_ip} and Q = ∏π_p^{e_p}. Hence arg z_i = const + 2 arg P_i (doubled arguments). g_ij = ∏p^{min(a_ip,a_jp)}
  divides P̄_iP_j, so t_ij = [P_i,P_j]/g_ij ∈ ℤ. It is nonzero because z_i≠z_j forces P̄_iP_j ∉ ℝ. Also
  |t_ij| ≤ √(n_in_j)|sin(θ/2)|/g_ij ≤ √d_ij·Δ/2, with d_ij = n_in_j/g_ij² = ∏p^{|a_ip−a_jp|}. With Δ ≤ ½N^{−1/4} this
  is (III.5).
* Plücker: g12g34·x − g13g24·y + g14g23·z = 0 with x = t12t34, y = t13t24, z = t14t23. After division by the common gcd the
  surviving exponent at p is h3−h2 on the matching that pairs the two lowest (and the two highest) values. I re-derived
  this: that matching gets h1+h3 and the other two get h1+h2. So A, B, C are pairwise coprime and ABC | N. For
  squarefree N, B is exactly the product of the primes whose cut splits the quartet 13|24, and similarly for A and C.
* "Each triple forces at most one subset": after exposing A and C, B = (Ax+Cz)/y is a single integer. By unique
  factorization its prime support is a single set of primes. This holds also with exponents ≤ E. CORRECT.
* Pair separation: 1 ≤ |t_ij| gives log d_ij ≥ W/2 + log 16 (W = log N). For squarefree N, d_ij is the product
  of the primes separating i and j. CORRECT.

## 2. Pair-separation mean, the middle-class probability, convexity: CORRECT

* Σ_{i<j} log d_ij = Σ_p (log p)·u_p with u_p = k_p(M−k_p). Pair separation gives a weighted mean of u that is
  ≥ M(M−1)/4. All log p lie in [D log r, D log r + log 2] and u ≤ M²/4, so the unweighted mean is ≥ M(M−1)/4 − O(M²/log r). CORRECT.
* For a uniform ordered quartet of distinct marks, a cut with k marks on one side puts {1,3} on one side and {2,4}
  on the other in 2·k(k−1)(M−k)(M−k−1) of the (M)_4 cases. The identity (k−1)(M−k−1) = u−M+1 gives
  2u(u−M+1)/(M)_4. CORRECT.
* f(u) = u(u−M+1) is convex and increasing for u ≥ (M−1)/2 ≤ M(M−1)/4. Jensen therefore gives E f ≥ f(M(M−1)/4) − o(1)
  = M(M−1)²(M−4)/16 − o(1). Multiplying by 2/(M)_4 gives F_M = (M−1)(M−4)/(8(M−2)(M−3)). CORRECT (M ≥ 5 for positivity).

## 3. E log V ≤ 3W/(2(M−1)): CORRECT

From (III.5), log|t_ij| ≤ ½log d_ij − W/4 − log 4. The mean of log d_ij over pairs is ≤ W·⌊M²/4⌋/C(M,2) ≤ WM/(2(M−1)).
So the mean of log|t| is ≤ W/(4(M−1)), and each of the six pairs of a uniform ordered quartet is a uniform pair. V ≥ 1
because |t| ≥ 1.

## 4. Inequality (III.7) and its arithmetic: CORRECT (reproduced exactly, `thresholds.py`)

E[sH − log V] ≥ r log r·[(D−1)F_M − 3D/(2(M−1)) − o(1)], with H = log((K−r)/r) = (D−1+o(1))log r and W = (D+o(1))r log r.
K ~ r^D/(2D log r) by the PNT in progressions, so K ≫ r needs D>1.

| case | (D−1)F_M | 3D/(2(M−1)) | holds |
|---|---|---|---|
| D=4, M=18 | 0.371875 | 0.352941 | yes |
| D=4, M=17 | 0.371429 | 0.375 | **no** |
| D→∞, M=14 / 13 | F=0.12311 vs 3/26=0.11538 / F=0.12273 vs 0.125 | | yes / no |
| D=3/2, M=38 / 37 | 0.0624008 vs 0.0608108 / 0.0623950 vs 0.0625 | | yes / no |

Monotonicity in M was checked beyond the threshold, so "at most 17 / 13 / 37" is the right reading.

## 5. Union bound, exposure and the adversarial profile: CORRECT WITH CORRECTIONS

* The count of nonzero signed triples with |xyz| ≤ T is ≤ 4·Σ_{|x||y|≤T} 2T/(|x||y|) ≤ 8T(1+log T)². CORRECT.
* binom(K−r+s, s) ≥ ((K−r+s)/s)^s ≥ ((K−r)/r)^s = e^{sH}. CORRECT. With T = exp(sH − c r log r) each (allocation, s) term is
  ≤ 8e^{−cr log r}(1+rH)². Over 4^r allocations the total is exp(−(c−o(1)) r log r). For s=0, T<1 and no triple exists. CORRECT.
* **Adversarial profile.** The adversary may choose the marks, cuts and quartet after seeing the primes. This is covered,
  because the bad event is contained in a union over data that do not depend on the sample: (map from sample positions
  to {A,B,C,unused}, triple (x,y,z) with |xyz| ≤ T(s)). T depends only on s, r, K, which are deterministic. Each event in the
  union has the stated probability. **Correction:** the positions must be a *uniformly random ordering* of the drawn
  set, i.e. draw an ordered sequence of distinct primes; the law of the set is unchanged. Then, conditioned on the
  values at the non-B positions, the B-values form a uniform s-subset of the remaining K−r+s primes. If "order the
  sampled primes" means increasing order, the conditional law is constrained and the bound 1/binom is not justified as
  written. This is a one-word repair. (The notes file says "sampled prime slots", which is consistent with the repair.)
* The exceptional fraction is exp(−c_{D,M,C} r log r) for r ≥ r_0(D,M,C). For each fixed C the exceptional set
  depends on C. This is not uniform in C, but the statement claims only "for each fixed C". A larger cluster contains an
  M-subcluster, and cuts that become constant are allowed (u=0), so excluding M = 18 suffices. CORRECT.

## 6. Actual arcs, several unit classes, arbitrary C: CORRECT

Multiplying P_i by v ∈ {1, 1+i, i, i(1+i)} (v/v̄ = 1, i, −1, −i; checked) absorbs the unit. g_ij is odd, so it still divides
the new P̄'_iP'_j. Nonvanishing persists, and |t'| ≤ √2·√2·√d·Δ/2 = √d·C·N^{−1/4} for arc angle Δ = C·N^{−1/4}.
Then log d_ij ≥ W/2 − 2log C, the mean-cut estimate moves by O_C(1/W), and E log V moves by O_C(1). Both are absorbed by the
strict margin. The Plücker identity holds with the same g's, hence with the same A, B, C. No subdivision of the arc is needed.

## 7. Which circles are covered: CORRECT WITH CORRECTIONS (scope)

The theorem is about clusters whose **normalized** norm (after removing the Gaussian gcd of the retained points,
which absorbs powers of 2 and inert primes) is a good parent N of the family. Covered:

* x²+y² = N itself, for N in the good set, with any arcs of length C√R and any units.
* x²+y² = n with n = N·|G|², for clusters all of whose points are divisible by a common Gaussian integer G. Examples:
  n = 2^a·N·m² with m composed of inert primes, or a split factor common to all retained points. Here the normalized arc
  is shorter by (|G|)^{1/2}, so it is covered with the same C.

**Not covered:** any cluster on n = N·q (or N·q² etc.) in which a new split prime q is a non-constant cut, and more
generally any n whose normalized norm is not exactly a family member. The family itself (r primes all in [r^D,2r^D])
is a density-zero, rigid set of integers. Squarefreeness is used only to make each prime a single binary cut. The
bounded-exponent extension removes it for exponents ≤ E chosen adversarially on a random support, but only on that
support.

## 8. Extensions (D>1; 37 at D=3/2; exponents ≤ E): CORRECT

For D>1 we have F_M → 1/8 while the right side → 0, so a constant exists for every D>1. The exponent extension was checked:
the B-type 2|2 layers of a prime number h3−h2 ≤ E, and they all induce the same matching. The expected B-exponent is
≥ (F_M−o(1))·S over S = Σe_p layers, so the expected number of distinct B-primes is ≥ (F_M−o(1))S/E. W ≤ (D+o(1))S log r,
there are ≤ (1+3E)^r allocation states (which also encode the adversary's exponents), and T depends only on s'. The
condition (D−1)F_M/E > 3D/(2(M−1)) is correct.

## 9. Is it trivial? Yes for D>2 (a first-moment count over pairs does better)

**Proposition (referee; proved here).** Let N be a random parent of the family, C fixed, and ρ_M = ⌊M²/4⌋/C(M,2)
(ρ_2=1, ρ_3=ρ_4=2/3, ρ_5=ρ_6=3/5, …, → 1/2). If D > 4ρ_M, then except for a fraction exp(−c r log r), no arc of length
C√R carries M lattice points. In particular **D>4 ⇒ at most 1 point; 8/3<D≤4 ⇒ at most 2; any D>2 ⇒ O_D(1).**

*Proof.* Let two distinct points have separating prime set S, with |S| = k. Then ω = P_iP̄_j/g_ij is a Gaussian integer of norm
d = ∏_{p∈S}p, and z_i z̄_j = (unit)·N ω²/d. A small angle forces arg ω to lie within Δ/2 of (π/4)ℤ. ω is not real and not on a
diagonal, since d is odd and >1. Gaussian integers of norm ≤ (2P)^k in this set number ≤ c·(2P)^k·C·P^{−r/4}, with no
boundary term, because the off-axis coordinate is a nonzero integer. Each such ω fixes S's prime set. Expose the other
positions (in random order as in §5). The probability is ≤ c·C·(2P)^k P^{−r/4}/binom(K−r+k,k) ≤ c·C·(6kD log r)^k·r^{−Dr/4}.
Sum over the 2^r position sets. In an M-cluster some pair has k ≤ ρ_M r, because Σ_pairs k_ij = Σ_p u_p ≤ r⌊M²/4⌋. The
total is exp(−(D/4 − ρ_M − o(1)) r log r). ∎

Comparison (`thresholds.py`, smallest excluded M):

| D | 1.1 | 1.5 | 2.0 | 2.1 | 2.2 | 2.5 | 3 | 4 | 10 | 100 |
|---|---|---|---|---|---|---|---|---|---|---|
| Astra (III.7) | 134 | 38 | 26 | 25 | 24 | 22 | 20 | 18 | 15 | 14 |
| pair count | — | — | — | 21 | 11 | 5 | 3 | 3 | 2 | 2 |

So the quartet argument is the better tool only for 1 < D ≲ 2.08, and it is the only tool here for 1<D≤2. The headline
"17 at D=4, 13 for large D" should be replaced by "2 at D=4, 1 for D>4 (elementary); the quartet argument gives
the first constant bound for 1<D≤2, e.g. 37 at D=3/2". Both arguments lose the same entropy factor k! ≈ r^k relative to
the Poisson heuristic, which predicts ≤ 1 point for every D>0 (expected close pairs ≈ 4^r N^{−1/4} → 0).

## 10. Novelty against the record

grep for "almost all", "typical", "random (prime|parent|support)", "first moment", "density" in `docs/claude_review`
finds no statement about random parents of this shape. The related items are:

* Bourgain–Rudnick's almost-all separation (overview §6, `lit18`). It concerns density-one sets of n and does not cover this
  density-zero family.
* The almost-all sector caps of round 10 (Thm B, Prop C). These concern analytic inputs, not parents.
* The round 21 first-moment threshold for *constructing* fakes.

Astra's own III.7 "almost every parent excludes every assignment" is the Paley-specific precursor in this same round.
The quartet-exposure bound is new to the record. The phenomenon "random parents have tiny clusters" is not deep:
the pair count above gives it for D>2.

## 11. Numerical sanity experiment (`experiment.py`, independent code)

I enumerated **all** 4·2^r lattice points of random parents, using exact Gaussian products. Points were sorted by angle, and
arc membership was tested by the exact integer chord test |z_i−z_j|⁴ ≤ C⁴N (a superset of arc ≤ C√R).

| r | D | K | parents | max points on an arc, C=1 / 10 / 100 |
|---|---|---|---|---|
| 10 | 3 | 67 | 20 | 1 / 1 / 2 |
| 12 | 3 | 107 | 20 | 1 / 1 / 1 |
| 14 | 3 | 160 | 10 | 1 / 1 / 1 |
| 16 | 2.2 | 33 | 10 | 1 / 1 / 1 |
| 14, 16 | 4 | 1747, 2837 | 10, 6 | 1 / 1 / 1 |
| 18 | 2, 3 | 24, 314 | 4, 4 | 1 / 1 / 1 |

Close pairs occur in clumps. The only close pairs (32 of them, at r=10, C=100) came from a single near-axis ω shared by many
point pairs. This explains why the total pair counts sit below the naive Poisson count in mean-of-small-samples. The
experiment is consistent with the theorem and with the stronger elementary bound (typical value 1). As expected at
these sizes, it cannot distinguish 17 from 2.

## Verdicts

1. Definitions, (III.1), A,B,C, uniqueness of the forced subset: CORRECT.
2. Pair mean, middle-class probability, F_M by convexity: CORRECT.
3. E log V ≤ 3W/(2(M−1)): CORRECT.
4. (III.7) arithmetic (18 yes / 17 no at D=4; 14/13 large D; 38/37 at D=3/2): CORRECT.
5. Union bound and exposure against an adaptive profile: CORRECT WITH CORRECTIONS (random ordering of the sample is
   required; the exceptional set depends on C and M).
6. Passage to actual arcs, several units, arbitrary C: CORRECT.
7. Coverage: CORRECT WITH CORRECTIONS (only clusters whose normalized norm is a good family member; new non-constant
   split cuts are not covered; the family has density zero).
8. Extensions D>1 and e_p ≤ E: CORRECT.
9. Triviality: the D ≥ 4 headline (17, 13) is superseded by an elementary pair first-moment bound (2 at D=4, 1 at D>4,
   constant for all D>2). The nontrivial content is 1<D≤2.
10. Significance: peripheral. No change to the growth rate, the uniform bound, Questions A/B or route 5.
