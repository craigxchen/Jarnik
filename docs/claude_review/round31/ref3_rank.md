# Referee report, round 31, block "rank": astra_free_III/answer.md, sections III.4 and III.5

Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref3_rank/` (`fourcycle.py`, `paley_rank.py`). Astra's scripts `iii_tri_balanced.py` and `iii_rank_extension.py` were also rerun; both reproduce their stated output.

## Summary verdict

Every proof in III.4 and III.5 that I checked line by line is correct. The four-cycle example is genuine, and I verified it with independent code. There are only minor notational corrections. All results in the block are **obstructions to methods**. None of them changes the growth rate, the uniform bound, or the status of Questions A, B or route 5. The main content is a sharp, clean reason why the coordinator's routes (a) and (b) fail:

- **(a), rank.** Any multiplicative coset holding the extension data of s third vertices has rank at least about 3s/16 to s/8. The bound is linear in s, so every uniform unit-equation count, being exponential in the rank, is vacuous here.
- **(b), C4 exclusion.** Excluding monochromatic four-cycles can neither be forced by the edge-density count nor justified at the four-point level.

These results are **peripheral to the goal**, but they are useful for closing off routes. The record already rated single unit equations and ESS as "vacuous" and "exponential in rank" (overview line ~1281; lit18/tools.md row 14 and §2.5). What is new is the precise, cluster-specific linear rank lower bound and the explicit C4 example.

No unit-equation theorem is *used* in any proof. Evertse–Schlickewei–Schmidt (ESS) is cited only as the benchmark that the linear bound beats. The cited statement is correct: for an n-variable equation a_1x_1+...+a_nx_n=1 with x in a rank-r subgroup Γ of (K^*)^n, char K = 0, the number of nondegenerate solutions is at most exp((6n)^{3n}(r+1)), uniformly in the coefficients (Ann. of Math. 155 (2002) 807–836). For n = 2 the Beukers–Schlickewei bound is 2^{8(r+1)}; some secondary sources quote 2^{9(r+1)}, and the difference is immaterial here. Sources: [Annals page](https://annals.math.princeton.edu/2002/155-3/p04), [arXiv math/0409604](https://arxiv.org/pdf/math/0409604).

## III.4

### (A) Exact edge version: s ≤ 8r_0+4 — CORRECT (in fact strict: s < 8r_0+4)

Checked line by line.

1. **Local formula for the divisor split.** At prime p, v_p(a_x) = (|a_i−a_j|+|a_i−a_x|−|a_j−a_x|)/2. A case check in the three positions of a_x relative to [u,v] gives min(a_x,v)−min(a_x,u) = z_xp when a_i is the lower anchor, and L_p − z_xp otherwise. This is affine in z_xp, with the map fixed by the anchors.
2. **Projection.** c_x(p) = (2z_xp−L_p)√(log p/(W L_p)) is exactly the inner product of the threshold-sign row f_x with the normalized block indicator: the sum over the block of f_x is 2z−L. So C = FΠ^T with Π a partial isometry.
3. **Operator norm.** ‖F‖² ≤ 2 holds because the row Gram G is PSD with unit diagonal, and G_xy = 1−2w(a_x−a_y)/W < 0 by pair separation. Hence I−G is nonnegative with top eigenvalue ≤ 1, so by Perron–Frobenius its spectral radius is ≤ 1 and spec G ⊂ [0,2]. The same holds for every subset of rows and survives centering.
4. **Lower bound on spread.** Draw x, y independently and uniformly. Then E(distance on retained layers)/W ≥ (1−1/s)/2 − (1−D)/2, since each discarded layer has disagreement probability 2q(1−q) ≤ 1/2. Cauchy–Schwarz gives |c_x−c_y|² ≥ 4X²/D, Jensen gives the expectation, and (D−1/s)²/D ≥ D−2/s. The pair identity then gives a total centered squared norm ≥ sD/2 − 1.
5. **Upper bound from the coset.** A rank-r_0 coset in K^* (any field: valuations at primes above p are e(P|p)v_p) has affine valuation image of dimension ≤ r_0, so the centered c_x lie in a linear space V with dim V ≤ r_0. Hence the total is ≤ ‖C‖²·dim V ≤ 2r_0.
6. **Conclusion.** D = w(a_i−a_j)/W > 1/2 by pair separation, so s < (4r_0+2)/D < 8r_0+4. The a_x/b_x = a_x²/d_ij variant follows in the same way.

The proof uses **only pair separation w(a_x−a_y) > W/2**. As Astra says, it therefore holds verbatim in the angle (fake) model. It is a combinatorial statement about profiles and is F1-transparent.

### (B) Three-anchor version (III.4): s ≤ (16/3+o(1))r_0 + O(1) — CORRECT WITH A NOTATIONAL CORRECTION

- D > 3/4 − 4/√m follows from the identity d_ij+d_jk+d_ik = 2·range, together with e_p log p ≤ 4W/m and the discarded primes p < √m (II.4).
- Selecting the extreme-anchor components recovers z_xp with error E_xp = difference of two cofactor valuations. |E_xp| ≤ 2B for p ≥ √m.
- The squared error is ≤ (4/W)Σ E² log p ≤ 24B² log m / W.
- W ≥ log((m−1)!), since m affinely independent rows need m−1 primes.
- The perturbation step: √(sD/2−1) ≤ √(2r_0) + √(sε) with sε = O_B(1) for s ≤ m.

The corrected statement is s ≤ (4/D)r_0 + O_B(√r_0) + O_B(1) = (16/3 + O(m^{-1/2}) + o(1))r_0 + O_B(1). **The additive O(1) depends on B** (Astra writes o_B(1) in front of r_0 and an unqualified O(1)). This is harmless.

### (C) One rational coordinate: s ≤ (8+o_B(1))r_0 + O_B(1) — CORRECT

- min(a_x,a_i)−min(a_x,a_k) = ±(clipped coordinate), with the sign fixed per prime.
- D > 1/2 − o(1).
- Varying rational multipliers of height ≤ m^{B'} add errors of the same type.

### (D) Moving offsets in (II.1) — CORRECT

- The term (b_x/a_x)(v_x/u_x)² has valuation −2v_p(a_x) + const + 2v_p(v_x/u_x). This is the clipped coordinate up to a fixed affine map, with errors from offsets ≤ m^B.
- The one-coordinate proof applies and gives factor 4/D → 8.

### (E) Consequence for ESS — CORRECT; peripheral

- Mean log|t_ij| ≤ W/(4(m−1)) − log 4 = O_b(log m).
- By Markov, at most 1% of edges have |t| > m^B.
- A random triple and a random x lose at most 3% to bad edges, so some triple has ≫ m good extensions.
- Then r_0 ≥ (3/16−o(1))s − O_B(1) ≫ m. Equivalently, a rank-o(log m) coset holds o(log m)+O(1) extensions.

ESS, Beukers–Schlickewei and Amoroso–Viada are all exponential in the rank, so they say nothing here. The record's general caution is that unit-equation counts are known to be near-sharp exponentially in rank in general (Konyagin–Soundararajan lower bounds, cited in lit18). The new point is that this rank must be linear in the number of solutions **in every grouping**, not only for the atom group. The caveat Astra states is right: the result concerns this specific unit equation in cross ratios, not other additive systems.

### (F) COMPUTED anchor-rank check — CORRECT (reproduced independently)

- I built the ordinary Paley double-core rows myself: column pair (a, ±), bit [χ(a−x)=1] off the diagonal and ± at x=a. This agrees with ii_paley.double_core for q ≡ 3 mod 4.
- For q = 7, 11, 19, every one of the 35, 165 and 969 anchor triples gives projected affine rank 3, 7 and 15 = q−4, computed with exact rational rank (sympy).
- I also checked the exact edge inequality numerically on all pairs (equal-weight model): sD/2−1 ≤ total centered norm ≤ 2·rank and ‖C‖² ≤ 2 hold, with minimum slack 7.6, 14.5 and 28.5.
- Astra's own script reproduces its stated output.

### (G) Star identity and the fixed-weight caveat (iii_star_rank.md, iii_star_targets.md) — mostly CORRECT; one part NOT CHECKED

- **Coset to marks: CORRECT.** If the word vectors of k-tuples lie in a rank-r coset (mod torsion), then all marks used lie among ≤ k(r+1) marks. The proof uses injectivity of zero-sum formal vectors → valuations (real affine independence, II.4) and an affine basis of ≤ r+1 tuples. The comparison with ESS is trivially correct.
- **Triangle squareclasses: CORRECT.** A coset of a rank-r subgroup of Q_{>0}^{*6} meets ≤ 2^r square-class images, and distinct triangles have distinct parity images by II.4. So ESS gives no saving.
- **Place-uniformity example: CORRECT.** In the absolute normalization over Q(i), the infinite place with forms X, Y gives 2j·1. The finite forms X±iY over all places dividing 4j²+1 give (4j²+1)^{-1/2} each by the product formula. The product is 2j/(4j²+1) < 1/(2j) = H(2j,1)^{-1}. This gives M distinct lines with fixed forms and a growing place set.
  - This is the standard dependence of Subspace counts on |S|, already recorded in the record (lit18 row 13: Evertse's (2^{60n²}δ^{-7n})^s). It is not new information.
- **Fixed-center bound: CORRECT.** At most 576/δ⁴ oriented edges with gain w(v)−w(v−b) ≥ δW, and O(δ^{-5}) on a ray. I checked:
  - the matching bound 16/δ², from rounding cost 2W/√h, the inequality |d|−|d−b| ≤ bd/e for d ∈ {0,±e}, and imbalance 2W/√h;
  - the star bound 9/δ², using that |z|+|z−e| is minimal on [0,e];
  - the impossibility of both orientations, from w(v−b)+w(−v−b) ≥ 2w(v);
  - the vertex-cover count;
  - the ray net: the gain is 1-Lipschitz in b, successful centers have w(b) ≤ 2W, and 4/δ net points suffice.
- **NOT CHECKED in detail:** the "twisted height" translation (the archimedean ½log|st| ≥ W/8 step and the centered local adjustments), the fixed-higher-words Walsh argument, and the Auerbach-basis √m parameter bound. They are peripheral, and I give no verdict on them.

## III.5

### (H) Density threshold — CORRECT

- **Multiplicity bound.** If the largest multiplicity is f, the j-th smallest label is ≥ j/f. Then Σ log ≥ log E! − E log f ≥ E log E − E − E log f, so f ≥ E e^{−1−H}.
- **Mean log label.** H ≤ W/(4(m−1)) + O(1), because log t ≤ L_ij/2 − log 4 and mean L ≤ W/(2(m−1)).
- **At W = 3m log m.** f ≥ m^{5/4−o(1)}. This forces C_10 by Bondy–Simonovits (ex(n,C_{2k}) = O(k n^{1+1/k}); for k = 5 this is n^{6/5}). It does not force C_4 (projective-plane incidence graphs give ~n^{3/2} edges) and not even C_8 (n^{5/4}, borderline).
- **C_10 calibration.** No monochromatic C_10 would need H ≥ (4/5)log m − O(1), i.e. W ≥ (16/5)m log m − O(m). That would give the constant 5/8 in place of 2/3. This is purely hypothetical: nothing in the record supports a C_10 exclusion.
- **Budget remark.** "For b > 8 the budget permits all labels distinct" is a correct budget statement (2 log m < (b/4) log m), not a construction.

### (I) Parity and two equal offsets — CORRECT

- In one unit class, z_i = g'Ā and z_j = g'A with A primitive of odd norm. So z_j − z_i = 2it·g' and v_{1+i} = 2 + 2v_2(t).
- With a fixed valuation k, the vertices of each edge share a class mod (1+i)^k and differ mod (1+i)^{k+1}. Each component therefore alternates between two lifts, so the graph is bipartite.
- **Two equal offsets t from one anchor.** N(gcd(A_j, A_k)) divides Im(Ā_j A_k) = t(s_k−s_j). It is odd and coprime to t, because gcd(s_j,t) = 1. The reduced αβ̄ is primitive: I checked all four ways a prime could divide both it and its conjugate. Its imaginary part is t(s_j−s_k)/g, and s_j ≡ s_k mod 2 because s ≢ t mod 2. So the third edge has 2t | cofactor.
- **The (II.1) specialization.** With u = v = t: c = (a+b±2s)/h² with h = |w|/t an (even) integer. From ac − t² = ((a±s)/h)², integrality gives h | a ± s. Correct.
- Consistency with the four-cycle example: anchors 0 and 1 each have two edges of cofactor 2, and the opposite edges 12 and 03 have cofactors 560 = 4·140 and 7344 = 4·1836.

### (J) The balanced four-cycle — CORRECT (verified exactly with independent code)

Independent checks in `fourcycle.py` (sympy for the quartics; my own Gaussian arithmetic for the points):

- **Quartics.** All six B_ij are monic quartics in n with constant imaginary parts (+19966500, +16110900, +14871600, −3855600, −5094900, −1239300), matching the stated absolute values.
- **Period.** lcm(|b_01|,|b_02|,|b_13|,|b_23|) = 86434978500. At n0 = 79658575767 the gcds are 9983250, 8055450, 2547450, 619650, with b/g = 2 and Re/g odd. Both depend only on n mod b (Re/g mod 2 is determined by Re mod 2g = b), so verification at one residue suffices, as claimed.
- **Points, at five points of the progression** (t = 0, 1, 7, 10^6, 10^20 periods):
  - The common Gaussian gcd has constant norm 859172194859062500. After dividing it out, the gcd of the four points is a unit, the four points are distinct, and all have norm N (135 digits at t = 0).
  - N ≡ 1 mod 4, and all prime factors below 2·10^5 are ≡ 1 mod 4.
  - For each edge I computed the primitive half-angle A = s+it directly from z_i z̄_j and checked that z_i Ā = z_j A **exactly**: the ratio is A/Ā with no unit, so the four points are in one unit class.
  - Offsets are (2,2,7344,560,2,2) at every tested n. All are even, so the edges are primary-class. The K_{2,2} edges are 01, 02, 13, 23.
- **Arc.** Angular width × N^{2/7} = 300.13, constant, so the width is ≍ N^{−2/7}. The arc constant (arc length)/√R, with R = √N, is 4.6·10^{−3} at t = 0 and decreases to 0 along the progression (< 1/128 throughout). These are genuine clusters on (sub)critical arcs.
- **Balance.** The seven columns are the seven nonempty subsets of marks {1,2,3}, so each nonconstant four-mark cut occurs once, with weight log N(H_k) − O(1) = 2 log n + O(1). The bounded-factor argument (primes common to two linear values divide fixed root differences; remove min(valuation, fixed bound)) is correct.
  - **Quantitative caveat:** the O(1) is large. As a crude upper bound, stripping all powers of the 22 "bounded" primes costs 18–31 nats of log N. At the certificate residue itself, the column weights after this stripping range from 28.6 to 45.8 against log N/7 = 44.3, so that instance is far from balanced. Balance holds only asymptotically along the progression; at t = 10^20 the spread is 123.7–138.0 against 136.6.
  - The theorem as stated ("arbitrarily large N, weights (log N)/7 + O(1)") is correct.
- **Smallest cofactor.** In the primary convention a primitive A with A/Ā in the primary class is ± a product of primary factors, so it has even imaginary part, and 2 is the minimal cofactor. Correct.
- **Calibration.** On five marks with 16 unoriented cuts, two of them (the constant cut and the cut isolating the fifth mark) restrict to the constant cut on four marks, so a balanced four-subtuple carries 14/16 = 7/8 of W. The parent width e^{−W/4} becomes H^{−2/7} for the four-point norm H = e^{7W/8}. The example's N^{−2/7} matches. Correct.

**Significance of the four-cycle:** it shows that "no monochromatic C_4" cannot be derived from four-point consequences. It does not bear on whether a large parent can host many such cycles. At b = 3 the forced density m^{5/4} does not force C_4 anyway, so C_4 exclusion was never the operative target. The example is peripheral; its value is that it removes a tempting false lemma.

## Novelty and significance (blunt)

- **Novelty.**
  - The linear rank lower bound (8r_0+4 exact; 16/3 and 8 asymptotic) is new relative to the record. The record contains only the generic statement that unit-equation bounds are exponential in the rank.
  - The parity / 2t lemma and the C_4 example are new.
  - The place-uniformity example is a restatement of known |S|-dependence.
- **Significance.**
  - None of these results changes the growth rate (still (2/3+ε) log R/loglog R), the uniform bound, or the status of Questions A or B or route 5.
  - They are clean negative results that retire the coordinator's suggestion (a): organise a unit equation so that the rank is o(log m). For the cross-ratio / triangle unit equations at fixed anchors, this is now **provably impossible**: the rank is ≥ cs.
  - They also show that (b), cofactor pigeonhole, can at best force long even cycles (C_10), whose exclusion would only give the constant 5/8 and is itself unsupported.
  - All of this is F1-transparent (it holds for angle fakes with pair separation).
