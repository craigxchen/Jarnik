# Referee report, block "char2": Part VII, section VII.6 (characteristic two)

Referee scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_char2/`
(`gf.py` is my own F_{2^k}, F_2[x]-multivariate and F_4[t] arithmetic; there is no sympy. `symbolic.py`, `fibres.py`, `genuine.py`,
`genuine_big.py`, `dim4.py`, `dim4count.py` and `gcdtest.py` each have a matching `*.out`).

**Overall verdict on VII.6: CORRECT.** One correction about which differential matters, one simplification, and one
strengthening are below; the strengthening is supported numerically but not proved. The claim stands: the
bounded-term ("sparse") theorem of VII.2 fails over F_4(t)/F_2(t) for every delta < 1/30, with one term budget for all
members. The covering-family inequality (VII.2) is the step that fails, and it is in fact plainly false in
characteristic 2. Novelty against the record is modest. Significance for the integer goal is low: this is a
calibration result.

---

## (a) Algebra: subspace polynomials, common norm, contact. CORRECT

* c = sum x_j lies in x_j + H_j, so c is not in H_j, and V = H_j ⊔ (c+H_j). Hence F_j(Y)F_j(Y+c) = P_V(Y)/(Y(Y+c)).
  The recurrence P_{H+<v>} = P_H^2 + P_H(v)P_H gives additivity (as in Part V).
* **Checked symbolically** over F_2[x_1..x_5,Y] (`symbolic.out`). Each F_j = P_j/Y has Y-exponents in {15,7,3,1,0}.
  max_{i<j} deg_Y(F_i - F_j) = **7** exactly, against 15. The gcd of the F_j is 1, because the five H_j meet in 0.
* Split model with generic constants: h = 8 and W = d_+ + d_- = 15 + 15 = 30, since 1/F_j is proportional to F_j(t+c). So
  h/W = 4/15 = 1/4 + 1/60, i.e. h/D = 8/15 = 1/2 + 1/30.

## (b) Dominance via Moore minors. CORRECT

* **Symbolic checks.** y_j := F_j(0) = (product of the 15 nonzero elements of H_j) equals the 4x4 Moore determinant
  deleting column j, checked by Leibniz expansion over F_2 for all j. sum_j y_j x_j^{2^a} = 0 holds identically for
  a = 0,1,2,3. As a sanity check, at a = 4 the sum is the 5x5 Moore determinant, which is the product of the 31 forms.
* **Her inverse-Frobenius argument is correct.** Taking 2^a-th roots of row a gives the linear system
  sum_j y_j^{2^{-a}} x_j = 0. Raising its rows to the 8th power and reversing them gives the Moore matrix of y. Its 4x4
  minors are nonzero because the y_j are F_2-independent; her GL_5(F_2) column-operation argument for this is
  correct, since Frobenius commutes with F_2-linear column operations. So the kernel is the line through x, and
  [y] determines [x].
* **Independent exhaustive check over F_32** (`fibres_dim5.out`). There are 322,560 projective points with x_1 = 1 and
  {x_j} an F_2-basis of F_32. They give **322,560 distinct images [y]**, so every fibre has size 1. The map is injective
  on rational points, as the argument predicts.
* Normalising to c = 1 is legitimate: all minors have degree 15, and the c = 1 slice meets every scaling orbit with
  c ≠ 0.

## (c) Differential rank. CORRECT, with a correction of emphasis

* **Symbolic check.** sum_j x_j^{2^a} ∂y_j/∂x_k = 0 holds identically for a = 1,2,3 and all k. At random points of
  F_{2^20} the affine Jacobian of y has rank 2, and the projective rank of d[y] is **1**. This is exactly her claim.
* **Correction.** The covering argument (VII.2) needs the rank of the *whole* family (x, Y) -> [F_1(Y):...:F_5(Y)],
  source derivative included, not just the Y = 0 slice. I computed it at 4 random points of F_{2^20}^6. The
  projective rank is **2**: rank 2 from the x-directions alone, and rank 1 from the Y-direction alone. This is still
  < 4, so her conclusion holds, but "rank at most one" should be stated for the slice only.
* **No conflict with the record.** The R31-V ruling ("Part V is not the inseparable failure; every ratio has nonzero
  t-derivative", reconcile (d)) concerns single curves. Here, too, each curve is separable (Y-rank 1). What is
  inseparable is the evaluation map of the family onto T: its image has dimension 4 but its rank is 2. The two
  statements are about different maps and are compatible.

## (d) Genuine conjugate-norm model. CORRECT; the contact bound can be sharpened

* Setup: c = 1, x_j = t^{n_j} (j ≤ 4), x_5 = 1 + sum, Y = t^L + θ with θ^2 + θ + 1 = 0. Then conj(Y) = Y + 1 = Y + c. The F_j have
  F_2[t] coefficients, so conj z_j = F_j(Y+1) and the norm is common. Distinct positive n_j give F_2-independence: a
  sum involving x_5 has constant term 1. The degree is exactly 15L for L > B.
* **Sharpening.** The coefficient e_8 of Y^7 has t-degree ≤ 8B, so max deg(z_i - z_j) ≤ 7L + 8B and **h ≥ 8L - 8B**.
  Her bound 8L - 16B is valid but weaker.
* **Exact F_4[t] reproductions** (`genuine.out`, `genuine_big.out`). For each instance I checked that all norms
  z·z̄ are equal and lie in F_2[t], that the marks are distinct, the full gcd over F_4[t], D after gcd removal, and
  h = D_raw - max deg(z_i - z_j), which is unchanged by the gcd.

| n | L | raw deg | gcd deg | D | h | h/D | 8L-8B | 8L-16B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| (1,2,3,4) | 100 | 1500 | 7 | 1493 | 776 | 0.5198 | 768 | 736 |
| (1,2,3,4) | 300 | 4500 | 8 | 4492 | 2376 | 0.5289 | 2368 | 2336 |
| (1,2,3,4) | 700 | 10500 | 7 | 10493 | 5576 | **0.5314** | 5568 | 5536 |
| (1,4,6,9) | 60 | 900 | 15 | 885 | 428 | 0.4836 | 408 | 336 |

  h/D tends to 8/15 = 0.5333 from below at rate O(B/L), as claimed. For small L/B the ratio is below 1/2, so the
  statement genuinely needs "L large after n is fixed". That dependence is legitimate in a counterexample.
* **Primitivity.** The tuple gcd is small but **nonzero** (degree 7 to 15), and she removes it. That matches the record's
  notion (gcd of the tuple = 1; expert_note §3, "Setting"). Individual marks are **not** coprime to their conjugates:
  deg gcd(z_j, z̄_j) ranges from 6 to 28. So single marks are divisible by nonconstant elements of F_2[t]. The
  record's definition allows this, but "primitive" should not be read as "each mark primitive".
* Monomials per mark do not depend on L (170 to 234 for n = (1,2,3,4)), far below 6^15.

## (e) Term budget and density. CORRECT

* **Term budget.** Each factor t^L + θ + v has at most 6 monomials (at most 4 powers t^{n_k}, one constant, t^L).
  So b = 6^15 for numerator and denominator, uniformly in n, L and the gcd. This is crude but finite and uniform; no
  template-dependent budget is involved.
* **Density.** A nonzero P on the torus, even over the algebraic closure of F_2(t), pulls back to a nonzero
  Φ(x_1..x_4, Y). This uses dominance over F̄_2 and base change. Products of infinite sets {t^n} are Zariski dense, so
  one fixes n with Φ(t^n, Y) ≢ 0 (avoiding diagonals and denominators). Then all large L avoid its finitely many
  zeros and poles and also exceed the contact threshold.
* **Quantifiers.** The requirement is "for every proper Z and every W_0 there is an escaping member with
  h ≥ (1/4+ε)W". n may depend on Z and L on (n, δ). There is no "eventually along a subsequence" step here.
* **Simplification.** VII.2 is stated for arbitrary q_j ∈ C(t), split model included. Its characteristic-2 analogue
  is already refuted with **b = 5**: take q_j = F_j(t^k)/F_1(t^k) with generic constant x ∈ F̄_2^5. Then W = 30k,
  h = 8k, and density follows from (b). The genuine F_4/F_2 version is needed only for the A_5-type statement with
  the norm structure.

## (f) Which step fails. CORRECT WITH CORRECTIONS

* **The covering inequality (VII.2) is false in characteristic 2.** The split family above (finite type, W_0 = 30,
  h_0 = 8, dominating T) has 4h_0 + 1 = 33 > 30. The failing line is "all derivatives together span the four-space
  by generic smoothness". The family rank is 2, so the 4x4 log-derivative determinant vanishes identically. Every
  other line of that proof (pole and zero counting) is characteristic-free, and p | h would only add zeros at
  infinity.
* **The gcd estimate is characteristic-free and true in characteristic p.** Each of its steps (Dirichlet
  approximation, the resultant in the ideal (A, C) after the weight substitution, bounded support box,
  injectivity of the weights, the root count O(L/Q)) uses no characteristic.
  * gcd(t^a - 1, t^b - 1) = t^{gcd(a,b)} - 1 over every field.
  * Sparse high powers such as (1+Z)^{p^k} = 1 + Z^{p^k} are absorbed into the exponent vector during the fixed-torus
    passage, so they do not violate the "fixed box" hypothesis.
  * The one char-p caveat is that coprimality must be checked in characteristic p, e.g. 1 + Z_1^2 + Z_2^2 = (1+Z_1+Z_2)^2
    over F_2.
  * Sanity test (`gcdtest.out`): 7 fixed-box pairs over F_2, n = (L, aL+b), L = 64..2048. The maximum gcd degree stays
    in [0, 16] and does not grow.
* The passage to a fixed torus and the frozen-family step are also characteristic-free. So in VII.2 the
  characteristic enters **only** through (VII.2). This is consistent with VII.2 being true in characteristic 0, but it
  does not verify VII.2's characteristic-0 steps; those belong to the gcd and uniform blocks.

## (g) Ternary analogue (VI.4). Bounded terms YES for a fixed template; density NO; δ far below 1/30

* For fixed (q, k, forms, signs), (VI.5) gives F_j = P_V + 2c^{q-1}A - 2c^{(q-1)/2}A^{(q+1)/2}, with A linearized
  (k terms).
  * At q = 27, k = 3: A^{14} = A^9·A^3·A^2 has 54 X-exponents, so F_j has ≤ 58 X-monomials.
  * Its coefficients lie in F_9[t] with degree ≤ (k-1)n.
  * After X = t^L, numerator and denominator of every ratio F_j/F_1 have a number of monomials bounded independently
    of L.
* Contact: h/D -> (n-h)/(2(n-g)) = 1/2 + (g-h)/(2(n-g)). For the recorded five-mark template (n = 19683, h = 729,
  g = 1301) this is **1/2 + 0.01556**, below 1/2 + 1/30.
* **Not dense.** Every ratio curve is t -> R(t, t^L) for one fixed rational map R on the (t, X)-plane. So all members
  lie in a proper closed subset of T of dimension ≤ 2, and they do not refute a sparse theorem.
* The record's dense F_3 refutation (Weil-balanced subtuples, q -> ∞) has term counts that grow with q.
* A dense ternary bounded-term family would need generic generators x_b in place of t^b, with k ≥ 4 (n = 27^4) for five
  marks. It would also need a characteristic-3 dominance proof. **Not established.** Even then the orientation
  bound g ≥ n/16 gives δ ≈ 0.0136 < 1/30. So VII.6's "no term-budget assertion for the ternary examples" is the right
  caution.

## Referee extension (CONJECTURE, strongly supported): the same mechanism reaches every δ < 1/14

* Construction: take V = span(1, x_2, x_3, x_4), four-dimensional, and five of the eight hyperplanes H not containing
  1. Set F_H = P_H(Y)/Y, of degree 7.
  * deg_Y(F_H - F_H') ≤ 3 (symbolic, `dim4.out`), so h/D -> 4/7 = 1/2 + 1/14 in the genuine model.
  * The split-model term budget is 4.
  * It is still below the characteristic-free pair bound 3/5.
* **Exact F_4[t] instance** (n = (3,4,11), L = 200): D = 1394 and h = 774, so **h/D = 0.5552 > 8/15**, with norms
  checked (`genuine.out`).
* **Density is not proved.** The Jacobian has projective rank 2 here as well, so the Jacobian test cannot decide
  dominance.
  * Evidence: (x_2,x_3,x_4,Y) -> T is **injective on all F_32-points (322,560) and all F_64-points (9,999,360)** for
    two hyperplane selections (`fibres_dim4.out`, `dim4count.out`).
  * Images/q^3 grows from 9.8 to 38.1, which fits image dimension 4.
  * The Lachaud–Rolland-type bound for a degree ≤ 343 threefold does not exclude dimension 3 at q = 64. So this is
    evidence, not proof.
* If dominance holds, the characteristic-2 sparse failure extends to every δ < 1/14. That matches the "rank-4
  subtuple" value δ_S -> 1/14 noted in ref5_reconcile, where those subtuples were set aside as possibly special. A
  generic-injectivity proof analogous to her Moore-kernel argument would settle it.

## Standing cautions

* **Primitive height.** D is computed after removing the full tuple gcd, exactly. W = 2D also needs the conjugate
  tuple to be gcd-free, which holds by symmetry.
* **Template-dependent thresholds.** The L threshold depends on n. That is legitimate in a counterexample. The budget
  b is uniform.
* **Uniform exceptional set.** The negation is proved against every proper closed Z, even one over the algebraic
  closure of F_2(t).
* **Fixed-S Subspace, drifting phase.** Not used.

## Novelty and significance

* **Already in the record (R31 V/VI):**
  * F_2(t) and F_3(t) analogues of A_5(δ<1/30) fail densely, via Part V/VI.4 clusters with m -> ∞ and unbounded
    complexity.
  * (IV.2) fails in characteristic 2.
  * Subspace polynomials and the Artin–Schreier translation mechanism.
* **New in VII.6:**
  * A fixed-size five-mark family, with V generic and Moore minors, that is dense with a **uniform term budget**.
    This separates the sparse theorem's characteristic-0 proof from characteristic 2.
  * The identification of the exact failing step: the inseparable evaluation map, i.e. failure of generic smoothness.
* Moore determinants and subspace (Ore) polynomials are classical. The novelty is in the application and is modest.
* **Significance for the goal: low.** It shows that the separability hypothesis in (VII.2) is essential and that the
  bounded-term restriction cannot be made characteristic-free. It says nothing about the unrestricted
  characteristic-0 question or about integers.
