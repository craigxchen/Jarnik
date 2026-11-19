# Judge report, round 31: Part VII of the Question B thread (Astra, "Five marks over a characteristic-zero function field")

Source: `astra_lemma_VII/answer.md` from line 1565 (`# Part VII`) to the end, with `vii_pell_check.py`, `reply7.txt`, `final7.txt`.
Referee reports in this directory: `ref7L_gcd.md`, `ref7L_uniform.md`, `ref7L_circle.md`, `ref7L_pell.md`, `ref7L_char2.md`,
`ref7L_search.md`, `ref7L_meaning.md`. Judge script: `README_VIIL/pell5.py` (output `README_VIIL/pell5.out`). I also re-ran the
search referee's `ref7L_search/verify_d11.py`.
Claims of my own are labelled **[J]** and are unrefereed.

## 0. Verdict (five lines)

1. Nothing in Part VII is false. No referee returned FALSE on a claim of Astra. The two GAP verdicts (ternary term budget; D = 9 real model) are on statements Astra explicitly does **not** make, and I uphold them as gaps in those non-claims.
2. The sparse theorem VII.2 is a **theorem**, after one wording correction (frozen family with **exact** supports). I checked both steps Astra asked to have prioritised myself [J]: the moving-coefficient gcd estimate and the one-exponent-vector passage to a common exceptional set. One proper Q-closed Z_{b,eps} for all degrees, coefficients, exponents and common factors is legitimate. Every "eventually" is either along a countable diagonal subsequence in a proof by contradiction or over a finite set fixed in advance.
3. (VII.4) 20h <= 9D + 18r is a theorem. The VII.4 structure facts at 2h > D are theorems. The Pell families of VII.5 are correct: I reproduced D = 9, h = 4; D = 17, h = 8; D = 25, h = 12 exactly with my own code [J], with h = 4n for **every** choice of distinct seeds (correction). The characteristic-2 failure VII.6 is correct; the rank statement is corrected from "<= 1" to "2 for the family, still < 4".
4. New and correct relative to the record: VII.2 (as a statement; its gcd lemma is a special case of Guo-Sun-Wang and its torus passage is standard Zannier-type reduction); (VII.4); the four-independence lemma; dense five-mark Pell families at W = 4h + 2; and the consequence Astra did not draw, alpha^ess_e(-K) <= 4 over Q. From the referees: the first characteristic-0 five-mark configurations with 2h > D (rigid, D = 11, 13); H_5(7) = 3 exactly; (VII.2) is false in characteristic 2.
5. For Question B over Q, the growth rate (2/3+eps) log R/loglog R and the uniform bound M(C): **no change**. The unrestricted characteristic-0 five-point question is open, as Astra says. The only Diophantine output is a calibration: the known window for alpha^ess is now [10/3, 4] (previously the upper end was open), and Question B asks for > 15/4. Close the thread as internal work.

## 1. Per block: referee-corrected statements

### 1.1 The gcd estimate (referee: gcd; attribution by meaning)

**Proved (statement as corrected).**
> Let n = n_l in Z^r with L = ||n||_inf -> infinity, such that for every fixed nonzero u in Z^r, n_l . u != 0 for all large l. For each l let F_1, ..., F_k be Laurent polynomials with exponents in a fixed box B. Their coefficients may depend on l, but they must have no common nonunit factor in C[Z_1^(+-1), ..., Z_r^(+-1)]. Then the total common vanishing on C^*, with multiplicities, of F_1(t^n), ..., F_k(t^n) is o(L). The o(L) is uniform in the coefficients. It has no rate, and can decay as slowly as ~1/log L.

* The Dirichlet step holds with ||v||_inf <= b' (sharper than b'+1).
* Characteristic-free.
* The cited Amoroso-Leroux-Sombra slides are decorative.
* **Attribution (meaning).** The lemma is the case S = {0, infinity}, height-0 coefficients of Guo-Sun-Wang, JNT 231 (2022), Thm 8 (= GNSW arXiv 2106.15881, Thm 17). The meaning referee read it in extracted text; I did not check it myself.

**Computed (referee).**
* Resultant checks with sympy for r = 2, 3 (nonzero, single v-weight, in the ideal).
* Gcd degrees for planted resonances decay as stated.
* The hypothesis is sharp: n = (l, l+1, l+2) gives gcd degree L.

**[J] Judge check of the step.** I went through it line by line. The points that carry the weight:
1. **Finite data for fixed Q.** For fixed Q, b' <= Q^(r-1), so v ranges over a finite set. The X-degree d = O(b'), so the support of R lies in a box 2d·B of bounded size. Hence the set E_Q of possible exponent differences of R is finite and **independent of the coefficients**. The hypothesis "n.u != 0 eventually" applied to finitely many u gives one threshold l_0(Q, B, (n_l)). Uniformity over coefficients is genuine.
2. **The v-weight bound.** R is v-homogeneous: R(Z S^v) = S^(const) R. So for exponents u, u' of R, n.(u-u') = (n - L v/b').(u-u'). Here ||n - Lv/b'||_inf <= L/(b'Q) and ||u-u'||_1 = O(b'), so the b' cancels and the bound is O(L/Q) with a constant depending only on (B, r). This is why Q -> infinity is allowed.
3. **R(t^n) controls the gcd with multiplicities.** R = alpha·A + gamma·C in the Laurent ring, and A, C are fixed combinations of the F_i. Every common zero of the F_i(t^n) of multiplicity m is a zero of R(t^n) of multiplicity >= m. R(t^n) != 0 because its monomials have distinct n-weights. Its zero count on C^* is (max weight) - (min weight) = O(L/Q).
4. **Coprimality survives.** v-substitution is an automorphism of K[Z^(+-1), X^(+-1)]. Stripping X^a removes the only unit. Gauss's lemma passes to C(Z)[X], so the resultant is nonzero. If d = 0, A is itself v-homogeneous and R = A works.

Verdict: **correct.** No cost is discarded. The gcd is the full Laurent gcd of the list, and the bound is on what remains after it is divided out.

**Withdrawn / Gap.** None.

### 1.2 Height, covering inequality, torus passage, one exceptional set (referee: uniform)

**Proved.**
* **(VII.1).** W = d_+ + d_- = sum_v (max_i a_iv - min_i a_iv). W = 2D in the conjugate-norm model, with D the **primitive** degree. h >= (1/2+delta)D is equivalent to h >= (1/4+delta/2)W. The exceptional set may be taken in the constant torus: real translations t -> t+s and t -> t^k preserve the class and the image, and the latter multiplies W and h by k.
* **(VII.2), sharpened (referee).** Take a finite-type characteristic-0 family of maps with W <= W_0, all equal to e mod t^(-h_0), whose evaluation map dominates T. Then
  4h_0 + 1 <= #(charged places of a generic member) <= W_0.
  * Proof: log-Wronskian of the source derivative and three parameter derivatives; simple poles only at charged places.
  * On the Pell family the determinant has ord_inf = 4h+1 exactly and W simple poles, so the inequality is sharp within one.
  * This generalises round 29 Theorem D to arbitrary multiplicities and the split model.
* **Passage to a fixed torus and (VII.3).** W = d_F + d_G - o(d_F + d_G), conditional on §1.1, which holds.
  * ||n|| is comparable to d_F + d_G.
  * o(·) cannot be improved to O(1): n = (mk', (2m+1)k').

**Corrected statement of the "one exceptional set" step (referee: uniform).**
> Fix n_* from the tail. Let P be the locally closed, Q-constructible set of coefficient lists (c, lambda) with **exactly** the frozen supports: every coefficient on the support of each F_i and G_i is nonzero, F_i - F_1 has exact support S_i1, and F_i G_i = F_1 G_1 as Laurent identities, with Z -> lambda Z, lambda in (C^*)^r. Every member, evaluated at Z = lambda t^(n_*), has W <= d_F^* + d_G^* and contact >= h_*. Moreover 4h_* > d_F^* + d_G^*.

So by (VII.2) no component of P dominates T, and the closure Y of its evaluation image is a proper Q-closed set. At t = 1, Y contains Phi_c(lambda) for all lambda, hence Phi_{c_l}(u^(n_l)) for every later term l, hence every later ratio curve.

If the supports are only bounded by the frozen ones, the top coefficient of F_1 can vanish and the contact drops. Astra's phrase "scaling a monomial by a nonzero scalar preserves them" shows that exact supports are meant.

**[J] Judge check of the passage from one n_* to a common exceptional set.**
1. **Contradiction scheme.** Negate the theorem with a Q-defined Z. Enumerate the nonzero P in Q[Y^(+-1)]. At step l pick a qualifying map with W >= l escaping Z(P_1) ∪ ... ∪ Z(P_l). This is a legitimate countable choice.
2. **Diagonal subsequences.** There are at most 8b distinct exponents, so finitely many combinatorial placements. For each of the countably many u in Z^r, "n_l.u = 0 eventually" or "!= 0 eventually" can be forced by a diagonal subsequence. The first set is a saturated subgroup Lambda, finitely generated, so its relations hold simultaneously on one tail. Fixed supports are then exponent-injective on a tail. Every "eventually" is over finitely many fixed objects or along this one diagonal subsequence: **legitimate**.
3. **Inputs to (VII.2).** In the quotient torus the lists F, G are formal presentations of q and 1/q, so F_iG_i = F_1G_1 holds formally. Each later term l is a member c_l of P. n_* is on the tail, so its contact h_* is the weight difference read from the supports, and the tail inequality 4h_* > d_F^* + d_G^* comes from h >= (1/4+eps)W and (VII.3).
4. **The contradiction.** Y is proper and Q-closed, so Y ⊂ Z(P_{l_0}) for some index l_0. Every curve with l > l_0 escapes Z(P_{l_0}), yet lies in Y. Contradiction.
5. **Uniformity.** Z_{b,eps} is a finite union of hypersurfaces Z(P_1..P_{l_0}), and W_0 is the corresponding index. Both depend on (b, eps) only. The threshold is not template-dependent: it is fixed before any tuple is chosen. Fixed-degree exceptional sets per exponent vector would **not** have given this, and Astra says so.

Verdict: **correct**, with the exact-support wording.

**Against the standing cautions.**
* Fixed-S Subspace: not used. S = {0, infinity} is intrinsic to the monomial torus.
* Gcd/denominator costs: kept. W is the primitive height and the resonant common factors are paid in (VII.3).
* Drifting phase: coefficients are kept exactly in the family.
* Template threshold: see point 5 above.

**Computed (referee).** Pell at fixed n composed with t^k is bounded-term (b <= 8n+2 monomials per numerator and denominator), dense at every height, with h/W = 4n/(16n+2) -> 1/4. So **the threshold 1/4 in VII.2 is optimal as b grows**: for every c < 1/4 the analogue with h >= cW fails for large b.

**Withdrawn / Gap.** None.

### 1.3 VII.3 / VII.4 circle obstruction and surviving relations (referee: circle)

**Proved.**
* **Triple lemma.** For any circle with real centre and any three distinct primitive marks: 2h <= d + 2u, where d is the reduced degree of the triple and u its off-circle degree per mark. It is sharp: with w = (3+4i)/5, the marks G1G2, conj(G1G2), G1 conj(G2), G_1 = t^d - w, G_2 = t^d + w, have r = 0 and 2h = D for d = 1, 2, 3.
* **(VII.4).** For five distinct primitive marks and any circle with real centre, 20h <= 9D + 18r, where r is half the number of roots of N off the circle, with multiplicity. Hence h >= (1/2+delta)D forces r/D >= 1/18 + 10delta/9. If r = 0: five marks have h <= 9D/20 and three marks have h <= D/2, attained.
* **Lines are genuinely excluded.** The analogue for a line is **false**, checked exactly by the referee:
  * three primitive marks with all roots on the imaginary axis, D = 5, h = 3;
  * five marks with D = 6, h = 3, so 20h = 60 > 54 = 9D.
* **VII.4.**
  * For 2h > D the relation space is conjugation-stable.
  * There are no 2- or 3-term relations.
  * Every four marks are independent.
  * A five-term relation is unique, real, of full support, and forces only 2|3 layers.
  * 2h > D is necessary: in the D = 4 cluster the five-term relation is not real.
* **Ramification test.** If there is a real relation a, then either the anchored ratios are multiplicatively dependent or >= 3h+1 finite places have sum a_i v(q_i) != 0. This is weak: for squarefree N it is implied by the pair bound h <= 3D/5.

**Computed (referee).**
* Slack checks of (VII.4) on the Pell families, all 56 five-subsets of the round-29 cube cluster, the D = 4 cluster and 25 random genuine tuples. The minimum slack off-circle is 32.
* r = 0 searches find no solutions at (D, h) = (4,2), (5,3), (6,3).

**Not claimed / no integer analogue.** The proof uses a circle inversion swapping infinity with a finite place. Z has no such automorphism.

**Withdrawn / Gap.** None.

### 1.4 VII.5 Pell families and pair-degree bounds (referee: pell; cross-checked by search and meaning)

**Proved (corrected statement).**
> A(t) = [[t^2+1, t],[t, 1]], beta = t + 2i, G_j = X_j + iY_j with (X_j, Y_j)^T = A^n (1, s_j)^T. Marks: z_0 = beta prod conj(G_l), and z_jk = conj(beta) G_j G_k prod_{l != j,k} conj(G_l). For r = 4 and edges 12, 13, 23, 14 this gives five marks; for r = 3 and edges 12, 13, 23, four marks.
> For **every** choice of pairwise distinct real seeds s_j the marks are primitive, distinct and of one norm, with D = 2rn+1 and h = 4n.
* Astra's condition "some edge has s_j + s_k != 0" only affects whether an anchor pair reaches contact 4n+1. The pair (z_12, z_13) is pinned to contact exactly 4n by G_s conj(G_u) - conj(G_s) G_u = 2i(s-u).
* Primitivity also holds for all distinct real seeds, not just generic ones.
* For fixed n, the union over seeds of the ratio curves is Zariski dense; already their values at t = 1 are. No single curve is dense.
* Five marks: D = 8n+1, h = 4n, so 4h + 2 = W. This is the maximum allowed by (VII.2) for a dominating family, with slack 1. Four marks: D = 6n+1, h = 4n, every delta < 1/6. Three marks: G_j = t^d - a_j + i gives D = 3d, h = 2d, dense at the endpoint delta = 1/6.
* Pair-degree bounds: h/D <= 2/3, 2/3, 3/5 for 3, 4, 5 marks. Equality h = 2D/3 is impossible at four marks (characteristic != 2).
* Within multiplicative Pell-block templates, an exhaustive search found no five-mark example with 2h > D. The best dense template has 2h - D = -1, and the best non-dense one 2h = D (3-dimensional family).

**[J] Computed, own exact Q(i)[t] arithmetic (`README_VIIL/pell5.py`, no sympy).** Seeds (1/2, -2, 7/3, 4) and (-1/3, 5/2, 6, -7). The script asserts equal norms z·conj(z), a full gcd of degree 0 (Euclid over Q(i)), distinctness and a common degree.

| seeds | n | D | h | all ten deg(z_i - z_j) | N squarefree | seed->phase Jacobian at t=1, smallest s.v. |
|---|---|---|---|---|---|---|
| (1/2,-2,7/3,4) | 1 | 9 | 4 | 5 | yes | 2.5e-2 |
| (1/2,-2,7/3,4) | 2 | 17 | 8 | 9 | yes | 3.7e-3 |
| (1/2,-2,7/3,4) | 3 | 25 | 12 | 13 | yes | 5.4e-4 |
| (-1/3,5/2,6,-7) | 1 | 9 | 4 | 5 | yes | 2.0e-2 |
| (-1/3,5/2,6,-7) | 2 | 17 | 8 | 9 | yes | 3.0e-3 |
| (-1/3,5/2,6,-7) | 3 | 25 | 12 | 13 | yes | 4.3e-4 |

* All marks are monic, and **all ten** pair contacts equal 4n. This confirms the referee's correction.
* The floating-point Jacobian has numerical rank 4, consistent with density.

**Integer specialisation (pell, meaning; specialisation lemma of VII.7).**
* At integers T the Gaussian gcd is bounded: it divides a fixed Bezout constant. Observed values: {160, 2080, 4160}.
* R ≍ T^D, and the five points lie on arcs ≍ R^(1/2 + 1/(2D)).
* For every n these give Zariski-dense primitive five-tuples with exponent alpha <= 4 + 1/(2n). **Hence alpha^ess_e(-K) <= 4 over Q, i.e. A_5(delta) is false for every delta < 0.** The profile is unbalanced (1|4 cuts), so nothing follows for (L_5°).

**Withdrawn / Gap.** None.

### 1.5 VII.6 characteristic 2 (referee: char2)

**Proved (corrected).**
* **Identities.** F_j(Y) F_j(Y+c) = P_V(Y)/(Y(Y+c)). deg_Y(F_i - F_j) <= 7 < 15. The full gcd is 1.
* **Dominance via Moore minors.** The values at Y = 0 are the 4x4 Moore minors, and y determines [x]. The referee checked this exhaustively over F_32: 322,560 bases, all fibres of size 1.
* **Genuine model.** c = 1, x_j = t^(n_j), Y = t^L + theta over F_4(t)/F_2(t).
  * Contact: h >= 8L - 8B (sharpened from 8L - 16B).
  * D <= 15L and h/D -> 8/15.
  * Each mark is a product of 15 factors with <= 6 monomials, so the term budget is b = 6^15, uniform in L, n and the gcd.
  * The tuples are Zariski dense over the algebraic closure of F_2(t), for every delta < 1/30.
* **Rank (corrected).** "Rank <= 1" holds only for the differential restricted to the Y = 0 slice. The rank that matters for (VII.2), that of the whole family (x, Y) with the source derivative included, is **2**. That is still < 4.
* **(VII.2) is false in characteristic 2.** The fixed split family has W_0 = 30 and h_0 = 8 and dominates, yet 4h_0 + 1 = 33 > 30. Generic smoothness is the only characteristic-dependent step of VII.2: the gcd estimate, the torus passage and the frozen family are characteristic-free.
* **Simplification (referee).** The split family F_j(t^k)/F_1(t^k) with generic constants already refutes the char-2 analogue of VII.2, with b = 5 and h/W = 4/15.

**Gap (on a non-claim).** The ternary analogue over F_9(t)/F_3(t) with a uniform term budget at delta < 1/30 is **not established**. For a fixed template the budget is uniform, but delta -> 0.01556 < 1/30 and the family is not dense. Astra asserts nothing here ("I make no term-budget assertion"), so no claim of hers is affected. I uphold the referee's GAP as a description of the ternary question.

**Not checked (referee's own).** The analogue with a 4-dimensional V and five hyperplanes would give h/D -> 4/7, i.e. failure for every delta < 1/14, with term budget 4. It is injective on F_32 and F_64, but dominance is unproved, since the Jacobian has rank 2.

### 1.6 Search for genuine char-0 five-mark tuples above D/2 (referee: search)

**Computed exactly.**
* H_5(7) = 3: there are no five marks at D = 7 with h = 4. Hence M(7,0) = 4. This comes from Groebner bases over Q in Singular with coincidences saturated; all types give the ideal (1), and a positive control passes.
* H_5(D) <= floor(D/2) for D <= 8 and D = 12. h = 7 is impossible at D = 11. Both by exact type combinatorics.

**Computed numerically (not interval-certified).**
* Real (conjugation-model) five-mark configurations exist at D = 11, h = 6 (4 of 18 cut types) and at D = 13, h = 7 (16 of 24 types). They are rigid (expected dimension 0 modulo affine maps).
* [J] I re-ran `ref7L_search/verify_d11.py` on the 40-digit roots of the D = 11 type-12 example:
  * the coefficients of t^10 .. t^6 of all ten z_i - z_j are <= 1.3e-39;
  * the t^5 coefficients are of order 0.1-1, so h = 6;
  * the norm discrepancy is 2.4e-79;
  * min |Im root| = 0.133.
  This reproduces the claim to the precision of the data. Distinctness and primitivity rest on the referee's checker.

**Gap.** D = 9, h = 5 in the conjugation model. The split model has 218 verified non-real solutions; the conjugation model shows no real solution in 3000 starts. The exact count did not finish. "H_5(9) = 4 in the conjugation model" is **evidence only**.

**Not checked.** Emptiness of the types with negative expected dimension (even D at h = D/2 + 1, K5 at D = 10, D = 14 at h = 8) rests on random-start Newton only.

**Significance.** Astra's strict inequality h > D/2 is not vacuous in characteristic 0. Any proof of the characteristic-0 A_5 must allow rigid curves with h = (D+1)/2 in the exceptional set. For fixed delta, such configurations exist only for D <= 1/(2delta), together with their t -> t^k compositions, which have the same image.

### 1.7 VII.7 and meaning (referee: meaning)

**Proved.**
* **Two-place Chebyshev version.** Ratios with numerators and denominators sums of <= b Chebyshev polynomials T_j(at + b_0) satisfy the conclusion of VII.2. Pulling back by at + b_0 = (u + 1/u)/2 doubles the height and gives contact h at u = 0 and u = infinity, so a dominating frozen family would need 8h <= 2W + o(W). This is a small extension. It has no integer meaning: integer values of T_j lie in varying quadratic orders, where VI.1's estimate is false.
* **Specialisation lemma.** Suppose ratio curves over Q(i)[t] with h/D >= 1/2 + delta_0 are Zariski dense. Then they refute A_5(delta) over Q for every 0 < delta < delta_0. Families over C alone do not specialise. The quantifiers are correct: one escaping curve per Z suffices.
* **Published function-field tools are vacuous here.** In the W = 2D normalisation, the published results pay chi_S, and chi_S ≍ W:
  * Mason-Stothers gives h <= #roots(N) - 1;
  * Brownawell-Masser, Corvaja-Zannier 2008, Capuano-Turchet, GNSW Thms 1-4 and the Guo-Sun-Wang / Levin-Wang gcd theorems are all vacuous;
  * Noguchi-Winkelmann-Yamanoi needs entire curves;
  * Wang's Subspace Theorem with S = {infinity} gives only 4h <= 5D, with an exceptional set that is a dense union of hyperplanes.
  The GNSW Theorem 4 statement was read in the PDF (text extraction). Astra's statement that the unrestricted question is open stands.

**Corrected "what transfers" paragraph.**
> The Dirichlet plus weight-homogeneous-resultant gcd step transfers verbatim to a rank-one unit group, i.e. one fixed integer base T with coefficients of height o(L log T): gcd(F(T^n), G(T^n)) divides R(T^n), which has size T^(O(L/Q)). It fails for rank >= 2, where the Subspace Theorem is needed (Bugeaud-Corvaja-Zannier, Levin 2019). The covering-family degree count has no arithmetic analogue.

**Corrected comparison with VI.1.**
* "Arbitrary complex coefficients" in VII.2 are height-0 constants. Their number-field analogue is fixed coefficients, so VII.2 is **not** Diophantinely stronger than VI.1.
* Its one real gain is that there is no hypothesis on raw common factors (no honest-block hypothesis).

## 2. What is new relative to the record, and correct

1. **VII.2, the sparse theorem in characteristic 0.**
   * It is the first exceptional set in the record that is uniform across all degrees, coefficients, exponent vectors and cancelled common factors for unbounded-degree families.
   * The new idea is the single-n_* frozen-family test.
   * Ingredients are standard: the gcd lemma is a Guo-Sun-Wang special case, the torus passage is Zannier/Fuchs-Zannier lacunary reduction, and the covering inequality is round 29 Theorem D generalised.
   * Its threshold 1/4 is optimal as b grows (referee).
2. **(VII.2) in place-count form.** 4h_0 + 1 <= #charged places, sharp within one on Pell, and false in characteristic 2 (W = 30, h = 8).
3. **(VII.4)**, the triple lemma, and the line counterexamples (referee). (VII.4) is a clean reflection bound but function-field-only.
4. **VII.4 structure at 2h > D.** Every four marks are independent, and a five-term relation is real, of full support, with only 2|3 layers.
5. **Dense polynomial Pell families.**
   * Five marks at W = 4h + 2, which is Theorem D's maximum.
   * Four marks for every delta < 1/6.
   * Three marks at the endpoint delta = 1/6.
   * The pair-degree bounds, and the impossibility of h = 2D/3 at four marks.
6. **alpha^ess_e(-K) <= 4 over Q** (pell and meaning referees; Astra did not draw it). The record had the window (15/4, 4) "open" at the upper end (expert note l.338) and dense exponent 4 as open (ref6L_subspace). The literature status was not checked. It is the expected answer from Vojta and McKinnon-type heuristics, so it is a calibration, not a surprise.
7. **Characteristic 2 with a uniform term budget** (b = 6^15, or b = 5 in the split model) via Moore minors. The record already had dense characteristic-2 and characteristic-3 refutations, but with unbounded complexity.
8. **From the referees.**
   * The first char-0 five-mark configurations with 2h > D (D = 11, 13; rigid, numerically verified).
   * H_5(7) = 3, hence M(7,0) = 4, exactly.
   * Complete type tables for D <= 15.

Not new: the covering argument itself (round 29), the gcd lemma (literature), the torus reduction (literature), and three-mark endpoint families (classical Jarnik configuration).

## 3. Rulings

**(a) VII.2 is a THEOREM.** Exact statement:
> Fix b >= 1 and eps > 0. There exist W_0 and a proper Zariski-closed Z_{b,eps} ⊂ T = G_m^4, defined over Q (both ineffective), with the following property. Let q = (q_2, ..., q_5) in C(t)^*, with q_j(infinity) = 1. Suppose each q_j is a quotient of two Laurent polynomials in t - t_0 (one common affine source change allowed), each with at most b monomials. Exponents, coefficients, roots and common factors are arbitrary, and the presentations are chosen separately for each j, before cancellation. If the intrinsic (primitive) height W(q) = d_+ + d_- >= W_0 and the contact h = min_j ord_inf(q_j - 1) >= (1/4 + eps) W, then q(P^1) ⊂ Z_{b,eps}.
>
> For genuine tuples, W = 2D (D the primitive degree), and this gives the characteristic-0 function-field A_5(delta) on the b-sparse class with eps = delta/2, for every delta > 0.

Proof hygiene:
* The frozen family must have exact supports (§1.2).
* The quantifiers are correct, and Z depends on (b, eps) only.
* The proof uses characteristic 0 exactly once, in (VII.2), and the statement is false in characteristic 2 (§1.5).
* The exponent 1/4 is optimal as b -> infinity.

What it does not cover: generic tuples. A generic polynomial of degree d has no b-term presentation when d+1 > 2b, and Astra says so. It has no integer content beyond VI.1.

**(b) (VII.4) is a THEOREM:**
> For five distinct primitive marks in the conjugate-norm model over C[t] and any Euclidean circle with real centre, 20h <= 9D + 18r, where r is half the number of roots of N off that circle, with multiplicity.

* It is false with "circle" replaced by "line" (§1.3).
* It bites only when >= 17/18 of the roots lie on one such circle.
* It has no integer analogue.

**(c) Pell families: correct, reproduced [J].**
* Five marks: D = 8n+1, h = 4n, for all distinct real seeds, dense for each fixed n.
* Four marks: D = 6n+1, h = 4n.
* Three marks: D = 3d, h = 2d.

What the integer specialisations give:
* Zariski-dense primitive five-tuples of Gaussian integers on arcs R^(1/2 + 1/(2(8n+1)) + o(1)) for every n. So **alpha^ess_e(-K) <= 4 and A_5(delta) is false for every delta < 0**.
* Four points: A_4(delta) is false for every delta < 1/6, consistent with VI.4.

They give nothing for A_5(delta > 0), for (L_5°), or for the exact exponent 4 on a dense set. Pell or multiplicative-block templates are capped at 2h <= D for five marks.

**(d) Characteristic 2: correct, with the rank correction.**
> Over F_4(t)/F_2(t), there are genuine five-mark sequences with uniform term budget b = 6^15, Zariski dense, with h/D -> 8/15, so the char-2 analogue of VII.2 fails for every delta < 1/30. The covering inequality (VII.2) fails in char 2 (W_0 = 30, h_0 = 8, 33 > 30), because the family's evaluation differential has rank 2 < 4.

The ternary analogue with a uniform budget is open (GAP, not claimed).

**(e) Status of the unrestricted characteristic-0 five-point question: OPEN.**
* Proved facts:
  * at most 24 marks (Part IV);
  * for each fixed degree D, Theorem D/(VII.2) puts all qualifying tuples with h > (2D-1)/4 into a proper closed Z_D, but the union over D may be dense;
  * true on the b-sparse class (VII.2);
  * true when most roots of N lie on one circle with real centre ((VII.4));
  * false for four marks below delta = 1/6.
* Data:
  * rigid real five-mark configurations exist strictly above D/2 (D = 11, h = 6; D = 13, h = 7, so delta = 1/22, 1/26);
  * none with h >= (D+3)/2 up to D = 14 (numerical);
  * none at D = 7, h = 4 (exact).
* [J] heuristic: the expected dimension 2D - 4h + 2 is negative once h >= (1/2+delta)D and D > 1/(2delta). That is what one expects if the char-0 A_5(delta) holds for every delta > 0, but excess intersection (as in Pell structure) is exactly what such counts miss. It is not evidence for a proof.
* No published function-field theorem decides the question (§1.7).

**(f) Consequences over Q. Blunt: nothing moves.**
* Question B (A_5(delta) for some delta < 1/30): **undecided, no progress.** Part VII contains no argument that applies to generic Gaussian integers. Its only integer output is the upper calibration alpha^ess <= 4, which goes in the wrong direction for Question B.
* Growth rate (2/3+eps) log R/loglog R: **unchanged.**
* Uniform bound M(C): **unchanged.**
* The characteristic-0 function-field question is itself an open function-field Vojta instance. It shares the integer problem's difficulty: generic complexity, no bounded presentation. It is not a stepping stone.

## 4. Expert-note updates (recommendations only)

1. §3, line 338: "The open window alpha^ess in (15/4, 4) is untouched." Replace with: "Unconditionally 10/3 <= alpha^ess_e(-K) <= 4. The upper bound comes from integer specialisations of the polynomial Pell five-mark families of round 31 Part VII.5, which are Zariski-dense with exponent 4 + 1/(2n) for every n (A_5(delta) false for every delta < 0). Question B asks for alpha^ess > 15/4; Vojta predicts exactly 4. Exponent exactly 4 on a Zariski-dense set remains open, and nothing is known for the balanced (L_5°)."
2. §3, function-field calibration bullet. Add: "In characteristic 0 (round 31 Part VII): the function-field A_5(delta) holds for every delta > 0 on quotients of b-sparse Laurent polynomials, with one ineffective Q-closed Z_{b,eps} across all degrees (exponent 1/4 optimal as b grows). It also holds when >= 17/18 of the roots of N lie on one circle with real centre (20h <= 9D + 18r). It is open in general. Rigid five-mark configurations with h = (D+1)/2 exist (D = 11, 13; numerical). The sparse theorem is false over F_4(t)/F_2(t) with a uniform term budget, because the covering inequality 4h+1 <= W fails in characteristic 2 (W = 30, h = 8)."
3. Round 29 Theorem D entry. Add the place-count form 4h_0 + 1 <= #charged places, its sharpness within one on the Pell family, and "characteristic 0 only".
4. Record M(7,0) = 4 as **exact** (replacing "evidence"). Keep "M(9,0) = 4 (evidence)" as evidence.
5. Literature table. Add Guo-Sun-Wang JNT 2022 Thm 8 (the sparse gcd lemma at S = {0, infinity}), and the remark that every chi_S-paying function-field tool is vacuous in the W = 2D normalisation.

## 5. Open sub-questions, ranked by honest probability of progress

1. **Bookkeeping** (P ≈ 0.95). Record the items of §4. *First step:* the five edits above.
2. **Characteristic 2 at delta < 1/14** (P ≈ 0.5; calibration value only). Prove dominance of the 4-dimensional-V, five-hyperplane map. *First step:* a Moore-kernel argument as in VII.6 for the 4-dimensional V, starting from the referee's injectivity data over F_64.
3. **Certify the rigid char-0 configurations and test them for a common subvariety** (P ≈ 0.4 for certification; ≈ 0.1 that it says anything structural). *First step:* run a Krawczyk certification on the D = 11 type-12 root. Then evaluate the ratio curves at D = 11 and D = 13 on the H_S hypersurfaces and subtori of bridge §4.2.
4. **Single-base integer version of VI.1 without the honest-block hypothesis** (P ≈ 0.3; low value). *First step:* write VI.1 with ratios that are quotients of <= b terms a·T^n (fixed T), using the rank-one resultant transfer of §1.7.
5. **Dense exponent exactly 4 over Q** (P ≈ 0.1). *First step:* note [J] that (VII.2) forbids W = 4h on any single dominating finite-type family over Q(i)[t] (it needs 4h + 1 <= W). So exponent exactly 4 densely can only come from a non-polynomial-family source, or from a countable union of non-dominating curves. First decide whether W = 4h + 1 (exponent 4 + 1/h) is attainable by a dominating five-mark family; Pell has W = 4h + 2.
6. **Unrestricted char-0 five-point question** (P ≈ 0.03 internally). *First step:* none internally. Pose it as outreach to function-field Vojta specialists, in the form of §3(e).
7. **Question B itself via this thread** (P ≈ 0). No mechanism in Part VII reaches generic Gaussian integers.

**Recommendation: close the Part VII / characteristic-0 A_5 thread** after item 1. Items 2 and 3 are optional calibrations and should be done only if cheap. If the thread continues at all, the one question worth asking is item 6, posed to outside specialists:
> Over C(t), with conjugation on coefficients: is there delta > 0 and one proper closed Z ⊂ T_5 containing every ratio curve of five distinct primitive marks of one norm with h >= (1/2+delta)D?

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| gcd estimate | 1 | correct, uniform in coefficients; a special case of Guo-Sun-Wang; internal tool only |
| uniform (VII.1-VII.2 sparse theorem) | 2 | correct, genuinely uniform theorem, new statement, optimal exponent; sparse class only, no integer content beyond VI.1 |
| circle ((VII.4), VII.4) | 1 | correct, clean, new; bites only for near-circular root sets; no integer analogue |
| Pell (VII.5) | 2 | correct, reproduced [J]; sharp for Theorem D; gives alpha^ess <= 4 over Q, a real but upper-side calibration |
| char2 (VII.6) | 1 | correct with rank correction; isolates characteristic dependence in (VII.2) |
| search (referee data) | 1.5 | new: char-0 configurations above D/2 exist (rigid); H_5(7) = 3 exact; decides nothing |
| meaning (VII.7) | 1 | correct; attribution and over-reading corrected; the alpha^ess <= 4 consequence drawn |
| **Part VII overall** | **1.5** | sound, honest, nothing false; one exact-support wording fix; no movement on Question B, the growth rate or M(C) |

Literature: GNSW Thms 1-4 and 17 were read by the meaning referee in extracted PDF text. The Guo-Sun-Wang attribution and the Zannier/Fuchs-Zannier provenance are the referees'; I did not fetch any source myself. My own checks are the gcd estimate and the one-n_* passage (§1.1, §1.2), the exact Pell reproduction (§1.4) and a re-run of the D = 11 verification (§1.6).
