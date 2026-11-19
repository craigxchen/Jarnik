# Judge report, round 31: Part VI of the Question B thread (Astra, "Genuine five-point tuples at the 10/3 endpoint")

Source: `astra_lemma_VI/answer.md` from line 1319 (`# Part VI`), with `vi_pell_four.py`.
Referee reports in this directory: `ref6L_subspace.md`, `ref6L_count.md`, `ref6L_kodaira.md`, `ref6L_pell.md`,
`ref6L_descent.md`, `ref6L_meaning.md`. Judge scripts: `README_VIL/pell4.py` (output `README_VIL/pell4.out`).
Claims of my own are labelled **[J]** and are unrefereed.

## 0. Verdict (five lines)

1. Nothing in Part VI is false. The one FALSE verdict (meaning referee) is aimed at a degree-uniform weakening of VI.1. Astra does not claim it and refutes it herself with the v_n example. I uphold the verdict as a statement about that weakening.
2. VI.1 is a theorem. If K, S and b are fixed and the coefficients have subpower height, then genuine five-tuples whose primitive pair blocks are sums of at most b S-units have essential exponent >= 4 (> 15/4). The exceptional set is ineffective. The proof respects the standing Subspace caution, because S is fixed by hypothesis and not tied to the parent. I checked the delicate step (formal Jacobian = 0 implies a fixed polynomial relation on the sequence) myself [J].
3. The rest of Part VI is also correct. VI.5 is a uniform, effective rigidity theorem for a two-real-parameter toy class. VI.3 shows the cubic cover is Mbar_{0,7}'s triple cover with kappa = 0, which closes the general-type route. VI.4 makes A_4(delta) false for every delta < 1/6; I reproduced this with exact arithmetic [J]. VI.6 proves: non-density in (*) is equivalent to A_5(delta) for some delta < 1/10.
4. Novelty: VI.1 is the first exponent->=4 theorem in the record for genuine five-tuples with free parent primes. The other new items are A_4 for delta < 1/6, the kappa = 0 identification and the VI.6 converse. The X^(3/10) count and the extension test are the record's endpoint chart and rigidity, restated.
5. Nothing changes for Question B, (*), the growth rate log R/loglog R, or M(C). (*) is still undecided, as Astra says, and by VI.6 it is as hard as beating Liouville/Cilleruelo-Cordoba at all. Close the (*) sub-thread. At most one narrow follow-up is worth trying (§5).

## 1. Per block: referee-corrected statements

### 1.1 VI.1, Subspace obstruction (referee: subspace; cross-checked by meaning)

**Proved.**

* **(VI.3) Fixed-place cancellation estimate.** Fix a number field K and a finite set S of its places. Take m sequences x_j in K^* with h(x_j) = O(L), sum_{v not in S} |log|x_j|_v| = o(L), and liminf h(x_j/x_k)/L > 0 for j != k. Then for every v in S, log|sum_j x_j|_v = max_j log|x_j|_v + o(L). In particular no fixed linear relation recurs on a subsequence.
  This is an Evertse (1984) / Corvaja-Zannier (2005, Prop. 1) type nondegenerate S-unit-sum estimate, extended to almost-S-units. It uses only the projective number-field Subspace Theorem with fixed K, S and fixed linear forms. The meaning referee read the CZ statement in the PDF. The subspace referee could not reach arXiv and checked against the standard form from memory.
* **Size reduction (corrected wording).** Read "remove zero subsums" as "take a representation of minimal length". Then no proper subsum vanishes, including subsums that contain -G_e, and all units have height O(log H).
* **VI.1 (exact statement).** Fix K ⊃ Q(i), a finite S containing the infinite places, b >= 1 and eps > 0. Then there exist a proper Zariski-closed Z ⊂ T_5, rho > 0 and H_0, all depending only on (K, S, b, eps) and all ineffective, such that the following holds. Take a primitive five-tuple of norm H >= H_0 whose anchored ratios satisfy q_i = eta_i prod_{e in E} (G_e/conj G_e)^{b_ie}, where:
  * there is at most one Gaussian block per cut;
  * prod_e |G_e|^2 <= H, where H is the **primitive** norm;
  * the eta_i have norm 1 and h(eta_i) <= rho log H;
  * each G_e is a sum of at most b terms a·u with u in O_{K,S}^* and h(a) <= rho log H.

  If Delta <= H^(-1/4-eps), the tuple lies in Z.
* **Corollary (meaning referee's Cor. U).** Suppose the five-tuples are built from at most m_0 template blocks J_s. Each J_s is a sum of at most b S-units with free exponents, and the blocks are pairwise coprime up to H^(o(1)). Suppose further the distances are <= R''^(1/2-eps). Then the tuples lie in an ineffective proper closed Z. This covers golden, silver and Pell (trace) blocks. Both hypotheses are needed: bounded m_0, and subpower overlap (gcd(f_45, f_135) = f_45 shows overlap can be total).
* **One-chord test.** Fix a finite set of primes. In every separating sextet of a putative (*) sequence, some real or imaginary coordinate of a block has prime-to-set part >= P^c eventually. The subspace referee verified the parity step for all 64 x 64 sign patterns.
* **Varying orders break (VI.3).** With v_n = (n + sqrt(n^2-4))/2, the relation v_{n+1} + v_{n+1}^{-1} - v_n - v_n^{-1} - 1 = 0 is exact, since the two traces are n+1 and n. It has no proper vanishing subsum: a subsum containing v_{n+1} must contain -v_n, and v_{n+1} - v_n = 1 + v_n^{-1} - v_{n+1}^{-1} > 1 then forces all five terms [J]. Yet h(v_n) ~ (1/2) log n, so (VI.3) fails with field degree <= 4.

**Corrected scope sentence.** Astra's scope paragraph says traces and Pell constructions are covered "with the parent primes unrestricted". It should read: covered **provided the gcd-reduced (primitive) blocks are still bounded fixed-field S-unit sums, or the raw Gaussian gcd is H^(o(1))**. Not covered:
* varying fields or orders (shown by v_n);
* an unbounded number of summands;
* S growing with H;
* generic Gaussian blocks. These form all but an X^(o(1)) proportion of Gaussian integers.

**Computed (referee).** The golden 8-point family has two-term S-unit blocks in Q(i, sqrt5). Its best 5-subsets have Delta·H^(1/4) = 3.70245917364 at n = 100, 160, 220, i.e. exponent exactly 4 on a curve. This is consistent with VI.1, which excludes only Delta <= H^(-1/4-eps) off Z.

**Withdrawn.** Nothing withdrawn from Astra. The meaning referee's FALSE is on the hypothetical statement "VI.3 or Subspace uniform over fields of bounded degree". Ruling: that statement is false. The v_n relation is a counterexample at degree 4. With b = 4 and degree <= 8, such a uniform statement would cover every Gaussian integer and prove far too much. Astra's §VI.1 and headline paragraph 4 already state that the degree cannot replace fixedness. No claim of Astra is affected.

**[J] Judge check of the delicate step of VI.1.** The step: the formal logarithmic Jacobian vanishes, therefore some fixed P in K[Y] vanishes on the numerical tail. Line by line:

1. **The field k.** It is the fraction field of K[small generators]/I, where I is the ideal of polynomials vanishing eventually on a diagonal subsequence. The polynomials over K in finitely many generators form a countable set, so a diagonal subsequence on which each is eventually zero or eventually nonzero exists. I is prime: if fg is eventually zero and f, g are both eventually nonzero, their product is eventually nonzero, a contradiction. k is finitely generated over K and has characteristic 0, so a separating transcendence basis exists. Its dual derivations extend to k(X) by D X_j = 0, together with X_j ∂/∂X_j.
2. **Spanning.** These derivations span Der_K(k(X)) over k(X). All 4x4 minors of the matrix (D log Q_i) vanish, so dQ_2 ∧ ... ∧ dQ_5 = 0 in Ω_{k(X)/K}. In characteristic 0 this is equivalent to algebraic dependence over K (Jacobian criterion). So there is a nonzero P in K[Y_2..Y_5] with P(Q) = 0 in k(X).
3. **From formal to numerical.** Clear denominators: N = P(Q)·prod(B~_e)^a·(denominators of the eta) = 0 in k[X^{±1}]. Each coefficient of N lies in k and comes from a polynomial in the small generators. N = 0 means each such coefficient is in I, i.e. is eventually zero numerically. There are finitely many coefficients. The units satisfy U = (Lambda-monomial)·X^nu exactly, by the basis change, so evaluation is a ring homomorphism into sequences modulo eventual zero. Hence N(numerical) = 0 eventually. The cleared factors are conj G_e != 0 and nonzero eta denominators, so P(q_n) = 0 eventually.
4. **Contradiction.** P is a nonzero Laurent polynomial over the normal closure, so it appears at some finite index of the enumeration. The sequence escapes it from that index on, and the subsequence inherits this.

The only point needing care is that (VI.4) be applied only to f with coefficients in k, whose nonzero coefficients evaluate eventually nonzero with height o(L). The construction guarantees this. The degree bookkeeping also checks:
* deg_v(D Q/Q) <= deg_v(Q-1) <= -(1/4+eps) e_v at the infinite places;
* deg_v of prod B B~ <= e_v, from the primitive-norm hypothesis;
* deg_v F <= -4 eps e_v at the infinite places and <= 0 at the finite places of S;
* summing over S against sum_v <w_v, nu> = 0 for any nu in the support gives a contradiction.

Verdict on the step: **correct**.

### 1.2 VI.2, counting (referee: count)

**Proved.**
* **(a)** The map from the ten block arguments to the four anchored arguments is -2 times the anchored K5 incidence matrix: rank 4, Smith form diag(2,2,2,2), surjective on tori. The Haar fraction for an arc of angle delta < pi is 5(delta/2pi)^4, with delta = 2 arcsin(Delta/2). Standard.
* **(c)** In X <= H <= 2X there are at most X^(3/10+o(1)) ordered primitive ratio classes with the ten-pair profile and width H^(-3/10+o(1)). The proof's reduced word C is the record's X + iL (cotangent.md D1, Theorem B). A one-line proof in that chart: L_0 = H^(o(1)), A ~ L H^(3/10), and X_j - A divides A^2 + L^2. In that form the profile hypothesis is not needed. The bound gives no closed set. It sits X^(1/2) above the heuristic, and the easy bound is X^(7/10).
* **Threshold lemma (referee, new).** A ten-pair profile needs H >= (5·13·17·29·37)^2 ≈ 1.4e12, because each split prime supplies at most two pair cuts.
* **Extension test.** |q_i - q_j| >= e^(-d_ij/2) and sum d_ij = 6W - 2W_sing. At the endpoint W_sing = o(W) and all w_ij = W/10 + o(W). This is the record's (cotangent.md B3), re-derived.

**Heuristic.** (b) predicts P^(20-4s) tuples at width P^(-s). That is P^(-4) per dyadic window at the endpoint, which is summable, so the heuristic predicts finitely many (*) tuples. This is the record's N^(-1/5) heuristic (cotangent.md l.260-264).

**Computed (referee, independent `enum_chart.py`).**
* The enumerator reproduces the record's run exactly: 407 classes, min K5 = 4.9475.
* No new cluster for A <= 10^6 with L = 6, 12 and K5 <= 8, i.e. H up to about 1e19.
* With K5 <= 20, L <= 30 and A <= 10^5: 1928 classes up to H = 4.6e15. None has all ten pair cuts; at most 8 distinct pair cuts occur. The best width exponent is 0.2299 (target 0.3).
* These searches are not exhaustive in L.

**Gap / Withdrawn.** None. Novelty is overstated: the "PROVED inverse bound" is the record's chart restated.

### 1.3 VI.3, Kodaira dimension and curve exclusion (referee: kodaira)

**Proved.**
* The cover w^3 = F = prod_{0<=a<b<=5} (y_a - y_b) (degree 15 on P^4) is the cofactor-cube-class cover of II.3/III.9.
* Its base model B is Kapranov's Mbar_{0,7}, with -K_B = 5h - E_3 - 2E_4 - 3E_5 = (1/3)·B_2.
* The branch divisor is B_2, all 21 two-point boundary divisors, and lies in |-3K_B|. K_Y = π^*(-K_B), and the singularities are canonical.
* h^0(B, jM) = 1 if 3 | j and 0 otherwise, so h^0(mK_Y) = 1 for all m >= 1 and kappa(Y) = 0. Equivalently, kappa(Mbar_{0,7}, -K) = 0. The referee confirmed this by an independent rank computation mod 1000003 (kernels 0,0,1,0,0 for j = 1..5).
* Y is not uniruled. No dominant algebraic family of rational curves has a constant (or algebraically varying) cube class.

**Corrected statement of the curve exclusion.**
* For each D_0, the fixed rational carriers of primitive-height degree <= D_0 that carry infinitely many (*)-points lie in a proper closed Z_{D_0}. Z_{D_0} **depends on D_0**. The union over D_0 may be Zariski dense: this is the template-dependent-threshold caution.
* Genus >= 1 carriers do not exist (Faltings; Mordell-Weil plus Roth in genus 1).
* Each rational carrier gives arcs R^(2/5) with O(1) cut errors. All collisions sit at one rational or real quadratic point, only pair layers occur, and 3 | div J.

Sporadic points and curves of growing degree are untouched.

**Withdrawn / Gap.** None. kappa = 0 is borderline, so no Lang/Vojta conclusion is available on this cover.

### 1.4 VI.4, four points (referee: pell)

**Proved.**
* **Pell family (1,1,1).** A = [[5,2],[2,1]], G_j = A^n v_j. The tuples are Zariski dense with Delta ≍ H^(-1/3), i.e. arc ≍ R^(1/3), and escape every proper closed Z with constants O_Z. Sharpness follows from triangle area >= 1/2 = abc/4R.
* **Mixed rates (6,7,7).** Two-sided Delta ≍ H^(-3/10) against primitive norm λ^(40n). The block-1 overlap is bounded **for all seeds** when the rates are 6 and 7; genericity is needed only for the equal-rate pair (corrected).
* **Fibonacci one-unit-class variant.** G_j = F_{k+1} + iF_k, k = 6(n + d_j), β = -1-2i.
* **Consequence.** A_4(delta) is false for every delta < 1/6. For delta >= 1/6 it holds vacuously by Cilleruelo-Cordoba. So k = 4 is settled. The record had only delta < 1/14 (expert note §3, gauss Thm 4.1); that result is for the balanced (L_4°) version, which remains a separate and stronger statement.
* **Extension statements.** An endpoint five-tuple needs all ten primitive triples three-point-sharp simultaneously. Four-point R^(1/3)-sharp sequences and the (6,7,7) family cannot be sub-tuples. This follows simply from the exponent: endpoint quartets have H_4^(-3/10).
* **Consistency with Cilleruelo-Granville.** The smallest limiting constant found is ((320+256√2)/7)^(1/3) = 4.60154. This exceeds the three-point constant (40 + (40/3)√10)^(1/3) = 4.34737 and the second threshold (40+20√5)^(1/3) = 4.39202. There is no conflict.

**[J] Computed (my own exact Gaussian arithmetic, `README_VIL/pell4.py`).** Seeds (3,1), (1,2), (2,-1). The script asserts:
* the four raw norms are equal;
* exact Euclidean Gaussian gcd division;
* all four quotient norms equal H;
* the four points are distinct.

| rates | n | N(gcd) | log10 H | log Δ / log H | chord/R^(1/3) | chord/R^(2/5) |
|---|---|---|---|---|---|---|
| (1,1,1) | 10 | 2 | 47.4 | -0.3101 | 12.62049 | 0.33 |
| (1,1,1) | 15 | 10 | 69.7 | -0.3209 | 7.38051 | 0.035 |
| (1,1,1) | 20 | 2 | 93.4 | -0.3215 | 12.62049 | 0.0098 |
| (1,1,1) | 30 | 10 | 138.6 | -0.3271 | 7.38051 | 0.00018 |
| (6,7,7) | 4 | 2 | 124.0 | -0.3032 | 5.4e3 | 0.39616 |
| (6,7,7) | 6 | 10 | 184.5 | -0.3033 | 3.5e5 | 0.24444 |
| (6,7,7) | 12 | 10 | 368.2 | -0.3017 | 4.6e11 | 0.24444 |

* The gcd is bounded (2 or 10, periodic in n). Here Δ = max chord / sqrt(H).
* For (1,1,1), chord/R^(1/3) takes the two exact constants 12.62049 and 7.38051, both above 4.34737. The local exponent tends to -1/3 like a constant over log H.
* For (6,7,7), chord/R^(2/5) is stable at 0.39616 / 0.24444, i.e. exponent -3/10.

This reproduces Astra's and the referee's claims with independent code.

**Novelty.** Four-point R^(1/3) families are Cilleruelo-Granville's (2007 §5). New are the escape quantifier and the resulting A_4 statement. The R^(2/5) four-point family adds nothing: CG's literal four-point question is answered by their own R^(1/3) families.

### 1.5 VI.5, moving affine pencil (referee: descent)

**Proved (corrected statement).**
* **Class.** A sequence in the ten-pair profile with blocks G_e = a_e B + c_e. Here B is an arbitrary Gaussian integer, and the moving coefficients a_e, c_e lie in Q(i) with h = o(log P), P = H^(1/20).
* **Conclusion.** For fixed 0 <= nu < 1, all large terms with Delta <= P^(-6+nu) (i.e. P^(-6) <= Delta <= P^(-5-(1-nu))) lie in a fixed closed set of dimension <= 2. That set is the closure of the d = 1 constant-cofactor template curves. It is absolute, possibly empty, and uncomputed.
* **Threshold and exponent.** The threshold is effective. The essential exponent on the class is >= 4.
* **Tight points.** Only nu < 2 is needed for (VI.6); nu < 1 is used in the line case, and that is the sharp point. "Degree at most four" in the Liouville step should read degree <= 2 (harmless).
* **Size of the class.** It is a codimension-18 slice, about P^(2+o(1)) of the P^20 configurations.
* **Provenance.** The uniformity step (rank <= 2) and the two-moment input are the record's III.6 and II.5. The new step is the polynomial-division descent that forces B onto a subpower-height rational line.

**Computed (referee).**
* Exact identities for t_ij = Im(κ A(B) conj D(B)).
* The II.5 sixteen-point system rechecked by Groebner (degree 16; at most 4 < 10 distinct values).
* Exhaustive box search: 976,152 strict k = 2 templates (nu = 3), none with k = 3, none with exponent < 4.

**Gap / Withdrawn.** None.

### 1.6 VI.6, meaning (referee: meaning)

**Proved.**
* **Forward direction.** A Zariski-dense sequence in (*) refutes A_5(delta) for every delta < 1/10. It does not refute (L_5°), since singleton cuts are o(W) and the profile is unbalanced, nor the uniform bound.
* **Converse.** If A_5(delta) fails for every delta in (0, 1/10), a diagonal sequence lies in (*). Two-sided width and the ten-pair profile are forced by the extension inequalities; the LP optimum max min d_ij = 3W/5 is unique. Hence **essential non-existence in (*) ⇔ A_5(delta) for some delta in (0, 1/10)**.
* **Relation to Question B.** Question B (delta < 1/30) implies this. "Strictly weaker" should read **"logically weaker; no separation known"**.

## 2. What is new relative to the record and correct

1. **VI.1.** The first theorem in the record giving essential exponent >= 4 for genuine five-tuples with unrestricted parent primes, on the class of bounded fixed-field S-unit blocks. It covers all profiles, so on that class it gives both A_5(delta > 0) and (L_5°). The record's only previous fixed-S results had S = the parent primes (round18 Thm 2.5, useless). The method is standard: Evertse / Corvaja-Zannier Subspace, a formal Laurent model and a logarithmic Jacobian. I could not find the theorem itself in the literature, and neither could the referees.
2. **A_4(delta) is false for every delta < 1/6** (VI.4 escape quantifier plus Cilleruelo-Cordoba vacuity). This is the full non-vacuous range; the record had delta < 1/14.
3. **The cubic cover attached to the cofactor cube class.** It is the triple cover of Mbar_{0,7} branched on B_2 in |-3K| and has kappa = 0. The general-type route on this cover is closed.
4. **Rigidity of fixed rational carriers.** One rational or real quadratic collision point, pair layers only, 3 | div J, and non-existence of genus >= 1 carriers. This extends III.6 in hypotheses but has a weaker conclusion.
5. **VI.6 converse:** (*) is equivalent to A_5(delta) for some delta < 1/10. This upgrades README_VL §3(b) to an equivalence.
6. **The VI.5 division descent** (new step, toy class).
7. Minor new items: the one-chord test, the v_n example, the referee's threshold lemma H >= 1.4e12, and the new negative search data.

Not new: the X^(3/10) count (record chart), the volume heuristic (record), the endpoint rigidity w_ij = W/10 (cotangent B3, README_VL 1.1), and the four-point families themselves (CG).

## 3. Rulings

**(a) VI.1 is a THEOREM.**
* Exact hypotheses: as in §1.1. K ⊃ Q(i), S finite containing the infinite places, and b fixed independently of the tuples. Coefficients a_ej and phases eta_i of height <= rho log H. The bound prod |G_e|^2 <= H is taken with the **primitive** norm H.
* Conclusion: an ineffective proper closed Z(K, S, b, eps) containing every tuple with Delta <= H^(-1/4-eps), H >= H_0.

Against the standing cautions:
* **Subspace exceptional sets.** Respected. S is the fixed set of places of the additive ingredients, not the parent primes. Z depends on (K, S, b, eps) only, so the claim "for all parents" is legitimate inside the class. Varying S is excluded and must not be read in.
* **Gcd and denominator costs.** Respected, with the scope correction: constructions with large raw gcd are covered only if the reduced blocks stay bounded S-unit sums.
* **Drifting phase.** Respected: eta_i is kept exactly in the coefficient field.
* **Template-dependent threshold.** Not applicable, since H_0 depends on (K, S, b, eps) only.

The ineffectivity is stated honestly. What VI.1 covers: fixed-order Pell, trace and sparse-digit constructions with honest blocks. What it does not cover: varying fields, an unbounded number of summands, growing S, or generic Gaussian integers. Generic Gaussian integers are what Question B is about.

**(b) VI.5 is a theorem**, with the corrected statement of §1.5. It is uniform, effective and uses no Subspace Theorem. It applies to a codimension-18 toy class inside the ten-pair profile only.

**(c) Kodaira and curve exclusion.**
* kappa(Y) = 0, h^0(mK) = 1: **proved** and independently verified.
* Non-uniruledness and the exclusion of dominant rational families: **proved**.
* Curve exclusion: **proved for each fixed degree bound D_0**, with Z_{D_0} depending on D_0. It is not a single closed set for all curves, and it excludes nothing sporadic.

**(d) Four-point families.**
* Correct, and reproduced exactly by me [J] (§1.4).
* The families are known (CG 2007). The Zariski-dense escape, and hence "A_4(delta) false for all delta < 1/6", is new to the record.
* Consistent with CG's three-point theorem: every limiting constant found exceeds (40 + (40/3)√10)^(1/3) = 4.347.
* Caveat: the referee could not reach the CG 2009 full text, so the claimed finiteness of families per t is unchecked.

**(e) Count and data.**
* The X^(3/10+o(1)) bound is correct but is the record's chart restated. It is X^(1/2) from the heuristic and produces no closed set.
* The referee's data are honest negatives: nothing with the ten-pair profile up to H = 4.6e15 (K5 <= 20, L <= 30), best exponent 0.2299 against 0.3.
* The data do not test (*): no ten-pair tuple can exist below 1.4e12, and the searches are not exhaustive in L.

**(f) Consequences.**
* **Question B: no change.** (*) is undecided. By VI.6, deciding (*) means beating Liouville/Cilleruelo-Cordoba at five points (some delta < 1/10). That is not easier than Question B, only logically weaker.
* **Growth rate (2/3+eps) log R/loglog R: no change.**
* **Uniform bound M(C): no change.**
* The only real gain is a no-go for constructions: no counterexample to A_5 or (L_5°) can come from bounded fixed-order unit templates with honest blocks. In addition, k = 4 is fully settled.

## 4. Expert-note updates (recommendations only; do not edit the record from here)

1. §3, the bullet "The 4-point version is false". Add: "In fact A_4(delta) is false for every delta < 1/6, the full non-vacuous range. Pell families G_j = A^{n} v_j with A = [[5,2],[2,1]] give Zariski-dense four-tuples at arc O_Z(R^(1/3)) (round 31 Part VI.4, referee-checked and reproduced). The 1/14 statement remains the best for the balanced (L_4°)."
2. §3, after the Vojta paragraph. Add: "**Unconditional on a thin class** (round 31 VI.1, proved). If each primitive block is a sum of at most b S-units of a fixed field K, with fixed S and subpower coefficients, then the essential exponent is >= 4 (ineffective Z(K, S, b, eps)). This covers fixed-order Pell, trace and golden templates with honest blocks. It fails for varying orders (v_n = (n+√(n^2-4))/2) and says nothing about generic tuples."
3. §3(b) / §5 item 5. Add: "A degree-uniform or growing-S version of VI.1 is false or meaningless (v_n example). Fixed-S results remain S-dependent."
4. §3(b). Replace any wording that treats the ten-pair profile (*) as a shortcut with: "Non-density in (*) is **equivalent** to A_5(delta) for some delta < 1/10 (round 31 VI.6). It is logically weaker than Question B, but no separation is known."
5. §5 (what any proof must avoid). Add: "General-type arguments on the cofactor cube cover: the cover is the triple cover of Mbar_{0,7} branched on B_2, with kappa = 0 (round 31 VI.3)."

## 5. Open sub-questions, ranked by honest probability of progress

1. **Record the closed items** (P ≈ 0.9; bookkeeping): A_4 for delta < 1/6, kappa = 0, VI.1 with its scope, and the VI.6 equivalence. *First step:* the expert-note edits of §4.
2. **(Q-res) Resonant templates.** Extend VI.1 to blocks that are unbounded products of a·λ^r + b·λ'^r in one fixed quadratic order (cyclotomic quotients), or find a counterexample there (P ≈ 0.1 for a decisive answer). *First step:* test the record's cyclotomic five-row family numerically for Zariski-dense arcs R^(1/2-delta) at five points.
3. **Merge the VI.3 carrier rigidity with III.6** (P ≈ 0.3, low value), to extend the dimension <= 2 bound to non-strict carriers and quadratic contact. *First step:* write the template statement with the "one collision point, 3 | div J" input.
4. **Certified solve of the strict four-moment system** of VI.5 (P ≈ 0.3, toy value). If it is empty, Z = ∅ on the class. *First step:* run msolve on the II.8 k <= 4 system over a real line.
5. **C port of the endpoint-chart search into H >= 1e20** (P ≈ 0.05 of finding anything; negative data only). *First step:* port `enum_chart.py`, with K5 <= 30.
6. **Decide (*) itself** (P ≈ 0.01). *First step:* none recommended. By VI.6 this is A_5 at some delta < 1/10, and no available tool reaches generic Gaussian blocks.

**Recommendation.** Close the integer-model / (*) sub-thread of Question B after item 1. Item 2 is the one optional narrow follow-up, and it should be done only as a numerical test first.

## 6. Scores (0-10, on the goal)

| block | score | reason |
|---|---|---|
| VI.1 subspace | 3 | correct, new to the record, legitimate fixed-S theorem; exponent 4 on a thin X^(o(1)) class; a no-go for constructions, not a step to a proof |
| VI.2 count | 0.5 | correct; record chart and heuristic restated; data negative and below the relevant height |
| VI.3 Kodaira | 1 | correct, new identification, kappa = 0 closes a route; curve exclusion only per degree bound |
| VI.4 four points | 1.5 | correct, reproduced [J]; settles k = 4 fully (A_4 false for delta < 1/6); no bearing on k = 5 beyond caution |
| VI.5 descent | 0.5 | correct, uniform and effective, but a codimension-18 toy class; half from the record |
| VI.6 meaning | 1.5 | correct equivalence; shows (*) is the hard core, not a shortcut |
| **Part VI overall** | **2** | sound throughout, honest and well calibrated; nothing false; no movement on Question B, (*), the growth rate or M(C) |

Literature: the Subspace / Evertse / Corvaja-Zannier statements were checked by the meaning referee against the CZ PDF and by the subspace referee from memory (arXiv was blocked in that session). Faltings, Roth and Mordell-Weil are used in their standard forms. I did not fetch any source myself. My own checks are the Jacobian-to-relation step (§1.1) and the exact Pell reproduction (§1.4).
