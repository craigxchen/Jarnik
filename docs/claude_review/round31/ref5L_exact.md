# Referee report, round 31, Part V (Question B thread), block "exact"

Scope: the positive statements of Part V (`astra_lemma_V/answer.md`, lines 1174-1317) and what they mean.
* (a) The rigidity theorem of §V.4: defects at most P^(2-eps), plus three isotropic endpoint edges that do not share
  a vertex, force exact norms.
* (b) The remark in §V.4 / summary item 11 that a dense nine-block family contradicts A_5. Also: what is known for
  5 to 9 isotropic blocks.
* (c) Which conjunction of hypotheses is the actual lever.
* (d) An assessment for Question B.

V.1 to V.2 (full isotropy, the size-free construction) belong to block "full" (`ref5L_full.md`). The V.3 families
belong to block "partial". I recompute the defects of those families only, as (a) requires.

Scripts are in `ref5L_exact/`. All of them are my own code.

| script | what it does |
|---|---|
| `defects.py` | Rebuilds the four-star and the five-block families of V.3 over Q(i)[T]. Rechecks all minors and the isotropy pattern, and computes the **degree of every defect**. |
| `qstar.py` | Tests whether the interpolation identity (IV.5) holds on the four-star, i.e. whether an interpolating Q_* exists. |
| `int_defect.py` | Measures log\|delta_i\|/log P on Astra's 191-digit integer instance (unsheared columns), and tests the linear identity exactly. |
| `graphs.py` | Enumerates all isotropic edge sets of K5 and checks the peeling conditions (§1.4). |
| `ffsearch.c` | An **exhaustive** function-field test of the theorem and of its extension over F_11[T] and F_13[T] (§1.5). |
| `peel_proof_check.py` | A trivial Gauss-lemma check. |

## 0. Verdict in seven lines

1. **(a) The theorem is CORRECT.** I checked every step and found no gap. It uses only item 1 of IV.4 (Q_* exists,
   and D·E is integral with D = prod|t_ij|), not the rank-one dichotomy. Width, index one, squarefreeness and
   primality are never used. It works for any rank of E, so it reaches past the rank-one range (below P^(1-eps)) up
   to P^(2-eps). An exhaustive F_13[T] toy search finds non-exact five-block nonstar models at defect exactly "P^2",
   so the exponent 2 is sharp there.
2. **The same peeling argument gives more.**
   * Below P^(1-eps), *one* isotropic incidence forces exact norms.
   * Below P^(3-eps), the whole IV.5 range, a "3-2-1" isotropic configuration forces exact norms. Every set of
     **8** isotropic blocks is one, and so is K4 on four marks (6 blocks).
   * The thresholds in P are 1, 2, 3 for 1, at least 5, and at least 8 isotropic blocks.
3. **The hypothesis fails massively on the V.3 families, as it must.** Both the four-star and the five-block curve
   have defects of exact degree 7 in T, i.e. delta ~ P^(7/2) = n_i^(7/8): relative error N^(-1/20). That exceeds
   P^3, so even the IV.5 linear identity fails on them. I verified this symbolically and on the 191-digit instance.
   Part V's isotropy is bought by giving up almost all norm accuracy: the IV.2 models had defect P^(1+o(1)), these
   have P^(3.5).
4. **(b) The nine-block remark is CORRECT, with this exact constant.** A Zariski-dense family with k isotropic blocks,
   |V_i|^2/n_i -> 1, subpower cofactors and width N^(-3/10+o(1)) yields genuine Gaussian 5-tuples at exponent
   alpha <= (20-k)/3. These violate A_5(delta) for every delta <= (k-8)/(2(20-k)).
   * k = 9 gives alpha <= 11/3 and delta <= 1/22 > 1/30.
   * k = 8 gives alpha <= 4, delta = 0, and contradicts nothing.
   * So 9 is the threshold, and it rests on a conjecture, not on a theorem.
5. **(c) Isotropy and size are complementary, and every lever proved so far is a conjunction of the two.** The
   proved trade-off is: defect P^(s-eps) together with peeling degrees (s, s-1, s-2). Each lever *characterises*
   genuine tuples. None gives an inequality for genuine tuples, since they satisfy every hypothesis with delta = 0.
6. **(d) Nothing moved on Question B.** The cleanest restatement is "fully isotropic integer 2x5 matrices with
   subpower multipliers c_i". It is literally the Gaussian factorisation V_i = G_i U_i with |U_i|^2 = N^(o(1)), so it
   is not a new or easier object. These reformulations capture only the ten-pair, Liouville-endpoint slice
   (alpha = 10/3 + o(1)), which is strictly weaker than Question B's 15/4.
7. **Score for Part V on the goal: 1/10.** It is correct, and the rigidity theorem is new and clean (it is the first
   statement in which *partial* isotropy pays). There is no progress on Question B.

## 1. (a) The rigidity theorem

### 1.1 Statement as verified

The statement as I verified it, with every hypothesis that is actually used:

> **Theorem (V.4 rigidity).**
> * *Setup.* Let J_e (e in E(K5)) be pairwise coprime positive integers with J_e = P^(1+o(1)), in the endpoint
>   orientation: n_i = prod_{e containing i} J_e and g_ij = J_ij.
> * *Columns.* Let V_1..V_5 in Z^2 be primitive, with J_ij | [V_i,V_j] = t_ij J_ij, 0 < |t_ij| <= P^(o(1)), and
>   delta_i = |V_i|^2 - n_i.
> * *Hypotheses.* Suppose max|delta_i| <= P^(2-eps) for a fixed eps > 0. Let S be the set of isotropic edges,
>   e = {i,j} with J_e | |V_i|^2 (equivalently J_e | |V_j|^2, see below). Suppose S contains three edges that do
>   not all pass through one mark.
> * *Conclusion.* For P large, delta_1 = ... = delta_5 = 0.
>
> In particular this holds whenever |S| >= 5.

### 1.2 Line-by-line check

1. **Propagation (V.1).** Mod p | J_e, the columns V_i and V_j are proportional by a unit, since J_e | [V_i,V_j] and
   both columns are primitive. So isotropy at one endpoint gives it at the other. "Isotropic edge" is therefore
   well defined. Since J_e | n_i, isotropy at (i,e) is the same as J_e | delta_i. **Correct.**
2. **E exists (IV.4 item 1, refereed in `ref4L_dichotomy.md`).**
   * For h = 4 and kappa = o(1), (IV.5) has right side at most P^(1+(2-eps)-4+o(1)) = P^(-1-eps+o(1)).
   * The left side is 0 or at least P^(-o(1)).
   * Hence the linear identity holds. Q_* exists with Q_*(V_i) = n_i, and D·Q_* is integral, D = H_t =
     prod|t_ab| = P^(o(1)).
   * E = I - Q_* satisfies E(V_i) = delta_i, and DE is in Z[X,Y]_2.
   * The range actually available is P^(3-eps). Astra says so. **Correct.**
3. **A vertex on two edges.**
   * Three edges in K5 cannot be pairwise disjoint (that would need six vertices), so some mark a lies on two of
     them.
   * If the triple is not a star at a, the third edge jk misses a.
   * J_ab·J_ac (coprime factors) divides delta_a. Their product is at least P^(2-o(1)), which exceeds
     |delta_a| <= P^(2-eps). So delta_a = 0. **Correct.**
4. **Factorisation.**
   * DE(V_a) = 0 with V_a primitive, so [V_a,X] divides DE in Q[X,Y].
   * Gauss's lemma (content of [V_a,X] is 1) makes L integral. I tested the content identity on 2000 random cases
     (`peel_proof_check.py`). **Correct.**
5. **L(V_j) = 0.**
   * D·delta_j = [V_a,V_j]·L(V_j) = t_aj J_aj L(V_j), and J_jk | delta_j.
   * J_jk is coprime to J_aj, so J_jk / gcd(J_jk, t_aj) divides L(V_j). That divisor is at least P^(1-o(1)).
   * But |L(V_j)| <= P^(o(1)) P^(2-eps) / P^(1-o(1)) = P^(1-eps+o(1)). So L(V_j) = 0.
   * If delta_j = 0 the same conclusion holds trivially, because [V_a,V_j] != 0. **Correct.**
6. **Conclusion.**
   * L(V_j) = 0 with V_j primitive gives L = m[V_j,X], m in Z.
   * At a third mark k: |D delta_k| = |m| |t_ak J_ak t_jk J_jk| >= |m| P^(2-o(1)).
   * Since |delta_k| <= P^(2-eps), m = 0, so DE = 0 and every delta vanishes. **Correct.**
7. **"Five edges contain a nonstar triple."**
   * Pairwise intersecting edge sets in K5 are stars or triangles. A set of at least 3 edges whose triples are all
     stars is therefore a star, which has at most 4 edges.
   * Hence 5 or more edges always contain a nonstar triple. **Correct.**
   * Exhaustive check (`graphs.py`, s = 2 row): 100 of 120 three-edge sets qualify (the 20 failures are the stars
     K_{1,3}), 205 of 210 four-edge sets (the 5 failures are the stars K_{1,4}), and all sets of 5 or more edges.

### 1.3 Roles of IV.4 and the rank-one branch

* **Only IV.4 item 1 is used** (the linear identity and the integrality of DE). The rank-one dichotomy (item 2) needs
  defect below P^(1-eps). The theorem lives mostly in P^1..P^2, where E can have rank two, and the argument never
  needs the rank.
* **Isotropy is precisely what kills the rank-one branch** (the analyst's Theorem R1, a genuine cluster for
  Z[sqrt k]). Suppose E = lambda·l^2 and two isotropic edges meet at a. Then delta_a = 0, so l = ±[V_a,·], and every
  other |delta_j| = |lambda| t_aj^2 J_aj^2 >= P^(2-o(1)).
  * So below P^(2-eps), rank-one is impossible once two isotropic edges meet.
  * Below P^(1-eps), one incidence suffices.
  * This answers the Part IV open branch in the presence of isotropy, and only there.
* **What survives at a star** (2 to 4 isotropic edges, all at a, defect below P^(2-eps)):
  * E = [V_a,X]·L(X)/D with L not proportional to [V_a,·]. That is a *split* (rational-isotropic) indefinite
    residual form, with delta_a = 0 and delta_j = t_aj J_aj L(V_j)/D.
  * This is neither constructed nor excluded. Neither Part V nor I have an example.

### 1.4 Extension (mine): the peeling theorem, valid in the whole IV.5 range

Same hypotheses, but with max|delta_i| <= P^(s-eps) for some 0 < s <= 3. Suppose there are distinct marks a, j, k
with:
* deg_S(a) >= s,
* deg_{S-a}(j) >= s-1 (isotropic edges at j avoiding a),
* deg_{S-{a,j}}(k) >= s-2.

Then all delta_i = 0.

*Proof.*
1. **E exists**, because s <= 3 (§1.2 step 2).
2. **delta_a = 0.** The product of deg_S(a) blocks is at least P^(s-o(1)) and divides delta_a.
3. **L(V_j) = 0.** DE = [V_a,X]L. The isotropic blocks at j avoiding a have product at least P^(s-1-o(1)) and are
   coprime to J_aj, so they divide L(V_j) up to a factor dividing t_aj. But |L(V_j)| <= P^(s-1-eps+o(1)).
4. **m = 0.** L = m[V_j,X]. The isotropic blocks at k avoiding a and j divide m t_ak t_jk, so up to P^(o(1)) they
   divide m. But |m| <= D|delta_k| / |t_ak t_jk J_ak J_jk| <= P^(s-2-eps+o(1)). So m = 0 and DE = 0. QED

What each level s gives:

| s | condition | consequence |
|---|---|---|
| 1 | one isotropic edge | below P^(1-eps), one incidence forces exact norms. It also follows from IV.4 plus the referee's "one zero forces all". |
| 2 | Astra's nonstar triple | Astra's theorem. |
| 3 | 3-2-1 configuration | `graphs.py`: 125/210 six-edge sets, 110/120 seven-edge sets, and **all sets of at least 8 edges** qualify (K4 on four marks does). The only failing 7-edge graphs have degrees (2,2,2,4,4), the complement of a triangle. |

Above P^3 the interpolation identity itself can fail, and nothing of this kind remains.

### 1.5 Numerical tests

**The V.3 families violate the hypothesis, as they must** (`defects.py`, my own reconstruction).
* **Four-star** at the certificate shape (w = (-21/16, 45/16, -153/32, -39/32)):
  * isotropic blocks exactly 12, 13, 14, 15, i.e. 8 incidences;
  * delta_1 = 0 identically (|G|^2 = n_1);
  * delta_2..delta_5 have **degree exactly 7** in T, while n_i has degree 8.
* **Five-block curve** (T-chart, all ten cofactors constant): isotropic blocks exactly 12, 13, 14, 15, 23; delta_1 = 0;
  delta_2..delta_5 of degree 7.
* With P ~ T^2 (J_e quadratic in T), the defects are P^(7/2) = n^(7/8). The relative error is T^(-1) = N^(-1/20).
  * "|V_i|^2/n_i -> 1" is therefore true, but only at rate N^(-1/20).
  * This is far coarser than IV.2's N^(-3/10), and above the whole IV.5 range P^(3-eps).
  * `qstar.py`: the quadratic Q through marks 1-3 misses n_4 and n_5 by polynomials of degree 6. So **no**
    interpolating Q_* exists, and the families lie outside the IV.4 framework altogether.
* **191-digit instance** (`int_defect.py`): log|delta_i|/log P = exact, 3.912, 3.909, 3.850, 3.888 (pre-asymptotic
  dilation constants; the asymptotic value is 3.5). The linear identity Σ n_i V_i / prod_j m_ij = 0 is exactly
  **false**.
* So the answer to (a)'s second question is **yes, it fails, by a factor P^(3/2)**. The five-block curve contains the
  nonstar triangle 12, 13, 23, so the theorem *requires* its defect to be at least P^(2-o(1)). It is P^(3.5).

**Exhaustive function-field test** (`ffsearch.c`). The proof is form-agnostic and works verbatim over k[T], with
"P^(s-eps)" read as "deg delta <= s-1".

*The model:*
* k = F_p, form I = XY (split);
* columns (z_i, z'_i) of monic quadratics, ten distinct linear blocks T - a_e;
* constant nonzero cofactors and primitive columns.

The search enumerates *all* such tuples with every deg delta_i <= 2 (the IV.5 range), up to translation and scaling of
T. It then checks propagation and the three predictions T1 (deg 0, one edge), T2 (deg <= 1, nonstar triple) and
T3 (deg <= 2, 3-2-1).

*Results:*
* **p = 11** (4 min 43 s): 28 800 tuples, all with deg delta = 2 exactly.
  * The isotropic sets: 5 280 with one edge, and 23 520 with two edges (15 840 sharing a vertex, 7 680 disjoint).
  * 0 propagation failures and 0 violations.
* **p = 13**: 293 760 tuples, all with deg delta = 2, including **480 with five isotropic blocks** (K_{1,4} plus one
  edge, a nonstar), and 0 violations. In this toy the P^2 endpoint of the theorem is therefore **sharp**. Details
  are in §6.

*Caveat:* in this toy no tuple reaches deg delta <= 1, so T2 and T3 are tested only vacuously there. The run does
confirm the following:
* the propagation lemma;
* that two isotropic edges, even sharing a vertex, do not force exactness at the P^(3-eps) level, consistent with
  the peeling threshold.

## 2. (b) The nine-block remark, and the table for 5 to 9 blocks

### 2.1 (V.2) checked

The bound H <= (1+o(1)) N prod_{e not in E} J_e holds. The argument:
* Factor the oriented Gaussian prime at each known isotropic incidence out of V_i; this is integral by V.1.
* The signed ranges of these factors cost prod_{e in E} J_e (range 1 at each such prime).
* Each residual range is at most Σ_i |s_ip(U)|, and prod_p p^(Σ_i |s_ip|) <= prod_i |U_i|^2 =
  (1+o(1)) prod_{e not in E} J_e^2.
* Ranges are subadditive. **Correct.** It needs |V_i|^2 <= (1+o(1)) n_i; a constant factor C costs only C^5.

### 2.2 The exponent

* With k isotropic blocks and equal weights, H <= N^((20-k)/10+o(1)).
* The genuine points have the *same* ratio point, hence the same Δ = N^(-3/10+o(1)), and Gaussian gcd one.
* So the approximation exponent satisfies alpha = log H / log(1/Δ) <= (20-k)/3.
* A_5(delta) forbids, off Z, Δ <= H^(-1/4-delta/2) (mutual distances <= R''^(1/2-delta), with H = R''^2). This is
  violated iff delta <= (k-8)/(2(20-k)).

| k | alpha bound | A_5(delta) violated for delta <= | consequence |
|---|---|---|---|
| 9 | 11/3 | **1/22** | Since 1/22 > 1/30, this contradicts A_5(delta) for every delta < 1/30, and Vojta's prediction alpha^ess >= 4. |
| 8 | 4 | 0 | Contradicts nothing (it meets the Vojta value). |
| 10 | 10/3 | 1/10 | Exact genuine tuples (c_i = 1). |

So the remark is **CORRECT**, and the sharp form is k >= 9. Two caveats:
* It contradicts *conjectures* (A_5, Question B, V_E(5)), not a theorem. Liouville alpha >= 10/3 is respected for
  every k.
* The genuine tuple has a ten-pair (plus one) profile, which is not balanced. A_5 applies and (L_5°) does not.
  (Block "full" makes the same point.)

### 2.3 The table requested

Endpoint orientation, width N^(-3/10+o(1)), subpower cofactors, escape from every fixed closed set unless stated.
P = N^(1/10).

| isotropic blocks k | (i) existence proved | (ii) forced genuine (exact norms) if max defect <= ... | what a dense family at 10/3 would give | (iii) open |
|---|---|---|---|---|
| 0 | dense, defect P^(1+o(1)) (IV.2) | none (the rank-one dichotomy only, below P^(1-eps)) | nothing | — |
| 1-3 | dense: subsets of the four-star, defect P^(3.5) | P^(1-eps) for k >= 1; P^(2-eps) for k = 3 only if the triple is a nonstar | genuine approximants with alpha <= (20-k)/3 >= 17/3 | dense families with defect below P^2 (stars), or below P^1 (none: forced) |
| 4 | **dense** (four-star K_{1,4}), defect P^(3.5) | P^(1-eps); P^(2-eps) unless the 4 edges form a star; a K_{1,4} is not forced below P^3 | alpha <= 16/3 | star K_{1,4} with defect below P^(2-eps): the split residual E = [V_a,X]L/D is not excluded |
| 5 | **on one rational curve** (star plus edge 23), defect P^(3.5); no dense family | **P^(2-eps)** (Astra) | alpha <= 5 | dense k = 5 at any defect; any k = 5 in the window [P^2, P^3) |
| 6, 7 | none | P^(2-eps); P^(3-eps) for 3-2-1 graphs (K4 included; all 7-edge graphs except the complement of a triangle) | alpha <= 14/3, 13/3 | existence of any kind |
| 8 | none | **P^(3-eps)** (peeling) | alpha <= 4 (the Vojta value; no contradiction) | existence of any kind |
| 9 | none | P^(3-eps); at any defect o(n) the transfer gives approximants violating A_5 | alpha <= 11/3: **contradicts A_5(delta <= 1/22), hence Question B as conjectured** | existence (expected impossible) |
| 10 | size-free only (V.2, multipliers unbounded) | any defect <= n_i (c_i = 1) | genuine at alpha = 10/3 (A_5 for delta < 1/10) | dense genuine tuples at alpha = 10/3 + o(1): open, the ten-pair slice of A_5 |

"Forced genuineness" excludes nothing unconditionally. It says that such models *are* genuine ten-pair tuples at
alpha = 10/3, and whether those exist densely is itself open (it is the delta < 1/10 slice of A_5).

## 3. (c) Which conjunction is the lever

Isotropy alone is not a lever: V.2 is correct per block "full", with arbitrarily small width and no size bound. Nor
are size alone (IV.2), divisibilities alone, or exact values of another form (R3). Every forcing statement in
Parts III-V is a conjunction of **size** (archimedean) and **local isotropy** (non-archimedean), with a
quantitative trade-off. All of them assume the ten private pair blocks of equal weight, primitive columns, the ten
divisibilities and cofactors of N^(o(1)).

| | size input | isotropy input | conclusion |
|---|---|---|---|
| G0 | delta_i = 0 | (automatic) | genuine. Part III; the dictionary in reverse. |
| G1 | \|V_i\|^2 < 2n_i (any relative error below 1) | all 20 incidences | c_i = 1, genuine. |
| G1' | \|V_i\|^2 <= N^(o(1)) n_i, i.e. Δ >= N^(-3/10-o(1)) | all 20 | genuine at H = N^(1+o(1)), with the same ratio point and the same exponent up to o(1) (V.1) |
| G2 | max\|delta\| <= P^(1-eps) | **one** incidence | exact (IV.4 plus §1.4, s = 1) |
| G3 | max\|delta\| <= P^(2-eps) | a nonstar triple (any 5 blocks) | exact (Astra) |
| G4 | max\|delta\| <= P^(3-eps) | a 3-2-1 set (any 8 blocks; K4) | exact (§1.4) |
| G5 | \|V_i\|^2/n_i -> 1 | 9 blocks | genuine *approximants* at alpha <= 11/3. This is a transfer, not an equivalence. |

On G1': under full isotropy |W_i| >= 1 gives Δ <= N^(-3/10+o(1)) automatically, so the "two-sided" hypothesis is
exactly the archimedean size bound in disguise, not a third ingredient.

The pattern: each unit of log P removed from the size requirement is paid for by one more layer of isotropic
incidences. This holds in the peeling order (s, s-1, s-2), up to the IV.5 ceiling P^3. Beyond P^3 nothing is known.

**None of these yields an inequality for genuine tuples.** Every one is a characterisation (G0-G4) or a one-way
transfer *into* the genuine set (G5). Genuine tuples satisfy all the hypotheses with delta = 0 and full isotropy,
so the theorems say only that "fake" models must sacrifice size or isotropy. They do not bound genuine widths. As
the brief predicted, these are characterisations.

## 4. (d) Honest assessment

* **Has anything moved?** No. The growth rate, the uniform bound and Question B are unchanged. The new rigidity
  theorem is a correct, clean, elementary statement about integer models. It locates where fakes must fail:
  * any five isotropic blocks force defect at least P^(2-o(1));
  * any eight isotropic blocks force defect at least P^(3-o(1));
  * one incidence forces defect at least P^(1-o(1)).

  It produces no exponent for genuine tuples.
* **Is there a cleaner formulation?** Formally yes. The ten-pair, Liouville-endpoint slice of Question B is
  equivalent to: *no Zariski-dense family of primitive integer 2x5 matrices with ten private pair blocks, subpower
  cofactors and isotropy at a nonstar triple of blocks, with defect below P^(2-eps).* But this is a dressing of
  "genuine tuples". The full-isotropy plus subpower-multiplier version is literally V_i = G_i U_i with |U_i|^2 =
  N^(o(1)), the Gaussian factorisation. Nothing is easier in it.
  * The isotropy congruences are exactly the information that the Gaussian prime pi_e divides V_i.
  * The size constraint is the same archimedean constraint as before.
  * The only "new" feature is a flexibility in the trade-off (partial isotropy against defect). That would help
    only if one could prove that genuine-like objects in some auxiliary construction automatically have small
    defect and some isotropy. No such auxiliary construction is in the record.
* **Scope caution.** All of Parts III-V address the ten-pair cut profile at exponent 10/3 + o(1). Question B
  (delta < 1/30) also needs exponents up to 15/4 and the fifteen-cut profiles. Excluding dense genuine ten-pair tuples
  at 10/3 is strictly weaker than Question B.
* **Novelty (grep of `docs/claude_review` and `round31`).**
  * New: the forcing of genuineness by *partial* isotropy plus defect (V.4), and the nine-block transfer (V.2).
  * The nearest earlier item is Part III's "bounded mean log slack in a fully isotropic integral model forces
    genuine points on a positive fraction of marks" (README_III §1.4).
  * The factorisation V_i = G_i U_i and the range formula are older (per block "full").
* **Score: 1/10** on the goal. The positive mathematics is correct. The Part V families are less accurate in norm
  than Part IV's (n^(7/8) against n^(1/4+)), which shows how hard isotropy and size are to combine. That is evidence,
  not proof, that the conjunction is the right lever.

## 5. Recommendations

1. In the record, state the theorem with its actual hypotheses (§1.1): no width, index or primality is used; only
   IV.4 item 1. Add the peeling extension (§1.4), with thresholds P^1, P^2, P^3 for 1, at least 5 and at least 8
   blocks.
2. Record that the V.3 families have defect n^(7/8) (relative N^(-1/20)) and fail the IV.5 identity, so "|V_i|^2/n_i
   -> 1" must always carry its rate.
3. Record the nine-block constant exactly: alpha <= (20-k)/3, contradicting A_5(delta <= (k-8)/(2(20-k))).
4. The one open, finite-looking question this pass isolates: **can an isotropic star (2 to 4 edges at one mark)
   coexist with defect below P^(2-eps)**, via a split residual E = [V_a,X]L/D? A proof that it cannot would push the
   forcing threshold to 2 blocks at P^2. Either answer is calibration only.

## 6. Exhaustive F_13 run (completed)

`./ffsearch 13 ex13.txt` (about 25 min) found **293 760** tuples, all with deg delta = 2 exactly. There were 0
propagation failures and 0 violations of T1, T2 or T3.

| isotropic sets | count |
|---|---|
| no edge | 32 400 |
| one edge | 115 200 |
| two edges | 114 960 (stars or disjoint) |
| three edges | 29 280, **all nonstar** |
| four edges | 1 440, nonstar |
| five edges | 480, nonstar; all of shape "K_{1,4} plus one edge", exactly the graph of Astra's five-block curve |

`check_ff_example.py` rechecked one five-edge tuple independently in sympy:
* ten distinct linear blocks, all minors of degree 1, primitive columns;
* isotropic blocks 12, 13, 14, 15, 23;
* delta_1 = 0, and deg delta_2..delta_5 = 2.

**Reading.** In this toy, defect "P^2" (degree 2) is compatible with five isotropic blocks and a nonstar triangle,
while degree at most 1 (P^(2-eps)) never occurs. So the exponent 2 in Astra's theorem is **sharp in the
function-field toy**: the theorem cannot be improved to P^(2+eps) with 5 blocks. This is a toy statement (split form,
monic (2,2) columns, constant cofactors). It is not a statement over Z, but it is the honest numerical test
available. The 3-2-1 configurations (at least 6 edges) never occur, so T3 is untested beyond vacuity.
