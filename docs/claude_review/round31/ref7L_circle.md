# Referee report ref7L_circle: Part VII, sections VII.3 (two-real-point / circle bound) and VII.4 (linear dependences)

Scripts: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_circle/`
(`circ.py` checker, `fam_pell.py`, `fam_record.py`, `rand_tuples.py`, `vii4_checks.py`, `circle_search.py`,
`line_search.py`, `line_verify.py`; outputs `c52.txt`, `c53.txt`, `c63a.txt`, `c63b.txt`, `line5.txt`).
All polynomial identities, gcds, degrees and contacts use exact arithmetic over Q(i) with sympy. Root positions come from
mpmath at 60 digits on the squarefree parts of exact `sqf_list` factorizations, so multiplicities are exact. "On the
circle" uses a tolerance of 1e-25.

## Verdicts in brief

| claim | verdict |
|---|---|
| (VII.4) 20h <= 9D + 18r for five distinct primitive marks, any circle with real centre | **CORRECT** (checked line by line; numerically consistent; triple lemma sharp) |
| triple lemma 2h <= d + 2u | **CORRECT**, sharp at r = 0 (exact equality family below) |
| summation: sum d_T <= 9D, sum u_T <= 9r | **CORRECT** |
| "all roots on one circle => even three marks need 2h <= D" | **CORRECT**, sharp |
| "a line is not included" | **CORRECT, and the line analogue is FALSE** (exact counterexamples below) |
| VII.4: conjugation-stability of the relation space for 2h > D | **CORRECT**; the hypothesis is necessary (fails at 2h = D, exact example) |
| VII.4: no three-term relations (h > 0); every four marks independent (2h > D) | **CORRECT** |
| VII.4: a five-term relation is real, has full support, and forces only 2\|3 layers | **CORRECT** |
| VII.4: ramification test (multiplicative dependence, or >= 3h+1 places with sum a_i v(q_i) != 0) | **CORRECT**; its force is limited (see below) |
| integer analogue | none known; the natural analogue is the record's open "twin"/residue-class problem |

## (a) VII.3: the bound (VII.4)

### Triple inequality, line by line

Setting: F_1, F_2, F_3 is the triple with its common gcd removed. It has degree d, the leading coefficients are equal
(q_j(inf) = 1), and the triple contact is ord_inf(F_i/F_j - 1) >= h, because ratios are unchanged. Use a real affine
coordinate x = (t-c)/rho. This commutes with coefficient conjugation and sends the circle to |x| = 1.

1. **Reversal.** Define F^*(x) = x^d Fbar(1/x). Then F_i^* - F_j^* = x^d (Fbar_i - Fbar_j)(1/x) has order
   >= d - deg(F_i - F_j) >= h at 0. Also F_i^*(0) = conj(lc) = 1. **OK.**
2. **Factorization.** Conjugation preserves a circle with real centre, so roots of Fbar_i are conjugates of roots of
   F_i, on-circle to on-circle. Hence U_i Ubar_i = N_off, which gives deg U_i = u for all i. For |alpha| = 1,
   (x - alpha)^* = 1 - alphabar x = -alphabar (x - alpha). So K_i/K_i^* is a product of the factors -alpha, a constant
   of modulus 1. **OK.** (alpha = 0 is never on the circle.)
3. **Degree of p.** N_off^* = U_i^* (Ubar_i)^*, since reversal is multiplicative. So
   H_i N_off^* = c_i U_i (Ubar_i)^*, a polynomial of degree <= 2u. The point p = [H_i] therefore has degree <= 2u,
   possibly less after removing a gcd. **OK.**
4. **Order at infinity.** Subtracting column 1 gives Delta = (P_2-P_1)(F_3-F_1) - (P_3-P_1)(F_2-F_1), of degree
   <= 2u + d - h. Equivalently, Delta has order >= h at infinity as a section of degree d + 2u. **OK.** Unreduced
   polynomial coordinates for p do no harm: the bound only gets weaker.
5. **Order at zero.** This is the step I was asked to check. H_i = c_i U_i/U_i^* is regular at 0, because
   U_i^*(0) = 1 for monic U_i. Write F_i = H_i F_i^* = H_i (F_1^* + O(x^h)). With M = N_off^* (M(0) = 1) and
   lambda = F_1^*/M (regular at 0), this gives F_i - lambda P_i = O(x^h). The row operation f -> f - lambda p leaves the
   determinant unchanged, and every entry of the new row is O(x^h). So ord_0 Delta >= h. **The claimed vanishing order
   at zero is right.** "Agree projectively to order h" is exactly this.
6. **Conclusion.** Delta is a polynomial of degree <= d + 2u - h divisible by x^h. If 2h > d + 2u, then Delta = 0.
   **OK.**
7. **Case: H_i not all equal.** Then e and p are independent over C(x), so F_i = A + B P_i. If two of the H_i were equal,
   the corresponding F_i would be equal, so the H_i are distinct and B != 0. At generic real x:
   * |F_i(x)|^2 = N'(x) is a common value.
   * |H_i(x)| = |F_i(x)| / (|x|^d |F_i(1/x)|) is also common.
   * So the three distinct points F_i(x) lie on the circle about 0 of positive radius, and on the circle about A(x) of
     radius |B(x)| |H(x)| (> 0, since the points are distinct).
   * Three distinct points determine their circle, so A(x) = 0 at infinitely many real x, hence A = 0.
   * Then F_i^* = F_i/H_i = B for all i, so the F_i are equal. Contradiction. **OK.**
8. **Case: all H_i equal (= H).** Then F_i - F_j = H (F_i^* - F_j^*) has order >= h at 0, since H is regular there, and
   degree <= d - h. So 2h > d forces F_i = F_j. **OK** (this case gives the stronger 2h <= d).

Hence 2h <= d_T + 2u_T for every triple and every circle with real centre. Independent check: the per-triple minimum
over all candidate circles of d_T + 2u_T - 2h_T is >= 0 in every tested triple. The tested triples are all triples of
every family listed below (several hundred).

**Sharpness of the triple lemma (new exact family).** Take w = (3+4i)/5 (|w| = 1), G_1 = t^d - w, G_2 = t^d + w and
the marks G_1G_2, conj(G_1G_2), G_1 conj(G_2). These are primitive, with all roots on the unit circle (r = 0), D = 2d
and h = d. So 2h = D exactly, for d = 1, 2, 3 (`vii4_checks.py`, triple slack 0). Astra's three-mark statement
"2h <= D when all roots lie on such a circle" is attained.

### Summation

At a non-real root alpha (real roots are impossible: at real alpha, v_alpha(z_j) = v_alpha(N)/2 for all j, which is a
common factor), let e_j = v_alpha(z_j), so v_alphabar(z_j) = m - e_j. Primitivity gives min e = 0 and max e = m. The
per-mark contribution to D is m, which equals the number of thresholds. After removing the triple's gcd at alpha and
alphabar, the triple's contribution is max_T e - min_T e, the number of thresholds that separate T.

A 2|3 cut fails to separate only the 3-side triple, so it separates 9 triples. A 1|4 cut fails to separate the four
triples inside the 4-side, so it separates 6. Off-circle pairs are counted the same way. So sum_T d_T <= 9D and
sum_T u_T <= 9r, and summing 2h <= d_T + 2u_T over the 10 triples gives 20h <= 9D + 18r. **CORRECT.**

The corollary r/D >= 1/18 + 10 delta/9 for h >= (1/2 + delta)D is arithmetic. The computed sum_T d_T agrees with
9 D_{2|3} + 6 D_{1|4}: for example 72 = 8·9 for Pell5 at n = 1, where the beta layer is a 1|4 cut, and 36 = 9·4 for the
D = 4 cluster, which has only 2|3 cuts.

### Standing cautions

* The height is the primitive one: gcds are removed per triple, and per-triple reduced degrees are used in the sum.
  Using the five-tuple's D per triple would give the false count 10D.
* No phase is dropped: the c_i are carried and cancel in the circle argument.
* No uniformity issue arises: the statement is per tuple with explicit constants.

### Numerical tests of (VII.4)

These come from `fam_pell.py`, `fam_record.py` and `rand_tuples.py`. Slack = 9D + 18 r_min - 20h, where r_min is
minimized over all circles with real centre.

| family | D | h | r_min | slack | min triple slack |
|---|---|---|---|---|---|
| Pell 5 marks n=1,2,3 | 9,17,25 | 4,8,12 | 7,15,23 | 127,263,399 | 6,16,26 |
| Pell 4 marks n=1,2,3 (triples only) | 7,13,19 | 4,8,12 | 5,11,17 | – | 6,16,26 |
| three-mark family, a=(1,3,8), d=1,3,6 | 3,9,18 | 2,6,12 | 1,6,12 | – | 1,9,18 |
| three-mark family, a=(1,-1,3) | 3,9,18 | 2,6,12 | 1,3,6 | – | 1,3,6 |
| cube cluster (round 29 Thm C), all 56 five-subsets | 8 | 4 | 4 | min 64 | 4 |
| D=4 six-point cluster (roots a_k +- ic), all 6 five-subsets, 3 parameter choices | 4 | 2 | 2 | 32 | 1 |
| 25 random genuine five-tuples with h = 2, 3 imposed exactly | 4..8 | 2,3 | 2..6 | min 32 | min 1 |

The Pell families are far from equality, since their roots are in general position (r ≈ D). The record clusters are also
far from equality.

**Searches with all roots on one circle (r = 0).** Here (VII.4) predicts 20h <= 9D. In `circle_search.py` the marks
are e^{i sigma_jk theta_k}, and contact >= h is the system sum_k (sigma_jk - sigma_1k) sin(m theta_k) = 0 for m < h.
All 15^D-type column codes were tried, with random-start least squares, then nondegeneracy, distinctness, primitivity
and the true contact were checked.

* D=4, h>=2 (forbidden: 40 > 36): 0 solutions among 1435 codes.
* D=5, h>=3 (forbidden): 0 solutions among 7693 codes.
* D=5, h>=2 (allowed, sanity check): 43 solutions, max contact 2. Slack is 45 - 40 = 5, the closest to equality found
  for five marks.
* D=6, h>=3 (forbidden: 60 > 54): 0 solutions among 30450 rank>=3 codes (2 random starts each)

The search is heuristic: random starts can miss solutions. It is consistent with the proof, and the sanity run shows it
does find on-circle five-tuples when they are allowed.

**What (VII.4) says for special parents.**
* Primitive parents never have real roots, so "all roots real" is impossible.
* If all roots of N are on one circle with real centre, for example N(t) = t^{2D} - 2cos(phi) t^D + 1 or other
  unit-circle/cyclotomic-type parents, then h <= 9D/20 for five marks and h <= D/2 for any three marks. So no such
  parent carries even a genuine triple with 2h > D, and the triple bound is attained.
* Astra's VII.7 Chebyshev setting is not of this type: T_j(at+b) has real roots, and the parents there are not
  restricted.

**Line analogue: FALSE (exact).** The conjugation-invariant lines are the real axis and the vertical lines Re t = c.
* Real axis: primitive parents have no real roots, so r = D and the statement is empty.
* Vertical lines, the limits of circles with real centre: the analogue fails.
  * `line_search.py` / `line_verify.py`: the marks z_j = prod (t - i s) with s-multisets (-3,-1,-1,1,2),
    (-2,-1,-1,-1,3), (-3,-2,1,1,1) are three distinct primitive marks of D = 5 with all roots on the imaginary axis and
    pairwise difference degree 2, so h = 3 and 2h = 6 > D = 5. The polynomials (exact):
    t^5+2i t^4+6t^3+8i t^2+5t+6i, t^5+2i t^4+6t^3+20i t^2-19t-6i, t^5+2i t^4+6t^3+4i t^2+13t-6i.
  * Five marks: s-multisets (-5,-1,-1,2,3,4), (-4,-3,1,1,2,5), (-3,-2,-1,-1,4,5), (-4,-2,-1,1,3,5), (-5,-2,1,1,3,4)
    give D = 6, h = 3, so 20h = 60 > 54 = 9D with every root on one line.
* Mechanism: for the imaginary axis the reflection is t -> -tbar. Composed with conjugation this is t -> -t, which
  fixes infinity, so the symmetry gives no second contact point. For a circle, x -> 1/x swaps infinity and the centre,
  and that is the whole proof.
* So "a line is not included" is not a gap in the statement but a necessary exclusion. The constant in (VII.4) is
  therefore not continuous under the degeneration circle -> line: r jumps.

## (b) VII.4: linear dependences

* **Conjugation-stability (2h > D).**
  * Sum a_i q_i = 0 gives sum a_i = 0 at infinity.
  * Conjugating and using zbar = N/z gives sum abar_i/q_i = 0.
  * 1/q = 2 - q + O(t^{-2h}), so sum abar_i z_i = O(t^{D-2h}), which is 0. **CORRECT.**
  * The hypothesis is necessary. In the D=4 six-point cluster (2h = D) the five-mark relation is unique, and applying
    its conjugate gives the nonzero constant -112i (`vii4_checks.py`). So at the boundary the relation is not real.
* **Two- and three-term relations.**
  * Two-term: equal leading coefficients force equal marks.
  * Three-term: q_3(x) lies on the unit circle and on a circle centred at -a_1/a_3 != 0 for every real x, so it takes
    <= 2 values and is constant, hence 1.
  * **CORRECT**, needs only h > 0.
* **Four marks are independent (2h > D).**
  * A four-mark relation has full support, is unique, and is real up to a scalar.
  * The curve {sum a X = 0, sum a/X = 0} has quadratic term sum a_i u_i^2 on the plane. It is nondegenerate on
    {u_1 = 0, sum a u = 0}: the ambient diagonal form is nondegenerate, the hyperplane is the B-orthogonal complement of
    the vector 1, B(1,1) = sum a = 0, so the radical is span(1), and u_1 = 0 is a complement to it.
  * So there is an ordinary node, the germ of the parametrized curve lies on one smooth branch, and the tangent form
    ell satisfies ell(u(t)) = O(|u|^2) = O(t^{-2h}).
  * z_1 · ell(u) = sum b_i (z_i - z_1) = O(t^{D-2h}), which is 0.
  * ell is not proportional to the plane equation on the 2-plane, so this is a second relation, which yields a
    <= 3-term relation: contradiction. **CORRECT.**
  * Checked: the Pell four-mark families (2h > D) have rank 4 for n = 1, 2, 3, and the Pell five-mark families have
    rank 5.
* **Full support and real; only 2|3 layers.**
  * At each root, a full-support relation needs the minimal exponent attained twice.
  * The reciprocal relation sum a_i / z_i = 0, from sum a_i zbar_i = 0, needs the maximum attained twice.
  * A 1|4 threshold would make one extreme unique. **CORRECT.**
* **Ramification test.**
  * Real a gives sum a_i (q_i - 1/q_i) = 0, i.e. sum a_i sinh(ell_i) = 0, so sum a_i ell_i = O(t^{-3h}).
  * omega = d(sum a ell) has order >= 3h - 1 at infinity and only simple poles, with residues sum a_i v(q_i).
  * If omega != 0, deg div = -2 gives >= 3h + 1 poles.
  * If omega = 0, a real kernel vector of the integral valuation matrix exists, hence a rational one, and the character
    is constant, equal to 1. **CORRECT.**
  * Force of the test: since #places <= 2D, it only constrains parents with many repeated roots. For squarefree N,
    3h + 1 <= 2D is weaker than the pair-degree bound h <= 3D/5. Astra says this herself.
* None of VII.4 has been tested on an actual five-mark tuple with 2h > D, because none is known. The five-point problem
  is exactly what is open. The statements are conditional structure results.

## (c) Integer analogue

Honestly, none.
* The proof's engine is the anti-holomorphic inversion in a circle with real centre. It swaps infinity with a finite
  place (the centre), so a contact of order h at infinity automatically produces a second contact of order h at a
  finite place, for the self-reciprocal part of the parent.
* Z has no automorphism exchanging the archimedean place with a finite one. "Roots of N on a circle" (factors K with K^*
  proportional to K) has no counterpart among Gaussian primes: the self-conjugate-up-to-unit primes (1+i and inert q)
  correspond to real roots, which are common factors, not to circle roots.
* After specialisation t = T (VII.7), a second contact at the finite place t = c would read z_i(T) ≡ z_j(T) mod
  (T - c)^h. That is lattice points that are both archimedean-close and congruent modulo a large modulus. This is the
  record's residue-class / "twin" problem (round 29, clothes lane: T(N,q), twin boundedness), which is open and is not
  supplied by any symmetry of the integer problem.
* So (VII.4) has no integer analogue, and it gives no input to A_5(delta) over Z.

## (d) Novelty and significance

* **Novelty.** grep shows "second contact" and "self-reciprocal" are absent. "Reversal", "reflection" and "unit circle"
  occur only in other contexts. The closest earlier statements are both Astra's own:
  * (i) the Part II exclusion "ten roots of a degree-one solution cannot lie on a generalized circle", in answer.md
    line 401;
  * (ii) the free thread's "ellipse obstruction" (astra_free_II–VI line 494; README_II "an axis-parallel ellipse with
    real centre gives at most two weight-m directions"), which uses the unit-root identity a_j = a_0 abar_{m-j}, i.e.
    self-reciprocity.

  (VII.4) is a genuinely new quantitative version: no sparsity, arbitrary multiplicities, an explicit off-circle
  defect, and a clean determinant/reversal proof. I know of no literature statement in this form. It is in the spirit
  of Mason–Stothers / Wronskian-type two-point arguments. The VII.4 linear-dependence facts are new to the record as
  far as grep shows ("conjugation-stable", "node" appear nowhere relevant). They are elementary.
* **Significance: low for the goal.**
  * (VII.4) rules out only parents with >= 17/18 of their roots on one circle with real centre, and the threshold
    worsens with delta. Generic parents, and every known construction (Pell, cube, D=4 cluster, random), have r ≈ D,
    where the bound is far from binding (slack >= 32 in all tested five-tuples off the circle).
  * It is a function-field-only mechanism with no integer analogue.
  * VII.4 gives structure (independence of four marks, reality and the 2|3 restriction on a putative relation) that a
    future proof could use, but it excludes nothing by itself.
  * Both sections are correct and honest about their scope. Neither moves A_5(delta).
