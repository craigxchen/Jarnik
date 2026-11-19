# Referee report, round 31, Part III, block "linear"

Scope: `astra_lemma_III/answer.md`, sections III.5, III.7, III.8, III.9 and III.10. These cover the strict linear-block ten-pair
system: the staged elimination, the two-height and symmetry tests, the tangent bound, the global attempt, the census, and
the one-list (barred-list) elimination.

My scripts are in `ref3_linear/`. None of them calls or imports Astra's code.

## 0. Verdict in brief

* **III.8, tangent bound: CORRECT.** I checked it line by line. The two constants it uses, the line-graph determinant
  −192 and the quadric-evaluation determinant −3, I recomputed (`constants.py`). This is the best new result in the block.
* **III.10, one-list test (0 / ≤1 / ≤9 completions at ranks 5 / 4 / 3, rank ≤ 2 impossible): CORRECT.** Line by line,
  including the injectivity argument for the moment map `F`. I reproduced the rank-4 example over F_101 with my own code
  (`onelist.py`). The twenty roots are distinct and the ten scales disagree. My scales differ from Astra's only by the
  free projective normalization of the jet: both ratios λ_01/λ_02 equal 84 mod 101.
* **III.5 census (no strict twisted points for p = 7, 11, 19, 23, 31): CORRECT.** I used two
  independent algorithms of my own and cross-validated them on positive controls.
* **III.5 pivot certificates: CORRECT.** I reproduced all 3 × 36 nonzero resultants (`pivotcert.py`). But they only give
  dimension ≤ 4 over the anchor base, which III.8 already beats.
* **III.5 elimination formulas, finiteness of the pivot fibres, and the two-real-part empty slice: CORRECT.** Checked line
  by line.
* **III.7 five-cycle case: CORRECT.** I re-derived it with my own Gröbner basis over Q(ζ_5), which is [1].
* **III.7 real winding restriction: CORRECT.**
* **III.7 relabelling freeness: CORRECT, conditional on the reflection slice.** The group theory is right. One case (d = 4
  with a swap) is skipped in the text but is easy to fill.
* **III.7 two-height slice and reflection slice: NOT REPRODUCED.** I checked the reductions where I could. I did not redo
  the degree-32 and degree-42 eliminations.
* **III.9: the parts I checked are CORRECT, and they are peripheral.** Checked: the saturation set-up, the unimodular pivot,
  the plane count, the interpolation degree, the consequences of a hypothetical template, avoidance of the quadric cones,
  the positive-definite Hankel test, and the genus-four cube cover. The Gröbner runs are explicitly incomplete and claim
  nothing.
* **III.10 degree bounds and the degree-cutoff equations: CORRECT at the level checked.** I checked the Bézout and
  pole/order counting. I did not check every constant. These results are peripheral.
* **Decision.** I cannot decide the system. It stays open over Q(i), over R and over C.
* **New observation (mine).** At p = 7 the census is already empty using only the first **four** moments. At p = 11 the
  four-moment system has exactly 220 strict twisted solutions; both methods agree on this number.

**Significance: peripheral.** The block decides nothing about the growth rate, the uniform bound, Question A, Question B or
route 5. Even a strict Q(i) linear template would give one curve of approximants at exponent 10/3. Such a curve only
enlarges the exceptional set (expert note, Question B). III.6 already shows that templates of any fixed degree sweep at
most a surface.

## 1. The system and the definitions I used

There are five marks with stars `P_i = ∏_{j≠i}(x − r_ij)` and ten linear blocks. On the split chart the second list
`r'_e` is independent; on the real form it is complex conjugation. The moment equations (II.8) are

  `Σ_{e∋i}(r_e^k − r'_e^k)` independent of i, k = 1..5,

that is, `E(r^k − r'^k) = 0` with E the 4×10 matrix of star differences. **Strict** means all twenty values `r_e, r'_e` are
distinct. As Astra notes (III.9), this already forces every `t_ij ≠ 0`. If some `t_ij = 0`, the two six-factor words are
equal, so their root multisets coincide, which is a collision.

**Twisted F_p-point** (p ≡ 3 mod 4, F_{p²} = F_p[i]): roots `r_e ∈ F_{p²}`, with `r'_e = r_e^p` (Frobenius = reduced
conjugation). Strict then means the 20 values are distinct. Equivalently, the ten roots are nonreal (y ≠ 0) and lie in ten
distinct conjugate pairs. The moment condition is `Im Σ_{e∋i} r_e^k` (the y-coordinate in F_p) independent of i for k = 1..5.

This is exactly what Astra's `p3_anchor_census.cpp` tests: the encoding is x + p·y, `conj_` is y ↦ −y, `moments` stores the
y-coordinate of z^k, and every root is required to have index ≥ p, i.e. y ≠ 0.

## 2. III.8: tangent dimension ≤ 3 (≤ 2 at full second-moment rank). CORRECT.

I went through every step. Here k is algebraically closed, of characteristic 0 or p > 5, and the gauge fixes
`δr_01 = δr'_01 = 0`.

1. **Orders at infinity.** At a solution `log u = O(x^{-6})`, with the x^{-6} coefficient a nonzero multiple of `2i t_{5i}`.
   So `v_0 = O(x^{-7})`, with exact order 7 in some component (this uses 6 invertible). A Zariski tangent vector preserves
   the coefficients of x^{-1}..x^{-5}, so `v_δ = O(x^{-6})`. This is exactly the linearized moment system.
2. **Dependence of v_0, v_δ, v_δ'.** Project modulo `b_e`. The 3×3 determinant has at most simple poles at the 18 roots
   of the other edges, because two singular columns have proportional residues. The residues at the roots of e itself
   die. The determinant vanishes to order ≥ 7+6+6 = 19 at infinity, so it is identically 0. Some `b_e` lies outside any
   3-dimensional span, which gives the dependence. Correct.
3. **No common zero of v_0 off the roots.** Correct. The steps:
   * `(s+it)/(s−it)` has derivative `−2it s'/(s−it)²`, so all `s'_ij(x_0) = 0`.
   * Then `D'_ij(x_0) = 0`. The line-graph adjacency matrix (det −192, recomputed) gives `n'_e(x_0) = 0` for every e, so all
     blocks are centred at x_0.
   * At x_0 the two words take equal values, contradicting `2it ≠ 0`.
4. **Injectivity of δ ↦ W_δ = D·(v_0∧v_δ).** `W_δ` is a polynomial of degree ≤ 20 − 13 = 7. If `v_δ = f v_0`, then:
   * by step 3, f has no pole off the roots;
   * f is regular at the roots, since there `v_0` has a simple pole with residue `b_e ≠ 0`;
   * the orders at infinity give deg f ≤ 1.

   Comparing residues gives `f(r_e) = −δr_e` and `f(r'_e) = −δr'_e`. The gauge then forces f = 0 and δ = 0. So the kernel
   is exactly the infinitesimal affine maps. Correct.
5. **Rank-one module.** The W_δ generate a rank-one torsion-free k[x]-submodule of k[x]^6, which is free with generator G
   of degree g. The multipliers have degree ≤ 7 − g, so dim T ≤ 8 − g. Correct.

   **Incidence:** `G(r_e) ∈ b_e∧V`, because G is a k[x]-combination of the W_δ. So `L_ef = vol(G∧b_e∧b_f)` is divisible
   by `n_e n_f`, and g ≥ 4. Correct.
6. **g = 4.** Correct. The steps:
   * `L = diag(n) K diag(n)` with K constant of rank 2 (G is decomposable, `L = BᵀΩ_G B`).
   * The rows `w^j·diag(n)` lie in rowspace(B) over k(x), hence coefficientwise over k. This gives six linear forms with
     `ℓ_ja(b_e) = w^j_e [x^a]n_e`.
   * The 2×2 minors are quadrics vanishing at the ten `b_e`. The 10×10 evaluation matrix has det −3 (recomputed), so the
     minors vanish identically.
   * A linear space of rank-≤1 2×3 matrices has a fixed row or a fixed column (standard). A fixed left factor makes
     w¹ ∥ w² (via the monic coefficient). A fixed right factor makes two blocks equal.

   Hence g ≥ 5 and dim T ≤ 3.
7. **g = 5.** Correct. `K(x) = K_0 + xK_1` is a decomposable pencil, so it has a common vector w. Then w and `(w_e x_e)`
   lie in rowspace(B). With w = λB, λ is a nonzero left null vector of `B·diag(x_e)·Z`. Full rank 4 therefore forces
   g ≥ 6, i.e. dim T ≤ 2.

Each step uses only the invertibility of 2, 3 and 5 (for the moment–jet equivalence and the order count) and works at an
arbitrary point, singular or not. **The statement is a correct upper bound on the local dimension of the strict scheme
modulo affine maps, in char 0 and every p > 5.**

Context: the naive count is 18 unknowns after the gauge against 20 equations, so the expected dimension is −2. The bound
"≤ 3" is a genuine theorem, but far from the expected answer. It does not exclude finitely many points or curves.

## 3. Census. Own reproduction.

### Method A: `census.cpp` (no gauge, no elimination)

1. Enumerate all oriented 4-stars: four distinct conjugate pairs, one orientation each. That is C(P,4)·16 stars, with
   P = (p² − p)/2 pairs.
2. Bucket the stars by their signature `(Im p_1, …, Im p_K)` ∈ F_p^K.
3. In each bucket, find 5-sets of stars that pairwise meet in exactly one element, have ten distinct shared elements, and
   are pairwise non-conjugate. This is exactly a strict labelling of K_5 up to S_5.

### Method B: `gauged.cpp`

1. Real affine maps act simply transitively on nonreal elements, so assume one edge root is i.
2. Stars 0 and 1 are the two stars containing i. Bucket the i-stars by signature.
3. For each ordered same-signature pair and each of the 6 matchings, look up the inner pairs (r23, r24, r34) in a hash
   table keyed by two-element signature sums.
4. This counts unlabelled solutions containing i. Each unlabelled solution set is counted twice, once per ordering of
   stars 0 and 1.

| p | K | Method A: unlabelled strict solutions | Method B: solutions containing i | consistency (A·10/(p²−p) = B) |
|---|---|---|---|---|
| 7 | 3 (control) | 1638 | 390 | 1638·10/42 = 390 ✓ |
| 7 | 4 | **0** | 0 | ✓ |
| 7 | 5 | **0** (95,760 stars, 16,660 classes, max 14, as in Part II) | 0 | ✓ |
| 11 | 3 (control) | – | 25,290 | – |
| 11 | 4 | **220** | 20 | 220·10/110 = 20 ✓ |
| 11 | 5 | **0** (5,456,880 stars, 160,930 classes) | 0 | ✓ |
| 19 | 5 | – | **0** | – |
| 23 | 5 | – | **0** (10.5M i-stars, 68.3M ordered star pairs, 1m52s) | – |
| 31 | 5 | – | **0** (16 signature passes, 603.5M ordered star pairs, 19m) | – |

`verify_sol.py` rechecks the printed control solutions directly. They are genuine:
* the K = 3 solution at p = 7 satisfies k = 1, 2, 3 and fails k = 4;
* the K = 4 solution at p = 11 satisfies k = 1–4 and fails k = 5.

So the search code does find solutions when they exist. Methods A and B agree wherever both ran. Astra's census, by a third
algorithm (anchor 0–4 and pivot b), also reports zero for 7, 11, 19, 23, 31.

**New (mine): at p = 7 even the four-moment system has no strict twisted point.** Part II's F_49 statement therefore holds
with k ≤ 4. At p = 11 four moments do not suffice: there are 220 unlabelled solutions, i.e. 20 per gauge.

### What the census implies

**Over Q(i).** Let r_e ∈ Q(i) be a strict solution and p ∈ {7, 11, 19, 23, 31}, an inert prime.
* Every real rational affine change that makes all roots p-integral reduces the system to a twisted F_p-solution. The
  moment equations have integer coefficients, and conjugation becomes Frobenius.
* The census then forces a collision mod p among the 20 values. Such a collision may be `r_e ≡ r̄_e`, i.e. a root that
  becomes real.
* Equivalently: for each of these primes, the 20 points `r_e, r̄_e` are **not p-adically equidistant**. Take the
  normalization that makes the minimal distance a unit; some other distance is then a non-unit.

This restricts the p-adic cluster picture at five small primes. It does not exclude a Q(i)-solution. A rational solution
has only finitely many bad primes, but nothing prevents these five from being among them.

**Over C (and over R).** The census implies nothing. Twisted F_{p²}-rational points are not geometric points. If the
geometric scheme were finite and nonempty, its points would be defined over some number field. They could fail to reduce
to F_{p²}-rational strict points at all five primes. And no Lang–Weil argument applies at such small p.

Astra's own sentence "These point counts alone do not decide geometric emptiness or dimension" is accurate.

### III.5 pivot certificates (reproduced, `pivotcert.py`)

I wrote my own F_{p²} arithmetic and computed the following:
* the bidegree-(2,2) divided differences W_j;
* the homogeneous resultants E_12 and E_23, both of exact degree 8 for all matchings;
* the homogeneous 16×16 resultant `Res_b(E_12, E_23)`.

For all 36 (v, v̄)-matchings at the three stated base points (p = 7, 11, 19) this resultant is nonzero. I also checked
that the base points have equal imaginary power sums in degrees 1–5 and 12 distinct roots, and that `ŪV − UV̄` is a nonzero
constant.

The deduction is correct: the product over matchings is a regular function on the 5-dimensional anchor base, nonzero mod p
and hence nonzero in char 0. Strict solutions therefore lie over a proper hypersurface, so the strict scheme has dimension
≤ 4. That bound is **weaker than III.8**. The certificates only say that the fibre over a generic anchor is empty.

## 4. III.10: one-list elimination. CORRECT.

* **Jet reformulation.** The moment conditions hold if and only if `P'_i/P_i ≡ 1 + λa(x) mod x^{-6}`, with a common jet.
  Then `P'_i = P_i + λR_i`, with `R_i = (aP_i)_{≥0}`, and the x^{-1} condition `Σ_k p_ik a_{k+1} = 0` holds. Both
  directions are correct.
* **Rank ≤ 2 is impossible.** The other three quartics would vanish at the shared root u_ij.
* **The polynomials C_ij.** Each `C_ij` has degree ≤ 1. In a solution it equals `±2it_ij·λ^{-1}(x − v_ij)` up to sign
  conventions, so its degree is exactly 1. The converse via common λ is correct.
* **Rank 4.** The barred list is unique.
* **Rank 3.** The cross-product forms have degree 4 + 5 = 9.
* **Injectivity of the differential of F(v) = (B v^k)_{k ≤ 5} on distinct inputs.** Correct. For s ≤ 5 the polynomial C
  has negative degree. For 6 ≤ s ≤ 10, quadrics evaluated on C(x) span at most 2s − 11 < s values, against rank s from det
  −3.

  Hence strict fibres of F are finite. If all cross-products vanished on the kernel line, an open set of directions would
  give infinitely many barred lists with the same F. So ≤ 9 completions. Correct.
* **Computed rank-4 example (F_101).** I reproduced it independently.
  * The rank is 4, the ten C_ij are linear, and the 20 roots are distinct.
  * The barred list is (97, 62, 63, 95, 50, 44, 43, 22, 16, 73).
  * The scales are (11, 35, 27, 42, 92, 55, 80, 2, 85, 18), which are not equal.
  * Ratio λ_01/λ_02 = 84 = 17/100 mod 101, so this matches Astra's normalization.
* **Component-degree bound 120⁴.** Correct; this is the standard refined Bézout induction.
* **Degree-cutoff and fifteen-cut cutoff.** The pole/order counts are right:
  * ≤ 2D simple poles, against order ≥ 4k + 1 > 2D at infinity;
  * coefficient degrees 6 or 8;
  * graph degree 7 or 9;
  * resulting degree O(6^{24d_0}), resp. O((d_0+1)^{15}·8^{30d_0+4}).

  I did not audit every constant. **Peripheral**: it only controls tuples that come from polynomial families with these
  contact properties. Astra says so.

## 5. III.5 elimination and the empty slice: CORRECT

* **The cubic identity** `(X − b)C_j = A_jL_j + (1 − A_j)R_j`. It comes from `t_04 n_04 P_j = t_j4 n_j4 P_0 + t_0j n_0j P_4`
  after cancelling `(X − a)(X − u_j)(X − v_j)`. Then `A_j = R_j(b)/(R_j(b) − L_j(b))`, W_j is symmetric of bidegree
  ≤ (2,2), and E_jk has degree 8. All verified.
* **Nonvanishing of E_jk.** Correct. A joint degree d | 3 forces d = 3. Then the two Möbius-related cubic functions share
  the zero v_l, the pole u_l and the value 1 at ∞, so they are equal. That contradicts their different zero sets.
  Separability holds since p > 3.
* **One or two real parts.** Correct.
  * Two centres: each colour has constant sums and cube sums. A nonempty colour misses no vertex (by the K_4 and 3xyz
    argument). Leaves are excluded (by 3y_1y_2b and the (y_1+y_2)(y_2+y_3)(y_3+y_1) identity). The odd-cycle argument
    finishes.
  * One centre: the six purely imaginary factors give a real value, so t = 0.

## 6. III.7

* **Two heights.**
  * Checked: the first moment forces a regular graph; the one-colour case reduces to equal power sums of abscissas (as in
    II.5).
  * Not re-run: the remaining elimination (degree-32 P(x), gcd with its reversal equal to x − 2, plus four primes).
  * Formulation that would do it: with y eliminated linearly through the invertible circulant 1 + S², the unknowns are
    t_2, t_3, t_4, α, β = B², c, c′, c″, giving 10 equations in 8 unknowns. I set this up, but it is beyond sympy without
    Singular or msolve. Installing those was not permitted in this session.
* **Reflection slice.** The abscissa reduction (coboundary formula, D-operator) matches III.6. The degree-42 elimination
  and the golden-ratio point were **not reproduced**.
* **Relabelling freeness.** Correct, given that slice and the five-cycle case.
  * Orders: ord σ = ord μ.
  * No swap: gcd(d, 6) = 1, so d = 5.
  * Swap with odd d: μ^d = −1, impossible.
  * Swap with d = 4: −λ² has order ≤ 2. **This case is omitted in the text, but trivial.**
  * Swap with d = 2: λ = 1 forces a fixed edge with u_e = v_e; λ = −1 is a double transposition, i.e. the reflection
    slice.
* **Five-cycle case.** I re-derived the four bracket equations `G_k(A_1^k − B_1^k) + (A_2^k − B_2^k) = 0`, k = 1..4. The
  lex Gröbner basis over Q[z]/Φ_5, with A_1 − B_1 = 2 (legitimate by homogeneity), is **[1]** (`fivecycle_check.py`). This
  confirms Astra's cyclotomic contradiction for both primitive roots.
* **Real winding restriction.** Correct. Each factor has degree ±1, the ratio `(s+it)/(s−it)` with s real monic of even
  degree has winding 0, so the graph is regular of degree 0, 2 or 4, and degree 2 means C_5.

## 7. III.9

* **Saturation set-up.**
  * The gauge (u_01, v_01) = (0, 1) is a global slice.
  * The pivot minor of `[E; e_01]` on the first five lexicographic edges has det 1 (recomputed).
  * There are 14 free variables and 16 equations of degrees 2, 3, 4, 5 (four of each).
  * Strictness implies t ≠ 0.

  All correct. The FLINT and Buchberger runs are incomplete and claim nothing.
* **Plane count O(P^10 log P).** Correct, via a reduced basis and dyadic shells; the true order is P^10 (Schmidt).
* **Interpolation degree** `O(1 + (AK)^{1/4}L^{(1+ε)/4})` and the condition α + 10κ + ℓ < 4: correct but conditional.
* **Consequences of a rational template.** Columns 12, 13, 14, 25 have det −1 (recomputed); the CRT index-one
  construction is plausible. Avoiding the 15 quadric cones follows from a double pole against a regular side, which is
  correct. All of this is conditional on a template that nobody has.
* **Positive-definite Hankel test.** Correct.
  * F = U/V with U, V monic cubics having roots in the upper half-plane.
  * U·V̄ is the reduced word, so Im F = t/|V|² > 0.
  * The Schwarz integral holds on the lower half-plane. The expansion to order z^{-5} works because θ = O(x^{-6}): the
    remainder is O(|z|^{-6} log|z|).
  * The real parts of the coefficients are (m_{i,k} − m_{j,k})/k.

  A neat necessary inequality, but it excludes nothing yet.
* **Cubic cover** `z³ = kx∏(x − q_i)`: there are 6 totally ramified points, so g = 4 and Faltings gives finiteness. Correct;
  peripheral.
* **Degree-two heuristic.** It is labelled HEURISTIC and stays so.

## 8. What is decided and what is not

**Proved.** A strict linear template, over any algebraically closed field of characteristic 0 or > 5, has:
* modulo affine maps, local dimension ≤ 3, and ≤ 2 where `B·diag(x)·Z` has rank 4;
* at most 9 barred completions for any unbarred list;
* no one- or two-real-part, two-height, reflection or five-cycle symmetry. The two-height and reflection cases rest on
  Astra's unreproduced eliminations.
* projection onto a proper hypersurface of the anchor base.

**Computed and reproduced.** There are no strict twisted points over F_{p²} for p = 7, 11, 19, 23, 31, and at p = 7 none
even with four moments.

**Open.**
* Geometric emptiness of the strict scheme. Expected, since the system is overdetermined by 2, but unproved.
* Existence of a real point (the real form).
* Existence of a Q(i)-point.

My own numerical probe gives evidence only. In 200 Levenberg–Marquardt runs on the complex split system, 162 converged, and
every converged point is degenerate (relative minimal root separation ≤ 7·10⁻¹¹; `lm_split.py`). This agrees with the
record's 450 real-form runs (`growth/cotangent.md` §3).

The natural decisive computation is the III.9 saturation in Singular or msolve over a few primes and Q, or the eight-variable
rank-4 chart of III.10. I could not run either here.

## 9. Novelty and significance

* **Novelty.** The record (`growth/cotangent.md` §3) states the linear system and has numerical searches; Part II added
  the F_49 census, which `ref2_window.md` verified. Everything in III.5–III.10 is new this round:
  * the tangent bound;
  * the barred-list elimination with ≤ 1/≤ 9 completions;
  * the pivot elimination;
  * the censuses at p ≥ 11;
  * the symmetry slices;
  * the Hankel inequality.
* **Significance: peripheral.**
  * The linear-block system is one test case, of fixed degree d = 1, inside the polynomial-template programme. Deciding it
    either way changes neither the (2/3 + ε) log R / log log R bound, nor the uniform M(C) question, nor Question A,
    Question B or route 5.
  * A solution would give a single curve, and by III.6 ≤ 2-dimensional families. These only enlarge Z. The expert note
    says a refutation of (L_5°) needs a Zariski-dense supply.
  * Emptiness would kill only d = 1.
  * The tangent bound is a clean structural result, but it is an upper bound far above the expected (empty) answer.

## Addendum: runs at p = 23 and p = 31

Method B (`gauged.cpp`) finished with **0** strict twisted solutions at p = 23 (`out_p23_K5.txt`) and at p = 31
(`out_p31_K5.txt`, run in 16 hash passes to stay under 3.5 GB). Every prime Astra listed has now been reproduced
independently.

Output files in `ref3_linear/`: `out_p11_K5.txt`, `out_p11_K4.txt`, `out_p23_K5.txt`, `out_p31_K5.txt`. The p = 7, 11 (K ≤ 4)
and p = 19 counts in the table were printed to the terminal by the commands documented at the top of each source file.
