# Referee report on round16/reduction.md

Setting throughout: `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime `d = n - 1 >= 316 ln 4r`,
`H = U^2` (record) or `H' = e^{m'W/40}` (t1rigour Remark 3.3), any profile, any norm. All scripts are my own and are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_reduction_checks/`.
I did not import the author's `lib.py`.

## 0. Verdict

**The result survives, but it needs corrections. No growth improvement and no uniform bound are claimed, and none
follows.** `M <= (2/3 + eps) log R/loglog R` stands.

I re-derived every proof. Each theorem the reduction uses is correct in the direction it is used:
* Thm 2.1, in the direction "probabilistic form ⇒ `L^2` form ⇒ sup form";
* Prop 3.1 and Cor 3.2;
* Thm 4.1;
* Thms 5.1, 5.2 and 5.4, Lemma 5.5, Cor 5.6, Prop 5.7 and Thm 5.8;
* Prop 7.1.

My own exact computations confirm the cube formula, the product formula (also with a non-slab factor) and the
numbers in Prop 7.1. The defects are in the framing statements:

* **(M1, major: overclaim in a "proved" item.)** The text says "all probabilistic versions are the same statement up
  to `e^{O(d)}`". This includes any probability threshold, "most `k`", and the cube measure with threshold `mu(A)`.
  It is **not proved**: Thm 2.1(3)-(4) gives only one direction. It is **false in general**. I give exact toys in
  which `MIN^pair_inf` holds while `AC_p` (for every `p > 1/7`) and `AC_cube` fail at every level `delta >= 0` (§2). For
  regime boxes the converse is open. Only the following are equivalent up to `e^{O(d)}`:
  * the sup form;
  * the `L^2(mu_A)` form;
  * the `p <= 1/U` ("some `k`") small-ball form.

  The "most `k`" and cube forms are sufficient conditions that may be strictly stronger. This actually strengthens
  the author's main message: anti-concentration is the wrong frame.
* **(M2, major: quantifier error in the minimal open statement.)** §8 states `(R_H) <== (R_inf) + (P2)`. This is false
  as written. (P2) loses a factor `e^{-o(n ln n)}`, and `(R_inf)` as defined has no margin to absorb it. A matched
  margin is needed (§4).
* **(M3, minor-to-major: unproved "negligibility" of `e^{O(d)}`.)** The claims "the four statements are
  interchangeable for the growth question" and "`(R_H)` *is* MIN^pair for the growth question" are heuristics, not
  theorems. The Step-0 level `eps_0` is fixed by the cluster. The "room" `e^{(3/4) n ln n}` is the gap to the
  Chebyshev ceiling, not slack in any implication. "Shift `W` by `O(d)`" also moves the gap region `G`. The lossless
  directions are the ones that matter for a proof, and they are correct.
* **Minor:**
  * "(P1) holds wherever the certificates of §5 apply" is unproved. The certificates prove the *conclusion*
    (`kappa_H` lower bounds on `G_1`), not (P1). I supply a proof of (P1) for the uncut cube.
  * `r_0 = e^{o(d)}` needs `ln r = o(d)`.
  * The hypercontractivity constant is quoted from memory and could not be confirmed.
  * `n = 2901` is not in `d >= 316 ln 4r`.
  * Influences of `0.16`-`0.32` are called "`-> 0`".
  * Remark 4.2 has a factor-`m'` slip.
  * The few-cuts constant `-5.6` is asymptotic, not "explicit room".
  * "Regime slabs are thick" is argued only for balanced sign columns.

I found no fake-pattern-box counterexample to any proved statement. Every proved positive case (cube, one cut, few
cuts, thin cuts) needs twin classes or thin slabs, and gap-admissible fakes have neither (Prop 5.7 re-derived; §7).

## 1. Thm 2.1 / Cor 2.2 (the reduction to Cheb_A and Psi_A)

I re-derived each item of Thm 2.1 for arbitrary finite `A` in the ternary cube.

* **(1) Duality.** This is `l^1`-`l^inf` duality: `sup_{c in V} |<c,g>|/|c|_1 = dist_inf(g, V^perp)`, with
  `V^perp = Poly_{<m}|_A`. The normalisation `kappa_c = |Phi_c|/(|c|_1 d^m)` matches the record (round12 README §1).
  **Correct.**
* **(2) Sup versus `L^2`.**
  * `Psi_nu <= Cheb_A^2` holds because `supp nu` lies in `A`.
  * `Cheb_{supp nu}^2 <= Psi_nu / min nu` follows from `E_nu Q^2 >= min nu * max_{supp} Q^2`.
  * For `mu_A`, these give `Psi_A <= Cheb_A^2 <= U Psi_A`.

  **Correct.**
* **(3) First part.** `SB_nu <= 1 - p` gives `nu(|Q| > delta) >= p` for every admissible `Q`, hence `E Q^2 > p delta^2`.
  The inequality is strict because `Q^2 > delta^2` on a set of positive mass, and the minimum is attained. **Correct.**
* **(3) Second part.** `Cheb_A > delta <=> SB_{mu_A}(u; delta) < 1` is **correct**. The compactness argument is not even
  needed: `mu_A`-probabilities are multiples of `1/U`, so `Cheb_A > delta` gives `SB <= 1 - 1/U`.
* **(4) and (5).** **Correct.** In (5), `min_{f in V^perp} |g - f|_2^2 = |P_V g|^2 = sum_i <c_i, g_u>^2`, and the zero set is
  `TC_V`.
* **Remark 2.3 (Gaussian identity).** `Psi_{N(0,I)}(u) = m! |u|_2^{2m}` is correct. The `m`-th Wiener chaos component
  of `|u|^m g^m` is `|u|^m He_m(g)`, and `E He_m^2 = m!`.
* **Remark 2.4.** **Correct.** Injectivity of the full statistic on sum-zero vectors when the columns together with `1`
  span `R^n` is immediate.

**Cor 2.2.** The lossless sufficient directions are all correct:
* `AC_2 ==> MIN^pair_inf`;
* `AC_cube ==> MIN^pair_inf`;
* `AC_p ==> AC_2`, with the strictness discussed above.

The converse `MIN^pair_inf ==> AC_2` at level `eps_0 3^{-d/2}` is also correct. **The final sentence is not proved:**
"the four statements are interchangeable for the growth question (each at level `eps_0` implies the others at level
`eps_0 e^{-O(d)}`)". Specifically:
* (a) `MIN^pair_inf ==> AC_p` (fixed `p`) and `MIN^pair_inf ==> AC_cube` are not implied by anything in Thm 2.1, at any
  level (§2).
* (b) Even for `AC_2`, "at level `eps_0 e^{-O(d)}`" is not something the growth argument can use. Step 0 gives
  `kappa <= eps_0` with `eps_0 = C d e^{-W/4}/(m'+1)` fixed by the cluster. A statement at the lower level
  `eps_0 e^{-O(d)}` does not contradict it.

  The proposed fix "replace `W` by `W - 4.4 d`" also changes the gap region. `eps_xy = d_xy - W/2 + 2 ln C'` grows
  by `2.2 d`, so `G` grows. One would need MIN^pair for *modified* cluster data, and no theorem says regime data are
  closed under this shift.

  So the `e^{O(d)}` losses are "negligible" only in a heuristic sense: the room `e^{(3/4 - o(1)) n ln n}` measures how
  far the target lies below the Chebyshev ceiling. **Correct reading:** `AC_2` (and `(R_H)`) is a sufficient condition
  that is *at most* `e^{O(d)}` stronger than MIN^pair. It is not shown to be equivalent to it for the growth question.

**Exceptional set and separation.** For a collision pair `(x, y)` of `L`: if `v in L`, then
`|u_x - u_y| = |(u - v)_x - (u - v)_y| <= 2 |u - v|_inf`. So `dist_inf(u, L) >= (1/2) e^{-min_coll eps/2}` on `G`. The
minimum over the `P(L)` collision pairs is at most their average, which is at most `N_2 abar/P(L)`. So the per-pair
radius is at least the averaged radius of nearstrata L3.3. **Correct.** The statement "the exceptional set is exactly
`R^d \ G`" holds by definition, since MIN^pair is a minimum over `G_1`.

## 2. Probabilistic versions and "scale is absolute"

**Claim under test** (summary item 1, Cor 2.2, first row of the §8 status table): "All probabilistic versions (any
probability threshold, 'all or most `k`', uniform on `A` or the cube with threshold `mu(A)`) are the same statement up
to `e^{O(d)}`."

**Refutation in general (exact; `check_prob_forms.py`, `check_prob_forms_out.txt`).**

Setting: `d = 2` (`n = 3`, `u_3 = 0`), `m = 3`. Take the slab box `A = {k' in {-1,0,1}^2 : |k_1 + k_2| <= 1}`: one sign
column `beta = (1,1,-1)`, `|<k,beta>| <= 2`, 7 points. Take `u = (1/2, 1)`. Then `|u|_inf = 1`, and the pairwise gaps
are `1/2, 1/2, 1`, so `u` lies in `G_1` whenever `eps_12, eps_13 >= 2 ln 2`.
* **`MIN^pair_inf` holds.** `Psi^Sigma_A(u) = 3/8 > 0`. Exactly on an integer grid, `TC_V(A,3)` lies in the braid
  arrangement `u_1 u_2 (u_1 - u_2) = 0`. Since `J(A)_3 = span{k_1 k_2 (k_1 + k_2)}`, `Cheb_A > 0` on the compact set `G_1`,
  so `MIN^pair_inf` holds for all small `eps_0`.
* **`AC_cube` fails at every `delta >= 0`.** Take `B` = the cube minus `{(1,1), (-1,0)}`, also 7 points. It has
  `delta_{(1,1)} + (1/2) delta_{(-1,0)} in I(B)`, with leading form `(1/4) k_1 k_2 (k_1 + 2 k_2)`. This is proportional to
  the top component `6 k_1 k_2 (k_1 + 2 k_2)` of `l_u^3`. Exactly, `Psi^Sigma_B(u) = 0`. So some `P` of degree `< 3` has
  `l_u^3 = P` on `B`, and `SB_mu(u; 0) >= 7/9 = mu(A)`. `AC_cube` demands `< mu(A)`.
* **`AC_p` fails for every `p > 1/7` at every `delta`.** `A` minus `(-1,0)` (6 points) has `Psi = 0` at `u`, so
  `SB_{mu_A}(u; 0) >= 6/7`.
* **Same phenomenon at top degree.** For the uncut cube `{-1,0,1}^2` at `m = 4`, `Psi_cube(u) = 1`, while the cube minus
  one point has `Psi = 0` for every `u`.

**What is true.**
* The sup form is the same statement as the `L^2(mu_A)` form, up to `U <= 3^d`.
* Both are the same as the small-ball form with threshold `p = 1/U` ("some `k`", level `sqrt U` times larger), up to
  `e^{O(d)}`.
* Every "most `k`" form (`AC_p` with `p > 1/U` fixed) and the cube form `AC_cube` *imply* `AC_2` or `MIN^pair_inf`.
  The converses fail in toys and are open for regime boxes. In the regime, `m' <= m_max` gives
  `r_0 <= 3^d/8 <= U/4`, so the trivial mechanism above (a subset of `A` with interpolation degree `< m'`) needs mass
  `<= 1/4` and cannot break `AC_{1/2}`. Whether a structured half-mass subset `B` of `A` with
  `TC(B, m') ∩ G ≠ ∅` exists is undecided.
* For `AC_cube`, the competing sets `B` range over *all* subsets of the cube of mass `>= mu(A)`. These include
  permuted and fake slab boxes of the same mass. So `AC_cube` is a statement about all such sets, not about `A`.

**"Scale is absolute" and "shape plays no role".** Both are correct for the sup and `L^2` forms:
`max_A |Q*| >= ||Q*||_{L^2(mu_A)}`. The correction above sharpens the author's conclusion. Probabilistic phrasings
are at best equivalent (threshold `1/U`) and otherwise possibly *stronger* than needed.

## 3. Prop 3.1-3.3 (slab ideal, q_A, restriction ratio)

* **Prop 3.1. Correct.**
  * `gr(phi) = 0` iff `phi|_A in Poly_{<m}|_A` iff `phi` is the degree-`m` leading form of an element of `I(A)`, so
    `ker gr = J(A)_m`. Surjectivity is clear.
  * For degree-compatible orders, `in(I(A)) = in(J(A))`.
  * The normal-form formula is standard. The slab polynomials `prod_{|t|<=j}(k(S_p) - t)` and the relations
    `k_x^3 - k_x` lie in `I(A)` and enter only through `J(A)_m`.
* **Cor 3.2. Correct.** The arithmetic `a + 2b < W/2 - 2 ln(C d^{m'+1}/(m'+1))` checks. The integrality remark (iii) is
  correctly labelled as useless, and its bound is consistent with Hadamard plus `lambda_max <= tr`.
* **Prop 3.3. Correct as a sufficient condition.**
  * `rho_A = mu(A) Psi_A/Psi_cube` and `Psi_A >= rho_A Psi_cube`.
  * The factor bound `d^{2m}/Psi_cube <= exp(2m ln(ed/m) + 1.46 m + m abar_n)` is re-derived in §5 below.
  * At `m'` this factor is `e^{O(m'(ln n + abar_n))}`. With `abar = W/(2d) + 2 ln C'`, `m' abar = O(W sqrt(ln r/d)) = o(W)`,
    so the threshold is `e^{-(1/2 - o(1))W}`.
  * The informal phrase "equivalently a restriction-ratio statement" holds only with the `u`-dependent threshold
    `mu(A) eps_0^2 d^{2m}/Psi_cube(u)`. The uniform `e^{-(1/2 - o(1))W}` version is sufficient, and is necessary only up
    to `e^{o(W)}`. This is the same caveat as M3.

## 4. Thm 4.1 (passage H = infinity -> V_H; Bombieri-Vaaler) and the §8 splitting

**Thm 4.1. Correct, modulo the quoted theorem.**

* **The quoted theorem.** The Bombieri-Vaaler statement used is: an `r_0 x U` integer matrix of rank `r_0` has
  `U - r_0` independent integer kernel vectors with `prod |x_i|_inf <= sqrt(det(M M^T))/g`.
  * My WebSearch found secondary restatements that agree with this shape: the sup norm, `sqrt(det(AA^T))`, and `g` the
    gcd of the maximal minors.
  * They disagree on the theorem number ("Thm 1" or "Thm 2"). The citation should be "Bombieri-Vaaler 1983, main
    theorem", with the number unconfirmed.
* **The argument.**
  * Hadamard gives `det(M_0 M_0^T) <= U^{r_0}`, since the rows are `{-1,0,1}`-valued of length `U`.
  * Each nonzero integer vector has `|x_i|_inf >= 1`. So fewer than `r_0 ln U/(2 ln H)` of the `x_i` exceed `H`, and
    the rest lie in `V_H`.
  * `ker(Phi|_V) = V(A, m+1)` holds because `c` annihilates all degree-`m` forms.
  * `Phi(V) = J(A)_m^perp` holds because `Phi(V)^perp` under the apolar pairing is `{f : f|_A in Poly_{<m}|_A} = J(A)_m`.

  All three **check**.
* **Minor 1: `r_0 = e^{o(d)}`.** The bound `binom(m'+d-1, d) = e^{O(w ln d)}` is `e^{o(d)}` only if `ln r = o(d)`. At
  the regime edge `ln 4r = d/316`, `m' ≈ 0.636 d ≈ m_max`, and `r_0` is of order `3^d/8`. The relative codimension is
  then `O(ln U/ln H)`, not `e^{-Omega(d)}`.

  In the growth range this is harmless. Split primes are `>= 5` and `sum_p e_p log p <= W`, so `r <= W/ln 5`. With
  `W = O(n ln n)` this gives `ln r = O(ln n)`. The text should say so. MIN^pair "for all `W > W_1`" allows larger `r`.
* **Minor 2: Remark 4.2.** A circuit of size `s` has Hadamard height `<= (s-1)^{(s-1)/2}`. This lies below
  `H' = e^{m'W/40}` for `s` up to about `0.15 m' n` at `W ~ 3 n ln n`, not `0.15 n`. The slip is harmless, since
  circuits have `s >= 3^{m'/2}` anyway. I re-derived that bound: a standard monomial `k^alpha` with `alpha <= 2` has
  `prod(alpha_x + 1) >= 3^{|alpha|/2}` divisors.

**(P1) "holds in the uncut cube, in one-cut boxes and wherever the certificates of §5 apply". Partly unproved.**
* (a) As written, `Phi(V_H) = Phi(V)` cannot hold literally: `V_H` is a `Q`-span, and `Phi(V)` is a real space. One
  should use the real span. `kappa_H` is unchanged, by continuity.
* (b) **The uncut cube: (P1) is true, and here is the proof missing from the text.**
  * For `c in V`, `Phi_c = sum_{alpha <= 2, |alpha| = m} a_alpha u^alpha` with `a_alpha = (m!/alpha!) sum_k c_k k^alpha`. Terms
    with some `alpha_x >= 3` reduce to degree `< m` and vanish against `c`. So `|a_alpha| <= (m!/alpha!) |c|_1`.
  * The product certificate `c_alpha` takes a second difference `(1,-2,1)` on `alpha_x = 2` and a first difference
    `(-1,0,1)` on `alpha_x = 1`. It has `|c_alpha|_1 = 4^{n_2} 2^{n_1}` and `Phi_{c_alpha} = m! 2^{n_1} u^alpha`.
  * Set `c' = sum_alpha a_alpha c_alpha/(m! 2^{n_1(alpha)})`. Then `Phi_{c'} = Phi_c` and `|c'|_1 <= K |c|_1` with
    `K = sum_alpha 2^{n_2(alpha)} = e^{O(m ln d)}`, which is negligible.
* (c) **One cut and few cuts: unproved.** Cor 5.6 gives `kappa_H` lower bounds *on `G_1`* from pair-zonotope
  certificates. It does not span `Phi(V)` with controlled norms. The correct statement: in these boxes the
  *conclusion* of the passage (MIN^pair with `kappa_H`) holds directly. The §8 status row "automatic for product
  certificates" should be read the same way.

**The §8 splitting (major, M2).** §8 defines:
* `(R_inf)`: `rho_A(u) > eps_0^2 d^{2m'}/Psi_cube(u)` on `G_1`. This gives exactly `Psi_A(u) > d^{2m'} eps_0^2`.
* `(P2)`: `Psi_{A,H'} >= e^{-o(n ln n)} Psi_A` on `G_1`.

Together they give only `Psi_{A,H'} > e^{-o(n ln n)} d^{2m'} eps_0^2`. **That is not `(R_H)`**, and it does not
contradict Step 0, whose level is fixed (M3). The informal gloss "an `e^{-(1/2 - o(1))W}` fraction" hides the
needed margin. The splitting must use a matched margin. For a function `lambda(n) >= 0`:

```text
(R_inf^lambda)  rho_A(u) > e^{lambda} eps_0^2 d^{2m'}/Psi_cube(u)  on G_1,
(P2^lambda)     Psi_{A,H'}(u) >= e^{-lambda} Psi_A(u)               on G_1.
```

Then `(R_inf^lambda) + (P2^lambda) ==> (R_H)`, which is now correct. For example, with `lambda = delta W` for a fixed
`delta > 0`, `(R_inf^lambda)` is the restriction ratio `> e^{-(1/2 - delta)W + O(m'(ln n + abar_n))}`.

Also, `kappa_H d^m >= Psi_{A,H}^{1/2}` and `Psi_{A,H} = Psi_A - |P_F g_u|^2/U` both re-derive correctly, and
`(R_H) ==> MIN^pair` is lossless.

## 5. Proved cases (§5)

**Thm 5.1 (uncut cube).**
* **(a) Formula. Correct.** I checked it exactly with my own least-squares code over `Q`
  (`check_cube_product_ref.py`): `d = 2, 3, 4`, every `m <= 2d + 1`, random rational `u`. There were 63 cases and 0
  mismatches, including `Psi = 0` for `m > 2d`.
* **(b) Lower bound. Correct.**
  * The kept term is `(m!/2^{m/2})^2 (2/9)^{m/2} prod_X u_x^4`.
  * On `G`, `prod_X u_x^4 >= e^{-2 sum_X eps_xn}`, and `sum_X eps_xn <= (m/2) abar_n` (the smallest `m/2` of `d` values).
  * `ln 2 + (1/2) ln 4.5 = 1.445 < 1.46`.
  * *Nit:* for odd `m` the kept term `alpha = 2*1_X + e_y` carries an extra factor `u_y^2 >= e^{-eps_yn}`, which the
    remark "factor 3/2" omits. This is harmless.
* **(c) MIN^pair for the cube box. Correct.**
  * The certificate is the `alpha = 2` zonotope, with `|c|_1 = 2^m` and `Phi = m! prod u_x^2`.
  * `W/4 = (d/2)(abar - 2 ln C')`.
  * The condition `abar >= (2/(1-theta))(theta ln(2e/theta) + ln C') + O(ln d/d)` re-derives.
  * The base-point averaging is correct: `(1/n) sum_y abar_y = (2/(nd)) sum eps <= abar`.

**Thm 5.2 (products). Correct.** The decomposition
`Poly_{<=m}|_{XxY} = (+)_{a+b<=m} H_a(X) (x) H_b(Y)` and the top-component computation re-derive. I checked it exactly
with `X` the one-cut block (`s = 3`, `L = 1`) and **`Y` an irregular 6-point subset of `{-1,0,1}^2`** (not a slab box,
which is more general than the author's check). There were 21 cases (`m <= 7`) and 0 mismatches. *Nit:* "disjoint
signed cuts give products" needs representatives of the sides that avoid the base coordinate `n`. Otherwise the
blocks are coupled through `k_n = -sum k'`.

**Thm 5.4, Lemma 5.5, Cor 5.6.**
* **Thm 5.4. Correct.** Disjoint pairs give ternary coordinates, and `max_t |sum t_i Delta_p| = sum |Delta_p|`.
* **Lemma 5.5. Correct.** For odd `|C|`, the inclusion probability is `1/|C| <= 1/(|C| - 1)`.
* **Cor 5.6(a). Correct.** Here `s_* >= (n-2)/2`, and the exponent is `(m/2) abar n(n-1)/((c_0-1)(n-2))`, about
  `m abar` for a balanced cut. The condition `m <= (1 - eps) n/2`, `abar >= C(eps)(1 + ln C')` re-derives.
* **Cor 5.6(b). Correct asymptotically.** The `(1 + O(2^q/n))` factor and the threshold
  `2^q < sqrt(d/(8 ln 4r))/16`, i.e. `q < (1/2) log_2(d/ln 4r) - 5.5`, check. Two caveats:
  * The margin from `-5.5` to `-5.6` is a factor `2^{-0.1}`. It then needs `abar` above roughly `30 (ln C' + ...)` to
    beat `d ln C'`.
  * So "explicit room" is too strong: the statement holds for `n >= n_0(C')`, with `n_0` not computed.
* **Common caveat.** These are model statements: by Prop 5.7, twin classes never occur in actual or gap-admissible
  data.

**Prop 5.7. Correct.** I re-derived it with the record's definitions (nearstrata §1:
`a_{p,x} - a_{p,y} = (e_p/2) Delta_p(x,y)`, `d_xy = sum_p log p |a_{p,x} - a_{p,y}|`, `sum_p e_p log p <= W`):
* `sum_i d_{x_i y_i} <= (w/2) W`;
* `sum_i d_{x_i y_i} >= |M|(W/2 - 2 ln C')`;
* hence `|M| < 2w` once `W > 8 ln C'`.

Theorem 5.4 needs `m <= 2|M| < 4w <= m'`, so it never applies. The exclusion covers pair zonotopes with any
multiplicities (with `alpha_i = 1`, `|M| >= m` is needed). Other product certificates are excluded only by the cited
nearstrata Remark 2.10 and fatpoint Thm 3.1, which I did not re-referee.

**Thm 5.8 (thin cuts). Correct modulo the quoted constant.**
* `E[Q^2 1_{A^c}] <= ||Q||_4^2 eta^{1/2}` is Cauchy-Schwarz. The Efron-Stein degree is at most the polynomial degree.
* I could not confirm O'Donnell Thm 10.21 with the constant `(sqrt(q-1) lambda^{1/q-1/2})^k`. A secondary source found
  by search quotes it in the weaker form `(10q/lambda)^k`.
* With the weaker form the threshold becomes `eta <= 120^{-4m}/4`. Every conclusion survives: thin cuts still give
  `rho_A >= 1/2`, and regime slabs are thick a fortiori.
* Sanity check on the constant: indicators of atoms give a ratio of `3^{m/8}`, far below `3^{3m/4}`.
* *Gap:* "regime slabs are thick" is argued only for balanced sign columns. For actual or gap-admissible data, a
  short argument gives a thick column:
  * Pair separation and `sum e_p log p <= W` give some column with `sum_{x<y} |Delta_p| >= N_2 (1 - 4 ln C'/W)`.
  * Cauchy-Schwarz then gives `Var<k, beta_p> = Theta(n)` for a suitable base point.

  A rigorous lower bound on the tail mass at `w = O(sqrt(n ln r))` needs a large-deviation lower bound, not
  Berry-Esseen. Since nothing depends on it, the claim should be labelled *evidence/heuristic* unless supplied.

**§5.5.** The bullets are correct. `|<k, beta>| <= |k|_1` for `beta in [-1,1]^n`, and `|k|_1` is even, so
profile-independent subsets lie in `S_w`. The standard monomials of `S_w` have support `<= 2j`, so the interpolation
degree is `<= 4j <= 2w < m'`. The sentence "so every regime proof must use a non-product, profile-specific certificate
or subset" is a heuristic summary, not a theorem.

## 6. Prop 7.1 (literature vacuity) and literature checks

**Prop 7.1. Numbers confirmed exactly** (`check_slab_conc_ref.py`; trinomial big integers, my own code):

| `n`, `r` | `w` | `j` | `2j+1` | `m'` | `Pr(Q_S != 0)`, `s = floor(d/2)` | `s = d` | Hoeffding (`s = d/2`) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2901, `n` | 465.97 | 232 | 465 | 1864 | `10^{-13.18}` | `10^{-6.91}` | `10^{-7.83}` |
| 5000, `n` | 629.33 | 314 | 629 | 2518 | `10^{-13.92}` | `10^{-7.30}` | `10^{-8.32}` |
| 6000, `n^2` | 949.50 | 474 | 949 | 3798 | `10^{-25.70}` | `10^{-13.22}` | `10^{-16.04}` |

Remarks.
* (i) `n = 2901` violates `d >= 316 ln 4r` (`316 ln 11604 = 2957.5`). It is the record's *exact* `m_max >= 4w`
  threshold (2900 for `r = n`, round12 README §2.1), so it should be called that. `n = 5000` and `n = 6000` are in the
  regime.
* (ii) At these parameters the normalised influence bound `(2j+1)/s` is `0.16`-`0.32`, not small. It tends to 0 only
  asymptotically (`~ 2w/d`).
* (iii) **Scope of the "consequence".** It refutes exactly the bounds of the form `F(D, eps, tau)` (degree, variance,
  max influence).
  * Meka-Nguyen-Vu is a bound in the **rank**, of the form `rank^{-1/(4D+1)}` (memory, as labelled). It is vacuous
    here for a different reason: `rank^{1/(4D+1)} ≈ 1` at `D >= 465`.
  * More importantly, `Q_S = l_{1_S}^{2j+1} + (lower)`, and `u = 1_S` is a collision point, outside `G`. So Prop 7.1
    does **not** show that anti-concentration fails for the relevant class `{l_u^{m'} - P : u in G_1}`.
  * What it shows is that tools which cannot see the restriction `u in G` are vacuous. This is the correct and useful
    reading. The write-up's "necessarily vacuous at the needed degree" should be qualified to "for bounds depending
    only on degree, variance and influences".
* (iv) In the cube-threshold form, the needed small-ball bound is only `< mu(A)`, a probability near `1/2`, not small
  (§2). Prop 7.1 shows `F(D, 0, tau)` must be `>= 1 - 10^{-13}`. So such tools cannot even give `< 1/2`. This is
  consistent.

**Citations located by my own search:**
* Kwan-Sauermann, *Resolution of the quadratic Littlewood-Offord problem*: confirmed as Compositio Math. 161(12)
  (2025) 3089-3139.
* Jin-Kwan-Sauermann-Wang, arXiv 2505.23335, *Algebraic aspects of the polynomial Littlewood-Offord problem*:
  confirmed. Its abstract says the MNV bound can be substantially improved unless `f` is close to a polynomial with
  special algebraic structure.
* Bombieri-Vaaler: the statement shape is confirmed; the theorem number is not.
* O'Donnell Thm 10.21: the constant is not confirmed (see §5).
* The remaining items are labelled [memory] by the author. I did not contradict them.
* The author's characterisation of each tool's hypotheses (degree, product measure, Gaussian for Kane, bounded
  degree for inverse Littlewood-Offord) matches what I know. None applies at degree `Theta(sqrt(n ln n))` on `mu_A`.
  I agree with the net verdict of §7.

## 7. Fake-pattern-box checks

* **Proved statements.** All of §2-§4 hold for arbitrary finite `A`. §5's positive cases need twin classes (Cor 5.6)
  or `mu(A^c) <= 3^{-3m}/4` (Thm 5.8). Gap-admissible fakes satisfy the same pair inequalities and the same
  `sum e_p log p <= W`, so Prop 5.7 excludes twin classes for them. Fakes also have thick slabs (§5 gap). So no
  proved positive case applies to a gap-admissible fake. **No fake contradicts a proved statement.**
* **The fake mechanism against `AC_cube`.** This is a further reason the cube form is not equivalent; it is a
  heuristic, not a theorem. `AC_cube` at the *actual* box quantifies over every set `B` in the cube with
  `mu(B) >= mu(A)`. This includes fake pattern boxes and relabelled boxes. Suppose such a `B` has a nondegenerate real
  tangent-cone point (fake.md Thm 2.1) whose pairwise gaps exceed the actual radii `e^{-eps_xy/2}` after relabelling.
  Then `AC_cube` fails for the actual cluster, whether or not MIN^pair holds. (Typical actual radii are
  `e^{-abar/2} ~ n^{-3/4}` at `W ~ 3 n ln n`.) The §2 toy is exactly this, with `B` a 7-point set of the same mass.
  `AC_cube` should therefore be dropped from the list of equivalent forms.
* **`(R_H)`.** The author's §8 "fake test" is correct. `(R_H)` is a statement about `(J(A), V_{H'}, G)` alone. By fake
  Thm 2.1 it fails for fake data with a nondegenerate real point of `TC_{H'}` in `G` satisfying (2.2). So any proof of
  `(R_H)` must separate actual from fake boxes through `J(A)` and `V_{H'}`. I found no missed fake that falsifies a
  claim of the write-up.

## 8. Corrected statement

Setting: `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`, regime `d >= 316 ln 4r` (with `ln r = o(d)`
where `r_0 = e^{o(d)}` is used), `W > W_1`, `H' = e^{m'W/40}`, any profile, any norm.

1. **(Proved.) Sup and `L^2` forms.** For every finite `A` in the ternary cube, every `m` and every `u`:
   * `Cheb_A(u) = d^m kappa_V(u)`;
   * `Psi_A(u) <= Cheb_A(u)^2 <= U Psi_A(u)`;
   * `Psi_A = (1/U) sum_i Phi_{c_i}^2` is a sum-of-squares form with zero set `TC_V`.

   For the cluster box, `AC_2` (`min_{G_1} Psi_A > d^{2m'} eps_0^2`) **implies** `MIN^pair_inf` losslessly, and is
   implied by `MIN^pair_inf` at level `eps_0 3^{-d/2}`.
2. **(Proved, one direction only.) Probabilistic forms.** Every "most `k`" form (`AC_p`, fixed `p > 1/U`) and the cube
   form `AC_cube` imply `AC_2`. The converses are **false in general**: there are exact toys in which `MIN^pair_inf`
   holds and `AC_p` (`p > 1/7`) and `AC_cube` fail at every level. They are **open** for regime boxes. Only the
   threshold `p = 1/U` ("some `k`") is equivalent, up to `e^{O(d)}`. The relevant quantity is the `L^2(mu_A)` size of
   the top component, not anti-concentration. The `e^{O(d)}` losses are *not* shown to be negligible for the growth
   argument, because the Step-0 level `eps_0` is fixed.
3. **(Proved.) Graded form and restriction ratio.**
   * `Psi_A(u) = q_A(l_u^m mod J(A)_m)`, with `q_A` positive definite on `R[k']_m/J(A)_m`.
   * `AC_2 <=> rho_A(u) > mu(A) eps_0^2 d^{2m'}/Psi_cube(u)` on `G_1`.
   * A sufficient condition is `rho_A > e^{-W/2 + O(m'(ln n + abar_n))}`.
4. **(Proved, Bombieri-Vaaler quoted.) Codimension.** `codim(V_H in V(A,m)) <= r_0 ln U/(2 ln H)`. This is relatively
   `e^{-Omega(d)}` when `ln r = o(d)`, which holds in the growth range because `r <= W/ln 5`.
5. **The passage.** (P1) holds for the uncut cube with `K = e^{O(m ln d)}` (proof in §4 above). It is unproved for one
   cut and few cuts, where MIN^pair with `kappa_H` is proved directly. (P1) and (P2) are open in the regime.
6. **(Proved.) Model cases.**
   * Cube: exact `Psi`, and MIN^pair for `m <= theta d`.
   * Products: exact product formula.
   * One cut: MIN^pair for every `L`.
   * Few cuts: MIN^pair for `q <= (1/2) log_2(d/ln 4r) - 5.6` and `n >= n_0(C')`.
   * Thin cuts: `rho_A >= 1/2` for `mu(A^c) <= 3^{-3m}/4`. The constant depends on the quoted hypercontractivity; with
     the weaker form it is `120^{-4m}/4`.
   * None applies to actual or gap-admissible regime data (Prop 5.7). "Regime slabs are thick" is proved only for
     balanced sign columns.
7. **(Proved, exact numbers.) Prop 7.1.** Anti-concentration bounds depending only on degree, variance and influences
   are vacuous at degree `>= 2j + 1` on the cube. This does not refute anti-concentration for the class
   `{l_u^{m'} - P : u in G_1}`.
8. **Minimal open statement (corrected).** `(R_H)`: `min_{G_1} Psi_{A,H'} > d^{2m'} eps_0^2`. It implies MIN^pair
   losslessly, and MIN^pair implies it at level `eps_0 3^{-d/2}`.
   * It follows from `(R_inf^lambda) + (P2^lambda)` **with a matched margin `lambda`** (§4). It does not follow from
     `(R_inf) + (P2)` as stated in §8 of the write-up.
   * It is open for every regime box.
   * It must distinguish actual from fake pattern boxes through `J(A)` and `V_{H'}` alone.
9. **Growth/uniformity: no movement.** `M <= (2/3 + eps) log R/loglog R` stands.

## 9. Scripts

All are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_reduction_checks/`, and each writes `*_out.txt`.

| script | what | type |
| --- | --- | --- |
| `rlib.py` | my own exact `Psi^Sigma_A(u)` over `Q` (python-flint): independent monomial columns, Gram solve | library |
| `check_cube_product_ref.py` | Thm 5.1(a): 63 cases, 0 mismatches. Thm 5.2 with an irregular non-slab factor: 21 cases, 0 mismatches | **exact** |
| `check_prob_forms.py` | toy refutation of the probabilistic-equivalence claim (`A` one-cut, `d = 2`, `m = 3`, `u = (1/2, 1)`): `Psi_A = 3/8`, `Psi_B = 0` with `|B| = |A| = 7`; a 6-point subset of `A` with `Psi = 0`; `TC(A,3)` grid check; cube minus a point at `m = 4` | **exact** |
| `check_slab_conc_ref.py` | Prop 7.1 tail probabilities at `n = 2901, 5000` (`r = n`) and `n = 6000` (`r = n^2`); regime-condition check | **exact** (big integers), printed in log10 |

Working notes: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_reduction_notes.md`.
