# Conditional route: Vojta's conjecture implies the uniform `C sqrt(R)` bound

Route: conditional theorems. Candidates examined: (a) Bombieri–Lang / Caporaso–Harris–Mazur,
(b) Vojta's conjecture, (c) `abc` and the n-term `abc` conjecture. The three verdicts follow.

* **(b) works, and needs no extra uniformity.** One instance of Vojta's Main Conjecture
  implies the full uniform theorem. The instance is for rational points, with `D = 0`, on one
  fixed smooth projective rational 7-fold over `Q`. That 7-fold is a resolution of the blow-up
  of the toric variety `X_7 = {x_1^2+y_1^2 = ... = x_7^2+y_7^2}` in `P^13` along its diagonal
  line `Y`. The conclusion is in fact stronger than the target. There is an **absolute**
  constant `M_0` such that for every `C > 0` and every `R >= c_0^2 C^42`, every arc of length
  at most `C sqrt(R)` contains at most `M_0` lattice points. With more points, the same
  statement holds for arcs of length `C R^alpha`, for every `alpha < 1`. That is the full
  Cilleruelo–Granville `R^(1-eps)` conjecture.
  No uniformity in a finite place set `S` is needed, and no uniformity across a family. The
  exceptional set is the standard one, which depends only on the variety, `A` and `eps`.
  The earlier "claimed complete proof" needed exceptional sets uniform in `S`. Its route went
  through `S`-units with `S` equal to the primes of `N`. The present route never introduces `S`.
* **Unconditional by-products.**
  1. A Bezout lemma gives an exact equivalence. The uniform theorem holds if and only if, for
     each `C`, the short-arc `M`-tuples are not Zariski dense in `X_M` for some `M`.
  2. The same framework, fed with the **Ru–Vojta theorem** (a theorem based on the Subspace
     Theorem) with `S = {infinity}`, gives an ineffective re-proof of the qualitative
     Cilleruelo–Córdoba theorem: arcs of length `C R^alpha` with `alpha < 1/2` hold boundedly
     many points. The input is a rigorous lower bound `beta(H,E) >= 2n/(n+1)^(1+1/n) -> 2`.
  3. The framework stalls exactly at `alpha = 1/2`, because reaching it needs
     `beta(H,E) > 2` on some `X_M`. Numerically `beta` exceeds the generic value by only about
     1% for `M <= 8`, so `beta(X_7)` is about `1.25`.
* **(c) `abc` gives only a restricted theorem.** Under `abc`, every pair of lattice points at
  distance `<= C R^alpha` has a "difference norm" `n` with `rad(n) >= c_delta R^(1-alpha-delta)/C`.
  Hence a circle with `rad(R^2) <= R^(1-alpha-eta)` has at most one lattice point on any arc of
  length `C R^alpha`, once `R >= R_0(C, eta, alpha)`. A rigorous proposition shows that every
  `abc`-type inequality is vacuous on squarefree near-extremal block relations. This covers
  pairs, Ptolemy relations, n-term `abc` and `abc` over `Q(i)`. So `abc` cannot give a growth
  improvement through these relations.
* **(a) gives no implication.** Bombieri–Lang alone gives nothing that I could find. `X_M` and
  its blow-ups are rational, so `K` is not big. The arc condition is an archimedean
  approximation condition, which Lang's conjecture does not see. Vojta's conjecture implies
  Bombieri–Lang, so hypothesis (b) is stronger, but it is standard.

Novelty check: I grepped the notes for "Vojta", "blow", "Zariski dense", "canonical class",
"toric", "abc" and "rad(". The notes use Ru–Vojta only in the `S`-unit/GCD formulation on
`P^2`, where `S` is the primes dividing `N`. See `claimed_complete_proof_salvage_audit.md`,
`codex_quantitative_audit_result.md` and `codex_uniform_bound_research.md` §8. The notes
mention `abc` only as a non-applicability remark for polynomial families
(`polynomial_value_sharing_abc_applicability_gap.md`) and as an informal label in
`balanced_clique_residual_factorization.md`. None of Theorems 3.2, 5.1, 6.2 or 8.2 below
appears there. I did not do a thorough literature search. A bounded web search found no
published derivation of the `C sqrt(R)` arc bound from Bombieri–Lang or Vojta.

Everything labelled *Theorem*, *Lemma* or *Proposition* below is proved. Items labelled
*Numerics* or *Heuristic* are evidence only.

---------------------------------------------------------------------------------------------

## 1. Setting

For `M >= 2` let

```text
X_M = { [x_1:y_1:...:x_M:y_M] in P^(2M-1) : x_1^2+y_1^2 = x_j^2+y_j^2  (j = 2..M) },
Y   = { x_j = x_1, y_j = y_1  (j = 2..M) }  (a line, Y = P^1),   Y subset X_M,
H   = hyperplane class.
```

Write `z_j = x_j + i y_j`. Over `Q(i)` put `u_j = x_j + i y_j` and `v_j = x_j - i y_j`. Then
`X_M` becomes the split variety `{u_1 v_1 = ... = u_M v_M}`.

For `C > 0` and `alpha in [0,1)`, let `T^(M)_{C,alpha}` be the set of points
`[z_1 : ... : z_M]` in `X_M(Q)` with the following properties. The `z_j` are pairwise distinct
Gaussian integers of one common norm `N = R^2`, and `|z_i - z_j| <= C R^alpha` for all `i, j`.
Write `T^(M)_C = T^(M)_{C,1/2}`.

Every set `S` of lattice points on an arc of length `<= C R^alpha` of `x^2+y^2 = R^2`
(arbitrary position, `R^2` in `Z`, `R` not necessarily an integer) has pairwise chords
`<= C R^alpha`. So every ordered `M`-tuple of distinct points of `S` lies in
`T^(M)_{C,alpha}`.

Heights: `h(P) = log max|coord|` for a primitive integral representative. The proximity to
`Y` at the archimedean place is
`lambda_Y(P) = log ||P||_inf - log max_{j, coord} |coord(z_j - z_1)|`.
This is a Weil function for the subscheme `Y`. Since `Y` is a linear subspace contained in
`X_M`, its ideal in `X_M` is generated by the linear forms `x_j - x_1` and `y_j - y_1`.
Weil functions are well defined up to `O(1)`.

## 2. Geometry of `X_M`

**Lemma 2.1.** The following hold for every `M >= 2`.

1. `X_M` is a geometrically integral complete intersection of `M-1` quadrics, of dimension
   `M`, with `omega_{X_M} = O(-2)`.
2. Over `Q(i)`, `Sing(X_M)` is the set of points at which at least two of the pairs
   `(u_j, v_j)` vanish. It has codimension 3, so `X_M` is normal.
3. `X_M` is a normal projective toric variety over `Q-bar`. The torus is
   `{(lambda_j, mu_j) : lambda_j mu_j = lambda_1 mu_1}/G_m`, acting by
   `u_j -> lambda_j u_j` and `v_j -> mu_j v_j`. `X_M` is Gorenstein, hence has canonical
   singularities.
4. `Y` lies in the smooth locus. Every point of `T^(M)_{C,alpha}` lies in the open torus,
   because all `u_j v_j = N` are nonzero. In particular these points are smooth points of
   `X_M` and are not on `Y`.

*Proof.*

1. In the affine cone `Chat = {u_j v_j = u_1 v_1}` of `A^(2M)`, put `t = u_1 v_1`.
   - The part with `t != 0` is isomorphic to `G_m x G_m^M`, with coordinates
     `(t, u_1, ..., u_M)` and `v_j = t/u_j`. It is irreducible of dimension `M+1`.
   - The part with `t = 0` is a union of `2^M` linear spaces of dimension `M`.
   - Every component of a complete intersection of codimension `M-1` in `A^(2M)` has
     dimension at least `M+1`. So no component lies in `{t = 0}`, and `Chat` is irreducible.
   - It is reduced, because it is Cohen–Macaulay and generically smooth.
   - Adjunction for complete intersections gives `omega = O(-2M + 2(M-1)) = O(-2)`.
2. The Jacobian of `f_j = u_1 v_1 - u_j v_j` (`j >= 2`) has rows of the form
   `(v_1, u_1 | 0 ... | -v_j, -u_j | 0 ...)`.
   - If at most one pair vanishes, the rows are independent, because every row has a private
     nonzero block or the block `(v_1, u_1) != 0`.
   - If two pairs vanish, two rows become dependent.
   - On the locus where two pairs vanish, `t = 0`. So each of the other `M-2` pairs has one
     zero coordinate. This locus has affine dimension `M-2`, which is codimension 3 in `Chat`.
   - Normality follows from Serre's criterion: `S_2` holds by Cohen–Macaulayness, and `R_1`
     holds because the singular locus has codimension at least 2.
   - Exact check (A) in `check_vojta_geometry.py` confirms the ranks for `M = 3..7`.
3. Two points of the open set `{all u_j, v_j != 0}` differ by
   `lambda_j = u'_j/u_j, mu_j = v'_j/v_j`. These satisfy `lambda_j mu_j = t'/t` for all `j`.
   So the open set is one torus orbit.
   - A normal variety with a dense torus orbit is toric.
   - Normal Gorenstein toric singularities are canonical. Toric singularities are klt when
     `K` is Q-Cartier, and integral discrepancies `> -1` are `>= 0`; see Cox–Little–Schenck,
     *Toric Varieties*, §11.4.
   - Locally, near a point where exactly the pairs `1..k` vanish, `X_M` is
     `{u_1v_1 = ... = u_kv_k} x (smooth)`. This is because `v_j = t/u_j` eliminates the
     other pairs.
4. On `Y` every pair equals one nonzero pair `(u, v)`. So `Y` is in the smooth locus by
   item 2. `[]`

**The model.** Let `pi : Bl = Bl_Y X_M -> X_M`, with exceptional divisor `E`. `Y` is smooth of
codimension `M-1` inside the smooth locus. Therefore `Bl` is smooth near `E`, and

```text
K_Bl = pi^* K_{X_M} + (M-2) E = -2 pi^*H + (M-2) E.
```

Let `phi : X^#_M -> Bl` be a resolution of singularities defined over `Q` (Hironaka). Choose it
to be an isomorphism over the smooth locus of `Bl`, which contains `E`. Since `Bl` has the
canonical singularities of `X_M` away from `E`,

```text
K_{X^#} = phi^*( -2 pi^*H + (M-2) E ) + sum a_i F_i,   a_i >= 0,      (2.1)
```

where the `F_i` lie over `Sing(X_M)`. Write `H, E` again for their pull-backs to `X^#`.

## 3. Bezout reduction: Zariski non-density is equivalent to the uniform theorem

**Lemma 3.1.** Fix `N != 0` and let `Q_N` be the conic `{x^2 + y^2 = N}` in `C^2`
(or `{xy = N}`). Let `F` be a polynomial in `(x_1, y_1, ..., x_M, y_M)` of degree at most `D`
in each pair. Let `S` be a subset of `Q_N` with `|S| >= 2D + M`. If
`F(s_1, ..., s_M) = 0` for all pairwise distinct `s_1, ..., s_M` in `S`, then `F` vanishes on
`Q_N^M`.

*Proof.* Induct on `j` to show that `F(w_1, .., w_j, s_{j+1}, .., s_M) = 0` for all `w` in
`Q_N` and all pairwise distinct `s` in `S`.

* The case `j = 0` is the hypothesis.
* For the step, fix `w_1..w_j` and distinct `s_{j+2}..s_M`. The polynomial
  `G(w) = F(w_1, .., w_j, w, s_{j+2}, .., s_M)` has degree at most `D`.
* `G` vanishes on `S \ {s_{j+2}, .., s_M}`, which has at least `2D + 1` points of the
  irreducible conic `Q_N`.
* By Bezout, `G` vanishes on `Q_N`. `[]`

**Proposition 3.2.** Let `M >= 1`, `C > 0` and `alpha in [0,1)`.

* **(i)** Let `W` be a proper Zariski-closed subset of `X_M`. Let `D` be the degree of a form
  `F` that vanishes on `W` but not on `X_M`. Suppose all ordered `M`-tuples of distinct
  lattice points of a short arc lie in `W`. Then that arc has at most `2D + M - 1` lattice
  points.
* **(ii)** Consequently the following are equivalent:
  - the uniform theorem (for exponent `alpha`), and
  - for every `C` there is an `M` such that `T^(M)_{C,alpha}` is not Zariski dense in `X_M`.

*Proof.*

* **(i)** If the arc carried `2D + M` points, Lemma 3.1 would make `F` vanish on `Q_N^M`.
  `F` is homogeneous, so `F(lambda p) = lambda^D F(p)`. Hence `F` also vanishes on
  `Q_{lambda^2 N}^M` for every `lambda` in `C^*`. That is the dense set `{t != 0}` of the
  irreducible cone `Chat`. So `F` would vanish on `X_M`, a contradiction.
* **(ii)** The direction "non-density implies the uniform theorem" is (i), applied with `W`
  equal to the closure of `T^(M)_{C,alpha}`. Conversely, if every such arc has at most
  `M(C)` points, then `T^(M(C)+1)_{C,alpha}` is empty. `[]`

Remark. Caporaso–Harris–Mazur pass from non-density on fibred powers to uniform bounds. Here
the family of circles `x^2 + y^2 = N` is isotrivial: all members are related by scaling. So
every fibred power is the single fixed variety `X_M`, and Lemma 3.1 replaces the CHM machinery.
This is why no uniformity in families is needed below.

## 4. Height and proximity of short-arc tuples

**Lemma 4.1.** Let `z_1, .., z_M` be pairwise distinct Gaussian integers with `|z_j| = R`
and `|z_i - z_j| <= C R^alpha`. Let `d` be the gcd of the `2M` coordinates, `R' = R/d`, and
`P = [z_1 : ... : z_M]`. Then:

* (a) `R' >= R^(1-alpha)/C`;
* (b) `|h(P) - log R'| <= (1/2) log 2`;
* (c) `lambda_Y(P) >= (1-alpha) log R' - log C - (1/2) log 2`.

*Proof.*

* (a) Two distinct points satisfy `z_i - z_j = d w` with `w` a nonzero Gaussian integer. So
  `d <= |z_i - z_j| <= C R^alpha`.
* (b) The primitive coordinates satisfy `|z_j/d| = R'`. So `max(|x|, |y|)` lies in
  `[R'/sqrt 2, R']`.
* (c) `|coord(z_j/d - z_1/d)| <= C R^alpha / d = C R'^alpha d^(alpha-1) <= C R'^alpha`, and
  `||P|| >= R'/sqrt 2`. `[]`

*Exact check (B).* The 7-subtuples of the eight-point translated Pell family were checked for
`n = 61, 121, ..., 361`, with `log R` up to 4181. The worst slack of (c) is `0.347`, which is
`(1/2) log 2` as expected. The quantity `(5 lambda_Y - 2h)/h` (the lower bound for
`h_K/h` used below) rises from `0.377` to `0.479`, towards `1/2`. See
`log_vojta_geometry.txt`.

## 5. The conditional theorem

**Hypothesis `V(M, eps)`.** This is Vojta's Main Conjecture (Vojta, LNM 1239, Conj. 3.4.3;
Bombieri–Gubler, Ch. 14) for the smooth projective variety `X^#_M` over `Q`, with `D = 0`,
the big divisor `A = H` (the pull-back of the hyperplane class) and the given `eps`:

> There are a proper Zariski-closed `Z` in `X^#_M` and a constant `c` such that
> `h_{K_{X^#}}(P) <= eps h_A(P) + c` for every `P` in `X^#_M(Q) \ Z`.

If one prefers `A` ample: an ample `A` satisfies `h_A <= m h_H + O(1)` off a proper closed set,
because `m H - A` is big for large `m`. Add that set to `Z` and use `eps/m`. Since `D = 0`, no
set of places `S` enters.

**Theorem 5.1.** Assume `V(7, 1/4)`. Then there are an absolute integer `M_0` and an
(ineffective) constant `c_0 >= 1` with the following property. For every `C > 0`, every real
`R` with `R^2` in `Z` and `R > c_0^2 C^42`, and every arc of length `<= C sqrt(R)` on
`x^2 + y^2 = R^2` in any position, the arc contains at most `M_0` lattice points. For
`R <= c_0^2 C^42`, the trivial bound `<= C sqrt(R) + 1` applies. So
`M(C) <= max(M_0, c_0 C^22 + 1) < infinity` for every `C`.

*Proof.* Let `P` be in `T^(7)_C`, and let `P^#` be its unique lift to `X^#_7`. The lift is
unique because `P` is a smooth point off `Y` (Lemma 2.1, item 4). By (2.1) and functoriality
of Weil heights,

```text
h_{K}(P^#) = -2 h(P) + 5 h_E(P^#) + sum a_i h_{F_i}(P^#) + O(1).
```

Each `F_i` and `E` is effective and does not contain `P^#`. So `h_{F_i}(P^#) >= -O(1)`. The
local Weil functions of `E` are bounded below at every place and are `>= 0` at almost all
places, so `h_E(P^#) >= lambda_{E,inf}(P^#) - O(1)`. Also
`lambda_{E,inf} = lambda_{Y,inf} o pi + O(1)`, because `pi^{-1} I_Y O_Bl = O(-E)` (Silverman,
Math. Ann. 279 (1987), on arithmetic distance functions). With Lemma 4.1 (`alpha = 1/2`):

```text
h_K(P^#) >= -2 log R' + 5( (1/2) log R' - log C ) - O(1) = (1/2) log R' - 5 log C - O(1).
```

If `P^#` is not in `Z`, then `V(7, 1/4)` gives `h_K <= (1/4) log R' + O(1)`. Hence
`R' <= c_0 C^20`, where `c_0` depends only on `X^#`, the chosen Weil functions and `c`.

Now take any arc with `R > c_0^2 C^42` and any 7 distinct lattice points on it. By
Lemma 4.1(a), `R' >= sqrt(R)/C > c_0 C^20`. So the 7-tuple lies in `W = phi pi (Z)`, a proper
closed subset of `X_7` that does not depend on `C`. Proposition 3.2(i) then gives at most
`M_0 := 2 deg F_W + 6` points, where `F_W` is a fixed form vanishing on `W` but not on `X_7`.

For small `R`, consecutive lattice points on an arc are at arc distance at least 1, so an arc
of length `L` has at most `L + 1` points. `[]`

**Theorem 5.2 (general exponent).** Let `alpha` be in `[0,1)` and put
`M_alpha = 3 + floor(2/(1-alpha))`. This gives `M_{1/2} = 7`, and `M_alpha -> infinity` as
`alpha -> 1`. Put `eps_alpha = ((M_alpha - 2)(1-alpha) - 2)/2 > 0`. Then `V(M_alpha, eps_alpha)`
implies the following. There is an absolute `M_0(alpha)` such that for every `C` and
`R >= R_1(alpha, C)`, every arc of length `<= C R^alpha` has at most `M_0(alpha)` lattice
points.

In particular, Vojta's conjecture implies the Cilleruelo–Granville conjecture for arcs
`R^(1-eps)` (and for `C R^(1-eps)`). The proof is the same: `h_K >= ((M-2)(1-alpha) - 2) log R'
- (M-2) log C - O(1)`, and `R' >= R^(1-alpha)/C` by Lemma 4.1(a).

**What uniformity is used.** Only the standard dependence of Vojta's exceptional set `Z` on
`(X^#, A, eps)`.

* No finite set of places `S` appears, because `D = 0`. The Gaussian arithmetic of `N` (its
  primes, gcds and valuation profiles) never enters. It is all absorbed into the fact that the
  tuple is a rational point of height `log R'` on one fixed variety.
* The earlier route (`claimed_complete_proof_salvage_audit.md`) put `u_j = z_j/z_0` as
  `S`-units on `P^2`, with `S` equal to the primes of `N`. That route needed an exceptional set
  independent of `S`, which is not part of any standard conjecture. Its geometric gain was a
  threshold of 3 points instead of 7.
* The constants `c_0` and `deg F_W` are ineffective, as Vojta's conjecture is.

**Corollary 5.3 (hyperbola, same proof).** Assume `V(7, 1/4)` for the split model
`X^split_7 = {x_1y_1 = ... = x_7y_7}`. Lemma 2.1 holds verbatim with `u = x`, `v = y`, and
Lemma 3.1 holds for the conics `xy = n`. Then the Erdős–Rosenfeld question has a positive
answer: there is an absolute `K` such that, for every `c > 0` and every `n >= n_0(c)`, `n` has
at most `K` divisors in `[sqrt(n) - c n^(1/4), sqrt(n) + c n^(1/4)]`.

Also, `V(3 + floor(1/eps), .)` implies Ruzsa's conjecture: at most `K_eps` divisors in
`sqrt(n) +- n^(1/2-eps)`.

*Proof sketch.* Consider the points `(d, n/d)`. Pairwise distances are
`<= 3c n^(1/4) = O(H^(1/2))` with `H = sqrt(n)(1 + o(1))`. The content `g` of a 7-tuple
satisfies `g <= 3c n^(1/4)`, so `H/g -> infinity`. Then argue as in Theorem 5.1.

I have not checked whether this corollary appears in the literature. The problem status is
taken from Chan, arXiv:1406.2230, and the Erdős-problems forum, via search snippets only.

### 5.1 Consistency checks: the known infinite families lie on rational quartics through Y

Under `V(7, .)` the infinitely many known short-arc tuples must lie in `Z`. They do lie on
curves.

* **Pell family (exact checks (C)).** In the eight-point translated Pell family, every point is
  `K_p` times a product of four `Z[i]`-linear forms in `(x_n, y_n)`. So all eight points are
  the values at `(x_n : y_n)` of one map `Gamma : P^1 -> X_8` of degree 4. The check reproduces
  the `n = 61` points.
  - At `(±sqrt5 : 1)` all eight words coincide, so `Gamma` meets `Y` there. The identity behind
    this is `((1 + sqrt5) + 2i)^2 = (2 + 2 sqrt5)(1 + 2i)`, which makes `F G-bar^2` real.
  - The Pell solutions are Roth-optimal approximations to `sqrt5`, with exponent 2 on `P^1`.
    Dividing by the degree 4 gives `lambda_Y/h -> 2/4 = 1/2`.
* **Cilleruelo–Granville family (exact check (D)).** The six-point family
  `prod_{j=1}^4 (a + j + i sigma_j)` (with `sum sigma = 0`) is a rational quartic.
  - The coefficients of `a^4` and `a^3` agree across all six points, and the coefficient of
    `a^2` does not. So the curve has contact order 2 with `Y` at the rational point
    `a = infinity`.
  - The exponent is `2 x 1/4 = 1/2`. Measured `lambda_Y/h` rises from `0.436` to `0.479`.

So the exponent `1/2` is attained along rational curves exactly as
(contact × Roth exponent)/degree. This matches McKinnon's heuristic that best approximations
lie on rational curves.

*Heuristic, not used.* `X_7` has about `H^2` rational points of height `<= H`. The tube of
radius `C H^(-1/2)` around `Y` has real codimension 6. So the expected number of "sporadic"
7-tuples in that tube is `H^(2-3)`, which is summable. Hypothesis `V(7, .)` would be false if
`T^(7)_C` were Zariski dense for some `C`, for example through infinitely many unit templates
with uniformly bounded arc constant. I know of no such family.

## 6. Why the route is not unconditional: the `beta` criterion

The strongest unconditional theorem of the same shape is the Ru–Vojta theorem
(Amer. J. Math. 142 (2020)). It is a consequence of the Subspace Theorem: for an effective
Cartier divisor `D` and a big `L`,
`beta(L, D) m_{D,S}(P) <= (1 + eps) h_L(P)` off a proper closed set, with
`beta(L, D) = liminf_N sum_{m >= 1} h^0(NL - mD) / (N h^0(NL))`.

**Proposition 6.1.** Let `beta_M = beta(H, E)` on `X^#_M`. If `(1 - alpha) beta_M > 1`, then
(unconditionally and ineffectively) arcs of length `C R^alpha` carry at most `M_0` lattice
points for `R >= R_1(C)`.

*Proof.* Apply Ru–Vojta with `S = {infinity}`, the single divisor `E` and `L = H`. Use
Lemma 4.1(c) and Proposition 3.2 exactly as in Theorem 5.1. The set `S` is fixed, so the
exceptional set is fixed. `[]`

**Proposition 6.2 (computation of `beta_M`, rigorous lower bound).** Put `n = M - 1`. Then

```text
beta_M  >=  2n / (n+1)^(1+1/n)   ->  2  (n -> infinity).
```

The values of this bound are 0.5, 0.770, 0.945, 1.070, 1.165, 1.239, 1.300 for
`M = 2..8`.

*Proof.*

1. **Basis.** `X_M` is a complete intersection, so it is projectively normal.
   `H^0(O(N))` has the monomial basis `t^k prod u_j^(e_j+) v_j^(e_j-)` with
   `2k + |e|_1 = N`, where `e` ranges over `Z^M`.
2. **Restriction near `Y`.** Near `Y` write `u_j = u w_j`, `v_j = v/w_j` with `w_1 = 1`. The
   order of vanishing along `Y` of a section is the minimum, over the weight
   `s = sum e_j` (with `s = N` mod 2), of the order at `w = 1` of a Laurent polynomial
   supported on
   `A_s = {e' in Z^n : |s - sum e'| + |e'|_1 <= N}`.
3. **Count.** Order `>= m` means the coefficients annihilate all polynomials of degree `< m`
   on `A_s`. Hence `h^0(NH - mE) = sum_s (|A_s| - h_{A_s}(m-1))`, where `h_A(d)` is the
   dimension of the space of polynomials of degree `<= d` restricted to `A`.
4. **Generic bound.** `h_A(d) <= min(|A|, binom(d+n, n))`. Therefore
   `sum_d (|A| - h_A(d)) >= (n/(n+1)) (n!|A|)^(1/n) |A| (1 - o(1))`.
5. **Volumes.** `|A_s| ~ vol(K_sigma) N^n`, where `sigma = s/N`, and
   `int_{-1}^{1} vol(K_sigma) dsigma = vol(B_1^(n+1)) = 2^(n+1)/(n+1)!`. The power-mean
   inequality `E[v^(1+1/n)] >= E[v]^(1+1/n)` then gives the bound. `[]`

**Corollary 6.3.** For every `alpha < 1/2`, choose `M` with `2n/(n+1)^(1+1/n) > 1/(1-alpha)`.
Proposition 6.1 then gives bounded counts on arcs `C R^alpha`. This is an ineffective re-proof
of the qualitative Cilleruelo–Córdoba theorem, which is a check that the framework is sound.

*Numerics (evidence, not proof).* Exact Hilbert functions of the `A_s` were computed mod two
primes (`beta_vanishing.py`, `beta_ratio.py`).

| `M` | largest `N` computed | `beta_N` | generic value at same `N` | ratio |
|---|---|---|---|---|
| 2 | 16 | 0.5000 | 0.5000 | 1.000 |
| 3 | 24 | 0.7725 | 0.7646 | 1.010 |
| 4 | 12 | 0.9275 | 0.9165 | 1.012 |
| 5 | 6 | 0.9869 | 0.9756 | 1.012 |
| 6 | 5 | 1.0252 | 1.0159 | 1.009 |
| 7 | 4 | 1.0368 | 1.0224 | 1.014 |
| 8 | 3 | 1.0208 | 1.0208 | 1.000 |

So the true `beta_M` is about 1% above the generic lower bound, and `beta_7` is about `1.25`.

**Where it stalls.** The comparison is between two thresholds.

* Vojta's canonical class gives non-density for proximity exponent `> 2/(M-2)`, which is
  `< 1/2` once `M >= 7`.
* Ru–Vojta gives non-density only for exponent `> 1/beta_M`, which is about `1/2 + Theta(log M / M)`.

The arc problem sits exactly at exponent `1/2`. This matches the notes' recurring finding that
"everything is critical at 1/2". The notes' unweighted optimum `(m+1)/(2m)` has the same shape.

**Open sub-question.** Is `beta(H, E) > 2` on `Bl_Y X_M` for some `M`? If so, Proposition 6.1
would prove the uniform theorem unconditionally, and even boundedness on arcs `R^(1/2+delta)`.

* The question is equivalent to an averaged volume computation for blow-ups of toric varieties
  at a general point of the torus. This is a Nagata-type problem.
* The generic bound tends to 2 only from below.
* The measured lattice excess (about 1%) does not show growth with `M` up to 8, but `M` would
  have to be in the hundreds before the generic bound exceeds `1.98`.

## 7. Bombieri–Lang / CHM (candidate a)

* `X_M` is rational. It is birational over `Q` to `(P^1)^M`, with coordinates the direction of
  `z_1` and the half-angle parameters of the unit-circle points `z_j/z_1`. (The separate map
  `z_j -> (x_j : y_j)` to `(P^1)^M` has degree `2^(M-1)`.)
* `K_{X^#} = -2H + (M-2)E + (effective)` is not big. A general complete-intersection curve
  `Gamma` avoids `E` and the `F_i`, which lie over sets of codimension at least 2. Such
  curves form a covering family with `K . Gamma = -2 H . Gamma < 0`, so `K` is not
  pseudo-effective. So Lang's conjecture says nothing about `X^#`.
* The arc condition is an archimedean approximation condition. Bombieri–Lang concerns rational
  points without such conditions. Turning approximation into rational points on a cover of
  general type introduces twists that vary with the point, which needs uniformity across
  twists. I found no standard statement that supplies this.
* In the balanced extraction regime, the arc condition does become exact equations
  `Im(eps_ij U_ij(H)) = b_ij` with bounded `b`. That is an integral-points problem, and the
  relevant conjecture is Lang–Vojta for integral points, not Bombieri–Lang.
  - Those varieties contain rational curves (the polynomial families of §5.1 and item 533).
  - They contain twisted tori (the Pell templates).
  - Unbalanced profiles give no exact equation.
* Since Vojta implies Bombieri–Lang, Theorem 5.1 is the natural conditional statement. No
  conditional theorem from Bombieri–Lang alone is claimed.

## 8. `abc` (candidate c)

`abc`: for every `delta > 0` there is `K_delta` with `c <= K_delta rad(abc)^(1+delta)` for
coprime positive integers with `a + b = c`.

**Proposition 8.1 (pair identity; `abc` only for the last line).**

*Setup.* Let `z != w` be Gaussian integers with `|z| = |w| = R`, let `g = gcd(z, w)`,
`z = g u` and `n = |u|^2 = R^2/|g|^2`.

*Identity (unconditional).*

* `w = g eta u-bar` for a unit `eta`, `u` is conjugate-primitive, and `u = A + Bi` with
  `gcd(A, B) = 1`.
* Put `(s, t) = (B, A)`, `(A, B)`, `(A - B, A + B)` or `(A + B, A - B)` according as
  `eta = 1, -1, i, -i`.
* Then `s^2 + t^2 = n` or `2n`, `gcd(s, t) = 1`, `s != 0`, and
  `|z - w| = |g| |u - eta u-bar| >= |g| sqrt2 |s|`.

*Under `abc`.* If moreover `|z - w| <= C R^alpha`, then for every `delta > 0`

```text
rad(n) >= K_delta^(-1/(1+delta)) n^(-delta/(1+delta)) R^(1-alpha) / (2 sqrt2 C)
       >= c_delta C^(-1) R^(1-alpha-2 delta).
```

*Proof.* Apply `abc` to `s^2 + t^2 = c_0`, where `c_0` is `n` or `2n`:
`n <= c_0 <= K_delta (2|st| rad n)^(1+delta)`, with `|t| <= sqrt(2n)` and
`|s| <= C R^(alpha-1) sqrt(n)`. `[]`

The identity was verified exactly on 1582 close pairs from 60 random circles (`check_abc_pairs.py`, (E)).

**Corollary 8.2 (restricted theorem under `abc`).** For every `alpha < 1`, `eta > 0` and
`C > 0` there is `R_0` with the following property. If `R >= R_0` and
`rad(R^2) <= R^(1-alpha-eta)`, then every arc of length `<= C R^alpha` on `x^2 + y^2 = R^2`
contains at most one lattice point. For `alpha = 1/2` the condition is
`rad(R^2) <= R^(1/2-eta)`.

*Proof.* `n` divides `R^2`, so `rad(n) <= R^(1-alpha-eta)`. Combine with Proposition 8.1 for
`delta = eta/4`. `[]`

This is a statement about circles with a highly powerful norm, which is a different class from
the notes' unconditional squarefull result (constant `3 + eps`, `inverse_valuation_profile.md`).

The corollary is purely asymptotic, because the `abc` constants are uncontrolled.

* Machin's identity `239^2 + 1 = 2 * 13^4` gives the pair `28561` and `28560 - 239i` on
  `x^2 + y^2 = 13^8`, at distance `< 1.42 sqrt(R)`, with `rad(n) = 13 = R^(1/4)`. Its `abc`
  quality is `1.254`.
* Circles `5^a 13^b <= 10^12` still carry close pairs with `log rad(n)/log R` as low as
  `0.21` (check (E')).

**Proposition 8.3 (`abc`-type inequalities are vacuous on squarefree blocks).** Let
`x_i = n_i y_i` (`i = 1..k`) satisfy:

* `sum x_i = 0` and `gcd(x_i) = 1`;
* the `n_i` are pairwise coprime and squarefree;
* `|y_{i*}| <= prod_{j != i*} n_j`, where `i*` is an index maximising `|x_i|`.

Then `max|x_i| <= rad(prod x_i)`. Consequently every inequality `max|x_i| <= K rad^kappa`
with `K, kappa >= 1` holds automatically. This covers `abc` (any `delta`), the n-term `abc` of
Browkin–Brzeziński (exponent `2k - 5 + delta`), and their `Q(i)` versions.

*Proof.* `rad(prod x_i) >= prod n_i >= n_{i*}|y_{i*}|`.

Over `Q(i)`, a split block `pi^1` contributes `N(pi) = p` to the radical norm, which is the
same as over `Q`. `[]`

*Application.* The Ptolemy relation `n_A a + n_B b = n_C c` and its five-point versions
(two_thirds_sketch.md) have block norms `e^(w(1+o(1)))` and cofactors `e^(o(w))` near the
extremal profile. So the hypothesis holds whenever the blocks are squarefree. The notes'
inverse results put the near-extremal weight on exponent-one primes, which is exactly where
Proposition 8.3 applies. High prime powers are the only regime in which `abc` bites
(Corollary 8.2), and Gaussian structure does not help.

*Exact data (F).* The 15 Ptolemy quadruples of the six-point family were checked at
`a = 20, 50, 100, 300`. The relation `a + b = c` holds exactly, with
`log rad(abc)/log max >= 2.2` and maximal `abc` quality `0.45`. On 1582 random close pairs, the
maximal quality of `s^2 + t^2 = c_0` is `0.998`.

## 9. Missing lemma (unconditional) and why the method stalls

**Missing lemma (exactly equivalent to the target, by Proposition 3.2).** For each `C`, the
set `T^(M)_C` of `M`-tuples of distinct lattice points of one origin-centred circle with
pairwise distances `<= C sqrt(R)` is not Zariski dense in the fixed toric variety `X_M`, for
some `M`. A fixed `M = 7` with `W` independent of `C` gives the absolute form.

Equivalently, rational points of `X_7` approximating the line `Y` to exponent
`>= 1/2 - o(1)` are degenerate. Vojta's conjecture predicts this for exponents above `2/5`.

**Why it stalls.**

* The available unconditional theorems on a fixed variety are Subspace, Evertse–Ferretti and
  Ru–Vojta with `S = {infinity}`. They use only the section ring, through `beta(H, E)`. That
  constant is provably at least `2 - O(log M/M)` and numerically only about 1% higher, so it
  reaches only exponents `> 1/2`.
* Using the `S`-unit structure of `N` lowers the threshold, but it makes the exceptional sets
  depend on `S`. That is the failure recorded in the salvage audit.
* Closing the gap requires either a canonical-class input (Vojta), or a proof that
  `beta(H, E) > 2` on some `Bl_Y X_M`, or a genuinely arithmetic input that the
  section-ring methods do not see.

## Files (all under `scratchpad/growth/conditional_checks/`)

All scripts run with `python3` from that folder, and each writes the `log_*.txt` file listed
with it.

* `check_vojta_geometry.py`: exact checks (A) to (D), with log `log_vojta_geometry.txt`.
  - (A) singular locus and the smoothness of `Y`;
  - (B) Lemma 4.1 on the Pell 7-subtuples;
  - (C) the Pell family as a quartic through `Y` at `(±sqrt5 : 1)`;
  - (D) the Cilleruelo–Granville family and its contact order.
* `check_abc_pairs.py`: (E) the pair identity, radicals and qualities; (E') powerful circles;
  (F) Ptolemy and `abc`. Log: `log_abc_pairs.txt`.
* `beta_vanishing.py`, `beta_generic.py`, `beta_ratio.py`: the vanishing-order numerics of §6.
  Logs: `log_beta_generic.txt`, and `log_beta_ratio.txt` (pairs `n,N` = `2,24 ... 7,3`).
  The `M = 2` row of the table was computed separately with `beta_vanishing.py 1 4 8 16`.
* `gauss.py`: helper functions for Gaussian integers.
