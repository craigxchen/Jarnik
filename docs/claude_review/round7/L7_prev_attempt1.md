# (L_7), round 7: the cub13 family is exact, Q-density of its members is not the missing step, and an addability criterion

Task: decide `(L_7)` of `threshold/L8.md` §1 and `growth/ptolemy.md` §7, starting from the round-4
material (`round4/L7.md`, `round4/referee_L7.md`). The round-4 referee found NUMERICALLY a covering
3-dimensional family of twisted cubics in the Kapranov model of `Mbar_{0,7}` meeting 13 of the 20
Kapranov planes (`cub13`), and located the remaining gap in the Zariski density of Q-rational members.
Scripts are in `L7_checks/r7/` (Section 8).

Conventions (as in L8 and round 4).

* *Theorem*, *Proposition*, *Lemma* and *Corollary* mean proved here, with all quantifiers, or cited
  with the exact place where the proof is written.
* *Data* means an exact finite computation (Python integers, `Fraction`, exact arithmetic in a number
  field). *NUMERICAL* means floating point and is never used as a proof.
* `X = Mbar_{0,7}`, `beta` the balanced class (`beta.D_S = 1` for all 56 boundary divisors),
  `pair(g)`, `triple(g)` the sums of `D_S.g` over the 21 splits of type `2|5`, resp. 35 of type `3|4`.
  `-K.g = pair(g)/3` and `(B+8K).g = triple(g) - (5/3) pair(g)` (`round4/L7.md` §4).
* A class `g` is *rich* if `(B+8K).g >= 0`, i.e. `triple(g) >= (5/3) pair(g)`. `beta` is rich with
  equality.

---------------------------------------------------------------------------------------------

## 0. Verdict

**`(L_7)` is still not decided.** The round produced one exact existence theorem, one positive
criterion, and a correction of the framing of the gap.

1. **Theorem 1 (exact; upgrades the round-4 referee's numerics).** The cub13 family exists. It has an
   explicit parametrisation by four contact positions `(r25, r34, r35, w)` subject to one equation
   `Res = 0` (Section 2). A member defined over a quintic field `K` is certified by exact arithmetic
   in `K`: it is a rational curve of class `cub13`, meeting 17 boundary divisors in 19 distinct
   points, each transversally and on exactly one divisor. The parameter point is a smooth point of
   `{Res = 0}`, and the differential of `P^1 x W -> X` there has rank 4. So the family is
   3-dimensional and dominates `X`.
   * **Corollary 2.** Over `C`, free rational curves of class `q beta` (for example
     `q = 131040`) cover `X`. This is `referee_L7.md` §3.3, now unconditional. So the geometric
     (function-field) analogue of `(L_7)` is false. The component of the space of stable maps that
     contains these curves has smooth **Q-points** (Galois-stable combs, Proposition 2.4).
2. **Correction of the framing (Section 4).** Q-rational members of the *cub13* family, even a
   Zariski-dense set of them, do not by any known argument give Q-rational curves of class
   `~ q beta`. What is missing is a **smoothing over Q**. The classes reachable from a fixed seed by
   every Q-rational gluing operation found here lie asymptotically in the cone of *addable* classes
   (Proposition 4.5). The seed, cub13 or otherwise, drops out.
3. **Theorem 3 (addability criterion, new).** Call a Q-family `V` of rational curves *addable* if
   its universal curve `U_V -> X` has a rational section over `Q`. Examples are the fibres of any
   Q-rational map `X ⇢ Z^3` with rational generic fibre. If some rich family is addable, then
   `(L_7)` is false, under explicit genericity hypotheses (G1)–(G3), which are open conditions. The
   proof iterates sections of ruled surfaces `C x P^1`, then applies Theorem 2(a) of round 4.
4. **What is and is not addable (Proposition 4.4, Data and NUMERICAL).**
   * Every addable class known here has `(B+8K).g <= -(7/9) B.g`: the forgetful fibres, the plane
     projections in the 7 Kapranov models, and the toric 1-parameter subgroups. So their cone misses
     `beta`.
   * A 3-dimensional rich family (such as cub13) is addable only if it has exactly one member
     through a general point. cub13 has at least 8 (NUMERICAL, referee), so it is not addable.
   * The Q-rational 4-dimensional family `cub12` is explicit (Proposition 3.1), but it is not rich:
     `(B+8K) = -3`.
5. **Arithmetic of the cub13 family (Data, NUMERICAL).**
   * `W` is birational to the hypersurface `G = 0` in `A^4`, with `G` of degrees `(7,5,5,5)`.
   * Every natural projection of `W` to `A^3` has degree `>= 4`, and so does the contact map to
     `E_126`.
   * An exhaustive search over the `45^4 = 4,100,625` parameter points of height `<= 6` found 44 rational
     zeros of the residual. Every one is degenerate: it lies on a spurious factor, where two contacts
     collide. So no Q-rational cub13 curve is known.
6. **Barriers (Theorem 5).**
   * (a) Product charts built from a group law (torus or vector group) have hidden-contact-free
     rulings, so `beta` is never in their ruling cone.
   * (b) The L6 common-argument (Galois–Boolean) construction with four moving points has genus
     `>= 5`, so it has finitely many Q-points by Faltings. This is already in `threshold/L8.md` §5
     (item 112). Added here: balance forces the four square classes to be independent.
   * (c) Single-point lifting by interpolation needs `25q` linear conditions against `16q + 1`
     parameters. This restates round 4 §8.

**Net effect.**
* Deciding `(L_7)` is now a purely arithmetic question about Q-points. A proof of `(L_7)` would at
  least have to explain why the Q-points of a geometrically irreducible variety *that has smooth
  Q-points* (a component of the space of approximately balanced rational curves) never produce
  balanced configurations off a fixed `Z`. No mechanism for that is known.
* For a refutation, the problem is now the existence of **one rich addable family**. Concrete
  candidates, with explicit equations, are in Section 6: rational quartics through all 10 coordinate
  planes and 7 hidden planes, `(B+8K) = +2`, dimension 4.
* Nothing here changes the status of the uniform theorem or of growth (Section 7).

---------------------------------------------------------------------------------------------

## 1. Setting: what a refutation must produce

`(L_7)`, its negation, the anchored configurations and the dictionary are as in `threshold/L8.md`
§1 and `round4/L7.md` §1. Three refereed facts are used.

* **Theorem 2(a) of round 4** (refereed: correct). Let `F` be a set of morphisms
  `f: P^1_Q -> X` such that each `f`:
  * is birational onto an image not contained in the boundary;
  * maps no point into two distinct boundary divisors;
  * satisfies `|deg f^*D_S - q_f| <= eps_f q_f` for all 56 `S`.

  If for every `eps > 0` the union of the images of `{f : eps_f <= eps}` is Zariski dense, then
  `(L_7)` is false.
* **Round-4 referee §2.2.** A dominant product chart `(P^1)^R ⇢ X` whose ruling cone contains `beta`
  produces such an `F`. So Theorem 2(b) is a special case of 2(a).
* **Residues** never need to be controlled (round 4, Theorem 1 and Corollary 1). Admissibility is
  automatic at the good primes of a curve (proof of 2(a)).

So **a refutation of `(L_7)` is exactly a supply of Q-rational, approximately balanced rational
curves whose union is dense.** The curves must be defined over `Q` and isomorphic to `P^1_Q`. Curves
over number fields do not help: a curve not defined over `Q` has only finitely many Q-points.

**Notation for the Kapranov model used throughout.**
* `P^4` with `p_1, ..., p_5 = e_0, ..., e_4` and `p_6 = (1,1,1,1,1)`.
* The point `(c_1 : ... : c_5)` of the open torus is the 7-point configuration
  `(y_1, ..., y_7) = (c_1, ..., c_5, 0, ∞)`, with label 6 at `0` and label 7 at `∞`. This is the
  Kapranov `psi_7`-model:
  * the coordinate hyperplane `{c_j = 0}` is `D_{j6}`;
  * the exceptional divisor over `p_j` is `D_{j7}`;
  * `{c_i = c_j}` is `D_{ij}`;
  * `E_{p_6} = D_{67}`.
* A curve `c(t) = (c_1(t), ..., c_5(t))` has a contact with `D_T` (`T ⊆ {1..6}`, `2 <= |T| <= 5`)
  at `t_0` exactly when `T` is a level set of `i -> c_i(t_0)`, where `c_6 := 0`. The contact is at a
  codimension-2 stratum when two level sets have size `>= 2`. This is how all contacts are computed
  below.
* The Kapranov dictionary is used only for *naming*. Theorem 1 is a statement about a curve
  `t -> (c_1(t), ..., c_5(t), 0, ∞)` in `Mbar_{0,7}` and its contact counts with the 56 `D_T`. The
  resulting `D_S`-vector is the one of `referee_L7.md` §3.1, under `D_{abc7} = D_{[6] \ abc}`.
* The 13 planes of cub13 (`referee_L7.md` §3.1) correspond to the following level sets `T`.

  | plane | `T` | plane | `T` | plane | `T` |
  |---|---|---|---|---|---|
  | `123` | `{4,5,6}` | `134` | `{2,5,6}` | `246` | `{1,3,5}` |
  | `124` | `{3,5,6}` | `135` | `{2,4,6}` | `256` | `{1,3,4}` |
  | `125` | `{3,4,6}` | `145` | `{2,3,6}` | `346` | `{1,2,5}` |
  | `126` | `{3,4,5}` | `136` | `{2,4,5}` | `356` | `{1,2,4}` |
  | | | | | `456` | `{1,2,3}` |

  The remaining contacts of the class are pair contacts: `D_{16}` three times, and `D_{14}`,
  `D_{15}`, `D_{23}` once each, for 19 points in all. The stabiliser of the class in `S_7` is
  `<(23), (45)>`.

---------------------------------------------------------------------------------------------

## 2. Theorem 1: the cub13 family, exactly

### 2.1 Construction

Write the curve as `c(t) = (h, X_2, X_3, X_4, X_5)`, all of degree at most 3.

* **The six coordinate planes `Π_{1bc}` (`b,c <= 5`).** Each forces two of `X_2..X_5` to have a
  common root. With every coordinate edge met once, this gives
  `X_d = lam_d prod_{e != d, e in {2..5}} (t - r_{de})`, with `r_{de} = r_{ed}`.
* **Gauge.** Put the contact with `Π_{126}` at `t = ∞`, which means `lam_3 = lam_4 = lam_5 = 1`.
  Put the contact with `Π_{123}` at `t = 0` (`r_45 = 0`) and the contact with `Π_{135}` at `t = 1`
  (`r_24 = 1`).
* **`Π_{136}` (`X_2 = X_4 = X_5` at `v`).** Since `X_4 = X_5` at both `∞` and `v`, `v` is the other
  root of a quadratic with a known root:
  `v = (r25 r35 - r34)/(r25 + r35 - 1 - r34)`. Then `lam_2 = X_4(v)/((v - r23)(v - 1)(v - r25))`.
* **`Π_{456}` (`h = X_2 = X_3` at `w`).** Take `w` as a free parameter and solve `X_2(w) = X_3(w)`
  for `r_23`. The equation is Möbius in `r_23`, so

  ```text
  r23 = (r25 r35 - r34)(r25 r34 r35 - r25 r34 - 2 r25 r35 w + r25 r35 + r25 w^2 - r34 r35
         - r34 w^2 + 2 r34 w + r35 w^2 - w^2) / ((r25 - 1)(r34 - w)(r35 - w)(r25 - r34 + r35 - 1)).
  ```
* **The four planes `Π_{246}, Π_{256}, Π_{346}, Π_{356}`** (`h = X_3 = X_5`, `h = X_3 = X_4`,
  `h = X_2 = X_5`, `h = X_2 = X_4`). Their contact parameters are the remaining roots of `X_3 - X_5`,
  `X_3 - X_4`, `X_2 - X_5` and `X_2 - X_4`. In each case the other roots are already known:
  `{r35, ∞}`, `{r34, ∞}`, `{r25, v}`, `{1, v}`. So all four are rational functions of the
  parameters. `h` is the Lagrange cubic through the four points `(t_k, X_{d_k}(t_k))`.
* **The last condition** is that `h` also passes through `(w, X_2(w))`:

  ```text
  Res(r25, r34, r35, w) := h(w) - X_2(w) = 0.
  ```

Everything above is a rational function of `p = (r25, r34, r35, w)` over `Q`
(`cub13exact.build`; `kap1.py`, `kap2.py` give the same formulas in sympy). The four parameters are
contact positions in an intrinsic gauge:
* `r25 = t_134`, `r34 = t_125`, `r35 = t_124`, `w = t_456`;
* fixed by the gauge: `t_126 = ∞`, `t_123 = 0`, `t_135 = 1`.

So distinct parameter points give distinct curves.

### 2.2 Statement

**Theorem 1.** Let `m(w) = 329715 w^5 - 8401893 w^4 + 44833850 w^3 - 89897718 w^2 + 73951827 w
- 20297381`, irreducible over `Q`, let `K = Q[w]/(m)`, and let `alpha` be the class of `w`. At
`p_0 = (3, -2, 5/3, alpha)` the following hold.

1. `Res(p_0) = 0` in `K`.
2. The curve `c(t)` (after the Möbius change `t = (3s - 1)/(2s + 5)`) has reduced parametrisation.
   For each of the 56 sets `T`, the polynomial whose roots are the points where `T` is exactly a
   level set has the degree the class `cub13` prescribes:
   * degree 1 for the 13 sets of the table in Section 1 and for `{1,4}`, `{1,5}`, `{2,3}`;
   * degree 3 for `{1,6}`;
   * degree 0 for the other 39.

   These 17 polynomials are squarefree and pairwise coprime. Every root is a simple root of each
   difference `c_i - c_j` with `i, j in T`. There is no contact at `s = ∞`. Hence `c` is a rational
   curve in `X` of class `cub13` that meets the boundary in 19 distinct points, each transversally
   and on a single divisor.
3. `grad Res(p_0) != 0`. So `{Res = 0}` is a smooth hypersurface near `p_0`, of dimension 3.
4. The differential of `(s, p) -> c(s; p)` on `P^1 x {Res = 0}`, at `s_0 = 2/7`, has rank 4.

**Corollary 1.1.**
* Let `W_0` be the irreducible component, over `Q-bar`, of the zero set of the numerator of `Res`
  that contains `p_0`, and let `W` be the union of its `Gal(Q-bar/Q)`-conjugates. Then `W` is
  defined over `Q` and is a 3-dimensional family of rational curves of class `cub13`.
* Each component of `W` dominates `X` through `P^1 x W -> X`, so its general member is a free curve.
* If the genuine factor `G` of the numerator has multidegree exactly `(7,5,5,5)` (Data on coordinate
  lines together with a NUMERICAL genuineness test, Section 3), then `W` is geometrically
  irreducible.

*Proof of Theorem 1.* Each item is an exact computation in `K`, implemented in `verify_cub13.py`:
* `nf.py` provides arithmetic in `K` by Euclid's algorithm modulo `m`, polynomials over `K` with
  gcd, and forward-mode automatic differentiation over `K`.
* Item 2 is decided by gcds in `K[s]`. Let `P_ij = c_i - c_j`. For each `T`, take
  `g_T = gcd_{i,j in T} P_ij`, then remove every root shared with some `P_{i_0 k}` with `k ∉ T`.
* Items 3 and 4 use exact derivatives propagated through every formula of 2.1, including the
  Lagrange interpolation.
* Irreducibility of `m` is `sympy.factor_list` over `Q`.
* For items 2–4 the statement over `K` is a statement about all five embeddings at once, because
  the computations are identities in `K`. `[]`

*Proof of Corollary 1.1.*
* By item 3, the zero set of `Res` is a smooth 3-fold near `p_0`. The class is constant on
  `W_0`, since curve classes are deformation invariant. The single-divisor transversal property
  holds at `p_0` by item 2, and it is open. So it holds on a Zariski open dense subset of `W_0`.
* By item 4 and the implicit function theorem over `C`, `P^1 x W_0 -> X` is dominant. A general
  member of a dominating family of rational curves is free (Kollár, *Rational curves on algebraic
  varieties*, II.3.11).
* `W` is `Gal`-stable, hence defined over `Q`. The conjugate components dominate as well, because
  Galois conjugation preserves dominance.
* Geometric irreducibility: conjugate components have equal multidegrees. So the number of
  components divides each of `7, 5, 5, 5`, and therefore equals 1. `[]`

*Run time.* `verify_cub13.py minpoly_3_-2_5o3.json` takes 2 s. It prints `True` for items 1–4.

### 2.3 Consequence over C

**Corollary 2.** Over `C`, free rational curves of class `131040 beta` exist and cover `X`.

*Proof.* This is the argument of `referee_L7.md` §3.3, whose only non-rigorous input was the
existence of a covering cub13 family. Corollary 1.1 now supplies that family.
* Exactly, `26 * 5040 * beta = 70 sum_{sigma in S_7} sigma_* cub13 + 21 sum_sigma sigma_* f`, where
  `f` is a forgetful fibre (`referee_L7_checks/cubic_check.py`).
* All the families involved cover `X`, and their general members are free.
* Glue `70*5040` translates of general cub13 curves and `21*5040` general fibres into a tree,
  attaching each new component at a general point of one already placed. This gives a tree of free
  rational curves.
* Such a tree smooths to a free rational curve of the total class, and a general deformation avoids
  any given codimension-2 set (Kollár II.7.6, II.3.7; Debarre, *Higher-dimensional algebraic
  geometry*, §4.2).
* So the general member meets every `D_S` only at single-divisor points, and the family covers
  `X`. `[]`

Corollary 2 is the geometric (function-field) form of "`(L_7)` is false": over `C(t)`, balanced
configurations off any proper closed `Z` exist. Therefore any proof of `(L_7)` must use arithmetic of
`Q`. Liouville and auxiliary-polynomial arguments are already excluded by `threshold/L8.md`
Theorem A.

### 2.4 Remark: the class-`q beta` curves have smooth Q-points

**Proposition 2.4.** For every `k >= 1` there is a class `g_k` with `g_k / (131040 k) -> beta` and a
component `M_k` of `Mbar_{0,0}(X, g_k)` such that:
* `M_k` is defined over `Q`;
* `M_k` has a smooth Q-point;
* its general member (over `C`) is a free rational curve meeting the boundary only at
  single-divisor points.

*Proof.*
1. **Teeth.** Members of `W` over number fields come in Galois orbits; Theorem 1 gives an orbit of
   5. Other lines `(a, b, c, w)`, with `(a, b, c)` rational, give infinitely many further orbits.
   Choose `14k` distinct such orbits of general (free) members. Their `S_7`-translates, together
   with forgetful fibres over Q-points, form Galois-stable finite sets of free curves.
2. **Handle.** Choose `s_0 in Q` and the Galois-stable sets of points `C_i(s_0)` together with
   their translates. By Lagrange interpolation with Q-rational node polynomials, choose a free
   Q-rational curve `H` (a rational curve in the Kapranov `P^4` of large degree) through all of
   these points.
3. **The comb.** Attach the teeth to `H` at these points: `70 k` cub13 curves of type `sigma`
   (`14k` orbits) and `21 k` fibres of type `sigma`, for each `sigma`, as in Corollary 2. The comb is
   `Gal`-stable, so it is a Q-point of the Kontsevich space.
   * It lies on a unique geometric component, which is therefore `Gal`-stable, so defined over
     `Q` and geometrically irreducible.
   * If the orbits of Theorem 1 cannot be varied freely enough, a tooth may be repeated in a
     different Galois orbit of points of `H`.
4. **Smoothness.** `H^1(C, f^*T_X) = 0` for a comb of free curves (Kollár II.7), so the comb is a
   smooth point.
5. **The class.** `g_k = [H] + 131040 k beta`, and `[H]` is fixed while `k -> infinity`. `[]`

So `(L_7)`, if true, would assert that the Q-points of a geometrically irreducible variety with
smooth Q-points (`M_k`, of dimension `-K.g_k + 1 >= 7 * 131040 k`) never give curves off `Z`. That
cannot be excluded in general, since varieties of general type behave this way. But `M_k`
parametrises rational curves on a rational variety, and no mechanism forcing such sparsity is known.
**The gap is the density of `M_k(Q)`.** Section 4 explains why Q-points of the cub13 family do not
close it.

---------------------------------------------------------------------------------------------

## 3. Arithmetic of the cub13 family

**Proposition 3.1 (cub12 is Q-rational).** Drop the condition at `Π_456`. The parameters
`(r23, r25, r34, r35)`, with the same gauge and the same formulas for `v`, `lam_2`, the four
contact parameters and the Lagrange cubic `h`, parametrise a 4-dimensional family `V` of rational
curves, Q-rationally. `V` dominates `X`. Its class `cub12` meets the 12 remaining planes once each:
* `D_23 = 2`, `D_12 = D_13 = D_14 = D_15 = 1`, `D_16 = 3`;
* `pair = 9`, `triple = 12`, `(B+8K) = -3` (`fibcone.py`).

`W` is the hypersurface `{exists t : h(t) = X_2(t) = X_3(t)}` in `V`.

*Proof.*
* Every step in 2.1 except the last is a rational function of `(r23, r25, r34, r35)`.
* `verify_cub12.py` checks exactly, at the Q-points `(r25, r34, r35, w) = (3, -2, 5/3, 7/2)` and
  `(2, 5/7, -3, 1/3)` (with `r23` from `w`), that the curve has class `cub12` and 21 distinct
  transversal single-divisor contacts. These are honest Q-rational curves.
* Dominance follows from Theorem 1, item 4, since `W ⊂ V`.
* The class is constant on `V` and the contact property is open, as in Corollary 1.1. `[]`

**Data and NUMERICAL on `W`.**

| quantity | value | type / script |
|---|---|---|
| degrees of the genuine factor `G` in `(r25, r34, r35, w)` | `(7, 5, 5, 5)` | Data (`kap3.py`), genuineness NUMERICAL (`kap5.py`) |
| all other factors of the numerator of `Res` on these lines | degenerate (colliding contacts) | NUMERICAL (`kap5.py`) |
| degree of `W -> A^3`, `(r25, r34, phi)`, for `phi` in `{r35, w, v, r23, lam2, tD23, t246, t256, t346, t356}` | between 4 and 8 | resultants exact, genuineness NUMERICAL (`projdeg.py`) |
| degree of the contact map `W -> E_126` | `>= 4` (4 fibre points found) | NUMERICAL (`fibre126.py`) |
| curves of class cub13 through a general point | `>= 8` | NUMERICAL (`referee_L7.md` §3.2) |
| `deg(ev_{s_0}: V -> X)` for cub12 at fixed intrinsic parameter `s_0` | `>= 4` | NUMERICAL (`cub12ev.py`) |
| Q-points of height `<= 6` with `Res = 0` | 44 found, **all degenerate** | Data (`fastsearch.py`, `verify_cub13.py`) |

* **The search** runs over all `(r25, r34, r35, w)` with entries `a/b`, `|a|, b <= 6`, excluding
  `0` and `1`: `45^4 = 4,100,625` points, of which 3,562,865 are not poles of `Res`. Each is tested
  by exact `Fraction` arithmetic.
* **The 44 zeros.** Each was checked by the exact verifier over `Q`. Every one fails item 2. Typically:
  * the `D_23` contact coincides with the `Π_123` contact (the curve passes through
    `D_23 ∩ D_456`), or
  * a root of `h` falls on a plane contact.

  These points lie on the spurious factors of the numerator. Those factors are rational
  hypersurfaces, which is why they have small rational points.
* **No Q-rational cub13 curve is known.**
* No natural projection of `W` is birational, so the rationality of `W` over `Q`, or even the density
  of `W(Q)`, remains open. This question is moot for `(L_7)` by Section 4.

---------------------------------------------------------------------------------------------

## 4. Adding families over Q: an addability criterion

### 4.1 Definitions

**Q-family.** A *Q-family of rational curves* `V` consists of the following:
* an irreducible Q-variety `B_V`;
* `U_V = P^1 x B_V`;
* a rational map `ev: U_V ⇢ X`, defined over `Q`, such that for general `b` the map
  `ev(., b): P^1 -> X` is birational onto a curve that meets the boundary only at single-divisor
  points.

Its class `[V]` is the class of a general member. Every family in this note is of this form: `W`,
`V = cub12`, forgetful fibres, and so on.

**Addable.** `V` is *addable* if `ev` is dominant and has a rational section over `Q`. That is a
rational map `sigma: X ⇢ U_V`, defined over `Q`, with `ev o sigma = id`. Equivalently, over the
generic point `eta` of `X` some member passes through `eta` and is defined over `Q(X)`, together with
the position of `eta` on it.

Two remarks.
* If `V` is addable, so is `tau(V)` for every `tau in S_7`, with section `tau o sigma o tau^{-1}`.
* A fibration is addable. If `phi: X ⇢ Z` is a Q-rational map to a 3-fold with geometrically
  rational generic fibre `X_eta`, and `X_eta(Q(Z)) != ∅`, then `X_eta = P^1_{Q(Z)}`. Taking `B_V = Z`
  (birationally) and `U_V = X`, the identity is a section.

### 4.2 Lemma 4.2: adding one family

**Lemma 4.2.** Let `V` be addable with section `sigma`. Let `f: P^1_Q -> X` be birational onto `C`,
with only single-divisor boundary contacts. For `c in C` let `M_c` be the member
`ev(P^1 x {pr_B sigma(c)})`, and `s_c` the position of `c` on it. Assume:
* (G1) `C` avoids the indeterminacy locus of `sigma`, and the map
  `iota: C x P^1 ⇢ X`, `(c, s) -> M_c(s)`, is defined at every point of the section
  `Sigma = {(c, s_c)}`;
* (G2) only finitely many `M_c` (`c in C`) meet a codimension-2 boundary stratum, and none is
  contained in one;
* (G3) `iota` is birational onto its image `S_C`.

Then for every `m >= 1` there are morphisms `f_m: P^1_Q -> X` with the following properties:
* each `f_m` is birational onto its image;
* each has only single-divisor contacts;
* each has class `[C] + m [V]`;
* their images are Zariski dense in the surface `S_C`.

*Proof.*
1. **The surface.** `Y = C x P^1 = P^1 x P^1` over `Q`, with fibre class `F = {c} x P^1` and
   horizontal class `h`. The section `Sigma` is the graph of `c -> s_c`, of class `h + delta F`
   for some `delta >= 0`. Then `iota(Sigma) = C` and `iota({c} x P^1) = M_c`.
2. **Resolution.** `iota` is undefined at finitely many points, none of them on `Sigma` (G1). Blow
   them up to get `iota~: S~ -> X`. Since the blown-up points are off `Sigma`, `Sigma~` is the
   total transform of `Sigma`.
3. **Sections.** For `m >= 1`, `D := Sigma~ + mF` is the pullback of `h + (delta + m) F`. On
   `P^1 x P^1` that class is base-point-free and very ample, since `delta + m >= 1`. So `|D|` is
   base-point-free and maps `S~` birationally to its image.
   * Its general member is the strict transform of a general curve in
     `|h + (delta + m) F| ≅ P^N_Q`. Such a curve is the graph of a map `C -> P^1` of degree
     `delta + m`, so it is a section, isomorphic to `P^1_Q`, and it avoids any given finite set of
     points.
   * Call it `sigma'`. Then `iota~_* [sigma'] = iota~_* [Sigma~] + m iota~_* F = [C] + m [V]`.
4. **Contacts.** Outside the exceptional curves, `iota~^{-1}(D_S ∩ D_{S'})` meets the strict
   transforms of fibres only where `M_c` meets a codimension-2 stratum. By (G2) this is a finite
   set of points. A general `sigma'` avoids this finite set. It is also disjoint from the
   exceptional curves, being the strict transform of a curve in `Y` that misses the blown-up
   points. So `iota~ o sigma'` meets the boundary only at single-divisor points.
5. **Birationality and density.** Birationality onto the image follows from (G3). The Q-points of
   `P^N_Q` are dense, and the curves sweep out `S~`, so their images are dense in `S_C`. `[]`

### 4.3 Theorem 3

**Theorem 3 (addability criterion).** Suppose there is an addable Q-family `V` with
`(B+8K).[V] >= 0`. Suppose also that at every step of the iteration below, the curves of `F_{j-1}`
satisfying (G1)–(G3), for the family being added, have Zariski-dense union. Then `(L_7)` is false.

*Proof.*
1. **The coefficients.** Since `N_1(X)^{S_7}` is spanned by the duals of the pair and triple
   boundary sums (Keel–McKernan), `sym(g)` has `D_S`-degree `pair(g)/21` on pair divisors and
   `triple(g)/35` on triple divisors. Put
   `a = 35/triple(V)` and `b = (21 - a pair(V))/6`. Then `b >= 0` because `V` is rich, and
   `beta = a sym([V]) + b sym(f)`.
2. **Seed.** Let `F_0` be the set of forgetful fibres `pi_7^{-1}(z)`, `z in Mbar_{0,6}(Q)` general.
   These are Q-curves with 6 single contacts, and their union is dense.
3. **Induction.** List the `2*5040` addable families `tau(V)`, `tau(f)` (`tau in S_7`). Give
   `tau(V)` the multiplicity `m_tau = max(1, round(a L / 5040))` and `tau(f)` the multiplicity
   `n_tau = max(1, round(b L / 5040))`.
   * Step `j` applies Lemma 4.2 to the curves of `F_{j-1}` that satisfy (G1)–(G3). By hypothesis
     their union is dense.
   * The new curves are dense in each `S_C`, and `S_C ⊇ C`. So `F_j`, the set of all new curves,
     again has dense union, and its members are Q-rational with single-divisor contacts.
4. **Conclusion.** After the last step every curve has class `L beta + R`. The remainder `R` is
   bounded independently of `L`: it collects the seed and at most one unit of rounding per family.
   So `eps_f = O(1/L)`. As `L -> infinity`, Theorem 2(a) of round 4 applies. `[]`

*Remarks.*
1. Theorem 3 needs **one** rich addable family. Its `S_7`-translates and the forgetful fibres do the
   rest.
   * The genericity hypotheses (G1)–(G3) fail only when a curve meets a fixed codimension-2 locus
     or is itself special. They are expected to hold.
   * They are part of the statement because no rich addable family is at hand on which to check
     them.
   * For the forgetful fibres they can be checked directly: `iota` is then the inverse of a
     forgetful map restricted to `pi^{-1}(pi(C))`.
2. The proof of `(b) ⊂ (a)` in `referee_L7.md` §2.2 is a different mechanism: curves inside a
   product chart. Theorem 3 does not need a global product structure. Different fibrations can be
   combined one step at a time.
3. Toric "enrichment" is a special case of Lemma 4.2, in which the curve is multiplied by `f(t)^{e_A}`
   in a torus chart. The 1-parameter-subgroup orbits `gamma_{a,b,A}` are the fibres of the fibration
   `T_ab -> T_ab/G_m^{e_A}`.

### 4.4 What is addable

**Proposition 4.4.**
1. **Fibrations.** The forgetful fibres `f_i`, the fibres of the plane projections of every
   Kapranov model, and the toric orbits `gamma_{a,b,A}` are addable. Data (`fibcone.py`):

   | addable class | pair | triple | `-K` | `(B+8K)` | `(B+8K)/B` |
   |---|---|---|---|---|---|
   | `f_i` (forgetful fibre) | 6 | 0 | 2 | -10 | -5/3 |
   | line through a general point (projection from a point) | 15 | 0 | 5 | -25 | -5/3 |
   | line meeting 3 Kapranov planes (projection from 3 planes) | 6 | 3 | 2 | -7 | -7/9 |
   | line meeting `L_ab` and `Π_cde` | 6 | 2 | 2 | -8 | -1 |
   | toric `gamma_{a,b,A}`, `abs(A) = 2` | 6 | 2 | 2 | -8 | -1 |

   So every known addable class satisfies `(B+8K).g <= -(7/9) B.g`, and the closed cone they
   generate does not contain `beta`.
2. **3-dimensional families.** A Q-family `V` with `dim B_V = 3` that dominates `X` is addable if and
   only if exactly one member passes through a general point. In that case `X` is fibred by its
   members.
3. **cub13 is not addable** if, as the round-4 referee found NUMERICALLY, at least 8 members pass
   through a general point. Even a Zariski-dense set of Q-rational cub13 curves would not change
   this.
4. **cub12.** This Q-rational 4-dimensional family is not rich: `(B+8K) = -3`. For a fixed intrinsic
   parameter `s_0`, `ev_{s_0}: B_V -> X` has degree `>= 4` (NUMERICAL, `cub12ev.py`), so this
   obvious candidate section does not exist.

*Proof.*
1. For a fibration the identity is the section (4.1). The table is computed from the Kapranov data
   `(e; m_a, m_ab, m_abc)`. A plane projection `P^4 ⇢ (P^1)^3` (three planes) or `P^2 x P^1` (a line
   and a plane) has as generic fibre a line in `P^4`, defined over the function field of the base. A
   line over any field is `P^1`. The toric orbits are the fibres of the quotient
   `T_ab ⇢ T_ab / G_m`, which has a section because `G_m`-torsors are trivial (Hilbert 90).
   Finally, `(B+8K).g <= -(7/9) B.g` for every generator gives the same inequality on the closed
   cone, while `(B+8K).beta = 0 < B.beta`.
2. `ev: P^1 x B_V -> X` is generically finite, of degree `N`. A rational section exists only if some
   component of `U_V` has degree 1 over `X`, and `U_V` is irreducible.
3. This follows from 2.
4. `(B+8K)` is computed in `fibcone.py`. `[]`

### 4.5 Q-points of a fixed seed do not matter

**Proposition 4.5.** Any curve produced by the operations above has class
`[seed] + sum_j m_j [V_j]`, with `V_j` addable. Such curves approach the ray of `beta`
(`eps -> 0`) only if `beta` lies in the closed cone spanned by the `[V_j]`, unless the seed itself is
already approximately balanced.

*Proof.* The seed class is fixed while `sum m_j -> infinity`. Divide by `sum m_j` and pass to the
limit. `[]`

**Consequence for the round-4 gap.** The round-4 referee wrote "the only missing ingredient for a
refutation via (a) is Zariski density of Q-rational members". That is correct for the members of
the *class-`q beta` component* `M_k` of Proposition 2.4. It is not correct for the cub13 family:
* Q-points of `W` can only serve as seeds or as teeth;
* teeth enter only through smoothing, which is not available over `Q`;
* seeds drop out by Proposition 4.5.

Two mechanisms are known that turn Q-rational small curves into Q-rational large ones:
* Theorem 2(b) needs a Q-rational product chart whose ruling cone contains `beta`;
* Theorem 3 needs a rich family with a rational section of `ev`.

Neither follows from Q-points of `W`. Through a general point, cub13 has 8 members, permuted
transitively by monodromy (referee §3.2), so `W` admits no section. The chart `P^1 x W` would give
only one cub13 ruling, plus whatever classes the coordinates of `W` produce, if `W` were Q-rational,
which is unknown.

---------------------------------------------------------------------------------------------

## 5. Barriers for natural Q-constructions

**Theorem 5.**

**(a) Group-law product charts.** Fix labels `a != b`, and let `T = T_ab` be the torus chart: `y_a = 0`,
`y_b = ∞`, coordinates `y_i` for `i ∉ {a,b}`, modulo scaling. Call `D_S` *hidden* if
`S ⊆ [7] \ {a,b}` and `abs(S) >= 3`; there are 16 such `S`, including `S = [7] \ {a,b}`, which gives
`D_ab`. Let `phi_1, ..., phi_R: P^1 ⇢ T` be rational curves, and
`Psi(y) = prod_k phi_k(y_k)` (coordinatewise product). Suppose that for each `k` the product map
`prod_{j != k} phi_j` has Zariski-dense image in `T`. Then every ruling class of `Psi` has
`D_S`-degree 0 for all 16 hidden `S`. So `beta` is not in the ruling cone, and Theorem 2(b) cannot
be applied to `Psi`. The same holds for the additive analogue: the affine chart with `y_b = ∞`,
translation by `G_a^6` modulo `Aff`, and hidden meaning `S ⊆ [7] \ {b}` with `abs(S) >= 3`.

**(b) The Galois–Boolean mechanism at `k = 7` has genus `>= 5`.** Take distinct anchors
`x_1, x_2, x_3` and an unmarked common argument `x_4`, and degree-2 rational functions `f_a, f_b,
f_c, f_d` from one pencil with the common arguments `x_1, ..., x_4`, as in `threshold/L8.md` §3.1.
Let `E` be a component of the normalisation of `{f_a(a) = f_b(b) = f_c(c) = f_d(d)}`, and
`Psi: E -> X`, `(a,b,c,d) -> (x_1, x_2, x_3, a, b, c, d)`. If `Psi(E)` meets each of the 15
divisors `D_{x_1 ∪ B}` (`∅ != B ⊆ {a,b,c,d}`), then `E` has genus `>= 5`. So `E(Q)` is finite
(Faltings), and such curves cannot refute `(L_7)`.

**(c) Lifting a seventh point by interpolation is overdetermined** (round 4, §8, in curve form).
Let `C_6 ⊂ Mbar_{0,6}` have class `~ 2 q beta_6`. A seventh coordinate `x_7`, of degree `~ 8q` in
an anchor frame, must join half of the contacts of `C_6`. That is `~ 25 q` linear conditions on
`16 q + 2` coefficients. So balanced lifts exist only over a locus of `C_6` of codimension
`~ 9 q - 1`, and no linear-algebra construction produces them for general `C_6`.

*Proof of (a).*
* `Delta_S = {y_i = y_j for all i, j in S}` is a subtorus of codimension `abs(S) - 1 >= 2`, and it
  is the image in `T` of the hidden divisor `D_S`.
* Fix a ruling direction `k` and let `Gamma = phi_k(P^1) ∩ T`.
  * The set `{tau : tau Gamma ∩ Delta_S != ∅} = Delta_S Gamma^{-1}` has dimension at most 3, so it
    is a proper closed subset.
  * At the finitely many points where `tau Gamma` leaves `T`, the limit point is determined by the
    leading coefficients multiplied by the corresponding coordinates of `tau`. It lies on a hidden
    `D_S` only on a proper closed set of `tau`.
* For general `tau` in the dense image of `prod_{j != k} phi_j`, the translate `tau Gamma` therefore
  meets no hidden divisor.
* The ruling class is the class of the general translate, so it has hidden degree 0. But `beta`
  has hidden degree 1 on all 16 hidden divisors, and a nonnegative combination of classes with
  hidden degree 0 has hidden degree 0.
* The additive case is identical, with subgroups of `G_a^6/Aff` in place of subtori. `[]`

*Proof of (b).*
* Let `T` be the common value. Each `m in {a,b,c,d}` gives a double cover `P^1_m -> P^1_T`. A
  mover equals an anchor `x_j` only over `T = y_j := f(x_j)`, because the pencil members agree only
  at the common arguments (L8 Lemma B.1(2)).
* A contact with `D_{x_1 ∪ B}` lies over `y_1`, at a point of `E` whose *pattern* (the set of movers
  equal to `x_1`) is `B`. So all 15 nonempty patterns occur over `y_1`.
* `y_1` is not a branch value of any `f_m`. If it were, both preimages of `y_1` under `f_m` would
  equal `x_1`, every pattern over `y_1` would contain `m`, and the pattern `{m'}` with `m' != m` could
  not occur.
* `E -> P^1_T` is Galois with group `G ⊆ (Z/2)^4` (`E` is a component of the fibre product of the
  four double covers). Since `y_1` is unramified, the fibre over it is a `G`-torsor, so it has
  `abs(G)` points. Hence `abs(G) >= 15`, so `G = (Z/2)^4`, and the four quadratic subextensions are
  independent.
* Over `Q-bar(T)`, square classes unramified outside a set of `n` points form an `F_2`-space of
  dimension `n - 1` (even subsets). So the set of branch values has `n >= 5` elements.
* All inertia groups are of order 2. Riemann–Hurwitz gives `2g - 2 = 16 * (-2) + 8n`, that is
  `g = 1 + 4(n - 4) >= 5`. `[]`

*Remarks.*
1. **Novelty of (b).** The genus bound for four independent quadratic covers is
   `threshold/L8.md` §5, item 112 (`seven_row_quadratic_quotient_obstruction.md`). New here:
   * balance forces independence, because dependent square classes give `abs(G) <= 8` and so at
     most 8 of the 15 patterns;
   * the branch value `y_1` is forced to be unramified.

   The L6 trick (`k = 6`: three movers, `G = (Z/2)^3`, `n = 4`, genus 1) therefore cannot extend in
   this form. The round-4 discussion (`round4/L7.md` Corollary 6) covered only class exactly `beta`,
   by parameter count.
2. (a) explains the failure of every attempt that combines families by the group law of a torus
   chart. Hidden contacts survive a product only at parameters where both factors are synchronised
   on the same subtorus. The toric enrichment of Remark 3 after Theorem 3 preserves the existing
   hidden contacts but never creates new ones.

---------------------------------------------------------------------------------------------

## 6. What a refutation now needs, and the first candidates

By Theorem 3, `(L_7)` is refuted by **one rich addable Q-family**, subject to (G1)–(G3). By
Theorem 5 and Proposition 4.4, no family found so far qualifies.

**6.1 Necessary conditions.** Let `V` be rich and addable.
* If `dim B_V = 3`, then `X` is fibred by the members of `V`: exactly one member passes through a
  general point (Proposition 4.4(2)).
* If `dim B_V >= 4`, the members through the generic point `eta` form a variety `Phi_eta` of
  dimension `dim B_V - 3` over `Q(X)`, and it must have a `Q(X)`-point. For `dim B_V = 4`,
  `Phi_eta` is a curve.
  * One sufficient condition is that `Phi_eta` has genus 0 and a canonical `Q(X)`-point, for
    example a degenerate member through `eta` that is unique of its type.
  * The simplest sufficient condition is that `ev_{s_0}: B_V -> X` is birational for some fixed
    `s_0`. It fails for cub12 (degree `>= 4`).
* Every rich covering class `g` has `psi_i`-degree `e_i >= 3` in every Kapranov model.
  * In the `psi_i`-model, with multiplicities `m_a`, `m_ab`, `m_abc` at points, lines and planes,
    `-K.g = 5e - 3 sum m_a - 2 sum m_ab - sum m_abc` and `triple(g) = sum m_ab + sum m_abc`. Hence
    `triple(g) <= 5e - (-K.g)`.
  * Richness means `triple(g) >= 5(-K.g)`, so `e >= (6/5)(-K.g) >= 12/5`, using `-K.g >= 2` for
    covering families.
* A full component with `dim B_V = 4` has `-K.[V] = 3`.

**6.2 Rich families of small degree** (Kapranov `psi_7`-model, incidences with planes only).
* A degree-`e` rational curve meeting `M` planes has `pair = 15e - 3M` and `triple = M`.
* It is rich iff `M >= 25e/6`, and the family has expected dimension `5e + 1 - M`. So:
  * `e = 3`: `M = 13`, dimension 3. cub13 is one such choice of 13 planes; other admissible
    13-sets give further 3-dimensional rich classes. All have `-K = 2`, so by 6.1 they are addable
    only if they are fibrations.
  * `e = 4`: `M = 17`, dimension 4. This is the first rich family of dimension `>= 4`.
  * `e = 5`: `M = 21` (dimension 5) or `M = 22` (dimension 4).
  * Larger `e` give rich families up to dimension about `5e/6 + 1`.

**The quartic family `Q4`** (explicit; existence not yet established).
* **Class.** All 10 coordinate planes and the 7 hidden planes with level sets
  `{134, 135, 145, 234, 235, 245, 345}` (all triples of `[5]` not containing `{1,2}`).
  Then `pair = 9` (`D_12 = 3`; `D_13, D_14, D_15, D_23, D_24, D_25 = 1`), `triple = 17`, and
  `(B+8K) = +2`.
* **Shape.** Every coordinate is a product over the edges of `K_5`:
  `c_m = lam_m prod_{e != m}(t - r_{me})`.
* **Gauge and Vieta structure.** Put `t_345 = ∞`, so `lam_3 = lam_4 = lam_5 = 1`. Then
  `(c_3 - c_4)/(t - r_34)` is a quadratic whose roots are exactly `t_134`, `t_234`, because
  `D_34 = 0`. Likewise `(c_3 - c_5)/(t - r_35)` has roots `t_135`, `t_235`, and
  `(c_4 - c_5)/(t - r_45)` has roots `t_145`, `t_245`.
* **Remaining conditions.** Label 1 must satisfy `c_1 = c_3` at `t_134` and `t_135`, and
  `c_1 = c_4` at `t_145`. Label 2 must satisfy the same at the conjugate roots.
* **Count.** 16 unknowns after gauge, 12 conditions, dimension 4.
* **Projection from `p_1`.** This reduces `Q4` to a 6-dimensional family of quartics in `P^3`. The
  positions `t_1xy` are then Vieta-rational, and two linear-in-`lam_1` lifting conditions remain,
  exactly as for cub13 (Section 2.1).
* **Numerics.** Damped Newton from 700 random starts (`q4.py`, `q4b.py`) and Newton homotopies
  (`q4h.py`) converge only to degenerate strata, where parameters collide (`q4bdbg*.py`). This is
  the same failure mode as round 4 at `k = 7`, and it is NOT evidence of emptiness.
* **Open questions for `Q4`.**
  * (i) Does it exist? Its expected dimension is 4, and its class has all `D_S >= 0`.
  * (ii) Is it Q-rational?
  * (iii) Does `Phi_eta` have a `Q(X)`-point?
* **What a positive answer gives.** If (i) and (iii) hold, then Theorem 3 refutes `(L_7)`, since
  `beta = (35/17) sym(Q4) + (7/17) sym(f)`. Check: `(35/17)(9, 17) + (7/17)(6, 0) = (21, 35)`.

**6.3 Rich fibrations.** A rich class `g` that is the fibre class of a Q-rational map `X ⇢ Z^3` would
suffice by itself. Fibre classes need not be minimal: lines through a general point have `-K = 5`.
Such a map is a congruence of order one of rich rational curves, that is, exactly one member
through a general point. All the classical ones in the Kapranov models are poor
(Proposition 4.4(1)). A Cremona transformation `chi: X ⇢ P^4` under which all 56 boundary divisors
become hypersurfaces of equal degree would make the lines of `P^4` balanced. That is the `k = 7`
analogue of `threshold/L8.md` Proposition C, the integral chart for `k = 5`. Whether such `chi`
exists is open.

---------------------------------------------------------------------------------------------

## 7. What the outcome means for the uniform theorem and for growth

* **No change of status.**
  * Proved: `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`).
  * Conditional: Vojta for one fixed variety gives a uniform bound (`growth/conditional.md`).
  * `(L_7)` remains open.
* **What `(L_7)` would do.** It is the `k = 7` instance of the missing lemma of the Ptolemy route
  (`ptolemy.md` §7, `threshold/L8.md` §1). If true, it would give the uniform theorem (via the
  Ptolemy steps), and moreover at a `k` where the `D = 0` Vojta inequality is vacuous
  (`K.beta = -7`).
* **If `(L_7)` is false**, which every piece of evidence here suggests, the route must work at
  `k >= 8`. There:
  * Liouville-type arguments are impossible (L8 Theorem A);
  * the needed input is `V_8`, a Vojta-type inequality.
* **The new information on which side `(L_7)` falls.**
  1. The geometric obstruction is gone. Over `C`, free balanced rational curves exist and cover
     `X` (Corollary 2), and the relevant component of curves has smooth Q-points
     (Proposition 2.4). A proof of `(L_7)` would therefore have to be purely arithmetic: Q-points of
     a variety with smooth Q-points avoiding the balanced locus. It would also have to use the
     Ptolemy structure; Theorem A of L8 rules out auxiliary polynomials.
  2. A disproof is reduced to one concrete algebro-geometric object over `Q`: a rich addable family
     (Theorem 3), with the first candidate written down (`Q4`, Section 6.2).
* **Growth.** `(L_7)` concerns whether `H_Sigma >= c w` for balanced anchored configurations. A
  refutation would show that the Ptolemy/S-unit system alone, at `k = 7`, gives no growing lower
  bound `H_Sigma >= f(w)`, `f -> infinity`. Corollary 2 already proves this for the function-field
  analogue. For the circle problem, the Gaussian data that the Ptolemy system forgets (isotropy
  `rho_T^2 = -1`, norms, singleton blocks; `ptolemy.md` §7 item 4) would then be indispensable at
  `k = 7`. This is consistent with the round-5/6 fake-circle barriers, which show that
  residue-level data alone stops at `W = O(M log M)`.

---------------------------------------------------------------------------------------------

## 8. Files (`scratchpad/round7/L7_checks/r7/`)

| file | content | type |
|---|---|---|
| `kap1.py`, `kap2.py` | sympy derivation of the parametrisation of Section 2.1 (`kap2.pkl`) | exact |
| `cub13exact.py` | the construction as plain field arithmetic (works over `Q`, a number field, or AD) | exact |
| `nf.py` | number field `Q[x]/(m)`, polynomials over it (gcd, derivative), forward-mode AD | exact |
| `minpoly.py`, `minpoly_3_-2_5o3.json` | the quintic `m(w)` on the line `(3, -2, 5/3, w)` | exact |
| `verify_cub13.py` | **Theorem 1**: residual, contact structure via gcds, smoothness, rank of differential (over `K` or over `Q`) | exact |
| `verify_cub12.py` | Proposition 3.1: exact class `cub12` of two Q-rational members | exact |
| `stab.py` | stabiliser `<(23),(45)>` of the cub13 class in `S_7` | exact |
| `kap3.py`, `kap5.py`, `kapcheck.py`, `kapcurve.py` | factor degrees of the numerator of `Res` on coordinate lines; NUMERICAL genuineness | exact + NUMERICAL |
| `slice2.py`, `projdeg.py` | bivariate slices; degrees of projections `W -> A^3` | exact elimination + NUMERICAL |
| `fibre126.py` | contact map `W -> E_126` has degree `>= 4` | NUMERICAL |
| `fastsearch.py`, `fs6_*.log`, `fs6_zeros.txt`, `classify_zeros.py` | height-6 search for Q-points of `W`; 44 zeros, all degenerate | exact |
| `cub12ev.py` | `deg ev_{s_0} >= 4` for the Q-rational family cub12 | NUMERICAL |
| `fibcone.py` | classes and `(B+8K)` of the addable fibration classes, cub12, cub13, `beta` | exact |
| `levelsets.py`, `q4.py`, `q4b.py`, `q4h.py`, `q4bdbg*.py` | attempted numerical construction of `Q4` (only degenerate solutions found) | NUMERICAL |
| `genus.py`, `genus_slice.py` | monodromy genus of slices of `W` (did not finish in 25 min; no result claimed) | NUMERICAL |
| `NOTES.md` | working notes | — |

Run each script from `L7_checks/r7/` with `python3`. `verify_cub13.py minpoly_3_-2_5o3.json`
takes about 2 s. `fastsearch.py 6 p 4` (`p = 0..3`) takes about 22 min per part.

**References.**
* Kollár, *Rational curves on algebraic varieties*: II.3.7, II.3.11 (free curves), II.7 (combs and
  smoothing).
* Debarre, *Higher-dimensional algebraic geometry*, §4.2.
* Keel–McKernan, *Contractible extremal rays on Mbar_{0,n}* (`S_n`-invariant classes).
* Kapranov, *Veronese curves and Grothendieck–Knudsen moduli space*.
* Faltings (finiteness of rational points on curves of genus `>= 2`).
* Project: `round4/L7.md` (Theorems 1, 2; Proposition 3), `round4/referee_L7.md` (§§2.2, 3.2,
  3.3), `threshold/L8.md` (Theorems A, B; §5), `growth/ptolemy.md` §7, research note item 112
  (`seven_row_quadratic_quotient_obstruction.md`).
