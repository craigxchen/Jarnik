# Round 12 — adversarial fake test for the T = 1 programme (route 5)

## 0. Summary of results

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. The task was the adversarial
fake test of the `T = 1` programme (R11 §8.3 item 4, §8.4). Results:

1. **The fake test reduces exactly to the tangent cone of the fake pattern set (Theorem 2.1, proved).** Fake
   positions are free up to the gap inequalities, and `A`, `m`, `V_H` depend only on the pattern set (and `r`). Every
   nondegenerate real point `u*` of `TC_H(A,m)` of a strictly admissible fake pattern set (skew double core, round-5
   Paley twin profile) is, to within `theta e^{-W/4}`, the direction of genuine gap-admissible fakes, at every large
   prime window; these fakes satisfy Step 0's conclusion `kappa_c(hat u) <= eps_0` for all `c in V_H` and all pair
   inequalities, so they refute ZE-TC_eps + arithmetic separation for that profile. No Diophantine hypothesis on `u*`
   is needed. At the barrier scale `W = O(n log n)` the same holds if `u*` has polynomial Diophantine quality (2.2)
   (contains nearstrata's `Gamma(u*) = O(log n)`). Conversely fakes have no other freedom: if quantitative ZE-TC holds
   on the fake pattern sets, fakes simply violate Step 0. **So the fake test for route 5 at `T = 1` is neither
   easier nor harder than ZE-TC on the barrier pattern sets.**
2. **Fakes on `Z_P` (Theorem 2.2, proved): regular tangent-cone points lift.** If `c -> Phi_c` is injective on
   `V(A,m)` (e.g. `m = rho`), `dim V <= d - 1`, and the forms have a regular real common zero whose cone germ is not
   in a hyperplane, then `Z_P` contains a real-analytic piece of nondegenerate points near `e`, and gap-admissible
   fakes lie **exactly on `Z_P`** (sublevel/Pyartli argument for the infinitely many characters).
3. **Barrier transfer realised in two fake models (computer-assisted, certified).** Skew double core, `T = 1`,
   `L = 2`, `m = rho = 8`:
   * `q = 7` (`n = 8`, 1 tangent form): a certified regular zero with all 28 gaps `>= 0.2496` (the maximum
     possible is `0.2857`), non-flat (Prop 3.1). Gap-admissible sdc7 fakes on `Z_P` exist (Cor 2.3). This settles the
     point R11 compute §6bis left open.
   * `q = 11` (`n = 12`, `dim V(A,8) = 5`, exact): Krawczyk-certified regular collision-free common zeros with gaps
     `>= 0.0100` and independent restricted Hessians (Prop 3.2); fakes on `Z_P` again (Cor 3.3).
   Both are underdetermined (`H_A(rho) < d`), LP-subcritical, with ad hoc `L`. In the regime `dim V_H >= |A|/2` and
   Theorem 2.2's mechanism is unavailable.
4. **Overdetermined orders behave as ZE-TC predicts (evidence; one exact case).**
   * pal7 at `L = 3`, `m = 11` (the round-5 Paley fake profile): its box *is* `M8'` up to a point permutation, so by
     R11 Thm 5.2 every real tangent-cone point and every point of `Bs` has a collision; **no gap-admissible pal7 fake
     lies on `Z_P`** (Prop 3.4, proved from refereed data).
   * sdc7 (`m = 7`, 43 forms), pal7 (`m = 7`, 95), sdc11 (`m = 7`, 2157), pal11 (`m = 8`, **436 forms, faithful**:
     `zeta <= 7 < 8 = rho`): the minimum of `kappa` over directions with all gaps `>= s` is positive for every `s > 0`
     tried and decays like `s^{4..8}`; minimisers sit next to cheap strata. For pal11, the round-5 fake profile at `n = 12` and the
     only faithful overdetermined fake model computed, the minimisers over `{sep >= s}` are again arithmetic
     progressions next to the 2-value stratum, with `kappa2* = 6e-11, 2e-8, 1e-7, 3e-6` at `s = 0.02, 0.05, 0.08, 0.12`.
5. **Growth of `dim V` in the skew-double-core family and shrinking of the cone (§6).** `dim V(A, rho) = 1, 5` for
   `q = 7, 11`, and `H_A(rho - 1) = 43, 2157`. The largest gap of a real tangent-cone point drops from `87%` of the
   maximum (`q = 7`) to `5.5-11%` (`q = 11`) at the top degree, and to `~0` once `H_A(m) >> d`. The cone shrinks
   *towards* the cheap strata but, at the underdetermined order, not onto them (collision-free certified points).
6. **Restricted ZE-TC_eps theorems (§5).** Explicit lower bounds exist from zonotopes (Thm 5.2), sections and
   products (Lemma 5.1), and transfer along Sylvester towers (Thm 5.4, conditional on a finite certificate not
   completed). **Prop 5.3:** zonotope lower bounds never apply in the regime, to fakes or to actual clusters (the
   R11 compute Prop 2.2 argument uses only the gap inequalities). No restricted theorem reaches an LP-supercritical box.

| task item | outcome | where |
|---|---|---|
| (1) `T = 1` boxes, `V(A,m)`, forms, real `TC` for barrier profiles | sdc (`n = 8, 12`), pal (`n = 8, 12`), `M8`, `M8'`: Hilbert functions, `zeta` bounds, forms; certified `TC` points | §3, §6 |
| (2) adversarial search | `kappa` landscapes over gap-constrained directions; fakes with `kappa = 0` exist exactly when `TC` has nondegenerate points (Thm 2.1); found in the underdetermined sdc models, not in any overdetermined model | §2, §4 |
| (2') gap-admissible fake ON `Z_P` | **exists** for sdc7, sdc11 models (`m = rho`); **impossible** for pal7 at `L = 3` | Cor 2.3, 3.3; Prop 3.4 |
| (3) restricted ZE-TC_eps theorems | zonotopal (explicit, blind in regime), section/product transfer, Sylvester tower (conditional) | §5 |
| (4) sdc `dim V` growth, cone shrinking | `1 -> 5` at top degree, `43 -> 2157` one below; cone shrinks, not onto cheap spans at the top degree | §6 |
| growth / uniform bound | **not proved** | — |

**Fake-circle test (mandatory).** The input fakes violate is still [D] Thm 4.1 (rigidity), entering through Step 0.
Theorem 2.1 shows that any *qualitative* failure of ZE-TC on a fake pattern set transfers the fake-circle barrier to
route 5, and Theorem 2.2 that in underdetermined systems it even transfers to the nonvanishing step itself (fakes on
`Z_P`); both phenomena are realised in models (sdc7, sdc11) and excluded in one overdetermined fake model (pal7,
`L = 3`).

## 1. Setting, notation, conventions

**Labels.** *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted with an exact source.
*Exact*: integer or exact-rational computation (python-flint, sympy). *Mod-p*: ranks modulo a prime `p ~ 2^58`
(two primes each time); a mod-`p` rank is `<=` the rational rank. *Certified*: exact rational arithmetic plus an
explicit error bound. *Evidence*: floating point or a finite search that no proof uses. Scripts:
`fake_checks/` (`x.py` writes `x_out.txt`). [D] = `round10/determinant.md`; R11 = `round11/README.md`;
`exact.md`, `construct.md` are `round5/`. The companion round-12 write-ups `zetc.md`, `t1rigour.md`,
`nearstrata.md` (same folder) are unrefereed; I use nothing from them except where re-proved or re-computed here.

**T = 1 slab-boxes.** `n` points, base point = row `0`, `d = n - 1`, `k' = (k_1, ..., k_d)`, `k_0 = -sum k'`. A
sign column with pattern `S` (rows with entry `-1`) gives `<k, beta> = +-2 k(P)`, `P = S` or its complement,
whichever omits the base row. For a set of patterns and a half-width `L`:

```text
A = A(profile, L) = { k' in {-1,0,1}^d : |k(P)| <= L for every pattern P }.
```

In the regime `L = floor(w/2)`, `w = sqrt(8 d ln 4r)` (R11 §5.3); in the computable models `L` is chosen ad hoc
(slabs do not cut at regime `w` for `n <= 12`). Repeated columns (copies) give the same slab, so `A` depends on the
*set* of patterns and on `L` only; the multiplicities enter only through `r` in `w`.
`V(A,m)`, `V_H`, `Phi_c(u) = sum_k c_k <k',u>^m`, `TC_H(A,m) = {u in R^d : Phi_c(u) = 0 for all c in V_H}` are as in
R11 §1. `H_A(j) = dim (S/J)_j` is the Hilbert function of the leading-form ideal `J = gr I(A)`; `rho(A)` its top
degree. The span of the tangent forms `{Phi_c : c in V(A,m)}` has dimension `H_A(m)` (the degree-`m` part of the
inverse system; R11 §2 R2), while `dim V(A,m) = sum_{j >= m} H_A(j)`. For `|u|_inf = 1`,
`kappa_c(u) = |Phi_c(u)|/(|c|_1 d^m)` (here `D/2 = d`), `kappa(u) = sup_{c in V_H} kappa_c(u)`.

* *Two-norm proxy.* `kappa2(u) := ||(I - P_{<m}) y_u||_2 / d^m`, `y_u = (<k',u>^m)_{k in A}`, `P_{<m}` the orthogonal
  projection onto `Poly_{<m}|_A`. Since `d^m kappa(u)` is the best uniform approximation error of `y_u` by
  `Poly_{<m}|_A` (when `V_H = V(A,m)`; nearstrata Prop 2.1, an LP duality that I re-derive in §4.1), one has
  `kappa2(u)/sqrt|A| <= kappa(u) <= kappa2(u)`.
* *Separation.* For `u in R^d` put `z = (0, u)` (all `n` positions) and `sep(u) = min_{x<y} |z_x - z_y|`; for
  `|u|_inf = 1` the largest possible value is `2/(n-1)` (equally spaced points), attained only up to reflection.
  `Gamma(u) = (2/binom(n,2)) sum_{x<y} ln^+ (1/|z_x - z_y|)` is nearstrata's pair defect.

**Fakes (angle models).** As in `exact.md` §1: distinct primes `p_j = 1 mod 4` (`j <= r`), `e_j = 1`, rows
`a_x in {0,1}^r`, positions `delta in [0, Delta]^n`, `Delta = C e^{-W/4}`, `W = sum_j log p_j`, windings
`k in Z^n`. The *pattern data* are the set `T` of patterns, the multiplicities `mu_S >= 1` (number of primes
carrying pattern `S`), and `k`. `L_T = {lambda in Q^n : sum lambda = 0, lambda(S) in Z for all S in T}`;
`lambda` is *relevant* if `lambda.k in Z`. With weights `omega_S = mu_S`:

```text
sigma_lambda(mu) = sum_{S in T} mu_S (|lambda(S)| - 1/2)          (exact.md Def. 2.1, normalised slack).
```

A profile is *gap-admissible* (an angle model) iff the two inequalities of `exact.md` Prop 1.1 hold for every
`lambda in L_T \ {0}`; for that proposition one needs (AI) (rows affinely independent) and
`prod_j (1 + 2/(sqrt p_j - 1)) < 1.75`. Gap-admissibility is equivalent to all pair/character inequalities
(`outside.md` Lemma 0.1); it is the defining property of the fake circles.

**Margin hypothesis (Mar).** (AI) holds on `T`; there is `rho_1 > 0` with `sigma_lambda(mu) >= rho_1 ||lambda||_1`
for every relevant `lambda in L_T \ {0}`; and every irrelevant `lambda` has `||lambda.k|| >= 1/den_T` (distance to the
nearest integer; `den_T` = exponent of `L_T/Z^n_0`).

**The barrier pattern data** (all satisfy (Mar)):
* **sdc_q** (skew double core, `exact.md` Thm 4.1): `n = q + 1`, `q = 3 mod 4` prime, `2q` balanced cuts
  `S_a, S'_a = S_a xor {0,a}`, `mu = 1`, `k = e_1`. Margin `sigma_c >= ||c||_1/2` on integral `c`; the
  fractional characters (`L_T/Z^n_0 = Z/(n/2)`) are irrelevant with `||lambda.k|| >= 2/n`.
* **pal_q** (round-5 Paley-I fake profile `P(H_q, id, b, a*)`, `construct.md` §2, `b = 8`, `k = 0`): `q` labels with
  `b` copies each plus one twin (flipped in one row) per row; `2q + 1` distinct patterns; all rational characters
  are integral, and `construct.md` Thm 1 gives `G(c) >= 4` on pairs and `G(c) >= 2||c||_1` on all other
  characters. Since `v_j(c) = (S^T c)_j/2`, `sigma_c = G(c)/2`, so `sigma_c >= ||c||_1` for every `c`: `rho_1 = 1`.
* **syl_k** (Sylvester-Hadamard `n = 2^k`, the profile of R11 J1; `mu = 1`, `k = 0`): not a strictly admissible fake
  profile by itself (Plotkin-tight with zero pair margin, `sigma = 0` on pairs); it is included as the
  "actual-type" comparison profile of R11 §2 (models `M8`, `M8'`).

Pattern sets are built exactly as in `round11/compute_checks/model8.py` (`fake_checks/flib.py`); for `q = 7` the
boxes reproduce R11's `|A| = 757` (sdc7, `L = 2`) and 1037 / 1877 (`M8`, `M8'`).

## 2. What a fake can do: transfer theorems

The point of this section: for fakes, **positions are free** subject only to the gap inequalities, and the box
`A`, the order `m` and `V_H` depend only on the pattern set (and `r`). So the fake test of R11 §8.3 item 4 reduces
to properties of the real tangent cone (Theorem 2.1) and, in the underdetermined case, of `Z_P` itself
(Theorem 2.2). Everything in this section is proved.

### 2.1 Tangent-cone points are directions of fakes

**Theorem 2.1 (fake directions).** Let pattern data `(T, mu, k)` on `n` rows satisfy (Mar), let `r = sum mu_S`,
and let `u* in R^d` (relative coordinates, base row `0`) be **nondegenerate** (`<v, u*> != 0` for all
`0 != v in Z^d`), `|u*|_inf = 1`. Let `C in (0,1]` and `theta > 0`. Then there is `P_0` such that for every
`P >= P_0` and every assignment of `r` distinct primes `p_j = 1 mod 4` in `[P, P e^{1/(4r)}]` to the patterns
(`mu_S` primes for `S`), there is a gap-admissible angle model with these primes, rows and windings whose relative
direction `hat u = u/|u|_inf` (`u_x = delta_x - delta_0`) satisfies `|hat u - u*|_2 <= theta e^{-W/4}`.

Consequently, for every finite `A` in the ternary cube and every `c` supported on `A` with `c in V(A,m)`:
`kappa_c(hat u) <= kappa_c(u*) + 2 m theta e^{-W/4}`. **In particular, if `u* in TC_H(A,m)` then the fakes have
`sup_{c in V_H} kappa_c(hat u) <= 2 m theta e^{-W/4}`, which is below Step 0's `eps_0 = C d e^{-W/4}/(m+1)` once
`theta < C d/(2m(m+1))`.**

*Proof.* Write `L = L_T`, `ell_0 = min{||lambda||_1 : 0 != lambda in L} > 0`, and `N(ell) = #{lambda in L :
||lambda||_1 <= ell}`, which is `<= (2 den_T ell + 1)^n`. Extend `u*` by `u*_0 = 0` to `R^n` and let
`tilde u = (u* - min u*)/(max u* - min u*) in [0,1]^n`; for zero-sum `lambda`, `lambda.tilde u = lambda.u*/s*`,
`s* in [1, 2]` the spread.

*Weights.* With `w_j = log p_j in [tau, tau + 1/(4r)]`, `tau = log P`, and `v(lambda) = sum_x lambda_x a_x in Z^r`:
`w(v(lambda)) - W/2 = sum_j w_j (|v_j(lambda)| - 1/2) >= tau sigma_lambda(mu) - 1/8` (terms with `v_j = 0` lose at
most `1/(8r)` each), and `w(v(lambda)) >= tau c_T ||lambda||_1` for a constant `c_T > 0` (the map
`lambda -> v(lambda)` is linear and injective by (AI)).

*Positions.* Fix `eta in (0, 1/4]` (chosen below), draw `xi` uniformly from `[-eta, eta]^n`, fix `t in [Delta/2, Delta]`,
and put `delta = (t/(1 + 2 eta)) (tilde u + xi + eta 1) in [0, Delta]^n`. For zero-sum `lambda`,
`lambda.delta = (t/(1+2eta)) (lambda.tilde u + lambda.xi)`, and `lambda.xi` has a density bounded by
`1/(2 eta |lambda|_inf) <= n/(2 eta ||lambda||_1)`.

*Sufficient condition.* Both inequalities of `exact.md` Prop 1.1 follow from

```text
dist( lambda.delta/2 + (pi/4) lambda.k , (pi/4) Z )  >=  arcsin( e^{-w(v(lambda))/2} )      for all 0 != lambda in L,   (2.1)
```

since `(pi/2)Z` and `pi/4 + (pi/2)Z` lie in `(pi/4)Z`. We check (2.1) in three classes.

1. *Relevant `lambda` with `||lambda||_1 <= Lambda`* (`Lambda` fixed below). Then `(pi/4)lambda.k in (pi/4)Z`, and
   for `|lambda.delta| <= pi/4` the left side is `|lambda.delta|/2`. Nondegeneracy gives `g_lambda := |lambda.tilde u| > 0`.
   If `eta ||lambda||_1 <= g_lambda/2` (true for all these finitely many `lambda` once `eta Lambda <= g_Lambda/2`,
   `g_Lambda = min g_lambda`), then `|lambda.delta| >= (Delta/2)(g_lambda/2)/(3/2) = Delta g_lambda/6` for every `xi`.
   Since `sin x >= 0.97 x` on `[0, pi/8]` and `e^{-w/2} <= e^{-W/4} e^{-tau sigma_lambda/2 + 1/16}`, (2.1) holds as soon as
   `C g_lambda/12 >= 1.07 e^{-tau rho_1 ||lambda||_1/2}`, i.e. as soon as
   `tau rho_1 ||lambda||_1 >= 2 ln(13/(C g_lambda))`. These are conditions on `tau` only.
2. *Irrelevant `lambda` with `||lambda||_1 <= Lambda_1 := 1/(n Delta)`.* Then `|lambda.delta|/2 <= 1/(2n)` and
   `(pi/4)||lambda.k|| >= pi/(4 den_T)`, so the left side of (2.1) is `>= pi/(4 den_T) - 1/(2n) >= 1/(4 den_T)`
   (`den_T <= n/2` for sdc, `= 1` for pal). The right side is `<= 2 e^{-tau c_T ell_0/2}`; (2.1) holds once
   `tau c_T ell_0 >= 2 ln(8 den_T)`.
3. *All remaining `lambda`* (relevant with `||lambda||_1 > Lambda`, irrelevant with `||lambda||_1 > Lambda_1`). The bad
   set of (2.1) for `lambda` is `{lambda.delta in B_lambda}`, `B_lambda` a union of intervals of length
   `4 arcsin(e^{-w/2}) <= 8 e^{-w/2}` centred on a progression of step `pi/2`; at most `||lambda||_1 Delta + 1` of
   them meet the range `|lambda.delta| <= ||lambda||_1 Delta/2`. With the density bound,
   `Prob(bad_lambda) <= (||lambda||_1 Delta + 1) 8 e^{-w/2} (1+2eta) n/(t eta ||lambda||_1)
   <= 24 n (1 + 1/(Delta ||lambda||_1)) e^{-w/2}/eta`. For relevant `lambda`, `e^{-w/2} <= e^{-W/4} e^{-tau rho_1 ||lambda||_1/2 + 1/16}`
   and `Delta = C e^{-W/4}`, so `Prob <= (50 n/(eta C)) e^{-tau rho_1 ||lambda||_1/2}` (using `||lambda||_1 >= ell_0`,
   absorbed for `tau` large). For irrelevant `lambda` with `||lambda||_1 > 1/(n Delta)`,
   `Prob <= 48 n^2 e^{-tau c_T ||lambda||_1/2}/eta`.

*Choice of parameters.* Take `eta = theta' e^{-W/4}` with `theta' = theta/(4 sqrt n)`, and `Lambda = r/rho_1 + 1`.
Since `W <= r(tau + 1)`, `1/eta <= e^{r(tau+1)/4}/theta'`. The sum of the class-3 probabilities is at most

```text
(50 n/(theta' C)) e^{r(tau+1)/4} sum_{ell > Lambda} dN(ell) e^{-tau rho_1 ell/2}  +  (48 n^2/theta') e^{r(tau+1)/4} sum_{ell > 1/(n Delta)} dN(ell) e^{-tau c_T ell/2},
```

and since `rho_1 Lambda/2 > r/4 + rho_1/2` and `N(ell)` grows polynomially, the first sum is
`O(e^{-tau rho_1/2 + O(log tau)})`, while the second is doubly exponentially small (`1/(n Delta) = e^{W/4}/(nC)`).
So for `tau >= tau_0` the total is `< 1`. Classes 1 and 2 hold for `tau >= tau_0'` (with `eta Lambda <= g_Lambda/2`,
true once `theta' e^{-W/4} Lambda <= g_Lambda/2`). Also `prod_j (1 + 2/(sqrt p_j - 1)) <= e^{2r/(sqrt P - 1)} < 1.75`
for `P` large. So for `P >= P_0` some `xi` makes (2.1) hold for every `lambda`, and `exact.md` Prop 1.1 gives a
gap-admissible realisation.

*Direction.* `u_x = delta_x - delta_0 = (t/(1+2eta)) s*^{-1} (u*_x + s*(xi_x - xi_0))`, so
`|u/(t s*^{-1}/(1+2eta)) - u*|_inf <= 2 s* eta <= 4 eta` and, after normalising by the sup norm,
`|hat u - u*|_2 <= 2 sqrt(d) 4 eta <= theta e^{-W/4}`.

*kappa.* For `|y|_inf <= 1` on the segment `[u*, hat u]`, `|<k',y>| <= d`, and
`|grad Phi_c(y)|_2 <= m |c|_1 d^{m-1} sqrt d`, so `|Phi_c(hat u) - Phi_c(u*)| <= m |c|_1 d^{m-1} sqrt d |hat u - u*|_2`
and `kappa_c(hat u) <= kappa_c(u*) + m |hat u - u*|_2/sqrt d`. QED

**Corollary 2.1' (what is refuted).** Let `(A, m)` be the `T = 1` box of fake pattern data satisfying (Mar), and
suppose `TC_H(A,m)` has a nondegenerate real point `u*`. Cheap families have total width `< m`, so there are finitely
many cheap spans `L_F`, and `u*` (lying on no rational hyperplane) has distance `delta* > 0` from their union. ZE-TC_eps
asserts that every `hat u` with `kappa(hat u) <= eps` lies within `rho_F ~ (eps 2^m)^{1/(m-d_F)}` (polynomial factors
aside) of some `L_F`. With `eps = eps_0(W) -> 0` the radii `rho_F -> 0`, while the fakes of Theorem 2.1 have
`kappa(hat u) <= eps_0` and `dist(hat u, L_F) >= delta* - theta e^{-W/4}`. So ZE-TC_eps fails for this box once `W` is
large, **with gap-admissible witnesses**; and since the fakes satisfy all pair inequalities, every arithmetic
separation lemma that uses only pair inequalities and the Plotkin average (as in `nearstrata.md`, `t1rigour.md`)
holds for them too. Hence ZE-TC_eps is false for this box, the fakes satisfy every other input of the route
(Step 0's conclusion, pair inequalities, separation), and route 5 as formulated cannot exclude clusters with this
pattern set: the fake-circle barrier transfers to it.

**Remarks.**
1. *No Diophantine hypothesis is needed.* Class 3 is handled by the random perturbation `xi`; nondegeneracy of `u*`
   is used only for the finitely many characters of class 1.
2. *Scale.* Theorem 2.1 lets `P -> infinity`, so its fakes have `W/(n log n)` unbounded. The proof is quantitative:
   the conditions are `tau rho_1 ||lambda||_1 >= 2 ln(13/(C g_lambda))` and `theta' e^{-W/4} Lambda <= g_Lambda/2`
   (class 1), `tau >= O(log n)` (class 2), and the class-3 sum, which is `< 1` once
   `tau(r/4 + rho_1/2) >= n ln(4 n Lambda) + O(log(n/(theta' C)))`, i.e. for `tau = O(log n)` when `r >= n`
   (true for sdc, `r = 2n - 2`, and pal, `r = 9n - 8`); here `Lambda = r/rho_1 + 1`. If the tangent-cone point satisfies

   ```text
   |lambda.u*| >= n^{-a ||lambda||_1}      for all relevant 0 != lambda in L_T with ||lambda||_1 <= r/rho_1 + 1,          (2.2)
   ```

   then class 1 holds for `tau >= (2a ln n + ln(13 s*/C))/rho_1` and the condition on `eta` for
   `r tau/4 >= a Lambda ln n + O(log n)`, i.e. `tau = O((a/rho_1 + 1) log n)` suffices, i.e. **fakes at `W = O(r log n)`, the barrier scale, with
   `kappa(hat u) <= eps_0`.** For pairs (2.2) is polynomial separation, `|u*_x - u*_y| >= n^{-2a}`; in nearstrata's
   language (2.2) contains `Gamma(u*) <= 4a ln n`. So a sufficient condition for a barrier-scale fake is: *does `TC_H(A,m)`
   contain a nondegenerate point with polynomial Diophantine quality (2.2)?*
3. *Residues.* The same proof works in the class `A_1(Q_0)` of `construct.md` (residue axiom (R)): the requirement
   (5.1) there has the extra factor `L^{||lambda||_1/2}` inside the arcsin, which is absorbed by `tau rho_1 ||lambda||_1/2`
   once `tau rho_1 > ln L`.
4. *What Theorem 2.1 does not give:* a fake **on** `Z_P`. A direction in the tangent cone is the Step-0 conclusion;
   membership in `Z_P` needs higher-order terms (§2.2).

### 2.2 Fakes on Z_P from regular tangent-cone points (underdetermined systems)

Since `A = -A`, `V_H = V_H^+ (+) V_H^-` (even and odd measures), and on real `delta` the equations of `Z_P` are
`F_c(delta) := sum_k c_k cos<k,delta> = 0` (`c` even) and `sum_k c_k sin<k,delta> = 0` (`c` odd) (R11 compute R3). So a
real basis `c_1, ..., c_s` of `V_H` consisting of even and odd vectors gives exactly `s = dim V_H` real equations.

**Theorem 2.2 (lifting).** Let `V' = V(A,m)` (any subspace `V'` of `V(A,m)` with `V' = -V'` and `V' >= V_H` works the same
way; since `Z(V') subset Z_P` for `V' >= V_H`, membership in `Z(V(A,m))` implies membership in `Z_P`). Suppose
`c -> Phi_c` is injective on `V'` (this holds whenever `m = rho(A)`: a `c in V(A,rho)` with `Phi_c == 0` annihilates
`Poly_{<= rho}|_A = R^A`), let `c_1, ..., c_s` be a real basis of `V'` of even and odd vectors, and let `u_0 in R^d`,
`|u_0|_inf = 1`, be a common zero of `Phi_{c_1}, ..., Phi_{c_s}` such that

* (R) the differentials `dPhi_{c_i}(u_0)` are linearly independent (so `s <= d - 1`, as `u_0` itself is in their
  common kernel by Euler's identity), and
* (F) the germ at `u_0` of the cone `Sigma = {Phi_{c_i} = 0, i <= s}` (a real-analytic manifold of dimension
  `d - s >= 1` near `u_0`) is not contained in any linear hyperplane.

Then there are `t_0 > 0` and a real-analytic map `psi : [0, t_0) x B -> R^d` (`B` an open ball in `R^{d-s-1}`)
with `psi(0, .)` a chart of `Sigma` near `u_0` (inside `{|u|_inf = 1}`), such that `t psi(t, y) in Z(V') subset Z_P` for all
`0 < t < t_0`, `y in B`. If the pattern data satisfy (Mar), then for every `C in (0,1]` and every `P >= P_0` (prime
window as in Theorem 2.1) there are **gap-admissible angle models with `delta in Z(V') subset Z_P`** (exactly), with
`|hat u - psi(0,y)| = O(e^{-W/4})` for some `y in B`. Their positions are nondegenerate.

*Proof.* (1) *Lifting.* For `c` of parity `(-1)^m` (all `c in V'` have it, since otherwise `Phi_c == 0`), the
Taylor expansion of `F_c(t u)` has no terms of order `< m` and

```text
G_c(t, u) := (-1)^{floor(m/2)} t^{-m} F_c(t u) = Phi_c(u)/m! + t^2 R_c(t, u),
```

with `R_c` real-analytic (power series in `t^2` with coefficients `sum_k c_k <k,u>^{m+2j}`). The map
`(t, u) -> (G_{c_i}(t,u))_i` is analytic near `(0, u_0)`, and its `u`-differential at `(0,u_0)` restricted to the
hyperplane `{u : u_{x_0} = 1}` (`x_0` a coordinate with `|u_{0,x_0}| = 1`; WLOG `u_{0,x_0} = 1`) has rank `s`, because
the `dPhi_{c_i}(u_0)` are independent and vanish on `u_0` (Euler), which is transverse to that hyperplane. The
implicit function theorem gives `psi` analytic in `(t, y)`; `t psi(t,y)` is a zero of every `F_{c_i}`, hence of every
`F_c`, `c in V'`, i.e. lies in `Z(V') subset Z_P`.

(2) *A sublevel lemma.* Put `chi(t, rho, y) = rho psi(t rho, y)` for `rho in [1/2, 1]` (so `t chi = t rho psi(t rho,y)` is
again in `Z_P`); for `t = 0`, `chi(0, ., .)` parametrises an open piece `Sigma'` of the cone `Sigma`, of dimension
`d - s`. By (F) and analyticity, at every point of the compact closure of the parameter box `Q = [1/2,1] x B'`
(`B'` a smaller closed ball) the derivatives of all orders of `chi(0, .)` span `R^d` (if `theta` were orthogonal to all
of them at one point, `theta.chi(0,.)` would be constant on the connected box, and, `chi(0,.)` being homogeneous in
`rho`, zero: `Sigma'` would lie in `theta^perp`). By compactness there are `k` and `gamma > 0` with
`max_{1 <= |alpha| <= k} |d^alpha (theta.chi(0, .))| >= gamma` on `Q` for every unit `theta`, and by continuity of the
derivatives in `t` the same holds with `gamma/2` for `chi(t, .)`, `t <= t_1`. Pyartli's lemma (the sublevel estimate
for functions with a non-vanishing derivative of order `<= k`, uniform in the family; see Bernik-Kleinbock-Margulis,
IMRN 2001, Lemma 3.3, or Kleinbock-Margulis, Ann. Math. 1998, Prop 3.4 for the analytic case) then gives `K_0` with

```text
| { (rho, y) in Q : |theta.chi(t, rho, y) - b| < epsilon } |  <=  K_0 epsilon^{1/k}      for all unit theta, real b, t <= t_1.   (2.3)
```

(3) *Gap-admissibility.* Now repeat the proof of Theorem 2.1 with `delta = t chi(t, rho, y)` (translated by a constant
vector to lie in `[0, Delta]^n`, which does not change `lambda.delta`; `|t chi|_inf <= t` so the spread is `<= 2t`),
`t in [Delta/4, Delta/2]` fixed, and `(rho, y)` uniform on `Q` in place of `xi`. Consider the three classes.
Class 1 (relevant, `||lambda||_1 <= Lambda`): by (F) the analytic function `(rho,y) -> lambda.chi(0,rho,y)` is not
identically zero, so `|lambda.chi(0, .)| >= g_Lambda > 0` on a sub-box `Q_1 subset Q` of positive measure for all these
finitely many `lambda` simultaneously (the zero sets are closed and nowhere dense; take a small box avoiding them),
and `|lambda.chi(t, .) - lambda.chi(0, .)| = O(t)`. So class 1 holds on `Q_1` for `tau` large exactly as before.
Classes 2-3 are handled by (2.3) applied with `theta = lambda/|lambda|_2`: the bad set of `lambda` has measure
`<= K_0 (count of centres) (8 e^{-w/2}/(t |lambda|_2))^{1/k}`, and the sums over `lambda` are as in Theorem 2.1 with
`e^{-w/2}` replaced by `e^{-w/(2k)}` and without the factor `1/eta` (the parameter box is fixed): `1/t <= 4/Delta =
4e^{W/4}/C` enters only as `e^{W/(4k)}`, which cancels against `e^{-w/(2k)} <= e^{-W/(4k) - tau sigma_lambda/(2k) + 1/(32k)}`
(relevant `lambda`; the remaining factor `e^{-tau rho_1 ||lambda||_1/(2k)}` is summable and `-> 0` as `tau -> infinity`),
and the irrelevant `lambda` with `||lambda||_1 > 1/(n Delta)` contribute a doubly exponentially small amount. So for `P >= P_0` a set of `(rho, y)` of positive measure in `Q_1` gives
(2.1) for all `lambda`. (2.1) for all integral `lambda` implies `lambda.delta != 0`, i.e. nondegeneracy. QED

**Corollary 2.3 (barrier transfer in the sdc7 model; proved, computer-assisted).** For the skew double core at
`q = 7` (`n = 8`), `T = 1`, `L = 2`, `m = 8` (R11 compute.md §6bis: `dim V = 1`, `V_H = V`), there are gap-admissible
angle models with the sdc7 pattern data whose positions lie **exactly on `Z_P`**, at every large prime window.
*Proof.* `m = rho = 8` (exact, §3), so injectivity holds; `s = 1`. Proposition 3.1 (`sdc7_cert.py`, exact + certified) exhibits a
zero `u*` of `Phi` with `grad Phi(u*) != 0` (R) and nonzero second fundamental form (F: a hypersurface germ contained
in a hyperplane is an open piece of it and has zero second fundamental form). Apply Theorem 2.2. QED

So the "barrier-transfer scenario" of round10 README §8.5 / R11 §7 item 4 is **realised by genuine
gap-admissible fakes** in the underdetermined sdc7 model (R11 left open whether positions on that hypersurface
could be made gap-admissible; they can, at every large prime window). This is a model statement (`L = 2` is not
[D]'s `L`, and the model is LP-subcritical); it confirms R11's reading that ZE can only be an overdetermination
statement.

## 3. V(A,m), tangent forms and real tangent cones for the barrier profiles

### 3.1 The models (exact / mod-p data)

All boxes are `T = 1`; `L` is ad hoc (at `n <= 12` the regime `w` gives no cut). Hilbert functions: exact over `Q`
(`model.py`, `|A| <= 2500`) or mod `p` with an isotypic decomposition under a cyclic coordinate symmetry
(`cyc.py`: shifts `x -> x+1` of `F_q`; `cycgen.py`: multipliers `x -> a x`), two primes each, validated against the
exact values at `n = 8`. `zeta` bounds: rigorous LP dual certificates (`lpzeta_sym.py`, `lpzeta_full.py`) and
verified zonotopes found by search (`zono_search*.py`).

| model | `n` | `L` | `|A|` | `H_A(0..rho)` | `rho` | `zeta` | faithful `m` | `H_A(m)` forms / `dim V(A,m)` | status |
|---|---|---|---|---|---|---|---|---|---|
| sdc7 | 8 | 2 | 757 | 1 7 28 77 161 251 188 43 1 | 8 | 6..7 | 8 (7 if `zeta = 6`) | 1 / 1 (m=8); 43 / 44 (m=7) | exact |
| sdc7 | 8 | 3 | 1765 | ... 301 153 35 | 10 | `>= 10` (zonotope, R11) | none | | exact |
| pal7 | 8 | 2 | 861 | 1 7 28 77 161 251 241 95 | 7 | 6..7 | 7 iff `zeta = 6` | 95 / 95 | exact |
| pal7 = `M8'` | 8 | 3 | 1877 | ... 195 63 7 | 11 | `<= 10` (R11) | 11 | 7 / 7 | exact; §3.4 |
| syl8 = `M8` | 8 | 2 | 1037 | ... 287 175 42 | 8 | `<= 7` (R11) | 8 | 42 / 42 | exact |
| sdc11 | 12 | 2 | 9373 | 1 11 66 275 880 2243 3735 2157 5 | 8 | 6..7 | 8 (7 if `zeta = 6`) | 5 / 5 (m=8); 2157 / 2162 (m=7) | mod p (2 primes); `dim V(A,8) = 5` exact |
| pal11 | 12 | 2 | 12093 | 1 11 66 275 880 2253 4664 3507 436 | 8 | 6..7 | 8 | **436 / 436** | mod p, `Z_5` blocks |
| sdc11, pal11 | 12 | 1 | 133 | 1 11 66 55 | 3 | | | | exact |

(pal11: blocks `s = 0, 1, 2` of the multiplier group `x -> 3x` computed mod `p`; blocks `3, 4` are Galois
conjugate to block 1 over `Q` and taken equal; the totals then sum to `|A|`. Zonotopes of length 6 were found in
sdc7, pal7, sdc11 by search; the rigorous LP bounds are `zeta <= 7` (sdc7, pal7: all ternary functionals, 42 resp. 3
automorphisms), `zeta <= 22/3` (sdc11, pal11: slab-orbit functionals).)

Reading of the table.
* For sdc the top degree is thin: `H_A(rho) = 1` (sdc7), `5` (sdc11), so at `m = rho` these fake-profile systems are
  **underdetermined** (`H_A(rho) < d - 1`). This is not general: `M8'` has 7, `M8` 42, pal7 (`L = 2`) 95 and pal11
  **436** forms at the top degree, all `>= d`.
* One degree lower the systems are strongly overdetermined (`H_A(rho - 1) = 43` for `d = 7`, `2157` for `d = 11`),
  but there `m = rho - 1 <= zeta` may hold (the LP bound `zeta <= 7` is not attained by any zonotope found so far),
  in which case the tangent cone is settled by a zonotope (R11 compute Prop 2.1; all of it lies in finitely many
  rational hyperplanes).
* No model with `n <= 12` is LP-supercritical (consistent with `zetc.md` §5 for Hadamard profiles): all slabs have
  `h_A(beta) = 4L >= m`.

### 3.2 sdc7 at m = rho = 8: the single tangent octic has comfortable well-separated zeros (certified)

`dim V(A,8) = 1`, `c` even, `|c|_inf = 42`, `|c|_1 = 1008`, support 687 (recomputed exactly, `sdc7_tc.py`; equal to
R11's vector up to sign). `Phi = sum_k c_k <k',u>^8`.

**Proposition 3.1 (certified, `sdc7_cert.py`).** There is a real zero `u*` of `Phi` with `|u*|_inf = 1 +- 1e-11`,
`grad Phi(u*) != 0`, nonzero second fundamental form of `{Phi = 0}` at `u*`, and

```text
sep(u*) >= 0.249565       (all 28 pair gaps of the 8 positions (0, u*), normalised),
```

against the absolute maximum `2/7 = 0.285714` of `sep` over all directions. Positions (rounded):
`u* = (0.50069, -0.90211, 0.75044, 1, 0.24963, -0.40127, -0.65146)`.

*Certificate.* A rational point `u_r` (denominator `10^14`) and a rational direction `g`; `p(t) = Phi(u_r + t g)` is
computed exactly in `Q[t]`; `p(-h) < 0 < p(h)` at `h = 1e-11`, and Sturm sequences show that `p` has exactly one root and
`p'` none in `[-h, h]`; so `u* = u_r + t* g` is a simple zero with `grad Phi(u*).g = p'(t*) != 0`. The gradient and
Hessian at `u_r` are computed exactly; with explicit bounds for the second and third derivatives on the segment
(`|k'.x| <= 7.007`, `|c|_1 = 1008`), the restricted Hessian `P H P` (`P` the projection onto `grad^perp`) at `u*`
differs from its value at `u_r` by at most `1.0e2` in Frobenius norm, while `||P_r H_r P_r||_F = 8.30e4`.
(Numerically the separation optimum on `{Phi = 0}` is `0.25`, attained at a degenerate boundary point with many
gaps equal to `0.25`; `u*` is a generic nearby point.)

**Corollary 2.3 is now proved** (with (R), (F) from Proposition 3.1): gap-admissible sdc7 fakes on `Z_P` exist at every
large prime window. Since `Gamma(u*) = 0.981` (equally spaced points: `0.921`) and every gap is `>= 0.2495`, `u*` satisfies (2.2) for pairs with
`a = O(1)`; the separation demanded of a fake (pair slack `O(log n)`) is far weaker than what `u*` offers. In this
model the tangent cone **does not shrink towards the collision strata**: it contains points almost as well
separated as any configuration of 8 points can be.

### 3.3 sdc11 at m = rho = 8: five octics, certified collision-free zeros, but the cone has shrunk

`dim V(A,8) = 5`, all five basis vectors even and shift-invariant (exact, §6). Numerically (`sdc11_tc.py`,
`sdc11_m8_land.py`; evidence) the real common zero set of the five octics behaves very differently from sdc7:
SLSQP maximisation of the separation over the zero set (60 starts) ends at zeros clustered within `~3e-3` of
3-value strata (cheap spans of cost `<= 6`), and the minimum of `sum_i (Phi_i/(|c_i|_1 d^8))^2` over `{sep >= s}`
is `~1e-20..1e-16` (zeros) for `s <= 0.01` but positive (`4e-14, 1e-11, 4e-10, 5e-9, 3e-8` at
`s = 0.02, 0.05, 0.08, 0.12, 0.16`, against typical values `~3e-8`) above; Newton from the `s = 0.02, 0.05` minimisers
falls back to zeros with separation `7e-4`, `3e-3`.

**Proposition 3.2 (certified, `krawczyk.py`; outputs `krawczyk_sdc11_out.txt`, `krawczyk_sdc11_b_out.txt`).** The
five sdc11 tangent octics (`L = 2`, `m = 8`) have real common zeros `u*`, `u**` with `u_0 = 1 = |u|_inf`, regular
(Jacobian of rank 5), with nonzero, linearly independent restricted Hessians (condition (F) of Theorem 2.2), and

```text
sep(u*) >= 0.005287  (Gamma = 3.03),        sep(u**) >= 0.010008  (Gamma = 4.03)      [max possible sep: 2/11 = 0.1818].
```

*Certificate.* Six coordinates are fixed to rationals (denominator `10^25`), the remaining five are enclosed by a
Krawczyk test in arb ball arithmetic (400 bits) on a box of radius `1e-30` (`K(X) - x~` has width `< 1e-40`), which
proves existence and uniqueness of a zero in the box and invertibility of the `5 x 5` Jacobian on it. For (F): an
orthonormal basis `B` of the approximate tangent space satisfies `|J(u*) B| <= 3.3e-17` (arb bound), the `5 x 36`
matrix of the restricted Hessians `B^T H_i B` has smallest singular value `0.0716` (resp. `0.0106`), and the
perturbation bound for passing to the exact tangent space is `<= 8.1e-13` (resp. `1.8e-13`).

**Corollary 3.3.** Theorem 2.2 applies to sdc11 (`m = rho = 8`, injectivity, `s = 5 <= d - 1`, (R), (F)): **there are
gap-admissible sdc11 fakes on `Z_P`** at every large prime window. So the barrier transfer of Corollary 2.3 is
not an `n = 8` accident; it persists at `n = 12`, in a system whose real tangent cone has nevertheless shrunk from
"almost maximal separation" (sdc7: `0.2496` of `0.2857`) to "between `5.5%` (certified) and `~11%` (no zero found beyond) of the maximal separation" (sdc11).

### 3.4 The round-5 Paley fake profile at q = 7: the overdetermined order is toric (proved)

**Proposition 3.4.** For the Paley-I fake profile pal7 (round-5 profile `P(H_7, id, b, a*)`, all `b`) the `T = 1`
box at `L = 3` equals the box of the Paley core alone (the twin slabs do not cut; exact set equality,
`|A| = 1877`), and a permutation of the 8 points fixing the base point (`(0,1,2,4,3,6,7,5)`, found by search and
checked) maps the Paley-core patterns onto the Sylvester patterns, hence this box onto `M8'`. Consequently, by R11
compute Thm 5.2(a) (refereed) and J4: at `m = 11` (faithful, `zeta <= 10 < 11 = rho`, `dim V = 7`), every point of
`Bs(A, 11)` in the complex torus has a collision `w_x = 1` or `w_x = w_y`, and every real common zero of the tangent
forms lies on a collision hyperplane. **So no gap-admissible pal7 fake lies on `Z_P`, and every gap-admissible pal7
direction has `kappa(hat u) > 0`** (qualitative ZE-TC^braid for a genuine fake profile at an overdetermined order).
*Proof.* A point permutation fixing the base acts on `k'` by a coordinate permutation, maps slabs to slabs and
collision hyperplanes `u_x = 0`, `u_x = u_y` to collision hyperplanes. QED
(At `L = 2` the pal7 box has `rho = 7` with `H_A(7) = 95` and LP bound `zeta <= 7`, so `m = 7` may be zonotopal; §4.)

## 4. Adversarial search: gap-admissible fakes with small sup_c kappa_c

### 4.1 What is searched

By Theorem 2.1 a fake direction is any point of the gap-admissible region, and only its tangent-cone data matter.
At model size the natural proxy for gap-admissibility is the normalised separation `sep(u)` (pair inequalities)
and `Gamma(u)`; the quantity searched is

```text
kappa*(s) = min { kappa(u) : |u|_inf = 1, sep(u) >= s }        (s in [0, 2/d]),
```

with `kappa` replaced by the proxy `kappa2`. *Duality (re-derived):* for `y in R^A` and a subspace `V` with
annihilator `V^perp`, `sup_{0 != c in V} |<c,y>|/|c|_1 = min_{g in V^perp} |y - g|_inf` (`<=`: `<c,y> = <c,y-g>`;
`>=`: Hahn-Banach, the dual of `(R^A/V^perp, quotient sup-norm)` is `(V, l^1)`). With `V = V(A,m)`,
`V^perp = Poly_{<m}|_A`, `y = y_u`: `d^m kappa_V(u) = min_g |y_u - g|_inf`, which lies between
`|y_u - P y_u|_2/sqrt|A|` and `|y_u - P y_u|_2`. Since `V_H subset V`, `kappa_H <= kappa_V <= kappa2`, so the reported
`kappa2` values are **upper bounds for `kappa_H`** at the reported points (rigorous up to floating point), and the
minima over `s` are evidence only (local optimisation, SLSQP, random multistart; `tcnum.py`, and the
symmetry-reduced `tcsym.py`, which block-diagonalises the projection over a cyclic coordinate symmetry and agrees
with `tcnum.py` on sdc7).

### 4.2 Results

`kappa2*(s)` (upper bounds at the points found; the minimum itself is evidence):

| model | `m` | forms | `s = 0` | `0.02` | `0.05` | `0.08` | `0.10` | `0.15` | `0.20` | max `sep` on `TC` found |
|---|---|---|---|---|---|---|---|---|---|---|
| sdc7 | 8 | 1 | `2e-16` | | | | `5e-16` | | `6e-16` | **0.2496 (certified, Prop 3.1)** |
| sdc7 | 7 | 43 | `1e-15` (collision) | | `8.3e-8` | | `2.1e-6` | `1.8e-5` | `4.6e-5` | 0 |
| pal7 | 7 | 95 | `1e-15` (collision) | | `1.3e-7` | | `3.7e-6` | | | 0 |
| sdc11 | 8 | 5 | `1e-20` (zeros up to `s = 0.01`) | `4e-14`\* | `1e-11`\* | `4e-10`\* | | `5e-9`\* (`s=.12`) | `3e-8`\* (`s=.16`) | **0.0100 (certified, Prop 3.2)**; none found `>= 0.02` |
| sdc11 | 7 | 2157 | `6e-12` (sep `1e-4`) | `6.5e-9` | `3.9e-7` | `6.6e-6` | | | | 0 |
| pal11 | 8 | 436 | `2e-15` (collision) | `6.0e-11` | `1.9e-8` | `1.1e-7` | | `2.8e-6` (`s=.12`) | `3.0e-5` (`s=.16`) | 0 |

(\*: sdc11 at `m = 8` reports `(sum_i (Phi_i/(|c_i|_1 d^8))^2)^{1/2}` over the 5 basis forms, typical size `3e-8`.)
Separation scale: `2/d` is the largest possible `sep` (`0.2857` at `n = 8`, `0.1818` at `n = 12`).

* **Underdetermined orders (sdc at `m = rho`).** At `n = 8` (one form) the tangent cone contains points that are
  almost as well separated as possible; at `n = 12` (five forms) it still contains certified collision-free points,
  but only up to `sep ~ 0.01-0.02`, i.e. `5.5-11%` of the maximum. These are honest fake directions (Theorem 2.1), and by Theorem 2.2 even honest fakes on
  `Z_P`. No restricted ZE-TC_eps theorem can hold there.
* **Overdetermined orders.** In every overdetermined model (`H_A(m) >= 6d`) the minimisers over `{sep >= s}` are
  arithmetic progressions of step exactly `s` together with one far point (e.g. sdc11, `m = 7`, `s = 0.05`: positions
  `-0.1, -0.05, 0, 0.05, ..., 0.4, 1`): the optimiser runs into the boundary of the separated region next to the
  cheap stratum "all but one point coincide with the base point" (R11 Prop 5.1 type; its span is `R e_x`, cost `2 < m`).
  `kappa2*` decays roughly like `s^{4..8}` as `s -> 0`, and no zero with `sep > 1e-4` was found. This is the
  behaviour ZE-TC^braid with a polynomial Lojasiewicz exponent predicts, at model size. In the regime the fake
  separation is `n^{-O(1)}` and the required bound is `kappa > eps_0 ~ e^{-W/4}`; a decay `s^{O(m)}` would be ample
  (R11 §5.3), but nothing here controls the exponent's growth with `n`.
* *Caveat* for sdc7, pal7, sdc11 at `m = 7`: `zeta = 7` is not excluded (§3.1); if a zonotope of length 7 exists, the
  collapse there is explained by one character (Lemma 5.1(1)). pal11 at `m = 8` is faithful (`zeta <= 7`).

## 5. Restricted ZE-TC_eps theorems (lower bounds for sup_c kappa_c)

Throughout, `|u|_inf = 1`, `d^m kappa(u) >= sup_{c in X} |Phi_c(u)|/|c|_1` for any subset `X` of `V_H`.

**Lemma 5.1 (three explicit lower-bound mechanisms; proved).**
1. *Zonotopes.* If `k_0 + sum_i [0, L_i] v_i` lies in `A` and `0 <= alpha_i <= L_i`, `sum alpha_i = m`, then the
   finite-difference vector `c` of `w^{k_0} prod_i (w^{v_i} - 1)^{alpha_i}` lies in `V(A,m)`, `|c|_1 <= 2^m`
   (so `c in V_H` once `H >= 2^m`), and `Phi_c(u) = m! prod_i <v_i, u>^{alpha_i}`. Hence
   `kappa(u) >= (m!/(2d)^m) prod_i |<v_i,u>|^{alpha_i}`.
2. *Sections.* If `V in Z^{d x j}` has rank `j`, `k_0 in Z^d` and `A' = {t in Z^j : k_0 + V t in A}`, then every `c' in V(A', m)` pushes forward to `c = psi_* c' in V(A,m)`
   with `|c|_1 = |c'|_1` and `Phi_c(u) = Phi'_{c'}(V^T u)`. Hence, writing `v = V^T u` and `d' = max_{t in A'} |t|_1`,
   `kappa_A(u) >= (d'/d)^m |v|_inf^m kappa_{A'}(v/|v|_inf)`.
3. *Products.* If `A = A_1 x A_2` (disjoint coordinates) and `m = m_1 + m_2`, then `c_1 (x) c_2 in V(A, m)` for
   `c_i in V(A_i, m_i)`, `|c_1 (x) c_2|_1 = |c_1|_1 |c_2|_1` and `Phi_{c_1 (x) c_2}(u_1,u_2) = binom(m, m_1)
   Phi_{c_1}(u_1) Phi_{c_2}(u_2)`.

*Proof.* (1) is R11 compute Prop 2.1 plus the observation that the degree-`m` part of
`prod_i (e^{<v_i,z>} - 1)^{alpha_i}` is `prod_i <v_i,z>^{alpha_i}`, while the degree-`m` part of `sum_k c_k e^{<k,z>}`
is `Phi_c(z)/m!`. (2) `sum_t c'_t (<k_0,u> + <t, v>)^m`: the terms with `<k_0,u>^{i}`, `i >= 1`, are moments of order
`< m` of `c'`, which vanish; moments of order `< m` of `c` are moments of `c'` of polynomials of degree `< m` in `t`. The
normalisation: `|<t,v>| <= d' |v|_inf` on `A'`. (3) Expand `(<k,u_1> + <l,u_2>)^m`; `c_1` kills `<k,u_1>^i` for
`i < m_1`, `c_2` kills `<l,u_2>^{m-i}` for `i > m_1`. QED

**Theorem 5.2 (zonotopal profiles: an explicit ZE-TC_eps; proved).** Suppose `A` contains a zonotope with
generators `v_1, ..., v_s`, lengths `L_i`, `sum L_i >= m`. Then every `u` with `|u|_inf = 1` satisfies
`kappa(u) >= (m!/(2d)^m) min_i |<v_i, u>|^m`, and in particular every gap-admissible fake (or actual) direction
satisfies `kappa(hat u) >= (m!/(2d)^m) (|u|_inf^{-1} 2 e^{-X_max/2})^m`, `X_max = max_i X_{v_i}`
(`X_v = sum_p log p |<v, a_p>|`), by the gap lemma `|<v,delta>| >= 2 e^{-X_v/2}` (`exact.md` Prop 1.1 (E)/(O) form).

**Proposition 5.3 (Theorem 5.2 never applies in the regime, for fakes or actual clusters; proved).** Let a
gap-admissible angle model (or an actual cluster) satisfy `W >= 16 ln n + 8 log^+ C + 14` and let `(A, m)` be its own
`T = 1` slab-box at any order `m >= 4w` (`h_A(beta_p) <= 2w <= m/2`). Then `zeta(A) < m`.
*Proof.* R11 compute Prop 2.2 verbatim: a zonotope with `sum L_i >= m` forces some `X_{v_i} <= W/4`, while the gap
inequality (which fakes satisfy by definition) gives `X_v >= W/2 - 2 log(C |v|_1/(2 sqrt 2))` for every `v` with
`<v,delta> != 0`, and `|v_i|_1 <= 2D` with `D = 2(n-1)` at `T = 1`. QED

So zonotope (character) lower bounds, the only mechanism with explicit constants known in general, are exactly as
blind to fakes as to actual clusters: they exist only for boxes that no gap-admissible configuration in the regime
can produce. The section and product mechanisms of Lemma 5.1 transport lower bounds from small models.

**Theorem 5.4 (Sylvester towers: restricted ZE-TC_eps from a section certificate; proved conditional on a
finite certificate, see below).** Let `L in {2, 3}`, `(m, A') = (8, M8)` resp. `(11, M8')`, and suppose
`kappa_{A'}(v) >= f(sigma)` for every `v in R^7`, `|v|_inf = 1`, whose 8 positions `(0, v)` have all pair gaps
`>= sigma` and, for `L = 2`, avoid the 27 tangent hyperplanes of R11 Thm 5.2(b) by `>= sigma`. Then for every `k >= 3`
and every `u` in the Sylvester box `Syl(2^k, L)` coordinates with `|u|_inf = 1`: for every subgroup `U <= F_2^k`
of order 8 containing a coordinate `x` with `|u_x| = 1`,

```text
kappa_{Syl(2^k,L)}(u)  >=  (7/(2^k - 1))^m  f( sep(u|_U) ).
```

*Proof.* zetc.md Lemma 4.1 (the section `{supp k subset U}` of `Syl(2^k,L)` is `Syl(8,L)`; re-checked: `hat k(a)`
depends only on `a|_U`) and Lemma 5.1(2) with `|u|_U|_inf = 1`. QED

Only `f(sigma) > 0` for all `sigma > 0` is known (from R11 Thm 5.2 for `M8'`, J4); an explicit `f` needs either a
Nullstellensatz certificate for the cubic/quartic factors of the `M8'` forms (`(prod collisions)^N in (G_1..G_7)` with
explicit cofactors) or a branch-and-bound over the separated region, and **neither was completed**: a branch-and-bound
with interval Taylor forms needs boxes of side `~sigma` near the boundary of the separated region, i.e. `~sigma^{-6}`
boxes, which is out of reach for the separations a fake has (`sigma ~ n^{-2}`). Numerical values of `f` are in §4
(evidence). Theorem 5.4 is in any case LP-subcritical (LP-supercritical needs `m >= 8L`; here `8 < 16`, `11 < 24`).

## 6. The skew double core family: dim V growth at T = 1 and the real tangent cone

Computable members: `q = 7` (`n = 8`) and `q = 11` (`n = 12`); `q = 19` has `3^19` cube points and is out of reach
(the box enumeration is feasible, but the Hilbert function is not, even with the `Z_19` symmetry). At `L = 1` both
are trivial (`rho = 3`). At `L = 3`, `q = 7` contains a zonotope of length 10 = `rho` (R11, `sdc7_L3_zeta.py`), so it is
unfaithful; `q = 11`, `L = 3` (`|A| = 71633`) was not computed. So the comparable data are at `L = 2`:

| `q` | `n` | `|A|` | `H_A(0..rho)` | `rho` | `zeta` | `dim V(A,8) = H(8)` | `H(7)` / `dim V(A,7)` |
|---|---|---|---|---|---|---|---|
| 7 | 8 | 757 | 1 7 28 77 161 251 188 43 1 | 8 | 6..7 | **1** | 43 / 44 |
| 11 | 12 | 9373 | 1 11 66 275 880 2243 3735 2157 5 | 8 | 6..7 | **5** | 2157 / 2162 |

(`q = 7` exact over `Q`; `q = 11` mod two primes, with `dim_Q V(A,8) = 5` proved exactly by the invariant-subspace
computation `sdc11_m8.py`; `zeta`: verified zonotope of length 6 and rigorous LP bound 7.)

**Answers.**
1. *Growth of `dim V`.* At the top degree `m = rho = 8` the system stays tiny and underdetermined: `1` form in
   `d = 7` variables, `5` forms in `d = 11` variables. One degree lower it is strongly overdetermined and grows fast
   (`43 -> 2157`, i.e. `~6 d -> ~196 d`). For sdc the top degree is thin, but this is not a general feature (the Paley twin profile pal11 has `H_A(rho) = 436 > d` at its top degree `rho = 8`, §4), so "`dim V` small" is a property of the order `m = rho` of this family, not of `T = 1` boxes; in the regime
   `m = ceil(4w) << m_max ~ 0.64 n`, deep below the top degree.
2. *Does the real tangent cone shrink to cheap spans once `dim V` is large?* Measured by the largest separation of a
   real tangent-cone point (normalised by the maximum `2/d`):

   | system | forms / `d` | max `sep` on `TC` / `(2/d)` | Theorem 2.2 (fakes on `Z_P`) |
   |---|---|---|---|
   | sdc7, `m = 8` | 1 / 7 | `>= 0.873` (certified) | applies |
   | sdc11, `m = 8` | 5 / 11 | `>= 0.055` (certified), `< ~0.11` (evidence) | applies |
   | sdc7, `m = 7` | 43 / 7 | `~0` (no zero with `sep > 0` found) | `H > d`: not applicable |
   | sdc11, `m = 7` | 2157 / 11 | `~0` (`< 1e-3`; evidence) | not applicable |

   So the cone shrinks quickly with the number of forms, and once `H_A(m) >> d` the numerical tangent cone sits on
   the collision/cheap strata: the near-minimisers of `kappa2` over `{sep >= s}` are arithmetic progressions next to
   the stratum "all but one point at the base point" (span `R e_x`, cost 2), with `kappa2*(s) ~ s^{4..8}`.
   **But "shrinks to cheap spans" is not what happens at the underdetermined order**: at sdc11, `m = 8`, the
   certified zeros are collision-free, with `Gamma = 3-4.5`, so qualitative ZE-TC^deg fails there, while the
   separation available on the cone has dropped by an order of magnitude. Whether the overdetermined orders contain
   *any* nondegenerate real point (ZE-TC^deg) is not decided by these computations; caveat `zeta = 7` (§3.1).

## 7. Consequences for route 5 and open points

1. **The fake test of the `T = 1` programme is exactly the tangent-cone hypothesis on fake pattern sets.**
   Theorem 2.1 shows that fakes have no freedom beyond choosing a direction: every nondegenerate point of
   `TC_H(A,m)` of a fake profile is the direction (to within `e^{-W/4}`) of genuine gap-admissible fakes satisfying
   Step 0's conclusion; conversely a fake direction is an arbitrary point of the gap-admissible region. So:
   * route 5 as formulated (Step 0 + ZE-TC_eps + arithmetic separation) is **capped** as soon as the `T = 1` box of
     *one* strictly admissible fake pattern set in the regime has a nondegenerate real tangent-cone point (at large
     `W`), or one with the polynomial Diophantine quality (2.2) (at the barrier scale `W = O(n log n)`);
   * conversely, if quantitative ZE-TC (nearstrata's MIN\*, or ZE-TC_eps relative to rational strata) holds on all
     profiles including the fake ones, the fakes are harmless: they then simply violate Step 0, which only actual
     clusters must satisfy.
   The adversarial fake search therefore cannot be easier or harder than deciding ZE-TC on the barrier pattern sets
   (Paley twin profile, skew double core); this is now a theorem, not a heuristic.
2. **Barrier transfer is real but, so far, only in underdetermined systems.** In the fake-profile models at the
   top degree (sdc7, sdc11; Props 3.1-3.2) the real tangent cone has regular collision-free points (very well
   separated at `n = 8`, moderately at `n = 12`), and Theorem 2.2 places genuine gap-admissible fakes **on `Z_P`**.
   Both are models with `H_A(m) < d`, at the top degree, LP-subcritical and at ad hoc `L`. In the regime, `m = ceil(4w)`
   is far below the top degree `m_max ~ 0.64 n` and `dim V_H >= |A|/2`, so the transversality mechanism of
   Theorem 2.2 is unavailable (a regular common zero needs `dim V_H <= d - 1`).
3. **Overdetermined models behave as ZE-TC predicts (evidence).** One degree below the top, the same fake profiles
   have `43` (`n = 8`) and `2157` (`n = 12`) tangent forms, and the numerical tangent-cone landscape (§4) collapses onto
   the collision strata: the minimum of `kappa` over directions with all gaps `>= s` is positive for every `s > 0`
   tried and decays like a power of `s`. The caveat is that at these orders `zeta = m` is not excluded (LP bound
   `zeta <= 7`), in which case the collapse is explained by a single character. The round-5 Paley fake profile at
   `n = 12` (pal11, `L = 2`) has a **faithful** overdetermined top degree (`zeta <= 7 < 8 = rho`, 436 octics in 11
   variables) and shows the same collapse, with no caveat; at `n = 8` (pal7, `L = 3` = `M8'`) the collapse is a theorem
   (Prop 3.4).
4. **Restricted positive theorems** (§5): explicit lower bounds come only from zonotopes (blind in the regime, for
   fakes and actual clusters alike, Prop 5.3), and from sections and products of small certified models (Sylvester
   towers, LP-subcritical). An explicit `f(sigma)` for `M8'` is the first missing finite computation.

**Open points (in order of usefulness).**
* Decide whether `TC_H(A, ceil(4w))` of a regime-size fake pattern set has nondegenerate points. Theorem 2.1 makes
  this the whole fake test. The computable proxies are overdetermined faithful fake-profile models: pal11 at
  `L = 2`, `m = 8` (436 octics) is the first one; deciding exactly whether its real tangent cone lies in the braid
  arrangement (as for `M8'`, by the R11 cone-reduction/saturation method) is the most useful next computation, then
  sdc11 at `L = 3` (`|A| = 71633`).
* An explicit Lojasiewicz function `f(sigma)` for one overdetermined zonotope-free model (Theorem 5.4's input), via a
  Nullstellensatz certificate.
* Whether the underdetermined top-degree phenomenon persists: `H_A(rho) = 1, 5` for sdc at `n = 8, 12`; if
  `H_A(rho) < d` for all sdc orders, Theorem 2.2 yields fakes on `Z_P` for a whole family of models at `m = rho`.

## Appendix: scripts and outputs (`fake_checks/`)

| script | content | status |
|---|---|---|
| `flib.py` | pattern sets (sdc, pal, Sylvester; as in R11 `model8.py`), box enumeration, reduced monomials | library |
| `sizes.py` | `|A|` for all models (`sizes_out.txt`) | exact |
| `model.py` | Hilbert functions (`model_small_out.txt`), exact `V(A,m)` bases | exact (`|A| <= 2500`), mod `p` otherwise |
| `cyc.py`, `cycgen.py` | isotypic Hilbert functions under `x -> x+1` / `x -> a x` (`cyc_sdc11_L2_out.txt`, `cycgen_pal11_L2_out.txt`) | mod `p`, two primes / Galois |
| `lpzeta_sym.py`, `lpzeta_sym2.py`, `lpzeta_full.py` | rigorous LP dual bounds on `zeta` (`lpzeta_*_out.txt`) | exact rational certificate |
| `zono_search.py`, `zono_search2.py` | verified zonotopes (length 6) | exact verification |
| `sdc7_tc.py`, `sdc7_cert.py` | sdc7 form, separation optimum, Prop 3.1 certificate (`sdc7_cert_out.txt`) | exact + certified |
| `sdc11_m8.py` | exact `V(A,8)` (dimension 5) for sdc11 via the invariant subspace, verified on all 17469 monomials | exact |
| `sdc11_tc.py`, `sdc11_m8_land.py`, `sdc11_signs.py`, `sdc11_coll.py`, `sdc11_psd.py` | sdc11 octics: zero search, landscape, sign structure, no collision factor, no semidefinite combination | evidence |
| `krawczyk.py` | Prop 3.2: Krawczyk certificate, separation bound, condition (F) (`krawczyk_sdc11_out.txt`, `krawczyk_sdc11_b_out.txt`) | certified (arb) |
| `tcnum.py`, `tcsym.py` | `kappa2` landscape over `{sep >= s}` (`tcnum_*_out.txt`, `tcsym_*_out.txt`) | evidence |

Working notes: `fake_notes.md`.
