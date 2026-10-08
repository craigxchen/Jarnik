# Referee report on `growth/conditional.md` (Vojta ⟹ uniform `C sqrt(R)` bound)

Referee: adversarial audit, 2026-10-08. All checks are in
`scratchpad/growth/referee_conditional_checks/`. The author's scripts were re-run in its `rerun/`
subfolder.

## Verdict

**The main result survives as stated.** Theorem 5.1 is a correct conditional theorem. It
assumes one instance of Vojta's Main Conjecture: rational points, `D = 0`, on a fixed smooth
projective model `X^#_7` of `Bl_Y X_7`. It concludes that there is an absolute `M_0` such that
every arc of length `<= C sqrt(R)` has at most `M_0` lattice points once `R > c_0^2 C^42`.

Theorem 5.2 (all `alpha < 1`), Corollary 5.3 (divisors near `sqrt n`), Proposition 3.2
(non-density equivalence), Propositions 6.1–6.2 (Ru–Vojta with `S = {infinity}`) and §8 (`abc`)
are also correct.

* I found no fatal or major error.
* The issues are minor: a few overstated sentences about the `beta` numerics, a loose
  converse in Proposition 3.2(ii), a corollary given only as a sketch (its details check out),
  and one summary sentence on Bombieri–Lang.
* Novelty: the result is new relative to the project notes. In the literature I found no prior
  derivation, but the method (Vojta on a blow-up gives approximation exponents) is standard.

## 1. Line-by-line audit of the load-bearing chain

**Lemma 2.1 (geometry). Correct.**

* The affine cone `{u_1v_1 = ... = u_Mv_M}` has a dense part `{t != 0}` isomorphic to
  `G_m^(M+1)`. The part `{t = 0}` is a union of `2^M` linear spaces of dimension `M`.
  Every component of a codimension `M-1` complete intersection has dimension `>= M+1`, so the
  cone is irreducible. It is reduced because it is Cohen–Macaulay and generically smooth.
* `omega = O(-2M + 2(M-1)) = O(-2)`.
* The Jacobian drops rank exactly when `>= 2` pairs `(u_j, v_j)` vanish. That locus has
  codimension 3, so `X_M` is normal by Serre's criterion.
* The torus `{lambda_j mu_j = const}` acts with a single orbit on `{all u_j v_j != 0}`.
* Normal Gorenstein toric implies canonical, since discrepancies are integers `> -1`.
* On `Y` every pair is `(x+iy, x-iy) != 0`, so `Y` lies in the smooth locus.
* Hence `K_Bl = pi^*K + (codim Y - 1) E = -2H + (M-2)E`, and (2.1) holds with `a_i >= 0`.
* `phi` is an isomorphism over the smooth locus of `Bl`. That locus contains `E` and the whole
  torus, so short-arc tuples lift uniquely and avoid `E` and all the `F_i`.

**Lemma 3.1 and Proposition 3.2(i) (Bezout). Correct.**

* `G` has degree `<= D` in one point `w` and vanishes on `>= 2D+1` points of the irreducible
  conic `Q_N` (`N != 0`). Bezout then forces `G` to vanish on `Q_N`.
* Homogeneity of `F` moves the vanishing to every `Q_{lambda^2 N}`. Their union is the dense
  open set `{t != 0}` of the irreducible cone.
* This is the `M`-variable version of the notes' "grid lemma" (`claimed_complete_proof_salvage_audit.md`
  §2.5). The real gain is the isotriviality remark: a single fixed variety `X_M` serves every `N`.

**Lemma 4.1. Correct.**

* `d <= |z_i - z_j|` because `z_i - z_j = d w` with `w` a nonzero Gaussian integer.
* `||P|| >= R'/sqrt 2` in any arc position.
* `C R^alpha / d = C R'^alpha d^(alpha-1) <= C R'^alpha`.
* `d` is the **rational** content. A common non-rational Gaussian factor is invisible to the
  projective point and is never divided out. So no Gaussian-gcd cost is discarded.

**Theorem 5.1. Correct.**

* `h_E >= lambda_{E,inf} - O(1)` holds off `E`. The local Weil functions of an effective
  divisor are bounded below at every place and are `>= 0` outside a finite set of places fixed
  by a model.
* `lambda_E = lambda_Y o pi + O(1)`, because `pi^*Y = E` as subschemes.
* `h_{F_i} >= -O(1)` off `F_i`.
* Arithmetic: `-2 log R' + 5((1/2) log R' - log C) = (1/2) log R' - 5 log C`.
* `V(7, 1/4)` then gives `R' <= c_0 C^20`. With `R' >= sqrt(R)/C`, this yields the threshold
  `R > c_0^2 C^42`.
* For smaller `R` the count is `<= C sqrt R + 1 <= c_0 C^22 + 1`.
* `W = phi pi(Z)` is proper, closed and independent of `C`, and `M_0 = 2 deg F_W + 6`.
* The constant can be loosened: any fixed `eps < 1/2` works, giving `R' <= c C^(5/(1/2-eps))`.

**Pitfall checklist (the user's list).**

* *Template-dependent thresholds treated as uniform:* none. `c_0`, `Z` and `W` depend only on
  `(X^#, A, eps)` and the chosen heights.
* *Discarded gcd/denominator costs:* none. See Lemma 4.1 above.
* *Drifting phases set to zero:* none. Every estimate is invariant under rotating the arc.
* *Non-uniform exceptional sets:* none. With `D = 0`, no place set `S` occurs. Vojta's
  exceptional set is the standard one for a single fixed variety.
* *Extraction steps assumed:* none. The argument applies to every 7-subtuple of the arc
  directly.

**Theorem 5.2. Correct.** `M_alpha = 3 + floor(2/(1-alpha))` gives
`(M_alpha - 2)(1-alpha) = (1 + floor(x))(1-alpha) > x(1-alpha) = 2` with `x = 2/(1-alpha)`, so
`eps_alpha > 0`.

**Corollary 5.3 (only sketched in the note). Correct; details checked.** Take divisors
`d_j in [sqrt n - c n^(1/4), sqrt n + c n^(1/4)]`. Then:

* `|n/d_j - sqrt n| <= c n^(1/4)(1 + o(1))`.
* All coordinate differences are `<= 2c n^(1/4)(1 + o(1))`, and the content satisfies
  `g <= 2c n^(1/4)`.
* So `h = (1/2) log n - log g + o(1)`, which lies between `(1/4) log n - log 2c` and
  `(1/2) log n`.
* `lambda_Y >= (1/4) log n - log 2c - o(1) >= h/2 - log 2c - o(1)`, and `h -> infinity`.
* For Ruzsa's window `n^(1/2-eps)`: `lambda_Y >= eps log n - O(1) >= 2 eps h - O(1)`, so
  `M - 2 > 1/eps` is needed, which is `M = 3 + floor(1/eps)` as stated.
* The conics `xy = n` are irreducible, so the Bezout lemma applies.

**Proposition 6.1 and Corollary 6.3. Correct.** Ru–Vojta (Amer. J. Math. 142, 2020) with one
Cartier divisor `E`, `L = H` big and `S = {infinity}` gives a fixed exceptional set. The
re-proof of Cilleruelo–Córdoba is qualitative and ineffective, so it is weaker than the
classical theorem, which gives explicit counts. It is a soundness check, not a new result.

**Proposition 6.2. Correct as a lower bound.**

* Basis `t^k u^(e+) v^(e-)`: distinct `(k, e)` give distinct torus characters
  `s^(k + sum e-) prod lambda_j^(e_j)`.
* The local form near `Y` is `(uv)^((N-s)/2) u^s w^(e')`. Distinct `s` give distinct powers of
  `v` in the chart `u = 1`, so `ord_Y = min_s ord_1(P_s)`.
* `h_A(d) <= binom(d+n, n)`.
* The volume identity holds via the unimodular shear `(sigma, x) -> (sigma - sum x, x)`.
* Jensen's inequality gives `beta_M >= 2n/(n+1)^(1+1/n)`. I recomputed the values:
  0.5, 0.7698, 0.9449, 1.070, 1.1647, 1.2395, 1.3002 for `M = 2..8`.

**Proposition 8.1, Corollary 8.2 and Proposition 8.3. Correct.**

* `w = g eta u-bar` follows from `gcd(u, u') = 1` and equal norms.
* `u` conjugate-primitive gives `gcd(A, B) = 1` and `A + B` odd, so `gcd(s, t) = 1` in all four
  cases.
* The `abc` step gives `rad(n) >= c_delta C^(-1) R^(1-alpha-2delta)`. Corollary 8.2 follows
  from `n | R^2`.
* Proposition 8.3 is a triviality, but it is correct.

## 2. Exact computations

| check | what | result |
|---|---|---|
| rerun | `check_vojta_geometry.py`, `check_abc_pairs.py` | logs reproduced byte-for-byte |
| R1 `ref_beta_direct.py` | `h^0(NH - mE)` computed **directly** in the real circle model (chart `z_1 = (1,y)`, `z_j = z_1 c(s_j)`, `c(s) = ((1-s^2) + 2si)/(1+s^2)`), independent of the author's `A_s` reduction | agrees with the author's `A_s` totals for `(M,N)` = (2,4), (2,8), (3,2..5), (4,2), (4,3), (5,2), mod two primes. So the reduction behind Proposition 6.2 and the `beta` table is right. |
| R2 `ref_curve_intersections.py`, `ref_pell_only.py`, `ref_cg6_m9_sample.py` | exact `E.Gamma~` (homogeneous gcd of the difference forms over `Q(i)`) and `K.Gamma~ = -2d + (M-2) E.Gamma~` for the known families | see §3. Every subtuple was enumerated except CG `d=6` on `X_9`, where a timeout cut the full run; there, 3000 random 9-subtuples were checked, all with `E.G = 2`, `K.G = 2`. |
| R3 `ref_circles.py` | Lemma 4.1(a),(c), as exact integer inequalities, on 24,488 short-arc subtuples (`C` = 1, 2, 4, 8; sizes 2–7) from 259 random circles (split, inert and dyadic factors) | 0 failures. The `Phi = 5 lambda_Y - 2h` bound holds on all 7-subtuples (min slack 1.78). |
| R4 `ref_circles.py` | Proposition 8.1 identity on 2,864,384 pairs (own Gaussian gcd) | 0 failures |
| Machin | `239^2 + 1 = 2*13^4`; `28561` and `28560 - 239i` on `13^8` | true; distance is exactly `sqrt2 * sqrt R`; `n = 13^4`, `rad n = R^(1/4)` |

## 3. Is the hypothesis consistent with what is known?

A conditional theorem is empty if its hypothesis is false. I looked for ways to refute
`V(7, 1/4)` (and `V(M_alpha, eps_alpha)`); none succeeded.

1. **Oganesyan, "Lattice points on small arcs" (arXiv:2107.09991).** The abstract claims
   unbounded counts on arcs `R^alpha` for every `alpha in (1/2, 1)`. Combined with Theorem 5.2
   that would refute Vojta's conjecture. The paper was **withdrawn** on 22 Aug 2021 because of
   a mistake in Lemma 2. The notes' `endpoint_recent_literature_check.md` already records this.
   No contradiction remains.

2. **Exact intersection numbers of the known short-arc curves** (R2):

   | curve | points | `d` | `E.Gamma~` | `K.Gamma~` on `X_6` / `X_7` / `X_8` / `X_9` |
   |---|---|---|---|---|
   | CG `prod_{k<=4}(a+k+i sigma_k)`, `sum sigma = 0` | 6 | 4 | 2 (at `a = infinity`) | 0 / – / – / – |
   | CG `d = 5`, `sum sigma = 1` (arcs `~R^(3/5)`) | 10 | 5 | 2 | – / **0** / 2 / – |
   | CG `d = 6`, `sum sigma = 0` (arcs `~R^(2/3)`) | 20 | 6 | 2 | – / −2 / 0 / 2 |
   | eight-point translated Pell quartic | 8 | 4 | 2 (at `X^2 = 5`) | – / **2** / 4 / – |

   * The Pell quartic is `K`-positive in `X_7`, with `K.Gamma~/d = 1/2 > 1/4`. So `V(7, 1/4)`
     forces it into `Z`, as §5.1 says. It is one curve, so there is no conflict.
   * The value `M_alpha = 3 + floor(2/(1-alpha))` is exactly the first `M` at which the CG
     curve with `alpha = 1 - 2/d` becomes `K`-positive. The `M_alpha` threshold of Theorem 5.2
     therefore cannot be lowered by this method.
   * The same holds for `M = 7` at `alpha = 1/2`: the six-point CG family has `K.Gamma~ = 0` in
     `X_6`, so no `D = 0` Vojta instance on `X^#_6` can give the result through this
     inequality.
   * The 10-point `d = 5` family sits exactly at the `V(7, .)` threshold `alpha = 3/5`.

3. **The known `K`-positive families are not Zariski dense (proved here).**
   * Take any family whose points are values `z_sigma = prod_k beta_k^(sigma_k)`, with
     `beta^(-) = conj(beta)` and `sigma` in a set `Sigma` of `M` vectors in `{+-1}^d`.
   * Projectively, `[z_sigma] = [e^(i sigma.theta)]_sigma` with `theta_k = arg beta_k`, because
     the moduli are common.
   * So the Zariski closure of the whole family has dimension `<= rank(Sigma) <= min(d, M)`.
   * For contact-2 CG-type curves, `K`-positivity needs `M >= d + 3 > rank(Sigma)`. For the
     Pell type, `rank <= 5 < 7`.
   * Hence none of the known families can contradict `V(7, 1/4)`. This is consistency
     evidence, not evidence for the conjecture.
   * A dimension count over constrained product families suggests the same for any contact
     order, but that part is heuristic.

## 4. Issues

1. *(minor)* §6 "Where it stalls", §9 and the summary say the Ru–Vojta route "reaches only
   exponents > 1/2", that the gap is `Theta(log M/M)`, and that `beta_7 ≈ 1.25`.
   * Only the lower bound `beta_M >= 2 - O(log M/M)` is proved.
   * The numerics stop at `N <= 4` for `M = 7`. There `beta_N ≈ 1.04`, far from the
     asymptotic regime, so the ~1% ratio does not determine `beta_7`.
   * Whether `beta_M > 2` for some `M` is open, as the note itself says under "Open
     sub-question". These sentences should be labelled heuristic.
2. *(minor)* Proposition 3.2(ii), converse. `T^(M)_{C,alpha}` is defined by pairwise chords,
   not arcs. Points with pairwise chords `<= l` lie on an arc of length `<= (pi/2) l` (when
   `l <= 2R`). So the right statement is `T^(M((pi/2)C)+1)_C = empty`. The equivalence is
   unaffected.
3. *(minor)* Corollary 5.3 is only sketched in the note. The details are supplied in §1 above.
4. *(minor)* The summary says "Bombieri–Lang alone gives **no** implication". The note itself
   says "nothing that I could find". Only the weaker wording is supported.
5. *(minor, context)* Corollary 8.2 is unconditional, and even effective, when `R^2` has a
   bounded number of prime factors, via Baker's bound on linear forms in `log(pi_p/pi_p-bar)`.
   Its genuine `abc` content is therefore the powerful case with many primes.
   Proposition 8.3 is a one-line triviality. It shows only that these particular relations are
   useless for `abc`, not that `abc` is useless in general.
6. *(scope, not an error)* The hypothesis `V(7, 1/4)` is much stronger than the target:
   * the same hypothesis (with any `eps < 1/2`) gives an absolute `M_0`;
   * `V(7, eps)` with small `eps` also handles every `alpha < 3/5`.

   It is a Roth-type statement for approximating a curve on a 7-fold, far beyond current
   methods (Subspace and Ru–Vojta reach only `beta_M`). The reduction is clean, but it does
   not bring the unconditional problem closer. By Proposition 3.2, the missing lemma is
   **exactly** equivalent to the target.

## 5. Novelty

* **Project notes.** "Vojta" appears in 24 notes. Every one is the Ru–Vojta, Autissier or
  Subspace (GCD) machinery, either on `Bl_{[1:1:1]} P^2` with `S` = primes of `N` or in the
  Segre/conductor routes, and all of them hit the `S`-uniformity wall
  (`research_reset_inverse_concentration.md`, `uniform_closed_set_vs_finite_remainder.md`).
  * No note uses Vojta's Main Conjecture with `D = 0`, the fibred power `X_M`, its canonical
    class, or an `abc`-conditional statement. I grepped for: Vojta's/main conjecture,
    Bombieri–Lang, Caporaso, fibred power, isotrivial, toric, abc conjecture, "under abc".
  * The `S`-free mechanism is new to the project. `X_M` encodes the common norm in its
    geometry, so no `S`-unit structure is needed. The `M = 7` threshold and the Bezout
    uniformization are also new here.
* **Literature** (WebSearch only; direct fetches of arXiv, UCL and UAM were blocked by the
  proxy).
  * I found no derivation of the `C sqrt R` arc bound, the Cilleruelo–Granville `R^(1-eps)`
    conjecture or Erdős–Rosenfeld/Ruzsa (Erdős #887/#886) from Vojta, Bombieri–Lang or `abc`.
  * Cilleruelo–Granville (math/0608109) use Bombieri–Lang and `abc` only for sumsets of squares
    and squares in arithmetic progressions.
  * Chan's papers (1303.2069, 1406.2230) are unconditional, and his `abc` papers treat nearby
    but different problems.
  * The **method** is standard: McKinnon (2007) and McKinnon–Roth (Invent. Math. 2015) derive
    approximation exponents from Vojta on blow-ups; Silverman (2005) and Yasufuku work with
    Vojta for blow-ups. The note should cite them.
  * The contribution is the application: the choice of `X_M`, the canonical-class computation,
    the isotriviality and Bezout uniformization, and the sharp `M_alpha`. That is modest but
    real.

## 6. Corrected statement (wording only; the mathematics stands)

> **Theorem.** Assume Vojta's Main Conjecture for rational points with `D = 0` and `A = H` on a
> smooth projective model `X^#_7` (over `Q`, isomorphic over the smooth locus) of
> `Bl_Y X_7`, for one fixed `eps < 1/2` (e.g. `eps = 1/4`). Then there are an absolute `M_0`
> and an ineffective `c_0` such that every arc of length `<= C sqrt R` on `x^2+y^2 = R^2`
> (`R^2 in Z`, any position) has `<= M_0` lattice points once `R > c_0^2 C^42` (for
> `eps = 1/4`). Hence `M(C) <= max(M_0, c_0 C^22 + 1)`.
>
> `V(M_alpha, eps_alpha)` with `M_alpha = 3 + floor(2/(1-alpha))` gives the same for arcs
> `C R^alpha`, for every `alpha < 1`. The value `M_alpha` is sharp for this method: the
> Cilleruelo–Granville curves of degree `d = 2/(1-alpha)` have `K.Gamma~ = 0` on
> `X_{M_alpha - 1}`. The split model gives the Erdős–Rosenfeld and Ruzsa statements.
>
> **Unconditional.**
> * Proposition 3.2: the uniform theorem holds iff for each `C` some `T^(M)_C` is not Zariski
>   dense.
> * Ru–Vojta with `S = {infinity}` and `beta_M >= 2n/(n+1)^(1+1/n)` reproves the qualitative
>   Cilleruelo–Córdoba theorem for `alpha < 1/2`.
> * Whether `beta_M > 2` for some `M`, which would give the unconditional theorem, is open.
>   The `beta` numerics are small-`N` evidence only.
>
> **Under `abc`:** Proposition 8.1 and Corollary 8.2 as stated. When `omega(R^2)` is bounded,
> Corollary 8.2 holds unconditionally by Baker.
