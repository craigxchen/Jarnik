# Referee report on `round12/fake.md` (adversarial fake test for the T = 1 programme)

Referee scripts: `referee_fake_checks/` (independent library `rlib.py`, built from the definitions in
`round5/exact.md` Thm 4.1, `round5/construct.md` §2 and R11 compute §3; no code from `fake_checks/` is imported).
Notes: `referee_fake_notes.md`.

## 0. Verdict

**The theorems survive, with small repairs. The computational reading of §4, §6 and §7 does not.** No growth
bound and no uniform bound are claimed, and none follows.

What survives (re-derived or recomputed independently):
* **Thm 2.1** (fake directions). Correct up to a factor-2 slip in the direction estimate.
  * It is the formal version of R11 §8.3 item 4.
  * Its "no Diophantine hypothesis" half holds only as `P -> infinity`, where the fakes are **not**
    barrier-scale.
  * The barrier-scale half (Remark 2) needs (2.2) for all relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`,
    not only for pairs.
  * "Neither easier nor harder than ZE-TC" is overstated: the two directions need different hypotheses.
* **Thm 2.2** (lifting). Correct, modulo a cited uniform sublevel (Pyartli / KM) lemma that is not proved.
* **Prop 3.1** (sdc7) and **Prop 3.2** (sdc11). Independently certified with my own exact bases and my own
  certificates. Cor 2.3 and Cor 3.3 hold.
  * I checked that they are not an artefact of the base-point choice: for sdc7 and sdc11 the top-degree system
    is the same for every base point.
* **Prop 3.4.** Correct, but only for the base point `inf`. At that base the box is the Paley-core / `M8'` box, in
  which the twins are invisible. So it is R11 Thm 5.2 relabelled, not a statement about the fake mechanism.
* **(Mar)** for sdc and pal; the zonotope / section / product lemmas; Prop 5.3. Correct.
* **pal11** (`L = 2`, `m = 8`) is faithful. `zeta <= 22/3` is independently certified.

What does not survive:
1. **Three of the four "overdetermined" models are unfaithful (exact, major).** sdc7, pal7 and sdc11 at `L = 2`
   contain lattice zonotopes of length `7 = m`, all with braid generators. So at `m = 7` their real tangent
   cones lie in the braid arrangement **by a single character**, and `kappa >= 5040 sep^7/(|c|_1 d^7)`.
   * The "collapse onto collision strata, `kappa*(s) ~ s^{4..8}`" of §4.2, §6 and §7 item 3 is this character.
   * By Prop 5.3 / R11 Prop 2.2 it carries no information about the regime.
   * The write-up's caveat ("not attained by any zonotope found so far") is false.
2. **The sdc11 real tangent cone has not shrunk (certified, major).** I certified (Krawczyk, (F) checked) four new
   regular collision-free common zeros of the five octics, the best with `sep >= 0.1062`. This refutes:
   * "none found `>= 0.02`";
   * "max sep `< ~11%` of the maximum";
   * "the separation available on the cone has dropped by an order of magnitude";
   * "the cone shrinks quickly with the number of forms".
3. **The landscape optimiser is unreliable, and the failures were omitted (major for the evidence).**
   * In sdc7 at `m = 8` the author's own run reports positive minima `2e-13, 2e-10, 5e-7` at `s = 0.02, 0.05, 0.15`,
     where the certified truth is 0. These rows are left out of the §4.2 table.
   * The "arithmetic-progression minimisers next to cheap strata" appear in models with zeros as well.
4. **Normalisation error (minor).** With the base pinned at 0 and `|u|_inf = 1`, the maximal separation is
   `1/ceil(d/2)` (`1/4`, `1/6`), not `2/d`. All the percentages change.

**Net effect on route 5.**
* The only faithful overdetermined fake model is pal11 at `m = 8`, and for it there is only unreliable SLSQP
  evidence.
* So the computations give **no** evidence, in either direction, that overdetermination kills nondegenerate
  tangent-cone points of fake profiles.
* What remains is Theorem 2.1, a clean (and essentially expected) reduction of the fake test to ZE-TC on fake
  pattern sets, plus Theorem 2.2 and its two underdetermined realisations.

## 1. Theorem 2.1 (fake directions) and Corollary 2.1'

**Verdict: correct up to a constant slip. It formalises R11 §8.3 item 4, and its two halves need different
hypotheses.**

I re-derived the proof line by line against `round5/exact.md` Prop 1.1.

1. *Sufficient condition (2.1).* The distance to `(pi/4)Z` is the minimum of `d_A` and `d_D`. So (2.1) implies both
   inequalities of Prop 1.1, the diagonal one with room `sqrt 2`. Correct.
2. *Weights.* `w(v(lambda)) - W/2 >= tau sigma_lambda - 1/8`. The loss is `<= w_j/2 - tau/2 <= 1/(8r)` per vanishing
   coordinate, at most `r` of them. Correct.
3. *Class 1.* `|lambda.xi| <= eta ||lambda||_1`, so `|lambda.delta| >= Delta g_lambda/6`. Only finitely many
   `lambda`, and nondegeneracy is used only here. Correct.
4. *Class 2.* With `den_T <= n/2` the left side is `>= (pi/4 - 1/8)/den_T`. Correct.
5. *Class 3.* `lambda.xi` has density `<= 1/(2 eta |lambda|_inf)`, and the interval count is `||lambda||_1 Delta + 1`.
   The factor `1/eta = e^{W/4}/theta'` is beaten by `e^{-tau rho_1 ||lambda||_1/2}` exactly when
   `rho_1 ||lambda||_1 > r/2`. That is why `Lambda = r/rho_1 + 1`. Correct. The same mechanism appears in
   `exact.md` Thm 5.1 step 4.
6. *Slip (minor).* `hat u - u*` has sup norm `<= 8 eta (1 + O(eta))`, so its 2-norm is
   `<= 8 sqrt(d) eta = 2 theta sqrt(d/n) e^{-W/4}`. The proof claims `<= theta e^{-W/4}`, which is off by a factor of
   2. Replace `theta'` by `theta/(8 sqrt n)`. Nothing downstream changes.

**(Mar) for the two fake profiles (exact, `marcheck.py`).** I built `sdc_q` from Paley-I (`exact.md` Thm 4.1) and
`pal_q` with `b = 8` (`construct.md` §2) independently. Both give `r = 2q` resp. `9n - 8`, and (AI) holds.
* `sdc_q`: `G_T = Z/4` (`q = 7`) and `Z/6` (`q = 11`), generated by `lambda* = (2/n)1 - 2e_0`. Every nonzero
  multiple has `lambda.e_1 notin Z`, so it is irrelevant.
* `pal_q`: `G_T` is trivial.
* Minimum of `sigma_c/||c||_1` over all integral `c` with `||c||_1 <= 6` (`q = 7`) or `<= 4` (`q = 11`):
  `1/2` for sdc (attained at pairs) and `3/2` for pal. This agrees with `rho_1 = 1/2` and `rho_1 = 1` (the proofs
  are `exact.md` Thm 4.1(iii) and `construct.md` Thm 1).

**Quantifiers.**
* As stated (`P -> infinity`), the fakes have `W/(n log n)` unbounded. They do not contradict the target
  `n = O(log W)`, so they are **not barrier-scale fakes**.
* What they do show is logical: if `TC_H(A,m)` of the box of a (Mar)-pattern set has a nondegenerate point, then
  ZE-TC_eps is false for that box. That much follows from `u*` alone, without fakes. The fakes add that Step 0's
  conclusion and every pair-inequality separation lemma are consistent with the failure. So route 5 cannot
  exclude that pattern set at any `W`. This part of Cor 2.1' is correct.
* The barrier-scale version is Remark 2, which needs the Diophantine condition (2.2). I re-derived its
  bookkeeping:
  * The class-3 sum needs `tau (r/4 + rho_1/2) >~ n ln(n Lambda)`. For sdc (`r = 2n - 2`, `rho_1 = 1/2`) this
    gives `tau >~ 4 ln n`, i.e. `W ~ 8 n ln n`.
  * Class 1 under (2.2) needs `W >= 4 a Lambda ln n + O(log n) = 16 a n ln n` for sdc.
  * So Remark 2 is correct. But (2.2) is required for **all relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`**
    (`~4n` for sdc), not just for pairs. The sentence "(2.2) contains `Gamma(u*) = O(log n)`" is only the pair
    part.
* §0 item 1 says "no Diophantine hypothesis is needed". That is true only for `P -> infinity`. At the barrier
  scale (2.2) is needed, and the write-up says so in Remark 2. §0 should say it too.

**"Neither easier nor harder than ZE-TC" (§0 item 1, §7 item 1): overstated.** The two directions have different
hypotheses:
* (fake ⇒ ZE-TC fails) needs a *nondegenerate* point for `P -> infinity`, or a point satisfying (2.2) at the
  barrier scale;
* (ZE-TC ⇒ fakes harmless) needs the *quantitative* ZE-TC_eps plus a separation lemma.

Between these there is a gap: a box whose real tangent cone has nondegenerate points, but only ones of very
poor Diophantine quality. In that case the qualitative test fails while the barrier-scale transfer is unproved.
The honest statement is: *the fake test is implied by qualitative ZE-TC failure, and at the barrier scale by
failure at a point satisfying (2.2)*.

## 2. Theorem 2.2 (lifting) and Corollaries 2.3, 3.3

**Verdict: correct, modulo one cited lemma.** The cited lemma is a uniform sublevel estimate, quoted, not proved.

1. *Parity.* "All `c in V'` have parity `(-1)^m`, since otherwise `Phi_c == 0`" uses the injectivity hypothesis:
   * an even `c in V(A,m)` with `m` odd has vanishing `m`-th moments, so `Phi_c == 0`;
   * injectivity then forces `c = 0`.

   So the sentence is correct, but only because injectivity is assumed. It should say so.
2. *Lifting.* `G_c = Phi_c/m! + t^2 R_c` with `R_c` entire in `t^2`. The IFT on the slice `{u_{x_0} = 1}` is
   correct: Euler's identity puts `u_0` in the common kernel of the differentials, and `u_0` is transverse to the
   slice. If several coordinates of `u_0` have modulus 1, the chart lies in `{u_{x_0} = 1}`, not in the
   sup-norm sphere. This is harmless.
3. *Spanning of derivatives.* The argument via homogeneity in `rho` and analytic continuation is correct. Note
   that it uses (F) at `u_0` and the fact that the image of `Q` contains `u_0` (`rho = 1`).
4. *Sublevel estimate (2.3).* This is the one non-elementary step. It needs the Pyartli estimate **uniformly**
   * in `theta` (compactness of the sphere),
   * in `b` (the derivative condition is translation invariant), and
   * in the family `chi(t, .)`, `t <= t_1`.

   The standard route works: cover `Q` by finitely many boxes; on each, one fixed directional derivative of one
   fixed order is `>= gamma'` for all `theta` near a fixed `theta_0` and all small `t`; then apply Fubini along
   lines. But this is not written out; the write-up cites BKM / KM, whose statements are for a single map. Label:
   *proved modulo a standard, cited uniformity statement*.
5. *Gap-admissibility.*
   * The interval count, the scaling `(e^{-w/2}/t)^{1/k}` and the cancellation of `1/t <= 4e^{W/4}/C` against
     `e^{-W/4}` are correct.
   * Class 1 is actually not needed here: the class-3 bound already tends to 0 for every fixed `lambda` as
     `tau -> infinity`.
   * Nondegeneracy of the lifted fakes follows from (2.1). Correct.

**Corollary 2.3 (sdc7) independently reproduced.**
* `sdc7_forms.py`: exact kernel; `dim V(A,8) = 1`, `|c|_inf = 42`, `|c|_1 = 1008`, support 687, even;
  `rank Poly_{<=8}|_A = |A| = 757` (so `rho = 8`).
* `sdc7_cert_ref.py`: my own certificate (exact rationals, a sign change of `Phi` on a segment of length
  `2e-22`, explicit Lipschitz bounds). It gives:
  * a real zero `u*` with `u*_3 = 1 = |u*|_inf` and `sep(u*) >= 0.249564`;
  * `|grad Phi(u*)|_2 >= 1.54e4 - 7e-12`;
  * `||P H P||_F >= 8.30e4 - 5e-10`.

  Prop 3.1 stands.

**Corollary 3.3 (sdc11) independently reproduced, and much stronger than stated** (§3 below).

**Base-point check (my adversarial test; negative).** [D] Thm 4.1 holds for every labelling, so an actual cluster
lies in `Z_P(A_b)` for **every** choice of base point `b`. I tested whether the "fakes on `Z_P`" of Cor 2.3 / 3.3
are an artefact of using one base point. They are not:
* **sdc7.** The eight base-point boxes are different sets of size 757, but they carry the **same** measure `c`
  (`multibase7b.py`, exact). Its support has `|k_x| <= 1` in every coordinate, base included. So `Z(V(A_b,8))` is
  independent of `b`.
* **sdc11.**
  * All five basis measures of `V(A_inf, 8)` are supported on `{k_inf = 0}` (`supp_out.txt`), so they lie in
    every `A_b`.
  * `dim V(C,8) = 5` mod `p` for the base-free core `C = A ∩ [-1,1]^12` (9263 points; `sdc11_bases.py`).
  * The base-`b != inf` boxes are larger (9553 points), but `dim V(A_1,8) = 5` mod `p`. That is an upper bound
    for the rational dimension, and the five rational measures lie in `V(A_1,8)`. So `V(A_1,8) = V(A_inf,8)`.
  * The cyclic shift fixes `inf` and acts transitively on the field points, so this holds for every `b`.

  So Cor 2.3 and Cor 3.3 survive the multi-base strengthening.
* **Side remark (structure the write-up missed).** In sdc11 the whole top-degree system is supported on
  `{k_inf = 0}`. So `Z(V(A,8))` is invariant under moving the base point `delta_inf`, and the five octics are
  forms in the 11 field-point positions only. The "12-point" fakes on `Z_P` are therefore 11-point
  configurations plus one **unconstrained** point. This changes no statement, but it should be said: in this
  model the base point is a free (cylinder) direction of `Z_P`.

## 3. Computations: models, Hilbert functions, certificates

**3.1 Hilbert functions (reproduced).** `hilb.py`: ranks of reduced monomials (`k_x^3 = k_x`), two primes, base
row 0 (= `inf`) and base row 7. All `n = 8` rows of the §3.1 table agree with mine at base `inf`:
* sdc7: L2 and L3;
* pal7: L2 `1 7 28 77 161 251 241 95`, and L3 `... 195 63 7`;
* syl8 = `M8` / `M8'`.

`sdc11_bases.py` gives `dim V(A,8) = 5` for sdc11 L2 mod `p`, and `sdc11_inv.py` gives an exact `Z_11`-invariant
basis of dimension 5, verified on all 17469 reduced monomials of degree `<= 7`. So `dim_Q V(A,8) = 5`, as claimed.
I did not recompute `H_A(7) = 2157` for sdc11 or `H_A(8) = 436` for pal11. Note that a mod-`p` value of
`dim V` is only an **upper bound** for the rational dimension. "pal11: 436 forms" is therefore `<= 436` until a
rational lower bound is supplied. `rho(pal11) <= 8` is rigorous, because full rank mod `p` implies full rank.

**3.2 Base-point dependence (quantifier error in Prop 3.4; minor-to-major).** The `T = 1` box depends on the base
point, because the base coordinate `k_b = -sum k'` is not bounded by 1. For pal7 (`iso.py`, `hilb.py`):

| base | `|A|` at L2 | `|A|` at L3 | `H_A` at L3 |
|---|---|---|---|
| `inf` (row 0) | 861 | 1877 | `... 195 63 7` (`rho = 11`; this is `M8'`) |
| rows 1, 4, 6, 7 | 807 | 1781 | `... 383 311 157 33` (`rho = 10`) |
| rows 2, 3, 5 | 801 | 1757 | not computed |

* The identification "pal7 at `L = 3` is `M8'`" holds only at base `inf`. There the twin slabs do not cut, so
  the box **is the Paley-core box**. My search finds 1344 point permutations mapping the Paley-core cuts onto
  the Sylvester cuts, 168 of them fixing row 0. The author's `(0,1,2,4,3,6,7,5)` is among them.
* Prop 3.4 is therefore correct as "no gap-admissible pal7 fake lies on `Z_P(A_inf, 11)`". An actual cluster lies
  on `Z_P(A_b)` for every `b`, so this is a valid exclusion.
* But it says nothing about the fake mechanism. At this base and `L` the twins, which are what makes the Paley
  profile a fake (integral characters), are **invisible**. The box is that of `H_8`, and Prop 3.4 is R11 Thm 5.2
  relabelled.
* The phrase "qualitative ZE-TC^braid for a genuine fake profile at an overdetermined order" (§3.4, §7 item 3)
  should read: *the Paley-core / Sylvester box `M8'`, which happens to be the base-`inf` `L = 3` box of pal7*.
* It is also LP-subcritical: the slab width is 12, while `m = 11`.

**3.3 The normalising maximum of `sep` is wrong (minor, but it propagates).**
* The write-up uses `2/(n-1)` (`0.2857` at `n = 8`, `0.1818` at `n = 12`) as "the largest possible `sep`" at
  `|u|_inf = 1`. That is the maximum for `n` free points in an interval of length 2.
* Here the base point is pinned at 0 and `|u|_inf = 1`. One side of 0 then carries `>= ceil(d/2)` points inside
  an interval of length 1, so `sep <= 1/ceil(d/2)`. Equally spaced points attain this.
* So the maximum is **`1/4` at `n = 8` and `1/6` at `n = 12`**. Consistently, the author's own `tcnum` run at
  `s = 0.25` returns the progression `(-0.75, ..., 1)`.
* The percentages in §0 item 5, §3.3 and §6 change accordingly:
  * sdc7: `0.2496/0.25 = 99.8%`, not 87%;
  * sdc11, certified by the author: `6.0%`, not 5.5%.

**3.4 sdc11: the real tangent cone has NOT shrunk (refutes §3.3, §4.2 row 4, §6 item 2).** The author reports
the following (evidence) and builds the "cone shrinks quickly with the number of forms" reading on it:
* SLSQP maximisation of the separation ends near 3-value strata;
* "none found `>= 0.02`";
* "max sep on TC `< ~0.11` of the maximum".

I maximised `t` subject to `F(u) = 0` (the 5 normalised octics), `u_0 = 1`, `|u_x| <= 1` and ordered gaps
`>= t`, with SLSQP from 470 random starts (`sdc11_maxsep.py`, seeds 0-3). Many runs end at zeros with separation
`0.04-0.106`. I then certified the four best with my own Krawczyk test (arb, 300 bits, box radius `1e-30`, five
unknowns, six coordinates fixed to rationals with denominator `10^25`), using **my own exact basis**
(`sdc11_cert_ref.py`). Output: `sdc11_cert_ref_search_out.txt`.

| certified zero | `sep >=` | `sep / (1/6)` | `sigma_min(J)` | `sigma_min` of the restricted-Hessian matrix (F) |
|---|---|---|---|---|
| Z1 | **0.106237** | **63.7%** | 6.2 | 9.3 |
| Z2 | 0.085680 | 51.4% | 2.5 | 3.4 |
| Z3 | 0.077125 | 46.3% | 2.8 | 2.7 |
| Z4 | 0.075862 | 45.5% | 0.69 | 1.3 |

Z1 has `u = (1, -1, 0.157939, -0.449516, -0.555753, -0.661989, 0.757405, -0.108272, 0.476649, 0.370412, 0.264176)`
and `Gamma = 1.32`. The author's equally spaced reference value `Gamma = 0.921` is for `n = 8`.

* Every point is regular, with independent restricted Hessians, so Theorem 2.2 applies to each.
* Errors in (F) are of order `1e-20`, against margins `>= 1`. They are computed in 80-digit arithmetic at the
  enclosed point, with a crude Lipschitz bound of `2e9 * 1e-30`.
* I also re-certified the author's two zeros with my own basis: `sep >= 0.005287` and `>= 0.010009`.

Consequences:
* The sdc11 tangent cone contains regular, collision-free points with **64%** of the largest possible
  separation, against **99.8%** for sdc7. "The separation available on the cone has dropped by an order of
  magnitude" (§6) is false.
* Within the computed models, the underdetermined top-degree cone does **not** shrink from `n = 8` to `n = 12`
  in any way the data can see.
* The §4.2 table row "sdc11, `m = 8`" reports positive `kappa*` at `s = 0.02, 0.05, 0.08`. These are optimiser
  failures: `kappa*(s) = 0` for every `s <= 0.106`.

## 4. The landscape evidence (§4, §6)

**4.1 Three of the four "overdetermined" models are unfaithful: they contain zonotopes of length 7 (exact; major).**
The write-up's caveat says "`zeta = 7` is not excluded ... the LP bound `zeta <= 7` is not attained by any zonotope
found so far". In fact it is attained, and in all three `m = 7` models the zonotope uses **braid generators only**
(coordinates `e_x` and differences `e_x - e_y`). `zono_dfs.py` / `zono_small.py` find them in 1-10 s.
`zono_verify.py` / `zono_verify11.py` check exactly that all 128 subset sums lie in `A`, that the
finite-difference measure `c` annihilates every reduced monomial of degree `<= 6`, and that
`Phi_c = 7! prod_i <v_i,u>` (identity at 10 random integer points; both sides are polynomials of degree 7, and
the identity is the Prop 2.1 / Lemma 5.1(1) computation).

| model (`T = 1`, base `inf`) | `k_0` | generators `v_1..v_7` | `|c|_1` |
|---|---|---|---|
| sdc7, `L = 2` | `(-1,0,-1,1,0,0,0)` | `e7, e6-e7, e3-e4, e3, e2-e4, e1-e2, e1-e6` | 126 |
| pal7, `L = 2` | `(0,-1,0,0,0,-1,0)` | `e7, e6-e7, e6, e5, e4, e2-e4, e2-e5` | 96 |
| sdc11, `L = 2` | `(-1,0,-1,0,...,0,1)` | `e1, e2, e1-e2, e3-e6, e3-e11, e6-e8, e8-e11` | 96 |

(Coordinates are numbered `1..d` = rows `1..d`.)

Consequences (exact, all quantifiers: these `(A, m = 7)`, `V_H = V(A,7)`, every real `u`):
* **ZE-TC^braid holds at `m = 7` in sdc7, pal7 and sdc11.** One character certifies it: every real common zero of
  the tangent forms has a collision.
* **A quantitative bound holds with one character:**
  `kappa(u) >= 5040 sep(u)^7/(|c|_1 d^7)`, i.e. `>= 4.9e-5 sep^7` for sdc7.
  This explains the "decay `s^{4..8}`" of §4.2 for these rows.
* By R11 compute Prop 2.2 and the write-up's own Prop 5.3, character certificates never exist for regime
  boxes. So these rows of §4.2, rows 3-4 of the §6 table, and §7 item 3 ("one degree below the top, the same fake
  profiles ... collapse onto the collision strata") say **nothing** about the regime. They are exactly the
  zonotopal situation R11 retired.
* The "caveat" should be replaced by the theorem: `zeta = 7 = m` for sdc7, pal7 and sdc11 at `L = 2`.

**4.2 The one faithful overdetermined model (pal11, `L = 2`, `m = 8`).**
* `zeta <= 22/3` is independently confirmed. `lpzeta_ref.py` uses the slab functionals and all 12 coordinate
  functionals (base included, width 4), all generators in `{-2..2}^11` filtered by the widths, cutting planes,
  and a final exact check over all 48.8M generators with rational weights of denominator 6. So `zeta <= 7 < 8`:
  the model is faithful, as claimed.
* Its landscape (§4.2 last row) is the only overdetermined evidence that is not explained by a character.
* It is floating-point SLSQP evidence, and §4.3 shows that this optimiser misreports minima in models where the
  truth is known.
* `H_A(8) = 436` is a mod-`p` value, i.e. an upper bound for `dim_Q V(A,8)` (§3.1).

**4.3 The optimiser behind the landscapes is unreliable, and failures were omitted from the table (major).**
* *sdc7, `m = 8`.* The author's own output `fake_checks/tcnum_sdc7_L2_m8_out.txt` reports the following minima of
  `kappa2` over `{sep >= s}`:

  | `s` | 0.02 | 0.05 | 0.15 |
  |---|---|---|---|
  | reported `min kappa2` | `2.1e-13` | `2.4e-10` | `5.2e-7` |

  By Prop 3.1 (certified, and reproduced here) the true minimum is **0** for every `s <= 0.2496`. These three
  failures are the same order of magnitude as the "positive minima" quoted as evidence for the overdetermined
  models, and their minimisers are again arithmetic progressions, e.g. `(-0.3, -0.15, 0, ..., 0.6, 1)` at
  `s = 0.15`.
* The §4.2 table prints only the successful entries for this row (`s = 0, 0.10, 0.20`) and omits `s = 0.02, 0.05,
  0.15`.
* *sdc11, `m = 8`.* The table reports positive values at `s = 0.02, 0.05, 0.08`, whereas certified zeros with
  `sep >= 0.106` exist (§3.4).

So "the minimisers sit next to cheap strata" and "`kappa2*(s) ~ s^{4..8}`" describe where SLSQP stops, not where the
minima are. The progression-plus-far-point shape is the optimiser running into the boundary of `{sep >= s}`. It
does that in models with zeros too. **No conclusion about the shape of the real tangent cone can be drawn from
these landscapes.** At most the §4.2 values are upper bounds for `kappa_H` at the printed points (that part is
correct).

**4.4 §6 (growth of `dim V`, shrinking of the cone).**
* *`dim V` growth.* The values are correct as stated (`1 -> 5` at the top degree, exact; `43 -> 2157` one
  degree lower, mod `p`).
* *"The cone shrinks quickly with the number of forms"* is **refuted**:
  * at the top degree, sdc11 has certified regular collision-free points at 64% of the maximal separation
    (§3.4);
  * one degree lower, the collapse is a single braid character (§4.1).

  Neither row measures a "shrinking" driven by the number of forms. The only faithful overdetermined datum is
  pal11, which is SLSQP evidence.

## 5. Restricted theorems (§5)

* **Lemma 5.1.** Correct.
  * (1): `|c|_1 <= 2^m`, and the degree-`m` part is `m! prod <v_i,z>^{alpha_i}`.
  * (2): moments of order `< m` of a push-forward are moments of `c'` of polynomials of degree `< m`, and
    `|<t,v>| <= d' |v|_inf`.
  * (3): binomial expansion.
* **Thm 5.2.** Correct, with one constant slip (minor). The gap lemma gives `|<v,delta>| >= 2 e^{-X_v/2}` only for
  even `lambda.k` (case (E)); for odd `lambda.k` it gives `sqrt 2 e^{-X_v/2}` (case (O) of `exact.md` Prop 1.1). For
  the pal profiles (`k = 0`) the stated form is right.
* **Prop 5.3.** Correct. At `T = 1`, `|v_i|_1 <= 2d + 2d = 2D` with `D = 2(n-1)`, so R11's inequality chain goes
  through verbatim. The fakes are nondegenerate by (2.1), so the character inequality applies to them.
* **Thm 5.4.** A correct implication.
  * The section claim checks out: with the base point `0` in `U`, `k(S_a)` depends only on `a|_U`, and
    `a|_U = 0` gives the trivial slab `k(S_a) = 0`.
  * But it is **empty as stated**: no `f` with `f(sigma) > 0` is computed, and positivity alone gives no explicit
    `ZE-TC_eps`.
  * It is LP-subcritical. The write-up says all of this.

**Overall for §5:** the zonotope mechanism is the only one with constants, and it is blind in the regime (Prop
5.3). §4.1 above shows that it is also exactly what drives three of the four overdetermined data points of §4.

## 6. Cheap subtori, twins and duplicates

I looked for cheap rational subtori that the write-up missed, and for twin and duplicate-column effects.

1. **Cheap spans in the `n = 12` models at `m = 8` lie in the braid arrangement (proved).**
   * `A` contains every ternary `k'` with `|k'|_1 <= 2`: such `k'` pass every slab at `L = 2`, and the base
     coordinate is in `[-2,2]`. So `h_A(v) >= 2(|v|_(1) + |v|_(2))`, the sum of the two largest `|v_x|`.
   * Hence `h_A(v) <= 3` forces `v = +-e_x`, and `h_A(v) <= 5` forces `v` ternary.
   * A cheap family (total width `<= 7`) is therefore one of:
     * a single direction, with `|v|_(1) + |v|_(2) <= 3`, hence at most 4 distinct entries;
     * one coordinate plus one ternary direction, whose span takes at most `2 * 3 = 6` position values;
     * at most three coordinates, which put 8 points at 0.

     In each case there are fewer than 12 distinct position values.
   * So every cheap span has collisions. This confirms that the certified sdc11 zeros (all collision-free) lie
     on no cheap span, and that ZE-TC for pal11 at `m = 8` is the braid statement. The write-up's reading is right
     here, and there are no missed multi-column cheap subtori at these sizes.
2. **Twins.**
   * In pal profiles each twin pair `S, S xor {x}` makes `e_x` the difference of two slab functionals. At `L = 2`
     the family `{1_S, e_x}` costs `4 + 2 = 6 < 8`, so it is cheap. Its span has at most 4 position values, so
     it is a collision stratum (by item 1).
   * At `L = 3`, base `inf`, the twin slabs do not cut at all (§3.2). So the twin structure, which is what makes
     pal a fake, never enters the computed pal models in any way that the tangent cone could see.
3. **Duplicates.** Repeated columns give the same slab. This is noted correctly; multiplicities enter only through
   `r` in `w`.
4. **A structure the write-up missed (sdc11).** The whole of `V(A,8)` is supported on `{k_inf = 0}` (§2). So
   `Z(V(A,8))` is a cylinder along the base-point direction. The "12-point" fakes are 11-point configurations
   with one unconstrained point. Harmless for the statements, but relevant to how "barrier transfer at `n = 12`"
   is read: the transfer is really an `n = 11` phenomenon with a free point.
5. **Base points.** Checked in §2: for sdc7 and sdc11 the top-degree system does not depend on the base point.
   For pal7 the boxes do depend on it (§3.2).

## 7. Corrected statement

> **(a) Theorem 2.1, corrected.** Let pattern data `(T, mu, k)` on `n` rows satisfy (Mar), and let `u* in R^{n-1}`,
> `|u*|_inf = 1`, be nondegenerate.
> * For every `C in (0,1]` and `theta > 0` there is `P_0` with the following property. For every `P >= P_0` and
>   every assignment of distinct primes `= 1 mod 4` in `[P, P e^{1/(4r)}]`, there are gap-admissible angle models
>   with relative direction `|hat u - u*|_2 <= theta e^{-W/4}` (with `theta' = theta/(8 sqrt n)`).
> * For every `T = 1` box `(A,m)` and every `c in V(A,m)`: `kappa_c(hat u) <= kappa_c(u*) + m theta e^{-W/4}/sqrt d`.
> * If moreover `|lambda.u*| >= n^{-a ||lambda||_1}` for all relevant `lambda` with `||lambda||_1 <= r/rho_1 + 1`,
>   the same holds at `W = O((a/rho_1 + 1) r log n)`.
> * Consequently: if `TC_H(A,m)` of a (Mar) pattern set has a nondegenerate real point, then ZE-TC_eps fails for
>   that box with gap-admissible witnesses, at large `W`; at the barrier scale, if the point satisfies (2.2).
>
> **(b) Theorem 2.2** as stated, with the uniform sublevel estimate (2.3) labelled as cited. Gap-admissible fakes
> lying exactly on `Z_P` exist for:
> * sdc7, `T = 1`, `L = 2`, `m = 8`, every base point;
> * sdc11, `T = 1`, `L = 2`, `m = 8`, every base point.
>
> In sdc11 the real cone contains certified regular collision-free points with `sep >= 0.1062` (64% of the maximum
> `1/6`). In sdc11 the system does not involve the base point at all.
>
> **(c) Overdetermined models (exact).** sdc7, pal7 and sdc11 at `T = 1`, `L = 2`, `m = 7`: `zeta(A) = 7 = m` (braid
> zonotopes). Every real tangent-cone point has a collision, and `kappa(u) >= 5040 sep(u)^7/(|c|_1 d^7)`. These are
> character (zonotope) certificates, never available in the regime (Prop 5.3). pal7 at `L = 3`, base `inf`, is
> `M8'` (R11 Thm 5.2), the twins invisible. The only faithful overdetermined fake model computed is pal11 (`L = 2`,
> `m = 8`, `zeta <= 22/3`, `dim V <= 436` mod `p`), for which only floating-point landscape evidence exists, from an
> optimiser shown to misreport minima.
>
> **(d)** No growth bound and no uniform bound. Whether the `T = 1` box of a regime-size fake pattern set has
> nondegenerate real tangent-cone points is open. No computation so far bears on it.

## 8. Scripts

All scripts are in `referee_fake_checks/`; each `x.py` writes `x_out.txt` (or the file named in the table).

| script | content | status |
|---|---|---|
| `rlib.py` | patterns (Paley-I sdc, construct.md pal with twins, Sylvester), `T = 1` boxes for any base, reduced monomials | library |
| `marcheck.py` | (AI), `G_T` (Smith form), relevance of fractional characters, small-character margins | exact |
| `hilb.py` | Hilbert functions, `n = 8` models, base `inf` and base 7, two primes | mod `p` (2 primes) |
| `iso.py`, `basefree.py` | pal7 `L = 3` box = Paley-core box = `M8'` (1344 permutations); boxes for every base | exact |
| `sdc7_forms.py`, `multibase7.py`, `multibase7b.py` | exact `V(A_b, 8)` for all 8 bases; identical measure | exact |
| `sdc7_cert_ref.py` | independent certificate for Prop 3.1 | exact + explicit bounds |
| `sdc11_inv.py`, `sdc11_bases.py`, `supp_out.txt` | exact invariant `V(A,8)` (dim 5, all monomials checked); `dim V` mod `p` for core, base `inf`, base 1 (`sdc11_*_out.txt`); support on `k_inf = 0` | exact / mod `p` |
| `sdc11_cert_ref.py` | Krawczyk (arb), separation, (F) for the author's zeros and for the new zeros (`sdc11_cert_ref_out.txt`, `sdc11_cert_ref_search_out.txt`) | certified |
| `sdc11_maxsep.py` | SLSQP maximisation of `sep` on the sdc11 zero set (seeds 0-3; `sdc11_maxsep_*.json`) | evidence (search only; results certified separately) |
| `zono_dfs.py`, `zono_small.py`, `zono_verify.py`, `zono_verify11.py` | length-7 braid zonotopes in sdc7, pal7, sdc11 (`L = 2`) and exact verification | exact |
| `lpzeta_ref.py` | independent LP bound `zeta <= 22/3` for pal11 and sdc11, exact final check over all generators | exact (rational dual) |
| `pal11_dim.py` | intended: `dim V(A,8)` mod `p` for pal11 | **not completed** (memory-bound list conversion; killed). `H_A(8) = 436` for pal11 is therefore not independently rechecked. No conclusion above depends on it. |
