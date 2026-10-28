# Round 22 / lit: Diophantine approximation of RATIONAL points on varieties vs. (L_5^circ) and A_5(delta)

Status: complete (round 22, literature lens). Unrefereed. Scripts: `lit_checks/` next to this file; notes:
`lit_notes.md`.

## 0. Verdict summary

**No theorem in the literature implies `A_5(delta)` for any `delta < 1/30`, or `(L_5^circ)`.** Verdict: `no_new_result`
for the target; one translation proposition (Prop 1.3) and one exact computation (§1.3) are new here.

1. **Translation (§1.1, [D]).** With `X = Lbar_5^tw` (toric Fano 4-fold of the anisotropic torus `T_5`), `L = -K_X`,
   `H = H_{-K}([zeta]) = R''^2` and `dist ≍ diam/R''`, an A-point at `delta` is `dist(y,e) <= H^{-(1+2delta)/4}`, i.e.
   McKinnon-Roth exponent `gamma = 4/(1+2delta)`. The record needs **`alpha_ess,e(-K) > 15/4` on balanced points**
   (Liouville gives `>= 15/4` balanced, `>= 10/3` in general; Vojta `V_E(5)` predicts `>= 4 = dim X`; best approximations
   `alpha_e(-K) <= 2` sit on diagonal conics and are irrelevant).
2. **Cap (Prop 1.3, [D] from the refereed mean-width theorem).** Every section / Seshadri / effective-threshold
   argument at `e` — the unconditional engine of McKinnon-Roth's Liouville theorem (with stable base locus), of
   Seshadri-constant bounds, and (by assessment) of Huang's torsor lower bounds — gives at most `15/4` on balanced
   points, for every toric `Q`-divisor `L`: `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}`. Exact ranks: at levels `m <= 5` the best
   section of `m(-K)` has `ord_e/m = 10/3` (the Vandermonde), i.e. **McKinnon-Roth-Liouville at `e` reproduces
   Cilleruelo-Córdoba and nothing more**.
3. **McKinnon's conjecture and all its proved cases** (McKinnon 2007; McKinnon-Roth 2016 cubic surfaces; Huang ANT 2021
   split toric; McKinnon-Satriano TAMS 2021; Lehmann-McKinnon-Satriano 2024; Monahan-Satriano; McKinnon cubic
   hypersurfaces; Manon-McKinnon-Satriano wonderful varieties) concern the **best** approximation constant and where
   it is attained; they give no lower bound on `alpha_ess` beyond `alpha_e`. All are split or homogeneous; the
   conditional ones take their quantitative input from Vojta on `Bl_P X`, i.e. from `V_E(5)` itself.
4. **Roth-type theorems** (McKinnon-Roth 2015, Ru-Vojta, Grieve, Heier-Levin) are weaker than Liouville at a rational
   target.
5. **Counting / local distribution** (Huang IJNT 2017, Acta Arith 2019, ANT 2021, arXiv 2111.01509; Huang-de Saxcé
   2025; Batyrev-Tschinkel for anisotropic tori, which does apply to `X` and gives `#{H <= B} ~ cB(log B)^{14}`) work at
   the natural scale with relative errors and cannot exclude a sparse Zariski-dense set of A-points.
6. **lit18 dismissal re-examined (§7).** In the `T_5` formulation the "drifting phase" objection is void: the phase is
   quotiented (not zeroed), the gcd cost is accounted for by Thm 3.4', and the exceptional sets of all candidate
   theorems are fixed by `(X, L, e)`. The correct reason for inapplicability is item 2-5: the theorems are legitimate
   but too weak.
7. **The gap (§8).** What is needed is an unconditional instance of "`alpha_ess > tau_x(-K)` at a rational point".
   No such theorem was found for any variety; in every example checked (split `P^1 x P^1`, split dP6, twisted `k = 3, 4`)
   the essential constant equals the section/Liouville value. Twisted `k = 5` is predicted (by counting and Vojta) to be
   the first case where it does not.

## 1. The target in the literature's normalisation (exponent arithmetic)

### 1.1 Dictionary

Notation of the record (`round20/gauss.md` §1): `T_5 = (R^1_{Q(i)/Q} G_m)^5 / diag`, the norm-one similarity torus
(dimension `n = 4`, anisotropic over `Q`: complex conjugation acts by `-1` on the cocharacter lattice
`N = Z^5/Z1`, so `T_5` has **no** nontrivial `Q`-split subtorus); `X := Lbar_5^tw` its smooth projective toric
compactification (permutohedral fan, rays `1_A`, `0 != A ⊊ [5]`, Galois swapping `1_A <-> 1_{A^c}`); `e` the identity,
a `Q`-rational point of the open torus; `Y = Bl_e X`, exceptional divisor `E`.

The literature (McKinnon 2007; McKinnon-Roth 2015/2016; Huang; McKinnon-Satriano; Lehmann-McKinnon-Satriano)
uses the following normalisation [S: definitions as quoted in search snippets of these papers; standard]:
for `x ∈ X(k)`, `L` ample, a fixed archimedean distance `dist` on `X(R)`,

* `alpha_x(L) = inf { gamma : there is a sequence y_i -> x, y_i ∈ X(k), y_i != x, with dist(x,y_i)^gamma H_L(y_i) bounded }`
  (small `alpha` = well approximable);
* the **essential constant** (Pagelot's notion, used by Huang, Def. 4.7 of ANT 15 (2021), and by
  Lehmann-McKinnon-Satriano): `alpha_ess,x(L) = sup_U alpha_x(L; points of U)` over dense Zariski-open `U`,
  i.e. the best exponent achievable by a **Zariski-dense** sequence.

**Lemma 1.1 (dictionary; [D], from `gauss.md` Lemma 1.1).** Take `L = -K_X = sum_rho E_rho`, which is ample
(`X` is a smooth toric Fano variety: the permutohedral variety is Fano, `-K` corresponding to the polytope
`P = {u ∈ M_R : <u, 1_A> >= -1 for all A}`, `M = {u ∈ Z^5 : sum u = 0}`; its width in every ray direction is `2`,
see Prop 1.3). For a primitive 5-tuple `zeta` on `|zeta| = R''`, with `y = [zeta] ∈ T_5(Q)`:
1. `h_{-K}(y) = W' = log R''^2` exactly (all archimedean toric local heights vanish because `|zeta_j/zeta_1| = 1`);
2. `dist(y, e) ≍ max_j |zeta_j/zeta_1 - 1| ≍ diam(zeta)/R''` (for `y` near `e`; any two archimedean distances on
   `X(R)` are comparable near `e`).

Hence:

| statement of the record | approximation-constant form (w.r.t. `L = -K_X`, `H = R''^2`) |
|---|---|
| A-point at `delta`: `diam <= R''^{1/2-delta}` | `dist(y,e) <= H^{-(1+2delta)/4}`, i.e. `gamma = 4/(1+2 delta)` |
| `A_5(delta)` | no Zariski-dense set of `y ∈ T_5(Q)` with `dist(y,e) <= H(y)^{-(1+2delta)/4}`; implied by `alpha_ess,e(-K) > 4/(1+2delta)`, implies `alpha_ess,e(-K) >= 4/(1+2delta)` |
| needed: `A_5(delta)` for one `delta < 1/30` | `alpha_ess,e(-K) > 15/4` (strictly), possibly only on the balanced set |
| `(L_5^circ)` | the same, restricted to `eta_0`-balanced points (all 15 `log n_A ∈ [(1-eta_0)w,(1+eta_0)w]`) in one unit class, off `Z ⊂ T_5` |
| Liouville on balanced points (pair bound `|t_ab| >= 1`, record) | `alpha >= 15/4` (diam `>= 2R''^{7/15}`, dist `>= H^{-4/15}`) |
| Cilleruelo-Córdoba, all distinct 5-tuples | `alpha >= 10/3` (diam `>~ R''^{2/5}`, dist `>~ H^{-3/10}`) |
| counting heuristic / Vojta `V_E(5)` (gauss Thm 4.2) | `alpha_ess,e(-K) >= 4 = dim X` |
| best approximations (diagonal conics `C_A`, see 1.2) | `alpha_e(-K) <= 2` |

Arithmetic check of the rows. `dist = diam/R''` and `H = R''^2`: `diam = R''^{1/2-delta}` gives
`dist = R''^{-1/2-delta} = H^{-(1+2delta)/4}`, and `dist^gamma H = 1` at `gamma = 4/(1+2delta)`. At `delta = 1/30`:
`(1+1/15)/4 = 4/15`, `gamma = 15/4`. At `delta = 1/10` (diam `R''^{2/5}`): `(1+1/5)/4 = 3/10`, `gamma = 10/3`. At
`delta = 0`: `gamma = 4`. The Vojta computation of `gauss.md` Thm 4.2 is the general one: on `Y = Bl_e X`,
`K_Y = pi^*K_X + (n-1)E`, and for `L = -K_X`, `lambda_{E,inf} + h_{K_Y} >= (n/gamma - 1) h_L - O(1)` at points with
`dist^gamma H_L <= 1`, so `V_E` gives `alpha_ess,e(-K_X) >= n` at **any** rational point of **any** variety with `-K_X`
ample (here `n = 4`, `n/gamma - 1 = (1+2delta) - 1 = 2delta = kappa`, matching the record).

**So the needed statement is:** the essential approximation constant of the rational point `e` of the twisted
toric Fano 4-fold `X`, with respect to `-K_X`, restricted to balanced points, is `> 15/4` (any amount); Vojta predicts
`>= 4`. It is a statement about `alpha_ess`, **not** about `alpha` (the best approximations, `<= 2`, lie on curves and
are irrelevant).

### 1.2 The curves through `e` (where best approximations live)

`Q`-subtori of dimension 1: since Galois acts by `-1` on `N`, **every** rank-1 sublattice `Z v ⊂ N` is
Galois-stable, and the corresponding subtorus is `R^1_{Q(i)/Q}G_m` embedded by `w -> w^v` (`v ∈ Z^5/Z1`). Its closure
`C_v` is a conic with the rational point `e`, so `C_v ≅ P^1_Q`. [D] For `v` with sorted entries
`v_(1) <= ... <= v_(5)`, the support function of `-K` is `psi(v) = max v - min v`, so `-K . C_v = psi(v) + psi(-v) =
2(max v - min v)`. On `C_v`, `y = (P/Pbar)^v` with `P = a + bi`: `dist(y,e) ≍ |b|/|P|`, `H(y) = |P|^{2(max v - min v)}`
(for `P` primitive), so along `C_v` the best exponent is `gamma = 2(max v - min v)` (Liouville on `P^1`).
* `v = 1_A` (two distinct values): `gamma = 2`. These are the **diagonal conics** (tuples taking two values only);
  they give `alpha_e(-K) <= 2`, and they lie in the big diagonal, which any admissible `Z` contains.
* `v` with five distinct entries: `max v - min v >= 4`, `gamma >= 8`: these curves approximate `e` **worse** than
  the Dirichlet rate `gamma = 4`.

So McKinnon-type statements ("best approximations lie on a rational curve") are trivially consistent here and say
nothing about `alpha_ess`.

### 1.3 Every Seshadri / effective-threshold (Liouville-with-base-locus) bound is at most `15/4`

The unconditional mechanism behind all "Liouville-type" theorems in this literature (McKinnon-Roth 2016's
Liouville analogue with stable base locus; Seshadri-constant lower bounds; Huang's torsor lower bounds, see §5) is:
a nonzero section `s ∈ H^0(X, mL)` vanishing to order `r` at `x`, with `s(y) != 0`, gives by the product formula
`m h_L(y) >= -sum_v log||s(y)||_v >= r(-log dist(x,y)) - O(1)`, i.e. `dist(x,y) >= c H_L(y)^{-m/r}` for every `y`
off `{s = 0}`. So it proves `alpha_ess,x(L) >= r/m`, and the best such bound is the effective threshold
`tau_x(L) = sup { r/m : s ∈ H^0(mL), ord_x s >= r, s != 0 }` (`= sup{t : pi^*L - tE` big`}`).

**Proposition 1.3 ([D]; proof from refereed `threshold/construct.md` Thm B).** Let `X_k = Lbar_k^tw`, `e` the identity.
(a) For every Galois-invariant toric `Q`-divisor `L = sum_{cuts c} a_c (E_A + E_{A^c})` with `a_c > 0` and every
section `s ∈ H^0(mL)`: `ord_e s <= 2^{3-k} m sum_c a_c` (`= (m/4) sum_c a_c` at `k = 5`). In particular `tau_e(-K_{X_k}) <= 4 - 2^{3-k}`, i.e. `<= 15/4`
at `k = 5`.
(b) On `eta`-balanced points (all `log n_c ∈ [(1-eta)w,(1+eta)w]`), the section bound of (a) gives, for every such
`L` and `s`, at best `dist(y,e) >= c e^{-4w(1+O(eta))} = c H_{-K}(y)^{-4/15 + O(eta)}`: exactly the Liouville exponent
`15/4`, never more.
(c) For `L = -K_{X_5}`: the Vandermonde `V = prod_{i<j}(w_i - w_j)` is a section of `3(-K)` with `ord_e V = 10`, so
`tau_e(-K_{X_5}) >= 10/3`, and the section bound off `{V = 0}` (the big diagonal) is exactly Cilleruelo-Córdoba.

*Proof.* (a) Over `Q(i)` a section of `mL` is a Laurent polynomial `f = sum_{u ∈ S} c_u x^u` with `S ⊂ mQ_L ∩ M`, where
`Q_L = {u : <u,1_A> >= -a_{c(A)}}` is the polytope of `L`. Its width in direction `1_A` is
`max_S <u,1_A> - min_S <u,1_A> <= a_{c(A)} + a_{c(A)} = 2 a_{c(A)}` per unit of `m` (using `-1_A ≡ 1_{A^c}` in `N`),
and `0` for `A ∈ {∅,[k]}` (homogeneity). The order of vanishing at `e = (1,..,1)` is at most `reg(S)` (`upper.md`
Lemma 1.1, as used in `construct.md`), and `construct.md` Theorem B gives
`reg(S) <= 2^{1-k} sum_{A ⊆ [k]} w_A(S) <= 2^{1-k} · 2 · sum_c 2 m a_c = 2^{3-k} m sum_c a_c`, which is `(m/4) sum_c a_c`
at `k = 5`. For `-K`, `a_c = 1` for all `2^{k-1} - 1` cuts: `tau <= 2^{3-k}(2^{k-1}-1) = 4 - 2^{3-k}`.
(b) On balanced points `h_L(y) = sum_c a_c log n_c >= (1 - eta) w sum_c a_c` (gauss Lemma 1.1: each Galois pair
`E_A + E_{A^c}` has height exactly `log n_c`). If `s(y) != 0` the section bound is `-log dist(y,e) <= (m/r) h_L(y) + O(1)`,
and by (a) `m/r >= 4/sum_c a_c`, so its right side is `>= 4(1-eta)w - O(1)`. Hence it can only exclude balanced points
with `-log dist > 4(1-eta)w - O(1)`; it never excludes `dist = e^{-4(1-eta)w}`, while
`H_{-K}(y) = e^{h_{-K}} ∈ [e^{15(1-eta)w}, e^{15(1+eta)w}]`, so `e^{-4w} = H^{-4/15 + O(eta)}`: exponent `gamma = 15/4`.
(c) The width of the support of `V` (Minkowski sum of the segments `[e_i,e_j]`) in direction `1_A` is `|A||A^c| <= 6
= 2·3`, so `V ∈ H^0(3(-K))`; `ord_e V = C(5,2) = 10`. Off `V = 0` the section bound gives `dist >= c H^{-3/10}`,
i.e. diam `>~ R''^{2/5}`. `□`

**Exact computation ([E], `lit_checks/reg.py`, `reg_k5.out`).** The maximal `ord_e` of a section of `m(-K)` is
`reg(mP ∩ M)` (`ord >= r` iff the coefficient vector is orthogonal to all polynomials of degree `< r` on the support).
Ranks modulo the prime `2147483629` (an upper bound for the rational `reg`, since rank mod `p <=` rank over `Q`):

| `k` | `m = 1` | `2` | `3` | `4` | bound (a) |
|---|---|---|---|---|---|
| 3 (`|mP∩M|` = 7, 19, 37, 61) | 3 | 6 | 9 | 12 | 3 per `m` |
| 4 | 3 | 6 | 9 | 12 | 3.5 per `m` |
| 5 (`|mP∩M|` = 21, 131, 471, 1251, 2751) | 3 | 6 | 10 | 13 | 3.75 per `m`; and `m = 5`: `<= 16` |

So at levels `m <= 5`, `max ord/m <= 10/3` (at `m = 5`: `Poly_{<=16}` already spans `R^{5P∩M}` mod `p`, so `ord <= 16`, ratio `<= 3.2`; `lit_checks/reg_at.py`, `reg_k5_m5.out`) at `k = 5`, attained by the Vandermonde at `m = 3`: McKinnon-Roth's
Liouville-with-base-locus theorem for `L = -K` at `e` reproduces **Cilleruelo-Córdoba and nothing more** at these
levels; `tau_e(-K) ∈ [10/3, 15/4]`, the exact value being open (it would matter only for the record's
"intermediate exponent" problem, `gauss.md` §7 item 2, not for the target). The Ehrhart data
`1, 21, 131, 471, 1251` give `(-K_{X_5})^4 = 70` (fourth difference 70), so the dimension count alone only gives
`tau >= 70^{1/4} ≈ 2.89`.

**Consequence.** Any theorem whose proof is a section/Seshadri/effective-threshold argument at `e` gives at most
`alpha_ess >= 15/4` on balanced points and (at the levels computed) `10/3` on all points. The target needs `> 15/4`.
This is the literature form of the record's Lemma K1 + `construct.md` Cor C3 cap (row 20 of the barrier map).

## 2. McKinnon's conjecture (2007) and the approximation constant alpha_x(L)

**Source.** D. McKinnon, *A conjecture on rational approximations to rational points*, J. Algebraic Geom. 16 (2007),
no. 2 (AMS page `ams.org/journals/jag/2007-16-02/S1056-3911-06-00458-9`), arXiv `math/0504303` [S; pages 257-303 [R]].

**Statement** (Conjecture 2.7 there [S]). Let `V` be a variety over a number field `k`, `D` ample, `P ∈ V(k)`. If `P`
lies on a rational curve, then the best approximations to `P` can be chosen to lie on a rational curve: there is a
rational curve `C ∋ P` with `alpha_P(D) = alpha_P(D|_C)`. Later papers call the existence of such a curve the
**Coba** ("curve of best approximation") conjecture, and a **weak form**: `P` is better approximated by points of some
curve than by any Zariski-dense sequence, i.e. `alpha_P(D) < alpha_ess,P(D)` (Monahan-Satriano, LMS) [S].

**Proved cases** (all [S] unless marked):
* split rational surfaces of Picard rank `<= 4` (McKinnon 2007); cubic surfaces, points off lines (McKinnon-Roth 2016);
  blow-ups of Hirzebruch surfaces at special configurations ("San19", cited in McKinnon 2310.01640; author/journal
  not verified);
* split toric varieties, `P` in the open torus, under a geometric hypothesis (simplicial pseudo-effective cone)
  (Z. Huang, ANT 15 (2021) 461-512, §5);
* smooth cubic hypersurfaces with a `k`-rational line, every rational point; Coba verified (McKinnon, arXiv
  2310.01640);
* wonderful compactifications of adjoint groups of type `A_n`, unconditional Coba; weaker form for all adjoint types
  (Manon-McKinnon-Satriano, arXiv 2504.05213 [A]);
* conditional on Vojta: smooth split toric surfaces, and "better approximated on a rational curve" for
  `Q`-factorial terminal split toric varieties of any dimension (McKinnon-Satriano, Trans. AMS 374 (2021) 3557-3577,
  arXiv 2004.05212, with a 2023 erratum); all split smooth projective surfaces (Lehmann-McKinnon-Satriano, arXiv
  2403.02480 [A]; they also show that Vojta implies every point is "essentially bounded" [S]); split horospherical
  varieties, weak form (Monahan-Satriano, arXiv 2308.11847 [S]).

**Verdict for `e ∈ X = Lbar_5^tw`: does not imply `A_5(delta)` or `(L_5^circ)`, even if proved for `X`.**
1. The conjecture (and its weak form) concerns `alpha_e` = the **best** exponent and where it is attained. Here the
   best approximations to `e` are already known to lie on curves we exclude anyway (§1.2: diagonal conics, exponent
   `<= 2`). The record needs a **quantitative lower bound on `alpha_ess`** (`> 15/4`); the conjecture gives none.
   The weak form gives `alpha_ess > alpha_e`, i.e. a bound `> alpha_e` (at most `> 2` here): useless.
2. Removing the curve of best approximation and re-applying the conjecture to the complement does not iterate to a
   bound on `alpha_ess`: each step only says the next-best constant is attained on another curve, with no value.
3. In every conditional proof listed, the lower bound for approximation by Zariski-dense sequences comes **from
   Vojta's conjecture on `Bl_P X`** with `D = E`; that is exactly the record's `V_E(5)` (gauss Thm 4.2). The
   unconditional proofs (McKinnon 2007 surfaces, MR16 cubic surfaces, Huang toric, McKinnon cubic hypersurfaces,
   MMS wonderful varieties) obtain their lower bounds from Liouville/Seshadri/torsor-integrality arguments, which are
   capped at `15/4` on our balanced set by Prop 1.3 (see §3 and §5 for each mechanism).
4. All proved cases are split (or have many rational curves through `P` coming from split structure: lines on cubics,
   1-parameter subgroups of split tori, centred primitive collections). `X` is the compactification of an
   **anisotropic** torus (no split subtorus over `Q`); none of the listed results covers a non-split torus (search for
   "non-split / anisotropic torus" + "approximation constant" found nothing).

## 3. McKinnon-Roth (Invent. Math. 2015): Liouville / Roth type theorems, Seshadri constants

**Sources.** D. McKinnon, M. Roth, *Seshadri constants, diophantine approximation, and Roth's theorem for arbitrary
varieties*, Invent. Math. 200 (2015) 513-583, arXiv 1306.2976 [S]. D. McKinnon, M. Roth, *An analogue of Liouville's
theorem and an application to cubic surfaces*, Eur. J. Math. 2 (2016) 929-959, arXiv 1306.2977 [S].
**Full texts not read this round** (arXiv, ar5iv, queensu.ca, uwaterloo.ca unreachable from this session); exact
theorem numbers and constants below are from search snippets [S] or recall [R].

**Statements.**
* (Invent., Roth-type) For `x ∈ X(kbar)` and `L` ample, a lower bound for `alpha_x(L)` in terms of the Seshadri
  constant `epsilon_x(L)` (and of an "expected order of vanishing" quantity, the later `beta_x(L)`), generalising
  Roth's theorem; a secondary source quotes "`alpha_x(L) >= epsilon_x(L)/2`" (Thm 6.3) [S]; on `P^n` the method gives
  the constant `n/(n+1)` [S]; equality/near-equality forces `alpha` to be computed on a proper subvariety [R].
  No exceptional set appears in the inequality itself (it bounds the infimum over **all** sequences).
* (EJM, Liouville-type) A Liouville analogue on arbitrary varieties, plus an extension involving the stable
  (asymptotic) base locus; computation of `epsilon_x(L)` and `alpha_x(L)` for all nef `L` on cubic surfaces and all
  `x` not on a line, proving McKinnon's conjecture there [S].

**Hypotheses.** Arbitrary projective `X` over a number field, `L` ample (nef in the cubic computation); a single
place (archimedean distance), so `S = {inf}` is fixed; heights are Weil heights for `L`; no splitness assumption.
Exceptional sets: none in the Roth-type bound; the base-locus extension excludes the stable base locus of
`pi^*L - gamma E`, a fixed closed set depending only on `(X, L, x, gamma)` (so uniform in `R`; no constraint issue).

**Verdict: none implies `A_5(delta)`, `delta < 1/30`, or `(L_5^circ)`.**
1. Any lower bound for `alpha_e(L)` (infimum over all sequences) is `<= alpha_e(L)`, and for `L = -K` this is `<= 2`
   (diagonal conics, §1.2). In record units: diam `>~ R''^{0}`: vacuous.
2. The base-locus (Liouville) extension is precisely the section bound of §1.3. Its best output is
   `alpha_ess >= tau_e(L)`, and Prop 1.3 shows: on balanced points exactly `15/4` for every toric `Q`-divisor `L`
   (never more), and for `-K` on all points `10/3` at levels `m <= 5` (Cilleruelo-Córdoba). This is the record's
   `construct.md` cap re-derived in the McKinnon-Roth language.
3. The Roth-type gain (factor `1/2` vs `1`, or `n/(n+1)`) is a gain only for **irrational** targets; for a rational
   target it is weaker than Liouville (the authors note that Liouville beats Roth on `P^1` exactly at rational points
   of the base field [S]). `e` is rational. Expectation constants `beta_x(L) = ∫_0^∞ vol(L - tE)dt / vol(L)` (Ru-Vojta,
   Grieve arXiv 1705.01023, Heier-Levin) satisfy `beta <= tau`, so they are capped exactly as in item 2. (Already recorded in
   `gauss.md` §4.4 for Ru-Vojta; extended here to the McKinnon-Roth family.)

## 4. McKinnon-Roth on cubic surfaces; McKinnon-Satriano; other proved cases

| paper | exact content (as verifiable) | hypotheses | exceptional set | implies `A_5`/`(L_5^circ)`? |
|---|---|---|---|---|
| McKinnon-Roth, EJM 2 (2016) [S] | `epsilon_x(L)`, `alpha_x(L)` for nef `L` on cubic surfaces, `x` off lines; McKinnon conj. for cubic surfaces | dimension 2, `x` rational | none (computes the infimum) | **no**: dimension 2; infimum only; mechanism = Liouville + lines |
| McKinnon-Satriano, TAMS 374 (2021) [S] | Vojta ⇒ McKinnon conj. for smooth split toric surfaces; for `Q`-factorial terminal split toric varieties, `P` is better approximated on a rational curve | **split**; Vojta for the surfaces case | Vojta's | **no**: split; conclusion concerns `alpha` not `alpha_ess`; the quantitative input is Vojta |
| Lehmann-McKinnon-Satriano, arXiv 2403.02480 [A] | strategy deducing a slightly weaker conjecture from Vojta; full conjecture for split surfaces | Vojta; split surfaces | Vojta's | **no**: conditional; dimension 2; split |
| Monahan-Satriano, arXiv 2308.11847 [S] | assuming Vojta, for smooth projective split horospherical `X` and rational `P`, a curve through `P` whose points approximate better than any Zariski-dense sequence | Vojta; split horospherical (includes split toric) | Vojta's | **no**: conditional; split; weak form only |
| McKinnon, arXiv 2310.01640 [S] | `alpha_P` for every rational `P` on smooth cubic hypersurfaces with a rational line; Coba | cubic hypersurfaces | none | **no** (other varieties; infimum) |
| Manon-McKinnon-Satriano, arXiv 2504.05213 [A] | unconditional Coba for wonderful compactifications of adjoint type `A_n`; weaker form for all adjoint types | wonderful varieties (split group compactifications) | none | **no**: concerns best approximations; not a torus compactification of an anisotropic torus |
| Huang-Wilsch, arXiv 2509.03998 [A] | integral approximation constant near a boundary rational point of a log pair; conjecture (curves with `<= 2` points at infinity), verified in examples | integral points, target in the boundary | — | **no**: our `e` is in the open torus and our points are rational, not integral |

A 2025 math.GM preprint claiming general "uniform height inequalities" (arXiv 2509.05300, P. Tiebekabe [A]) was
noted and **discarded** as unreliable (category math.GM, programmatic abstract, no checkable statement).

## 5. Zhizhong Huang: local distribution / approximation constants on toric varieties

**Sources** ([S] = bibliographic data and content from search snippets; [A] = abstract read in the offline mirror):
* H1. *Distribution locale des points rationnels de hauteur bornée sur une surface de del Pezzo de degré 6*, Int. J.
  Number Theory 13 (2017), no. 7, 1895-1930, arXiv 1610.02542 [S]. A measure describing the local distribution of
  rational points near a **general rational point** of a (split) degree-6 del Pezzo surface, after removing the
  locally accumulating subvarieties [S]. Note: the split dP6 is exactly `Lbar_3` (permutohedral surface), i.e. the
  **split** analogue of our `k = 3`.
* H2. *Approximation diophantienne et distribution locale sur une surface torique*, Acta Arith. 189 (2019) 1-94,
  arXiv 1703.01772; part II, arXiv 1805.03920 [S]. Defines the essential constant at `Q`
  (`sup alpha(Q, V)` over dense constructible `V`) and **locally accumulating subvarieties** (`Z` with
  `alpha(Q, W) < alpha_ess` for every dense open `W ⊂ Z`); studies a toric surface `Y_4`, using Pell-Fermat equations [S].
* H3. *Rational approximations on toric varieties*, Algebra Number Theory 15 (2021) 461-512, arXiv 1809.02001 [S].
  Split toric varieties, general rational point `Q` of the open torus, Salberger's universal torsor method. Under a
  geometric hypothesis (simplicial pseudo-effective cone, per a citation in 2308.11847 [S]) the best approximations
  are achieved on rational curves through `Q` of minimal degree, identified with Batyrev's centred primitive
  collections; Thm 6.1 [S]: for `X` not `≅ P^n`, `alpha` is attained on a finite union of `P^{N_i}`, fibres of open toric
  subvarieties `≅ P^{N} x G_m^{n-N}`, pairwise meeting only at `Q`. Def. 4.7: essential constant (credited to Pagelot).
* H4. *Equidistribution of rational points and the geometric sieve for toric varieties*, arXiv 2111.01509 (v3, 2026)
  [A]: Manin-Peyre equidistribution on smooth proper **split** toric varieties, with effective error terms for counts in
  an arbitrary (fixed) adelic neighbourhood; effective Ekedahl sieve.
* H5. Huang-de Saxcé, *Local distribution of rational points in flag varieties*, arXiv 2510.09954 [A]: for a flag
  variety over `Q` and `x ∈ X(R)` random, or `x ∈ X(Q)` with a non-degeneracy condition, suitably rescaled
  approximations equidistribute.
* H6. Huang-Wilsch, arXiv 2509.03998 [A]: integral points near a boundary point (see §4 table).

**Hypotheses vs. our setting.** H1-H4 are **split** toric (H1, H2 surfaces; H3, H4 any dimension); H5 is homogeneous
for a semisimple group (flag varieties). `X = Lbar_5^tw` is the compactification of an **anisotropic** torus. Heights:
anticanonical (or ample toric) Weil heights, same as ours (Lemma 1.1). Exceptional sets: fixed accumulating
subvarieties determined by `(X, L, Q)`: uniform, no constraint issue.

**Verdict: none implies `A_5(delta)` (`delta < 1/30`) or `(L_5^circ)`.** Reasons, each sufficient:
1. *Split only.* The torsor method of H3/H4 parametrises `T(Q)` by integer torsor coordinates of a **split** torus;
   for `T_5` the analogue is the Gaussian parametrisation `zeta_j = eps prod pi^{a} pibar^{e-a}` with norm conditions
   (a torsor over `Z[i]` with Galois descent), which these papers do not treat.
2. *What the method can prove.* In H3 the lower bounds for approximations off the accumulating set come from the
   non-vanishing of integer-valued expressions in torsor coordinates, i.e. a section bound (§1.3). [Assessment, not a
   theorem: I have not seen H3's proofs.] In our dictionary the analogous integers are the residues `t_ab`, block norms
   and forced-divisibility quotients, and Lemma K1 + Prop 1.3(b) cap every such bound at `15/4` on balanced points. A
   twisted extension of H3 would therefore reprove Liouville (`15/4` balanced, `10/3` = Cilleruelo-Córdoba in general)
   unless it contains a genuinely new non-section mechanism.
3. *Natural-scale local distribution cannot see A-points.* H1, H2, H4, H5 describe points at the **critical** scale
   (or in fixed neighbourhoods) with errors that are `o(main term)`. A-points at `delta > 0` in the height range
   `(B/2, B]` lie in the ball of radius `r_B B^{-1/4}` around `e`, with `r_B = 2^{(1+2delta)/4} B^{-delta/2} -> 0`. A Zariski-dense
   set of A-points can be as sparse as one per dyadic range, so it changes no limit measure and no asymptotic with
   relative error; excluding it needs a count in balls of radius `B^{-1/4-delta/2}` with **absolute** error `< 1`, which
   no equidistribution theorem provides (effective errors in H4 are for fixed neighbourhoods and are of size
   `B^{1-kappa}`), and which is a restatement of `A_5(delta)`.
4. *In the split examples Liouville is sharp.* [D] On `P^1 x P^1` at `(1,1)` with `-K = O(2,2)`: `tau = 4` (section
   `(x-1)^{2m}(y-1)^{2m}`; bidegree bound) and the Zariski-dense set `(1+1/q, 1+1/s)`, `q/s ∈ [1/2,2]`, has `dist ≍ 1/q`,
   `H ≍ q^4`, and every point off `{x=1} ∪ {y=1}` has `dist >~ max(1/q,1/s) >= H^{-1/4}`: so `alpha_ess = 4 = tau`, which is
   **larger** than the counting value `n = 2` (worse approximation than Dirichlet): near a split rational point the
   points are *repelled*, and the section (Liouville) bound is the truth. On the split dP6 (`= Lbar_3`) `tau_e(-K) = 3` (table of §1.3, `k = 3`:
   `reg = 3m`), and the twisted `k = 3` case is also Liouville-sharp (record: `A_3(delta < 1/6)` false). Our `k = 5` twisted case is
   the opposite regime: the section bound `tau <= 15/4` is **weaker** than the counting value `4`, and the needed
   statement is that the truth is the counting value. No paper found treats a rational point at which `alpha_ess`
   is shown to exceed `tau_x(L)`.

## 6. Tori, twisted tori, Losev-Manin spaces; local accumulation counting

* **Batyrev-Tschinkel**, *Rational points of bounded height on compactifications of anisotropic tori*, IMRN 1995
  (cited in their preprint alg-geom/9510014 as "to appear in IMRN" [S]; issue/pages [R]: no. 12, 591-635), and
  *Manin's conjecture for toric varieties*, J. Algebraic Geom. 7 (1998) (arXiv alg-geom/9510014 [S]; volume/pages
  [R]). Asymptotic `#{y ∈ T(K): H_{-K}(y) <= B} ~ c B (log B)^{rho - 1}` for smooth projective toric varieties over number
  fields, including anisotropic tori, by Poisson summation on `T(A_K)/T(K)`; Peyre's constant; equidistribution in
  fixed adelic neighbourhoods follows from the method.
  * *Applies to our `X`* (anisotropic torus, `-K` height = `W'`). [D] Here `rho = rank Pic(X)^{Gal}`: from
    `0 -> M -> Z^{Sigma(1)} -> Pic(X_{Qbar}) -> 0` and `M^{Gal} = 0` (Galois acts by `-1`), `rank Pic(X) = #` Galois orbits
    of rays `= 15`. So `#{[zeta]: R''^2 <= B} ~ c B (log B)^{14}`, matching the record's description by 15 cut norms with
    product `N'`.
  * *Verdict:* **no implication.** Global counts and fixed-neighbourhood equidistribution cannot exclude a sparse
    Zariski-dense set near `e` (reason 3 of §5). Heuristically the natural-scale count near `e` is
    `≍ r^4 (log B)^{14}` and the A-point count `≍ B^{-2delta}(log B)^{14}`, which is the record's `e^{-2 delta W'}` heuristic
    with the log factor made explicit (**heuristic**: it assumes the local density at `e` is the Tamagawa density, which
    fails at split rational points, §5 item 4).
* **Twisted toric / Losev-Manin / permutohedral geometry.** Voskresenskii-Klyachko (toric Fano varieties of root
  systems, 1985 [R]) is the source of "the permutohedral variety and its twisted forms are Fano"; Gallardo et al.,
  *Higher-dimensional Losev-Manin spaces and their geometry* (arXiv 2308.07911 [A]) and the 2024-2026 permutohedral /
  hypergraph-associahedra items in the mirror are purely geometric. None contains Diophantine approximation.
* **Intrinsic Diophantine approximation on quadrics** (Kleinbock-Merrill, *Rational approximation on spheres*, Israel J.
  Math. 209 (2015) 293-322; Fishman-Kleinbock-Merrill-Simmons, JEMS 24 (2022) 1045-1101 [S]): Dirichlet/Khintchine
  theory for approximating points of a rational quadric by its rational points, via homogeneous dynamics of
  `SO(Q)`. Our `k = 2` slice is a conic (where this is classical), but `T_5`'s compactification is not a quadric and the
  target is rational; the relevant inequality for a rational target on a conic is the pair bound `|t_ab| >= 1` itself.
  **No implication.**
* **Counting near a rational point / "local accumulation".** The only systematic local-accumulation counting near a
  rational point found is Huang's (H1-H5) and McKinnon's surface examples; all split or homogeneous, and all at the
  natural scale. **Nothing found for anisotropic tori.**

## 7. Re-examination of the lit18 dismissal (fixed point vs drifting phase)

`lit18/tools.md` §2.8 dismissed Huang / McKinnon-Roth because, on `X_M` (split over `Q(i)`), the cluster approximates the
point `(e^{i phi} : ... : e^{-i phi})` of a curve `Y`, `phi` drifts with the cluster, and a fixed-point theorem at the
cluster's centre would replace a drifting phase by a constant.

**In the `T_5` formulation that reason no longer applies** ([D], each point checked against the owner's constraints):
1. *Phase.* `[zeta] = (zeta_j/zeta_1)_j` is invariant under every similarity `zeta -> gamma zeta`, `gamma ∈ Q(i)^*`
   (`gauss.md` §1.1). The common phase is **quotiented out**, not set to zero: the statement `A_5(delta)` quantifies
   over all tuples in all positions, and Theorem 3.4' passes from an arbitrary arc (any `theta_0`) to `A_5` via
   Lemma 3.3, which needs only that `Z ⊂ T_5` be similarity-invariant (it is, being a subset of `T_5`). The target `e`
   is genuinely fixed and rational.
2. *Gaussian gcd / denominators.* The height used is `h_{-K}([zeta]) = W' = log R''^2` of the **primitive** reduced
   tuple, an intrinsic function of the class (gauss Lemma 1.1); the passage from the cluster to the reduced tuple is
   Lemma 1.2/3.2 with its explicit loss (`W_∅ >= 2^{-k}W`, or the balanced selection). Nothing is discarded.
3. *Exceptional sets.* Every theorem in §§2-6 that could apply has an exceptional set depending only on
   `(X, L, e, gamma)` (stable base loci, accumulating subvarieties, or Vojta's `Z_eps` for the fixed variety `Y`,
   `S = {inf}`). These are fixed proper closed subsets of `X`, uniform in `R`: allowed. No fixed-`S` Subspace theorem
   with `S` = primes of `N` is involved.
4. *Thresholds.* Constants depend on `(X, L, e)` only: uniform in `R`, not template-dependent.

**So the dismissal must be replaced, not repeated.** The corrected verdict is that these theorems are legitimately
applicable to `e ∈ T_5` and **too weak**: their unconditional content is section/Seshadri/torsor-integrality
(capped at `15/4` on balanced points by Prop 1.3, i.e. exactly the record's Liouville), or concerns the best
approximations `alpha_e` (`<= 2`, on curves we discard anyway), or is proved only for split varieties, or is conditional
on Vojta (where it coincides with `V_E(5)`). The lit18 equivalence "uniform-in-`phi` analogue on `Y` = uniform theorem"
(`growth/conditional.md` Prop 3.2) is consistent with this: on `T_5` the uniformity in `phi` is built in, and what
remains is exactly `A_5` / `(L_5^circ)`.

## 8. What is missing: precise gap between proved results and (L_5^circ)

**Needed (exact form, `[D]` from §1.1).** There are `gamma > 15/4`, `eta_0 > 0`, a proper Zariski-closed `Z ⊂ T_5` and `B_0`
such that every `eta_0`-balanced `y ∈ T_5(Q) \ Z` in one unit class with `H_{-K}(y) >= B_0` satisfies
`dist(y, e)^gamma · H_{-K}(y) >= 1`. (With `gamma = 4/(1+2delta)` this gives `(L_5^circ)` for that `delta < 1/30` via
gauss Prop 3.1's exponent bookkeeping; without the balance restriction it is `A_5(delta)`.)

**What the literature gives at `e`** (all uniform, all legitimate):

| input | exponent obtained (w.r.t. `-K`, balanced points) | status |
|---|---|---|
| any section / Seshadri / effective-threshold / Liouville-with-base-locus argument (McKinnon-Roth 2016; Huang-type torsor integrality) | `<= 15/4`, exactly `15/4` at best | proved cap (Prop 1.3(b); = record row 20) |
| same, all distinct tuples, `-K` sections of level `m <= 5` | `10/3` (Vandermonde) = Cilleruelo-Córdoba | exact (§1.3) |
| Roth-type (McKinnon-Roth 2015, Ru-Vojta, Heier-Levin, Grieve) | `<= tau <= 15/4` (rational target) | proved cap |
| McKinnon conjecture / Coba / weak form, proved cases | statements about `alpha_e <= 2`; no bound on `alpha_ess` | inapplicable |
| Manin / equidistribution (Batyrev-Tschinkel anisotropic tori; Huang split) | none at sub-natural scale | inapplicable |
| Vojta `V_E(5)` | `>= 4` | conjecture |

**The gap** is therefore one specific phenomenon: *an unconditional proof that, at a rational point `x` of a variety
with `-K` ample, the essential approximation constant exceeds the effective threshold `tau_x(-K)` (on a Zariski-dense
set of a prescribed shape).* I found **no theorem of this kind for any variety** (searches: McKinnon-Roth family,
McKinnon conjecture and all proved cases, Huang's local distribution, Batyrev-Tschinkel, intrinsic approximation on
quadrics, offline mirror 2024-03..2026-10 for "approximation constant", "essential constant", "local distribution",
"best approximation", "McKinnon"). In every computed example found or checked here (split dP6, `P^1 x P^1`, twisted
`k = 3, 4`), `alpha_ess` equals the section/Liouville value; `k = 5` twisted is predicted to be the first case where it
does not. Two structural reasons why the existing machinery cannot reach it:
1. **Rational target.** Roth/Schmidt/Faltings-Wüstholz gains require the target to be irrational (Galois-conjugate
   components of `D`); at a rational `e`, with `E` geometrically irreducible, they are weaker than Liouville.
2. **Non-split torus, but no extra rational curves.** The anisotropic `T_5` has 1-dimensional `Q`-subtori (conics
   `C_v`), but with five distinct coordinates they approximate at exponent `>= 8`, worse than Dirichlet (§1.2). So
   the "curve" mechanisms (McKinnon) give nothing in the relevant range, and the dense family that kills `k = 4`
   (gauss Thm 4.1) has no analogue (Thm 4.3(ii), killer's dimension count).

Owner's constraints: the needed statement and all candidate tools above have exceptional sets and constants fixed by
`(X, -K, e)`; none assumes uniform exceptional sets from fixed-`S` Subspace theorems, discards gcd costs, zeroes a
phase, or uses a template-dependent threshold.

## 9. Papers worth fetching in full

Priority order (all full texts were unreachable this round: `arxiv.org`, `ar5iv.arxiv.org`, `ar5iv.labs.arxiv.org`,
`mast.queensu.ca`, `smf.emath.fr`, `uwspace.uwaterloo.ca`, `sites.google.com` do not resolve in this session; these
hosts would have to be allowed in the session's network settings):

1. **McKinnon-Roth, Eur. J. Math. 2 (2016), arXiv 1306.2977.** Exact Liouville and stable-base-locus statements; check
   whether anything beyond `tau_x(L)` is proved off a closed set (expected: no; Prop 1.3 would then be the complete
   translation).
2. **McKinnon-Roth, Invent. Math. 200 (2015), arXiv 1306.2976.** Exact form of the Roth-type inequality and of the
   "expected order of vanishing" constant; confirm `<= tau` at rational points.
3. **Huang, ANT 15 (2021), arXiv 1809.02001**, and **Acta Arith. 189 (2019), arXiv 1703.01772 (+ II, 1805.03920).**
   Check the mechanism of the lower bounds off accumulating subvarieties (section-type or not), and why Pell-Fermat
   equations appear on the toric surface `Y_4` (a quadratic-irrational phenomenon inside a split variety might be the
   closest analogue of our anisotropic situation).
4. **Batyrev-Tschinkel, IMRN 1995 (anisotropic tori)** and J. Algebraic Geom. 1998. The adelic Fourier analysis on
   `T(A)/T(Q)` is the only existing machinery written for anisotropic tori; check whether local Fourier transforms at
   `inf` near `e` could give counts in balls shrinking faster than the natural scale (needed: absolute error `< 1` in
   balls of radius `B^{-1/4-delta/2}`, §5 item 3).
5. **Huang, arXiv 2111.01509** (effective errors in adelic neighbourhoods, split) and **Huang-de Saxcé, arXiv
   2510.09954** (flag varieties, rational `x` under non-degeneracy): size and uniformity of error terms in shrinking
   neighbourhoods.
6. **Lehmann-McKinnon-Satriano, arXiv 2403.02480**; **McKinnon-Satriano, TAMS 374 (2021), arXiv 2004.05212 (+ erratum
   2023)**; **Monahan-Satriano, arXiv 2308.11847**: the exact definition of "essentially bounded" and of the
   essential constant; their Vojta computations at `Bl_P X` (to cross-check gauss Thm 4.2's `kappa = 2delta`).
7. **McKinnon, JAG 16 (2007), arXiv math/0504303**: Conjecture 2.7 verbatim and the examples where `alpha_ess` is
   computed.
8. Lower priority: McKinnon arXiv 2310.01640 (cubic hypersurfaces); Manon-McKinnon-Satriano arXiv 2504.05213;
   Huang IJNT 13 (2017), arXiv 1610.02542 (split dP6 = split `Lbar_3`; direct comparison with twisted `k = 3`);
   Kleinbock-Merrill (Israel J. Math. 2015) and Fishman-Kleinbock-Merrill-Simmons (JEMS 2022).

## Evidence labels
[V] = verified from source text seen in this round (none: no full text was reachable); [A] = abstract read in the
offline arXiv RSS mirror (`lit18/tools/db.json`, 2024-03-05..2026-10-27); [S] = from a web-search snippet only;
[R] = recalled, not re-verified; [D] = derived here with proof; [E] = exact computation (`lit_checks/`; ranks mod the
prime 2147483629, which give valid upper bounds for the rational maximal vanishing order; the lower bound `10` at
`m = 3` is the explicit Vandermonde).

## Sources consulted (URLs as returned by search; none fetched successfully)

* McKinnon, JAG 2007: https://www.ams.org/journals/jag/2007-16-02/S1056-3911-06-00458-9 ; https://arxiv.org/abs/math/0504303
* McKinnon-Roth, Invent. Math. 2015: https://arxiv.org/abs/1306.2976
* McKinnon-Roth, Eur. J. Math. 2016: https://link.springer.com/article/10.1007/s40879-016-0113-5 ; https://arxiv.org/pdf/1306.2977
* McKinnon-Satriano, TAMS 2021: https://www.ams.org/tran/2021-374-05/S0002-9947-2021-08318-X/ ; https://arxiv.org/abs/2004.05212
* Lehmann-McKinnon-Satriano: https://arxiv.org/pdf/2403.02480 ; Monahan-Satriano: https://arxiv.org/pdf/2308.11847
* McKinnon, cubic hypersurfaces: https://arxiv.org/pdf/2310.01640 ; Manon-McKinnon-Satriano: https://arxiv.org/html/2504.05213v1
* Huang, ANT 2021: https://msp.org/ant/2021/15-2/p05.xhtml ; https://arxiv.org/pdf/1809.02001
* Huang, IJNT 2017: https://arxiv.org/abs/1610.02542 ; Acta Arith. 2019: https://arxiv.org/pdf/1703.01772 ; II: https://arxiv.org/pdf/1805.03920
* Huang 2111.01509, Huang-de Saxcé 2510.09954, Huang-Wilsch 2509.03998 (mirror abstracts)
* Batyrev-Tschinkel: https://arxiv.org/abs/alg-geom/9510014 (cites the IMRN 1995 anisotropic-tori paper)
* Kleinbock-Merrill: https://arxiv.org/abs/1301.0989 ; Fishman-Kleinbock-Merrill-Simmons: https://ems.press/journals/jems/articles/3108905
* Grieve, expectations and Roth's theorem: https://arxiv.org/pdf/1705.01023
