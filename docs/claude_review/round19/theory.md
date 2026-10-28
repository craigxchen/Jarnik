# Round 19 — theory: can a pencil on Mbar_{0,8} have positive abc-gain?

Status: complete; unrefereed. Labels: **proved** = full proof below from standard, cited facts; **evidence** = finite computation not used in any proof.

## 0. Summary and verdict

**Verdict: the abc lane is closed by theory (proved, unrefereed).** Every unit pencil `f = lambda A/C` on
`Mbar_{0,8}` (A, C products of differences, of any degree `D`, any `lambda ∈ C^*`) has

```text
gain(f) = sum_Y ((m_Y - 3)/3) Y.beta  <=  -2^{8-k}  <=  -1,
```

where `k` is the number of the 8 points the unit `A/C` involves (Theorem 3.3). So round-18 Prop 2.6 (abc + a unit
pencil with `gain > 0` gives `(L_8^pc)`) never applies. The same holds for every composite `h(A/C)` (Corollary 3.4).
The round-18 judge's stop rule fires, in a stronger form: the bound covers every `D` and every `lambda` at once, so
no `D = 5` search is needed. Only the sign of the gain enters Prop 2.6. The bound `-1` is weaker than the observed
`-32` but is enough.

**The new tool (Theorem 2.4, proved).** For every non-constant rational function `f` on `X = Mbar_{0,8}`, with
fibres `F_p ≡ L` and `R_p^int` the `beta`-degree of the reduced non-boundary part of `F_p`:

```text
sum_{p ∈ P^1} (L.beta - R_p^int)  <=  2 L.beta + 31 .
```

This is a `beta`-weighted log Riemann–Hurwitz inequality, with equality for the cross-ratio. Proof:
1. `df` spans a line subsheaf `O(div df)` of `Omega^1_X(log ∂)`.
2. `Omega^1_X(log ∂)` has a filtration by pull-backs along forgetful maps. Its line-bundle quotients have
   `beta`-degrees `16, 24, 28, 30, 31` (Lemma 2.1, exact).
3. `beta` is strictly positive on effective divisors (`threshold/L8.md` Thm A(a), via Keel–McKernan, refereed).
4. So `(div df).beta <= 31`, and `div df` is computed from the fibres (Lemma 2.3).

The consequence is `gain(f) <= 31 - n'(f)`, where `n'(f)` counts the boundary divisors on which `f` is constant
(Corollary 2.5). For a unit pencil the zero and polar fibres contain `N(u) >= 32` boundary divisors. This is
Schwartz–Zippel for the multilinear quadric `x ↦ ord_{D_x}(u)` on `{0,1}^7`, sharp for cross-ratios (Lemma 3.2).
The margin is exactly `32 - 31 = 1`.

**Where an abc lane could still open (§4; open).** Only through non-unit pencils. That needs two things, neither
of which exists:
* (i) a gcd extension of Prop 2.6 to pencils with base locus in codimension-2 boundary strata (an unrefereed sketch
  of the obstruction is in §4.1);
* (ii) a pencil of positive gain.

Such a pencil must be constant on at most 30 of the 119 boundary divisors. It must not be composite with a unit, must
involve at least 5 points, and must have three non-reduced fibres whose `beta`-average multiplicities satisfy
`sum 1/mu_p < 1` (§4.2). None is known, and Theorem 2.4 alone allows `gain <= 31` for such pencils.

**Unchanged:** `(L_8)` is open. The uniform bound is conditional on Vojta. Growth: `M <= (2/3 + eps) log R/loglog R`.
No claim here concerns actual circles.

## 1. Setup: the gain functional, beta, and what must be proved

**Notation.** `X = Mbar_{0,8}` over `C` (everything below is geometric; a pencil defined over `Q` is treated
over `C`, which changes no number: Galois-conjugate components have equal multiplicity and equal `beta`-degree).
`∂ = sum_T D_T` is the boundary, 119 smooth components meeting with simple normal crossings (Knudsen, Keel).
`M = X \ ∂ = M_{0,8}`. `beta` is the balanced curve class: the linear functional on `N^1(X)_Q` with
`D_T.beta = 1` for all 119 `T`. It is well defined because the `D_T` span `N^1(X)` and every Keel relation has 16
divisors on each side (checked: `theory_checks/c1_classes.py`, 1680 relations, 0 violations). In Kapranov's model
(anchor 8, coordinates `x_1..x_7` mod translation and scaling): `H.beta = 31`, `E_T.beta = 1`, `K.beta = 10`.

**Pencils and fibres.** A *pencil* is a non-constant rational function `f` on `X`. For `p ∈ P^1` let `F_p = f^*(p)`
be the fibre divisor: `F_∞` is the polar divisor of `f` and `F_p = div(f - p) + F_∞` for `p ∈ C`. All `F_p` are
effective, pairwise without common components, and linearly equivalent to one class `L` (the moving part; the base
locus has codimension `>= 2` and plays no role below). Write `F_p = sum_Y m_Y Y` over prime divisors and

```text
R_p      = (F_p)_red . beta                     (reduced fibre)
R_p^int  = (sum of the components Y ⊄ ∂ of F_p, reduced) . beta
n_p      = #{T : D_T ⊂ supp F_p},      n(f) = n_0 + n_1 + n_∞,
n'(f)    = #{T : f|_{D_T} is constant}  (>= n(f)).
```

So `R_p = R_p^int + n_p`. The round-18 gain (theory.md §2.4, referee_theory.md §2.6) is

```text
gain(f) = L.beta - (R_0 + R_1 + R_∞) = sum_{Y ⊂ F_0+F_1+F_∞} ((m_Y - 3)/3) Y.beta .
```

For a unit pencil `f = lambda A/C` this is exactly round 18's `L.beta - F_red.beta` (its `F_red` counts the `D_T`
with `max(nu_T(A),nu_T(B),nu_T(C)) > mu_T` and the interior factors `Y` of `B` with weight
`31 deg Y - sum nu_T(Y) = Y.beta`); the identity `gain = 10 - K^log_{X~/P^1}.beta~` of round 18 holds as well.

**Units.** `O^*(M)/C^*` is free of rank 20. In the Kapranov chart `M` is the complement of the 21 hyperplanes
`x_i = x_j` in `P^5 = P(C^7/C·1)`, so every unit is `u = c prod_{1<=i<j<=7} (x_i - x_j)^{e_ij}` with `e ∈ Z^21`,
`sum e_ij = 0` (degree 0), and

```text
ord_{D_T}(u) = e(T) := sum_{i<j, i,j ∈ T} e_ij        (T ⊂ [7], 2 <= |T| <= 6; these are the 119 divisors).   (1.1)
```

(For `|T| = 2` this is the exponent of `x_i - x_j`; for `|T| >= 3` it is the order of vanishing of the product along
the centre `L_T` of `E_T`, which is round 18's `nu_T`.) A *unit pencil* is `f = lambda u`, `lambda ∈ C^*`, `u` a
non-constant unit; for it `F_0`, `F_∞ ⊂ ∂` and `n_0 + n_∞ = N(u) := #{T : e(T) ≠ 0}`.

**What has to be proved (stop rule of round18/README.md §7.3).** Round 18 Prop 2.6: under abc, a unit pencil with
`gain > 0` gives `(L_8^pc)`. The judge's stop rule asks for `gain <= -32` on the untested finite families
(`D = 5`, all `lambda`; `D = 4`, `lambda ≠ ±1`). For closing the lane what matters is the sign: Prop 2.6 needs
`gain > 0`. This write-up proves `gain <= -1` for **every** unit pencil, of every degree, for every `lambda ∈ C^*`
(Theorem 3.3), so the hypothesis of Prop 2.6 is never satisfied.

**Standard facts used** (all classical; each is quoted with its source):
* (F1) `X` is smooth projective and `∂` is snc (Knudsen; Keel 1992).
* (F2) The forgetful map `pi: Mbar_{0,n+1} -> Mbar_{0,n}` is the universal curve (Knudsen). Near a node of a fibre it
  is, in suitable local analytic coordinates, `(u, v, y) -> (t = uv, y)`, where `{t = 0}` is the boundary divisor of
  `Mbar_{0,n}` whose general point has that node; `{u = 0}` and `{v = 0}` are the two boundary divisors of
  `Mbar_{0,n+1}` over it. Near a marking it is `(x, y) -> y` with the section `{x = 0}` a boundary divisor
  `D_{i,n+1}`; elsewhere it is smooth with `∂_{n+1} = pi^{-1}(∂_n)` locally.
* (F3) `beta` lies in the interior of the movable cone: `E.beta > 0` for every non-zero effective divisor `E`.
  This is `threshold/L8.md` Theorem A(a) (refereed, `threshold/referee_L8.md` item 3), proved there from
  Keel–McKernan, *Contractible extremal rays on Mbar_{0,n}*, Thm 1.3 (the effective cone of `Mbar_{0,n}/S_n` is
  spanned by `B_2, ..., B_{floor(n/2)}`), using only that `beta` is `S_8`-invariant and `D_T.beta > 0`.
* (F4) `K_{Mbar_{0,k}} = sum_T (|T|(k-|T|)/(k-1) - 2) D_T` (Keel–McKernan, Pandharipande), as used in round 18.

## 2. The log-cotangent filtration and a beta-weighted Riemann–Hurwitz inequality

### 2.1 The filtration

Fix an ordering of the 8 marked points. For `4 <= k <= 8` let `pi_k : X -> Mbar_{0,k}` forget the last `8 - k`
points (`pi_8 = id`) and put

```text
F_k := pi_k^* Omega^1_{Mbar_{0,k}}(log ∂_k)  ⊂  F_8 = Omega^1_X(log ∂),        F_3 := 0.
```

**Lemma 2.1 (proved).** For `4 <= k <= 8`, `F_{k-1} ⊂ F_k` is a subbundle, `F_k` is locally free of rank `k - 3`,
and `Lambda_k := F_k / F_{k-1}` is a line bundle with

```text
c_1(Lambda_k) = pi_k^*(K_k + ∂_k) - pi_{k-1}^*(K_{k-1} + ∂_{k-1})     (with K_3 + ∂_3 := 0),
Lambda_4.beta = 16,   Lambda_k.beta = 32 - 2^{8-k} (k = 5..8):   16, 24, 28, 30, 31.
```

*Proof.* (i) Let `pi : Y' = Mbar_{0,m+1} -> Y = Mbar_{0,m}` forget the last point. The boundary of `Y'` is
`pi^{-1}(∂_Y)` (reduced: `pi^* D_S = D_S + D_{S ∪ {m+1}}`) together with the `m` sections `D_{i,m+1}`. In the three
local charts of (F2):
* near a node, with coordinates `(u, v, y)` on `Y'` and `(t = uv, y)` on `Y`: `Omega^1_{Y'}(log ∂)` has the frame
  `du/u, dv/v, dlog y_j` or `dy_j` (according as `{y_j = 0}` is or is not boundary), and the image of
  `pi^* Omega^1_Y(log ∂)` is spanned by `dt/t = du/u + dv/v` and the same `y`-forms;
* near a marking, `(x, y) -> y` with the section `{x = 0}`: frame `dx/x` plus the `y`-forms;
* elsewhere `pi` is smooth and `∂_{Y'} = pi^{-1}∂_Y` locally: frame `dx` plus the `y`-forms.

In each chart the image of `pi^*Omega^1_Y(log ∂)` is spanned by part of a frame. So
`0 -> pi^*Omega^1_Y(log ∂) -> Omega^1_{Y'}(log ∂) -> Q -> 0` is exact and locally split, with `Q` invertible
(generated by `du/u`, `dx/x`, `dx` respectively). (This is the log-smoothness of the universal curve; `Q` is the
relative log dualizing sheaf `omega_pi(sum sigma_i) = psi_{m+1}`.)
(ii) Pulling the sequence for `Mbar_{0,k} -> Mbar_{0,k-1}` back along `pi_k` keeps it exact and locally split; its
first term is `F_{k-1}` and its middle term is `F_k`. Induction from `F_4 = pi_4^* O_{P^1}(-2 + 3)` gives the ranks.
(iii) `c_1(Omega^1(log ∂)) = K + ∂` on each `Mbar_{0,k}`; take determinants.
(iv) For a boundary divisor `D_S` of `Mbar_{0,k}`, `pi_k^* D_S` is the reduced sum of the `2^{8-k}` divisors `D_T`
of `X` with `T ∩ [k] = S` up to complement, so `(pi_k^* D).beta = 2^{8-k} D.beta_k`, `beta_k` the balanced class
of `Mbar_{0,k}`. With (F4), `(K_k + ∂_k).beta_k = (k-4) 2^{k-3} + 1`, hence
`(pi_k^*(K_k + ∂_k)).beta = 32(k-4) + 2^{8-k}` = 16, 40, 68, 98, 129, and the differences are as stated
(`theory_checks/c1_classes.py`, exact). `[]`

The same numbers come from `Lambda_k = pi_k^* psi_k`: `psi_k` on `Mbar_{0,k}` is `sum_{k ∈ T; a,b ∉ T} D_T`, of
`beta_k`-degree `2^{k-3} - 1`, and `2^{8-k}(2^{k-3} - 1) = 32 - 2^{8-k}`; `psi_j.beta = 31` is checked directly.

**Lemma 2.2 (line subsheaves of `Omega^1(log ∂)` have `beta`-degree `<= 31`; proved).** Let `Delta` be a divisor on
`X` and `O_X(Delta) -> Omega^1_X(log ∂)` an injective sheaf map. Let `k` be the least index with image in `F_k`.
Then `Lambda_k - Delta` is linearly equivalent to an effective divisor, and

```text
Delta.beta <= Lambda_k.beta <= 31.
```

*Proof.* By minimality the composite `O(Delta) -> F_k -> Lambda_k` is non-zero (its kernel would be contained in
`F_{k-1}`). A non-zero map of invertible sheaves on an integral scheme is a non-zero section of
`Lambda_k ⊗ O(-Delta)`, so `Lambda_k - Delta ~ E >= 0`, and `E.beta >= 0` by (F3). `[]`

### 2.2 The divisor of `df`

**Lemma 2.3 (proved).** Let `f` be a pencil and `theta = df`, a non-zero rational section of
`Omega^1_X(log ∂)`. For a prime divisor `Y` let `ord_Y theta` be the least order along `Y` of the coefficients of
`theta` in a local frame at the generic point of `Y`. Then

```text
div(theta) = sum_{c ∈ C} ( F_c - (F_c^int)_red )  -  F_∞  -  (F_∞^int)_red ,                     (2.1)
```

where `F^int` is the sum of the non-boundary components (with multiplicity), and the sum over `c` is finite.

*Proof.* Let `s` be a local equation of `Y` at its generic point and `y` complementary coordinates. Frames:
`ds, dy` if `Y ⊄ ∂`; `ds/s, dy` if `Y ⊂ ∂`.
* `Y ⊄ ∂`, `Y ⊂ F_c` (`c ∈ C`) with multiplicity `m`: `f - c = s^m g`, `g` a unit, so
  `df = s^{m-1}(m g ds + s dg)`, order `m - 1` (char 0).
* `Y ⊂ ∂`, `Y ⊂ F_c` (`c ∈ C`), multiplicity `m`: `df = s^m (m g ds/s + dg)`, order `m`.
* `Y ⊂ F_∞`, multiplicity `m`: `f = s^{-m} g`, `df = s^{-m}(-m g ds/s + dg) = s^{-m-1}(-m g ds + s dg)`: order `-m`
  if `Y ⊂ ∂`, `-(m + 1)` if not.
* `Y` in no fibre: `f` is regular at the generic point of `Y` and `f|_Y` is not constant (otherwise `Y ⊂ F_c`), so
  some `dy`-coefficient of `df` is non-zero on `Y`: order `0`.

Comparing coefficients with (2.1) gives equality. `[]`

Since `X` is normal and `Omega^1_X(log ∂)` is locally free, `1 -> theta` defines an injective map
`O_X(div theta) -> Omega^1_X(log ∂)` (it is regular in codimension one by Lemma 2.3, hence everywhere).

### 2.3 The inequality

For a pencil `f` let `k(f)` be the least `k` such that `f = g ∘ pi_S` for some forgetful map
`pi_S : X -> Mbar_{0,S}` with `|S| = k` (so `4 <= k(f) <= 8`).

**Theorem 2.4 (beta-weighted log Riemann–Hurwitz inequality; proved).** For every pencil `f` on `Mbar_{0,8}`:

```text
sum_{p ∈ P^1} ( L.beta - R_p^int )  <=  2 L.beta + Lambda_{k(f)}.beta  <=  2 L.beta + 31 .           (2.2)
```

Every term on the left is `>= n_p >= 0`, and only finitely many are non-zero.

*Proof.* Order the points so that `S` comes first. Then `theta = df = pi_S^* dg` lies in `F_{k(f)}`. Let
`k' <= k(f)` be the least index with `theta ∈ F_{k'}`. Lemma 2.2 with `Delta = div theta` gives
`(div theta).beta <= Lambda_{k'}.beta <= Lambda_{k(f)}.beta` (`Lambda_k.beta` increases with `k`). By (2.1) and
`F_p.beta = L.beta`,

```text
(div theta).beta = sum_{c ∈ C} (L.beta - R_c^int) - L.beta - R_∞^int ,
```

and adding `L.beta - R_∞^int` to both sides gives (2.2). Each term: `L.beta = F_p.beta = sum m_Y Y.beta >=
R_p^int + n_p`, because `Y.beta > 0` for every prime divisor (F3). For all but finitely many `p` the fibre `F_p` is
reduced (generic smoothness on the open set where `f` is a morphism) and contains no `D_T`, and then the term is 0.
`[]`

*Remarks.* (a) (2.2) is the analogue of the log Riemann–Hurwitz formula
`sum_p (deg f - #f^{-1}(p) ∩ C°) = 2 deg f + (2g - 2 + #(C \ C°))` for a curve `C ⊃ C°`. The role of the log
Euler characteristic is played by `Lambda_k.beta`. Indeed on `Mbar_{0,4} = P^1` with its three boundary points,
`2g - 2 + 3 = 1`, and `16 · 1 = Lambda_4.beta`.
(b) Equality holds for the cross-ratio `z` (`k = 4`): the fibres over `0, 1, ∞` are 16 boundary divisors each
and all other fibres are reduced and interior, so the left side is `3 · 16 = 48 = 2 · 16 + 16`. It also holds for
`4z(1 - z)`: `32 + 32 + 16 = 80 = 2 · 32 + 16`.
(c) The proof is not a curve argument. It uses the rank-5 bundle `Omega^1_X(log ∂)` and the positivity (F3) of
`beta` on all effective divisors. Round 18 showed that restriction to covering curves of class `d·beta` can give
only `gain <= (2g - 2)/d >= 10` (theory.md §2.4).

**Corollary 2.5 (proved).** For every pencil `f`,

```text
gain(f)  <=  Lambda_{k(f)}.beta  -  n(f)  -  sum_{p ∉ {0,1,∞}} (L.beta - R_p^int)  <=  31 - n'(f) .
```

*Proof.* `gain = L.beta - sum_{p ∈ {0,1,∞}} (R_p^int + n_p) = sum_{p ∈ {0,1,∞}} (L.beta - R_p^int) - 2 L.beta - n(f)`.
Bound the first sum by (2.2) minus the remaining terms. A boundary divisor on which `f` is constant with value
`c ∉ {0,1,∞}` lies in `F_c` and contributes at least 1 to the term at `p = c`. `[]`

## 3. Unit pencils: the abc lane is closed

**Lemma 3.1 (multilinear Schwartz–Zippel; proved).** A non-zero multilinear polynomial `P` of degree `<= d` in `n`
variables over a field is non-zero at `>= 2^{n-d}` points of `{0,1}^n`.

*Proof.* Induction on `n`; `n = 0` is clear. Write `P = x_n Q + R` with `Q, R` multilinear in `x_1..x_{n-1}`,
`deg Q <= d - 1`. If `Q = 0`, `P = R` and each of the two slices `x_n = 0, 1` has `>= 2^{n-1-d}` non-zeros. If
`Q ≠ 0` and `R = 0`, the slice `x_n = 1` is `Q`, with `>= 2^{(n-1)-(d-1)}` non-zeros. If `Q, R ≠ 0`: when
`Q + R ≠ 0` both slices give `>= 2^{n-1-d}`; when `Q + R = 0`, `R = -Q` has degree `<= d - 1` and the slice
`x_n = 0` alone gives `>= 2^{n-d}`. `[]`

**Lemma 3.2 (proved).** For every non-constant unit `u` of `M_{0,8}`, `N(u) = #{T : ord_{D_T}(u) ≠ 0} >= 32`.
Equality holds for the cross-ratios.

*Proof.* By (1.1), `ord_{D_T}(u) = q(1_T)` for the multilinear quadric `q(x) = sum_{i<j<=7} e_ij x_i x_j` on
`{0,1}^7`, with `e ≠ 0` because `u` is not constant. Lemma 3.1 gives `>= 2^{7-2} = 32` subsets `T ⊂ [7]` with
`q(1_T) ≠ 0`. None has `|T| <= 1` (`q = 0` there) or `T = [7]` (`q(1) = sum e_ij = 0`), so all are among the 119
divisors. For `(x_1 - x_2)/(x_1 - x_3)`, `q = x_1(x_2 - x_3)` has exactly 32 non-zeros
(`theory_checks/c4_support.py`). `[]`

**Theorem 3.3 (no unit pencil has positive gain; proved).** Let `u` be a non-constant unit of `M_{0,8}`, involving
exactly `k = k(u)` of the 8 points, and `lambda ∈ C^*`. Then

```text
gain(lambda u)  <=  Lambda_k.beta - N(u) - n_1  <=  -2^{8-k}  <=  -1      (and <= -16 if k = 4).
```

This holds for every degree `D` and every `lambda`, rational or not.

*Proof.* `F_0, F_∞ ⊂ ∂`, so `n(f) = N(u) + n_1`. Corollary 2.5 and Lemma 3.2 give
`gain <= Lambda_k.beta - 32`, and `Lambda_k.beta = 32 - 2^{8-k}` for `k >= 5`, `16` for `k = 4`. `[]`

**Corollary 3.4 (composites; proved).** If `f = h ∘ u` with `u` a non-constant unit and `h ∈ C(z)` non-constant,
then `gain(f) <= Lambda_{k(u)}.beta - N(u) <= -1`.

*Proof.* Let `G_t = u^*(t)`. Then `F_p(f) = sum_{t ∈ h^{-1}(p)} e_t G_t`, so `R_p(f) = sum_{t ∈ h^{-1}(p)} R_t(u)`
and `L_f = (deg h) L_u`. Put `S = h^{-1}{0,1,∞}` and `j = |S ∩ {0,∞}|`; Riemann–Hurwitz for `h` gives
`|S| >= deg h + 2`. Theorem 2.4 for `u`: the terms at `t = 0, ∞` equal `L_u.beta` (those fibres are boundary), so
`sum_{t ∈ C^*} (L_u.beta - R_t^int(u)) <= Lambda_k.beta`. Hence
`sum_{t ∈ S ∩ C^*} R_t(u) >= |S ∩ C^*| L_u.beta - Lambda_k.beta`, and with `R_0(u) = N_0`, `R_∞(u) = N_∞`
(the numbers of `D_T` with `e(T) > 0`, `< 0`):

```text
gain(f) <= (deg h - |S ∩ C^*|) L_u.beta + Lambda_k.beta - sum_{t ∈ S ∩ {0,∞}} N_t
        <= (j - 2) L_u.beta + Lambda_k.beta - sum_{t ∈ S ∩ {0,∞}} N_t .
```

Since `L_u.beta = sum_{e(T)>0} e(T) >= N_0` and likewise `>= N_∞`, each of `j = 0, 1, 2` gives
`gain <= Lambda_k.beta - N_0 - N_∞ = Lambda_k.beta - N(u)`. `[]`

**What this does to the round-18 lane.** Prop 2.6 of round 18 needs a unit pencil `lambda A/C` with `gain > 0`.
By Theorem 3.3 none exists, in any degree, for any `lambda ∈ Q^*` (indeed `C^*`). The data (all `gain <= -32` for
`D <= 4`, `lambda = ±1`, and `D <= 3`, all special `lambda`) are consistent with this and are now superseded as
evidence. The bound `-1` is weaker than the observed `-32`. It is enough, because only the sign enters
Prop 2.6. **The abc lane, as round 18 defined it, is closed by theory.** The judge's stop rule therefore fires, in a
stronger form: no `D = 5` (or any `D`) search is needed.

## 4. Where an abc lane could still open: non-unit pencils

Theorem 3.3 covers exactly the pencils to which round-18 Prop 2.6 applies. Two things would be needed to reopen an abc
lane, and both are open.

**4.1 A gcd input for non-unit pencils (open; the following is an unrefereed sketch of the obstruction).** If `A` or
`C` has a factor `Y ⊄ ∂`, the proof of Prop 2.6 breaks at the step `gcd(A(x), C(x)) = prod d_T^{mu_T}` on
`Bal^pc`. The step can survive only when the base locus `Bs(f)` avoids the open strata that `p`-adically carry
balanced pairwise-coprime points. For `p ∤ prod d_T`, `x mod p` lies in `M(F_p)`, so the condition is
`Bs(f) ∩ M = ∅`. For `p | d_T` (one `T` only, by pairwise coprimality), `x` is `p`-adically close to `D_T°` and far
from the other `D_T'`, and the excess `v_p(gcd) - mu_T v_p(d_T)` is the `p`-adic valuation of the gcd of the local
equations of the strict transforms of `A` and `C` along `D_T`. That is bounded only if `Bs(f) ∩ D_T° = ∅`. So the
natural class is

```text
P° := { pencils f regular on X° := X \ ∪_{T ≠ T'} (D_T ∩ D_T') }       (base locus in codimension-2 strata).
```

This class contains the unit pencils and their composites `h ∘ u`. A version of Prop 2.6 for `P°`, with the same
gain formula (interior components of all three fibres counted with cofactor weight `Y.beta`), is plausible but
**not proved here**. Round 18's remark (theory.md §2.4, after Prop 2.6) makes the same observation without proof.

**4.2 The configuration a positive-gain pencil in `P°` (or any pencil) would need (proved necessary conditions).**
Write `E_p = L.beta - R_p = sum_{Y ⊂ F_p} (m_Y - 1) Y.beta` (the `beta`-excess of the fibre). By the definitions,
`gain > 0` if and only if `E_0 + E_1 + E_∞ > 2 L.beta`. A pencil with `gain(f) >= 1` must satisfy all of the following.
1. `n'(f) <= 30` (Corollary 2.5). So `f` is non-constant on at least 89 of the 119 boundary divisors, and the three
   fibres together contain at most 30 boundary divisors.
2. `f` is not of the form `h ∘ u`, `u` a unit (Corollary 3.4). In particular no two of `F_0, F_1, F_∞` are both
   supported on `∂`. (If two are, a Möbius transformation makes `f` a unit pencil.)
3. `k(f) >= 5`. For `k(f) = 4`, `f = h ∘ chi` with `chi` a cross-ratio. Then gain `= 16(deg h - #h^{-1}{0,1,∞})
   <= -32`, by Riemann–Hurwitz on `P^1` (the round-18 morphism computation, which applies verbatim).
   For `k(f) = k`, Corollary 2.5 requires `n'(f) <= 31 - 2^{8-k}`.
4. All three fibres are non-reduced, with `E_p > 2 L.beta - E_q - E_r > 0` for `{p,q,r} = {0,1,∞}` (since each
   `E_q < L.beta`). Equivalently, the `beta`-average multiplicities `mu_p = L.beta / R_p` satisfy
   `1/mu_0 + 1/mu_1 + 1/mu_∞ < 1`, a "hyperbolic" triple.
5. The total excess off `{0,1,∞}` is bounded: `sum_{p ∉ {0,1,∞}} (L.beta - R_p^int) <= 30 - n(f)`.

What is **not** known: whether any pencil on `Mbar_{0,8}` (in `P°` or not) has positive gain. Theorem 2.4 alone allows
`gain <= 31` when `n'(f) = 0`. A non-unit pencil with three heavily non-reduced fibres and almost no boundary in them
is the only configuration left, and the record contains no example of a pencil with three non-reduced fibres
other than composites `h ∘ chi`.

*Remark (Campana, not used).* `Mbar_{0,8}` is simply connected. For a pencil with connected fibres, the orbifold
fundamental group of `(P^1, gcd-multiplicities)` is a quotient of `pi_1(X) = 1`. So at most two fibres have
gcd-multiplicity `> 1`, and those two multiplicities are coprime. Condition 4 is compatible with this: a fibre of
gcd-multiplicity 1 can still have large `mu_p` if its reduced components have small `beta`-degree. This is why
Campana-type statements (which see the inf- or gcd-multiplicity) cannot replace Theorem 2.4, which sees the
`beta`-weighted multiplicity.

## 5. Routes tried that do not work, and why the log-cotangent route does

* **Covering curves** (round 18): a covering family of genus `g` and class `d·beta` gives only
  `gain <= (2g - 2)/d`, and `K.beta = 10 <= 2g - 2`. Free rational curves have `K.gamma <= -2`, so no positive
  combination of them is `beta`.
* **Complete intersections of nef classes.** Products of `psi`-classes are movable. They would certify (F3)
  directly, but `beta` is **not** in the cone they span. The `S_8`-symmetrized pairings with `(D_2, D_3, D_4)`-type
  divisors are `(6,0,0)`, `(120,15,0)`, `(90,0,16)`, `(720,180,96)` and `(360,180,144)` for the exponent types
  `4, 31, 22, 211, 1111`, against `(1,1,1)` for `beta`. The exact LP is infeasible (`theory_checks/c2_psi_cone.py`).
  `(K + ∂)^4` gives `(1379, 244, 150)`, also pair-heavy (`theory_checks/t2.py`, via the exact intersection engine
  `theory_checks/mbar.py`). So (F3) is used through Keel–McKernan, as in `threshold/L8.md` Thm A(a), not through an
  explicit curve family.
* **Campana / orbifold base:** see the remark in §4.2.
* **Log canonical bundle formula / weak positivity** (Kawamata, Viehweg, Fujino). The relative
  `K^log_{X~/P^1}` is weakly positive only when the general fibre has `kappa >= 0`. The general fibre of a pencil on
  the rational variety `X` need not, so these theorems give no lower bound on `K^log_rel.beta`. Theorem 2.4 gives
  `K^log_rel.beta~ = 10 - gain >= n'(f) - 21`.
* **Mason on lines** (round 18 item 4) sees only Kapranov degrees; the `E_T` are invisible to lines.

**Why the log-cotangent bundle works.** `Omega^1_X(log ∂)` contains `O(div df)` for every pencil. All five graded
pieces of its forgetful filtration have `beta`-degree `<= 31 = psi_j.beta = 2^5 - 1`, while the zero and polar fibres
of a non-constant unit contain at least `2^5 = 32` boundary divisors (Schwartz–Zippel on the Boolean cube).
The margin is exactly one.

## 6. Exact checks (`theory_checks/`)

All arithmetic is exact (Python integers, `Fraction`, python-flint factorisation over `Q`).
* `c1_classes.py` -> `c1_classes.out`:
  * `beta` respects all 1680 Keel relations (16 = 16);
  * `(K + ∂).beta_k = 1, 5, 17, 49, 129` for `k = 4..8`;
  * pulled-back degrees `16, 40, 68, 98, 129`, so `Lambda_k.beta = 16, 24, 28, 30, 31`;
  * `psi_j.beta = 31` directly;
  * `K.beta = 10 = psi.beta - 2·119`.
* `mbar.py`, `t_mbar.py`: an exact intersection engine on `Mbar_{0,n}` by boundary restriction. Tests:
  `D^2 = -1` and `psi^2 = 1` on `Mbar_{0,5}`; `int psi^a` multinomial on `Mbar_{0,6}`, `Mbar_{0,8}`; a
  `D_T.psi^a` value against the product formula.
* `c2_psi_cone.py` -> `c2_psi_cone.out`: `beta ∉ cone(psi^a)` (exact LP over all basic solutions); `t2.py`:
  `D_T.(K + ∂)^4`.
* `c4_support.py` -> `c4_support.out`:
  * `N = 32` for the cross-ratio (Lemma 3.2 is sharp);
  * exhaustive check on all 36330 `±1` vectors with at most 4 edges and sum 0, and 200000 random vectors: `min N >= 32`;
  * *evidence:* smallest `N` found for units involving exactly `k = 4..8` points is 32, 48, 40, 48, 50.
* `c3_gain_vs_n.py` (adapted from round-18 `ref_gain.py`; jsonl per orbit and chunk, resumable): every unit pencil
  `lambda A/C` with `D <= 2` (and `D = 3`, see below) and `lambda` in
  `{±1, ±2, ±1/2, ±3, ±1/3, ±4, ±1/4}`. It checks Corollary 2.5 in the form `gain + n(f) <= 31`.
  *Result (evidence; consistency test only, no proof uses it; `c3_summary.out`).*

  | `D` | `A`-orbits | `(A, C, lambda)` triples | max gain | max `gain + n(f)` |
  | --- | --- | --- | --- | --- |
  | 1 | 1 | 280 | -32 | 16 |
  | 2 | 3 | 8,260 | -32 | 16 |
  | 3 | 8 | 138,600 | -32 | 16 |

  Every triple satisfies `gain + n(f) <= 16 <= 31`. The maximum 16 is `Lambda_4.beta`. It is attained, for example,
  by `((x_1 - x_2)/(x_1 - x_3))^D` with `lambda = 1` (gain `-32`, `n = 48`), where (2.2) is an equality. Runs were
  split into chunks of 250 `C`'s per process, because python-flint leaks memory; `c3_D3_err.txt` lists the empty
  trailing chunks (`FAIL`), and every orbit was checked complete against its count. The maximum gain `-32` reproduces
  round 18 at `lambda = ±1`. It extends the referee's `D <= 3` special-`lambda` result to the 14 listed values of
  `lambda`.

## 7. Consequences

1. **The abc lane of round 18 is closed by theory (proved).** No unit pencil `lambda A/C` on `Mbar_{0,8}`, of any
   degree and for any `lambda ∈ C^*`, has positive gain (Theorem 3.3: `gain <= -2^{8-k} <= -1`). So round-18
   Prop 2.6 never applies, and abc yields no `(L_8^pc)` through unit pencils or their composites (Corollary 3.4).
   The judge's stop rule fires, more strongly than stated: no `D = 5` search is needed, and its outcome cannot
   matter. The nested-chain extension (§3.3 of the round-18 README, J2) is moot for this lane.
2. **What is not closed.** An abc lane through non-unit pencils would need both (i) a gcd extension of Prop 2.6 to
   the class `P°` (base locus in codimension-2 strata), which is not proved, and (ii) a pencil satisfying conditions
   1-5 of §4.2, of which none is known. Theorem 2.4 does not exclude (ii).
3. **Unchanged.** `(L_8)` is open. The uniform bound is known only conditionally (Vojta). The best growth bound is
   `M <= (2/3 + eps) log R/loglog R`. Nothing here bears on actual circles or on the rate.
4. **New tool, possibly reusable.** Theorem 2.4 is a `beta`-weighted log Riemann–Hurwitz inequality,
   `sum_p (L.beta - R_p^int) <= 2 L.beta + 31`, valid for every rational function on `Mbar_{0,8}`. It is
   sharp for the cross-ratio. It holds verbatim for any movable class `gamma` with `31` replaced by
   `max_k Lambda_k.gamma`.
