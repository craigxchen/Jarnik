# Referee report on round12/t1rigour.md (T = 1 rigour)

## 0. Verdict

**The restricted theorem survives.** No growth bound and no uniform bound is claimed, and none follows. Every
proved statement on which the conditional reduction (Theorem 6.8) rests checks out. That covers §2 (T = 1 slab box),
§3 (rigidity with `W_1`), §4 (Step 0), Props 5.1-5.3, Lemmas 6.1-6.2, Theorem 6.3, Cor 6.4, Prop 6.5, Lemma 6.6,
Theorem 6.7, Lemmas 6.9-6.10 and Theorem 6.11(b). I re-derived all of them and recomputed every numerical input
with my own scripts. One claimed result is not supported, and three statements need correcting:

* **Major (F1).** Prop 5.5 ("scales compatible for `n >= 5185`") rests on `o_F <= m' - 2js - 2e`, i.e.
  `rho_A(F) >= 2(js + e)`. That inequality is not proved, and it is **false**: near-twin sign columns violate it
  (exact computation, §5). The author's script actually uses the weaker exponent `m' - d_F`, not the `o_F` of their
  own ceiling (5.1). With the lower bound on `rho_A(F)` that can be proved, compatibility is established (double
  precision) only for `n >= ~2.55 * 10^5` (`r = n`, `C = 1`). The bulleted `n >= 5185` claim in §0 item 4 must be
  withdrawn or relabelled as conditional. No theorem uses Prop 5.5, so Theorem 6.8 is unaffected.
* **Minor (F2).** The exact model of Theorem 6.11(a) is not zonotope-free (`zeta(A) >= 8 >= m`). It is also outside
  the regime (`m < 4w`) and uses a `w` that does not come from the formula. It is a toy only. The regime claim rests on
  6.11(b), which I confirm, including that it is zonotope-free, overdetermined and LP-supercritical.
* **Minor (F3).** §7.1 "J1 at T = 1": the coordinate family is not the cheapest family whose span contains J1's
  `p`-adic point. The single direction `1_{V_c}` gives cost `<= 4j + 2 floor((lambda m' + 4j)/4)`. So with
  `lambda = 1` the point lies in a cheap span inside the window from `n = 2^13`, not `2^15`. The headline "not cheap at
  `m' = ceil(4w)`" survives for `lambda >= 1`, which is the case that matters.
* **Structural (F4).** ZE-TC^rat is trivially true at every rational zero, since `R u` is a rational subspace of the
  cone. (Q1)+(Q2) can be replaced by **one** hypothesis on a **finite** family: quantitative ZE-TC relative to the
  collision subspaces with at least `binom(n,2)/Q_0` coincident pairs. This hypothesis is implied by (Q1)+(Q2) and
  suffices for Theorem 6.8 (§6). Definition 5.4 requires `𝓛` finite, but (Q2) applies it to a possibly infinite
  family; the reformulation removes that defect as well.
* **Minor (F5).** Prop 6.5(c) ("iff `k >= 16`") is checked by the author only for `k <= 22`. The tail is a one-line
  estimate (supplied in §6).

I found no cheap subtorus that Theorem 6.3 misses. Its width lower bound needs only `S_w ⊂ A`, which holds for
every profile. Twin and near-twin columns do not create new cheap spans, but they do make the projection cost
`rho_A` strictly smaller than any width-type count. That is exactly where F1 breaks.

## 1. Setup and what was checked

I read t1rigour.md in full, together with [D] §§1-4, R11 README §§0, 5, 8, compute.md §§1-2, theory.md §§1-3 and the
relevant referee reports. I did not import the author's code. My scripts are in `referee_t1rigour_checks/` and each
writes `*_out.txt`:

| script | what it checks | type |
| --- | --- | --- |
| `r_params.py` | `m_max` for every `n` in `[3,12000]` by big-integer comparison `8 binom(m+d-1,d) <= 3^d` at `m` and `m+1`; `m_max > c_* d - kappa`; `m_max <= 2d`; last failures of `m_max >= 4w` (`m^2 >= 128 d ln 4r` at 50 digits) for `r = n, n ln n, n^2`; interval values of `g(c_*) - ln 3`, `phi(316)`, `6 ln 3/sqrt 8`; the `316 ln 4r` thresholds; `c_0` | exact / interval |
| `r_sylv.py` | `m_max` at `n = 2^12..2^16` (log-gamma, then big-integer confirmation); Sylvester `2jk` against the window; tail `k >= 17` | exact |
| `r_rho.py` | `rho` by exact integer rank of the evaluation matrix: near-twin bands, the projection example, Lemma 6.9 for `s <= 6` | exact |
| `r_toy.py` | Theorem 6.11(a): `|A|`, `dim V(A,m)` for `m = 6..9`, common zeros by exact rank augmentation, and an explicit zonotope in `A` | exact |
| `r_consistency.py` | Prop 5.5's `G`, with the author's exponent and with the proved one | double precision (evidence) |
| `r_j1T1.py` | the indicator identity behind F3 (all `(a,b,c,x)`, `k = 3..6`); cheapness windows of `(b, a+b, 1_{V_c})` | exact / integer |

Results:
* All exact numbers of §2 agree: last failures `2899 / 3622 / 5808`, closed-form thresholds `2966 / 5930`, and
  `m_max(3000) = 1928`, `m_max(4096) = 2631`. The values `m_max(2^13..2^16) = 5260, 10518, 21034, 42065` also agree.
* The interval constants agree: `g(c_*) - ln 3 = -8.84097e-5`, `phi(316) = 1.660318`, `kappa = 2.213663`,
  `c_0 = 0.6417941`.
* The toy dimensions `193 / 96 / 33 / 6` and the stated zero pattern agree.

## 2. Section 2 (T=1 slab-box theorem): correct

I re-derived each part.

* **(a)** In `X_p = sum_{x<=d} k_x(beta_{p,x} - beta_{p,n})` the terms are independent and mean zero, each with
  range `<= 4`. Hoeffding then gives `2 exp(-w^2/(8d)) = 1/(2r)`.
* **(b)** The single-term binomial bound is right, and so is concavity of `g` on both sides. The upper-bound
  direction uses `g(c) - g(c_0) >= g'(c)(c - c_0)` for concave `g` with `c >= c_0`, which is correct. The lower bound
  needs `d >= 2` for `m_max` to exist (`n = 2` has no admissible `m`); this is harmless.
* **(d)** The Siegel count is right, and the condition `L >= 3 ln 3` is exact.
* **(e)** `phi' > 0` for `t > 32/c_*^2`; `phi(316) = 1.66032`.
* **(f)** Trivial.
* The exact thresholds and the isolated success at `5807` (`r = n^2`) are reproduced.
* The labels are correct: `2920 / 3642 / 5830` is double-precision evidence, and no proof uses it.

## 3. Section 3 (dichotomy, rigidity): correct

* **Lemma 3.1.** `D - |k|_1` is even because `|k|_1 ≡ sum k ≡ 0 (mod 2)`. Also `G_k = N^{D/2} e^{i<k,delta>}`, and
  `|<k',u>| <= d Theta`. (I1) is [D] §4 line 301ff and is stated for arbitrary supports with `|k|_1 <= D`.
* **Theorem 3.1.** `|c|_1 <= U^3 <= 3^{3d}` and `cost <= wW`. With `m'/2 - w >= m'/4` and `ln m'! >= m' ln(m'/e)` this
  gives (3.2) exactly.
* **Corollary 3.2.** I spot-checked `W_1^ex(3000, 1899) = 53.27` by hand. `W_1(m_max) <= 52.70 + 8 ln C'` at
  `d = 2899`, and the bound decreases in `d`.
* **Remark 3.3.** The margin `m'W/5` and the codimension bound `(36.97 + 14.65 d)/(m'W)` are both correct.
* **Quantifiers.** The theorem holds for every profile and every order in `[ceil(4w), m_max]`, with `H = U^2`,
  `T = 1`, any base point and any norm. Squarefreeness, normalisation and nondegeneracy are indeed not used.
  Distinctness of the points enters only through `u != 0` (§4).

## 4. Section 4 (Step 0): correct

* Lemma 4.1 is theory.md Lemma 2.4 with `D_{hat u} <= d`.
* For `|u|_inf >= e^{-W/2}` one needs `d_{xn} <= sum_p e_p log p <= W`, which holds.
* The log-rooms `6178` (at `m_max`) and `6198` (at `ceil(4w)`, `n = 3000`) are reproduced by hand.

## 5. Section 5 (transversal vanishing): Props 5.1-5.3 correct; Prop 5.5 not supported (F1)

**Props 5.1-5.3 are correct.**
* *Prop 5.1.* If `i < m' - rho(Lambda)`, then `<k',u_L>^{m'-i}` agrees on `A` with a polynomial of degree
  `<= rho(Lambda)` in `y = F_*(k')`, so the total degree is `< m'`.
* *Prop 5.2.* Since `eta <= 1`, `sum_{i>=o} binom(m',i)(1+eta)^{m'-i} eta^i <= eta^o (2+eta)^{m'}`.
* *Prop 5.3.* The argument is a correct compactness argument. As stated, the constants are ineffective.

**F1: the inequality used in Prop 5.5 is false.** Prop 5.5 asserts that a family with parameters `(s, e)` has
`o_F <= m' - 2js - 2e`, where `o_F = m' - rho_A(F)` is the projection-cost exponent of (5.1). That is the claim
`rho_A(F) >= 2(js + e)`. Widths give only the opposite inequality, `rho_A(F) <= d_F`, together with
`d_F >= 2(js + e)`, so nothing in the paper proves it. The author's script `consistency.py` in fact uses the exponent
`m' - d_F`, the R11 exponent. That exponent gives a *weaker* valid ceiling, not the ceiling (5.1) that the text
describes.

*Counterexample (near-twin sign columns, any regime size).*
* **Setup.** Take two sign columns with sides `S_1`, `S_2`, neither containing the base point, with
  `|S_1 \ S_2| = |S_2 \ S_1| = delta < j` and `|S_1 ∩ S_2| >= 2j + 2 delta`.
* **The image is a band.** The image of `A` under `F = (1_{S_1}, 1_{S_2})` lies in the band
  `{(t_1,t_2) : |t_i| <= j, |t_1 - t_2| <= 2 delta}`. In the two-column model it equals the band (`r_rho_out.txt`:
  image = band for `(j,delta) = (2,1), (3,1), (4,2)`).
* **Upper bound on `rho`.** Newton interpolation across the levels of `t_1 - t_2`, ordered `0, 1, -1, ...`, gives
  `rho <= 2j + 2 delta`. Exact rank gives `rho = 2j + 2 delta` for all `1 <= delta <= j <= 6`.
* **Comparison with the claim.** For `(j, delta) = (6, 1)`, `rho = 14`. The family `F` has `(s, e) = (2, 0)`, so the
  claim is `>= 24`. The family `F' = (1_{S_1}, 1_{S_1} - 1_{S_2})` has the same span and `(s, e) = (1, 2 delta)`, so the
  claim is `>= 16`. Both claims fail.
* **The family is cheap in the regime.** At `m' = ceil(4w)` the family is cheap, since `d_F <= 4j <= 2w < m'`. Its
  transversal order is therefore `>= m' - 2j - 2 delta`, which is larger than either claimed value.

A further caution: `rho` is not monotone under projection. For `Lambda = {(0,0),(1,1),(2,0)}`, `rho(Lambda) = 1`
while its first-coordinate projection has `rho = 2` (`r_rho_out.txt`). So "one long direction forces
`rho_A(F) >= 2j`" also needs an argument. Here is one.

**Lemma R1 (proved lower bound for the projection cost).** For every family `F` and every profile,
`rho_A(F) >= 2 min(j, nu_F)`, where `nu_F = min(floor(n/2), n - (largest F-class))`. Consequently, at any order:
* if `F` has a long direction (`ell_i >= 1`), then `rho_A(F) >= 2j`;
* if all directions are sparse, then `rho_A(F) >= 2 min(j, |E|)`.

*Proof.*
1. Take `T = min(j, nu_F)` disjoint pairs `(x_i, y_i)` with `x_i` and `y_i` in different `F`-classes. Such a matching
   exists because a complete multipartite graph has a matching of size `min(floor(n/2), n - largest part)`.
2. The points `k_0 + sum_i t_i (e_{x_i} - e_{y_i})`, with `k_0 = -sum_i (e_{x_i} - e_{y_i})` and `t in {0,1,2}^T`, lie
   in `S_w`, hence in `A` (Lemma 6.1).
3. Their images form a lattice zonotope in `F_*(A)` with nonzero generators `g_i = F_*(e_{x_i} - e_{y_i})` and
   `sum L_i = 2T`. By compute.md Prop 2.1, `V(F_*(A), 2T) != 0`, i.e. `rho_A(F) >= 2T`.
4. *Long direction.* If `v_i` has `ell_i >= 1`, then its bottom `j` points and its top `j` points take disjoint
   values, so no class meets both sets and `n - (largest class) >= j`.
5. *Sparse directions.* If every direction is sparse, the core is one class of size `n - |E| > n/2`.

QED

**Consequence for Prop 5.5.** Rerunning `G` with the proved exponent gives
`o_F <= m' - 2j` (`s >= 1`) and `o_F <= m' - 2 min(j, |E|)` (`s = 0`). The radii and the `Q_{s,e}` are those of
Theorem 6.7. The results (`r_consistency_out.txt`, `r = n`, `C = 1`, double precision, step 50 above `2*10^4`):

| exponent used | last `n` with `W`-coefficient `<= 0` | last `n` with `G(3n ln n - 14n) <= 0` |
| --- | --- | --- |
| author (`m' - 2js - 2e`) | 5184 (reproduced) | 3661 (reproduced) |
| proved (Lemma R1) | ~255100 | ~250600 |

Asymptotically the binding family has `s = 4`, `o = m' - 2j ~ 3w` and `Q ~ 16`, so the condition is
`n > (48 + o(1)) w`, not `(8 + o(1)) w`.

**Corrected statement of Prop 5.5.** Scale compatibility of the Theorem 6.7 radii with the rigorous ceiling (5.1)
holds (evidence, double precision) for `n >= ~2.55 * 10^5` at `r = n`, `C = 1`. For smaller `n` it holds if
`rho_A(F) >= 2(js + e)` for the cheap families of the profile at hand. That condition fails for near-twin columns,
although for those the Lemma 6.6 radius is far better than Theorem 6.7's worst case. This affects no theorem, since
Theorem 6.8 is asymptotic in `n`. It does affect the round's reading that the `T = 1`, `ceil(4w)` reduction is
"not vacuous from `n >= 5185`".

## 6. Section 6 (cheap families, separation, refutation)

### 6.1 Theorem 6.3, Cor 6.4 and Prop 6.5 are correct; no cheap family is missed

**Lemma 6.1.** `|<k, beta_p>| <= |k|_1 <= 2j <= w` for every `beta_p in [-1,1]^n`, so `S_w ⊂ A` holds for every
profile. That makes (6.1) universal.

**Lemma 6.2.** An entry above the window sits among the top `j - 1` sorted positions and adds `>= 1` by
integrality. The `ell = 0` equality case is right.

**Theorem 6.3.** It follows by summation.

**Twins and the base point.** I looked specifically for ways a cheap family could escape the structure theorem:
* *Twin points and twin directions.* Pair directions `e_x - e_y` cost exactly 4 and appear in `E`. Twin points cannot
  occur in actual clusters (Lemma 1.1).
* *Duplicate columns.* They leave `A` unchanged.
* *Near-twin columns.* `(1_{S_1}, 1_{S_2})` costs `4j` and has 4 core classes, which is covered.
* *Base-point directions.* `e_n` has `tau_j = 1`, so the lower bound applies, although its true width can be large.

In every case the lower bound (6.1) needs only `S_w ⊂ A`, so nothing escapes.

**Cor 6.4.**
* `10j >= 5w - 10 >= 4w > m' - 1` for `w >= 10`, so `b <= 4`.
* `|E| <= (m'-1)/2 < 2w`.
* `w >= sqrt(2528) ln 4 = 69.7` and `w <= 0.159 d`.
* "Attained" is correct: four sign columns with both sides `>= j` cost `8j < m'`, unless `w` is an even integer.

**Prop 6.5.**
* (a) is pigeonhole.
* (b) holds because `x -> (1[x_i=0])_i` is injective.
* (c) is reproduced: `2jk = 26340 > m_max = 21034` at `2^15`, and `40896 < 42065` at `2^16`.
* **F5.** The author checks "iff `k >= 16`" only for `k <= 22`. For `k >= 17` use `2jk <= kw` and
  `m_max > 0.6417 d - 2.2137`. The inequality `k sqrt(8 d ln 4d) < 0.6417 d - 2.2137` holds for `17 <= k < 80`
  (`r_sylv_out.txt`). Beyond that, `k^{3/2} 2^{-k/2}` decreases, so it holds for all `k >= 17`.

### 6.2 Lemma 6.6, Theorem 6.7 and the reduction (Theorem 6.8): correct, with a better formulation

**Lemma 6.6 and Theorem 6.7.**
* Lemma 1.1 gives `Norm gcd(z_x, z_y) = N e^{-d_xy}` on every norm: inert primes and `1+i` are common to all points.
* The averaging over collision pairs is valid because `eps_xy >= 0`.
* The sup-norm factor `1/2` is right.
* Convexity gives `P >= n'(n'-q_F)/(2q_F)`.
* `Q_F <= 16 n(n-1)/((n-2w)(n-2w-16))` holds at `ceil(4w)`.
* Only distinctness of the points and the pair gap are used. "Any cluster, any norm" is right.

**Theorem 6.8.** The implication is valid. Two remarks change how it should be stated.

* **(F4a) Finiteness.** Definition 5.4 takes `𝓛` finite, but (Q2) applies it to "all rational subspaces contained
  in `TC_H`". That set can be infinite: a real cone can contain infinitely many rational lines.
* **(F4b) ZE-TC^rat is degenerate at rational points.** `TC_H` is a cone, so for every rational zero `u` the line
  `R u` is a rational subspace inside `TC_H`. ZE-TC^rat therefore says nothing about rational zeros; it only
  constrains irrational ones. (Q1) correspondingly says: **every rational point of `TC_H` has
  `>= binom(n,2)/Q_0` coincident pairs.**
* **Cleaner equivalent hypothesis.** Let
  `𝓛_coll(Q_0) = {Delta_pi : pi a partition of [n] with >= binom(n,2)/Q_0 within-block pairs}`, where
  `Delta_pi = {u : u_x = u_y` for `x ~_pi y}` and `u_n = 0`. This is a **finite** set of rational subspaces.
  * If (Q1) holds, every `L` in `𝓛_A` lies in the member `Delta_{pi(L)}` of `𝓛_coll`, since
    `P(Delta_{pi(L)}) = P(L)`. A point `sigma`-far from all of `𝓛_coll` is then `sigma`-far from all of `𝓛_A`. So
    (Q1)+(Q2) imply ZE-TC_eps`(A, ceil(4w); 𝓛_coll(Q_0), sigma, eps_0)`.
  * Conversely, Theorem 6.8 applied with `𝓛 = 𝓛_coll(Q_0)` (Lemma 6.6 gives
    `sigma_L >= (1/2)e^{-Q_0 abar/2}`) needs nothing else.
  * The single hypothesis to state is therefore **ZE-TC^coll_eps**. Every `hat u` with `|hat u|_inf = 1` at sup-distance
    `>= (1/2)e^{-Q_0 abar/2}` from every collision subspace with `>= binom(n,2)/Q_0` coincident pairs has
    `kappa(hat u) > eps_0`.
  * Its qualitative core is **ZE-TC^coll**: every real common zero of the `Phi_c`, `c in V_H(A, ceil(4w))`, has at
    least `binom(n,2)/Q_0` coincident pairs.
  * It is consistent with everything proved: cheap spans (Theorem 6.7), `R^{S'}` (Theorem 6.11) and the uncut cube
    (`TC = {|supp u| < m'/2}`) all have `>= c n^2` coincidences. It is a quantitative braid statement and is not
    touched by the ZE-braid refutation, because that needs `m' > 2j log_2 n` (Prop 6.5).

### 6.3 Lemmas 6.9-6.10 and Theorem 6.11

**Lemma 6.9.** The Newton argument is correct; I checked the slice induction, including the dropped levels. Exact
values `rho(Lambda_t(s)) = s - |t|` and `rho(A2(s,j)) = s + j` hold for `2 <= s <= 6` (`r_rho_out.txt`).

**Lemma 6.10.** Correct.

**Product formula for one cut.** Correct, reading `J` as the leading-form ideal. `in(J_1) + in(J_2)` and
`in(J(A_1 x A_2))` have quotients of the same dimension `|A_1||A_2|`, and the bidegree argument is right.

**Theorem 6.11(a) (F2).**
* The statement is true: `R^{S'} ⊂ TC(A,m)` for `m = 6, 7, 8`, and no family of width cost `< m` spans it. It is
  confirmed exactly.
* But the model is **not faithful**. `A` contains the lattice zonotope with generators `e_1 - e_2`, `e_3 - e_4`,
  `e_5`, `e_6` (each `L = 2`, `sum L = 8`) (`r_toy_out.txt`). So `zeta(A) >= 8 >= m`, and by compute.md Prop 2.1
  every `(A, m)` with `m <= 8` is settled by a character. Consistently, `R^{S'}` lies in the hyperplane `u_5 = 0` of
  that certificate.
* At `m = 9` the model is faithful, but there `R^{S'}` is a cheap span (cost `8 < 9`).
* Also `w in [2,4)` is not `sqrt(8 d ln 4r) = 8.16`, and `m < 4w`.
* So (a) is a toy illustration only.

**Theorem 6.11(b): confirmed, and the model is faithful.**
* *Sizes.* `w = 563.83`, `j = 281`, `m' = 2256`, `b = 4` and `m_max = 2631` are all exact.
* *In the tangent cone.* `s + j <= 2255 < m'`.
* *In no cheap span.* A cheap span containing `R^{S'_0}` would need `s <= |E| + q_F - 1 <= 1127`. (The author's
  `<= 1128` is slightly lossy but still sufficient.) The cheap spans form a finite set and `R^{S'_0}` is irreducible,
  so a generic point of `R^{S'_0}` is a common zero outside every cheap span.
* *Faithful.* I checked the three properties the brief asks for:
  * LP-supercritical: `m' = ceil(4w) <= m_max`.
  * Overdetermined: `dim V_H >= |A|/2` (Thm 2.1(d)).
  * **Zonotope-free.** For a lattice zonotope `Z ⊂ A` with generators `v_i`, sum `h_Z(beta_p) <= 2w` over the
    `n - 1` Sylvester columns. For `v in Z^n_0 \ 0`, `sum_{a != 0} |<v, chi_a>| >= n |v|_2^2/|v|_1 >= n` (Parseval).
    Hence `sum L_i <= 2w(n-1)/n < 2w < m'`, i.e. `zeta(A) < m'`.
* *Caveat.* "Admissible" means R11's SL\* admissibility: rank `n - 1` and pair separation `(n-1)/(2n)`. The actual
  pair gap `d_xy >= W/2 - 2 ln C` would need `ln C >~ tau/4`, with `tau >= 2198` forced by `q_0 = 5`. Since
  ZE-TC is a statement about profiles, the refutation stands.

The cheap-span form of ZE-TC is indeed false at the recommended order.

## 7. Section 7 (J1-J3)

**J2.** The Chebyshev argument is correct. I recomputed by hand: `n = 2137`, `D/2 = d n`,
`W_TC = 4.661e6 = 1.0207 n^2`, factor 242.3. Confirmed.

**J3.** Confirmed by `r_params.py`: `m_max(3000) = 1928`, thresholds `2900 / 3623 / 5809`, `W_1 = 51.75` (exact form).

**J1 at `T = n`.** Re-derived, with no change. (The author's exact `m(2^13) = 49388999` comes from their 60-digit
log-gamma plus a big-integer boundary check; I did not recompute it.)

**J1 at `T = 1` (F3).**
* The author costs the coordinate family `(b, a+b, e_x : x in V_c)` at `4j + n/2`. That figure is correct for that
  family.
* But the indicator identity `1[a.x=c_1] 1[b.x=c_2] = (1/4)(1 + eps_1 chi_a + eps_2 chi_b + eps_1 eps_2 chi_{a+b})`
  (exact check, `k = 3..6`, 0 failures) bounds the single direction `1_{V_c}` on the `T = 1` box `A_{Q,lambda}`:
  `|k(V_c)| <= (lambda m' + 4j)/4`.
* So `(b, a+b, 1_{V_c})` costs at most `4j + 2 floor((lambda m' + 4j)/4)`. Its span contains `a'_p`.
* With `lambda = 1` this family is cheap for `m' in [4943, 5260]` at `n = 2^13`, `[7223, 10518]` at `2^14` and
  `[10535, 21034]` at `2^15` (`r_j1T1_out.txt`).
* So J1's point lies in a cheap span inside the window **from `n = 2^13`**, not `2^15`, and "cheap iff
  `m' > n/2 + 4j`" holds only for the coordinate family.
* At `m' = ceil(4w)` it is not cheap for `lambda >= 1`, so the author's conclusion there stands. For
  `lambda < 1/2` it would be cheap, but that is not a real widening (`h(a_p) < m'`).

## 8. Issues list

1. **(major, F1)** The Prop 5.5 premise `o_F <= m' - 2js - 2e` (i.e. `rho_A(F) >= 2(js+e)`) is unproved and false
   for near-twin sign columns, where exactly `rho_A = 2j + 2 delta < 2j + 4 delta <= 4j`. The "`n >= 5185`" scale
   compatibility in §0 item 4, Prop 5.5 and §8 item 2 is not established. With the proved Lemma R1 it holds
   (evidence) from `n ~ 2.55e5` (`r = n`, `C = 1`). The asymptotic condition becomes `n > (48+o(1))w`. No theorem is
   affected.
2. **(minor, F2)** The Theorem 6.11(a) toy has `zeta(A) >= 8 >= m`, so it is not zonotope-free and has to be read as
   an illustration, outside the regime. Part (b) carries the claim and is faithful.
3. **(minor, F3)** §7.1 "J1 at T = 1": a cheaper family exists via `1_{V_c}`. Cheapness inside the window starts at
   `n = 2^13` (`lambda = 1`), and the stated "iff" holds only for the coordinate family. The conclusion at
   `m' = ceil(4w)` (`lambda >= 1`) stands.
4. **(minor/structural, F4)**
   * Definition 5.4 requires `𝓛` finite, while (Q2) uses a possibly infinite `𝓛_A`.
   * ZE-TC^rat is vacuous at rational zeros.
   * (Q1)+(Q2) should be replaced by the single finite-family hypothesis ZE-TC^coll_eps (§6.2). It is implied by
     (Q1)+(Q2) and suffices for Theorem 6.8.
5. **(minor, F5)** Prop 6.5(c) needs the tail `k >= 23`. It is supplied in §6.1.
6. **(cosmetic)**
   * Theorem 6.11(b) can use `s <= |E| + q_F - 1`, which gives the range `1128 <= s <= 1974`.
   * The heuristic sentence about random sign profiles (§6.3, end) should be labelled as heuristic.

No error was found in §§2-4, Props 5.1-5.3, Lemmas 6.1-6.2, 6.6, 6.9, 6.10, Theorems 6.3, 6.7, 6.8 (as an
implication), 6.11(b) or the J2/J3 rechecks. Every "exact" or "interval" number I recomputed agrees.

## 9. Corrected statement

**Status: restricted theorem (survives).** No growth bound and no uniform bound; `M <= (2/3+eps) log R/loglog R`
stands. Proved for every profile, every norm and any base point, at `T = 1` with `w = sqrt(8(n-1) ln 4r)`,
`j = floor(w/2)`, `H = U^2`:

1. **Slab box.** `|A| >= 3^{n-1}/2`; `m_max > 0.6417(n-1) - 2.2137`; `dim V_H(A,m') >= |A|/2` for all `m' <= m_max`
   (`n >= 5`); `m_max >= 4w` if `n - 1 >= 316 ln 4r`. Exact thresholds are `2900 / 3623 / 5809`.
2. **Rigidity.** For every `m'` in `[ceil(4w), m_max]` and `W > W_1 = 24(n-1) ln 3/m' + 8 ln(e(n-1)C'/m')`.
3. **Step 0.** `kappa_c(hat u) <= C(n-1)e^{-W/4}/(m'+1)`.
4. **Transversal vanishing.** Order `>= m' - rho_A(F)` along the span of any family with `rho_A(F) < m'`. Ceiling
   `min(2^{1-m'}, (2+eta)^{m'} eta^{m'-rho_A(F)})`. New lower bound `rho_A(F) >= 2 min(j, nu_F)` (Lemma R1).
5. **Structure of cheap families.** For any order: total window length `<= floor((m'-1)/(2j))`, `|E| <= (m'-1)/2`,
   `<= 2^b` core classes. At `ceil(4w)`: `<= 4` long directions and `<= 16` classes, and cheap spans lie in the
   braid arrangement. At Sylvester orders `m' > 2j log_2 n` this fails, within the window exactly for `n = 2^k`,
   `k >= 16`.
6. **Separation.** Every cluster on every norm satisfies `dist(hat u, L) >= (1/2)exp(-binom(n,2) abar/(2P(L)))`,
   and hence `>= (1/2) exp(-16(1+O(w/n)) abar/2)` from every cheap span at `ceil(4w)`.
7. **Refutation.** The cheap-span form of ZE-TC is false for an admissible, zonotope-free, overdetermined,
   LP-supercritical profile (Sylvester `2^12` plus one light column, `m' = ceil(4w) = 2256`).
8. **Conditional reduction.** Hypothesis ZE-TC^coll_eps`(A, ceil(4w); Q_0)` implies `n = O(log W) = O(loglog R)` on
   all norms. This holds for any function `Q_0(n) >= 1`, since Lemma 6.6 supplies the radius for whatever `Q_0` is
   chosen; scale consistency, however, needs `Q_0 = O(1)` (cf. §5). The hypothesis says: for all profiles with
   `n >= n_1` and `n - 1 >= 316 ln 4r`, and all `W > W_1`, `kappa(hat u) > eps_0` at every `hat u` at distance
   `>= (1/2)e^{-Q_0 abar/2}` from all collision subspaces with `>= binom(n,2)/Q_0` coincident pairs. Its qualitative core: every real common zero has `>= binom(n,2)/Q_0`
   coincident pairs. (Q1)+(Q2) of the write-up imply it.
9. **Withdrawn / downgraded.** Scale compatibility from `n >= 5185` is withdrawn. With proved exponents it holds
   (double-precision evidence) from `n ~ 2.55*10^5` (`r = n`, `C = 1`). "J1 never cheap at `T = 1`" holds at
   `m' = ceil(4w)` for genuine widenings (`lambda >= 1`). Inside the window, J1's point becomes cheap from
   `n = 2^13` (not `2^15`).
