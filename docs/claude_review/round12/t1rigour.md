# Round 12: route 5 at box size T = 1 — rigour and cheap families

## 0. Summary of results

**No growth improvement and no uniform bound is proved.** `M <= (2/3 + eps) log R/loglog R` stands. What is proved
makes route 5 at box size `T = 1` rigorous up to one hypothesis, and changes that hypothesis.

1. **Closed-form `T = 1` versions of [D] Thm 2.1 and Thm 4.1 (§§2-3, proved).** For every profile:
   `|A| >= 3^{n-1}/2`; `m_max >= 0.6417(n-1) - 2.2137` (an `O(1)` loss; `m_max = 0.641794 n + O(log n)`);
   `dim V_H >= |A|/2` for every order `m' <= m_max` once `n >= 5`; `m_max >= 4w` whenever `n - 1 >= 316 ln 4r`, i.e. for
   `n >= 2966` (`r <= n`) and `n >= 5930` (`r <= n^2`); exact thresholds `2900 / 3623 / 5809` (`r = n, n ln n, n^2`) on
   `[2, 12000]`. Rigidity holds for **every** order `m' in [ceil(4w), m_max]` once
   `W > W_1 = 24 (n-1) ln 3/m' + 8 ln(e(n-1)C'/m')` (`<= 52.71 + 8 ln C'` at `m_max` for `n >= 2900`, `= O(sqrt(n/log r))` at
   `ceil(4w)`). No squarefree, normalisation or nondegeneracy hypothesis.
2. **Step 0 (§4, proved).** The actual direction has `kappa_c(hat u) <= eps_0 = C(n-1) e^{-W/4}/(m'+1)` for all
   `c in V_H`; the Chebyshev ceiling is `2^{1-m'}`, so the log-room is `W/4 - O(m')`.
3. **Transversal vanishing (§5, proved).** Along the span of any family of projection cost `rho_A(F) < m'` (every
   cheap family qualifies, `rho_A(F) <= d_F`) every `Phi_c` vanishes to order `>= m' - rho_A(F)` transversally, and
   `kappa <= min(2^{1-m'}, (2+eta)^{m'} eta^{m' - rho_A(F)})` at distance `eta`. A local converse (Prop 5.3) holds when the
   transversal leading forms have no common zero, with ineffective constants.
4. **Cheap families at T = 1 (§6, proved).** Every direction has width `>= 2 tau_j` (trimmed spread, `j = floor(w/2)`),
   and `tau_j >= j * (window length) + (#points outside the window)`. Hence **every** cheap family (any profile, any
   order) has total window length `<= b = floor((m'-1)/(2j))`, at most `(m'-1)/2` exceptional points, and its span
   lies in a subspace where the other points take `<= 2^b` values (Theorem 6.3).
   * At `m' = ceil(4w)`: `b <= 4`, so **at most four long directions** and **at most 16 classes**, both attained; the
     judge's sketch (R11 §5.3) is proved, including directions taking many values (Corollary 6.4).
   * **Braid arrangement: exact dichotomy.** Cheap spans lie in the braid arrangement whenever
     `2^b < n - (m'-1)/2` (always at `m' = ceil(4w)` in the regime); at Sylvester orders `n = 2^k` this holds iff
     `m' <= 2jk`, and for `k >= 16` the `log_2 n` coordinate columns are cheap and collision-free for
     `m' in (2jk, m_max]` (Proposition 6.5).
   * **The arithmetic separation lemma is proved** at `m' = ceil(4w)`: for every cheap span,
     `dist(hat u, L_F) >= (1/2) e^{-Q_F abar/2}`, `Q_F = 16(1 + O(w/n))`, `abar = W/(2(n-1)) + 2 ln C'`, i.e. `n^{-O(1)}`
     at `W = O(n log n)`, for **any** cluster (Lemma 6.6, Theorem 6.7). The scales are compatible with the transversal
     ceiling for `n >= 5185` (`r = n`, all `W >= 3n ln n - 14n`; Prop 5.5).
   * **The cheap-span form of ZE-TC is false** (Theorem 6.11): the tangent cone also contains *projection-cheap*
     rational subspaces. Exact model at `n = 7`; and an admissible regime profile (Sylvester `2^12` plus one light
     unbalanced column of size `1129..1974`) at `m' = ceil(4w)`. The correct core is ZE-TC^rat (common real zeros lie
     in rational subspaces of the tangent cone).
5. **Reduction (Theorem 6.8, conditional).** (Q1) every rational subspace of `TC_H(A, ceil(4w))` has `>= n^2/(2Q_0)`
   collision pairs, plus (Q2) the quantitative ZE-TC_eps relative to those subspaces, imply `n = O(log W) =
   O(loglog R)` on **all** norms; no rank lemma or nondegenerate subcluster is needed. (Q1) is proved for cheap spans
   only; (Q2) is open.
6. **J1, J2, J3 rechecked (§7): all confirmed** (J1's exact `m(2^13) = 49388999` and cost `4.70779e7`; J2's factor
   `242.3`; J3's thresholds and `W_1 = 51.75`). New: J1's family is never cheap at `T = 1`, `m' = ceil(4w)`.

| item (task) | status | where |
| --- | --- | --- |
| (1) Thm 2.1(a)-(d), Thm 4.1 at `T = 1`, closed forms, `W_1`, every order in the window | **proved** | §§2-3 |
| (1) `m_max >= 4w` thresholds | **proved** (`316 ln 4r`; `2966`, `5930`) + **exact** (`2900/3623/5809`) | Thm 2.1(e) |
| (2) Step 0, explicit `eps_0` | **proved** | Thm 4.2 |
| (3) transversal vanishing + quantitative version | **proved** (upper bound; local lower bound ineffective) | Props 5.1-5.3 |
| (4) at most four long directions, all cheap families characterised | **proved** | Thm 6.3, Cor 6.4 |
| (4) cheap spans in the braid arrangement for `n >> w` | **true** iff `2^b < n - (m'-1)/2` (always at `ceil(4w)`); **false** at Sylvester orders for `m' > 2j log_2 n` | Prop 6.5 |
| arithmetic separation from cheap spans | **proved** (any norm) | Thm 6.7 |
| ZE-TC in cheap-span form | **false** (exact model; admissible regime profile) | Thm 6.11 |
| (5) J1 | **confirmed** | §7.1 |
| growth or uniform bound | **not proved** | — |

## 1. Setting and notation at T = 1

Conventions (as in rounds 10-11). *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted
with an exact source. *Exact*: integer or exact-rational computation. *Evidence*: floating point or a finite
experiment that no proof uses. Scripts are in `t1rigour_checks/` (each writes `<name>_out.txt`). [D] is
`round10/determinant.md`; R11 is `round11/README.md`.

**Clusters.** `N >= 2`, `W = log N`, `R = sqrt N`; `n >= 2` distinct `z_x in Z[i]` with `|z_x|^2 = N` on an arc of length
`<= C sqrt R`; `C' = max(C, 1)`. `delta_x` is the angle of `z_x` minus the arc midpoint; the angular width is
`Theta <= C N^{-1/4} = C e^{-W/4}`. For a split `p = pi conj(pi)`, `p^{e_p} || N`: `a_{p,x} = v_pi(z_x)`,
`beta_p = (2 a_p - e_p 1)/e_p in [-1,1]^n`, `r` = number of split primes, `d_xy = sum_p log p |a_{p,x} - a_{p,y}|`.
Nothing below assumes `N` squarefree, normalised or odd.

**The T = 1 slab box.** Put `d = n - 1` and fix the base point `n`. Frequencies are `k in Z^n_0 = {sum k = 0}`, written
`k' = (k_1, ..., k_d)`, `k_n = -sum k'`. For `r >= 1` put

```text
w   = sqrt(8 d ln(4r)),      j = floor(w/2),      D = 2d,
B_1 = {k in Z^n_0 : k' in {-1,0,1}^d},     A = {k in B_1 : |<k, beta_p>| <= w for every split p},
m_max = the largest integer m with 8 binom(m + d - 1, d) <= 3^d,
U = |A|,   H = U^2,   E(m') = binom(m' + d - 1, d) = dim Poly_{<m'}(R^d).
```

For an order `m' >= 1`, `V(A, m')` is the space of `c in Q^A` with `sum_k c_k P(k') = 0` for every polynomial `P` of
degree `< m'`, and `V_H(A, m')` is the `Q`-span of the integer `c in V(A,m')` with `|c|_inf <= H`. The order is any
integer in the **admissible window** `[ceil(4w), m_max]`; the recommended choice is `m' = ceil(4w)`.
`Phi_c(u) = sum_k c_k <k',u>^{m'}` (`u in R^d`), `kappa_c(u) = |Phi_c(u)| / (|c|_1 d^{m'})` (note `D/2 = d`), and
`kappa(u) = sup_{c in V_H(A,m') \ 0} kappa_c(u)`.

*Positions.* `u_x = delta_x - delta_n` (`x <= d`), so `<k, delta> = <k', u>`; `u != 0` because the points are distinct;
`hat u = u/|u|_inf`. In `u`-coordinates the base point has coordinate `0`; it is convenient to write `u_n := 0`.

*Directions and widths.* An integer direction is a functional `k -> <k, v>` on `Z^n_0` given by `v in Z^n` modulo
`Z 1`; in `k'`-coordinates it is `k' -> <k', v' - v_n 1>`. Its width on a set `S` is
`h_S(v) = max_S <k,v> - min_S <k,v>`. A family `F = (v_1, ..., v_s)` of nonconstant integer directions has
*(width) cost* `d_F = sum_i h_A(v_i)` and is *cheap* for `(A, m')` if `d_F < m'` (R11 §1). Its span
`L_F = span_R{v_i}` is read in `delta`-space modulo `R 1`, equivalently in `u`-space as `{(v'_i - v_{i,n} 1)}`-span.
`dist(hat u, L) = min_{y in L} |hat u - y|_inf` (sup norm in `u`-space).

*Regularity.* For a finite `Lambda` in `Z^s`, `rho(Lambda)` is the least `e` such that every function on `Lambda` is
the restriction of a polynomial of degree `<= e` (so `V(Lambda, m) != 0` iff `m <= rho(Lambda)`). It is monotone:
`Lambda` contained in `Lambda'` implies `rho(Lambda) <= rho(Lambda')`. For a family `F` let
`F_*: k' -> (<k',v_i>)_i` and define the *projection cost* `rho_A(F) = rho(F_*(A))`. Since `F_*(A)` lies in the grid
`prod_i {values of <k,v_i> on A}`, tensor Lagrange interpolation gives `rho_A(F) <= d_F`.

**Imported, all refereed.** (I1) [D] Thm 4.1's proof: forced divisibility `Norm F_c(z) >= N^D prod_p p^{-h_{supp c}(a_p)}`
for integer `c` with `F_c(z) != 0`, and `<k,a_p> = (e_p/2)<k,beta_p>` on `Z^n_0`. (I2) theory.md Lemma 2.4 (Taylor
remainder) and Props 2.1-2.2 (tangent cone; cheap spans lie in the base locus). (I3) compute.md Prop 2.1
(zonotope certificate). Everything else is proved below. Two elementary facts are re-proved because they are
used quantitatively:

**Lemma 1.1 (pair gap).** For `x != y`: `e^{-d_xy/2} <= |delta_x - delta_y| <= Theta <= C e^{-W/4}`. Consequently
`d_xy >= W/2 - 2 ln C`.

*Proof.* `g = gcd(z_x, z_y)` has `v_pi(g) = min(a_x, a_y)`, `v_{conj pi}(g) = e_p - max(a_x, a_y)` at each split `p`, and
the full power at inert primes and at `1+i`, so `Norm g = N e^{-d_xy}`. As `g | z_x - z_y != 0`,
`N e^{-d_xy} <= |z_x - z_y|^2 = 4N sin^2((delta_x - delta_y)/2) <= N (delta_x - delta_y)^2`. QED

**Lemma 1.2 (Plotkin average).** `sum_{x<y} d_xy <= (n^2/4) W`, so the mean of `d_xy` over the `binom(n,2)` pairs is
`<= W/2 + W/(2(n-1))`.

*Proof.* For values `a_x in [0, e]`, `sum_{x<y} |a_x - a_y|` is convex in each `a_x`, so it is maximal at a
`{0,e}`-valued vector, where it is `e |S| |S^c| <= e n^2/4`. Sum with weights `log p` and use `sum_p e_p log p <= W`. QED

## 2. [D] Thm 2.1 (a)-(d) at T = 1 with closed-form constants

**Theorem 2.1 (T = 1 slab box).** Let `n >= 2`, `r >= 1`, and let `beta_1, ..., beta_r in [-1,1]^n` and weights
`omega_p >= 0` with `sum_p omega_p <= 1` be arbitrary (no balance, Plotkin, rank or squarefree hypothesis). With the
notation of §1:

* **(a)** `U = |A| >= 3^d/2`.
* **(b)** `m_max >= c_* d - kappa` with `c_* = 0.6417` and `kappa = ln 8/ln(1 + 1/c_*) < 2.2137`. Also
  `m_max <= 2d` and `m_max <= c_0 d + 1 + ln(3d/8)/ln(3/2)`, where `c_0 = 0.6417941...` is the root of
  `(1+c) ln(1+c) - c ln c = ln 3`. So `m_max = c_0 n + O(log n)`, and the lower bound loses only `O(1)`.
* **(c)** Every `k in B_1` has `|k|_1 <= D = 2d`; `h_A(beta_p) <= 2w` and `h_A(a_p) <= e_p w` for every `p`; for a sign column
  (`e_p = 1`) `<k, a_p>` is an integer in `[-j, j]` on `A`, so `h_A(a_p) <= 2j`.
* **(d)** If `d >= 4`, then for **every** order `m' in [1, m_max]`: `E(m') <= U/4` and `dim V_H(A, m') >= U/2`, with
  `H = U^2`.
* **(e)** If `n - 1 >= 316 ln(4r)`, then `m_max >= 4w`. In particular `m_max >= 4w` for all `n >= 2966` when `r <= n`,
  and for all `n >= 5930` when `r <= n^2`.
* **(f)** For every order `m' in [ceil(4w), m_max]` and every `c in V_H(A, m')`:
  `sum_p omega_p h_{supp c}(beta_p) <= 2w <= m'/2 <= ord(c)/2`.

*Proof.* (a) Draw `k'` uniformly from `{-1,0,1}^d`. For fixed `p`, `X_p = <k, beta_p> = sum_{x<=d} k_x (beta_{p,x} - beta_{p,n})`
is a sum of independent mean-zero terms with ranges of length `2|beta_{p,x} - beta_{p,n}| <= 4`. Hoeffding:
`P(|X_p| > w) <= 2 exp(-2w^2/(16 d)) = 2 exp(-ln 4r) = 1/(2r)`. A union bound over the `r` columns gives
`P(k notin A) <= 1/2`.

(b) *Binomial bound.* For integers `N >= K >= 0`, `binom(N,K) (K/N)^K ((N-K)/N)^{N-K} <= 1` (it is one term of the
binomial expansion of `(K/N + (N-K)/N)^N`; `0^0 = 1`). With `N = m - 1 + d`, `K = d` and `c = (m-1)/d` this reads
`binom(m+d-1, d) <= exp(d g(c))`, `g(c) = (1+c) ln(1+c) - c ln c` (expand `(m-1+d) ln(m-1+d) - d ln d - (m-1) ln(m-1)`).
So `8 binom(m+d-1,d) <= 3^d` holds whenever `g((m-1)/d) <= ln 3 - ln 8/d`. Now `g'(c) = ln(1 + 1/c) > 0` and
`g''(c) = -1/(c(1+c)) < 0`, so for `c <= c_*` concavity gives `g(c) <= g(c_*) - g'(c_*)(c_* - c)`. Interval
arithmetic gives `g(c_*) - ln 3 in [-8.8410e-5, -8.8409e-5] < 0` (`interval_constants_out.txt`). Hence every
integer `m` with `m - 1 <= c_* d - kappa` (where `g'(c_*) kappa = ln 8`) satisfies the defining inequality, and
`m_max >= floor(c_* d - kappa) + 1 > c_* d - kappa`.
*Upper bounds.* If `m >= 2d + 1` then `binom(m+d-1,d) >= ((m+d-1)/d)^d >= 3^d`, so `m_max <= 2d`. The largest term of
the same binomial expansion is the `K`-th, and there are `N + 1` terms, so `binom(N,K) >= exp(d g(c))/(N+1)`. Then
`8 binom <= 3^d` forces `g(c) <= ln 3 + ln((m+d)/8)/d <= ln 3 + ln(3d/8)/d`. Concavity at `c >= c_0` gives
`g(c) >= ln 3 + g'(c)(c - c_0)`, and `g'(c) >= g'(2) = ln(3/2)` for `c <= 2`; so `c - c_0 <= ln(3d/8)/(d ln(3/2))`.

(c) For `x <= d`, `|k_x| <= 1`, and `|k_n| = |sum k'| <= d`, so `|k|_1 <= 2d`. The width bounds are the slabs;
`<k,a_p> = (e_p/2)<k,beta_p>` on `Z^n_0`, and for `e_p = 1` it is the integer `k(S_p) = sum_{x in S_p} k_x`, of absolute
value `<= w/2`, hence `<= j`.

(d) `E(m') <= E(m_max) <= 3^d/8 <= U/4` by (a). Let `M` be the `E(m') x U` moment matrix (rows `k'^alpha`,
`|alpha| < m'`); its entries are in `{-1,0,1}` and `V(A, m') = ker M`. The `(H+1)^U` vectors of `{0..H}^U` map into
`[-UH, UH]^{E(m')}`, so some fibre has `>= (H+1)^U/(2UH+1)^{E(m')}` elements; differences with one fixed element are
distinct integer vectors of `V(A,m')` with `|.|_inf <= H`, so they lie in `V_H`. A subspace of dimension `delta` meets
`[-H,H]^U cap Z^U` in at most `(2H+1)^delta` points (it projects injectively to some `delta` coordinates). With
`L = ln U`, `H = U^2`: `ln(H+1) >= 2L`, `ln(2H+1) <= 2L + ln 3`, `ln(2UH+1) <= 3L + ln 3`, so

```text
dim V_H  >=  (2L U - (U/4)(3L + ln 3)) / (2L + ln 3)  =  U (1.25 L - 0.25 ln 3)/(2L + ln 3),
```

which is `>= U/2` iff `L >= 3 ln 3`. By (a), `L >= d ln 3 - ln 2 >= 3 ln 3` once `d >= 4`.

(e) Put `lambda = ln 4r >= ln 4` and `t = d/lambda`. By (b) it suffices that
`lambda phi(t) >= kappa`, `phi(t) = c_* t - 4 sqrt(8t)`. `phi` is increasing for `t >= 32/c_*^2 = 77.7`, and interval
arithmetic gives `phi(316) in [1.66031, 1.66032]`, so for `t >= 316`: `lambda phi(t) >= (ln 4)(1.66031) = 2.3016 > kappa`.
For `r <= n`, `n - 1 - 316 ln(4n)` is increasing for `n >= 316` and is `> 0` at `n = 2966`; for `r <= n^2`,
`n - 1 - 316 ln(4n^2)` is increasing for `n >= 632` and `> 0` at `n = 5930` (`thresholds_T1_out.txt`).

(f) `supp c` lies in `A`, so `h_{supp c}(beta_p) <= h_A(beta_p) <= 2w`; `sum omega_p <= 1`; `4w <= m'`; and
`ord(c) >= m'` for `c in V(A, m')`. QED

*Exact checks* (`params_T1.py`, exact big-integer comparison `8 binom(m+d-1,d) <= 3^d` at the boundary, all
`n in [2, 12000]`; output `params_T1_out.txt`).
* The closed-form bound (b) holds at every `n` in the range (0 violations); the slack `m_max - (c_* d - kappa)` is
  between 4 and 7.5 for `n in [500, 12000]`. `m_max/n` decreases from `0.6440` (`n = 1000`) to `0.64208`
  (`n = 12000`), consistently with `c_0 = 0.641794`.
* `m_max >= 4w` holds for every `n in [2900, 12000]` (`r = n`), every `n in [3623, 12000]` (`r = n ln n`) and every
  `n in [5809, 12000]` (`r = n^2`); the last failures are `2899`, `3622`, `5808`. (At `r = n^2` the condition holds at
  `5807`, fails at `5808`, and holds from `5809`: `m_max` is an integer step function.) This confirms J3(d).
* The closed form is nearly sharp: solving `c_* d - kappa >= 4 sqrt(8 d ln 4r)` directly gives `2920` (`r = n`),
  `3642` (`r = n ln n`), `5830` (`r = n^2`), checked to `n = 10^9` on a grid (floating point); the `K = 316` criterion
  of (e) gives `2966` and `5930` with a two-line monotonicity proof.

**Remarks.** (i) Compared with [D] at `T = n` (`n - 1 >= 300 ln 4r`, exact thresholds `2137/2679/4278`), the `T = 1`
box needs slightly more points (`316 ln 4r`; exact `2900/3623/5809`) but every other constant improves: the
coefficient bound is `B = 1`, `D = 2d` instead of `2dn`, and `ln U <= d ln 3`.
(ii) Only the column *directions* enter `A`; the weights enter only through `sum omega_p <= 1` in (f).

## 3. [D] Thm 4.1 at T = 1: the dichotomy with explicit W_1

**Lemma 3.1 (dichotomy on any T = 1 support).** Let `S` be any finite subset of `B_1`, `m' >= 1`, and `c` an integer
vector on `S` with `c in V(S, m')`, `c != 0`. For an actual cluster put `F_c(z) = sum_k c_k G_k(z)`,
`G_k(z) = N^{(D - |k|_1)/2} prod_x z_x^{k_x^+} conj(z_x)^{k_x^-} in Z[i]` (`D = 2d >= |k|_1` is even). If `F_c(z) != 0` then

```text
cost(c) := sum_{p split} log p * h_{supp c}(a_p)  >=  m'W/2 - 2 ln|c|_1 - 2m' ln(d C') + 2 ln m'! .        (3.1)
```

*Proof.* `G_k(z) = N^{D/2} e^{i<k,delta>}` because `sum k = 0`. Taylor: `|e^{ix} - sum_{t<m'} (ix)^t/t!| <= |x|^{m'}/m'!`, the
terms `t < m'` cancel against `c in V(S,m')`, and `|<k,delta>| = |<k',u>| <= |k'|_1 |u|_inf <= d Theta`; so
`|F_c(z)| <= N^{D/2} |c|_1 (d Theta)^{m'}/m'!`. Forced divisibility (I1): `Norm F_c(z) >= N^D e^{-cost(c)}`. Compare
`Norm = |.|^2`, use `ln Theta <= ln C - W/4`, and replace `C` by `C'` (an arc of length `<= C sqrt R` with `C < 1` also
has length `<= sqrt R`). QED

**Theorem 3.1 (rigidity at T = 1).** Let an actual cluster (any norm) have `n >= 5`, `r` split primes, and
`m_max >= 4w` (for example `n - 1 >= 316 ln 4r`, Theorem 2.1(e)). Let `m'` be **any** integer in `[ceil(4w), m_max]`
and suppose

```text
W  >  W_1(n, C, m')  :=  24 d ln 3 / m'  +  8 ln(e d C'/m').                                             (3.2)
```

Then `sum_{k in A} c_k e^{i<k,delta>} = 0`, i.e. `F_c(z) = 0`, for every `c in V_H(A, m')`. Equivalently
`delta in Z_P(m') := {delta : sum_k c_k e^{i<k,delta>} = 0 for all c in V_H(A,m')}`; and `f_delta|_A` lies in
`V_H(A,m')^perp`, a subspace of dimension `<= |A|/2` containing `Poly_{<m'}|_A`.

*Proof.* Let `c` be an integer vector of `V(A,m')` with `|c|_inf <= H = U^2`. Then `|c|_1 <= U^3 <= 3^{3d}`, and
`cost(c) <= sum_p log p (e_p/2) h_A(beta_p) <= w sum_p e_p log p <= wW` (Theorem 2.1(c)). If `F_c(z) != 0`, (3.1) gives
`(m'/2 - w) W <= 6d ln 3 + 2m' ln(dC') - 2 ln m'!`. As `m' >= 4w`, `m'/2 - w >= m'/4`, and `ln m'! >= m' ln(m'/e)`; so
`W <= 24 d ln 3/m' + 8 ln(e d C'/m') = W_1`, contradicting (3.2). So `F_c(z) = 0` on the integer generators of `V_H`, and
by linearity on `V_H`. The dimension statement is Theorem 2.1(d). QED

**Corollary 3.2 (explicit sizes).** `W_1(n,C,m')` is decreasing in `m'`, so `W > W_1(n, C, ceil(4w))` suffices for every
admissible order. With `lambda = ln 4r`:

```text
W_1(n, C, ceil(4w))  <=  2.3306 sqrt(d/lambda) + 8 ln(e sqrt(d/(128 lambda))) + 8 ln C',
W_1(n, C, m_max)     <=  24 ln 3 / (c_* - kappa/d) + 8 ln(e/(c_* - kappa/d)) + 8 ln C'   ( -> 52.63 + 8 ln C' ).
```

(`24 d ln 3/(4w) = (6 ln 3/sqrt 8) sqrt(d/lambda)` and `6 ln 3/sqrt 8 in [2.33050, 2.33051]`; for `m_max` use Theorem
2.1(b) and `m_max <= 2d < e d`.) The sharper exact form of the same proof,
`W_1^ex = [6d ln 3 + 2m' ln(dC') - 2 ln m'!]/(m'/2 - w)`, is tabulated in `thresholds_T1_out.txt` (`r = n`, `C = 1`):

| `n` | `m' = ceil(4w)` | `m_max` | `W_1^ex(ceil 4w)` | `W_1(ceil 4w)` (3.2) | `W_1^ex(m_max)` | `W_1(m_max)` (3.2) | `3n ln n - 14n` |
|---|---|---|---|---|---|---|---|
| 3000 | 1899 | 1928 | 53.27 | 53.30 | 51.75 | 52.55 | 30057 |
| 10^4 | 3683 | 6421 | 87.56 | 87.57 | 36.87 | 52.60 | 1.36e5 |
| 10^6 | 44112 | 641800 | 630.68 | 630.69 | 27.25 | 52.63 | 2.7e7 |

This confirms J3's `W_1 = 51.75` (`n = 3000`, maximal `m`) and J3b's `53` and `631` (`m' = ceil(4w)`). The last column is
the judge's illustrative near-extremal size (round 11); no proof here uses it.

**Remark 3.3 (larger heights; `V_H` is almost all of `V`).** The same proof works with `H' = e^{m'W/40}` in place of
`U^2` provided `W > 10 d ln 3/m' + 10 ln(e d C'/m')` (then `2 ln|c|_1 <= 2d ln 3 + m'W/20` and the margin over
`cost <= m'W/4` is `m'W/5`). By theory.md Lemma 1.1 with `B = 1`,
`codim(V_{H'} in V(A,m')) <= 40 (U ln 2 + (U/4) d ln 3)/(m'W)`, a relative codimension `<= 14.6 d/(m'W) + 37/(m'W)`
(as `dim V(A,m') >= 3U/4`), which is `O(1/(w ln n))` at `m' = ceil(4w)`, `W ~ 3 n ln n`. So the tangent cone may be taken
with respect to `V_{H'}`, which differs from `V(A, m')` by a vanishing fraction of conditions.

**Remarks.** (i) The hypotheses of Theorem 3.1 are profile-only (`m_max >= 4w`) plus `W > W_1 = O(sqrt(n/log r) + log C')`.
Every cluster with `n >= n_0(C)` satisfies the second by `two_thirds/proof.md` (`W >= (3 - o(1)) n ln n`). No
squarefree, normalisation or nondegeneracy hypothesis enters.
(ii) [D] Prop 4.2 and Cor 4.3 carry over verbatim: `Z_P(m')` is a closed real-analytic Lebesgue-null set, so almost
every fake violates (4.2) at `T = 1`.

## 4. Step 0 at T = 1 with explicit eps

**Lemma 4.1 (Taylor at the actual point).** Let `c in V(A, m')` (real), `u in R^d \ 0` with
`sum_k c_k e^{i<k',u>} = 0`, and `hat u = u/|u|_inf`. Then `kappa_c(hat u) <= d |u|_inf/(m'+1)`.

*Proof.* theory.md Lemma 2.4 with `t = |u|_inf` and `D_{hat u} = max_A |<k', hat u>| <= |k'|_1 <= d`:
`t^{m'} |Phi_c(hat u)|/m'! <= |c|_1 (t d)^{m'+1}/(m'+1)!`. Divide by `t^{m'} |c|_1 d^{m'}/m'!`. QED

**Theorem 4.2 (Step 0 at T = 1).** Under the hypotheses of Theorem 3.1 (any admissible order `m'`), the actual
`hat u` satisfies, for every `c in V_H(A, m')`,

```text
kappa_c(hat u)  <=  eps_0  :=  C d e^{-W/4} / (m'+1),       and   e^{-W/2} <= |u|_inf <= C e^{-W/4}.         (4.1)
```

*Proof.* Theorem 3.1 gives `delta in Z_P(m')`, i.e. the hypothesis of Lemma 4.1 for every `c in V_H`. Then
`|u|_inf = max_x |delta_x - delta_n| <= Theta <= C e^{-W/4}`; the lower bound is Lemma 1.1 with `d_xn <= W`. QED

**How small eps_0 is (J2/J3 recheck).** The Chebyshev bound of J2 (re-proved in §7.2) gives
`kappa_c(hat u) <= 2^{1-m'}` for every `hat u` and `c`. So Step 0 carries information only if `eps_0 < 2^{1-m'}`, i.e.
`W > W_0 := 4(m'-1) ln 2 + 4 ln(C d/(m'+1))`. At `m' = ceil(4w)` this is `W_0 = 16 w ln 2 + O(log n) = O(sqrt(n log r))`;
at `m' = m_max` it is `~ 1.78 n`. In either case the log-room `ln(2^{1-m'}/eps_0) = W/4 - (m'-1) ln 2 - ln(Cd/(m'+1))`
is `(3/4 - o(1)) n ln n` at near-extremal `W ~ 3 n ln n`: `6178` at `n = 3000` (`m_max`, confirming J3) and `6198` at
`m' = ceil(4w) = 1899`. The quantitative tangent-cone hypothesis therefore only has to beat the ceiling by a factor
`e^{-(3/4 - o(1)) n ln n}`.

## 5. Transversal vanishing along cheap spans

Throughout, `F = (v_1, ..., v_s)` is a family of integer directions with **projection cost** `rho_A(F) < m'`
(§1). Every cheap family qualifies, since `rho_A(F) <= d_F`. Write `o_F = m' - rho_A(F) >= m' - d_F >= 1`, and
`L = L_F`. In `k'`-coordinates `v_i` becomes `v'_i = (v_{i,x} - v_{i,n})_{x<=d}`.

**Proposition 5.1 (exact transversal order).** For every `c in V(A,m')`, `u_L in L`, `v in R^d` and `eta in R`,

```text
Phi_c(u_L + eta v)  =  sum_{i = o_F}^{m'} binom(m', i) eta^i Psi_{c,i}(u_L, v),     Psi_{c,i}(u_L, v) = sum_k c_k <k',u_L>^{m'-i} <k',v>^i,
```

i.e. `Psi_{c,i}(u_L, .) == 0` for `i < o_F`: every `Phi_c` vanishes to order `>= o_F >= m' - d_F` transversally along `L`.
Moreover `Psi_{c,i}(u_L, v) = 0` for all `i` when `v in L`. The same holds with `Phi_c` replaced by the whole function
`t -> sum_k c_k e^{it<k',u_L + eta v>}`: the base locus contains `L` (R11 §1, theory.md Prop 2.2; here with `rho_A`).

*Proof.* Write `u_L = sum_i t_i v'_i` (`t` real) and `y = F_*(k') in Lambda = F_*(A)`; then `<k',u_L> = <t, y>`. If
`m' - i <= rho(Lambda)` the function `y -> <t,y>^{m'-i}` is itself a polynomial of degree `m' - i`; otherwise, by
definition of `rho`, it agrees on `Lambda` with a polynomial of degree `<= rho(Lambda)`. In both cases
`<k',u_L>^{m'-i} <k',v>^i` agrees on `A` with a polynomial in `k'` of degree `<= min(m'-i, rho(Lambda)) + i`, which is
`< m'` when `i < m' - rho(Lambda)`; `c` annihilates it. For `v in L`, `u_L + eta v in L` and the same argument with
`i = 0` gives `Phi_c == 0` on `L`. The base-locus statement is the same argument applied to `e^{i<t,y>}`. QED

**Proposition 5.2 (quantitative upper bound).** Let `|hat u|_inf = 1` and `eta = dist(hat u, L) <= 1`. For every
`c in V(A, m')`,

```text
kappa_c(hat u)  <=  min( 2^{1-m'},  (2 + eta)^{m'} eta^{o_F} ).                                         (5.1)
```

*Proof.* Pick `u_L in L` with `hat u - u_L = eta v`, `|v|_inf = 1`; then `|u_L|_inf <= 1 + eta`. Since `|<k',y>| <= d |y|_inf`
on `B_1`, `|Psi_{c,i}| <= |c|_1 d^{m'} (1+eta)^{m'-i}`. By Proposition 5.1,
`kappa_c <= sum_{i >= o_F} binom(m',i) (1+eta)^{m'-i} eta^i <= eta^{o_F} (2+eta)^{m'}`. The first bound is J2 (§7.2). QED

So any radius in ZE-TC_eps is forced to be at least of the size the judge predicted: a point at distance `eta` from
`L` and farther from all other strata has `kappa <= 3^{m'} eta^{o_F}`, so ZE-TC_eps can hold with radius `rho_F` only if
`rho_F >= (eps 3^{-m'})^{1/o_F}` at such points. (R11 §8.2 wrote `(eps 2^{m})^{1/(m - d_F)}` "up to polynomial
factors"; the rigorous version has `3^{-m'}` in place of `2^{m}`, a factor `e^{O(m'/o_F)}`, and the sharper exponent
`o_F = m' - rho_A(F) >= m' - d_F`.)

**Proposition 5.3 (local lower bound from the transversal leading forms).** Extend `kappa` homogeneously:
`kappa(y) = sup_c |Phi_c(y)|/(|c|_1 d^{m'})` for all `y in R^d`. Let `M` be a linear complement of `L` in `R^d` and
`u_L in L` with `|u_L|_inf = 1`. Suppose the forms `v -> Psi_{c,o_F}(u_L, v)`, `c in V_H`, have no common zero on
`M \ 0`. Then there are `gamma, eta_0 > 0` and a neighbourhood `N` of `u_L` in `L` such that
`kappa(u' + eta v) >= gamma eta^{o_F}` for all `u' in N`, `v in M` with `|v|_inf = 1`, and `0 <= eta <= eta_0`.

*Proof.* Fix a basis `c_1, ..., c_q` of `V_H`. `G(u', v) = max_l |Psi_{c_l,o_F}(u', v)|/(|c_l|_1 d^{m'})` is continuous
and positive on the compact set `{u_L} x S`, `S = {v in M : |v|_inf = 1}`, hence `>= 2 gamma_1 > 0` on `N x S` for a
closed neighbourhood `N` of `u_L` in `L` with `|u'|_inf <= 2` on `N`. By Proposition 5.1 (valid at every `u' in L`)
and `|Psi_{c,i}(u',v)| <= |c|_1 d^{m'} 2^{m'-i}`,
`kappa(u' + eta v) >= binom(m',o_F) eta^{o_F} 2 gamma_1 - sum_{i > o_F} binom(m',i) eta^i 2^{m'-i}
>= eta^{o_F} (2 binom(m',o_F) gamma_1 - eta 3^{m'})` for `eta <= 1`. Take `eta_0 = binom(m',o_F) gamma_1 3^{-m'}` and
`gamma = binom(m',o_F) gamma_1`. QED

If the forms `Psi_{c,o_F}(u_L, .)` do have a common zero in `M`, the same argument runs at the first order
`o > o_F` at which they do not, but only pointwise in `u_L` (the intermediate `Psi_{c,i}`, `o_F <= i < o`, need not
vanish at nearby `u'`). The constants `gamma, eta_0` are **not** effective in this generality. Making them effective
(any bound `kappa >= e^{-o(n log n)}` at the radii of §6.4 suffices, by §4) is exactly the open quantitative content
of ZE-TC_eps.

**Definition 5.4 (the quantitative hypothesis, made precise).** For a finite set `𝓛` of rational subspaces of
`R^d`, radii `sigma = (sigma_L)_{L in 𝓛}` and `eps > 0`:

> **ZE-TC_eps(A, m'; 𝓛, sigma, eps).** Every `hat u` with `|hat u|_inf = 1` and `dist(hat u, L) >= sigma_L` for all
> `L in 𝓛` has `kappa(hat u) > eps`.

R11's ZE-TC_eps is the case `𝓛` = spans of cheap families. Step 0 (Theorem 4.2) supplies `eps = eps_0`; §6.4 supplies
the radii `sigma_F` for cheap spans unconditionally (Theorem 6.7); §6.5 shows that `𝓛` = cheap spans is too small in
general and must be enlarged to all rational subspaces contained in `TC(A, m')`.

**Proposition 5.5 (scale consistency at m' = ceil(4w)).** For `L = L_F` with `F` cheap, Proposition 5.2 shows that
ZE-TC_eps (radius `sigma_F`, level `eps_0`) fails at every point at distance exactly `sigma_F` from `L_F` and at
distance `>= sigma_G` from the other spans, unless `3^{m'} sigma_F^{o_F} > eps_0`. Since `sigma_F <= 1/2` and, for a family
with parameters `(s, e)` (`s = sum ell_i`, `e = |E|`), `o_F <= m' - 2js - 2e` and `sigma_F >= (1/2) e^{-Q_{s,e} abar/2}`, With the radii of Theorem 6.7
this holds for every cheap `F` simultaneously as soon as

```text
G := min_{s,e} [ m' ln 3 - (m' - 2js - 2e)(ln 2 + Q_{s,e} abar/2) - ln(C d/(m'+1)) + W/4 ]  >  0,
Q_{s,e} = 2^s n(n-1)/((n-e)(n-e-2^s)),   abar = W/(2(n-1)) + 2 ln C',
```

the minimum over `0 <= s <= b = floor((m'-1)/(2j))` and `0 <= e <= (m'-1)/2 - js`. *Computation* (`consistency.py`,
double precision; `r = n`, `C = 1`): at `W = 3n ln n - 14n`, `G > 0` for every `n >= 3662`; the coefficient of `W` in `G`
(`1/4 - max (m' - 2js - 2e) Q/(4(n-1))`) is positive for every `n >= 5185`, so for `n >= 5185` the condition holds for
all `W >= 3n ln n - 14n`. Asymptotically the binding families have `s = 2, 3` and the condition is `n > (8 + o(1)) w`.
So on scale grounds alone the `T = 1`, `m' = ceil(4w)` reduction is **not** vacuous: the rigorous transversal ceiling
never forces ZE-TC_eps to fail at the radii that the arithmetic supplies.

## 6. Cheap families at T = 1

Standing assumption in this section: `2j <= n` (true whenever `m_max >= 4w`, since then `4w <= m_max <= 2d`).

### 6.1 A universal lower bound for widths

**Lemma 6.1 (small balanced vectors).** `S_w := {k in {-1,0,1}^n : sum k = 0, |k|_1 <= 2j}` is contained in `A`, for
**every** profile (any `beta_p in [-1,1]^n`, any `r`). Consequently, for every `v in R^n`,

```text
h_A(v)  >=  2 tau_j(v),        tau_j(v) := sum_{i=1}^{j} (v_(n+1-i) - v_(i))                              (6.1)
```

(`v_(1) <= ... <= v_(n)` the sorted entries; `tau_j` is the sum of the `j` largest minus the `j` smallest entries).

*Proof.* `k in S_w` has `|k_x| <= 1` for all `x`, so `k in B_1`, and `|<k,beta_p>| <= |k|_1 <= 2j <= w`. With `P`
(resp. `N`) the positions of the `j` largest (resp. smallest) entries of `v`, disjoint as `2j <= n`, the vectors
`+-(1_P - 1_N)` lie in `S_w` and give `<k,v> = +-tau_j(v)`. QED

**Lemma 6.2 (window decomposition).** For integer `v` put `ell(v) = v_(n+1-j) - v_(j)` (window length), window
`Win(v) = [v_(j), v_(n+1-j)]`, and `e(v) = #{x : v_x notin Win(v)}`. Then

```text
tau_j(v)  >=  j ell(v) + e(v).                                                                             (6.2)
```

If `ell(v) = 0`, with `mu = v_(j)` the modal value: `tau_j(v) = |v - mu 1|_1`, and if moreover `v_n = mu` then
`h_A(v) = 2 |v - mu 1|_1` **exactly**.

*Proof.* For `i <= j`, `v_(n+1-i) - v_(i) = ell + (v_(n+1-i) - v_(n+1-j)) + (v_(j) - v_(i))`, both brackets `>= 0`. An entry
above the window sits at a sorted position `> n+1-j`, i.e. among the top `j - 1`, and contributes `>= 1` to the first
bracket sum (integrality); symmetrically below. If `ell = 0` every non-modal entry is among the top or bottom `j - 1`
positions, and the bracket sums are exactly `sum_x |v_x - mu|`. If also `v_n = mu`, `<k,v> = sum_{x<n} k_x (v_x - mu)` with
`|k_x| <= 1`, so `h_A(v) <= h_{B_1}(v) <= 2|v - mu 1|_1`; with (6.1) this is equality. QED

### 6.2 Structure of all cheap families

**Theorem 6.3 (structure of cheap families; any profile, any order).** Let `F = (v_1, ..., v_s)` be integer
directions with cost `d_F = sum_i h_A(v_i) <= m' - 1`. Write `ell_i, e_i, Win_i` as in Lemma 6.2 (each `v_i` represented
in `Z^n`; these quantities do not depend on the representative modulo `Z 1`), `b = floor((m'-1)/(2j))`,
`E = {x : v_{i,x} notin Win_i for some i}` (exceptional points), core `= [n] \ E`, and the pattern map
`pi(x) = (v_{i,x})_i`. Then:

* (i) `sum_i (j ell_i + e_i) <= d_F/2 <= (m'-1)/2`. Hence `sum_i ell_i <= b`; at most `b` directions are **long**
  (`j ell_i + e_i >= j`, which includes every direction with `ell_i >= 1`); and `|E| <= sum_i e_i <= (m'-1)/2 - j sum_i ell_i`.
* (ii) On the core, `pi` takes values in `prod_i Win_i`, so the core splits into at most
  `q_F := prod_i (ell_i + 1) <= 2^{sum ell_i} <= 2^b` classes.
* (iii) `L_F` is contained in `Pi(F) := {delta : delta_x = delta_y whenever x, y are core points with pi(x) = pi(y)}`, and
  `dim L_F <= q_F - 1 + |E|`.
* (iv) Conversely, the following families are cheap: any family of sparse directions (`ell_i = 0`, `v_{i,n} = mu_i`)
  and sign columns with `sum_i 2|v_i - mu_i 1|_1 + 2j * #(columns) < m'`; more generally any family with
  `sum_{sparse} 2|v_i - mu_i 1|_1 + sum_{columns} e_p w < m'`. For a sign column whose two sides both have `>= j`
  points, `h_A(a_p) = 2j` exactly; for a sparse direction with base point at the mode, `h_A = 2|v - mu 1|_1` exactly.

*Proof.* (i) Lemmas 6.1-6.2 summed over `i`; `j ell_i + e_i >= j` when `ell_i >= 1`. (ii) `ell + 1 <= 2^ell`. (iii) Each
`delta in L_F` is `sum_i t_i v_i` modulo `R 1`. Functions constant on the `q_F` core classes and arbitrary on `E`, modulo
constants, form a space of dimension `q_F + |E| - 1`. (iv) Theorem 2.1(c) and Lemma 6.2; for the column,
`tau_j(1_{S_p}) = min(j, |S_p|, n - |S_p|) = j`. QED

So, up to the slack between `2 tau_j` and the true width (which is zero for sparse directions with modal base point
and for columns with both sides `>= j`), **every** cheap family is "at most `b` long directions with total window
length `<= b`, plus exceptional points of total weight `<= (m'-1)/2`", and its span lives in a subspace in which all
but `(m'-1)/2` points take at most `2^b` values. Directions taking many values are covered: only `ell + 1` of their
values can occur on the core, and each value outside the window costs at least one unit of `e`.

**Corollary 6.4 (the judge's sketch at m' = ceil(4w), R11 §5.3: proved).** Let `m' = ceil(4w)` and `w >= 10`. Then
`b <= 4`, and every cheap family has

* at most **four** long directions, with total window length `sum ell_i <= 4`;
* at most `2w` exceptional points (`|E| <= (m'-1)/2 < 2w`);
* at most `16` core classes;
* hence, if `n > 2w + 16`, its span lies in the braid arrangement (some pair has `delta_x = delta_y` on `L_F`). Under
  `n - 1 >= 316 ln 4r` one has `69 < sqrt(2528) ln 4 <= w <= sqrt(8/316) d < 0.16 d`, so `w >= 10` and `n > 2w + 16`
  hold throughout that regime.

Both "four" and "16" are attained: four sign columns whose sides all have `>= j` points cost `8j <= 4w <= m'`, which is
`< m'` unless `w` is an even integer, and four Sylvester columns `1[x_i = 0]` (`i <= 4`) realise all 16 classes.

*Proof.* `10 j >= 10(w/2 - 1) >= 4w > m' - 1` for `w >= 10`, so `b <= 4`. The rest is Theorem 6.3; pigeonhole on the
`>= n - 2w > 16` core points. QED

*Exact check* (`structure_check.py`, exact integer enumeration of `A` for 6 random small profiles, `n = 12..14`,
`w in [4.5, 7]`): (6.1) and (6.2) on 2,400 random integer directions, and conclusions (i)-(ii) on 990 random cheap
families; no failures. Small instances of identities proved above.

### 6.3 Braid arrangement: the exact dichotomy in the order

**Proposition 6.5 (braid dichotomy).**
* (a) For every profile and every order `m'` with `2^b < n - floor((m'-1)/2)` (`b = floor((m'-1)/(2j))`), every cheap
  span lies in the braid arrangement. In particular this holds whenever `m' <= 2j floor(log_2(n/2)) + 1`
  (as `m' <= m_max < 0.65 n` for `n >= 2900`, by Theorem 2.1(b)).
* (b) Let `n = 2^k` and take the Sylvester-Hadamard profile (points `F_2^k`, columns `1[a.x = 0]`, `a != 0`, `r = n - 1`,
  any base point). The `k` columns `a = e_1, ..., e_k` have `h_A <= 2j` and their span is collision-free (generic
  points of it have pairwise distinct `delta_x`). So for every admissible order `m' > 2jk` the braid statement fails.
* (c) *Exact:* such orders exist iff `k >= 16` (`2jk < m_max` first at `n = 2^16`: `2jk = 40896 < m_max = 42065`;
  at `2^15`: `26340 > 21034`), and at every `k >= 12` the window is nonempty (`m_max >= 4w`). At `n = 2^k`, (a) and (b)
  are complementary: braid holds for every cheap span iff `m' <= 2jk`.

*Proof.* (a) Pigeonhole on Theorem 6.3(ii) with `|E| <= floor((m'-1)/2)`; for the second sentence, `b <= log_2(n/2)`
gives `2^b <= n/2 < n - m'/2`. (b) Theorem 2.1(c), and `x -> (1[x_i = 0])_i` is injective on `F_2^k`. (c)
`sylvester_out.txt` (log-gamma `m_max` at 60 digits, confirmed by exact big-integer comparison for `n <= 2^16`). If
`m' <= 2jk` then `b <= k - 1`, `2^b <= n/2 < n - m'/2`, and (a) applies; if `m' > 2jk`, (b) applies. QED

So the referee's ZE-braid refutation (`referee_compute.md` §2.8) has an exact `T = 1` analogue: at Sylvester
orders the `log_2 n`-column subtori become cheap exactly when `m' > 2j log_2 n`, which is inside the admissible window
from `n = 2^16` on, and **never** at `m' = ceil(4w)` (where `b <= 4`). For random sign profiles about `2 log_2 n`
columns separate all points, so the same happens for `m' >~ 4j log_2 n`.

### 6.4 The arithmetic separation lemma (proved) and the reduction

**Lemma 6.6 (collision separation; any cluster, any norm).** For a linear subspace `L` of `delta`-space modulo `R 1`
let `P(L) = #{x < y : delta_x = delta_y for all delta in L}` (the number of *collision pairs* of `L`). If `P(L) >= 1`,
the actual position satisfies

```text
dist(hat u, L)  >=  sigma_L := (1/2) exp( - binom(n,2) abar / (2 P(L)) ),     abar = W/(2(n-1)) + 2 ln C'.     (6.3)
```

*Proof.* Put `eps_xy = d_xy - W/2 + 2 ln C'`; by Lemma 1.1, `eps_xy >= 0`, and by Lemma 1.2,
`sum_{x<y} eps_xy <= binom(n,2) abar`. Averaging over the `P(L)` collision pairs, some collision pair `(x,y)` has
`eps_xy <= binom(n,2) abar/P(L)`. By Lemma 1.1, `|hat u_x - hat u_y| = |delta_x - delta_y|/|u|_inf >= e^{-d_xy/2}/(C e^{-W/4})
= (C'/C) e^{-eps_xy/2} >= e^{-eps_xy/2}`. Every `y_L in L` (in `u`-coordinates, `u_n = 0`) has `(y_L)_x = (y_L)_y`, so
`|hat u - y_L|_inf >= |hat u_x - hat u_y|/2`. QED

**Theorem 6.7 (separation from all cheap spans).** Let `F` be cheap for `(A, m')`, with the notation of Theorem 6.3,
`n' = n - |E| > q_F`. Then `P(L_F) >= n'(n' - q_F)/(2 q_F)`, so

```text
dist(hat u, L_F)  >=  sigma_F := (1/2) exp(-Q_F abar/2),      Q_F = q_F n(n-1) / (n'(n' - q_F)).
```

At `m' = ceil(4w)` with `w >= 10` and `n > 2w + 16`: `Q_F <= 16 n(n-1)/((n-2w)(n-2w-16)) = 16(1 + O(w/n))`, so
`sigma_F >= (1/2) C'^{-16(1+O(w/n))} exp(-(4 + O(w/n)) W/(n-1))`. For `W = c n ln n` this is `n^{-4c(1+o(1))}`:
**a polynomial separation**, for every cheap span at once, with no nondegeneracy, rank or squarefree hypothesis.

*Proof.* The core classes have sizes `n_1, ..., n_q` (`q <= q_F`) summing to `n'`; by Theorem 6.3(iii) every pair inside a
class is a collision pair of `L_F`; by convexity `sum binom(n_c, 2) >= q binom(n'/q, 2) >= n'(n'-q_F)/(2q_F)`. Apply
Lemma 6.6. QED

This is the "arithmetic separation lemma" of R11 §8.2, now **proved** at `m' = ceil(4w)` for cheap spans: R11 needed a
pair orthogonal to `L_F` with small excess; Theorem 6.3 supplies `>= n^2/32` such pairs and Lemma 1.2 bounds their
average excess. The radius needed by ZE-TC_eps near `L_F` is `~ (eps_0 3^{-m'})^{1/o_F} = e^{-(W/4 - O(m'))/o_F}`, while
`sigma_F` is only polynomially small; Proposition 5.5 records when the two are compatible (`n >= 5185` at `r = n`).

**Theorem 6.8 (conditional reduction at T = 1, m' = ceil(4w); all norms).** Let a cluster have `n >= 5`,
`m_max >= 4w`, `m' = ceil(4w)`, `W > W_1(n, C, m')`. Let `𝓛` be any set of rational subspaces with `P(L) >= 1`, and
`sigma_L` as in (6.3). If ZE-TC_eps`(A, m'; 𝓛, sigma, eps_0)` (Definition 5.4) holds for the cluster's own `T = 1` slab
box, the cluster does not exist. Consequently: if, for some `n_1`, every profile on `n >= n_1` points with
`n - 1 >= 316 ln 4r` satisfies

* **(Q1)** every rational subspace `L` of `R^d` contained in `TC_H(A, m')` has `P(L) >= binom(n,2)/Q_0(n)`, and
* **(Q2)** ZE-TC_eps`(A, m'; 𝓛_A, sigma, eps_0)` with `𝓛_A` = all those subspaces, `sigma_L = (1/2) e^{-Q_0 abar/2}`, and
  `eps_0 = C d e^{-W/4}/(m'+1)`, for every `W > W_1`,

then every cluster with `W > W_1(n, C, ceil(4w))` has `n < max(n_1, 316 ln(4r) + 1)`, so `n = O(log W) = O(loglog R)` on
**all** norms (`r <= W/ln 5`).

*Proof.* Theorem 4.2 gives `kappa(hat u) <= eps_0`; Lemma 6.6 gives `dist(hat u, L) >= sigma_L` for every `L in 𝓛`; the
hypothesis then gives `kappa(hat u) > eps_0`. QED

Compared with [D] Thm 5.7 and R11 §8.2, this route needs no nondegenerate subcluster, hence neither [D] Lemma 3.1
(squarefree rank) nor Lemma 3.2: distinctness of the points is all that Lemma 6.6 uses. Theorem 6.7 proves (Q1) with
`Q_0 = 16(1 + O(w/n))` **for the subspaces that are spans of cheap families**. §6.5 shows that this class is too small.

### 6.5 Beyond widths: the cheap-span form of ZE-TC is false; the correct qualitative core

R11 §8.2 states the qualitative core as *ZE-TC(A,m'): every real common zero of the `Phi_c`, `c in V_H`, lies in the
real span of a cheap family*. The tangent cone, however, contains rational subspaces that are not cheap spans.

**Lemma 6.9 (ternary slices and single cuts).** For `s >= 1`, `|t| <= s` and `0 <= j <= s`:
* `rho(Lambda_t(s)) <= s - |t|`, where `Lambda_t(s) = {k in {-1,0,1}^s : sum k = t}`;
* `rho(A2(s,j)) <= s + j`, where `A2(s,j) = {k in {-1,0,1}^s : |sum k| <= j}`.

Both are equalities for `s >= 2` (exact for `s <= 6`, all `t`, `j`: `single_cut_out.txt`; for `A2` the lower bound
`s + j - [s-j odd]` also follows from the zonotope with `j` coordinate generators of length 2 and `floor((s-j)/2)` pair
generators `e_x - e_y` of length 2).

*Proof.* (Newton interpolation across level sets, as in compute.md Prop 2.4.) If `Lambda = Lambda_{t_0} ⊔ ... ⊔ Lambda_{t_q}`
are level sets of `sigma = sum k` and every function on `Lambda_{t_i}` is a polynomial of degree `<= e_i`, then every
function on `Lambda` is `sum_i prod_{l<i} (sigma - t_l) q_i` with `deg q_i <= e_i` (choose `q_i` successively; the
product is a nonzero constant on `Lambda_{t_i}`), so `rho(Lambda) <= max_i (i + e_i)`.
*Slices*, by induction on `s` (`s = 1` is trivial). By symmetry `t >= 0`. On `Lambda_t(s)`, `k_s = t - sigma'(k_1..k_{s-1})`,
so `Lambda_t(s)` is the graph over `Gamma = {k'' in {-1,0,1}^{s-1} : sigma' in {t-1, t, t+1}}` and
`rho(Lambda_t(s)) <= rho(Gamma)`. Order the levels of `Gamma` as `(t-1, t, t+1)` if `t >= 1` and `(0, 1, -1)` if `t = 0`
(dropping absent levels). The induction hypothesis gives `max_i(i + s - 1 - |t_i|) = s - t` in both cases.
*Single cut*: order the levels `0, 1, -1, 2, -2, ...`; the `i`-th has `|t_i| = ceil(i/2)`, so
`max_{i <= 2j} (i + s - ceil(i/2)) = s + j`. QED

**Lemma 6.10 (coordinate projections).** If `Y` is a set of coordinates in `[d]` and `rho(pi_Y(A)) < m'`, then `R^Y` is
contained in `Bs(A, m')` and in `TC(A, m')` (and in `TC_H`). In particular, if the profile has a sign column (`e_p = 1`)
whose side `S'` not containing the base point has `s` points, then `rho(pi_{S'}(A)) <= s + j`, so
`R^{S'} ⊂ TC(A, m')` whenever `s + j < m'`.

*Proof.* For `u in R^Y`, `e^{i<k',u>}` and `<k',u>^{m'}` are functions of `k_Y in pi_Y(A)`, hence polynomials of degree
`<= rho(pi_Y(A)) < m'` in `k_Y`. On `A` the slab of the column reads `|sum_{x in S'} k_x| <= j`, so `pi_{S'}(A)` lies in
`A2(s, j)`; `rho` is monotone; Lemma 6.9. QED

**Theorem 6.11 (the cheap-span form of ZE-TC is false).**
* (a) *Exact model.* `n = 7`, one sign column with side `S' = {1,2,3,4}` (base point outside), `j = 1` (any
  `w in [2,4)`), `A = {k' in {-1,0,1}^6 : |k_1+k_2+k_3+k_4| <= 1}` (`|A| = 459`, `rho(A) = 9`). For `m in {6, 7, 8}`,
  `R^{S'}` lies in `TC(A,m)` and in `Bs(A,m)`, but in no span of a family of width cost `< m`.
* (b) *Regime family.* Let `n = 2^12` and take the Sylvester-Hadamard profile with all `log p in [tau, tau(1+2/n)]`, plus
  one sign column with side `S'_0` of size `s`, `1129 <= s <= 1974`, not containing the base point, carried by a prime
  `q_0` with `log q_0 <= (3n-2) tau/(n(n-1))` (e.g. `q_0 = 5`, `tau` large). Then `r = n`, the profile has exponent rank
  `n - 1` and pair separation `>= 1/2 - 1/(2n)` (the admissibility of SL\*, R11 §5.1), `m_max = 2631 >= m' = ceil(4w)
  = 2256`, `j = 281`, and `R^{S'_0}` lies in `TC_H(A, m')` but in no cheap span. So ZE-TC in the cheap-span form fails for
  an admissible regime profile at the recommended order `m' = ceil(4w)`.

*Proof.* (a) Lemma 6.9 gives `rho(A2(4,1)) = 5 < m`, so Lemma 6.10 applies. (Exact cross-check, `single_cut_out.txt`:
for `m = 6..9`, every vector of an exact basis of `V(A,m)` (dimensions 193, 96, 33, 6) has `Phi_c(u) = 0` at
`u = (1,3,7,19,0,0)`, and the generic full-support `u = (1,3,7,19,43,101)` is not a common zero.) A family whose span
contains `R^{S'}` has restrictions to `S'` spanning `R^{S'}`, hence at least 4 members with nonzero restriction; `A` is
the product of `A2(4,1)` and the cube on the two other coordinates, so `h_A(v) >= h_{A2}(v|_{S'})`, and every nonzero
`v|_{S'}` has `h_{A2} >= 2` (`+-(e_x - e_y)` and `+-e_x` lie in `A2(4,1)`). So the cost is `>= 8 >= m`.
(b) `w = sqrt(8 * 4095 * ln 16384) = 563.8`, `j = 281`, `m' = 2256`, `b = 4`, `m_max = 2631` (`sylvester_out.txt`,
`params_T1_out.txt`). Lemma 6.10: `s + j <= 2255 < m'`. If a cheap span contained `R^{S'_0}`, all points of `S'_0` would be
pairwise in different classes and outside the class of the base point, so by Theorem 6.3
`s <= |E| + q_F <= max_{0<=t<=4} (1127 - 281 t + 2^t) = 1128`, a contradiction. Admissibility: the Sylvester columns
give rank `n - 1`; a pair not separated by the new column is separated by `n/2` columns of total weight
`>= (n/2) tau`, against `W_P <= (n-1)(1+2/n) tau + log q_0`, and the ratio is `>= (n-1)/(2n)` iff
`(3n-2) tau >= n(n-1) log q_0`; pairs separated by the new column only gain. QED

**The correct qualitative core.** Replace "span of a cheap family" by "rational subspace contained in the tangent
cone":

> **ZE-TC^rat(A, m').** Every real `u` with `Phi_c(u) = 0` for all `c in V_H(A,m')` lies in a rational linear subspace
> contained in `TC_H(A, m')`.

By Proposition 5.1 (with `rho_A`) and Lemma 6.10 such subspaces include all spans of *projection-cheap* families
(`rho_A(F) < m'`); these lie in the base locus and can be much larger than width-cheap spans (Theorem 6.11). Whether
every rational subspace of the tangent cone is projection-cheap was not decided.

*Known cases.* On the uncut cube the leading-form ideal is `(k_x^3)` and `TC = {|supp u| < m'/2}`, a union of coordinate
subspaces, so ZE-TC^rat holds. For a single cut (base point outside the cut side `S'`) the box is the product of the
cube on `[d] \ S'` and `A2(s, j)`, and

```text
TC(A, m') = {(u_1, u_2) : 2|supp u_1| + nu(u_2) < m'},     nu(u_2) = max{i : <., u_2>^i notin J(A2(s,j))},
```

so ZE-TC^rat for one cut reduces to the single-cut cube. *Proof.* `J(A_1) + J(A_2)` is contained in `J(A_1 x A_2)` and both
quotients have total dimension `|A_1||A_2|`, so they are equal. In `S_1/J_1 (x) S_2/J_2` the components
`l_1^i (x) l_2^{m'-i}` of `(l_1 + l_2)^{m'}` have distinct bidegrees, so the power vanishes iff each component does,
iff `nu_1 + nu_2 < m'`; on the cube `nu_1(u_1) = 2|supp u_1|`. The exact model of Theorem 6.11(a) agrees:
`u = (1,3,7,19,5,0)` has `nu <= 2 + 5` and is a common zero exactly for `m = 8, 9` (`single_cut_out.txt`). I did not
decide ZE-TC^rat for one cut beyond the exact model (open point 3).

**What (Q1) needs beyond Theorem 6.7.** In Theorem 6.11 the extra subspace `R^{S'_0}` has `P >= binom(n - 1 - s, 2)`
collision pairs, so Lemma 6.6 still separates the actual point from it polynomially. In general, (Q1) asks for a
collision lower bound for every rational subspace of `TC_H(A,m')`. Widths give it (Theorem 6.7); projection costs do
not obviously: the only general lower bound I have is `rho_A(F) >= 2 min(j, n - (largest class of F))` (zonotope with
`j` pair generators inside `S_w`), which says nothing at `m' ~ 8j`. Whether a collision-free rational subspace (all classes singletons)
can lie in `TC_H(A, ceil(4w))` is open; if one does, Lemma 6.6 gives nothing and route 5 at `T = 1` needs the full
character inequalities of the relations orthogonal to it.

## 7. Recheck of J1, J2, J3

### 7.1 J1 (SL\* is false at Sylvester orders, `T = n`): confirmed

Every step of R11 §5.1 was re-derived; the numerical inputs were recomputed independently (`sylvester.py`).
* *Column identity.* For each class `c` of `(a.x, b.x)` the coefficients `(alpha, beta, gamma, delta)` are
  `(0,0): (1,-1,-1,2)`, `(0,1): (-1,1,1,2)`, `(1,0): (1,1,-1,-2)`, `(1,1): (1,-1,1,-2)` (the last is the judge's displayed
  example). Checked on all `(a, b, c, x)`, `k = 3..7` (9,326,784 cases), 0 failures; since only the four classes
  enter, it holds for every `k`. Subtracting the value at a base point outside `V_c` gives
  `a'_p = beta a'_b + gamma a'_{a+b} + delta sum_{x in V_c} e_x`, `|V_c| = n/4`.
* *Lightness.* `|Q| tau <= W_P/8 <= (n-1)(1+2/n) tau/8`, so `|Q| <= (n-1)(n+2)/(8n)`; `Q \ {a}` meets at most `|Q| - 1`
  of the `(n-2)/2` disjoint pairs `{b, a+b}`, so a pair with both narrow slabs exists.
* *Cheapness on `A_{Q,lambda}`* (`T = n`): cost `<= 4 floor(w/2) + (n/4)(2n)`. At `n = 2^13`: exact
  `m = 49388999` (log-gamma at 60 digits, confirmed by exact big-integer comparison at `m` and `m + 1`) against cost
  `4.70779e7`: cheap, and this **does** need the exact `m` (`0.7 n(n-1) = 4.697e7 <` cost). At `n = 2^14`: cost
  `1.73724e8 < 0.7 n(n-1) = 1.87893e8` (and `< m = 197534058`). The bound `(n/2 + 2 sqrt(8(n-1) ln 4(n-1)))/(0.7(n-1))`
  for cost`/(0.7 n(n-1))` is `0.9246` at `2^14` and decreasing (both summands decrease), so `k >= 14` follows from
  [D] Thm 2.1(b).
* *Conclusion and admissibility.* theory.md Prop 3.5(ii) (refereed) applies to any `(S, m)`; rank `n - 1`; pair
  separation `n^2/(2(n-1)(n+2)) >= (n-1)/(2n)` (iff `3n >= 2`); `n - 1 >= 300 ln(4(n-1))` from `n = 2^12`. **J1 is
  correct as stated.**
* *At T = 1 (new).* The same family costs `4j + 2(n/4) = 4j + n/2` on the `T = 1` box, so it is cheap at order `m'`
  iff `m' > n/2 + 4j`. This never happens at `m' = ceil(4w)` in the regime (it would need `n/2 < 4w + 1 - 4j <= 2w + 5`, while
  `n - 1 >= 316 ln 4r` gives `w < 0.16 d`), and happens
  inside the window `(.., m_max]` from `n = 2^15` on (`19896 < 21034`; at `2^14`: `10600 > 10518`). So a `T = 1` version
  of SL\* at `m' = ceil(4w)` is not refuted by J1. (Theorem R3's proof transfers to `T = 1` with Lemma 3.1 in place of
  theory.md Lemma 1.2, since `cost <= wW + m' sum_Q lambda_q log q <= 3m'W/8`; whether a `T = 1` SL\* is plausible
  was not examined.)

### 7.2 J2 (the leading-order criterion is vacuous at `T = n`): confirmed

*Proof re-derived.* If `|<k',u>| <= R_0` on `supp c` and `c in V(S, m)`, the push-forward `mu` of `c` under
`k -> <k',u>` has vanishing moments of order `< m`, so `Phi_c(u) = int (x^m - p) d mu` for every `p` of degree `< m`;
the monic Chebyshev polynomial gives `|Phi_c(u)| <= |c|_1 2 (R_0/2)^m`, hence `kappa_c <= 2^{1-m}` (`R_0 <= D/2` for
`|u|_inf = 1`). With `|u|_inf >= e^{-W/2}` (Lemma 1.1, `d_xn <= W`), theory.md (2.1) needs `W > W_TC = 2(m-1) ln 2 +
2 ln((D/2)/(m+1))`. *Numbers* (`j2_check_out.txt`; exact `m` with big-integer boundary check): `n = 2137`:
`m = 3362393`, `W_TC = 4.661e6 = 1.0207 n^2`, against `3n ln n - 14n = 19236` (factor `242.3`); `n = 3000`: factor `305.6`.
**J2 is correct.** (Its last bullet is also right: only the sufficient criterion with the uniform remainder fails.)

### 7.3 J3 (route 5 at `T = 1`): confirmed, with closed forms and two sharpenings

* (a) Confirmed (Theorem 2.1(a)).
* (b) Confirmed and sharpened: `m_max >= c_* d - kappa` with `c_* = 0.6417 < c_0 = 0.641794...`, `kappa < 2.2137`, an
  `O(1)` loss instead of `O(log n)` (Theorem 2.1(b)); exact `m/n in [0.64212, 0.6440]` on `[1000, 8000]` (J3:
  `[0.6421, 0.6440]`).
* (c) Confirmed: `dim V_H >= U/2` iff `ln U >= 3 ln 3 = 3.296`, true for `n >= 5` (Theorem 2.1(d)).
* (d) Confirmed exactly: thresholds `2900` (`r = n`), `3623` (`r = n ln n`), `5809` (`r = n^2`, with an isolated success
  at `5807`), now on `[2, 12000]`; closed-form sufficient conditions `n - 1 >= 316 ln 4r`, `n >= 2966` (`r <= n`),
  `n >= 5930` (`r <= n^2`) (Theorem 2.1(e)).
* *Theorem 4.1 at T = 1*: confirmed with the closed form `W_1 = 24 d ln 3/m' + 8 ln(e d C'/m')` for **every** order in
  the window (Theorem 3.1); exact-form values `51.75` (`n = 3000`, `m_max`), `53.27` and `630.68` (`m' = ceil(4w)`,
  `n = 3000`, `10^6`); `W_1(m_max)` is decreasing in `n` on the table; the closed form gives `W_1(m_max) <= 52.71 + 8 ln C'` for
  `n >= 2900`, tending to `52.63 + 8 ln C'`.
* *Scale separation*: `W_TC = 1.386 m` and log-room `6178, 9101, 22363` at `n = 3000, 4000, 8000`: confirmed (§4).
* *Near the strata*: the transversal order is proved (Proposition 5.1, with the sharper `m' - rho_A(F)`); the
  quantitative form is Proposition 5.2; the "arithmetic supply" heuristic is now a theorem (Theorem 6.7); the "at
  most `2^4` column classes" sketch is Corollary 6.4, for **all** directions including those taking many values.
* *Correction to the programme of R11 §8.2.* The qualitative core "common zeros lie in spans of cheap families" is
  false in general (Theorem 6.11, including an admissible regime profile at `m' = ceil(4w)`); it must be stated for
  rational subspaces of the tangent cone (ZE-TC^rat), and the separation lemma must then cover those (open (Q1)).

## 8. Open points

1. **(Q1) Collisions of rational subspaces of the tangent cone.** At `m' = ceil(4w)`: does every rational subspace
   `L` contained in `TC_H(A, m')` have `>= c n^2` collision pairs (as cheap spans do, Theorem 6.7, and as the
   projection-cheap spans of Theorem 6.11 do)? Is there any collision-free rational subspace in `TC_H(A, ceil(4w))`?
   The width method cannot decide this; a lower bound for `rho_A(F)` in terms of the classes of `F` is what is needed.
2. **ZE-TC^rat and (Q2).** Decide ZE-TC^rat`(A, ceil(4w))` for regime profiles, and make Proposition 5.3 effective:
   `kappa(hat u) >= e^{-o(n log n)}` at polynomial distance `sigma` from all rational subspaces of the tangent cone. With
   (Q1) this gives `n = O(log W)` on all norms (Theorem 6.8). Proposition 5.5 shows the scales are compatible for
   `n >= 5185` (`r = n`).
3. **One and two cuts.** For one cut the problem reduces to the single-cut cube `A2(s,j)` (§6.5). Compute
   `nu(u) = max{i : <.,u>^i notin J(A2)}` on `A2(s,j)`. The bounds proved here give
   `nu(u) <= min(2|supp u|, s + j, min_t (2|{x : u_x != t}| + 2j))`; if equality holds, ZE-TC^rat holds for one cut.
4. **SL\* at T = 1.** J1's family is not cheap at `m' = ceil(4w)` (§7.1); a `T = 1` grouped-widening hypothesis is
   untested.
5. **Constants.** The `316 ln 4r` and `n >= 5185` thresholds come from Hoeffding with a union bound and from
   Markov-type averaging; neither is optimised.
