# Round 21: the low-height window `P_+(N) < M^(2-delta)` and the large sieve

Scratch write-up (round 21). Scripts and outputs: `prove_checks/` next to this file. Repository files are cited by
their path under `docs/claude_review/`; they were read only.

Labels, as in the overview: **proved** (all quantifiers, here); **quoted** (a published or record result, with its
source); **exact** (integer or fixed-point computation with an error bound); **evidence** (a computation no proof
uses); **heuristic**. Nothing here has been refereed.

## 0. Summary and verdict

**No growth theorem in the window `P_+(N) < M^(2-delta)`.** `M <= (2/3 + eps) log R/loglog R` is still the best bound
there. The large-sieve route fails, and this write-up proves where and why. It also proves a cheap restricted growth
theorem just below the window, and it reduces the window question to a fake-existence question in which the large
sieve plays no part.

| # | statement | status |
| --- | --- | --- |
| Thm 2.1, Cor 2.2 | `M < 4 pi(P_+(N'); 4, 1)(1 + 2M log^+ C/W')`; hence **`M = o(log R/loglog R)` for clusters with `P_+(N') = o(log R)`**, uniformly in position. Near-extremal clusters have `P_+(N') > (eps/12.56) log R` (Cor 2.3), so in units of `M` the window is `[c M log M, M^(2-delta)]` | proved (from `stability.md` Thm 3.3 and Rosser-Schoenfeld); restricted growth, **below** the window |
| Remark 2.4 | counting cannot pass exponent 1: near-extremal clusters use `Theta(M)` primes, and `pi(M^a) >> M` for `a > 1` | proved |
| Lemma 3.2, Prop 3.3 | the large sieve for any set of angles is equivalent, up to a factor `< 23`, to its arc counts at scale `1/K`; for a cluster the arc counts are the lattice count `eta Y/2 + 2.53 Y^(1/2) + 1.58`, non-trivial in the window by a factor `M^(delta/2)` | proved, explicit constants |
| Lemma 4.1, Prop 4.2, Cor 4.3 | the realisation set is invariant under the group `H = ker(B mod pi/2)`; whenever the coordinates are pairwise independent in `Cok = Z^s/B^T Z^(M-1)`, it contains points that are as spread as genuine prime angles (large-sieve constant within a factor 3). **The short-arc condition never forces concentration**, so the large sieve cannot contradict it. Pinned profiles (`rank B = s`) are handled by the monomial gaps, which are stronger there (Remark 4.4) | proved (Remark 4.4 is a remark with a profile-dependent threshold) |
| Prop 5.2, 5.3, Cor 5.4 | `round10/analytic.md` Theorem B extends into the window: spread grafts satisfy every robust large sieve up to `1 + O(M^(-eta/2) (log M)^(1/2))` even when `s >> X^(1/2)`, and robust prime sums and box counts for `gamma > 1/(1+eta)`. So **if** spread gap-admissible window fakes exist, the window is capped for `A_0` + robust analytic inputs | proved; the cap is **conditional** on window fakes |
| Prop 6.1, 6.2 | the real obstruction to window fakes is the monomial gaps: the first-moment (union-bound) construction needs `sum_j p_j^(-1/2) = O(1)`, which fails in the window by `M^(delta/2)` (the same threshold `s ~ P^(1/2)` as the large sieve); lattice angles are gap-admissible but collapse to actual clusters | proved |
| Prop 8.1 | for `A_0`, growth in the window `<=>` (modulo the spreading proviso) `W_min` of gap-admissible angle models with primes `<= M^(2-delta)` is `>> M log M` | proved (a reduction) |
| §7.1 | exhaustive scan of all 200-smooth split-only `n <= 10^18` (`4.07 * 10^7` norms): the window is populated by small clusters, the record constants for `k = 7, 9, 10` reappear, and for every `k >= 7` the best `C_k` per decade **increases** with `n` beyond `10^9`; the extremal rate cannot be tested at `M <= 11` | evidence |
| §7.2 | every record cluster has `sum p_j^(-1/2)` between `0.99` and `1.43`, i.e. sits where first-moment fake constructions fail; the 7- and 10-point records are pinned (`rank B = s`, `D_j = 24, 8, 6, 4, 24`) | exact |
| §7.3 | the Haar measure of the gap-admissible set behaves like `exp(-1.2 E)` (Poisson with mild correlation) up to `sum p^(-1/2) ~ 1.2`: positive, no collapse to lattice points, so window fakes plausibly exist (Heuristic 6.3) | evidence (direct and sequential Monte Carlo) |

**Where it fails (§8.1).** Not because the cluster uses too few of the available angles (it uses enough for the large
sieve to be non-trivial), but because (F1) the large sieve sees only how many angles lie in an arc, (F2) the
short-arc condition constrains the angles only through `B phi mod pi/2`, which is invariant under a group that spreads
every free coordinate, while the cut/chain structure lives entirely in the exponent profile, and (F3) once spread,
grafts pass every robust input even with `s >> X^(1/2)`. What makes the window special is (F4): gap-admissibility for
many small primes, which first-moment methods cannot establish.

**Next step, if pursued (§8.4):** a rigorous sequential construction of gap-admissible spread angle models with
primes in `[M^(3/2), 2M^(3/2)]` on the Paley double-core profile. Success would cap the window with Corollary 5.4; a
reasoned failure would be the first evidence that the window is genuinely different.

## 1. Setting, notation, imported results

**Clusters.** `N = R^2`, `W = log N`, and `z_1, ..., z_M` are distinct Gaussian integers of norm `N` on an arc of
length `<= C sqrt R` (so `|z_x - z_y| <= C sqrt R`). After division by `G = gcd(z_1, ..., z_M)`
(`two_thirds/proof.md` Lemma 1.1) we get `u_x = z_x/G` of norm `N' = prod_(j=1)^(omega) p_j^(e_j)`, all `p_j = 1 mod 4`,
`W' = log N' <= W`, and

```text
u_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj),   arg u_x = <2a_x - e, theta> + (pi/2) k_x,   theta_j = arg pi_j .
```

`s := omega(N')` is the number of cluster primes. `K' = K_0 + 4 log^+ C`, `K_0 < 13.9092`. `P_+(n)` is the largest
prime factor; `pi(x; 4, 1)` counts primes `p <= x`, `p = 1 mod 4`.

**The window.** For `delta in (0, 1)`, the window is the class of clusters with `P_+(N') <= M^(2 - delta)`. (Dividing by
`G` does not change the angles and does not increase the chord constant (`round8/stability.md` §1), so every
statement below is about `N'`; since `P_+(N') <= P_+(N)`, a hypothesis on `P_+(N)` implies the one on `P_+(N')`.)

**Imported results (all proved and refereed in the record).**

* (I1) `two_thirds/proof.md` Thm 7.5: `3 M log M <= W' + K' M`.
* (I2) `round8/stability.md` Thm 3.3(2) with (3.0): for `M >= 4`, every prime power `p^e || N'` satisfies
  `e log p <= 2 beta'`, where `beta' = max_(ceil(M/2) <= s' <= M) b'_C(s')`,
  `b'_C(s') = W'/s' + (4(s'-1)/s') log C - 8 F*(s')/s'^2`, and
  `beta' < min(2W'/M + 4 log^+ C, 2W'/M - 3 log(M/2) + K')`.
* (I3) `round8/stability.md` Cor 3.6: if `M >= eps log R/loglog R` and `R >= R_0(eps, C)`, then
  `P_+(N) <= M^(4/eps - 3 + o(1))`, and `M/4 (1+o(1))^(-1) < omega(N') <= Omega(N') <= (2/eps + o(1)) M`.
* (I4) `round10/analytic.md` Lemma 3.2 and Theorem B (graft barrier) and Lemma 5.2 (lattice cross-gap:
  `|sin(arg pi - arg alpha)| >= (Norm pi Norm alpha)^(-1/2)` for non-parallel `pi, alpha in Z[i]`).
* (I5) `overview.md` Lemma 5.1 (monomial gaps): for every `v in Z^s \ {0}`, with `w(v) = sum |v_j| log p_j`,
  `dist(<v, theta>, (pi/2)Z) >= arcsin(e^(-w(v)/2))` and `dist(<v, theta>, pi/4 + (pi/2)Z) >= arcsin(e^(-w(v)/2)/sqrt 2)`.
* (Q1) Rosser-Schoenfeld (as cited in `round8/stability.md` Prop 3.5): `pi(x) < 1.25506 x/log x` for `x > 1`.
* (Q2) Montgomery-Vaughan's large sieve (quoted; H. L. Montgomery and R. C. Vaughan, "The large sieve",
  Mathematika 20 (1973)): for real `x_1, ..., x_n` whose classes in `R/Z` are pairwise at distance `>= delta`, any
  integer `K_1` and `L >= 1`, and complex `b_n`,
  `sum_(k = K_1 + 1)^(K_1 + L) |sum_n b_n e(k x_n)|^2 <= (L - 1 + delta^(-1)) sum_n |b_n|^2`.

**The two angle classes.** For a set `S` of split primes and real `phi in R^S`, the *realisation set* of a profile
`(a_x)` is

```text
Real(Delta) = { phi : there are theta_0, distinct delta_x in [0, Delta], integers k_x with
                      <2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x for all x } .
```

Actual clusters have `theta in Real(C e^(-W/4))`. An *angle model* (fake) is a profile with a `phi in Real(Delta)`
satisfying (I5) for all `v` (`overview.md` Def 5.2).

## 2. The counting boundary of the window (proved)

The window has a lower edge, and below it growth is easy. This section proves that and shows why counting stops
there.

**Theorem 2.1 (prime count bounds the cluster).** Let `M >= 5` distinct lattice points lie on an arc of length
`<= C sqrt R` (any position, any `C > 0`). Then `N' > 1`, and with `s = omega(N')`, `Y = P_+(N')`:

1. `W' <= 2 beta' s < (4W'/M + 8 log^+ C) s`, and `s <= pi(Y; 4, 1)`.
2. `M < 4 pi(Y; 4, 1) (1 + 2 M log^+ C / W')`. In particular `M < 4 pi(P_+(N'); 4, 1)` when `C <= 1`.
3. If `3 log M > K'`, then `2 M log^+ C/W' <= 2 log^+ C/(3 log M - K')`.
4. (Sharper form.) Put `A = W'/(M log M)`. If `D := 4A log M - 6 log(M/2) + 2K' > 0`, then
   `M < s D/(A log M)`. As `M -> infinity` with `A` bounded, this reads `M <= ((4A - 6)/A + o(1)) s`; at the
   2/3-extremal rate `A -> 3` it is `M <= (2 + o(1)) pi(Y; 4, 1)`.

*Proof.* Five distinct points of one norm force `N' > 1`, since norm 1 carries only the four units. By (I2) each
prime power `p^e || N'` has `e log p <= 2 beta'`; summing over the `s` primes of `N'` gives
`W' = sum e_p log p <= 2 beta' s`, and (I2) gives `2 beta' < 4W'/M + 8 log^+ C` and
`2 beta' < 4W'/M - 6 log(M/2) + 2K' = D` (the last equality is the definition of `D` with `W' = A M log M`).
All primes of `N'` are `= 1 mod 4` and `<= Y`, so `s <= pi(Y; 4, 1)`. This is item 1. Dividing
`W' < s (4W'/M + 8 log^+ C)` by `W' > 0` gives `1 < s (4/M + 8 log^+ C/W')`, i.e. item 2. Item 3: by (I1),
`W' >= 3M log M - K' M > 0`, so `M/W' <= 1/(3 log M - K')`. Item 4: `W' <= 2 beta' s < D s` and `W' = A M log M`.
QED

**Corollary 2.2 (restricted growth below the window).** Fix `C > 0` and let `M_1(C)` be the least integer
`>= 5` with `3 log M_1 >= K' + 8 log^+ C`. Every cluster with `M >= M_1(C)` satisfies

```text
M < 5 pi(Y; 4, 1) < 6.28 Y/log Y,        Y = P_+(N'),
```

and hence `M < max( M_1(C), 12.56 max((log R)^(1/2), Y)/loglog R )` once `(log R)^(1/2) >= e`. Consequently,
**for every function `Y(R) = o(log R)`, the clusters with `P_+(N') <= Y(R)` (in particular those with
`P_+(N) <= Y(R)`) satisfy `M = o(log R/loglog R)`**, uniformly in the position of the arc.

*Proof.* For `M >= M_1(C)`, item 3 gives `2M log^+ C/W' <= 2 log^+ C/(8 log^+ C) = 1/4` (and `0` if `C <= 1`), so
item 2 gives `M < 5 pi(Y; 4, 1) <= 5 pi(Y) < 6.28 Y/log Y` by (Q1) (`5 * 1.25506 < 6.28`). The function
`y/log y` increases for `y >= e`, and `Y >= 5`. If `Y <= (log R)^(1/2)`, then
`Y/log Y <= (log R)^(1/2)/((1/2) loglog R)`; otherwise `log Y > (1/2) loglog R` and `Y/log Y < 2Y/loglog R`.
QED

**Corollary 2.3 (near-extremal clusters sit at or above the lower edge).** Let `0 < eps <= 1`. If
`M >= max(M_1(C), eps log R/loglog R)` and `(log R)^(1/2) > 12.56/eps`, then `P_+(N') > (eps/12.56) log R`.
With (I1) and (I3), every near-extremal cluster therefore has

```text
(eps/12.56) log R  <  P_+(N')  <=  M^(4/eps - 3 + o(1)),        and   log R >= (3/2) M log M - (K'/2) M .
```

*Proof.* Corollary 2.2 gives `eps log R/loglog R <= M < 12.56 max((log R)^(1/2), Y)/loglog R`, and
`eps log R > 12.56 (log R)^(1/2)` by hypothesis, so the maximum is `Y` and `Y > (eps/12.56) log R`. The last
inequality is (I1) with `W' <= W = 2 log R`. QED

So, in units of `M`, the window of the task is `P_+(N') in [c_eps M log M, M^(2-delta)]`; below `c_eps M log M`
nothing near-extremal exists, and Corollary 2.2 is a (cheap) restricted growth theorem there.

**Remark 2.4 (why counting stops at exponent 1).** Theorem 2.1 is a statement about how many primes a cluster
needs: `omega(N') > M/4 (1 - o(1))`. Near-extremal clusters also have `omega(N') <= (2/eps + o(1)) M` (I3), and
for every `a > 1`, `pi(M^a; 4, 1) ~ M^a/(2a log M)` is much larger than `M`. So in the window
`[M^(1+eta), M^(2-delta)]` there are `M^(1 + eta - o(1))` times more admissible primes than a cluster uses, and no
statement that depends on the primes of `N'` only through their number and sizes can exclude near-extremal
clusters there. Any window theorem must say **which** primes can occur together: through their angles (the large
sieve, §§3-6) or their residues (P3-type data, overview §6.3).

*Exact check* (`prove_checks/check_thm21.py`): Theorem 2.1(1)-(2) on every window of `M >= 5` consecutive points
with chord constant `<= 9` on the record circles of `growth/families.md`; see §7.

## 3. What the large sieve says about actual clusters in the window (proved)

Notation: `e(t) = e^(2 pi i t)`. For a Gaussian prime `pi_j` put `x_j = 2 theta_j/pi in R/Z`, so that the Hecke
character `lambda_k(pi_j) = e^(4ik theta_j) = e(k x_j)`. An arc of length `delta` in `x` is an arc of length
`pi delta/2` in `theta mod pi/2`.

**Lemma 3.1 (sector count).** Let `Y >= 1` and let `I` be an arc of `R/(pi/2)Z` of length `eta in (0, pi/2]`. The
number of associate classes of `alpha in Z[i] \ {0}` with `Norm alpha <= Y` and `arg alpha mod pi/2 in I` is at most

```text
eta Y/2 + ((2 + eta)/sqrt 2) Y^(1/2) + pi/2   <=   eta Y/2 + 2.53 Y^(1/2) + 1.58 .
```

*Proof.* Lift `I` to a half-open interval `J` of length `eta <= pi/2`. Each associate class has exactly one member
with argument in any half-open interval of length `pi/2`, so at most one with argument in `J`. These members lie in
the closed sector `K = {r e^(i psi) : 0 <= r <= Y^(1/2), psi in closure(J)}`, which is convex (`eta <= pi`), with
area `eta Y/2` and perimeter `(2 + eta) Y^(1/2)`. The squares `alpha + [-1/2, 1/2)^2` are pairwise disjoint and lie in
`K + B(0, 2^(-1/2))`, whose area is `eta Y/2 + (2 + eta) Y^(1/2) 2^(-1/2) + pi/2` by Steiner's formula. QED

**Lemma 3.2 (the large sieve is controlled by arc counts).** Let `x_1, ..., x_n in R/Z`, `delta in (0, 1)`, and
suppose every half-open arc of length `delta` contains at most `m >= 1` of the points (with multiplicity). Then for
every integer `K >= 1` and complex `b`,

```text
sum_(|k| <= K) | sum_n b_n e(k x_n) |^2  <=  2m (2K + delta^(-1)) ||b||_2^2 .
```

Conversely, if `t` of the points lie in a closed arc of length `1/(8K)`, then the optimal constant `Lambda` in
`sum_(|k| <= K) |sum b_n e(k x_n)|^2 <= Lambda ||b||^2` satisfies `Lambda >= (2K + 1) t/sqrt 2`.

*Proof.* Sort the points cyclically, `x_0 <= ... <= x_(n-1)` in `[0, 1)`. If `n <= m`, Cauchy-Schwarz over the `n`
singletons gives `sum_k |sum b e|^2 <= n (2K + 1) ||b||^2`, which is within the claim. Otherwise write
`n = tm + rho`, `t >= 1`, `0 <= rho < m`. Form the `m` classes `{c, c + m, ..., c + (t-1)m}` (`0 <= c < m`) and the
`rho` singletons `{tm}, ..., {n-1}`: fewer than `2m` classes. In a class, consecutive members `x_i, x_(i+m)` satisfy
`x_(i+m) - x_i >= delta`, since otherwise the half-open arc `[x_i, x_i + delta)` would contain the `m + 1` points
`x_i, ..., x_(i+m)`. For the cyclic gap from the last member `L = c + (t-1)m` forward to `x_c` (only when `t >= 2`):
the forward arc from `x_L` to `x_c` contains `x_L`, `x_c` and the `n - 1 - L + c = m + rho - 1` points with indices
`L+1, ..., n-1, 0, ..., c-1`, so `m + rho + 1 >= m + 1` points; if it were shorter than `delta`, the arc
`[x_L, x_L + delta)` would contain them all. So all cyclic gaps inside a class are `>= delta`, and each class is
`delta`-separated in `R/Z`. Apply (Q2) to each class with the `L = 2K + 1` frequencies `-K, ..., K`
(`L - 1 + delta^(-1) = 2K + delta^(-1)`), and Cauchy-Schwarz over the `< 2m` classes:
`|sum_c S_c(k)|^2 <= 2m sum_c |S_c(k)|^2`. Summing over `k` gives the claim.
Converse: put `b_n = 1` on the `t` points of the arc and `0` elsewhere. For `|k| <= K` and two such points,
`|k (x_n - x_n')| <= 1/8`, so `Re e(k(x_n - x_n')) >= cos(pi/4)`, and `|sum b_n e(k x_n)|^2 >= t^2/sqrt 2` for each of
the `2K + 1` values of `k`. Hence `(2K + 1) t^2/sqrt 2 <= Lambda t`. QED

So, up to the absolute factor `16 sqrt 2 < 23` (an arc of length `1/(2K)` is covered by four arcs of length `1/(8K)`,
so `Lambda <= 32 K m''` and `Lambda >= 2K m''/sqrt 2` with `m''` the maximal count in an arc of length `1/(8K)`), **the large sieve for a set of angles at frequency range `K` says exactly how many of them can lie in an
arc of length `~1/K`, and nothing else.**

**Proposition 3.3 (angular spreading of an actual cluster; proved).** Let a normalised cluster have prime set
`{p_1, ..., p_s}`, all `<= Y`, and fixed primes `pi_j | p_j`. Then:

1. For every arc `I` of `R/(pi/2)Z` of length `eta <= pi/2`:
   `#{j : theta_j in I} <= min(s, eta Y/2 + 2.53 Y^(1/2) + 1.58)`.
2. For every `K >= 1` and complex `b`:
   `sum_(|k| <= K) |sum_j b_j lambda_k(pi_j)|^2 <= (pi Y + 20.24 K Y^(1/2) + 12.64 K) ||b||^2`.
3. In the window, `Y <= M^(2-delta)` and `s > M/5` for `M >= M_1(C)` (Corollary 2.2), so `s/Y^(1/2) > M^(delta/2)/5`:
   item 2 at `K = Y^(1/2)` beats the trivial bound `(2K + 1) s ||b||^2` by a factor `>= c M^(delta/2)`, and item 1
   forbids more than `eta Y/2 + 2.53 Y^(1/2) + 1.58` cluster angles in an arc of length `eta`. This is the
   non-triviality of the large sieve that the task refers to.

*Proof.* The `pi_j` are pairwise non-associate (distinct norms), so item 1 is Lemma 3.1. Item 2: in `x`-coordinates
an arc of length `delta = 1/(2K)` is a `theta`-arc of length `pi/(4K)`, so it contains at most
`m = pi Y/(8K) + 2.53 Y^(1/2) + 1.58` points by item 1; Lemma 3.2 gives `2m(2K + 2K) = 8Km`. Item 3 is arithmetic. QED

This is all the large sieve (Gallagher, Montgomery-Vaughan, Hecke-character versions alike) can say about the
cluster's own prime angles: an upper bound on how many of them lie in an arc, at every scale `>= 1/Y`. Note that
item 1 is the lattice count, so it holds for the angles of **any** `s` non-associate Gaussian integers of norm
`<= Y`; it does not use primality, nor the cluster.

## 4. Spreading never conflicts with the short-arc condition (proved)

A growth proof through the large sieve would have to show: *the short-arc condition forces the cluster's prime
angles to concentrate* (more than `~ eta Y + Y^(1/2)` of them in some arc of length `eta`), contradicting
Proposition 3.3. This section proves that the short-arc condition never forces concentration, for every profile
whose coordinates are pairwise independent in the sense below. The profile, the cut slack and the chain lemma
(`two_thirds/proof.md` §4) are combinatorial data about the exponents `a_xj`; the angles enter the problem only
through `B phi mod pi/2`, and `Real(Delta)` is invariant under a large group that moves every coordinate.

**Setting.** Fix a profile on `s` primes, `T = R/(pi/2)Z`, and the integer `(M-1) x s` matrix `B` with rows
`2(a_x - a_1)`, `x = 2, ..., M`. Put

```text
H = { h in T^s : B h in (pi/2) Z^(M-1) },        Cok = Z^s / B^T Z^(M-1) .
```

`H` is a closed subgroup of `T^s`; its Pontryagin dual is `Cok` (a character `v in Z^s` of `T^s` is trivial on `H` iff
`v in B^T Z^(M-1)`, the annihilator of `H`). Let `D_j in {1, 2, ..., infinity}` be the order of `e_j` in `Cok`, and call
`i != j` *independent* if `<e_i> cap <e_j> = 0` in `Cok`.

**Lemma 4.1.** (a) `Real(Delta) + H = Real(Delta)`. (b) For Haar-random `h in H`, `h_j` is Haar-distributed on the
subgroup of `T` of order `D_j` (on `T` if `D_j = infinity`); hence for any arc `I` of length `eta` and any `c in T`,
`Prob(c + h_j in I) <= 2 eta/pi + 1/D_j`. (c) If `i, j` are independent, `h_i` and `h_j` are independent.

*Proof.* (a) If `phi` realises the profile with `(theta_0, delta_x, k_x)` and `h in H`, then
`<2a_x - e, phi + h> = <2a_x - e, phi> + <2a_1 - e, h> + <2(a_x - a_1), h>`, and the last term is in `(pi/2)Z`; so
`phi + h` realises it with `theta_0 + <2a_1 - e, h>`, the same `delta_x`, and shifted `k_x`. Apply this to `-h` too.
(b) The image of `H` under the character `h -> h_j` is the closed subgroup of `T` dual to `<e_j> subset Cok`, i.e. of
order `D_j`, and the push-forward of Haar measure under a continuous surjective homomorphism of compact groups is
Haar measure. A coset of the subgroup of order `D_j` has at most `eta D_j/(pi/2) + 1` points in `I`. (c) The image of
`H` under `(h_i, h_j)` is dual to `<e_i, e_j> = <e_i> (+) <e_j>`, so it is the product of the two images, and the
push-forward of Haar measure is the product measure. QED

**Proposition 4.2 (spread points of the realisation set).** Let `J` be a set of `n_J` coordinates that are pairwise
independent and have `D_j >= D_0` (`D_0 = infinity` allowed), and let `L >= 1` with `2^L <= D_0`. For every
`phi in Real(Delta)` there is `h in H` such that `phi' = phi + h in Real(Delta)` and, for every arc `I` of `T` of
length `eta >= (pi/2) 2^(-L)`,

```text
#{ j : phi'_j in I }  <=  (16/pi) s eta + 4 (L s)^(1/2) + (s - n_J) .
```

Consequently, for every `K <= 2^(L-1)` and complex `b`,
`sum_(|k| <= K) |sum_j b_j e^(4ik phi'_j)|^2 <= (32 s + 32 K (L s)^(1/2) + 8K (s - n_J)) ||b||^2`.

*Proof.* For `l = 1, ..., L` let `F_l` be the `2^l` grid arcs `[(pi/2) i 2^(-l), (pi/2)(i+1) 2^(-l))`. For `I in F_l`,
`N_I = #{j in J : phi_j + h_j in I}` is a sum of pairwise independent indicators (Lemma 4.1(c)), so
`Var N_I <= E N_I <= n_J (2^(-l) + 1/D_0) =: E_l` (Lemma 4.1(b)). By Chebyshev,
`Prob(N_I >= E_l + t_l) <= E_l/t_l^2`. With `t_l = (2 L n_J (1 + 2^l/D_0))^(1/2)` the union over `l <= L` and
`I in F_l` has probability `<= sum_l 2^l E_l/t_l^2 = sum_l 1/(2L) = 1/2 < 1`. Fix `h` outside this event. Since
`2^l <= D_0`, every grid arc of level `l` then has `N_I < 2 n_J 2^(-l) + 2 (L n_J)^(1/2)`. An arc of length
`eta in ((pi/2)2^(-l-1), (pi/2)2^(-l)]` (`l <= L`; for `eta > pi/4` use `l = 0` and the trivial bound `n_J`) meets
at most two grid arcs of level `l`, so it holds fewer than `4 n_J 2^(-l) + 4 (L n_J)^(1/2) < (16/pi) n_J eta + 4(L n_J)^(1/2)`
points of `J`, and at most `s - n_J` others. `phi' in Real(Delta)` by Lemma 4.1(a). The large-sieve bound is
Lemma 3.2 with `delta = 1/(2K)` in `x = 2 phi/pi` (a `T`-arc of length `pi/(4K) >= (pi/2)2^(-L)`), with
`m <= 4s/K + 4(Ls)^(1/2) + (s - n_J)` and `8Km` as in Proposition 3.3. QED

**Corollary 4.3 (the large sieve cannot contradict the short-arc condition).** Let a profile on `s <= pi(Y; 4, 1)`
primes `<= Y` have all `s` coordinates pairwise independent with `D_j = infinity`. (This holds iff the row space of
`B` over `Q` contains no nonzero vector supported on at most two coordinates, because `B^T Z^(M-1)` has finite index in
the integer points of that row space; it needs `rank B <= s - 2`.) Then for
every `Delta > 0`, `Real(Delta)` contains angle vectors `phi'` whose large-sieve constant at every `K <= Y/4` is at
most `32 s + 39 K (s log Y)^(1/2)`. For `Y >= e^13` this is at most `3 (pi Y + 20.24 K Y^(1/2))`, three times the bound
that the **actual** prime angles satisfy (Proposition 3.3(2)); and at every scale `eta >= pi/Y` the arc counts of
`phi'` are at most `(16/pi) s eta + 4.9 (s log Y)^(1/2)`, within an absolute factor of Proposition 3.3(1). Hence no
large-sieve or sector-count statement relaxed by a fixed factor (published large sieves carry non-optimal absolute
constants) is violated by **all** points of `Real(Delta)`: spreading alone can never certify that a profile has no
realisation by genuine angles.

*Proof.* Take `J` = all coordinates, `D_0 = infinity`, `L = floor(log_2 Y)`; then `L <= 1.4427 log Y` and
`2^(L-1) >= Y/4 >= K`. Proposition 4.2 gives the constant `32 s + 32 K (Ls)^(1/2) <= 32 s + 38.5 K (s log Y)^(1/2)`
and the arc counts `(16/pi) s eta + 4 (Ls)^(1/2) <= (16/pi) s eta + 4.9 (s log Y)^(1/2)` for `eta >= (pi/2) 2^(-L)`,
which covers `eta >= pi/Y` because `2^(-L) < 2/Y`.
With `s <= pi(Y) < 1.26 Y/log Y`: `32 s < 40.4 Y/log Y <= pi Y` once `log Y >= 12.9`, and
`39 K (s log Y)^(1/2) < 39 K (1.26 Y)^(1/2) < 44 K Y^(1/2) < 3 * 20.24 K Y^(1/2)`. QED

**Remark 4.4 (profiles with pinned coordinates).** If some `D_j` is small, the coordinate `phi_j` is confined to a
coset of a cyclic group of order `D_j`, up to the realisation tolerance. In the extreme case `rank B = s` (which is
typical when `s <= M - 1`), `H` is finite and, **provided** `Delta` is smaller than a constant `c(B) > 0` depending
on the profile (the distance from `(pi/2)Z^(M-1) \ B R^s` to `B R^s`, divided by `(M-1)^(1/2)`), every
`phi in Real(Delta)` lies within `||B^+||_(inf) Delta` of `H`. Then each `theta_j` lies within that distance of
`(pi/2)(1/D_j) Z`, and the monomial gap (I5) at `v = D_j e_j` forces
`D_j log p_j >= W/2 - 2 log(D_j ||B^+|| C) - O(1)`, i.e. `D_j >~ W/(2 log p_j)`. The large sieve in the same
situation only asks that the `s` angles occupy at least `s/(2.53 Y^(1/2) + 1.58)` distinct points of these cosets
(Proposition 3.3(1) at `eta -> 0`), a weaker requirement in the window (by (I1), `W/(2 log p_j) >= (3/(2a) - o(1)) M`,
against `s/Y^(1/2) = O(M^(1 - a/2))` when `s = O(M)` and `Y = M^a`, `1 < a < 2`). So in the pinned case the large sieve adds nothing
to the monomial gaps. This remark is not used below; its threshold `c(B)` is profile-dependent and is **not**
claimed uniform in `R`.

## 5. The graft barrier inside the window (proved, with a conditional corollary)

`round10/analytic.md` Theorem B needs the grafted primes at height `>= M^A` with `A > 2` because its large-sieve
estimate (Lemma 3.2(1) there) is the worst case `||D|| <= (8 s (2K + 1))^(1/2) ||a||`, which is attained only when
the `s` grafted angles are concentrated. In the window `s >> X^(1/2)`, so that estimate is useless. With
Lemma 3.2 here, the worst case is replaced by the arc counts, and the window ceases to be special for the large
sieve.

**Definition 5.1.** A set of `s` angles `psi_1, ..., psi_s in T = R/(pi/2)Z` is `L`-spread if (i) every arc of `T` of
length `eta >= (pi/2) 2^(-L)` contains at most `(16/pi) s eta + 4 (L s)^(1/2)` of them, and (ii) the `2s` angles
`psi_j, pi/2 - psi_j` are pairwise at distance `>= (pi/2) 2^(-L)` in `T`.

Genuine Gaussian prime angles of norm `<= X` satisfy (ii) with `(pi/2)2^(-L) <= 1/X` by the lattice cross-gap (I4)
(also for the pair `theta_p, pi/2 - theta_q`, which is `arg pi_p` against `arg(i conj pi_q)`). Proposition 4.2 produces
points of `Real(Delta)` satisfying (i).

**Proposition 5.2 (spread grafts satisfy the robust large sieve; proved).** Fix `c > 0`. Let `S` be a set of `s`
split primes in `[P_0, P]`, `phi in T^S`, and `theta^S` the graft (analytic.md Def 3.1: `theta^S(varpi_p) = phi_p`,
`theta^S(i conj varpi_p) = pi/2 - phi_p` for `p in S`, genuine elsewhere). Suppose that both `(phi_p)_(p in S)` and the
genuine `(theta_p)_(p in S)` are `L`-spread. Then for every `X >= 2`, `K >= 1` and every inequality
`sum_(|k|<=K) |sum_(Norm varpi <= X) a_varpi lambda_k(varpi)|^2 <= Lambda ||a||^2` with `Lambda >= c(X + K X^(1/2))`
that holds for the genuine angle function, the graft satisfies it with `Lambda` replaced by `(1 + t^(1/2))^2 Lambda`,

```text
t  <=  max( (256/c) (s/P_0 + (L s/P_0)^(1/2)) ,  16/(c P_0^(1/2)) ) .
```

*Proof.* If `X < P_0` nothing changes. Let `X >= P_0` and `D_k = sum_(varpi above S, Norm <= X) a_varpi (lambda^S_k - lambda_k)(varpi)`.
By Minkowski's inequality in `l^2(k)`, `(sum_k |sum a lambda^S_k|^2)^(1/2) <= Lambda^(1/2) ||a|| + (sum_k |D_k|^2)^(1/2)`, and
`(sum_k |D_k|^2)^(1/2) <= U_graft + U_gen`, where `U` is the large-sieve norm of the `2 s_X <= 2s` grafted (resp. genuine)
points above `S cap [1, X]`, in the coordinates `x = 2 psi/pi in R/Z`.
*Case `K <= 2^(L-1)`.* A `x`-arc of length `1/(2K)` is a `T`-arc of length `pi/(4K) >= (pi/2) 2^(-L)`, so by (i) applied
to `psi` and to `pi/2 - psi` it contains at most `m = 2((16/pi) s pi/(4K) + 4(Ls)^(1/2)) = 8s/K + 8(Ls)^(1/2)` points.
Lemma 3.2 with `delta = 1/(2K)` gives `U^2 <= 8 K m ||a||^2 = 64 (s + K (Ls)^(1/2)) ||a||^2` for each family, so
`(sum |D_k|^2)^(1/2) <= 16 (s + K(Ls)^(1/2))^(1/2) ||a||`, and `t = 256 (s + K (Ls)^(1/2))/Lambda`; use
`Lambda >= c(X + K X^(1/2))` and `X >= P_0`.
*Case `K > 2^(L-1)`.* By (ii) the `2 s_X` points of each family are `2^(-L)`-separated in `R/Z`, so (Q2) gives
`U^2 <= (2K + 2^L) ||a||^2 <= 4K ||a||^2`; then `(sum|D_k|^2)^(1/2) <= 4 K^(1/2) ||a||` and
`t = 16K/Lambda <= 16/(c X^(1/2))`. QED

**Proposition 5.3 (prime sums and boxes in the window; proved).** Theorem B's estimates for (R2) and (R3)
(`analytic.md` §3) hold verbatim for any graft on `s` primes `>= P_0`: the relative change is
`<= 4 s log X_1/(c X_1^gamma)` for (R2) and `<= 2 s/(c x_B^gamma)` for (R3), which tend to `0` when
`s log P_0 = o(P_0^gamma)`. In the window with `P_0 >= M^(1 + eta)` and `s = O(M)` this holds for every
`gamma > 1/(1 + eta)`. *Proof.* `analytic.md` Theorem B proof, items (R2), (R3), with `M^A` replaced by `P_0`; spreading
is not used. QED

**Corollary 5.4 (conditional cap for the window).** Fix `eta, delta in (0, 1)`, `c > 0`, `gamma > 1/(1 + eta)`, and
`kappa > 0`. Suppose that for an unbounded set of `M` there are angle models (class `A_0` of the overview: profile,
every monomial gap, Haar-generic properties) with `M` points, constant `C_0`, `W <= kappa M log M`, `s = O(M)` primes
all in `[M^(1+eta), M^(2-delta)]`, and with both their angles and the genuine angles of their primes
`L`-spread for some `L = O(log M)` with `(pi/2)2^(-L) <= M^(-2)`. Then no argument that uses only `A_0` together
with `(c, gamma)`-robust prime-supported ensemble inputs (Hecke large sieve over primes, Hecke prime sums, sector and
box counts, with constants fixed in advance) can prove `M = o(log R/loglog R)` for clusters in the window.

*Proof.* Apply Propositions 5.2 and 5.3 to the grafts: `s/P_0 <= O(M^(-eta))`, `(Ls/P_0)^(1/2) = O(M^(-eta/2) (log M)^(1/2))`,
`P_0^(-1/2) <= M^(-1/2)`, and `s log P_0/P_0^gamma -> 0`. So every robust input holds for the graft up to `1 + o(1)`, and
the angle models satisfy `A_0` at `W <= kappa M log M`. This is the logic of `analytic.md` Corollary B.1. QED

**What this says.** Inside the window the large sieve is no more powerful than in Theorem B's regime, **as soon as
the fake's angles are spread**; and spread points are always available in the realisation set (Proposition 4.2).
The hypothesis of Corollary 5.4 that is not known to hold is the **existence of gap-admissible angle models with all
primes `<= M^(2-delta)`**. That, not the large sieve, is what keeps the window open. Section 6 isolates why.

*Caveats.* (1) Corollary 5.4 is conditional on the existence of window fakes and is therefore not a cap. (2) The
precision of the robust inputs is fixed before the target constant (the round-9 R1 caveat applies unchanged).
(3) Proposition 5.2 needs the **genuine** angles of `S` to be spread; for an arbitrary `S` they only satisfy the
lattice bound of Proposition 3.3(1), which is not `o(Lambda)`. A fake construction would choose `S`, and a random
choice of `S` among the admissible primes is spread with high probability if Gaussian prime angles in short norm
intervals obey an upper-bound sieve in polar boxes of angle `>= P^(-1/2)` (an assumption, not checked here), so this is a
construction issue, not a new obstruction.

## 6. The actual obstruction: monomial gaps at many small primes

Every refereed fake is built the same way: fix a profile, take the angle vector Haar-generic on the realisation
torus, and prove gap-admissibility by a union bound over characters `v` (`round8/lw.md` Thm 6.1; `round5/construct.md`).
In the window that union bound fails, and the threshold at which it fails is the same `s ~ P^(1/2)` at which the
large sieve becomes non-trivial.

For `v in Z^s \ {0}` let `G_v` be the gap event of (I5) at `v`, and `Bad_v` its complement:
`dist(<v, phi>, (pi/2)Z) < arcsin(e^(-w(v)/2))` or `dist(<v, phi>, pi/4 + (pi/2)Z) < arcsin(e^(-w(v)/2)/sqrt 2)`.
`Bad_v = Bad_(-v)`.

**Proposition 6.1 (first-moment threshold; proved).** Let `S` be `s` split primes, all `<= P`.

1. If `phi` is Haar-distributed on `T^s` (or on any closed subgroup on which the character `v` is non-trivial and has
   image `T`), then `Prob(Bad_v) = (4/pi)(arcsin(e^(-w(v)/2)) + arcsin(e^(-w(v)/2)/sqrt 2)) >= 2.17 e^(-w(v)/2)` (for
   `w(v) >= log 5`, where the two excluded arcs are disjoint).
2. `sum_(v != 0) e^(-w(v)/2) = prod_j (1 + 2/(p_j^(1/2) - 1)) - 1 >= exp(sum_j 2/(p_j^(1/2) + 1)) - 1 >= exp(2s/(P^(1/2) + 1)) - 1`.
   Hence the expected number of violated gap events (counting `v, -v` once) of a Haar-random `phi in T^s` is
   `>= 1.08 (exp(2s/(P^(1/2) + 1)) - 1)`.
3. For Haar-random `phi in T^s` the single-prime gaps (`v = e_j`) are independent, and
   `Prob(all single-prime gaps hold) <= exp(-2.17 sum_j p_j^(-1/2)) <= exp(-2.17 s/P^(1/2))`.
4. In the window (`P <= M^(2-delta)`, `s > M/5` for `M >= M_1(C)` by Corollary 2.2), the expected number in item 2
   is `>= exp(0.4 M^(delta/2) - 0.4) - 2` and the probability in item 3 is `<= exp(-0.43 M^(delta/2))`. All refereed fakes
   instead have `sum_j p_j^(-1/2) <= 1/32` (`lw.md` 6.1: `s = 2(M-1)` primes `>= (64M)^2`).

*Proof.* 1. Under the stated hypothesis `<v, phi> mod pi/2` is uniform on `T` (push-forward of Haar under a
surjective character), and the two excluded sets are arcs of total length `2 arcsin(e^(-w/2)) + 2 arcsin(e^(-w/2)/sqrt2)`,
disjoint when `e^(-w/2) <= 5^(-1/2)` (their half-lengths then sum to at most `arctan(1/2) + arctan(1/3) = pi/4`); divide by `pi/2` and use
`arcsin x >= x`: `(4/pi)(1 + 2^(-1/2)) = 2.173...`. 2. `w` is additive over coordinates, so the sum factorises; for one
prime `sum_(m in Z) p^(-|m|/2) = 1 + 2/(p^(1/2) - 1)`; then `log(1 + x) >= x/(1 + x)` with `x = 2/(p^(1/2) - 1)` gives
`log(1 + x) >= 2/(p^(1/2) + 1)`. The expected number of violated events is `(1/2) sum_(v != 0) Prob(Bad_v)`
(every `w(v) >= log 5`), and `(1/2) 2.17 > 1.08`. 3. The `e_j`-events depend on disjoint coordinates; use item 1
and `1 - x <= e^(-x)`. 4. `2s/(P^(1/2) + 1) >= (2M/5)/(M^(1 - delta/2) + 1) >= 0.4 M^(delta/2) - 0.4`, `1.08(e^x - 1) >= e^x - 2`, and
`2.17 s/P^(1/2) >= 2.17 (M/5) M^(delta/2 - 1) > 0.43 M^(delta/2)`; the `lw.md` primes
satisfy `p^(-1/2) <= 1/(64M)`. QED

**Reading.** A first-moment argument needs `sum_j p_j^(-1/2) = O(1)`, i.e. `s = O(P^(1/2))`, i.e. `P >~ M^2`; with the
equal-weight windows of the known constructions this becomes `P >~ r^2 log r` (`round10/README.md` §3). The same
quantity `s/P^(1/2)` measures the non-triviality of the large sieve for the cluster's angles (Proposition 3.3(3)).
This is not a coincidence of notation: both say that `s` Gaussian-prime-like angles at height `P` are more numerous
than the `~ P^(1/2)` lattice directions that any arc of length `~P^(-1/2)` can hold. **But the large sieve is a
statement about single angles (Lemma 3.2), and it is satisfied by spread points of the realisation set
(Proposition 4.2), whereas the gap events are statements about all integer combinations `<v, phi>`.** In the
window, the obstruction to fakes is the latter.

**Proposition 6.2 (lattice angles collapse; proved, a restatement of `analytic.md` Thm D).** Gap-admissible angle
vectors of the form `phi_j = arg beta_j`, `beta_j in Z[i]` primitive, `Norm beta_j <= p_j`, with the `beta_j` and their
conjugates multiplicatively independent, exist in every window (take `beta_j = pi_j`, or any Gaussian primes of
smaller norm), but every realisation of a profile by such a vector is the shadow of an actual cluster on
`x^2 + y^2 = prod Norm(beta_j)^(e_j) <= N`. *Proof.* `analytic.md` Theorem D, second paragraph: `y_x = i^(-k_x) prod
beta_j^(a_xj) conj(beta_j)^(e_j - a_xj)` are distinct Gaussian integers of one norm `N'' <= N` on an arc of length
`<= C N''^(1/4)`. (The gap (I5) for such `phi` is the lattice bound `|Im A| >= 1` for the non-real Gaussian integer
`A = prod beta_j^(v_j^+) conj(beta_j)^(v_j^-)`.) QED

So a window fake must use **non-lattice** angles, and for those the only available existence technique (Haar plus
union bound) fails by Proposition 6.1. Whether non-lattice gap-admissible points of `Real(Delta)` exist with primes
`<= M^(2-delta)` is the precise open question.

**Heuristic 6.3 (labelled heuristic; evidence in §7.3).** Model the events `Bad_v` for Haar-random `phi` as a
Poisson process with total intensity `E = (1/2) sum_(v != 0) Prob(Bad_v)`. Then `Prob(phi is gap-admissible) ~ e^(-kappa E)`
with `kappa` of order 1 (the data of §7.3 give `kappa ~ 1.1-1.3` for `E <= 10`) and `E ~ exp(2 sum_j p_j^(-1/2))`: doubly exponentially small in `s/P^(1/2)`, but positive. On the realisation set the
events with `v` outside the relation lattice behave the same way, and those inside are the profile's margin
conditions, which the known profiles satisfy with room. If this model is right, gap-admissible window fakes exist
with Haar probability `exp(-exp(O(M^(1 - a/2))))` for primes near `M^a`, `1 < a < 2`, and the window is fake-free only
because first-moment methods cannot prove their existence. Nothing here proves this, and Theorem D shows that the
lattice points, which do exist, are not fakes.

## 7. Computation: is the window populated?

All of this is **evidence**; no proof above uses it. Scripts and outputs are in `prove_checks/`.

### 7.1 Exhaustive scan of smooth norms (`smoothscan.c`, `smoothscan2.c`)

*Method.* The scan of `growth/families.md` §6, restricted to smooth norms and pushed to `10^18`: every
`n <= X` whose prime factors are all `= 1 mod 4` and `<= B` (by the normalisation argument of families.md §6, the
optimal constant for `k` points is attained on such `n` with the content removed). The `d(n)` point classes are
angles mod `pi/2` stored as `u64` with `2^64 = pi/2`, built multiplicatively from 50-digit prime angles
(`mkangles.py`), sorted by radix sort; `C_k(n) = n^(1/4) * (least span of k cyclically consecutive angles)`, exactly as in
families.md. Fixed-point error per angle `< 10^(-15)` against spans `> 10^(-6)`. *Validation:* with `B = 400`,
`X = 10^10` the scan reproduces every entry of the referee's brute force `growth/referee_families_checks/brute_1e10.txt`
whose norm is 400-smooth (`k = 4, 5, 7, 8, 9, 10, 11, 12`; the brute-force entries for `k = 3, 6` have a prime factor
`> 400`: `9520005625 = 5^4 13 673 1741`, `1944545425 = 5^2 17 89 101 509`). Runs: `B = 200, X = 10^18`
(`40,715,511` norms, `1.86 * 10^10` point classes) and `B = 113, X = 10^18` (`4,656,319` norms), each a few minutes.
For each `k <= 24`, decade of `n` and largest prime `P_+(n)`, the best `C_k` is recorded.

*Best `C_k` per decade of `n`, over 200-smooth `n`: all, and restricted to `P_+(n) <= k^2`* (`sum2_B200_1e18.txt`):

| decade of n | k=7 all / P+<=k^2 | k=8 all / P+<=k^2 | k=9 all / P+<=k^2 | k=10 all / P+<=k^2 | k=11 all / P+<=k^2 |
| --- | --- | --- | --- | --- | --- |
| 1e6 | 6.07 / 6.60 | 7.53 / 7.64 | 9.26 / 9.26 | 9.91 / 9.91 | 12.88 / 12.88 |
| 1e8 | 5.27 / 6.52 | 7.68 / 9.35 | 9.71 / 10.50 | 10.33 / 12.65 | 16.46 / 16.46 |
| 1e9 | 6.50 / 7.48 | 7.03 / 8.27 | 7.88 / 9.35 | 8.24 / 10.20 | 16.17 / 16.17 |
| 1e10 | 7.12 / 10.16 | 8.32 / 11.16 | 10.20 / 13.61 | 11.31 / 15.22 | 16.50 / 16.50 |
| 1e11 | 6.61 / 10.69 | 7.52 / 10.09 | 11.56 / 11.56 | 12.97 / 15.19 | 18.76 / 18.76 |
| 1e12 | 7.34 / 13.45 | 8.90 / 9.06 | 11.51 / 11.51 | 11.68 / 11.68 | 16.14 / 16.14 |
| 1e13 | 7.12 / 13.86 | 7.92 / 14.71 | 12.29 / 17.21 | 16.08 / 16.56 | 19.70 / 19.70 |
| 1e14 | 7.87 / 20.72 | 9.67 / 19.37 | 15.03 / 23.86 | 18.56 / 26.13 | 23.64 / 23.64 |
| 1e15 | 8.03 / 21.80 | 9.84 / 23.75 | 15.70 / 23.20 | 18.42 / 24.23 | 27.63 / 27.63 |
| 1e16 | 8.40 / 32.60 | 10.78 / 23.92 | 16.34 / 16.34 | 17.37 / 17.37 | 33.21 / 33.21 |
| 1e17 | 8.89 / 43.30 | 10.84 / 30.86 | 21.99 / 24.43 | 25.98 / 25.98 | 38.15 / 40.32 |

*Reading (evidence).*
1. Up to `10^18` the 200-smooth norms reproduce the record constants `5.2699` (`k = 7`), `7.8804` (`k = 9`) and
   `8.2393` (`k = 10`, both at `n = 1176852625`); the golden `k = 8` record `5.8978` is not in this run because
   `537606725 = 5^2 17 61 89 233` has `233 > 200`. The scan finds **nothing better at larger `n`**: for every `k >= 7`
   the best constant per decade **increases** with `n` beyond `10^9` (for `k = 10`: `8.24` at `10^9`, `~17-19` at
   `10^14-10^16`, `26` at `10^17`). Larger smooth norms carry more points (`d(n)` up to the thousands), but not more
   points on short arcs.
2. **The window at `delta = 0` is populated** at these sizes, but only trivially: the best clusters have
   `P_+(n)` between `k^2` and `k^3` for `k <= 8` (e.g. `k = 7`: `113 = 2.3 k^2`), about `k^2` for `k = 9, 10`
   (`113 = 1.13 * 10^2`), and below `k^2` for `k >= 11` because all scanned primes are `<= 200 < 11^2 + 80`. Restricting
   to `P_+ <= k^2` costs a factor up to `~1.5` in `C` at `k <= 10` and nothing at `k >= 11`.
3. **The extremal rate cannot be tested at these sizes.** The 2/3 inequality (I1) is vacuous while
   `3 log M <= K_0 ~ 14`, i.e. for `M <= 100` at `C <= 1`; here `M <= 11` and `W <= 41`. In units of the growth
   dictionary, the best clusters have `A = W'/(M log M)` between `0.9` and `2.7` (`check_thm21_out.txt`), which is a
   statement about the `O(M)` terms, not about the rate.

So the data are consistent with the window being populated by small clusters and with **no** growth in it; they
cannot distinguish "near-extremal clusters live in the window" from "they do not".

### 7.2 Record clusters (`check_thm21.py`, `cok_structure.py`)

Theorem 2.1(1)-(2) hold on all 267 windows (`k >= 5`, chord constant `<= 9`) of the seven record norms of
`growth/families.md` (no violation; exact Gaussian gcds). For the record clusters:

| cluster | `N'` | `s` | `P_+` | `P_+/k^2` | `sum_j p_j^(-1/2)` | `rank B` | Smith invariants of `B^T` | orders `D_j` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10 points, `C = 8.2386` (chord) | `5^3 13^2 17 29 113` | 5 | 113 | 1.13 | 1.247 | 5 | `2,2,2,2,24` | `24, 8, 6, 4, 24` |
| 7 points, `C = 5.2695` | `5^2 13^2 17 29 113` | 5 | 113 | 2.31 | 1.247 | 5 | `2,2,2,2,24` | `24, 8, 6, 4, 24` |
| golden 8 points, `C = 5.8974` | `5^2 17 61 89 233` | 5 | 233 | 3.64 | 0.989 | 4 | `2,2,2,2,0` | all `inf` |
| 5 points, `C = 2.2584` | `5^3 13 17 29 37 73` | 6 | 73 | 2.92 | 1.434 | 4 | `2,2,2,2` | all `inf` |

Two observations. (a) **Every record cluster has `sum_j p_j^(-1/2) ~ 1` or more**, i.e. it lives where the
first-moment construction of fakes fails (Proposition 6.1). Actual small clusters are made of few small primes; the
refereed fakes need many large ones. (b) The 7- and 10-point records are **pinned** (`rank B = s`, finite `H`): each
`theta_j` lies near `(pi/2) Z/D_j` with `D_j in {24, 8, 6, 4, 24}` (Remark 4.4), e.g. `theta_29 = arctan(2/5) = 0.3805`
against `pi/8 = 0.3927`. The golden 8-point cluster has
`rank B = s - 1`, so every pair of coordinates is dependent (10 of 10) and Corollary 4.3 does not apply; the 5-point
cluster (`rank B = s - 2`) has exactly one dependent pair (`5, 37`) out of 15, so Proposition 4.2 applies with `J` any
five of its six coordinates that avoid the pair (`cok_structure_out.txt`).

### 7.3 Monomial gaps for Haar-random angles (`poisson_gaps.py`, `smc_gaps.py`)

Question: how large is the Haar measure of the gap-admissible set
`Good(WMAX) = {phi in T^s : (I5) holds for every v with w(v) <= WMAX}` when `sum_j p_j^(-1/2)` is of order 1, i.e. at
the edge of the first-moment regime (Proposition 6.1)? `E` is the Poisson intensity
`sum_(v up to sign, w(v) <= WMAX) Prob(Bad_v)`. Direct Monte Carlo (`2 * 10^5` samples) and sequential Monte Carlo
(`2 * 10^4` particles, coordinates added one at a time):

| primes | `sum p^(-1/2)` | `WMAX` | `#v` | `E` | `log P_good` | `log P_good / (-E)` |
| --- | --- | --- | --- | --- | --- | --- |
| 101, 109, 113, 137, 149, 157 | 0.537 | 24 | 1039 | 2.09 | direct `log 0.0971 = -2.33`; SMC `-2.35` | 1.12 |
| 29, 37, 41, 53, 61, 73, 89, 97 | 1.096 | 22 | 6786 | 8.52 | direct `log(2.0 * 10^(-5)) = -10.8 +- 0.7`; SMC `-11.26` | 1.32 |
| 53, 61, 73 | 0.382 | 22 | 115 | 1.26 | SMC `-1.48` | 1.18 |
| first 5 split primes `>= 53` | 0.590 | 22 | 759 | 2.46 | SMC `-2.91` | 1.18 |
| first 7 | 0.785 | 22 | 2381 | 4.11 | SMC `-4.96` | 1.21 |
| first 9 | 0.965 | 22 | 4945 | 6.18 | SMC `-7.50` | 1.22 |
| first 10 | 1.047 | 22 | 6587 | 7.33 | SMC `-8.94` | 1.22 |

| first 11 | 1.126 | 22 | 8619 | 8.60 | SMC `-10.50` | 1.22 |
| first 12 | 1.203 | 22 | 11131 | 9.96 | SMC `-12.16` | 1.22 |

Tail intensities beyond `WMAX` (upper bounds from the product formula, `tail_bound.py`, using
`Prob(Bad_v) <= 2.18 e^(-w/2)` for `w >= 10`): `0.012` (first row), `0.287` (second), and `0.001, 0.010, 0.073, 0.286,
0.480, 0.753, 1.115` for the rows with the first 3, 5, 7, 9, 10, 11, 12 split primes `>= 53`. So the last two rows
miss up to about 10% of the intensity.

The genuine angles `theta_p` pass every enumerated gap (with a relative tolerance `10^(-9)` for the equality cases
`|Im A_v| = 1`, e.g. `101 = 10^2 + 1`), as Lemma 5.1 requires. *Reading (evidence).* Over this range the measure of
the gap-admissible set behaves like `exp(-kappa E)` with `kappa ~ 1.2`: a Poisson law with mildly positive
correlation, positive and smooth in `E`, with no sign of a collapse to the lattice points. Since
`E ~ exp(2 sum p_j^(-1/2))` (Proposition 6.1(2)), this supports Heuristic 6.3 in the form
`Prob(Good) ~ exp(-kappa exp(2 sum_j p_j^(-1/2)))`, doubly exponentially small but positive. The experiment reaches
only `sum p^(-1/2) ~ 1.1`, truncates at `w(v) <= WMAX` (the tail intensity beyond `WMAX` is listed above), and
says nothing about the realisation set or about `M -> infinity`.

### 7.4 Sanity checks of §3-§4 (`check_ls.py`)

(a) Genuine split primes `p <= Y` (`Y = 2 * 10^4`: 1125 primes; `Y = 10^5`: 4783 primes, random subset of 2500 for the
eigenvalue): the largest arc count at scale `pi/(4K)` is far below the lattice bound of Proposition 3.3(1)
(`K = Y^(1/2)`: 17 against 415, and 31 against 926), and the optimal large-sieve constant (largest eigenvalue of the
Dirichlet-kernel Gram matrix) lies between the two bounds of Lemma 3.2 in all six cases (e.g. `Y = 10^5`, `K = 316`:
`4924 <= Lambda = 9069 <= 42976`). (b) Proposition 4.2 demo: a random 0/1 profile with `M = 40` rows on the first 60
split primes in `[137, 881]` (`rank B = 39`, `dim H = 21`); three random points of `phi_0 + H` (with `phi_0 ~ 10^(-9)`)
realise all 40 angles within `2 * 10^(-8)` of each other mod `pi/2`, while the 60 coordinates are spread: at most 21-24,
8-9 and 4 coordinates in arcs of length `pi/8`, `pi/32`, `pi/128`, against the bounds `164, 92, 83` of Proposition 4.2.
(The demo approximates Haar measure on `H` by a long linear flow.)

## 8. Where the attempt fails, precisely; what a window theorem would need

### 8.1 The attempted argument and its break points

The intended proof: (i) in the window the large sieve at `K ~ P^(1/2)` is non-trivial for the cluster's prime
angles; (ii) the short-arc condition, through the cut/chain structure, forces those angles to concentrate;
(iii) contradiction unless `W/(M log M) -> infinity`. Step (i) is true (Proposition 3.3(3)). The argument breaks at
(ii), for four proved reasons.

* **(F1) The large sieve is spreading and nothing else.** For any set of angles, the optimal large-sieve constant at
  frequency range `K` is determined, up to an absolute factor `< 23`, by the maximal number of angles in an arc of
  length `~1/K` (Lemma 3.2, both directions). For the cluster this is the lattice count
  `eta Y/2 + O(Y^(1/2))` (Proposition 3.3(1)), which holds for **any** non-associate Gaussian integers of norm
  `<= Y`, prime or not, in a cluster or not.
* **(F2) The short-arc condition does not see single angles.** It constrains `phi` only through `B phi mod pi/2`, and
  the realisation set is invariant under the group `H` (Lemma 4.1). When the coordinates are pairwise independent
  in `Cok`, Haar-random points of `phi + H` are spread at every scale down to `1/Y` with constants comparable to the
  genuine ones (Proposition 4.2, Corollary 4.3). The cut slack `Q`, the level collisions and the chain lemma are
  functions of the exponent profile alone (`two_thirds/proof.md` §§2-4: no angle enters `F*`), so they cannot supply
  the missing coupling. When coordinates are pinned (`rank B = s`, the typical case for `s <= M - 1`, and the case of
  the 7- and 10-point records, §7.2), each `theta_j` sits near `(pi/2)Z/D_j`, and the monomial gap already forces more
  than the large sieve does (Remark 4.4).
* **(F3) Consequently the window is not special for robust analytic inputs once fakes are spread.** Theorem B extends
  to grafts with `s >> X^(1/2)` (Propositions 5.2, 5.3), so `A_0` plus every robust prime-supported input is capped in
  the window **if** spread gap-admissible window fakes exist (Corollary 5.4).
* **(F4) What actually distinguishes the window is the monomial gaps.** The first-moment construction of fakes needs
  `sum_j p_j^(-1/2) = O(1)`, which fails in the window by a factor `M^(delta/2)` (Proposition 6.1); lattice angles are
  gap-admissible but collapse to actual clusters (Proposition 6.2). The threshold `s ~ P^(1/2)` is the same as for the
  large sieve, but the large sieve sees single angles while the gaps see every integer combination.

The failure is therefore **not** "the cluster uses too few of the available angles for the large sieve to bite"
(it uses enough: Proposition 3.3(3)). It is that the large sieve bites only on concentration, and nothing forces
concentration.

### 8.2 The reduction that remains

**Proposition 8.1 (proved; the window question for `A_0`).** Fix `delta in (0, 1)` and `C > 0`. Let
`W_min^win(M)` be the infimum of `W` over angle models (overview Def 5.2: profile, realisation, every monomial gap)
with `M` points, constant `C`, and all primes `<= M^(2-delta)`.
1. If `W_min^win(M)/(M log M) -> infinity`, then actual clusters with `P_+(N') <= M^(2-delta)` satisfy
   `M = o(log R/loglog R)`.
2. If `W_min^win(M) <= kappa M log M` along an unbounded set of `M`, attained by models whose angles and whose
   primes' genuine angles are `L`-spread with primes `>= M^(1+eta)`, then no argument from `A_0` and robust
   prime-supported analytic inputs (large sieve included) proves item 1's conclusion (Corollary 5.4).

*Proof.* 1. Every normalised actual cluster is an angle model with `phi = theta` and its own primes
(overview §5.1), and `W = log N >= log N'`; then `2 log R >= W >= f(M) M log M` with `f -> infinity`; if `M >= (log R)^(1/2)` this gives
`M <= 4 log R/(f(M) loglog R)`, and otherwise `M <= (log R)^(1/2)`; either way `M = o(log R/loglog R)`. 2. Corollary 5.4. QED

So, for every input class considered in this round, **a restricted growth theorem in the window is equivalent (up to
the spreading proviso) to the non-existence of gap-admissible angle models with small primes at
`W = O(M log M)`**: a question about real vectors satisfying the Diophantine inequalities (I5) for all
`v in Z^s`, with `sum_j p_j^(-1/2) -> infinity`. The large sieve does not enter it.

### 8.3 What each outcome would need, and the evidence

* *Growth in the window via `A_0`* needs `W_min^win/(M log M) -> infinity`: every point of `Real(Delta)` with small
  primes violates some monomial gap unless it is a lattice point. The evidence of §7.3 points the other way: on the
  full torus the gap-admissible set has Haar measure `~ exp(-1.2 E)`, positive and smooth in `E` up to
  `sum p^(-1/2) ~ 1.2`, with no sign of a collapse to lattice points. A proof would also have to explain why the
  equal-weight fakes of `lw.md` 6.1 exist just above the window (`P ~ M^2 log M`) but not below it.
* *A cap of the window* needs a construction of gap-admissible spread angle models with primes `<= M^(2-delta)`. The
  union bound is dead (Proposition 6.1); what would be needed is a product-of-conditional-probabilities argument on
  `Real(Delta)` (the sequential structure that `smc_gaps.py` exploits numerically), i.e. a correlation inequality for
  the events `Bad_v` restricted to a rational subtorus. Heuristic 6.3 and §7.3 predict that such models exist.
* *Residue inputs* (P3-all, overview §6.3) are independent of the window and remain open there as everywhere.

### 8.4 Assessment

* A restricted growth theorem in the window cannot come from the large sieve: (F1) and (F2) are proved, so the large
  sieve is equivalent to spreading and spreading is compatible with the short-arc condition; a formal cap of the
  window additionally needs window fakes (Corollary 5.4). Probability that
  some other analytic input changes this inside this program: `< 2%`.
* Restricted growth below the window (`P_+(N') = o(log R)`, Corollary 2.2): **proved**, elementary, and sharp for the
  counting method (Remark 2.4).
* Most useful next step, if the window is pursued: try to **construct** window fakes by a sequential (conditional
  first-moment) argument on the realisation torus, starting with the Paley double-core profile on primes in
  `[M^(3/2), 2M^(3/2)]` grouped into columns of nearly equal total weight. Success would move the window from §5.5 to
  the capped table (with Corollary 5.4); a clean failure, with a reason, would be the first sign that the window is
  genuinely different. Estimated chance of a rigorous construction in one round: about 20-30%.

## 9. Constraint audit

The four constraints of the project owner, item by item.

* *Uniform exceptional sets from fixed-`S` Subspace Theorems.* Not used anywhere. No Subspace-type input appears.
* *Gaussian gcd or denominator costs.* Theorem 2.1 works with the normalised `N'` and keeps the content: (I1), (I2)
  are stated for `W'`, and Corollary 2.2 passes from `P_+(N)` to `P_+(N') <= P_+(N)` only in the harmless direction.
  The exact check divides by the exact Gaussian gcd of each window. Proposition 6.2 keeps the norm
  `N'' = prod Norm(beta_j)^(e_j)` of the shadow cluster; it is not discarded.
* *A drifting common phase replaced by zero.* The common phase `theta_0` is kept as a free parameter throughout
  §4: Lemma 4.1(a) shows that `h in H` moves it to `theta_0 + <2a_1 - e, h>`, and nothing sets it to zero. The
  realisation set is defined with `theta_0` quantified existentially.
* *A template-dependent threshold treated as uniform in the radius.* Theorem 2.1 and Corollaries 2.2-2.3 have
  thresholds depending only on `C` (and `eps`). Propositions 3.3, 4.2, 5.2, 6.1 are finite statements with explicit
  constants and no thresholds. Remark 4.4 has a profile-dependent threshold `c(B)`; it is flagged there, used only in
  a remark, and not claimed uniform. Corollary 5.4 is conditional and takes its parameters (`eta, delta, c, gamma,
  kappa`) fixed in advance, as in the round-9 R1 caveat.

## 10. Files

All in `prove_checks/` (run from that directory; outputs next to the scripts).

| file | what | type |
| --- | --- | --- |
| `mkangles.py`, `primes400.txt` | split primes `<= 400`, `pi_p = a + bi`, `theta_p` as `u64` (`2^64 = pi/2`) from 50-digit `atan2` | exact |
| `smoothscan.c`, `smoothscan2.c` | exhaustive scan of `B`-smooth split-only `n <= X`, best `C_k` per (decade, `P_+`) | exact fixed point |
| `scan_B400_1e10.txt` | validation run against `growth/referee_families_checks/brute_1e10.txt` | evidence |
| `scan2_B200_1e18.txt`, `scan2_B113_1e18.txt`, `sum2_B200_1e18.txt`, `summarize2.py` | the §7.1 runs and table | evidence |
| `scan_B200_1e16.txt`, `scan_B113_1e18.txt`, `sum_*.txt`, `summarize.py`, `count_smooth.py` | first (unbinned) runs and norm counts | evidence |
| `glib.py` | exact Gaussian-integer helpers (gcd, points of a norm, windows) | exact |
| `check_thm21.py`, `check_thm21_out.txt` | Theorem 2.1 on the record clusters (exact gcds) | exact |
| `cok_structure.py`, `cok_structure_out.txt` | `rank B`, Smith invariants, orders `D_j` (§7.2) | exact |
| `poisson_gaps.py`, `poisson_gaps_out.txt` | direct Monte Carlo of the gap-admissible measure vs Poisson | evidence |
| `smc_gaps.py`, `run_smc.sh`, `smc_out.txt` | sequential Monte Carlo of the same, small probabilities | evidence |
| `check_ls.py`, `check_ls_out.txt` | Lemma 3.2 / Prop 3.3 on genuine primes; Prop 4.2 demo | evidence |
| `../notes.md` | working notes | — |
