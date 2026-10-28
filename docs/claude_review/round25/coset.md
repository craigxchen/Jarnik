# Round 25: Proposition A' on realisation cosets (fakes open point 2 / round-21 judge §7)

Scripts and outputs: `coset_checks/`. Labels as in the overview: **proved** / **exact** / **evidence** /
**heuristic**. Nothing here has been refereed. Notation follows `round21/fakes.md` §1.1 unless stated.

## 0. Summary and verdict

**Verdict: partial.** No growth improvement, and no uniform bound. Scenario P is **not** proved on realisation cosets
for `m >= cM` in the window, and scenario R is not proved either. The target was "Proposition A' on realisation
cosets, as `Lambda -> infinity`". It is carried out (Theorem 4.1), but it is the wrong tool, and the round replaces it
with sharper statements.

1. **Theorem 2.2 (proved; full torus).** For distinct split primes, `mu(G(p)) > 0` **iff** `5 ∉ p`, for every `Lambda`.
   The mechanism is slack, not sequential integration. Lemma 2.1: a point that is gap-admissible for slightly smaller
   weights is a Lebesgue density point of `G(p)`. Transplant vectors (`arg beta_j`, `N(beta_j)` = the preceding split
   prime) are such points. This answers `fakes.md` open point 4 (no prime set without 5 is rigid). It supersedes
   Proposition A' (`Lambda <~ (1/4) log P`) on the full torus. It also gives `mu(G(p)) >= exp(-C_tau r^2 log P)`
   (Corollary 2.3), which refutes the Poisson extrapolation `exp(-c(Pi - 1))` at large `Lambda`, and with it `fakes.md`
   open point 3 for the coset `T^r` (Remark 2.3(2)).
2. **Theorem 3.1 (proved; sandwich on every coset).** If `G(p') ∩ coset != ∅` for some `p' < p`, then
   `nu(G(p) ∩ coset) > 0`, which in turn implies `G(p) ∩ coset != ∅`. The Haar-positivity form of (Q) is therefore the
   existence form up to an arbitrarily small loss in the weights. Theorem A's smallness, `exp(-c M^(delta/2))`, never
   decides positivity.
3. **Theorem 4.1 / Corollary 4.2 (proved): Proposition A' on cosets**, with exact step sums in graph coordinates.
   Dependent primes add mass at a step with no integration of their own.
4. **Theorem 5.1 / Corollary 5.2 (proved): barrier.** In any basis and any order, the step masses sum to at least
   `(kappa/2)(Pi_J - 1)`. So every step-by-step union-bound certificate needs `Pi_J < 1 + 2m/kappa`, which fails in the
   window when `m >= cM`, for `M >= M_0(c, delta)` (explicit, e.g. `4.6e5` at `(c, delta) = (1, 1/2)`). The judge's
   proposed route **cannot** reach the window.
5. **Proposition 6.1 (proved): circularity of arithmetic seeds.** A slack arithmetic point on a realisation coset is an
   actual window cluster with `log N_beta < W` (and `>= 3M log M - O(M)`). This is J2 for the measure form of scenario P.
6. **§7: what decides (Q) now.** Realisation cosets are `Delta`-translates of base cosets. Proposition 7.1 and
   Corollary 7.3 (proved) show that a base coset carrying a non-character slack seed, together with valid characters,
   gives positive Haar measure on the realisation coset whenever the local measure near the seed exceeds the
   heavy-window mass `e^{-W/4 + O(r)}`. This proves scenario P conditionally for **small** `m` (about
   `m <= M^{(2-delta)/5}`, sketch). It fails by a power of `M` for `m >= cM`. Lemma 7.2 (proved): Gaussian seeds force
   every character to be heavy. The `A_0` window question becomes (Q_seed) (§7.3): a profile-design part and a
   local-measure part.
7. **Evidence (§8).** Lines far from seeds are covered: exact emptiness for 61 lines in `T^8`, including one whose
   characters hold. Lines through or near seeds survive (heuristic A versus B). In the window the Haar-generic part of
   `G` is predicted empty. What survives is structured.

**Status for the map.** The low-height window stays **open** for `A_0`. The measure question is closed:
positivity reduces to existence (Theorem 3.1), and union bounds provably cannot help (Corollary 5.2). The full-torus
picture is corrected (Theorem 2.2, Remark 2.3). The decisive item is (Q_seed). Neither of its answers gives growth (J2).

## 1. Setting

**1.1 Angle models, gap sets, weights.** `T = R/(pi/2)Z` with normalised Haar measure `mu_1`; `mu` (or `mu_{T^r}`) is
Haar probability measure on `T^r`. A profile is `(p, e, (a_x)_{x <= M})` with distinct split primes `p_1, ..., p_r`,
`u_x = 2a_x - e`, `lambda_x = u_x - u_1`. For **real** weights `p'_1, ..., p'_r >= 5` and `v in Z^r` put
`w'(v) = sum_j |v_j| log p'_j` and `y'_v = e^{-w'(v)/2}`; for the genuine primes we write `w(v)`, `y_v`. The
**window** of `v` is the open set

```text
W_v = { t in T : dist(t, (pi/2)Z) < arcsin y_v }  ∪  { t in T : dist(t, pi/4 + (pi/2)Z) < arcsin(y_v/sqrt 2) },
omega(v) := mu_1(W_v) = (4/pi)(arcsin y_v + arcsin(y_v/sqrt 2))            (y_v <= 5^(-1/2): the two parts are disjoint),
```

and `G(p') ⊂ T^r` is the closed set of `phi` with `<v, phi> ∉ W'_v` for every `v != 0` (windows built from `y'_v`).
For the primes `p` this is the gap-admissible set `G(p)` of `round21/fakes.md` §1.1 (overview Lemma 5.1). Two
elementary facts are used throughout:

* *Monotonicity.* If `5 <= p'_j <= p_j` for all `j`, then `y'_v >= y_v`, so `W'_v ⊃ W_v` and `G(p') ⊂ G(p)`.
* *Window size.* `kappa y_v <= omega(v) <= (2 + sqrt 2) y_v`, `kappa = (4/pi)(1 + 2^(-1/2)) = 2.1735`, because
  `y <= arcsin y <= (pi/2) y` on `[0, 1]`. For `y_v <= 13^(-1/2)` the upper constant improves to `kappa_1 = 2.203`
  (`round21/fakes.md` Prop. A').

**1.2 Realisation cosets.** `L = span_Z{lambda_x}`, `L^sat = L_R ∩ Z^r` (its saturation), `l = rank L`, `m = r - l`
(free dimension), `D = [L^sat : L]`. `H = {phi : <lambda, phi> = 0 (mod pi/2) for all lambda in L}`, `H_0` its
identity component. For fixed `delta in [0, Delta]^M`, `Real(delta)` is empty or a coset of `H`, hence a finite union
of cosets of `H_0` (`fakes.md` §1.1). Let `K = {b in Z^r : <b, lambda> = 0 for all lambda in L}` (a saturated lattice
of rank `m`), let `b_1, ..., b_m` be a `Z`-basis of `K` and `B = [b_1 | ... | b_m] in Z^(r x m)`.

> **Lemma 1.1 (coordinates on a coset).** `theta -> B theta` is an isomorphism of topological groups `T^m -> H_0`.
> Hence for every `g in T^r`, `theta -> g + B theta` carries Haar measure of `T^m` to the Haar probability measure
> `nu` of the coset `g + H_0`. For `v in Z^r` put `mu_v := B^T v in Z^m`. Then `<v, g + B theta> = <v, g> + <mu_v, theta>`,
> and `mu_v = 0` iff `v in L^sat`.

*Proof.* `K` is saturated, so `b_1, ..., b_m` extend to a basis of `Z^r` and there is `A in Z^(m x r)` with `AB = I_m`.
If `B theta = 0` in `T^r` then `theta = AB theta = 0` in `T^m`: the map is injective. Its image is a compact connected
subgroup of dimension `m`, contained in `H` because `<lambda, B theta> = <B^T lambda, theta> = 0` for `lambda in L`. As
`dim H = r - l = m` and `H_0` is the identity component, the image is `H_0`; a continuous bijective homomorphism of
compact groups is an isomorphism and carries Haar measure to Haar measure. Finally `B^T v = 0` iff `v` is orthogonal
to `K`, iff `v in K^perp ∩ Z^r = L_R ∩ Z^r`. QED

Monomials with `mu_v = 0` are the **characters** of the coset: `<v, phi>` is constant (`= <v, g>`) on `g + H_0`, and
their gap conditions are a property of the coset alone (the "character margins" of `fakes.md`; for `v in L` the value
is `sum_x c_x (delta_x - delta_1)`, for `v in L^sat \ L` it is a `D`-torsion value up to `O(|v|_1 Delta)`). All other
monomials are **non-characters**. The **full torus** is the case `L = 0`: `B = I_r`, `m = r`, no characters.

**1.3 Arithmetic (transplant) vectors.** For Gaussian integers `beta_j != 0` write `arg beta = (arg beta_j mod pi/2)_j
in T^r`, and for `v in Z^r` put `gamma_v = prod_{v_j > 0} beta_j^{v_j} prod_{v_j < 0} conj(beta_j)^{|v_j|}`, so that
`<v, arg beta> = arg gamma_v (mod pi/2)` and `|gamma_v|^2 = prod_j N(beta_j)^{|v_j|}`.

> **Lemma 1.2 (Gaussian-prime seeds).** Let `q_1, ..., q_r` be distinct split primes and `beta_j` a Gaussian prime of
> norm `q_j`. Then `arg beta in G(q)` (weights `q_j`).

*Proof.* (This is `fakes.md` Lemma T in the special case of prime norms; the short argument is repeated because §2
depends on it.) Let `v != 0` and `gamma_v = X + iY`. Then `dist(arg gamma_v, (pi/2)Z) = arcsin(min(|X|,|Y|)/|gamma_v|)` and
`dist(arg gamma_v, pi/4 + (pi/2)Z) = arcsin(min(|X - Y|, |X + Y|)/(sqrt 2 |gamma_v|))`. `X^2 + Y^2 = prod q_j^{|v_j|}` is odd,
so `X`, `Y` have different parity and `X ± Y` are odd, hence non-zero. If `Y = 0` (or `X = 0`), `conj(gamma_v) = ±gamma_v`.
Take `j` with `v_j != 0`, say `v_j > 0` (the other case is symmetric). Then `beta_j | gamma_v`, so
`conj(beta_j) | conj(gamma_v) = ±gamma_v`.
But the Gaussian primes dividing `gamma_v` are the `beta_i` (`v_i > 0`) and `conj(beta_i)` (`v_i < 0`), of pairwise
distinct norms except that `beta_j`, `conj(beta_j)` share the norm `q_j`; as `q_j` is an odd split prime, `conj(beta_j)` is
not associate to `beta_j`, and only `beta_j` occurs. Contradiction. So `min(|X|, |Y|) >= 1` and `min |X ± Y| >= 1`, i.e.
both distances are at least `arcsin(|gamma_v|^{-1})`, resp. `arcsin(|gamma_v|^{-1}/sqrt 2)`, with
`|gamma_v|^{-1} = e^{-w_q(v)/2}`. QED

## 2. The full torus: positivity for every prime set without 5

The mechanism is not sequential integration but **slack**: a point that is gap-admissible for smaller weights is a
Lebesgue density point of `G(p)`. The lemma is stated on an arbitrary coset because §3 uses it there; the full torus is
`B = I_r`.

> **Lemma 2.1 (slack points are density points).** Let `g + H_0` be a coset with coordinates `theta -> g + B theta`
> (Lemma 1.1), `b* = max_i sum_k |B_ik|` (so `|mu_v|_1 <= b* |v|_1`). Let `5 <= p'_j < p_j` for all `j`, put
> `tau = min_j (p_j/p'_j)^(1/2) > 1`, and let `psi = g + B theta_0 in G(p')`. For `0 < eps < pi/4` let `Q_eps` be the
> image of `theta_0 + [-eps, eps]^m`. Then
>
> ```text
> nu(G(p) ∩ Q_eps) >= (1 - F(eps)) nu(Q_eps),      nu(Q_eps) = (4 eps/pi)^m,
> F(eps) = c_2 m sum_{v in U(eps)/±} (1 + b* |v|_1)/(tau^{|v|_1} - 1),    c_2 = pi (1 + sqrt 2) = 7.585...,
> U(eps) = { v not in L^sat : |mu_v|_1 eps >= sigma_v },                sigma_v = (tau^{|v|_1} - 1) y_v / sqrt 2,
> ```
>
> and `F(eps) -> 0` as `eps -> 0`. In particular `nu(G(p) ∩ (g + H_0)) > 0`.

*Proof.* (i) *Slack.* Let `v != 0`. Since `psi in G(p')`, `dist(<v,psi>, (pi/2)Z) >= arcsin y'_v`, and
`y'_v = y_v prod_j (p_j/p'_j)^{|v_j|/2} >= tau^{|v|_1} y_v`, with `y'_v <= 5^(-1/2)`. As `arcsin b - arcsin a >= b - a` for
`0 <= a <= b <= 1`, the distance exceeds the window radius `arcsin y_v` by at least `(tau^{|v|_1} - 1) y_v >= sigma_v`. In the
same way `dist(<v,psi>, pi/4 + (pi/2)Z)` exceeds `arcsin(y_v/sqrt 2)` by at least `sigma_v`.

(ii) *Which monomials can fail in `Q_eps`.* For `theta = theta_0 + t`, `|t|_inf <= eps`, Lemma 1.1 gives
`<v, g + B theta> = <v, psi> + <mu_v, t>` with `|<mu_v, t>| <= |mu_v|_1 eps`. Both distance functions are 1-Lipschitz. If
`v in L^sat` then `mu_v = 0` and the condition holds at every point of `Q_eps` (it holds at `psi`, with `G(p') ⊂ G(p)`).
If `v not in L^sat` and `|mu_v|_1 eps < sigma_v`, it holds at every point of `Q_eps` by (i). So every point of
`Q_eps \ G(p)` lies in `E_v = {theta in Q_eps : <v, g + B theta> in W_v}` for some `v in U(eps)`; and `E_v = E_{-v}`.

(iii) *Size of `E_v` for `v in U(eps)`.* Under the uniform distribution of `t` on `[-eps, eps]^m`, `s = <mu_v, t>` is a sum
of independent terms one of which, `mu_{v,k} t_k` with `|mu_{v,k}| = |mu_v|_inf >= 1`, has density `<= 1/(2 eps |mu_v|_inf)`;
hence `s` has density `<= 1/(2 eps |mu_v|_inf)`. The values of `s` lie in an interval of length `2|mu_v|_1 eps`, which meets at
most `4|mu_v|_1 eps/pi + 2` translates (by `(pi/2)Z`) of each of the two components of the window, of lengths
`2 arcsin y_v <= pi y_v` and `2 arcsin(y_v/sqrt 2) <= pi y_v/sqrt 2`. Hence

```text
nu(E_v)/nu(Q_eps) <= pi (1 + 2^(-1/2)) y_v (4 |mu_v|_1 eps/pi + 2) / (2 eps |mu_v|_inf).
```

For `v in U(eps)`, `y_v <= sqrt 2 |mu_v|_1 eps/(tau^{|v|_1} - 1)`; inserting this, and using `|mu_v|_1 <= m |mu_v|_inf`,
`2|mu_v|_1 eps/pi <= b*|v|_1/2`, gives `nu(E_v)/nu(Q_eps) <= c_2 m (1 + b*|v|_1)/(tau^{|v|_1} - 1)`.

(iv) *Summation.* By (ii)-(iii), `nu(Q_eps \ G(p)) <= F(eps) nu(Q_eps)`. The series
`sum_{v != 0} (1 + b*|v|_1)/(tau^{|v|_1} - 1)` converges: `#{v : |v|_1 = n} <= 2^r C(n + r - 1, r - 1)` is polynomial in `n`
and `tau > 1`. Each `v not in L^sat` leaves `U(eps)` once `eps < sigma_v/|mu_v|_1`. By dominated convergence `F(eps) -> 0`. QED

> **Theorem 2.2 (positivity on the full torus; dichotomy).** Let `p_1, ..., p_r` be distinct split primes. Then
> `mu_{T^r}(G(p)) > 0` **if and only if** `5 ∉ {p_1, ..., p_r}`. There is no condition on `Lambda = sum_j p_j^(-1/2)`.

*Proof.* If `p_j = 5`, then `G(p) ⊂ {phi : phi_j in G_1(5)}` and `G_1(5)` has two points (`fakes.md` Proposition B,
refereed), so `mu(G(p)) = 0`. Conversely, let `s_1 = 5 < s_2 = 13 < s_3 = 17 < ...` be the split primes. If 5 is not
among the `p_j`, then `p_j = s_{k_j}` with `k_j >= 2`; put `q_j = s_{k_j - 1}`. The `q_j` are distinct split primes with
`5 <= q_j < p_j`. Let `beta_j` be a Gaussian prime of norm `q_j`. By Lemma 1.2, `arg beta in G(q)`. Lemma 2.1 on the full
torus (`L = 0`, `B = I`, `m = r`, `b* = 1`) with `p' = q` gives `mu(G(p)) > 0`. QED

> **Corollary 2.3 (quantitative form).** Suppose there are distinct split primes `q_j <= p_j/tau^2` (`tau > 1`), and let
> `P = max_j p_j`, `z = tau^(-1/2)`. Let `n_0` be the least integer with
> `(c_2 r/2)(tau/(tau - 1)) tau^(-n_0/2) ((1 + z)/(1 - z))^r (1 + 2rz/(1 - z^2)) <= 1/2`, and
> `eps_* = min(pi/5, (tau - 1) P^{-(n_0 - 1)/2}/(sqrt 2 n_0))`. Then `mu(G(p)) >= (1/2)(4 eps_*/pi)^r`. For fixed `tau`,
> `n_0 = (2/log tau) log((1 + z)/(1 - z)) r + O_tau(log r)`, so
> `log mu(G(p)) >= -(1 + o(1)) (log((1+z)/(1-z))/log tau) r^2 log P` (constant `2.54` at `tau = 2`).

*Proof.* In Lemma 2.1 (full torus), `(1 + n)/(tau^n - 1) <= (tau/(tau - 1))(1 + n) tau^{-n}` for `n >= 1`. If every `v` with
`|v|_1 < n_0` is safe, then by Rankin's trick with `z^2 = 1/tau`,
`F(eps) <= (c_2 r/2)(tau/(tau-1)) tau^{-n_0/2} sum_v (1 + |v|_1) z^{|v|_1}`, and
`sum_v (1 + |v|_1) z^{|v|_1} = f(z) + z f'(z)` with `f(z) = ((1+z)/(1-z))^r`, `f'(z) = f(z) 2r/(1 - z^2)`. So `F(eps) <= 1/2`.
Safety of `|v|_1 <= n_0 - 1` at `eps = eps_*`: `|v|_1 eps_* < (tau - 1) P^{-(n_0-1)/2}/sqrt 2 <= sigma_v`, because
`tau^{|v|_1} - 1 >= tau - 1` and `y_v >= P^{-|v|_1/2}`. The asymptotics of `n_0` are read off the defining inequality. QED

**Remarks.** (1) *What this corrects.* `fakes.md` open point 4 asked whether some set of split primes `>= 13` has
`mu(G(p)) = 0`; Theorem 2.2 says never. Proposition A' (`Lambda <~ (1/4) log P`) is superseded on the full torus: the
limit there was the union bound, not the set. Theorem A is untouched (it is an upper bound); positivity and Haar
scarcity coexist. (2) *The Poisson extrapolation is false at large `Lambda`.* `fakes.md` §5.1 extrapolated
`mu(G(p)) ≈ exp(-c(Pi - 1))`, `c ≈ 1.1`, and inferred `mu(G(p)) <= exp(-c exp(M^(delta/2)/4))` in the window "if the
law persists". Corollary 2.3 gives `mu(G(p)) >= exp(-C_tau r^2 log P)`. Take window prime sets of the kind fakes use:
`r >= M/4` primes in `[P/2, P]`, `P = M^(2-delta)`, with `r` at most the number of split primes in `[5, P/8]` (so
`tau = 2` seeds exist; true for `M/4 <= r <= M` once `M >= M_0(delta)`, by the prime number theorem in progressions).
Then `Pi >= exp(Lambda) >= exp(r P^(-1/2)) >= exp(M^(delta/2)/4)` (`fakes.md` (1.1)), while the exponent of
Corollary 2.3 is `(2.54 + o(1)) r^2 log P <= (5.1 + o(1)) M^2 log M`. So for `M >= M_1(delta)` the extrapolated value
`exp(-1.1(Pi - 1))` is smaller than a proved lower bound. Fakes' open point 3 (a doubly
exponential upper bound `nu(G ∩ coset) <= exp(-c(Pi - 1))`) is therefore **false in that generality**: it fails for the
coset `T^r` itself (`L = 0`). Whether it holds on realisation cosets of barrier-scale profiles, which need not pass
near arithmetic points (§6), is open. The numbers are astronomically
large: the crossover needs `Lambda >~ log r` (`coset_checks/quant_bound.py`, §8.3). (3) *Where the measure sits.*
Lemma 2.1 puts positive measure in neighbourhoods of arithmetic points. §8.1 shows (evidence) that for 8 primes
`13..73` the survival rate in cubes around the transplant vector rises to 1, from 0.07 at `eps = 1e-4` to 1.0 at
`eps = 1e-7`, while uniform points survive at rate 0 (`T = 24`). Whether the positive measure at large `Lambda` lies
only near arithmetic points is the heuristic question of §7.

## 3. Cosets: Haar positivity is existence with slack (the sandwich)

> **Theorem 3.1 (sandwich).** Let `P` be any profile, `delta in [0, Delta]^M`, and `C` a coset of `H_0` contained in
> `Real(delta)` (or any coset of `H_0`). For every choice of weights `5 <= p'_j < p_j`:
>
> ```text
> G(p') ∩ C != ∅   ==>   nu(G(p) ∩ C) > 0   ==>   G(p) ∩ C != ∅ .
> ```
>
> The same holds with `C` replaced by `Real(delta)` (a finite union of such cosets) and `nu` by its Haar measure.

*Proof.* The first implication is Lemma 2.1 applied at a point `psi in G(p') ∩ C` (characters included: they hold at
`psi` for the weights `p' <= p`, and are constant on `C`). The second is trivial. QED

**Consequences for (Q).** The judge's target asks for `nu(G ∩ coset) > 0`, scenario P. By Theorem 3.1 this is the
existence question for a gap-admissible point of `Real(delta)`, up to an arbitrarily small reduction of the weights:

* any construction of **one** angle model with slack, `phi in G(p') ∩ Real(delta)` with `p' < p`, gives positive Haar
  measure, so measure positivity is never the obstruction;
* conversely, a proof that `nu(G(p) ∩ Real(delta)) = 0` would have to exclude slack points, i.e. it is (essentially)
  an emptiness proof of scenario-R type.

Theorem A (`fakes.md`: `nu <= exp(-2.17 sum_J p_j^(-1/2))`) is compatible with all of this. In the window, a positive
measure is at most `exp(-c M^(delta/2))`, but its positivity is decided by existence alone. For capping only
existence matters, so "Haar positivity" was never a separate requirement: Theorem 3.1 makes this precise and
shows that it cannot be the hard part. What remains is the existence part of (Q). §6 shows that the slack points
supplied by arithmetic cannot be used on realisation cosets without presupposing actual clusters.

## 4. Proposition A' on cosets: sequential completion

This is the extension the judge asked for. It is proved in full. By §5, however, it cannot reach the window, and by
§2–3 it is not needed for positivity wherever a slack point is available.

Fix a coset `g + H_0`, a basis `B` of `K` and the coordinate order `theta_1, ..., theta_m`. For `k <= m` let
`F_k = {mu in Z^m : mu_k > 0, mu_{k'} = 0 for k' > k}`; every `mu != 0` lies in exactly one of `±F_1, ..., ±F_m`. For
`mu in Z^m \ {0}` define the **fibre window**

```text
B_mu = ∪_{v in Z^r : mu_v = mu} ( W_v - <v, g> )  ⊂ T,          U(mu) = mu_1(B_mu) <= min(1, sum_{v : mu_v = mu} omega(v)).
```

(`W_v` is symmetric, so the conditions for all `v` with `mu_v = ±mu` say exactly `<mu, theta> ∉ B_mu`.)

> **Theorem 4.1 (sequential completion on a coset).** Assume the character condition
> `(Ch)`: `<v, g> ∉ W_v` for every `v in L^sat \ {0}`. Put `U_k = sum_{mu in F_k} U(mu)`. If `U_k < 1` for every `k`, then
> `nu(G(p) ∩ (g + H_0)) >= prod_{k=1}^m (1 - U_k) > 0`.

*Proof.* By Lemma 1.1, `G(p) ∩ (g + H_0)` corresponds to the set of `theta in T^m` such that `(Ch)` holds and
`<mu, theta> ∉ B_mu` for all `mu in ∪_k F_k`. The conditions with `mu in F_k` involve only `theta_1, ..., theta_k`, and
`theta_k` enters with coefficient `mu_k >= 1`. For fixed `theta_{<k}` let `A_k(theta_{<k})` be the set of `theta_k`
satisfying all of them. Its complement is `∪_{mu in F_k} {t : mu_k t + <mu_{<k}, theta_{<k}> in B_mu}`. Since
`t -> mu_k t + c` maps Haar measure to Haar measure, each set has measure `U(mu)`, so `mu_1(A_k) >= 1 - U_k` for
every `theta_{<k}`. Integrate over `theta_m`, then `theta_{m-1}`, ..., then `theta_1` (Fubini; all sets are closed). QED

For `B = I` (the full torus, `L = 0`) and `F_k` ordered by increasing primes, `U_k <= kappa_1 sum_{v_k >= 1, v_{>k} = 0} y_v =
kappa_1 Pi_{<k}/(sqrt(p_k) - 1)`: this is Proposition A'. On a general coset the sums are explicit in **graph
coordinates**. Suppose the coset is a graph over a set `J = {j_1, ..., j_m}` of free primes: there are `C in Z^((r-m) x m)`
(rows `c_i`, `i ∉ J`) and constants `g_i` with `phi_i = <c_i, phi_J> + g_i` on `g + H_0`. (This is the case
`det B_J = ±1`. In general `pi_J : H_0 -> T^J` is an isogeny of degree `|det B_J|`, and Theorem 4.1 applies in any basis.)
With `theta = phi_J`, `mu_v = v_J + f(v_{J^c})`, `f(x) = C^T x`.

> **Corollary 4.2 (graph coordinates; exact step sums).** For `x in Z^(J^c)` put `y(x) = prod_{i ∉ J} p_i^{-|x_i|/2}` and
> `E_{>k}(x) = prod_{k' > k} p_{j_k'}^{-|f(x)_{k'}|/2}`. Let `Pi_j = 1 + 2/(sqrt(p_j) - 1)` and
> `Pi_{J,<k} = prod_{k' < k} Pi_{j_k'}`. Then
>
> ```text
> sum_{v : last(mu_v) = k} y_v = Pi_{J,<k} sum_{x in Z^(J^c)} y(x) E_{>k}(x) ( Pi_{j_k} - p_{j_k}^{-|f(x)_k|/2} ),
> ```
>
> and `U_k <= (c/2) x (this sum)`, with `c = kappa_1` if all primes are `>= 13` (`c = 2 + sqrt 2` in general). So
> `nu(G(p) ∩ (g + H_0)) > 0` whenever `(Ch)` holds and `(c/2)` times the right-hand side is `< 1` for every `k`.

*Proof.* Write `v = (v_J, x)`. Then `last(mu_v) = k` iff `v_{j_k'} = -f(x)_{k'}` for `k' > k` and `v_{j_k} != -f(x)_k`. For such
`v`, `y_v = y(x) E_{>k}(x) prod_{k' < k} p_{j_k'}^{-|v_{j_k'}|/2} p_{j_k}^{-|v_{j_k}|/2}`. Sum over free `v_{j_k'} in Z` (`k' < k`):
`Pi_{J,<k}`. Sum over `v_{j_k} in Z \ {-f(x)_k}`: `Pi_{j_k} - p_{j_k}^{-|f(x)_k|/2}`. The set `{v : last(mu_v) = k}` is symmetric
under `v -> -v`, and `F_k` takes one of each pair, which gives the factor `1/2`. Finally `U(mu) <= sum_{mu_v = mu} omega(v)`. QED

*Exact check.* `coset_checks/graph_formula_check.py`: for 5 primes `13..41`, `m = 3` and an integer `C`, brute-force
summation over a box agrees with the formula to relative error `5.6e-8`, the size of the box truncation.

**Reading.** The term `x = 0` is Proposition A' (`Pi_{j_k} - 1 = 2/(sqrt(p_{j_k}) - 1)`). A dependent prime `i` whose
row `c_i` has last non-zero entry at position `k` contributes `x = ±e_i` with weight `p_i^{-1/2}(Pi_{j_k} - p_{j_k}^{-|c_{ik}|/2})`.
This contribution carries **no** factor `p_{j_k}^{-1/2}`: dependent primes act as extra mass at step `k` with no
integration of their own. If all primes lie in `[P/2, P]` and at most `h` dependent primes have their last entry at each
step (and the rows `c_i` are bounded), the conditions of Corollary 4.2 hold for
`Lambda <= (1/4) log P - (1/2) log(1 + h) - O(1)` (*sketch*, the same computation as for Proposition A'; not used below). Hence Proposition A' extends to cosets with `Lambda -> infinity`, but only up to
`Lambda ~ (1/4) log P`. In the window `Lambda >= c M^(delta/2) >> log P`, and §5 shows that no basis and no ordering
helps there.

## 5. The barrier: no sequential union-bound certificate reaches the window

> **Theorem 5.1 (total step mass).** For every coset `g + H_0`, every basis `B` of `K` and every coordinate order, and
> every `J ⊂ {1..r}` whose unit vectors are linearly independent modulo `L_R`,
>
> ```text
> sum_{k=1}^m U_k  >=  (kappa/2) (Pi_J - 1),        Pi_J = prod_{j in J} (1 + 2/(sqrt(p_j) - 1)) >= exp( sum_{j in J} p_j^(-1/2) ).
> ```
>
> Hence the hypothesis of Theorem 4.1 (`U_k < 1` for all `k`) implies `Pi_J < 1 + 2m/kappa`. The same holds for any
> upper bounds `U'_k >= U_k` used in its place, in particular for Corollary 4.2 and for Proposition A'.

*Proof.* Let `v in Z^J \ {0}` (support in `J`). Then `v ∉ L_R`, so `mu_v != 0`. If `mu_v = mu_{v'}` for `v, v'` supported
on `J`, then `v - v' in L_R ∩ span(e_J) = 0`. So `v -> ±mu_v` is injective on `(Z^J \ {0})/±`. Replace `v` by `-v` if
necessary so that `mu_v in ∪_k F_k`. Then `U(mu_v) >= mu_1(W_v - <v, g>) = omega(v) >= kappa y_v`. Summing over distinct
`mu`: `sum_k U_k >= kappa sum_{(Z^J \ 0)/±} y_v = (kappa/2)(Pi_J - 1)`, since `sum_{v in Z^J} y_v = Pi_J`. Last,
`log(1 + 2/(s - 1)) >= 1/s` for `s = sqrt p >= sqrt 5` (`fakes.md` (1.1)). QED

> **Corollary 5.2 (the window is out of reach of every such certificate).** Let `0 < delta < 1`, `c > 0`, and consider
> any profile with all primes `<= M^(2-delta)` and free dimension `m >= cM`. If
> `exp(c M^(delta/2)) >= 1 + 2M^2/kappa`, then no basis and no ordering satisfies the hypothesis of Theorem 4.1, on any
> coset of `H_0`. This holds for all `M >= M_0(c, delta)`; for example
> `M_0 ≈ 4.5e3` for `(c, delta) = (1/4, 1)`, `4.6e5` for `(1, 1/2)`, `7.0e8` for `(1/4, 1/2)`, `4.9e17` for `(1/2, 1/4)`
> (`coset_checks/barrier_M0.py`).

*Proof.* Take `J` with `|J| = m` (`fakes.md` Theorem A(e)). Then `sum_{j in J} p_j^(-1/2) >= m M^(-1 + delta/2) >= c M^(delta/2)`,
so `Pi_J >= exp(c M^(delta/2))`. Also `m <= r < M^2` (there are fewer than `M^(2-delta)` primes below `M^(2-delta)`). If Theorem
4.1 applied, `Pi_J < 1 + 2m/kappa < 1 + 2M^2/kappa`, a contradiction. The function
`c M^(delta/2) - log(1 + 2M^2/kappa)` is increasing beyond its last zero, which gives the threshold. QED

**What the barrier says and does not say.** It is a statement about a *method*: positivity certified step by step by a
union bound, in any coordinates. Such certificates need `Pi_J = O(m)`, that is `Lambda_J <~ log m`. On the full torus
with primes in `[P/2, P]` this is `Lambda <~ (1/4) log P + O(log log P)`, so Proposition A' is optimal for its method up
to lower-order terms. The barrier does **not** say that positivity fails. On the full torus it holds for every `Lambda`
(Theorem 2.2), through slack points that no step-by-step union bound sees. On realisation cosets the question is open,
and §§6–7 explain why the full-torus mechanism does not transfer.

## 6. Arithmetic seeds on realisation cosets are actual clusters

Theorem 2.2 seeds Lemma 2.1 with arithmetic vectors. On a realisation coset an arithmetic seed is an actual cluster.
This is the referee's strengthening of `fakes.md` Lemma T (`referee_fakes.md` §4), combined with the `(2/3)`
theorem.

> **Proposition 6.1.** Let `(p, e, (a_x))` be a profile with `M >= 2` points, `W = sum_j e_j log p_j`,
> `Delta = C e^{-W/4}`, and let `beta_j in Z[i] \ {0}` with `N(beta_j) <= p_j`. If `psi = arg beta in G(p) ∩ Real(delta)` for
> some `delta in [0, Delta]^M`, then the Gaussian integers `z_x = prod_j beta_j^{a_xj} conj(beta_j)^{e_j - a_xj}` are, after
> multiplication by suitable units, `M` distinct lattice points of the circle `|z|^2 = N_beta := prod_j N(beta_j)^{e_j}`,
> lying on an arc of length `<= C_beta N_beta^{1/4}` with `C_beta = C (N_beta e^{-W})^{1/4} <= C`. Consequently
> `log N_beta >= 3M log M - (14 + 4 log^+ C) M`, and every prime factor of `N_beta` is `<= max_j p_j`.

*Proof.* `arg z_x = <u_x, (arg beta_j)_j>`, which is `<u_x, psi>` modulo `pi/2`. By realisation there are units `eps_x` with
`arg(eps_x z_x) in theta_0 + [0, Delta]` modulo `2 pi`. For `x != y` and `d = a_x - a_y != 0`,
`z_x/z_y = prod_j (beta_j/conj(beta_j))^{d_j} = gamma_{lambda}/prod_j N(beta_j)^{|d_j|}`, with `lambda = u_x - u_y = 2d`. If
`eps_x z_x = eps_y z_y`, then `arg gamma_lambda in (pi/2)Z`, i.e. `<lambda, psi> in (pi/2)Z`, contradicting `psi in G(p)`.
The `M` points lie on the circle of radius `R_beta = N_beta^{1/2}`, in an angle `Delta`, hence on an arc of length
`Delta N_beta^{1/2} = C (N_beta/e^W)^{1/4} N_beta^{1/4}`. Since `N_beta <= e^W`, `C_beta <= C`. Theorem 1.1 of the expert note
(`two_thirds/proof.md` Thm 7.5) gives `3M log M <= log N_beta + (14 + 4 log^+ C_beta) M`. Prime factors of `N_beta`
divide some `N(beta_j) <= p_j`. QED

**Consequence (the circularity, made precise).** Let `A > 0`. Suppose Lemma 2.1 is seeded with an arithmetic slack
point (`N(beta_j) < p_j`) on a realisation coset of a barrier-scale window profile (`W <= A M log M`, primes
`<= M^(2-delta)`). Then there is an **actual** window cluster of `M` points with `log N_beta < W <= A M log M` and arc
constant `<= C`. So the density mechanism yields window fakes from arithmetic seeds only when actual window clusters
at barrier scale already exist. If they exist for unbounded `M`, growth fails inside the window, and nothing is gained by
capping `A_0` there. This is the judge's J2, now for the **measure** form of scenario P as well. Scenario P therefore
needs a **non-arithmetic** slack point on `Real(delta)`. By Theorem 3.1 such a point is also necessary, up to an
arbitrarily small loss in the weights.

**Characters are what forbids `H_0` itself.** On the coset `H_0` through `0` (`delta = 0`, `g = 0`) every character takes
the value `0` and fails its gap condition, so `G(p) ∩ H_0 = ∅`. Yet `H_0` can carry points that satisfy every
**non-character** condition with slack (Lemma 7.2 with `g_0 = 0`). The pair margins
`|<lambda_x, phi>| >= arcsin y_{lambda_x}` force a realisation coset to be a translate of such a base coset by an `eta`
of size about `Delta`. §7 quantifies what survives this translation.

## 7. What decides (Q) now: seeded base cosets

A realisation coset is a `Delta`-translate of a component of `H`. As `delta -> 0`, `Real(delta)` tends to `H`, whose
components are the **base cosets** `C_0 = g_t + H_0` (`g_t` torsion). So the right object is not a generic coset but a
base coset and its small translates `C_eta = C_0 + eta`, `|eta|_inf <= rho`, with `rho` of order `Delta` times the norm of a
right inverse of `(lambda_x)` (a profile-dependent factor, which must not be treated as uniform).

### 7.1 Perturbation from a seeded base coset (proved)

> **Proposition 7.1.** Let `C_0 = g_0 + H_0`, and let `psi_0 in C_0` satisfy every **non-character** gap condition for
> weights `p'` (`5 <= p'_j < p_j`, `tau = min (p_j/p'_j)^(1/2)`); characters are not required at `psi_0`. Let `eta in T^r`,
> `C_eta = C_0 + eta`, and assume `(Ch_eta)`: `<v, g_0 + eta> ∉ W_v` for all `v in L^sat \ {0}`. Fix `1 < tau_1 < tau` and let
> `Heavy(eta) = {v ∉ L^sat : |<v, eta>| > (tau^{|v|_1} - tau_1^{|v|_1}) y_v/sqrt 2}`. Then, for `0 < eps < pi/4`,
>
> ```text
> nu(G(p) ∩ C_eta) >= (4 eps/pi)^m (1 - F_1(eps)) - sum_{v in Heavy(eta)/±} omega(v),
> ```
>
> where `F_1` is `F` of Lemma 2.1 with `tau_1` in place of `tau` and the sum restricted to `v ∉ Heavy(eta)`.

*Proof.* Put `psi = psi_0 + eta in C_eta`. For `v ∉ L^sat ∪ Heavy(eta)`, Lemma 2.1(i) at `psi_0` gives an excess of at least
`(tau^{|v|_1} - 1) y_v/sqrt 2` over both window radii. The shift `<v, eta>` uses at most `(tau^{|v|_1} - tau_1^{|v|_1}) y_v/sqrt 2`
of it, leaving `sigma_v(tau_1)`. Steps (ii)-(iv) of Lemma 2.1 used only this excess, so the points of `Q_eps` (around
`psi`) that fail such a `v` have relative measure at most `F_1(eps)`. Characters hold on all of `C_eta` by `(Ch_eta)`. For
`v in Heavy(eta)`, `mu_v != 0`, and the set of `theta in T^m` with `<v, g_0 + eta> + <mu_v, theta> in W_v` has Haar measure
exactly `omega(v)`. QED

> **Lemma 7.2 (Gaussian seeds on base cosets).** Let `B` be a basis of `K` and `q_1, ..., q_m` distinct split primes with
> `prod_k q_k^{|B_ik|} <= p_i/(2 tau^2)` for every `i <= r` (`tau > 1`). Let `beta_k` be a Gaussian prime of norm `q_k` and
> `g_0 in ((pi/4)Z)^r` (such points lie in `H`, because `L ⊂ 2Z^r`). Then `psi_0 = g_0 + B arg beta` satisfies every non-character gap condition for the weights
> `p'_i = p_i/tau^2`.
> **Conversely, characters on such a coset are forced to be heavy.** If `(Ch_eta)` holds with `|eta|_inf <= rho`, then
> every `v in L^sat \ {0}` has `arcsin(y_v/sqrt 2) <= |v|_1 rho`.

*Proof.* For `v ∉ L^sat`, `<v, psi_0> = <v, g_0> + arg gamma_{mu_v}` with `<v, g_0> in (pi/4)Z` and `gamma_mu` built from the
`beta_k` as in §1.3. By Lemma 1.2, `gamma_{mu_v} = X + iY` has `X, Y, X ± Y` non-zero. Adding an even multiple of `pi/4`
keeps the two distances; adding an odd one exchanges them. In both cases the distance to `(pi/2)Z`, and to
`pi/4 + (pi/2)Z`, is at least `arcsin(1/(sqrt 2 |gamma_{mu_v}|))`. Now
`log |gamma_{mu_v}|^2 = sum_k |(B^T v)_k| log q_k <= sum_i |v_i| sum_k |B_ik| log q_k <= w(v) - |v|_1 log(2 tau^2)`.
So `1/(sqrt 2 |gamma_{mu_v}|) >= tau^{|v|_1} y_v`, which is the gap for weights `p'`. For the converse: `<v, g_0> in (pi/4)Z`
is a window centre, so `dist(<v, g_0 + eta>, (pi/4)Z) <= |v|_1 rho`. `(Ch_eta)` then needs
`arcsin(y_v/sqrt 2) <= |v|_1 rho`. QED

> **Corollary 7.3 (explicit sufficient condition).** In Proposition 7.1 let `|eta|_inf <= rho`, `P = max_j p_j`,
> `z = tau_1^(-1/2)`, `f(z) = ((1+z)/(1-z))^r`, let `n_0` be the least integer with
> `(c_2 m/2)(tau_1/(tau_1 - 1)) tau_1^(-n_0/2) f(z)(1 + 2 b* r z/(1 - z^2)) <= 1/2`, and
> `eps_* = min(pi/5, (tau_1 - 1) P^{-(n_0-1)/2}/(sqrt 2 b* n_0))`. For any `s in (0, 1)`,
> `sum_{Heavy(eta)/±} omega(v) <= H_s(rho) := (2 + sqrt 2)(sqrt 2 rho/(tau - tau_1))^s sum_{v != 0} |v|_1^s y_v^{1-s}`.
> Hence `nu(G(p) ∩ C_eta) > 0` whenever `(1/2)(4 eps_*/pi)^m > H_s(rho)`.

*Proof.* As in Corollary 2.3, with the majorant `(1 + b*|v|_1)` (`sum_v (1 + b*|v|_1) z^{|v|_1} = f(z)(1 + 2b* r z/(1 - z^2))`) and
with safety of `|v|_1 < n_0` from `|mu_v|_1 eps_* <= b*|v|_1 eps_* < sigma_v(tau_1)`. For `v in Heavy(eta)`,
`y_v < sqrt 2 |v|_1 rho/(tau^{|v|_1} - tau_1^{|v|_1}) <= sqrt 2 |v|_1 rho/(tau - tau_1)`, so `y_v <= y_v^{1-s}(sqrt 2 |v|_1 rho/(tau - tau_1))^s`,
and `omega <= (2 + sqrt 2) y`. QED

**Reading (asymptotic sketch, not used elsewhere).** Note that the sums run over `v in Z^r`, so `n_0` grows with `r`, not
with `m`: `n_0 ≈ 4 r z/log tau_1 + O(log(m r b*)/log tau_1)` for large `tau_1`. Hence `log(1/first term) ≈ m n_0 (log P)/2`.
For a realisation coset `rho ≈ kappa_inv C e^{-W/4}`, and `log H_s(rho) ≈ -s W/4 + O_s(r)`. At barrier scale, with
`W ≈ r e-bar log P`, the condition becomes roughly `m < c e-bar sqrt(tau_1) log tau_1`. Gaussian seeds need `m` distinct
`q_k` and `tau_1^2 <~ P/(2 max q_k)`, so in the window the proved range is about `m <= M^{(2-delta)/5 - o(1)}` (with
`(1)`, `(2)` below as hypotheses). Therefore:

* *Scenario P, proved for small free dimension, conditionally on the profile.* Suppose a barrier-scale window profile
  satisfies Corollary 7.3's inequality: in particular `m` is small as above, (1) its characters hold on some
  `Real(delta)`, and (2) the base coset carries a seed as in Lemma 7.2. Then `G(p) ∩ Real(delta)` has positive Haar
  measure, and its points are not arithmetic. They lie at distance about `Delta` from `C_0`, and arithmetic points of
  `Real(delta)` are actual clusters (Proposition 6.1). Existence of such profiles in the window is **not** proved here.
  (1) is the row-3 profile-input problem, and by Lemma 7.2 it forces every character to be heavy.
* *For `m >= cM` the proved bound fails, by a power of `M`.* The deterministic local estimate, which needs every
  `v in Z^r` with `|v|_1 < n_0` to be safe, is the loss. If, near the seed, the unsafe and heavy windows were placed like
  independent random slabs, their expected share of `Q_eps` would be about their total mass (`e^{-W/4 + O(r)}` for the
  heavy ones), and positivity would follow for every `m` (*heuristic B*, not proved).

### 7.2 Two heuristics, and what the computations say (evidence)

* *Heuristic A (cosets far from seeds are empty).* For a coset that passes near no slack seed, treat the
  non-character windows as random slabs. The expected number of surviving cells of side `e^{-X/2}` is about
  `exp(m X/2 - S(X))`, with `S(X)` the window mass of monomials of weight `<= X`. On the full torus with primes in
  `[P/2, P]` this predicts that, once `Lambda >~ (1/4) log P + O(log log P)` (where Proposition A' stops), all of `G(p)`
  lies near seeds; `G(p)` itself is never empty there (Theorem 2.2). In the window `S ~ e^{2 Lambda} >> m Lambda log P`. *Evidence* (§8.2): for 8 primes `13..73`, all 60 random lines `g + R b` are
  covered by non-character windows of weight `<= 22`. This is a certified emptiness of `G ∩ line` for each of them, and 9
  of them pass the character test. The line through a transplant vector translated by `eta = 1e-3` is also covered
  (cover margin `4.4e-6`), although its characters hold.
* *Heuristic B (seeded base cosets survive).* Near a seed, light monomials are protected deterministically and only
  heavy ones are random. *Evidence* (§8.2): all 30 lines through the transplant vector keep surviving length
  (`>= 1.6e-5` at `T = 22`), and the survivors sit within `1.6e-4` of it. On the base line `R(1, ..., 1)` (Gaussian seed
  `arg(2+i)`, all `5 < p_i`) the non-character survivors have length `0.302`, and the translates keep `0.29`
  (`eta = 1e-3`) and `0.18` (`eta = 1e-2`) before dying at `eta = 0.1`. Caveat: in this toy `L^sat = (1,...,1)^perp` has
  light vectors, so its characters would fail. Lemma 7.2 shows that a genuine profile seeded in this way must have only
  heavy characters. The toy also overstates the survival because all `v` of a fibre `v + L^sat` share one centre on the
  base line, which is not so for heavy `L^sat`.

The heuristics are consistent. The Haar-generic part of `G` disappears in the window (A). What survives is
structured: neighbourhoods of arithmetic points on the full torus, and `Delta`-neighbourhoods of seeds on base cosets
(B). This replaces the Poisson-law reasoning behind scenario P in `fakes.md` §6.3, which Remark 2.3(2) shows is false at
large `Lambda` on the full torus.

### 7.3 The reduced form of (Q)

Theorems 3.1 and 5.1, Proposition 6.1 and §7.1 reduce the `A_0` window question, **for the seeding mechanism**, to:

> **(Q_seed).** For `M -> infinity`, are there barrier-scale profiles (`W = O(M log M)`), with all primes
> `<= M^(2-delta)` and `m >= cM`, such that (1) some `Real(delta)` satisfies all character conditions, with every
> `v in L^sat \ 0` heavy (`w(v) >= W/2 - O(log W)`, Lemma 7.2), and (2) the base coset carries a non-character slack seed
> (for example the Gaussian seeds of Lemma 7.2)? And, given such a profile, is the local measure near the seed larger
> than the heavy-window mass `e^{-W/4 + O(r)}`?

A yes to both parts gives scenario P (`A_0` capped in the window). The first part is a design problem on profiles
(row-3 inputs). The second is the gap between the deterministic bound of Corollary 7.3 (first term
`~ exp(-m n_0 (log P)/2)`, `n_0` growing with `r`) and heuristic B. A no to the first part, for every
profile, would leave only unseeded cosets, where heuristic A predicts emptiness (scenario R). Neither answer gives growth
(J2).

## 8. Computations

All scripts are in `coset_checks/`, written for this round (`common.py` holds the shared enumeration of monomials up to
sign and the exact window half-widths `(2/pi) arcsin`). They run under `prlimit --as=3500000000`. Normalised coordinates
`x = phi/(pi/2) in R/Z` are used. "Truncated at `T`" means that all monomials with `w(v) <= T` are imposed. A truncated
test is an **upper** bound for the true gap set, so "covered at `T`" certifies emptiness, while "survives at `T`" is
evidence only.

### 8.1 Density at a transplant vector (`density_check.py` -> `density_check_out.json`; evidence for Lemma 2.1)

Primes `13, 17, 29, 37, 41, 53, 61, 73`, transplant `q = 5, 13, 17, 29, 37, 41, 53, 61` (`tau^2 = min p/q = 73/61`),
`T = 24` (30,228 monomials up to sign). The transplant vector passes every test, with minimum margins `1.3e-5` (0-family)
and `2.6e-6` (`pi/4`-family). Survival rate of uniform points in the cube of half-side `eps` around it (4000 samples):

| `eps` | `1e-1` | `1e-2` | `1e-3` | `1e-4` | `1e-5` | `1e-6` | `1e-7` | `1e-8` |
|---|---|---|---|---|---|---|---|---|
| survival | 0 | 0 | 0.00025 | 0.0715 | 0.473 | 0.915 | 1.0 | 1.0 |

Uniform points on `T^8` survive at rate 0 (4000 samples). The density tends to 1 as Lemma 2.1 predicts, but only at
scales `<= 1e-6`, which reflects the small slack `tau = 1.094`.

### 8.2 One-dimensional cosets (`line_cosets.py`, `line_verify.py`, `line_batch.py`, `seeded_base.py`)

For a line `g + R b` (`b in Z^8` primitive), each monomial gives the frequency `n = <b, v>` and the shift `c = <v, g>`. If
`n = 0` it is a character, checked at `g`. Otherwise its windows are exactly `2(|n| + 2)` intervals on `[0, 1)`, and their
union is computed by sorting.

*(a) Lines with direction `b = (1, -1, 2, 0, 1, -2, 1, 1)`* (`line_cosets_out.json`). "char" is the minimum character
margin (negative means a character fails, which empties the line trivially); "len" is the surviving length.

| `g` | `T = 16` char / len | `T = 20` | `T = 24` | `T = 28` |
|---|---|---|---|---|
| transplant `x_beta` | `4.6e-4` / `2.9e-4` | `1.5e-4` / `3.0e-5` | `9.7e-6` / `2.0e-5` | `9.7e-6` / `1.7e-5` |
| `x_beta + 1e-5 u` | `4.7e-4` / `2.9e-4` | `1.6e-4` / `3.0e-5` | `1.3e-5` / `1.0e-5` | `1.3e-5` / `8.3e-6` |
| `x_beta + 1e-3 u` | `7.2e-4` / `1.6e-4` | `5.6e-5` / **0** | `5.6e-5` / **0** | `3.1e-6` / **0** |
| `x_beta + 1e-1 u` | `< 0` / 0 | `< 0` / 0 | `< 0` / 0 | `< 0` / 0 |
| uniform `g` (4 draws) | `< 0` / 0 | 1 of 4 with char `> 0`, all len 0 | `< 0` / 0 | `< 0` / 0 |

The non-character window mass is `S_nc = 12.7, 14.3, 15.0, 15.3` and the number of intervals is `1.1e4` to `6.9e5`. On the
`eta = 1e-3` line the characters hold (margin `5.6e-5`), yet the non-character windows cover the line at `T = 20`. The
smallest overlap of consecutive covering intervals is `4.4e-6`, far above double-precision error (`line_verify_out.json`).
So `G(p) ∩ (x_beta + 1e-3 u + R b) = ∅`: **exact** (a finite cover by open intervals with margin), not evidence. On the
transplant line, all survivors lie within `theta ∈ (-1.6e-4, 5e-5)` of the transplant point (mod 1).

*(b) Batch at `T = 22`* (`line_batch_out.json`): 60 lines with random `b in [-3, 3]^8` and uniform `g`: 9 pass every
character, **0** keep any non-character survivor. For all 60 the non-character windows cover the line with overlap
margin `>= 6.0e-4` (`batch_margins.py` -> `batch_margins_out.json`), so each of these lines has `G(p) ∩ line = ∅` (exact). 30 lines through `x_beta` with random `b`: all 30 survive, with minimum
length `1.6e-5`, and all pass the characters.

*(c) Seeded base line* (`seeded_base_out.json`; heuristic B, with the caveat of §7.2): `b = (1, ..., 1)`, `T = 22`,
characters ignored; `eta u` is a translate in a random direction with `u ⊥ b`, `|u|_inf = 1`.

| `eta` | 0 | `1e-7` | `1e-5` | `1e-4` | `1e-3` | `1e-2` | `1e-1` |
|---|---|---|---|---|---|---|---|
| surviving length (3 directions) | 0.3019 | 0.3019 | 0.3018 | 0.3010-0.3013 | 0.2909-0.2918 | 0.176-0.191 | 0 |

At `eta = 0` the value is `mu_1(G_1(13)) = 0.3019` (`fakes.md` §5.3). On the base line every fibre has a single centre, and
its lightest member is `n e_1`.

### 8.3 Constants (`quant_bound.py`, `barrier_M0.py`, `graph_formula_check.py`)

* Corollary 2.3 (`quant_bound_out.json`). `r = 8`, `P = 73`, `tau = 1.094`: `n_0 = 939`, `log mu >= -1.6e4`. `r = 2500`,
  `P = 10^6`, `tau = 2`: `n_0 = 12772`, `log mu >= -2.2e8`, where the Poisson extrapolation gives `-162`. The proved bound
  overtakes the extrapolation only when `Lambda >~ log r`, i.e. far beyond feasible sizes (Remark 2.3(2) is asymptotic).
* Corollary 5.2 thresholds (`barrier_M0_out.json`): `M_0 ≈ 4.5e3` `(c, delta) = (1/4, 1)`; `4.6e5` `(1, 1/2)`; `7.0e8`
  `(1/4, 1/2)`; `4.9e17` `(1/2, 1/4)`; `4.3e20` `(1/4, 1/4)`.
* Corollary 4.2 formula (`graph_formula_check_out.txt`): relative error `5.6e-8` against brute force.

## 9. Ledger, constraints, open points

### 9.1 Ledger

| statement | status |
|---|---|
| Lemma 1.1 (coordinates on `H_0`; characters `= L^sat`) | proved |
| Lemma 1.2 (Gaussian primes of distinct split norms give points of `G(q)`) | proved (special case of `fakes.md` Lemma T) |
| Lemma 2.1 (slack points are Lebesgue density points, on any coset) | proved |
| **Theorem 2.2** (`mu(G(p)) > 0` iff `5 ∉ p`; no condition on `Lambda`) | **proved** |
| Corollary 2.3 (`mu(G(p)) >= exp(-C_tau r^2 log P)`, explicit) | proved |
| Remark 2.3(2): the Poisson extrapolation, and `fakes.md` open point 3 for the coset `T^r`, are false at large `Lambda` | proved (asymptotic; crossover only at `Lambda >~ log r`) |
| **Theorem 3.1** (sandwich: slack existence => positive Haar measure => existence, on every coset) | **proved** |
| Theorem 4.1, Corollary 4.2 (Proposition A' on cosets; exact step sums in graph coordinates) | proved; formula checked exactly |
| **Theorem 5.1, Corollary 5.2** (no union-bound certificate, in any basis or order, reaches the window when `m >= cM`; explicit `M_0`) | **proved** |
| Proposition 6.1 (arithmetic points of `Real(delta)` are actual clusters with `log N_beta >= 3M log M - O(M)`) | proved (referee form of Lemma T plus the `(2/3)` theorem) |
| Proposition 7.1, Lemma 7.2, Corollary 7.3 (perturbation from a seeded base coset; Gaussian seeds force heavy characters) | proved |
| Reading of Corollary 7.3 (`m <~ M^{(2-delta)/5}` proved range, under hypotheses (1), (2)) | sketch |
| Heuristic A (unseeded cosets empty in the window), heuristic B (seeded base cosets survive) | heuristic |
| Emptiness of the specific lines of §8.2(a) at `eta = 1e-3` and of the 60 random lines of §8.2(b) | exact (finite covers with margin) |
| Density near transplants, survival near seeds | evidence |

### 9.2 Owner constraints

* *No uniform exceptional sets from fixed-`S` Subspace theorems.* None are used.
* *Gaussian gcd and denominator costs.* The torsion index `D = [L^sat : L]`, the torsion characters and the right inverse
  of `(lambda_x)` (`kappa_inv`) are kept as explicit hypotheses (`(Ch)`, `(Ch_eta)`, `rho`). Lemma 7.2 restricts seeds to
  `(pi/4)`-torsion base cosets and states what that forces.
* *No drifting common phase set to zero.* Realisation is used in the coset form `Real(delta)`, for every `delta`, with
  `theta_0` absorbed exactly as in `fakes.md` §1.1. The translate `eta` between a base coset and a realisation coset is
  kept and enters through `Heavy(eta)`.
* *No template-dependent threshold treated as uniform.* `M_0(c, delta)` in Corollary 5.2 depends only on `(c, delta)`
  and is computed. The conditions of Corollary 7.3 depend on the profile (`B`, `b*`, `kappa_inv`, `tau`). They are stated
  as profile hypotheses, and nothing uniform in `R` is claimed from them.

### 9.3 Open points, smallest first

1. *(Q_seed), design part.* Are there barrier-scale profiles with primes `<= M^(2-delta)`, `m >= cM`, all characters
   heavy and valid on some `Real(delta)`, whose base coset carries a Gaussian seed (Lemma 7.2: a basis of `K` and
   distinct `q_k` with `prod_k q_k^{|B_ik|} <= p_i/(2 tau^2)`)? Grouped profiles (blocks of identical columns) give
   sparse bases of `K`, but their characters then contain the representative profile, which is pinned. This connects
   (Q_seed) to the pinned case (`fakes.md` open point 5).
2. *(Q_seed), analytic part.* Close the gap between the deterministic local bound `exp(-m n_0 (log P)/2)` and the
   heavy-window mass `e^{-W/4 + O(r)}`. An equidistribution statement for the arguments `arg gamma_mu` of heavy Gaussian
   monomials, relative to the shifts `<v, eta>`, would suffice for heuristic B.
3. *Scenario R for unseeded cosets.* Prove heuristic A for an infinite family: for example, that a fixed line `g + R b`
   with non-seeded `g` is covered by the windows of weight `<= T(Lambda)`. §8.2 gives exact covers only at fixed small
   sizes.
4. *Base cosets with torsion of order `>= 3`.* Their shifts are cyclotomic, and Gaussian seeds lose their Liouville
   bound. Angles of Weil numbers in CM fields of degree `2d` satisfy the gap conditions with weights `q^d <= p` (a
   Liouville argument in the field, *remark, not checked in detail*). Whether such seeds exist at barrier scale is open.
5. The pinned case `m = 0` is untouched, as before.
