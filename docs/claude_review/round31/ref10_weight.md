# Referee report, Part X, block "weight": X.4 and the calibrations paragraph of X.5

Scripts (own code, exact integer arithmetic except the floating logarithms of the lcm table):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref10_weight/lcm_check.py`,
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref10_weight/tuples.py`,
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref10_weight/bounded_rank_example.py`.
Astra's scripts were not imported.

## 0. Verdict in five lines

1. (X.4) is correct, with constant one and with a slightly stronger form: at most `D(x)` labels have
   `log p_j < theta_3(x)`. It is a one-line consequence of (1.2) at unit one on `v = e_j - e_k`. **Proved.**
2. The restricted growth rate is correct as a logical implication. The lcm estimate `log L(x) = o(x)` is correct,
   with one wrong but harmless intermediate bound. Numerically `log L(x)/x` falls from 0.14 at `x = 100` to 0.073 at
   `x = 10^6`. **Proved, Computed.**
3. The bounded-rank hypothesis enters only through `D(K log m) = o(m)`. The corrected statement uses that condition
   directly (§1.2).
4. For genuine circles the hypothesis and the conclusion are statements about the prime support alone. The cluster
   plays no role. The class of genuine parents meeting the hypothesis is nonempty, but on it all but `o(m)` primes
   exceed `m^{K/2 - o(1)}` for every `K`, so the conclusion is immediate from the sizes of the primes. The natural
   parents (all split primes in `[P, 2P]`, the first `s` split primes) fail the hypothesis: their rank is maximal,
   equal to `pi_3(x)`, and tends to infinity. **Proved.**
5. Verdict: X.4 is a correct reformulation, not a restricted rate theorem with content for genuine circles. For
   fakes it is a necessary condition on witnesses: a fake with `W <= A m log m` must have nearly injective tuples,
   `D(K log m) >= s - (2A/K + o(1)) m`. Generic residues satisfy this, so it excludes no known kind of fake and
   does not touch P3-all\*.

## 1. Statements, corrected where necessary

Notation, as in `sharp.md` Def. 1.1. A member of `F^all_X(M, C)` has distinct primes `p_j = 1 mod 4` (the `s`
labels), angles `phi_j` and residues `r_j(q)` of norm `p_j` at every odd prime `q <= X`, `q != p_j`. The weight of a
monomial is `w(v) = sum_j |v_j| log p_j`. `W = sum_j e_j log p_j >= sum_j log p_j`. For an inert prime `l = 3 mod 4`,
`Z[i]/l = F_{l^2}`, and `omega_j(l) = r_j(l)/conj(r_j(l))` lies in the norm-one group `mu_{l+1}`, which is cyclic.
`T_j(x) = (omega_j(l))_{l <= x, l = 3 (4)}` is the tuple of label `j`, and `D(x)` is the number of distinct tuples
among the `s` labels. `theta_3(x) = sum_{l <= x, l = 3 (4)} log l` and `pi_3(x)` is the number of these primes.

### 1.1 (X.4) — Proved (CORRECT, with a stronger form)

**Statement.** Fix a member of `F^all_X(M, C)` and `x <= X`. Then
`#{j : log p_j < theta_3(x)} <= D(x)`, and hence `W >= sum_j log p_j >= (s - D(x)) theta_3(x)`.

**Proof (checked line by line).** Suppose `T_j(x) = T_k(x)` with `j != k`. Take `v = e_j - e_k` and `s = 0` in
(1.2). Then `prod_i (r_i/conj r_i)^{v_i} = omega_j/omega_k = 1 mod l` for every inert `l <= x`. Every such `l` is
`<= X` and does not divide `N` (the `p_i` are `1 mod 4`), so `a_(l,0)(v) >= 1`, and `m_(v,0) >= prod_{l <= x} l =
e^{theta_3(x)}`. Here `m_(v,0)` is the product over all colliding prime powers `q^a <= X`, `q` not dividing `N`.
Using only the inert primes `l <= x` is a lower bound, and higher powers and split primes can only increase it. The
constant is `nu_0 = 1`, and `|sin| <= 1` gives `e^{theta_3(x)} <= m_(v,0) <= e^{w(v)/2} = sqrt(p_j p_k)`, that is,
`(log p_j + log p_k)/2 >= theta_3(x)`. So two labels in one tuple class cannot both have `log p < theta_3(x)`, and
each class holds at most one label below the threshold. There are `D(x)` classes.

**Bookkeeping in Astra's form, checked.** In a class `c` with `|c| >= 2`, summing the pair inequality over its
`|c|(|c|-1)/2` pairs gives `(|c| - 1) sum_{j in c} log p_j >= |c|(|c| - 1) theta_3`, hence
`sum_{j in c} log p_j >= |c| theta_3`. Summing over the non-singleton classes gives
`sum_j log p_j >= (s - #singletons) theta_3 >= (s - D(x)) theta_3`. Correct.

**Requirements.** The cutoff must satisfy `X >= x`, or more exactly `X >` the largest inert prime `<= x`. Only the
unit-one tests on the differences `e_j - e_k` at those levels are used. The parent, witnesses, angles and cutoff are
fixed first, and the statement holds for every `x <= X`.

**Novelty.** Minimal. The statement is (1.2) with `|sin| <= 1`, as in the proof of `sharp.md` §1.2 (`m_(v,s) | D_s`).

### 1.2 The restricted growth rate — Proved (CORRECT WITH CORRECTIONS)

**Corrected statement (quoted in full).** *Fix `C`. Consider a sequence of members of `F^all_{X}(m, C)` with
`m -> infinity` and cutoffs `X = X(m) >= K log m` for every fixed `K` (for example `X = 3m`). Then for each fixed
`K > 0`*
`liminf W/(m log m) >= (K/2) liminf (s - D(K log m))/m`.
*In particular, if `D(K log m) = o(m)` for every fixed `K`, then `W/(m log m) -> infinity`. A sufficient condition
for `D(K log m) = o(m)` is Astra's: the subgroup of `prod_{l <= K log m} mu_{l+1}` generated by the `s` tuples is
generated by at most `d_K` elements, with `d_K` independent of the instance. Contrapositive: if
`W <= A m log m` along the sequence, then `D(K log m) >= s - (2A/K + o(1)) m` for every fixed `K`.*

**Checks.**
* *Group bound.* `prod_l mu_{l+1}` has exponent `L(x) = lcm_{l <= x, l = 3 (4)} (l + 1)`. An abelian group of
  exponent dividing `L` generated by `d` elements has at most `L^d` elements, and `D(x)` is at most the size of the
  subgroup generated by the tuples. Correct.
* *`log L(x) = o(x)`.* Correct, with one wrong intermediate bound.
  * The primes `r <= x/B` contribute `sum_{r <= x/B} floor(log(x+1)/log r) log r <= pi(x/B) log(x+1)
    = (1 + o(1)) x/B`. This is `O(x/B)`, as needed. It is not `psi(x/B)` as written in the answer: a prime
    `r <= x/B` can occur to a power up to `x + 1`, not only up to `x/B`. Wrong wording, right order. The work notes
    contain the same slip.
  * Prime powers `r^2 | l + 1` with `r > x/B` need `r <= sqrt(x+1) < x/B` once `x > B^2 + 1`, so there are none. The
    stated `O(sqrt x log^2 x)` is a harmless overestimate.
  * A prime `r > x/B` with `r | l + 1 <= x + 1` has `l + 1 = h r` with `h < B + 1`. The answer's `h <= 2B` also works.
    It occurs to the first power and contributes `<= log(x+1)`. For fixed `h`, the primes `r` with `hr - 1` prime
    have upper relative density `<= prod_{t <= z, t odd, t not dividing h} (1 - 1/(t-1))` among primes. This holds
    by PNT in the fixed progressions modulo `prod t`: for each `t`, the excluded class `r = h^{-1} mod t` is one of
    the `t - 1` invertible classes. The `O_z(1)` exceptions `hr - 1 = t` do not matter. The product tends to 0 as
    `z -> infinity`. So the count is `o(x/log x)` and the contribution is `o(x)`. Summing over `h <= 2B` and letting
    `B -> infinity` gives `limsup log L(x)/x <= c/B` for every `B`. **Proved.** (A Brun or Selberg upper-bound sieve
    gives the sharper `O_h(x/log^2 x)` count directly. Heuristically `log L(x)` has order `x loglog x/log x`.)
* *Numerics (Computed, exact lcm by factorisation, `lcm_check.py`).*

  | `x` | `log L(x)` | `log L(x)/x` | `theta_3(x)/x` | share of primes `r > x/10` | share of primes `r > x/100` |
  |---:|---:|---:|---:|---:|---:|
  | 10^2 | 14.4 | 0.1445 | 0.439 | 0.052 | 0.145 |
  | 10^3 | 108.5 | 0.1085 | 0.495 | 0.030 | 0.089 |
  | 10^4 | 964.9 | 0.0965 | 0.498 | 0.026 | 0.084 |
  | 10^5 | 8267.0 | 0.0827 | 0.4995 | 0.021 | 0.072 |
  | 10^6 | 72934.7 | 0.0729 | 0.5002 | 0.017 | 0.063 |

  The ratio decreases slowly, as claimed. The decay is slow: at these heights `L(x)^d` is about `e^{0.07 d x}`. So
  `D <= L(K log m)^{d_K} < m^{1-eps}` holds only once `0.07 K d_K < 1`, or asymptotically. The statement is
  asymptotic, so this is not an error.
* *Limit.* With `x = K log m`: `log L(x)^{d_K} = d_K o(K log m) = o_K(log m)`, so `D(x) <= m^{o(1)} = o(m)` for fixed
  `K`, `d_K`. Here `o(1) -> 0` because `x -> infinity` with `m`. By PNT mod 4, `theta_3(x) = (1/2 + o(1)) K log m`.
  The rank bound `m <= s + 1 + K_C` (Theorem 1.2 / VIII.1, with `K_C` depending only on `C`) gives
  `s >= m - O_C(1)`. Hence `W >= (m - o(m))(K/2 + o(1)) log m`. The cutoff needed is `K log m < 3m` for
  `m >= m_0(K)`. Correct.
* *Quantifiers.* `C` is fixed, then `K` and `d_K`, then `m >= m_0(C, K, d_K)`. The sequence fixes the parent,
  witnesses, angles and cutoff of each instance.

**Correction.** The hypothesis "bounded number of generators" is used only to get `D = o(m)`. The corrected
statement above, in terms of `D`, is strictly more general and shows what X.4 really constrains.

### 1.3 Meaning for genuine circles — the content is a property of the prime support

**(a) Proved.** For genuine factors `r_j = pi_j`, `omega_j(l) = omega_k(l)` means
`pi_j conj(pi_k) = conj(pi_j conj(pi_k)) mod l`, that is, `l | Im(pi_j conj pi_k)`. This integer is nonzero for
`p_j != p_k`, and `|Im| <= sqrt(p_j p_k)`. So the equality of tuples at all inert `l <= x` forces
`e^{theta_3(x)} <= sqrt(p_j p_k)` for **any** two split primes, with no cluster, no arc and no angles. Hence, for
genuine factors, (X.4) and its stronger form (`#{j : log p_j < theta_3(x)} <= D(x)`) hold for every finite set of
split primes with any choice of Gaussian factors. Astra says this herself ("principally forces the support itself to
have large weight"; work notes: "every actual dense configuration of distinct physical primes, irrespective of the
cluster, has growing residue rank"). The answer's wording undersells it: for genuine circles the cluster plays no
role at all.

**(b) Vacuity on natural parents. Proved, Computed.**
* *Proved.* If `theta_3(x) > (log p_j + log p_k)/2` for all pairs, then `D(x) = s` and (X.4) reads `W >= 0`. For the
  first `s` split primes (`p_s ~ 2 s log s`) this happens once `theta_3(x) > log p_s`, that is, at
  `x = (2 + o(1)) log s`. So for every `K > 2`, `D(K log s) = s` eventually and (X.4) is vacuous. The same holds for
  all split primes in `[P, 2P]` once `e^{theta_3(x)} > 2P`.
* *Computed (`tuples.py`; Gaussian factor `a + bi`, `a` odd, `b > 0` even; exact discrete logarithms in
  `mu_{l+1}`; exact rank by lattice determinants, `d_r(H) = log_r(|H_r|/|rH_r|)`).*

  | parent | `s` | first `x` with `D(x) = s` | rank of the tuple group at `x = 3, 7, ..., 67` |
  |---|---:|---:|---|
  | split primes in [100, 200] | 10 | 7 | 1, 2, ..., 10 = `pi_3(x)` |
  | split primes in [500, 1000] | 36 | 19 | `pi_3(x)` at every `x` |
  | split primes in [1000, 2000] | 67 | 19 | `pi_3(x)` at every `x` |
  | split primes in [5000, 10000] | 280 | 19 | `pi_3(x)` at every `x` |
  | first 100 split primes | 100 | 27 (`K ~ 5.9`) | `pi_3(x)` |
  | first 1000 split primes | 1000 | 27 (`K ~ 3.9`) | `pi_3(x)` |

  In every same-class pair found (up to 124,540 pairs per row), `prod_l l` divides `Im(pi_j conj pi_k) != 0` and
  `(prod l)^2 <= p_j p_k`, with 0 violations. The rank is always the largest possible value `pi_3(x)`: every
  subgroup of a product of `pi_3(x)` cyclic groups has rank at most `pi_3(x)`, and the 2-rank is already full. The
  2-component of `omega_j(l)` is the Legendre symbol `(p_j/l)`, since `omega^{(l+1)/2} = pi^{-(l^2-1)/2} =
  p^{-(l-1)/2}`. This identity was verified exactly in every case.
* *Proved (the hypothesis fails on natural parents).* Fix `x_0`. For all split primes in `[P, 2P]` and `P >= P_0(x_0)`,
  PNT in the fixed modulus `4 prod_{l <= x_0} l` gives every sign pattern `((p/l))_{l <= x_0}`. So the F_2-span of
  the Legendre vectors is full, and the rank at level `x_0` is `pi_3(x_0)`. Projection does not increase the minimal
  number of generators, so the rank at `x = K log m >= x_0` is at least `pi_3(x_0)`. Hence for every fixed `K` the
  rank is unbounded along `P -> infinity`, and the dyadic parent fails the hypothesis. The same holds for any parent
  family whose primes equidistribute in fixed residue classes.

**(c) Nonemptiness for genuine parents, and why it is empty of content. Proved, Computed.**
* *Example (Proved).* Let `Q_y = prod_{l <= y, l = 3 (4)} l`, and take any split primes `p = a^2 + (Q_y b)^2`
  (infinitely many by Hecke's prime ideal theorem in the ray class mod `Q_y`). Then `pi = a + Q_y b i` has
  `omega = 1 mod l` for every `l <= y`, so `D(y) = 1` and the rank is 0. Choosing `y = y(m)` with
  `y/log m -> infinity` gives genuine parents with `s -> infinity` that meet the hypothesis for every `K`, with
  `d_K = 0`. *Computed* (`bounded_rank_example.py`): `y = 23`, `Q = 100947`, 60 such primes, `D = 1`, rank 0,
  smallest prime `4.08 * 10^10 > Q^2`.
* *Proved.* Conversely, on any genuine parent meeting the hypothesis, (a) together with `D = m^{o(1)}` shows that all
  but `m^{o(1)}` primes satisfy `log p_j >= theta_3(K log m) = (K/2 + o(1)) log m`, for every `K`. So the conclusion
  `W/(m log m) -> infinity` follows from the sizes of the primes alone. No lattice-point geometry is used, and it is
  no growth-rate statement about clusters. Whether such parents carry clusters of `M ~ s` points on short arcs is
  irrelevant to the conclusion and was not examined.
* *Answer to the crucial question.* For genuine circles the hypotheses can be satisfied only on parents with
  superpolynomially large primes (relative to `m`), and there the conclusion is trivial. On the natural parents the
  hypothesis is false. **The restricted rate theorem has no content for genuine circles.**

**(d) For fakes. Proved, Heuristic.** For fakes the hypothesis concerns the witnesses `r_j`, which are free up to the
norm. The theorem says (contrapositive of §1.2) that a fake with `W <= A m log m` must have
`D(K log m) >= s - (2A/K + o(1)) m`. Its tuples must therefore be almost injective at all levels up to `K log m`.
Since `|prod_{l <= K log m} mu_{l+1}| = m^{K/2 + o(1)}` is much larger than `s`, generic or genuine-like residues
satisfy this. The IX.3 witnesses `r_j = i^{u_j} pi_j prod gamma_l^{V_jl}` carry the independent factors `pi_j`, and
their rank is not bounded. (Heuristic for the IX.3 witnesses: not computed here.) X.2's deliberately failing omitted
difference is the same mechanism used in the opposite direction (equal ratios at one prime force a failure), and it
is consistent with X.4.

### 1.4 Calibrations paragraph of X.5

* *Function-field transfer of (X.4) and of the exponent estimate.* The arithmetic is **Proved** in the setting of the
  constant-field extension `F_{p^2}(T)/F_p(T)`. An irreducible of degree `j` is inert exactly for odd `j`; its
  norm-one group is cyclic of order `p^j + 1`. `log lcm_{j <= k, j odd}(p^j + 1) <= sum_{j <= k} log(p^j + 1) =
  O_p(k^2)`. The sum of the log norms of the inert irreducibles of degree `<= k` is `>> p^k log p` by the
  irreducible-polynomial count. So `log L = o(theta)` holds, even polylogarithmically. Given the analogue of the
  unit-one gap (a determinant or degree bound), the implication transfers verbatim. I did not check that "the
  previously constructed field families already have `W/(m log m) -> infinity`", or the precise definition of the
  field analogue of the gap: **Not checked.** As calibration, the transfer shows only that the implication does not
  prove too much. In both settings the content is the hypothesis.
* *"In T5 the comparable norm-one ratios ... are absent."* I could not locate the definition of the "twin relaxation
  T5" that this sentence refers to (the record's T1–T5 in round 27 are pre-registered tests of a different kind).
  **Not checked.**
* *"Fixed-field examples do not satisfy the local hypotheses; groups at prime powers are generally not cyclic."* For
  Z[i] at odd `q`, the norm-one group mod `q^a` is cyclic of order `(q -+ 1) q^{a-1}`. In `F_q[T]` the 1-units mod
  `P^a` form a `p`-group of large rank, not cyclic for `a >= 2` in general. **Correct** as a remark. The order
  argument `d_a | 4nG` of X.1 does rely on cyclicity. That X.1 is refereed elsewhere.
* The remaining calibration sentences (X.1 tests true factors as equalities; IX.3 data use a large rational lift;
  X.2 has no cofinite rational lift at each label by IX.7) belong to blocks refereed elsewhere. **Not checked here.**

## 2. Method

* (X.4): line-by-line check against `sharp.md` (1.2), including `nu_0 = 1`, the definition of `m_(v,s)` as a
  product over all colliding prime powers `<= X` not dividing `N`, and `w(v) = sum |v_j| log p_j`.
* `lcm_check.py`: a smallest-prime-factor sieve to `10^6 + 2`. For each prime `l = 3 mod 4`, factor `l + 1` and keep
  the maximal exponents.
* `tuples.py`: for each prime, the Gaussian factor is found by search. `omega = pi * conj(pi)^{-1}` mod `l` in
  `F_l[i]`. The discrete logarithm uses a generator `h^{l-1}` of `mu_{l+1}`, checked to have exact order `l + 1`.
  `D` is the number of distinct tuples. The exact minimal number of generators is `max_r d_r(H)`, with
  `d_r(H) = log_r(det(rV + diag)/det(V + diag))` computed by integer Hermite elimination modulo `r^{max}`. The
  Legendre-symbol identity and the determinant divisibility are asserted for every label and every same-class pair.
* `bounded_rank_example.py`: primes `a^2 + (Q_y b)^2`.

## 3. Labels summary

| item | label |
|---|---|
| (X.4), constant, cutoff requirement, bookkeeping | Proved (correct; stronger form available) |
| `log L(x) = o(x)` | Proved (the `psi(x/B)` step is misworded; the order is right) |
| `log L(x)/x` table to `10^6` | Computed (decreasing: 0.145 to 0.073) |
| restricted rate under bounded rank | Proved; corrected to the `D(K log m) = o(m)` form |
| genuine factors: (X.4) is a statement about the prime support alone | Proved |
| natural genuine parents: rank `= pi_3(x)` (computed), unbounded (proved); (X.4) vacuous for `K > 2` | Proved, Computed |
| genuine parents meeting the hypothesis exist, but only with primes `>= m^{K/2 - o(1)}` for every `K` | Proved |
| "restricted rate theorem with content for genuine circles" | False as a description; it is a reformulation |
| function-field transfer of the implication | Proved for the exponent arithmetic; field families Not checked |
| T5 sentence; other calibration sentences | Not checked |
