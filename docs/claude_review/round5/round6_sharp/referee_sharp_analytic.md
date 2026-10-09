# Referee report on `round6/sharp.md`: analytic estimates and union bounds

Target: `round6/sharp.md` ("fake circles for residue data that factor through the primes (P3)"), together
with every step of `round5/lower.md` that it uses (Lemma S, Prop. 3.3, Lemma 5.2, Lemma 5.3, Prop. 1.2/1.3).
Lens: analytic estimates and union bounds. Background read: `round4/outside.md` §§0-2, `two_thirds/proof.md`
§§2-3, `round5/construct.md` §0.

Referee scripts are in `round6/referee_analytic_checks/`. Each runs from that directory with
`python3 <script>`, and its output is stored next to it as `*_out.txt` (or `rerun*_<script>.txt`).
* New in this pass: `a_mu.py`, `a_height.py`, `a_p65.py`, `a_thm51.py`, `a_thm64.py`, `a_appA.py`.
* Written in an earlier pass of this referee task and rerun now: `r_lemmaS.py`, `r_slack.py`, `r_fibre.py`,
  `r_gs.py`, `r_struct.py`, `r_lemma33_mc.py`, `r_thm51.py`, `r_thm64.py` (outputs `rerun3_*.txt`).
  `r_height.py` has an older output only; `a_height.py` replaces it.
* The author's scripts in `round6/sharp_checks/` were rerun unchanged (outputs `rerun2_check_*.txt`). All
  five print their PASSED lines, with the numbers quoted in `sharp.md` §8.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The main results survive.** That is Theorem 5.1 (P3-char, `W <= (b max(A,2) + o(1)) M log M`, `b >= 33`),
Corollary 5.2, and Theorem 6.4 (P3-all, `W <= (2Ab + o(1)) M log^2 M/loglog M`). I re-derived every margin,
moment bound and union bound and found them correct. There is one sub-step gap (the `mu_M` bound for pair
demand), and it is easily closed.

**One major error, in the sharpness part (§6).** Prop. 6.1 states the height condition correctly, for each
unit separately: `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`. But the Verdict (item 5), the proof of Prop. 6.5, and the
evidence script `check_p3all.py` use an **aggregated** form, `m_v = prod_q q^(max_s a_(q,s)) <= sqrt2 e^(w/2)`.
That form is **not** a necessary condition. Actual Gaussian primes violate it in 8976 of 8977 cases, and in
every one of 8000 two-prime monomials (`a_height.py`). For example, `pi = 3+2i` gives `m_v = 15 > sqrt 26`.
So "some column has `m_(e_j) > sqrt2 e^(w_j/2)`", as written, contradicts nothing.

The proposition survives after a per-unit repair (§5 below), with the same threshold. At accessible `M`,
the aggregated version overstates the failure range by a factor of about 1.8 (`a_p65.py`).

**Minor issues.**
* The scope of Prop. 6.5 (uniform random fibre only) is narrower than the Verdict's wording "per-modulus
  independent residue choices".
* Lemma 5.1b (= `lower.md` Lemma 5.2) has a restriction gap when a deleted prime is `<= X`. It does not
  affect the theorems, where all `p_j > X`.
* The status of Appendix A is overstated. Its use of the Green-Sanders citation is correct.
* `M_0` is astronomically large with the proved constants. This is not an error, but it is worth stating
  next to the "explicit" bound.

---------------------------------------------------------------------------------------------

## 1. Lemma S (`lower.md` Lemma S = `sharp.md` Lemma 4.2): re-derived, correct

I rederived each step myself. Notation: `T_a = sum_x c_x H(x,a)` over the `q_0` nonconstant columns,
`F = ||T||_1`, `n = ||c||_1`, `m_2 = ||c||_2^2`, `u_a = |T_a|/n`, `kappa = #{u_a > 1/2}`.

* **Defect identity.** Parseval over the full Hadamard matrix with `T_const = sum c = 0` gives
  `sum T_a^2 = M m_2`. Hence `sum |T_a|(n - |T_a|) = nF - M m_2`, and `m_2 >= n`.
* **Window (integrality).** `sum u - kappa = sum_(u<=1/2) u - sum_(u>1/2)(1-u)`. On the first range
  `u <= 2u(1-u)`, and on the second `1-u <= 2u(1-u)`. So `|sum u - kappa| <= 2 sum u(1-u) <= 2(F-M)/n`,
  which is `|F - kappa n| <= 2(F-M)`.
* **Inversion.** `c = H T/M`, including the zero constant coordinate. So `|c_x| <= F/M`, and amplitude
  `>= 2` gives `F >= 2M`.
* **Windows.** For ternary `c` with support `h >= 4` and `F < 17M/16`:
  `kappa h in (F - M/8, 3F - 2M) subset (7M/8, 19M/16)`.
* **`kappa >= 3`.** Then `h < 19M/48 < (M-2)/2` for `M > 9.6`, so Lemma 2.2 applies.
  * *Lemma 2.2 re-derived.* Reduce by the exact power of `q_0` dividing the finite part. Then
    `T_a = -G(a) mod q_0` off the support, with `G(Y) = sum c_x (Y-x)^d - c_inf`. The coefficient of `Y^k`
    (`1 <= k <= d`) is `binom(d,k) sum c_x (-x)^(d-k)`, and the Vandermonde argument on `s' <= d` nodes
    forces `G` nonconstant. So at most `s' + d` labels have zero transform.
  * Parseval allows at most `M/h` labels with `|T| = h`. Every other nonzero `T_a` is even, with
    `2 <= |T_a| <= h-2`.
  * Hence `F - M >= f(h) = (2(h-2)/h)(M/2 - h - M/h)`. This gives `f(4) = M/4 - 4` and
    `f(6) = (4/3)(M/3 - 6)`, and by concavity `f >= (3/2) min(3M/8 - 8, 5M/48 - 48/19)` on
    `[8, 19M/48]`. All three are `>= M/16` for `M >= 44`. The numerical minimum over `43 <= q <= 10^5` of
    `min_h f(h)/(M/16)` is 1.975 (`r_lemmaS.py` D1).
* **`kappa = 1`.** `1 - u_a < 2(M/16)/h` gives `|T_a| > h - M/8 > 3M/4`. Column peeling then gives
  `F >= M + |T_a| - ||T_a| - M| = min(2|T_a|, 2M) > 3M/2`. (`H_a` is zero-sum and `T(H_a) = M e_a`.)
* **`kappa = 2`.** `|T_a|, |T_b| > h - M/8 > 5M/16`. Half-sum peeling, with
  `g = (sH_a + s'H_b)/2` integral and zero-sum and `T(g) = (M/2)(s e_a + s' e_b)`, gives
  `F >= M + phi(|T_a|) + phi(|T_b|) > 5M/4`.

All cases are covered, because `kappa = 0` is impossible. The exceptional set is used exactly where `c` is
excluded.

*Exact checks.* `r_lemmaS.py`, rerun:
* Lemma 2.2: exhaustive at `q = 7, 11` (571,656 vectors) and random at `q = 43, 47` (including entries
  divisible by `q`). Minimum slack 1.
* Lemma 2.1 items 1-6: 12,000 vectors, no failures.
* At `q = 43`: exhaustive support 4 and 6, minimum `F = 60` and `68`, against `17M/16 = 46.75`. Perturbed
  columns and half-sums have minimum `F = 80`.
* The author's `check_lemmaS6.py`, rerun, agrees.

## 2. Prop. 4.3 (= `lower.md` Prop. 3.3, canonical assignment): correct

* **Slack.** `s_c = (1/2) sum_j w_j(|tau_j| - 1)`, where every `tau_j` is even. Terms with `tau_j != 0` are
  `>= tau`. At most `r` terms equal `-w_j >= -tau - 1/r`. Hence `s_c >= (tau/2)(2B - r) - 1/2`.
* **Identity (3.1).** It is exact for this profile: `M` flipped physical columns, one per row, in label
  `alpha(x)`. Each flip term is `>= -2|c_x|`.
* **Pairs.** `F = M` (rows differ in `M/2` nonconstant columns), and the two flip terms are `+-2`. This
  gives `b - 4`.
* **Columns.** At least `M-2` rows have `T_(alpha(z)) = 0`, giving `+2` each, and at most 2 rows give `-2`.
  This gives `b + 2M - 8`.
* **Half-sums.** The support has `M/2` rows, at most 4 of them in labels `{a, b}`. This gives
  `b + M - 16`.
* **Other ternary.** Lemma S gives `(b/16 - 2)M + b >= kappa_b M`.
* **Amplitude `h`.** `b(hM - M + 1) - 2hM - [hM(b-4)/2 + b] = bM(h/2 - 1) >= 0`.

*Exact check* (`r_slack.py`, rerun): `q = 43, 47, 59, 67`, `b = 33, 40`. The identity, twins, rank and
nonconstancy are verified. The pair bound `b - 4` is attained, and the column bound `b + 2M - 8` is
attained. Every family respects its bound, and the slack inequality has minimum gap 0.495.

## 3. Theorem 5.1 (P3-char): every margin and union bound re-derived; correct

### 3.1 Weights, PNT inputs and Theorem A.1

* `Z = max(3M^A, 21r^2)(log M)^2`. `[Z, 2Z]` holds about `Z/(2 log Z) >> r^2` primes `= 1 mod 4`. This is
  PNT in progressions; no short-interval input is needed. Pigeonhole over `r` intervals of length `Z/r`
  gives `w_j in [tau, tau + 1/r]` and `W <= r tau + 1`.
* `tau <= log 2 + log(21 b^2 M^(max(A,2))) + 2 loglog M`, as stated.
* `p_j > X`, so no level data are needed. `p_j >= 21 r^2 (log M)^2`, so `Pi - 1 = O(1/log M)`, and
  `outside.md` Theorem A.1(ii) applies. I checked its proof: the uniformity of `<v, phi>` for `v` outside `L`
  and the `(8/pi) arcsin` slab bound.
* `M' = (1 + o(1)) M` for the least prime Paley order `M' >= M` follows from PNT for primes `= 3 mod 4`.

### 3.2 Design lemmas

* **Lemma 3.1** (own implementation in `a_mu.py`). Up to `q < 2^23` the multiplicity is 4,
  `2p' <= m_q` and `q < 8p'` hold everywhere, and `Lambda'(d) <= 9 log d + 18.72 omega(d)` holds for all
  `d < 2^20`.
  * I recomputed the Rosser-Schoenfeld ratio `R(k)` for `k >= 20` (`max = 8.1047`, `<= 8.14` for
    `k > 2000`). So the multiplicity is at most 9.
* **Lemma 3.3(1).** On the small-`q` level-1 classes `kappa_x - kappa_y = (sigma_x - sigma_y) + 2(...)` has
  absolute value `<= 2p' - 1 < m_q`, so it vanishes mod `m_q` iff both parts vanish.
  * `Z/m_(q^a) = Z/m_q x Z/q^(a-1)` is compatible with reduction. Level `a*` is injective. Large-`q`
    injections never collide with unit 1.
  * The pair demand is `<= Lambda'(d) + log d`.
* **Gap (minor) in "`<= mu_M log M`".** The text derives `Lambda'(d) <= (9 + 26/loglog d) log d` and then
  replaces `d` by `M`.
  * `(10 + 26/loglog d) log d` is not increasing. It is 315, 124 and 104 at `d = 3, 4, 5`, all above
    `mu_44 log 44 = 111.8` (`a_mu.py` (3)).
  * The correct route: `log x/loglog x` increases for `x >= e^e`, and `omega(d) <= 2` for `d < 30`. So
    `omega(d) <= 1.3841 log M/loglog M` for every `d < M`, `M >= 16`. Hence
    `log m_(xy,1) <= 10 log M + 25.91 log M/loglog M <= mu_M log M`. The claim is true.
  * Exact worst-case (parity-blind) pair demand: `5.7-6.8 log M` for `M = 64 ... 65536` (`a_mu.py` (2)),
    matching `check_design.py`.
* **Lemma 3.3(2).** `q^(a*) < 32 M^2` (because `p'(q) q^(a*-2) < M` and `q < 8p'`). The exact value of
  `log m_struct` is `7.35-8.0 M`, against the bounds `pi(4M) log(32M^2) = 8.8-9.8 M` and `11 M`
  (`a_thm51.py` (1)).
* **Lemma 3.3(3).**
  * Large `q`, level 1: condition on `y_x` for `x != x_0`. `y_(x_0)` is then uniform on at least
    `m/2 - M + 1 > m/4` values, and `c_(x_0) y = t mod m/2` has at most `gcd <= c_min` solutions. So the
    bound is `4 c_min/m` per unit and `16 c_min/m <= 32 c_min/q` for some unit.
  * Lifts: one congruence mod `q`, with the unit forced to reduce compatibly. At most one unit per prime
    can occur, since `iota` has order 4.
  * Monte Carlo (`r_lemma33_mc.py`): the maximum empirical/bound ratio is 0.55.

### 3.3 Events and sufficiency (step 3)

* **Non-pairs.** `|sin((c.delta - arg eta)/2)| = sin dist(c.delta/2 - arg(eta)/2, pi Z)
  >= sin dist(c.delta/2, (pi/4)Z)`, because `arg(eta)/2 in (pi/4)Z`. Non-failure forces that distance
  to be `< pi/8` and `>= arcsin(m^gen e^(-w/2))`. Then `m_(c,eta) <= m^gen`. (G) on `L` follows since
  `m^gen >= 1`.
* **Pairs, unit 1.** `|c.delta/2| <= Delta/2 < pi/8`.
* **Pairs, unit `!= 1`.** The left side is `>= sin(pi/8)`. The right side is `<= 2^(-1/2) e^(-1)`, because
  `w >= W/2` (pair slack `> 0`).
* **Distinctness.** It follows from the pair events.
* **Lemma 5.1a (= `lower.md` 5.3).** Density `<= 2/Delta`, range `n Delta/2`, at most `2n Delta/pi + 2`
  bad intervals of length `<= pi eps`. With `e^(-w/2) = (Delta/C) e^(-s/2)` and `Delta <= C`, this gives
  `K_c = 4(n + pi/C)`. Monte Carlo ratio: 0.31.

### 3.4 The theta-moment union bound (step 4)

* `Pr[Fail_c] <= E[min(1, K_c m^gen e^(-s/2))] <= K_c^theta e^(-theta s/2) m_struct^theta E[m_rand^theta]`
  holds, since `delta` is independent of the design.
* **Per-prime moments.** For tail-dominated laws the exact maximum of `E[q^(theta a)]` is
  `1 + sum_a t_a (q^(theta a) - q^(theta(a-1)))`.
  * Small-`q` lifts: `<= 1 + q^(theta-1)(1 - q^(-theta))/(1 - q^(theta-1)) <= 1 + q^(theta-1)` for
    `theta <= 1/2`, inside the claimed `1 + 2q^(theta-1)`. Worst ratio 0.50 (`a_thm51.py` (2)).
  * Large `q`: `<= 1 + 32 c_min q^(theta-1)`, inside 64. Worst ratio 0.25.
  * Primes dividing `c_(x_0)`: the trivial bound `q^(theta A_q) <= X^theta` applies.
  * `Sigma_theta <= X^theta/theta`: worst ratio 0.24 (`a_thm51.py` (3)).
  * `S(c) = c_min o(M)` uniformly in `c_min`, because `log_2(c) log(c K)/c = O(log K)`.
  * Independence across `q` holds.
* **Families.**
  * Ternary: `3^M`, `s >= tau kappa_b M/2 - 1/2`.
  * Amplitude: `(2h+1)^M`, `s >= 7 tau h M - 1/2`, which needs `(b-4)/4 >= 7`, true for `b >= 32`.
  * Pairs: `M^(mu_M)` against `(b-4)tau/4 >= 14.5 log M`.
  * Off-unit pairs (Markov with `W >= 64 M log M`).
  * Columns and half-sums are inside the ternary family with slack `>= M`.
  * Since `v(c) != 0` for `c != 0`, the union covers all characters.
* **Characters with large common factors** (for example `c = h(e_x - e_y)` with `q | h`). They collide at
  all lifts, but this is exactly what the "`q | c_(x_0)`" clause pays for. Their slack `7 tau h M` grows
  with `h`.
* **The design needs uniform `Theta(M)` slack for all non-pairs.**
  * A structured design collides on additive quadruples at most small `q`: exact
    `log m ~ 1.4-1.7 M`, `r_struct.py`.
  * So the margins of `construct.md` Theorem 1 (`G >= 2n`) would not suffice here.
  * This is why `b >= 33` with Lemma S is used. The design choice is consistent.
* **Quantifiers.** `M_0 = M_0(A, b, C)`, and `C` enters only through `log(1/C)` factors. That is fine.
* **Size of `M_0` with the proved constants** (`a_thm51.py` (4), `a_mu.py` (4); `b = 33`, `C = 1`):
  * the ternary family closes only for `log M >~ 450`;
  * the pair family, with the proved `mu_M = 10 + 26/loglog M`, needs `mu_M < 12.5`, i.e.
    `log M >~ 3.5 * 10^4` (it closes accidentally near `log M = 10` and fails for
    `50 <= log M <= 3 * 10^4`).

  Caution: the earlier `r_thm51.py` (d) prints "pairs: 4" and "amplitude2: 4". Those are the *first*
  `log M` at which the family drops below 1/4, not a threshold beyond which it stays there, because the
  pair family is non-monotone in `M`. `a_mu.py` (4) and `a_thm51.py` (4) supersede it.

  The theorem is asymptotic, so this is not an error. But the "explicit" bound
  `W <= b(M-1)(max(A,2) log M + 2 loglog M + log(42b^2)) + 1` at prime Paley orders is valid only for
  `M >= M_0`, and with the proved constants `M_0 > e^(3*10^4)`. Say so next to the formula.

### 3.5 Completion and restriction

* **Step 5.** `rank S = M` (twins). Theorem A.1(ii) needs character-admissibility on `L`, which step 3
  gives, and `Pi < 1.75`. Correct.
* **Lemma 5.1b (= `lower.md` Lemma 5.2), minor gap.** Suppose a deleted (constant) column has
  `p_l <= X'`. Then the restricted model has `p_l` not dividing `N'`, so it needs point residues at `p_l`
  and the full character condition (1.1) there. The original model carried only level residues at `p_l`,
  with the same-level **pair** clause; non-pair characters were never constrained at `p_l`. So "inherited"
  fails there.
  * The same gap is in `lower.md` Lemma 5.2. There, multiplying by `beta` of norm `N_new/N` is not even
    defined modulo `p_l`.
  * It is harmless for Theorems 5.1 and 6.4, where every `p_j > X`.
  * Fix: add "all `p_j > X`" (or "no deleted prime is `<= X'`") to the hypotheses.

### 3.6 Corollary 5.2

`log R = W/2` and `loglog R = (1 + o(1)) log M` give `f(R) >= (2/K - o(1)) log R/loglog R`. Correct. The
bracket `[1/33, 2/3]` for `A <= 2` is correct.

## 4. Theorem 6.4 (P3-all): correct, including the rate

* **Lemma 6.2 re-derived.** `c'_x = eps_x(v_(f(x)) - v_(s(x))) = c_x mod G` and `||c'||_1 <= ||v||_1 < G/2`
  force `sum c' = 0`. With binary rows, both `v` and `c'^T a` have entries `< G/2` in absolute value, so
  they are equal.
* **Lemma 6.3 re-derived.**
  * The parametrisation (2) is a bijection, with `r - M + 1` free coordinates.
  * Reduction `K'_Q -> K'_(q^a)` is onto (twins make `A` onto over every `Z/n`).
  * The image of `<v, .>` has index `d`, and `d > 1` forces `v in L + dZ^r`. The `psi o A` factorisation
    is valid over `Z/d`, since `A` is onto.
  * Hence `Pr <= 2d/m <= 4||v||_1/m` per unit.
  * Exact check (`r_fibre.py`, rerun): `q_0 = 7, 11`, about 6,400 vectors each. Maximum `d/(2||v||_1)` is
    0.25, and maximum `G/(2||v||_1)` is 0.47.
* **Step 2.**
  * The fibre shifts all point classes by a common `2t`, so the character events are unchanged. `tau` is
    larger, and every family bound in §3.4 is monotone in `tau`.
  * (1.2) for `v in L`: non-pairs follow from `m^gen`. Pairs with `s = 0` follow from `Fail_xy`, which
    gives `|sin| >= m e^(-w/2)`, the `nu_0 = 1` form. Pairs with `s != 0` follow from `Off_xy`.
  * The text cites only `m^gen` here, which pairs do not have. This is a wording slip.
* **Steps 3-4.** The slab measure is `(2/pi) arcsin(y) <= y`, with the union over `s` bounded by
  `min(1, 4 m_v e^(-w/2))`.
  * Per-prime bounds: `(E_max - 1)/(58 kk q^(theta-1))` is at most 0.41 for `q > 24kk`. For `q <= 24kk`
    the stated bound is weaker than the trivial `X^theta = e^beta` (`a_thm64.py` (i)).
  * Summed over all primes `<= X <= 10^7`: at most 0.047 of the claimed `kk`-bracket (`a_thm64.py` (ii)).
  * Count `<= (2r)^kk`. The exponent per `kk` is exactly `-4`, so `E[Haar(bad)] <= 4/(e^4 - 1) + tiny < 1/4`.
  * Markov plus step 2 give a good design with positive probability.
* **Rate.** `tau_* = (2A + o(1)) log^2 M/loglog M` is correct, but the `o(1)` is
  `~ (2 logloglog M + 4)/loglog M`.
  * `beta > 0` (so `theta > 0`) only from `log M = 4612` (`A = 1`) and `5049` (`A = 2`).
  * `tau_*/(log^2 M/loglog M)` is 34.0, 10.1, 6.7, 4.1 and 3.0 at `log M = 10^4, 10^5, 10^6, 10^9, 10^15`
    (`A = 1`, limit 2; `a_thm64.py` (iv)).
  * The statement is purely asymptotic, which is fine, but it should be described that way.

## 5. Prop. 6.1 / 6.5 (height condition, sharpness): one major error, repairable

* **Prop. 6.1 is correct as stated, per unit.** For actual circles, `q^a | D_s` with
  `D_s in {Y, X-Y, X, X+Y}`, and `|D_s| = nu_s^(-1) |A_v| |sin(arg A_v - s pi/4)|`. I verified all four
  identities.
* **The error.** The aggregated `m_v = prod_q q^(max_s a_(q,s))` is only an upper bound, used for the
  completion sum. The inequality `m_v <= sqrt2 e^(w/2)` is **not** implied by membership: the four `D_s`
  are pairwise almost coprime, so `m_v` is about `|XY(X^2 - Y^2)| ~ e^(2w)` for actual data. Exact data
  (`a_height.py`):
  * 8976 of the 8977 primes `p = 1 mod 4` below `2 * 10^5` violate the aggregated inequality (moduli up to
    `10^12`, and also up to 3000). Only `p = 5` does not.
  * All 8000 two-prime monomials violate it.
  * Zero cases violate Prop. 6.1.

  Yet the Verdict item 5 calls `m_v <= sqrt2 e^(w(v)/2)` "the necessary height condition (Prop. 6.1)", the
  proof of Prop. 6.5 ends with "`m_(e_j) > sqrt2 e^(w_j/2)`. So the height condition fails", and
  `check_p3all.py` (docstring and Parts A/B) measures the aggregated quantity. As written, the proof of
  Prop. 6.5 establishes nothing.
* **Repair (checked).** Fix a unit `s`. For `q in Q_0` (`q = +-3 mod 8`, so `m_q = 4 mod 8`) the labels of
  non-twin columns are uniform on their parity class. So `e_j` collides with unit `i^s` at `q` with
  probability `2/m_q` if `parity(iota s) = nu_j(q)`, and 0 otherwise.
  * Put `lambda_(j,s) = sum 2/m_q` over those `q`. Then `lambda_(j,0) + lambda_(j,1) = lambda/2`, so some
    `s_j` has `lambda_(j,s_j) >= lambda/4`, a positive constant.
  * The Poisson lower bound `P(>= kk hits) >= (1 - o(1)) e^(-2 lambda') lambda'^kk/kk!` still equals
    `exp(-(1 - eta_0 + o(1)) log M)`, because `kk log kk ~ (1 - eta_0) log M` dominates and `lambda'` only
    enters `o(log M)`.
  * Independence over columns is unchanged. So w.h.p. some `(j, s)` has `m_(e_j,s) > nu_s^(-1) e^(w_j/2)`,
    a genuine violation of Prop. 6.1.
  * The threshold `(2 + o(1)) log^2 M/loglog M` is unchanged.
* **Quantitatively** (`a_p65.py`, exact DP over primes `<= 3M`, `r = 33M`): the `tau` at which the expected
  number of violators drops to 1, in units of `log^2 M/loglog M`:

  | `M` | per-unit (correct) | aggregated (as written) |
  |---|---|---|
  | `10^3` | 3.48 | 6.42 |
  | `10^4` | 3.36 | 6.16 |
  | `10^5` | 3.25 | 5.91 |
  | `10^6` | 3.17 | 5.70 |

  Both tend slowly to 2. So the "obstruction" quoted from aggregated data is about 1.8 times too strong at
  every accessible `M`. Monte Carlo maxima of the per-unit demand agree with the DP tails.
* **Understated.** Taking `Q_0 = {q = +-3 mod 8 : X^(1-eta_0) <= q <= X}` instead gives failure for
  `tau <= (2A - eps) log^2 M/loglog M`. So Theorem 6.4 is sharp **for every `A`** for the uniform fibre,
  not only at `X = 3M`.
* **Scope (minor).** Prop. 6.5 concerns the uniform random fibre of Def. 6.3. Its witnesses are non-twin
  columns, whose labels are free in their parity class (Lemma 6.3(2)). Restricting those labels away from
  the 4-torsion removes all non-twin single-prime witnesses at no cost.
  * The `2M` twin columns would still carry Poisson tails, and removing them needs a coordinated choice
    that nobody has analysed.
  * So "with per-modulus independent residue choices, the height condition fails w.h.p." (Verdict item 5)
    is broader than what is proved. "For the uniform random fibre" is accurate. §6.4's own sentence ("a
    property of the method, not of the class") is accurate.

## 6. Appendix A (Sylvester multi-block route): sketch only; citation used correctly

* **Status.** The appendix says "a referee pass is still owed". Verdict item 7 ("also works for III-small
  and for P3-char") and Appendix B (rows "union bounds converge with `t` absolute" and "Green-Sanders
  statement") read as established results. They should say *sketch*.
* **Citation use: correct.**
  * Green-Sanders (Ann. Math. 168 (2008) 1025-1054; GAFA 18 (2008) 144-162) is a theorem about **Boolean**
    functions, with `||f||_A = sum |fhat|` and `fhat = E f chi`, so `||1_coset||_A = 1`, checked in
    `a_appA.py` (2). The GAFA abstract gives `+-1` combinations of `2^(2^(O(M^4)))` **subgroup** indicators.
  * Lemma A.2 reduces integer-valued `f` with `||f||_A <= 3` to Boolean level sets, through the Banach-algebra
    product. The maximum of the product is exactly 225 (at `v = +-1`), and the coefficient count is
    `sum_(v != 0) |v| = 12` (`a_appA.py` (1)). Cosets versus subgroups costs a factor 2 in `F_2^n`.
  * Sanders arXiv:1610.07092 exists and gives `exp(M^(4+o(1)))` for the general (integer-valued /
    idempotent) case, but it is not needed.
  * The later improvement `exp(M^(3+o(1)))` in `F_2^n` (Sanders, "Boolean functions with small spectral
    norm, revisited") is not cited. That is harmless.
* **The "Counting/Assembly" union bound is under-specified.**
  * With the stated `Pr[f_i structured] <= |S_L|/Mult(mu_c)`, `Mult >= binom(M, ceil(s/4))` and a union over
    **characters**, the bound at the stated threshold `s = 4L log M`, `t = 5` **diverges**.
  * A union over **multiset classes** (`Mult(mu)` characters each, so the total is
    `|S_L|^t Mult^(1-t)`) converges (`a_appA.py` (4)): for `m = 40 ... 2000`, `L = 10, 1000` the natural
    logs are `+9.5e3 ... +1.6e9` (naive) against `-1.0e4 ... -4.1e9` (multiset).
  * The sketch should state the multiset union.
* **"`m_c <= e^(2.6M)` (all accessible moduli)".** This holds for III-small moduli (`m_Q < 2M`, so
  `q^a <= 2.5M` for `q >= 5`). For P3-char at `X = 3M^A` the accessible product is `e^(psi(3M^A))`. That is
  about `e^(3.1M)` for `A = 1`, and far more for `A > 1`, and the claimed absorption fails.
  * P3-char would need the theta-moment machinery of Theorem 5.1 (random design, `m_struct m_rand`).
  * The claim "works for P3-char" is therefore unsupported as written.
* **Verified pieces.**
  * Lemma A.1's counting. `#{|T| = n} > M/(2n)` follows from the ternary defect identity.
  * The coplanarity probability of 4 random distinct points is exactly `1/(M-3)` (exhaustive for
    `m = 3..6`, `a_appA.py` (3)).
  * `t >= 5` is needed for support 4.
  * The margin `Tw >= max(-M, E_t - 1 - 2n)` is consistent only if `Tw` counts each twin column as
    `|.| - 1`; the text should define it.

## 7. Smaller points

1. Verdict item 1: "This contains the task's class III-small" has the inclusion reversed. The table is
   right: `F^char_(3M) subset III-small`, so fakes in P3-char are fakes in III-small.
2. Prop. 6.5: there are `(b-2)M - b` non-twin columns, not "at least `(b-2)M`". This is immaterial.
3. Prop. 7.1 uses (1.1) for `c`, which requires `v(c) != 0`. That holds on saturated profiles (§4) but
   should be stated.
4. §8's evidence row for `check_p3all.py` describes aggregated single-column moduli as height-relevant;
   relabel it.
5. `check_realise.py` uses `X = 2M^2`, not `3M^A`. The text says so, and the check is consistent.
6. `check_design.py` Part 1 now runs to `q < 2^25`, and the text says `2^23`. Harmless.

## 8. Corrected statements

* **Theorem 5.1:** unchanged. In the pair step, replace the `d`-wise chain by
  `omega(d) <= 1.3841 log M/loglog M` for all `d < M` (`M >= 16`), which gives
  `log m_(xy,1) <= (10 + 26/loglog M) log M`. Add: "`M_0` is ineffective in practice; with the proved
  constants the pair family needs `loglog M > 10.4`."
* **Lemma 5.1b / `lower.md` Lemma 5.2:** add the hypothesis "every `p_j > X`" (or "no deleted column has
  `p_j <= X'`").
* **Theorem 6.4:** unchanged. Say "for `log M >= 4612` the parameters are defined; the `o(1)` is
  `O(logloglog M/loglog M)`".
* **Prop. 6.5:** *Fix `eps > 0` and `A >= 1`, and use the construction of Theorem 6.4 with the uniform
  random fibre, `X = 3M^A`, and `tau <= (2A - eps) log^2 M/loglog M`. Then with probability `-> 1` some
  single-prime monomial `e_j` and some unit `s` have `m_(e_j,s) > nu_s^(-1) e^(w_j/2)`. So Prop. 6.1 fails,
  and no `delta` or `phi` completes the data.* The proof is per unit, as in §5.
* **Verdict item 5:** "the necessary height condition `m_(v,s) <= nu_s^(-1) e^(w(v)/2)` for each unit `s`
  (Prop. 6.1) fails w.h.p. for the uniform random fibre below `tau = (2A - eps) log^2 M/loglog M`".
* **Verdict item 7 and Appendix B:** "Sketch (not refereed): the Sylvester route plausibly works for
  III-small. For P3-char it would additionally need the theta-moment design of Theorem 5.1."

## 9. Script index (this pass)

| script | checks | result |
|---|---|---|
| `a_mu.py` | Lemma 3.1 (multiplicity, `Lambda'` bound), exact pair demand, non-monotonicity of the `d`-wise `mu` chain, pair-family threshold with the proved `mu_M` | all claims hold; chain gap; pair family needs `log M >~ 3.5e4` |
| `a_height.py` | per-unit versus aggregated height for actual Gaussian primes and 2-prime monomials | per-unit: 0 violations; aggregated: 8976/8977 and 8000/8000 violations |
| `a_p65.py` | exact DP and Monte Carlo heavy tails, per-unit versus aggregated | per-unit threshold 3.48 down to 3.17; aggregated 6.42 down to 5.70 (`M = 10^3 ... 10^6`) |
| `a_thm51.py` | `m_struct`, per-prime theta-moments, `Sigma_theta`, family exponents with proved constants | all bounds hold (worst ratios 0.50, 0.25, 0.24); `M_0` astronomical |
| `a_thm64.py` | step-4 moment bounds, summed bracket, `(2r)^kk` count, `beta`/`tau_*` asymptotics | worst ratios 0.41 and 0.047; `beta > 0` from `log M = 4612` |
| `a_appA.py` | 225 / 12, coset norm, coplanarity `1/(M-3)`, Assembly union bound | as in §6 |
| `r_*.py` (rerun) | Lemma S, Prop. 4.3, Lemma 6.2/6.3, Lemma 3.3(3) Monte Carlo, Thm 5.1 / 6.4 numerics, Green-Sanders level sets | all pass; `r_thm51.py` (d) reports first crossings only (see §3.4) |
| `rerun2_check_*.txt` | the author's five scripts | all PASSED |
