# Referee report (pass 2, independent, adversarial) on round17/hlight.md: H_light for sdc_q

Scripts: `round17/referee_hlight_checks/pass2/`. Each `x.py` / `x.c` writes `x_out.txt`. Nothing is imported from the author's
`hlight_checks/`. The library `plib.py` builds `sdc_q` directly from the record's definition (round5/exact.md Thm 4.1 and
fake.md §1): Paley-I skew-Hadamard `K`, normalised by `DKD` so that `K(x, base) = 1`, `S_a = {x : K(x,a) = -1}`, and
`S'_a = S_a xor {inf, a}`. It does **not** use the author's `P_a = a - D`, `P'_a = a + D`.
An earlier referee pass on the same file is kept, unchanged, as `round17/referee_hlight_pass1.md`. I did not rely on it.
Where I agree with it, I agree because my own checks led to the same place.

## 0. Verdict

**The result survives.** Theorem 2.8 is correct as a theorem, with all quantifiers, for every prime `q = 3 mod 4`, `q >= 7`:
every integral `lambda` with `sum lambda = 0` that is not an integer multiple of a pair has `sigma_lambda >= LB(q)`.
Three things follow from it:
* the elementary `c_0 = 1/6` for every `q`;
* the Weil-based `c_0 >= 0.2349` at every regime prime `q >= 3251`;
* sharpness up to `O(q^{-1/2})`.

All three are correct. One proviso applies at `q = 7`, where the "in particular" does not follow from the displayed `LB`; this is issue M1. I re-derived every lemma by hand. Every number I recomputed agrees with the write-up, using my own code built from the record's definition.

There are no fatal or major issues and six minor ones (§4). They are a `q = 7` slip, an overstatement about base 4-terms, a solver result labelled "exact", the status of the Weil citation, an incomplete threshold list in §4, and the rounding behind the constants in the table. None of them changes a proved conclusion.

CIRCLE(q) is **not** decided, and the write-up correctly says so. There is no growth or uniformity claim, and none is made.

## 1. Definitions checked against the record

* **Pattern family.** In the definition, `DKD` makes `K(inf, a) = -1`. So `inf ∈ S_a`,
  `S_a = {inf} ∪ {x : chi(a-x) = -1}`, its complement is the author's `P_a = {a} ∪ {x : chi(a-x) = 1}`, and
  `S'_a = {a} ∪ {x : chi(a-x) = -1}` is the author's `P'_a`. All `2q` cuts are distinct and balanced (`|S| = (q+1)/2`). For
  `sum lambda = 0` we have `|lambda(S^c)| = |lambda(S)|`, so `sigma` is unchanged by the complement.
  * *Note:* my first build omitted the `DKD` normalisation. It gave unbalanced sets of sizes `(q-1)/2, (q+3)/2` and a formula mismatch. Pattern sets must be built with the normalisation. The write-up's sets are the normalised ones.
* **Lemma 1.1 formula vs the definition** (`c1_definition_out.txt`). There are 0 mismatches over 300 random integral `lambda` per `q`, for `q = 7..59`. The test includes `lambda_inf != 0`, `|F|` up to 2, and entries up to 5. Both orientation conventions of the Paley core give identical `sigma`. The pair slack is `sigma_xy = 1` for all pairs, including pairs with `inf`.
* **Relevance / integrality** (the domain of H_light). I give my own proof that the integral `lambda` are exactly the relevant ones in `L_T`.
  * If `lambda ∈ L_T`, then `lambda(P_a) + lambda(P'_a) = lambda_a - lambda_inf ∈ Z`. So `lambda = mu + t·1` with `mu` integral.
  * Since `|P_a| = (q+1)/2`, `lambda(P_a) ∈ Z` forces `t ∈ (2/(q+1))Z`.
  * With `k = e_1`, `lambda.k = mu_1 + t`. So `lambda` is relevant iff `t ∈ Z` iff `lambda` is integral.

  An exact Smith normal form check confirms `L_T/Z^n_0 = Z/((q+1)/2)` for `q = 7, 11, 19, 23`. This matches fake.md §1, `Z/(n/2)`.
  So Theorem 2.8, stated for all integral non-pair `lambda`, covers the whole domain of H_light in barrier Lemma 3.4: relevant, non-pair, and light.
* **Slack.** `sigma_lambda = sum_S (|lambda(S)| - 1/2)` over the `2q` weight-1 patterns, which equals `N(lambda) - q`. This agrees with barrier.md §1 and fake.md §1, with `mu = 1` and `r = 2q`.

## 2. Step-by-step re-derivation (by hand)

* **L1.1.** `1_QR = (1 - delta_0 + chi)/2` on differences gives `lambda(P_a) = (h + c)/2` and `lambda(P'_a) = (h - c)/2`. Also `|U| + |V| = max(|U+V|, |U-V|)`. Correct.
* **L1.2.**
  * (1) For `a != 0`, `sum_t chi(t(a-t)) = sum_{u != -1} chi(u) = -chi(-1) = 1`, so `chi*chi = 1 - q delta_0` and `c*chi = F - q f`.
  * (2) The inner sum equals `chi(-1)(chi*chi)(t'-t) = q delta - 1`. This gives `sum c^2 = q||f||^2 - F^2`, `sum h^2 = ||f||^2 + (q+2)F^2`, and, using `lambda_inf = -F`, `sum_S lambda(S)^2 = (q+1)||lambda||_2^2/2`.
  * (3) and (4) are immediate.

  Correct. The Remark is also correct: the matrix `S` is skew and `S S^T = qI`, which I checked by hand row by row, and `c(a) = mu_a - mu_inf`.
* **L2.1.**
  * (1) `N >= |sum h| = (q+1)|F|`.
  * (2) `q f(b) = F - sum_{a != b} c(a) chi(b-a)` gives `||c||_1 >= 2q - 1 + |c(b)|`, and `N >= ||c||_1`.
  * The remaining case `|F| <= 1`, `|f| <= 1` is ternary, because `lambda_inf = -F`.

  Correct.
* **L2.2.** `|t| = (2/n)t^2 + (2/n)|t|(n/2 - |t|)`. For integers `1 <= |t| <= n/2 - 1`, `|t|(n/2-|t|) >= n/2 - 1`. With `||lambda||_2^2 = n` (ternary) this gives `N = q + 1 + ...`. Correct.
* **L2.3.**
  * Euler's criterion `chi(u) = u^m` holds at `u = 0` because `m >= 1`.
  * The coefficient of `X^{m-j}` in `sum f(t)(X-t)^m` is `binom(m,j)(-1)^j M_j`, and `binom(m,j) != 0 (mod q)` because `m < q`.
  * Vandermonde on `T` with `|T| <= m + 1 = (q+1)/2` gives `j_0 <= |T| - 1 <= m`, so `P != 0` with `deg P = m - j_0 < q`.
  * Every integer zero of `f*chi` is a root of `P` mod `q`.

  Correct, with the hypotheses exactly as stated. Equality at pairs: `j_0 = 1` and Jacobsthal gives `(q-3)/2` zeros.
* **P2.4 (`|F| = 1`).**
  * For `n >= 6`: `U + V = h ∈ {1, 2}` while sums of entries in `{0, ±n/2}` lie in `{0, ±n/2, ±n}` with `n/2 >= 3`. Hence `nu >= q + 1 - n/2`.
  * `phi(n) = q + 2 - n/2 - 2(q+1)/n` is concave, `phi(6) = 2(q-2)/3`, `phi(q+1) = (q-1)/2`, and `n <= q + 1`.
  * For `n = 4`: only the 2 points with `h = 2` can be trivial, so `nu >= q - 3`.

  Correct.
* **P2.5 (`F = 0`).**
  * On `T`, `U + V = ±1`, so `nu >= n` there.
  * Off `T`, `U = -V = c/2`, `c` is even and `|c| <= n`.
  * `N_0 <= (q-3)/2` by L2.3 (`j_0 >= 1` because `M_0 = F = 0`).
  * `sum_{a ∉ T} c^2 <= qn - n`, so `4M + n^2E <= qn - n`. Hence `M(1 - 4/n^2) >= (q+3)/2 - n - (q-1)/n`.
  * `sigma >= n - 1 + (2n/(n+2))[(q+3)/2 - n - (q-1)/n] = g(n)`, with `g(4) = q/3` and `g(6) = (q-3)/2`.
  * `(n+2)(g - (q-3)/2) = [q(n-6) + 11n - 2n^2 + 6]/2`, which is `>= 0` for `q >= 2n+1`.

  I recomputed all of these, and all are correct. When the bracket is negative, `g(n) <= n - 1`, so (1) already implies (2).
* **P2.6 (dual pairs).** `c = (q/2)(delta_1 - delta_0) - f + e`, `c(0) = -(q+1)/2`, `c(1) = (q-1)/2`, and `||f||_1 = 1 + (q-1)/2`
  (Jacobsthal: `A - D = -1`, `A + D = q - 2`). So `sigma = (q-1)/2`. The character is ternary and not a pair for `q >= 7`. Correct.
* **L2.7.**
  * With `x_i = f_i eps_i`, I solved `|sum x_i| = a + b e_2(x) + c x_1x_2x_3x_4` from the three orbit values `4, 2, 0`. This gives `a = 3/2`, `b = 1/2`, `c = -1/2`.
  * `e_2(x) = sum_I (prod_{I} f) eps_I`, and `prod_T f = +1`, so `s_T = -1`. The coefficients are as claimed.
  * `|S_I| <= 1 + 2`: the two points of `T \ I` are removed, and the points of `I` contribute 0.
  * `S_T` is the full sum, because `g(a) = 0` on `T`.
  * Constants: `(3/2)(q-4) - 9 - (3/2)sqrt q + 4 = 1.5q - 11 - 1.5 sqrt q`.

  Correct. Weil is applied to `g = prod(X - t_i)`, which is monic, squarefree, of degree 4 and not a square, with `chi` of order 2. The `(deg g - 1)` in the write-up equals the true `d - 1`, where `d` is the number of distinct roots.
* **Thm 2.8 / Cor 2.9.** The case assembly is correct. Multiples of pairs: `N(l·pair) = |l|(q+1)`. For regime `q`, `LB(q)/(2q) = 1/4 - 5.5/q - 0.75/sqrt q` is increasing. Its exact value at `q = 3251` is `1528.97 = 0.235154 r`, which is `>= 0.2349 r`.

## 3. Independent exact computations

| check | result | file |
| --- | --- | --- |
| `sdc_q` from the record's definition; formula vs definition; pairs; `L_T` | balanced, `2q` distinct cuts; 0 mismatches (300 random `lambda` each `q = 7..59`); `sigma_xy = 1`; SNF `L_T/Z^n_0 = Z/((q+1)/2)` (`q = 7, 11, 19, 23`) | `c1_definition_out.txt` |
| exact global minimum, `q = 7` (all `lambda`, `|lambda_x| <= 2`), `q = 11` (all ternary); with Lemma 2.1 these are global | `3`, `5` (`= (q-1)/2`); pair multiples min 1 | `c2_exhaust_out.txt` |
| `q = 7`, `F = 0`, `n >= 6` (the M1 case) | min `sigma = 5` (`>= 7/3`) | `c2_exhaust_out.txt` |
| `q = 11`, every class `(F, n)`, exhaustive ternary | every class minimum `>=` its proved bound. Tight cases: `F = 0`, `n = 6` (`5 = n - 1`; the dual pair), `F = 0`, `n = 10` (`9 = n - 1`), `F = 1`, `n = 4` (`5 = (q-1)/2`) | `c2_exhaust_out.txt` |
| 4-term `F = 0`, own C (`r4b.c`, normal form `d_0 + d_1 - d_c - d_e`), `q = 7..1019` (16 primes) and `3251, 3259, 3271, 5479` | minima `3, 5, 7, 9, 13, 17, ..., 43 (103), 87 (199), 229 (499), 481 (1019)`, **`1571` (3251), `1575` (3259), `1581` (3271), `2667` (5479)**, identical to the write-up; 0 violations of `sigma >= q/2 - 11 - 1.5 sqrt q` (exact integer test); `max |S_T| <= 2 sqrt q + 1` everywhere (`149` vs `149.04` at 5479) | `r4b_small_out.txt`, `r4b_regime_out.txt` |
| normal form vs brute force (no symmetry used) | `7`, `9` at `q = 19, 23` (equal) | `c3_misc_out.txt` |
| Lemma 2.3 | 0 violations: 6884 random `f` (`q = 7..131`, `|T| <= (q+1)/2`, values including multiples of `q`); 51194 ternary `f` at `q = 11` (exhaustive) | `c3_misc_out.txt` |
| dual pairs | `sigma = (q-1)/2`, `||lambda||_1 = (q+1)/2`: from the definition at 12 primes `7..83`; by formula at `3251, 3259, 5479` | `c3_misc_out.txt` |
| base 4-terms `e_x + e_y - e_z - e_inf` | values exactly `{(q-1)/2, (q+3)/2}` for `q = 7..103` (M2) | `c3_misc_out.txt` |
| regime arithmetic (interval arithmetic) | first regime prime `3251` (`q - 316 ln 8q = 38.49`; `3203` fails by `4.81`); `w = 514.1998`, `j = 257`, `m' = 2057 < q`, `2m' < q` false (needs R16-J2 `m < N`); `q/2 - 11 - 1.5 sqrt q = 1528.97 = 0.235154 r` | `c4_regime_out.txt` |

All of these are exact integer computations, apart from the interval-certified regime arithmetic. The tables are evidence for nothing beyond what they state. The proof does not use them.

## 4. Issues

**M1 (minor; the `q = 7` slip).** At `q = 7`, `LB(7) = min(max(7/3, ...), (7-3)/2) = 2 < 7/3`. So "in particular H_light(1/6)
for every q" does not follow from the displayed `LB` at `q = 7`. The §3 table entry "proved LB/r = .1667" at `q = 7` disagrees with the displayed formula, which gives `2/14 = .1429`.
* The claim itself is true. At `q = 7`, `n >= 6 > (q+1)/2`, so Prop 2.5(1) gives `sigma >= n - 1 >= 5`.
* Exactly, `min sigma = 3` at `q = 7` (`c2_exhaust_out.txt`).
* Fix: state `LB(q)` for `q >= 11`, where `(q-3)/2 >= q/3` iff `q >= 9`, and add "`sigma >= 3` at `q = 7`". Alternatively, write the `n >= 6` bound as `max(n-1, g(n))` in place of `(q-3)/2`.

**M2 (minor; overstatement in §0).** §0 says the base 4-terms `e_x + e_y - e_z - e_inf` "have `sigma = (q-1)/2` exactly". As a
statement about the family this is false. The family takes exactly two values, `(q-1)/2` and `(q+3)/2`, for every
`q = 7..103` (`c3_misc_out.txt`). Prop 2.4's wording, "attained ... for a suitable `z`", is the correct one. The
sharpness claim does not depend on this, because the dual pairs give it.

**M3 (minor; a solver result labelled "exact").** §0 and §3 call the HiGHS ternary MILP minima for `q <= 43` "exact global
minima". A floating-point MILP with default tolerances is not an exact certificate. The conclusion is nevertheless right, and the MILP is not needed:
* For `19 <= q <= 43` the exact C enumeration of 4-terms gives a minimum below `(q-3)/2`. By Theorem 2.8's case analysis that is the global minimum.
* For `q = 7, 11` my exhaustive enumeration (§3), together with Lemma 2.1, gives `3` and `5` exactly.

Relabel the MILP as "solver-optimal, consistent with" these results.

**M4 (minor; citation status; the hypotheses are right).** Weil's bound is used with its true hypotheses: `chi` nontrivial of order 2, `g` monic,
not a square, and `d` = number of distinct roots (= `deg g = 4` here).
* [L] Citing papers located by search this round identify the source as Lidl–Niederreiter, *Finite Fields*, **Theorem 5.41**. In that statement `psi` has order `m > 1`, `f` is monic and not an `m`-th power, `d` is the number of distinct roots in a splitting field, and the bound is `(d-1) q^{1/2}`. The papers are [arXiv:1902.04736](https://arxiv.org/pdf/1902.04736) and [arXiv:2012.01483](https://arxiv.org/pdf/2012.01483). I did not see the book's own text, so the label stays [L] at citation level.
* [M] For a squarefree quartic, Hasse's bound for the genus-1 curve `y^2 = g` gives `|S_T + 1| <= 2 sqrt q`. That would sharpen Lemma 2.7 to `q/2 - 10.5 - sqrt q`. My enumeration finds `max |S_T| <= 2 sqrt q + 1` at every checked `q`, which is evidence only. This is optional.
* Independently of any citation, at the four regime primes computed, the bound `sigma >= q/2 - 11 - 1.5 sqrt q` is verified for **every**
  4-term character by exact enumeration (0 violations).

**M5 (minor; consequence statement in §4).**
* (a) §4 item 2 lists `tau >= max(tau_D, tau_2, 4)`. Barrier Thm 2.3 also requires `sigma*(tau) <= r + 1`,
  `||lambda||_1 Delta <= 1/2` for light `lambda`, and `n >= 11`. §0 says "thresholds of barrier Thm 2.3", which is right; item 2 should say the same. All of these are `O(log q)`, or automatic.
* (b) "`inf_theta kappa_H <= ...`" should read "some point of `C_1` has `kappa_H <= ...`". This is what Thm 3.2(i) takes as input. The two are equivalent only because `kappa_H` is continuous on the compact circle, so the infimum is attained, and that step is not said.

**M6 (minor; rounding).** The "proved `LB/r`" row and the regime constant use `ceil(sqrt q)` in place of `sqrt q`. Examples:
* `q = 199`: the table gives `.1667`, while the displayed formula gives `67.34/398 = .1692`;
* `q = 3251`: the table uses `1527.5`, while the true value is `1528.97`.

These values are conservative, so the statements stay true. But the table is not "LB(q)" as defined, and the write-up should say so.

## 5. Constants, regime, kappa_H/kappa_V, consequences (§4 of the write-up)

* **Regime.** `n = q + 1`, `d = n - 1 = q`, `r = 2q`, `4r = 8q`. The first prime `q = 3 mod 4` with `q >= 316 ln 8q` is `3251` (interval-certified). Its `w = sqrt(8 q ln 8q) = 514.1998`, `j = 257` and `m' = 2057 < q`.
  * `2m' < q` fails at `3251`, so Thm 3.2 / 3.5 need the R16-J2/J3 relaxation `m < N` (N = q odd), which is used here.
  * `m' < q` for every `q >= 3251`, since `q - 1 - 4w` is increasing and equals `1193.2` at `3251`.
* **Cap threshold.** The write-up carries `eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` and `W = 2q tau = r tau` over from barrier Thm 3.5 unchanged. This is right: Theorem 2.8 touches only the Diophantine side (Lemma 3.4), and neither `eps_0` nor `m'` nor the Thm 3.2 factor. Thm 3.2(i) multiplies the threshold by `sqrt(m'+1) e^{pi^2/(4q)}` and returns `eps_0/2`. Barrier Thm 2.3 adds `m' theta e^{-W/4}/sqrt d`, which stays `<= eps_0/2` for small `theta`. This gives the fakes' `kappa_H <= eps_0`.
* **`kappa_H` vs `kappa_V`.** The cap test, Thm 3.2 and Thm 2.3's conclusion are all about `kappa_H`, and the write-up uses `kappa_H` throughout. It makes no claim that would need `kappa_V`, so the one-sided inequality `kappa_H <= kappa_V` is never invoked in the wrong direction. Correct.
* **`tau_D`.** Barrier Lemma 3.4's formula with `r = 2q`, `n = q + 1` and `2/c_0 = 12` is reproduced exactly. With `c_0 = 0.2349` the factor is `8.5143`, which the write-up rounds up to `8.52`. That is fine.
  * I also checked the Lemma 3.4 measure count. It needs `tau c_0 q >= (q+1) ln(8q+5) + (q-1) ln(4q+2) + O(ln(1/C))`. That holds once `tau >= (1/c_0)(2 ln(8q+5) + ...)`, so the stated `tau_D` is conservative by a factor of about 2.
  * At `q = 3251` the `C`-free parts are `244` (`c_0 = 1/6`) and `173` (`c_0 = 0.2349`). Both are `O(log q)`, so the barrier scale `W = r tau` is preserved.
* **Lemma 3.4's other inputs.** Every input other than H_light is unchanged. Pairs give `sigma_xy = 1`, which I confirmed exactly. If `lambda` is constant on `F_q`, it is heavy. Relevant light `lambda` are integral, which I proved in §1.
* **Prop 2.5.** Condition 2 reads `o sigma_min(L) > r/2 + O(o n/K)` at `tau = K ln n`, and `sigma_min(L) >= LB(q)` because a collision-free `Lambda_L` contains no pair multiple. So:
  * `o > 1/(2c_0)`: `o >= 4` for `c_0 = 1/6`;
  * `o >= 3` for `c_0 >= 0.2349`, since `1/(2c_0) = 2.129`;
  * in both cases `K >= K_0(c_0, o)` is needed. The `tau_L` correction, of about `ln(2r+3)/c_0`, stands.

  Correct as stated.
* **§4 item 4.** `sigma* = r/2 + 1 + (positive) > r/2 > LB(q)`, so the 4-terms and dual pairs are light. Also `0.2349/2 = 0.11745`, which justifies `e^{-0.117 r tau}`. Correct.
* **Evidence vs proof.** `anneal.py` is labelled evidence. The exact minima are stated to be unused in any proof, which is true: Thm 2.8 uses only L1.1–L2.7. The only mislabel is M3.

## 6. Corrected statement

**Theorem (H_light for `sdc_q`; proved).** Let `q = 3 mod 4` be prime, `q >= 7`, and `r = 2q`. Let `lambda ∈ Z^{F_q ∪ {inf}}`, with `sum lambda = 0` and
`lambda` not an integer multiple of a pair `e_x - e_y`. Then:
1. `sigma_lambda >= q/3 = r/6`. This proof is elementary: Lemmas 2.1–2.3 and Props 2.4–2.5, where at `q = 7` the case `F = 0`, `n >= 6` uses `sigma >= n - 1 >= 5`. So H_light(1/6) holds for every such `q`.
2. For `q >= 11`, `sigma_lambda >= LB(q) = min(max(q/3, q/2 - 11 - 1.5 sqrt q), (q-3)/2)`. This uses Weil's bound for the squarefree monic quartic `prod(X - t_i)` with the quadratic character ([L] Lidl–Niederreiter Thm 5.41, located via citing papers).
3. Hence `sigma_lambda >= 0.23515 r` (`>= 0.2349 r`) at every regime prime `q >= 3251`, and `c_0(q) >= 1/4 - 5.5/q - 0.75/sqrt q`.
4. Sharpness: the dual pairs give `sigma = (q-1)/2`, so `c_0 <= (q-1)/(4q)`.
5. Some base 4-terms `e_x + e_y - e_z - e_inf` have `sigma = (q-1)/2`; the others have `(q+3)/2`.

**Exact data, not used in the proof.** The global minima of `sigma` over non-pair integral `lambda` are `3, 5` at `q = 7, 11` (exhaustive). For
`19 <= q <= 1091` and `q = 3251, 3259, 3271, 5479`, the global minimum is the exact 4-term minimum, for example `1571 = 0.2416 r` at `q = 3251`. The HiGHS runs are solver-optimal and only corroborate these values.

**Consequences.**
* Barrier Lemma 3.4 holds unconditionally, with `tau_D` as displayed in hlight §4.1.
* Barrier Thm 3.5, as corrected by R16 (`m' < q` relaxation, AGL(1,q) invariance), holds without hypotheses. The cap test is: suppose some point of `C_1` has `kappa_H <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` at a regime prime `q`, at `W = 2q tau` with `tau = O(log q)` at least `tau_D` and **all** thresholds of barrier Thm 2.3 (`tau_2`, `4`, `sigma* <= r+1`, `||lambda||_1 Delta <= 1/2`). Then route 5 is capped.
* Prop 2.5 then needs `o >= 3` at regime `q`, and `o >= 4` for all `q`, each with `K >= K_0`.

**What stays open.** CIRCLE(q) itself is not decided, and there is no growth or uniformity improvement.
