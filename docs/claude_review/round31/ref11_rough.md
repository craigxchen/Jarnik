# Referee report, Round 31 Part XI, block "rough": XI.1, the bound (XI.1), the exact recomputation, the twin and function-field tests

Scripts (my own; her script was not imported):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref11_rough/dmin.py`,
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref11_rough/checks.py`.

## 0. Results

1. **(XI.3), the dilated collision bound: CORRECT (Proved).** I checked every step: (XI.4), the Bezout step, the class counts, the lift from a collision to `l | B_ij`, Cauchy-Schwarz and the Mertens step. No step needs `v` to be coprime to `q`.
2. **(XI.1) on both families: CORRECT, with the quantifier stated (Proved).** The statement is: for every `C > 0` there is `P_0(C)` such that `M_max(q_P, C) <= (s/4)(1 + O_C(loglog P/log P))` for all `P >= P_0(C)`, on both families. So for every fixed `c < 1/12`, including `c = 1/24`, `M_max(q_P, C) <= (1/3 - c)s` for `P >= P_0(C, c)`.
3. **Novelty: none for genuine circles. It is a one-line corollary of the record (Proved).** At `lambda = 1`, (XI.3) is literally the first bullet of `growth/auxiliary.md` l.74: "`W >= 4 M log M - O_C(M)` if no prime `<= M` divides `N`, so `M <= (1/2 + o(1)) log R/loglog R`". `round31/ref8_extend.md` Prop 3.1 (l.178-186) states the same: "If no split prime `<= M` divides `N'`, this gives `W' >= 4M log M - O_C(M)`." On `q_P`, Theorem 1.1 gives `M <= (1/3 + o(1))s < P`. So no prime `<= M` divides `N' | q_P`, the record's coefficient 4 applies, and `M <= (1/4 + o(1))s`. Her own words "already implicit in the record" understate this: the bound is explicit there, and specialising it needs only the observation `M < P`. The previous judge's phrase "a proved bound of order s/3" (README_X §0, §3(e)) was therefore already out of date: the record gave `s/4` on `q_P`. The one genuinely new ingredient of XI.1 is the rational dilation `lambda`. It matters only for the fakes of XI.2, which are outside this block.
4. **Honest size.** This is a constant-factor improvement (`1/3` to `1/4`) on a family where the computed truth is `M = 1`. It is not a rate, and it does not use (T1) or the cocycle beyond the record's collision count. The (A1) question as posed (a fixed `c > 0`) is answered, but only because the record already answered it.
5. **Exact rows: REPRODUCED (Computed).** At `P = 200`: `s = 16`, `65,536` unit orbits, `d_min^2 = 15219819238844501369583650`, `floor(d_min/q^{1/4}) = 468`. At `P = 250`: `s = 20`, `1,048,576` orbits, `d_min^2 = 675416293412257838261903282453746`, `floor = 4556`. Both certify `M_max(q_P, 8) = 1`.
6. **Twin test: CORRECT (Proved arithmetic).** Function-field test: **CORRECT** for the undilated argument (Proved). One sentence about the record's binary and odd-characteristic families was not re-checked.

## 1. Definitions

- `q` is an odd squarefree product of rational primes `p ≡ 1 (mod 4)`, and `W = log q`.
- For each `p | q`, fix a Gaussian prime `pi_p` with `pi_p conj(pi_p) = p`.
- A mark is `A_i = eta_i prod_{p|q} f_ip`, with `f_ip ∈ {pi_p, conj pi_p}` and `eta_i ∈ {±1, ±i}`. Marks have norm exactly `q`.
- `H_ij = prod_{p: f_ip = f_jp} f_ip` and `Q_ij = prod_{p: f_ip ≠ f_jp} p`. Then `|H_ij|^2 = q/Q_ij`.
- `lambda = u/v ∈ Q^×`, written in lowest terms with `v > 0`.
- `D_ij = (A_j - A_i)/H_ij` and `B_ij = lambda D_ij`.
- Hypothesis (XI.2): for all `i < j`, `B_ij ∈ Z[i]` and `|lambda (A_j - A_i)| <= C q^{1/4}`.
- `q_P` is the product of all primes `p ≡ 1 (mod 4)` (family 1), or `p ≡ 5 (mod 8)` (family 2), in `[P, 2P]`, and `s = omega(q_P)`.

## 2. Check of XI.1

**(XI.4) Proved.**
- Let `u_p` be the number of rows choosing `pi_p`. The prime `p` separates exactly `u_p(m - u_p) <= m^2/4` unordered pairs.
- Hence `sum_{i<j} log Q_ij <= (m^2/4)W`, and the mean over `binom(m,2)` pairs is at most `mW/(2(m-1))`.
- `log|B_ij| = log|lambda(A_j - A_i)| - log|H_ij| <= log C + W/4 - W/2 + (1/2) log Q_ij`.
- Taking means gives `mean log|B_ij| <= log C - W/4 + mW/(4(m-1)) = log C + W/(4(m-1))`.
- Constant columns (`u_p ∈ {0, m}`) contribute nothing, so the normalisation by the global gcd is not needed.

**(XI.5) Proved.** The steps are as follows.
- *`D_ij` is integral.* Write `A_i = eta_i H_ij X`, `A_j = eta_j H_ij X'`. Then `D_ij = eta_j X' - eta_i X ∈ Z[i]`.
- *Bezout.* Pick rational integers `alpha, beta` with `alpha u + beta v = 1`. Then `D_ij/v = alpha B_ij + beta D_ij ∈ Z[i]`. So `v | D_ij`, and hence `v | H_ij D_ij = A_j - A_i`. This holds for the whole of `v` and uses no coprimality of `v` with `q`. Only the `l`-parts for odd `l <= m` are used later.
- *Primes of `v` dividing `q`.* If `p | v` and `p | q`, then all rows agree at `p`. Indeed, if `f_ip = pi_p ≠ f_jp`, then `pi_p | A_i` and `pi_p ∤ A_j`, so `p ∤ A_j - A_i`. Such primes are constant columns and cause no trouble.
- *Classes, case `a = 0`.* `#{(x, y) mod l : x^2 + y^2 ≡ q}` equals `l - (-1/l)`. This is `l - 1` for `l ≡ 1 (mod 4)` and `l + 1` for `l ≡ 3 (mod 4)`. I checked it by brute force for `l = 3, 5, 7, 11`.
- *Classes, case `a >= 1`.* All marks agree mod `l^a`, because `l^a | v`. The lifts of norm `q` mod `l^{a+1}` satisfy one affine condition `2xh + 2yk ≡ (q - x^2 - y^2)/l^a (mod l)`. This uses `2a >= a + 1`, and `l^a | q - x^2 - y^2` because the marks have norm exactly `q`. The pair `(x, y)` is not `≡ (0, 0)` because `l ∤ q`. So there are exactly `l` lifts. I checked this by brute force for `l ∈ {3, 5, 7, 11}`, `a ∈ {1, 2}` and three values of `q`: every class mod `l^a` of norm `q` has exactly `l` norm-`q` lifts.
- *Collision gives `l | B_ij`.* `H_ij` is invertible mod `l` because its norm divides `q`. So `l^{a+1} | D_ij`. Together with `v | D_ij` and `ord_l v = a`, this gives `lv | D_ij`. Hence `B_ij/l = u D_ij/(lv) ∈ Z[i]`.
- *Products.* Distinct collision primes are coprime rational integers, so their product divides `B_ij`. `B_ij ≠ 0` because the marks are distinct and `lambda ≠ 0`. Hence `log|B_ij| >= sum log l` over the collision primes of the pair.
- *Cauchy-Schwarz.* The number of classes is at most `l + 1`, so `sum n_h^2 >= m^2/(l+1)`. The colliding fraction of unordered pairs is `(sum n_h^2 - m)/(m(m-1)) >= (m/(m-1))/(l+1) - 1/(m-1)`. A negative lower bound for large `l` is harmless.
- *Mertens.* `sum_{l<=m} log l/(l+1) = log m + O(1)` and `theta(m) = O(m)`. Omitting `l = 2` (and the case `l | 2q` in the lift step) costs `O(1)`.
- *Combination.* `log m - O(1) <= log C + W/(4(m-1))`, which gives (XI.3).

The phrase "simultaneously for all four units" is harmless. In XI.1, `eta_i` varies per mark and `H_ij` carries no unit, so there is no unit factor in the count. Ordered and unordered pairs are used consistently.

**Why the constant is 4 and not 3 (Proved, record).** In Theorem 1.1, split primes `l <= M` that divide `N'` give no collisions: the marks are spread over levels there. In `ref8_extend` Prop 3.1, each split `p <= M` with `p ∤ N'` adds `2M log p/(p-1)`. The total `2M sum_{p≡1(4), p<=M} log p/(p-1) ~ M log M` lifts 3 to 4. On `q_P` every prime `<= M` is coprime to `q_P`, since `M < P`.

## 3. Check of the application (XI.1)

The following steps are Proved.
- *Size of `s`.* `s ~ P/(h log P)` with `h = 2` (family 1) or `h = 4` (family 2), by the prime number theorem in progressions.
- *Size of `W`.* `s log P <= W <= s log 2P`.
- *No large clusters.* Suppose a cluster has `m = floor(P/2)` marks. No prime `<= m` divides `q_P`, so (XI.3) gives `W >~ 2P log P`. This contradicts `W = O(P)` once `P >= P_0(C)`. The step is valid but unnecessary: Theorem 1.1 already gives `M <= (1/3 + o(1))s < P`.
- *The bound.* Suppose `M >= s/4`. Then `log M >= log P - loglog P - log(4h) + o(1)`. So `M - 1 <= s(log P + log 2)/(4(log P - loglog P - O_C(1)))`, which is `(s/4)(1 + O_C(loglog P/log P))`. The `+1` is absorbed because `s -> infinity`.
- *Comparison.* Theorem 1.1 gives `3M log M <= W + O_C(M)`, so `M <= (1/3 + o(1))s`. (XI.3) gives `(1/4 + o(1))s`. Theorem 1.2 gives only `M <= 1 + s + K_C`.
- *Consistency.* In terms of `R`: `log R/loglog R ~ s/2` on `q_P`, so `s/4` is exactly the record's `(1/2 + o(1)) log R/loglog R` for rough `N`.
- *Quantifier.* The statement holds for `P >= P_0(C)`. For small `P`, `s` can be `0` while `M` can be as large as `4`, so the threshold `P_0(C)` is needed.

## 4. Exact recomputation (Computed)

- *Method.* I enumerated all `2^s` conjugate choices and rotated each into `x > 0, y >= 0`. I sorted by the exact integer key `floor(y 2^K/x)` with `K = 4·bitlen(q) + 20`; distinct slopes differ by `> 1/q`, so the order is exact. I then verified strict angular order by cross products. Finally I took adjacent gaps plus the gap from the last point to `i` times the first.
- *Validation.* Brute force over all lattice points matches the code for `P = 10, 20, 30` (`d_min^2 = 2, 116, 2`).
- *Results.*

| family | P | s | orbits | d_min^2 | floor(d_min/q^{1/4}) | d_min^4 > 8^4 q |
|---|---:|---:|---:|---:|---:|---|
| 1 mod 4 | 200 | 16 | 65,536 | 15219819238844501369583650 | 468 | yes |
| 1 mod 4 | 250 | 20 | 1,048,576 | 675416293412257838261903282453746 | 4556 | yes |
| 5 mod 8 | 200 | 9 | 512 | 117088209047764 | 25 | yes |
| 5 mod 8 | 250 | 10 | 1024 | 2140461115510426 | 20 | yes |
| 5 mod 8 | 300 | 10 | 1024 | 421293725319704400 | 174 | yes |

Both of her rows agree exactly. The floor values were computed as exact integer fourth roots. The 5 mod 8 rows are mine, added as a check on the second family.

- *Certificate.* A closest pair is adjacent in angle, since a point strictly inside the shorter arc between them is closer to both. A chord is no longer than its arc. So `d_min > 8 q^{1/4}` gives `M_max(q_P, C) = 1` for all `C <= 8`. Proved.

## 5. Twin and function-field tests

**Twin test (Proved, arithmetic).**
- With common divisor of size `q/Q_ij` and interval length `C sqrt q`: `log|b_ij| <= log C + W/2 - W + log Q_ij`. The mean is at most `log C + W/(2(m-1))`.
- With the collision bound `log m - O(1)` (integers have `l` classes mod `l`), this gives `W >= 2(m-1) log m - O_C(m)`, i.e. `m <= (1/2 + o(1)) log q/loglog q`.
- This is exactly the size `m ~ (1/2) log q/loglog q` at which (T5) is built (reply10, (T5)). So the coefficient 2 is consistent with (T5) being sharp for the twin.
- The angle-fake remark is also consistent: VIII.3 has `W <= (23/3)M log M`, so `23/3 > 4`; its exponents are capped at two and its primes exceed `M`.

**Function-field test (Proved, undilated).**
- At a base prime of degree `d` in a constant quadratic extension, the nonzero norm fibre has `b^d - 1` elements (split, `d` even) or `b^d + 1` elements (inert, `d` odd).
- `|d I_d - b^d| <= sum_{e<=d/2} b^e = O_b(b^{d/2})`. So `sum_{d<=k} d I_d/(b^d + 1) = k + O_b(1)` and `sum_{d<=k} d I_d = O(b^k) = O_b(m)`.
- The degree version of (XI.4) gives `mean deg B <= c + deg q/(4(m-1))`.
- So `deg q >= 4(m-1) log_b m - O_b(m)`, provided no parent prime has degree `<= k = floor(log_b m)`.
- Not re-checked: the claim that the record's Parts V-VI binary families have degree of order `m^2 log m`. If true, it is consistent with the bound.

## 6. What this block does not cover

I did not referee XI.2 (the cocycle forcing `z_i = w + lambda A_i`, and the pinning `lambda = ±1`). The `s/4` bound for fakes depends on XI.2. Given XI.2, it follows from (XI.3) with `lambda`, which I checked. Minor: `final11.txt` links to `/scratchpad/astra/free/answer.md` rather than to `round31/astra_free_XI/answer.md`.
