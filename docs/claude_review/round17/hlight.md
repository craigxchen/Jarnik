# Round 17: H_light for sdc_q is a theorem (sharp up to O(q^{-1/2}))

Status: author's write-up, unrefereed. Scripts: `round17/hlight_checks/` (each `x.py` / `x.c` writes `x_out.txt`).
Labels as in the record: *proved* (all quantifiers), *exact* (integer / rational computation), *evidence*,
*[L]* literature located by search this round, *[M]* quoted from memory.

## 0. Summary

**Main result (proved, elementary).** For every prime `q = 3 mod 4`, `q >= 7`, and every integral
`lambda in Z^{F_q ∪ {inf}}` with `sum lambda = 0` that is not an integer multiple of a pair `e_x - e_y`, the
`sdc_q` slack satisfies

```text
sigma_lambda  >=  q/3  =  r/6            (r = 2q).
```

So **H_light(c_0) holds with c_0 = 1/6 for every q**, for all non-pair integral characters, light or not, at
every norm. The proof uses no character-sum estimate: a pattern-sum formula (Lemma 1.1), the second-moment
("2-design") identity, and a polynomial-degree bound for the zeros of the Paley convolution `f * chi` (Lemma 2.3,
`#zeros <= (q-3)/2`, sharp at pairs).

**Sharp form (proved; uses Weil's bound, [L]).** For every such `q` and `lambda`,

```text
sigma_lambda  >=  min( max(q/3, q/2 - 11 - 1.5 sqrt q),  (q-3)/2 ),
```

so at every regime prime (`q >= 3251`) `sigma_lambda >= 0.2349 r`, and `c_0(q) = 1/4 - O(q^{-1/2})`.
This is **sharp up to O(q^{-1/2})**: the explicit "dual pairs" (Prop 2.6) have `sigma = (q-1)/2 = (1/4 - 1/(4q)) r`
for every `q`, and the base 4-terms `e_x + e_y - e_z - e_inf` have `sigma = (q-1)/2` exactly (exact, `q <= 83`).

**Exact minima (exact).** With the reduction theorem the global minimum of `sigma` over **all** non-pair integral
`lambda` is a finite computation over 4-term characters once that minimum is below `(q-3)/2`:
* `q <= 43`: HiGHS MILP over all ternary `F = 0` characters (complete, no coefficient bound) reproduces the
  round-16 referee's values `sigma_min/r = 0.2143, 0.2273, 0.1842, 0.1957, 0.2097, 0.1977`; these are now exact
  global minima (the referee's MILP had `|lambda_x| <= 3`).
* all 93 primes `7 <= q <= 1091`, `q = 3 mod 4`, and the regime primes `q = 3251, 3259, 3271, 5479`, by exhaustive
  enumeration of 4-term characters (C, exact integers) plus the theorem: at `q = 3251` the exact global minimum is
  `sigma_min = 1571 = 0.2416 r` (attained by a 4-term character); see §3.

**Consequence (§4).** Barrier Lemma 3.4 and Theorem 3.5 (as corrected by the round-16 referee) hold
**unconditionally**: H_light is removed from the hypothesis of the circle test. With `c_0 = 1/6` the Diophantine
threshold is `tau_D = max(12(2 ln(8q+5) + ln(1000/C)), 4 ln(8q+5) + 8, 2 ln(1000 q (q+1)^2/C))`. The cap test is now:

> If `inf_theta kappa_H(hat u_theta) <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` for one regime `q` at `W = 2q tau`,
> `tau >= max(tau_D, thresholds of barrier Thm 2.3)`, then route 5 is capped.

**What is not done.** CIRCLE(q) itself (the near-zero question for `kappa_H` on the Fourier circle) is **not**
decided here; this round's task was its hypothesis H_light. Growth and uniformity are unchanged.


## 1. Setting and the pattern-sum formula

**Data** (barrier.md §1; referee_barrier §6, `sdc_sym.py`). `q = 3 mod 4` prime, rows `F_q ∪ {inf}` (`n = q + 1`;
`inf` is the base row, contained in no relative pattern), `chi` the quadratic character (`chi(0) = 0`,
`chi(-1) = -1`), `D = {0} ∪ QR`. The `r = 2q` relative patterns, all of weight 1, are

```text
P_a = {a} ∪ {x : chi(a - x) = 1} = a - D,      P'_a = {a} ∪ {x : chi(a - x) = -1} = a + D      (a in F_q).
```

For `lambda in Q^{n}`, `sum lambda = 0`: `sigma_lambda = sum_S (|lambda(S)| - 1/2) = N(lambda) - q`, where
`N(lambda) := sum_S |lambda(S)|` over the `2q` patterns.

**H_light(c_0)** (barrier Lemma 3.4, referee wording): every relevant integral `lambda` that is not a multiple of a
pair `e_x - e_y` and is light has `sigma_lambda >= c_0 r`. Here `k = e_1`, so every integral `lambda` is relevant,
and the only non-integral elements of `L_T` are irrelevant (fake.md §1; referee_barrier §3). We prove the
statement for **all** non-pair integral `lambda`, light or heavy.

**Notation.** Write `lambda = (f, lambda_inf)` with `f : F_q -> Z`, `F := sum_{x in F_q} f(x) = -lambda_inf`,
`h := F + f` (so `h(a) = lambda_a - lambda_inf`), `c := f * chi`, `c(a) = sum_t f(t) chi(a - t)`,
`T = supp f`, `n_lambda := ||lambda||_1`.

**Lemma 1.1 (pattern sums; proved).** For every `a in F_q`:
`lambda(P_a) = (h(a) + c(a))/2`, `lambda(P'_a) = (h(a) - c(a))/2`. Hence

```text
N(lambda) = sum_{a in F_q} max(|h(a)|, |c(a)|),       sigma_lambda = sum_a max(|F + f(a)|, |(f*chi)(a)|) - q .
```

*Proof.* With `1_QR = (1 - delta_0 + chi)/2` and `1_NR = (1 - delta_0 - chi)/2` (functions of the difference),
`lambda(P_a) = f(a) + sum_{s in QR} f(a - s) = f(a) + (F - f(a) + c(a))/2 = (h(a) + c(a))/2`, and the same with
`-chi` for `P'_a`. Finally `|(h+c)/2| + |(h-c)/2| = max(|h|, |c|)`. QED
(Exact check against the pattern definition: `check_formula.py`, 3200 random integral `lambda`, `q = 7..59`, 0
mismatches.)

**Lemma 1.2 (identities; proved).** For `q = 3 mod 4`:
1. `chi * chi = 1 - q delta_0`; hence `c * chi = F - q f`.
2. `sum_a c(a)^2 = q ||f||_2^2 - F^2`, `sum_a h(a)^2 = ||f||_2^2 + (q+2) F^2`, and
   `sum_S lambda(S)^2 = (q+1) ||lambda||_2^2 / 2`.
3. `c(a) = h(a) (mod 2)` for every `a`; `sum_a h(a) = (q+1) F`.
4. `|lambda(S)| <= n_lambda/2` for every pattern `S`.

*Proof.* (1) `(chi*chi)(a) = sum_t chi(t) chi(a - t)`; for `a = 0` it is `chi(-1)(q-1) = -(q-1)`; for `a != 0`,
`sum_t chi(t(a-t)) = sum_{t != 0} chi(a/t - 1) = sum_{u != -1} chi(u) = -chi(-1) = 1`. Then
`c * chi = f * (chi * chi) = F - q f`. (2) `sum_a c(a)^2 = sum_{t,t'} f(t) f(t') sum_a chi(a-t) chi(a-t')`, and the
inner sum is (substituting `b = a - t`) `chi(-1) (chi*chi)(t'-t) = -(1 - q delta_{tt'}) = q delta_{tt'} - 1`; the second identity is direct; the third
is `sum_a (h^2 + c^2)/2` with `||lambda||_2^2 = ||f||_2^2 + F^2`. (3) `c(a) = sum_{t != a} f(t) chi(a-t) = sum_{t != a} f(t)
= F - f(a) = F + f(a) (mod 2)`. (4) `sum lambda = 0`, so `lambda(S) <= sum_x max(lambda_x, 0) = n_lambda/2`, and likewise
for `-lambda(S)`. QED

**Remark (conference-matrix form).** With `S(a,t) = chi(a-t)`, `S(a,inf) = 1`, `S(inf,t) = -1`, `S(inf,inf) = 0`,
`S` is a skew conference matrix (`S S^T = qI`) and `c(a) = mu_a - mu_inf` for `mu = S lambda`. So
`N(lambda) = sum_a max(|lambda_a - lambda_inf|, |mu_a - mu_inf|)`: an uncertainty functional for the pair
`(lambda, S lambda)`. Pairs (sparse `lambda`) and "dual pairs" (sparse `S lambda`, Prop 2.6) are its two extremal
families. This is the structure behind the round-5 margin `sigma >= ||lambda||_1/2` (exact.md Thm 4.1(iii)).

## 2. The theorem

Throughout, `lambda` is integral, `sum lambda = 0`, `lambda != 0`, and `lambda` is not an integer multiple of a pair
`e_x - e_y` (`x, y in F_q ∪ {inf}`). Notation of §1.

### 2.1 Reduction to ternary characters

**Lemma 2.1 (proved).**
1. If `|F| >= 2`, then `sigma_lambda >= q + 2`.
2. If `|F| <= 1` and `||f||_inf >= 2`, then `sigma_lambda >= q - 1`.

*Proof.* (1) `N >= sum_a |h(a)| >= |sum_a h(a)| = (q+1)|F|` (Lemmas 1.1, 1.2(3)). (2) Take `b` with `|f(b)| >= 2`.
By Lemma 1.2(1), `q f(b) = F - sum_{a != b} c(a) chi(b - a)`, so `q |f(b)| <= |F| + ||c||_1 - |c(b)|` and
`N >= ||c||_1 >= 2q - 1`. QED

So it remains to treat **ternary** `lambda` (`lambda in {-1,0,1}^n`; then `F = -lambda_inf in {-1,0,1}`),
with `n_lambda = ||lambda||_1 = ||lambda||_2^2` even and `n_lambda >= 4` (`n_lambda = 2` is a pair).

### 2.2 Two tools

**Lemma 2.2 (defect identity; proved).** For ternary `lambda` with `n = n_lambda`, writing `t_S = lambda(S)`,

```text
sigma_lambda = 1 + (2/n) sum_S |t_S| (n/2 - |t_S|)  >=  1 + (1 - 2/n) nu(lambda),
nu(lambda) := #{ S : t_S notin {0, +n/2, -n/2} } .
```

*Proof.* `|t| = (2/n) t^2 + (2/n)|t|(n/2 - |t|)`; sum over `S` and use `sum_S t_S^2 = (q+1) n/2` (Lemma 1.2(2),
`||lambda||_2^2 = n`): `N = q + 1 + (2/n) sum_S |t_S|(n/2 - |t_S|)`. Each `t_S` is an integer with `|t_S| <= n/2`
(Lemma 1.2(4)); if `1 <= |t_S| <= n/2 - 1` then `|t_S|(n/2 - |t_S|) >= n/2 - 1`. QED

(This is the round-5 "defect identity" of exact.md Thm 4.1(iii), in pattern coordinates.)

**Lemma 2.3 (zeros of Paley convolutions; proved).** Let `f : F_q -> Z` with support `T`, `1 <= |T| <= (q+1)/2`,
and `q ∤ f(t)` for some `t in T`. Let `j_0 = min{ j >= 0 : sum_{t in T} f(t) t^j != 0 (mod q) }` (`0^0 = 1`). Then
`j_0 <= |T| - 1` and

```text
#{ a in F_q : (f*chi)(a) = 0 }  <=  (q-1)/2 - j_0 .
```

In particular, if `q | F` (e.g. `F = 0`), the Paley convolution `f * chi` has at most `(q-3)/2` zeros. Equality
holds for every pair `f = delta_x - delta_y` (exact, `check_lemmas.py`; it is Jacobsthal's count).

*Proof.* Let `m = (q-1)/2`. Euler's criterion gives `chi(u) = u^m` in `F_q` for every `u in F_q` (including
`u = 0`). Hence `(f*chi)(a) = P(a) (mod q)` with

```text
P(X) = sum_{t in T} f(t) (X - t)^m = sum_{j=0}^{m} binom(m,j) (-1)^j M_j X^{m-j}  in F_q[X],   M_j = sum_t f(t) t^j .
```

If `M_0 = ... = M_{|T|-1} = 0` in `F_q`, the Vandermonde matrix `(t^j)_{t in T, j < |T|}` (distinct `t`) is invertible
over `F_q`, so `f = 0 (mod q)` on `T`, contrary to the hypothesis. So `j_0 <= |T| - 1 <= m`. Since `0 <= j <= m < q`,
`binom(m, j) != 0 (mod q)`, so the coefficient of `X^{m - j_0}` is nonzero and the coefficients of `X^{m-j}`,
`j < j_0`, vanish: `P != 0` and `deg P = m - j_0`. Every zero of `f*chi` (an integer equal to 0) is a root of `P`, so
there are at most `m - j_0` of them. `M_0 = F`. QED

*Remark.* The bound uses only the size of the support, never the values: it is an uncertainty principle
(in the spirit of Chebotarev / Tao's prime-cyclic uncertainty principle [M]) for the Paley circulant, specific to `chi` being a power map. It is the step
that the round-16 judge expected to need Weil-type estimates; none are needed.

### 2.3 The ternary cases

**Proposition 2.4 (`|F| = 1`; proved).** Ternary non-pair `lambda` with `|F| = 1`:
`sigma_lambda >= (q-1)/2` if `n_lambda = 4`, and `sigma_lambda >= (q+1)/2` if `n_lambda >= 6`.
The value `(q-1)/2` is attained, e.g. by `e_0 + e_1 - e_z - e_inf` for a suitable `z` (exact, every `q <= 83`,
`exhaust_F1.py`).

*Proof.* By negation assume `lambda_inf = -1`, `F = 1`. Let `m_+`, `m_-` count the `+1`, `-1` entries of `f`;
`m_+ - m_- = 1`, `n = m_+ + m_- + 1`, so `m_- = n/2 - 1` and `n <= q + 1`. Here `h(a) = f(a) + 1 in {0,1,2}`, and
`h(a) != 0` for exactly `q - m_- = q + 1 - n/2` points `a`. At such `a` the two entries `U = (h+c)/2`, `V = (h-c)/2`
satisfy `U + V = h(a)`; both lie in `{0, +-n/2}` only if `|h(a)| = n/2` (`U = V = +-n/2` would give `|h| = n > 2`, and
opposite signs give `h = 0`). Since `|h| <= 2`:
* `n >= 6`: every such `a` contributes at least one nontrivial entry, `nu >= q + 1 - n/2`, and Lemma 2.2 gives
  `sigma >= 1 + phi(n)` with `phi(n) = (1 - 2/n)(q + 1 - n/2) = q + 2 - n/2 - 2(q+1)/n`, concave on `[6, q+1]`; its
  minimum is `min(phi(6), phi(q+1)) = min(2(q-2)/3, (q-1)/2) = (q-1)/2` for `q >= 5`.
* `n = 4`: only the `m_+ = 2` points with `h = 2` can be fully trivial, so `nu >= q - 1 - 2`, and
  `sigma >= 1 + (q-3)/2 = (q-1)/2`. QED

**Proposition 2.5 (`F = 0`; proved).** Ternary non-pair `lambda` with `F = 0` (so `lambda_inf = 0`,
`n := n_lambda = |T| >= 4`):
1. `sigma_lambda >= n - 1`;
2. if `n <= (q+1)/2`, also `sigma_lambda >= g(n) := (q(n-2) + 4n - n^2)/(n+2)`; `g(4) = q/3`, `g(6) = (q-3)/2`;
3. hence `sigma_lambda >= q/3` if `n = 4`, and `sigma_lambda >= (q-3)/2` if `n >= 6`.

*Proof.* Here `h = f`. *On `T`:* `U + V = f(a) = +-1`, so not both entries are in `{0, +-n/2}` (that would need
`n/2 = 1`): `nu >= n` from `T`. *Off `T`:* `h = 0`, `U = c/2 = -V`, `c(a)` even (Lemma 1.2(3)); both entries are
nontrivial iff `c(a) notin {0, +-n}`. Let `N_0, M, E` count the `a notin T` with `c(a) = 0`, `c(a) notin {0,+-n}`,
`|c(a)| = n`; `N_0 + M + E = q - n`.
(1) `nu >= n` and Lemma 2.2: `sigma >= 1 + (1 - 2/n) n = n - 1`.
(2) Lemma 2.3 (`F = 0`, `1 <= |T| <= (q+1)/2`, values `+-1`): `N_0 <= (q-3)/2`. `c` is odd on `T`, so
`sum_{a notin T} c(a)^2 <= q n - n` (Lemma 1.2(2)), hence `4M + n^2 E <= qn - n`. Eliminating `E`,
`M (1 - 4/n^2) >= q - n - (q-3)/2 - (q-1)/n`. Lemma 2.2 with `nu >= n + 2M`:

```text
sigma >= n - 1 + (2(n-2)/n) M >= n - 1 + (2n/(n+2)) ((q+3)/2 - n - (q-1)/n) = g(n).
```

(3) `n = 4`: `g(4) = q/3`. `n >= 6`: if `n >= (q-1)/2`, (1) gives `sigma >= (q-3)/2`. If `6 <= n <= (q-1)/2`,
`(n+2)(g(n) - (q-3)/2) = (1/2)[q(n-6) + 11n - 2n^2 + 6] >= (1/2)[(2n+1)(n-6) + 11n - 2n^2 + 6] = 0`
(using `q >= 2n + 1`). QED

**Proposition 2.6 (dual pairs: sharpness; proved).** For every `q = 3 mod 4` the ternary character
`f(x) = (chi(x) - chi(x-1))/2` (`x != 0, 1`), `f(0) = 0`, `f(1) = 1`, `lambda_inf = 0` has
`F = 0`, `||lambda||_1 = (q+1)/2` and `sigma_lambda = (q-1)/2` exactly.
*Proof.* With `g_b(x) = chi(x - b)`: `f = (g_0 - g_1)/2 + e`, `e = (delta_1 - delta_0)/2`. By Lemma 1.2(1),
`g_b * chi = 1 - q delta_b`, so `((g_0 - g_1)/2) * chi = (q/2)(delta_1 - delta_0)`, and `e * chi = -(g_0 - g_1)/2 =
-(f - e)`. So `c = (q/2)(delta_1 - delta_0) - f + e`: for `a != 0, 1`, `|c(a)| = |f(a)| = |h(a)|`; `c(0) = -(q+1)/2`,
`h(0) = 0`; `c(1) = (q-1)/2`, `h(1) = 1`. Lemma 1.1: `N = sum_{a != 0,1} |f(a)| + q = q + ||f||_1 - 1`.
`sum f = 0` and, by Jacobsthal (`sum_x chi(x)chi(x-1) = -1`), `#{x != 0,1 : chi(x) != chi(x-1)} = (q-1)/2`, so
`||f||_1 = (q+1)/2` and `sigma = (q-1)/2`. QED (Exact recomputation for 18 values `q = 7..3251`: `dualpair.py`.)

So the dual pairs are exact equality cases of Prop 2.5(1), and `c_0 <= (q-1)/(4q) < 1/4` for every `q`.

### 2.4 The 4-term characters: Weil

**Lemma 2.7 (proved, using Weil's bound [L]).** For ternary `lambda` with `F = 0`, `n_lambda = 4`:
`sigma_lambda >= q/2 - 11 - (3/2) sqrt q`.

*Proof.* `T = {t_1..t_4}`, `f(t_i) = f_i in {+-1}`, `sum f_i = 0`. For `a notin T`, `eps_i(a) := chi(a - t_i) in {+-1}`
and `c(a) = Phi(eps(a))`, `Phi(eps) := |sum f_i eps_i|`. Its Walsh expansion (exact, `weil4.py`) is
`Phi = 3/2 + (1/2) sum_{|I| = 2} s_I eps_I + (1/2) s_T eps_T` with signs `s_I in {+-1}` (depending on `f`). Since
`eps_I(a) = chi(prod_{i in I}(a - t_i))`,

```text
sum_{a notin T} |c(a)| = (3/2)(q - 4) + (1/2) sum_{|I|=2} s_I S_I + (1/2) s_T S_T,   S_I = sum_{a notin T} chi(prod_{i in I}(a - t_i)).
```

For `|I| = 2`, `sum_{a in F_q} chi((a - t)(a - t')) = -1` (`t != t'`; the substitution of Lemma 1.2(1)), and the two
points of `T \ I` change it by at most 2: `|S_I| <= 3`. For `I = T` the quartic `prod (X - t_i)` is squarefree, so
Weil's bound `|sum_a chi(g(a))| <= (deg g - 1) sqrt q` for `g` not a constant times a square gives `|S_T| <= 3 sqrt q`
([L]: the multiplicative-character form "`|sum chi(g)| <= (d-1) sqrt q` unless `g` is an `m`-th power" is quoted as
Theorem 6.2.2 of a textbook in arXiv:1910.09825; the classical source is Weil, cf. Lidl-Niederreiter *Finite Fields*,
theorem number unconfirmed). So `sum_{a notin T}|c(a)| >= 1.5q - 6 - 9 - 1.5 sqrt q`; on `T`, `|c| >= 1`. Hence
`N >= 1.5q - 11 - 1.5 sqrt q`. QED

### 2.5 Main theorem

**Theorem 2.8 (H_light for sdc_q; proved).** Let `q = 3 mod 4` be prime, `q >= 7`, `r = 2q`. For every integral
`lambda`, `sum lambda = 0`, not an integer multiple of a pair:

```text
sigma_lambda  >=  LB(q) := min( max(q/3, q/2 - 11 - 1.5 sqrt q),  (q-3)/2 ).
```

In particular **H_light(1/6) holds for every q** (elementary: Lemmas 2.1-2.3, Props 2.4-2.5), and
**H_light(c_0) holds with `c_0 = LB(q)/(2q) >= 0.2349` for every regime `q >= 3251`**; `c_0(q) -> 1/4`, and no
`c_0 > (q-1)/(4q)` is possible (Prop 2.6).

*Proof.* Lemma 2.1 (`>= q - 1`), Prop 2.4 (`>= (q-1)/2`), Prop 2.5 (`>= q/3` at `n = 4`, `>= (q-3)/2` at `n >= 6`),
Lemma 2.7 (`n = 4`). For `q >= 3251`, `q/2 - 11 - 1.5 sqrt q` is increasing and equals `1527.5...` at `q = 3251`
(rounded `sqrt q` upward in `c0table.py`, exact rationals), i.e. `>= 0.2349 r`; the minimum over the other classes is
`(q-3)/2 > 0.2349 r`. QED

**Corollary 2.9 (all characters, any norm).** Integer multiples `l(e_x - e_y)`, `|l| >= 2`, have
`sigma = |l|(q+1) - q >= q + 2`. Hence every nonzero integral `lambda` other than `+-(e_x - e_y)` has
`sigma_lambda >= LB(q)`, and the integral characters with `sigma < r/6` are exactly the pairs (`sigma = 1`).

## 3. Exact computations

All in `hlight_checks/`; `hlib.py` builds the patterns from the definition (independently of round-16 code).

| check | result | file |
| --- | --- | --- |
| Lemma 1.1 formula vs pattern sums | 0 mismatches, 3200 random integral `lambda`, `q = 7..59`; `chi*chi = 1 - q delta_0` | `check_formula_out.txt` |
| Lemma 1.2(2), Lemma 2.3 zero bound, Props 2.4/2.5 bounds | 0 violations on 19 200 random `f` (`q = 7..199`, supports `1..(q+1)/2`, values in `{+-1,+-2}`); zero bound attained by every pair | `check_lemmas_out.txt` |
| Walsh coefficients of `|sum f_i eps_i|`, `n = 4,6,8,10` | `n = 4`: `3/2`, `|coef| = 1/2` at `|I| = 2, 4` | `weil4_out.txt` |
| `|F| = 1`, `n = 4` (exhaustive up to `AGL(1,q)`) | min `sigma = (q-1)/2` for every `q <= 83` (Prop 2.4 sharp); `n = 6` (`q <= 47`): min `>= (q+1)/2` | `exhaust_F1_out.txt` |
| dual pairs | `sigma = (q-1)/2`, `||lambda||_1 = (q+1)/2`, 18 values of `q` up to 3251 | `dualpair_out.txt` |
| MILP (HiGHS) over **all** ternary `F = 0` `f`, `n >= 4`, `f(0) = 1`, `f(1) != 0` (WLOG by `AGL(1,q)` and negation) | optimal, `sigma_min = 3, 5, 7, 9, 13, 17` for `q = 7, 11, 19, 23, 31, 43` | `milp_tern_out.txt` |
| exhaustive 4-term `F = 0` (C, `min4.c`), all `7 <= q <= 1091` and `q = 3251, 3259, 3271, 5479` | 0 violations of `max(q/3, q/2 - 11 - 1.5 sqrt q)`; 4-term min `< (q-3)/2` for all `q >= 19` | `min4_small_out.txt`, `min4_regime_out.txt`, `summarize_min4_out.txt` |
| rigorous `LB(q)` table (exact rationals) | `c_0 = 1/6` (`q <= 199`), `0.2044` (499), `0.2211` (1019), `0.2349` (3251), `0.2387` (5479), `0.2476` (100003) | `c0table_out.txt` |

**Exact global minima.** By Theorem 2.8's case analysis every non-pair integral `lambda` outside the class
"ternary, `F = 0`, `n = 4`" has `sigma >= (q-3)/2`. So whenever the exact 4-term minimum is below `(q-3)/2` it is the
exact global minimum of `sigma` over all non-pair integral characters. This holds for every `19 <= q <= 1091` and for
the four regime primes; for `q = 7, 11` the MILP gives the global minimum `(q-1)/2`. Selected values:

| `q` | 7 | 11 | 19 | 23 | 31 | 43 | 103 | 199 | 499 | 1019 | **3251** | 3259 | 3271 | 5479 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `sigma_min` | 3 | 5 | 7 | 9 | 13 | 17 | 43 | 87 | 229 | 481 | **1571** | 1575 | 1581 | 2667 |
| `sigma_min / r` | .2143 | .2273 | .1842 | .1957 | .2097 | .1977 | .2087 | .2186 | .2295 | .2360 | **.2416** | .2416 | .2417 | .2434 |
| proved `LB/r` | .1667 | .1667 | .1667 | .1667 | .1667 | .1667 | .1667 | .1667 | .2044 | .2211 | .2349 | .2350 | .2350 | .2387 |

The minimum over `q >= 19` is `0.1842` at `q = 19`; the minima increase towards `1/4` at the `O(q^{-1/2})` rate of
Lemma 2.7 (the 4-term deficit `q/2 - sigma_min` is `54.5` at `q = 3251`, against the proved `98`).

## 4. Consequences for the barrier test (barrier.md §3.3, as corrected by referee_barrier)

1. **Lemma 3.4 is unconditional.** Its only hypothesis, H_light(`c_0`), holds with `c_0 = 1/6` for all `q`
   (Theorem 2.8). Its threshold becomes, with `r = 2q`, `n = q + 1`,

   ```text
   tau_D = max( 12 (2 ln(8q+5) + ln(1000/C)),  4 ln(8q+5) + 8,  2 ln(1000 q (q+1)^2 / C) )   (= O(log q + log 1/C)),
   ```

   and for regime `q` one may use `c_0 = 0.2349` (factor `2/c_0 = 8.52` instead of 12).
2. **Theorem 3.5 (with the referee's `m < N` relaxation and AGL(1,q)-invariance) holds without hypothesis.**
   Precisely: let `q >= 3251` be prime, `q = 3 mod 4`, `w = sqrt(8q ln 8q)`, `j = floor(w/2)`, `m' = ceil(4w)`
   (`m' < q`, R16-J2), and `tau >= max(tau_D, tau_2, 4)` with `tau = O(log q)` (barrier Thm 2.3 thresholds),
   `W = 2q tau`. If

   ```text
   inf_theta kappa_H(hat u_theta) <= eps_0 / (2 sqrt(m'+1) e^{pi^2/(4q)}),    u_theta(x) = cos(theta + 2 pi x/q), u_inf = 0,
   ```

   then every point of every Fourier circle `C_a` has `kappa_H <= eps_0/2`, half of `C_1` satisfies (2.2') for all
   light characters, and barrier Thm 2.3 produces gap-admissible `sdc_q` fakes at `W = 2q tau` with
   `kappa_H(hat u) <= eps_0`: **route 5 (MIN^pair / (R_H) form) is capped.** The test remains quantitative: an exact
   "circle not in `TC`" does not decide it (referee_barrier §4).
3. **Prop 2.5 (collision-free rational strata)** now has hypothesis 3 for free; with `sigma_min(L) >= LB(q)`,
   condition 2 needs `o > r/(2 LB(q))`: `o >= 4` for every `q` (`c_0 = 1/6`), `o >= 3` for regime `q`
   (`c_0 >= 0.2349 > 1/6`). The referee's corrections (`tau_L` enlarged by `~ ln(2r+3)/c_0`; `rho_1 >= 1/2`) stand.
4. **Light characters exist and are not exceptional.** Since `LB(q) < r/2 < sigma*(tau)`, the 4-terms
   (`sigma ~ 0.24 r`), the base 4-terms and the dual pairs (`sigma = (q-1)/2`) are light. Theorem 2.8 says they cost
   `e^{-tau sigma/2} <= e^{-0.117 r tau}` in (2.2'), which is what Lemma 3.4's measure count needs. The dual pairs
   (`||lambda||_1 = (q+1)/2`) show that the "intermediate norm range" of the round-16 plan is populated by light
   characters of the same slack as the 4-terms; the polynomial zero bound (Lemma 2.3) is what controls them.
5. **What CIRCLE(q) now is.** With H_light a theorem, the circle test is a single quantitative question about
   `kappa_H` on `C_1` of the regime `sdc_q` box, with no remaining hypothesis on the Diophantine side. This round
   does not decide it; route 5's status is unchanged otherwise (open, not capped, not advanced on growth).

## 5. Scope, `pal_q`, and what was not done

* **Growth / uniformity: no change.** Theorem 2.8 is a statement about the fake profile `sdc_q`; it removes a
  hypothesis from a conditional *cap*, and moves no bound for actual circles.
* **`pal_q` not treated.** `pal_q` is not translation-invariant (barrier §3.2), so H_light for it matters only for
  Prop 2.5. Its patterns are the Paley columns with multiplicity `b = 8` plus twins; the record already has
  `sigma >= ||lambda||_1` there (fake.md §1, `rho_1 = 1`). The method of §2 (pattern formula, defect identity, the
  zero bound of Lemma 2.3) should transfer, but I have not written or checked it.
* **Sylvester (`syl_k`)**: not a fake; additive quadruples have `sigma = 0` there (barrier §2.3), so no analogue of
  Theorem 2.8 holds.
* **Lemma 2.7 is the only use of Weil.** Without it the theorem still gives `c_0 = 1/6` for all `q`; with it
  `c_0 >= 0.2349` in the regime. The exact regime minima (`0.2416 r` at `q = 3251`) are computations, not used in
  any proof.

## 6. Open points

1. CIRCLE(q) itself (README-16 §7.1): upper bound (cap) or lower bound (non-zonotopal certificate) for
   `inf_{C_1} kappa_H` at the barrier scale. Untouched here.
2. The exact constant: is `min sigma / r = 1/4 - Theta(q^{-1/2})` attained only by 4-term characters for all
   `q >= 19`? (True for `19 <= q <= 1091` and the four regime primes computed.)
3. H_light for `pal_q` (for Prop 2.5 only).
4. Lemma 2.3 is a general uncertainty principle for `f -> f * chi`; its natural strengthening
   (`#zeros <= (q-1)/2 - j_0` together with control of the values `+-||f||_1`) might be useful for the
   lower-bound side of CIRCLE(q), where `<k', psi_1>` is a Fourier coefficient of a ternary `k`. Not explored.

## 7. Files (`round17/hlight_checks/`)

| file | content |
| --- | --- |
| `hlib.py` | patterns, `sigma` by definition and by Lemma 1.1, exact Paley convolution |
| `check_formula.py` | Lemma 1.1 and `chi*chi` check |
| `check_lemmas.py` | Lemma 1.2(2), Lemma 2.3 (zero bound, tightness), Props 2.4/2.5 bounds on random inputs |
| `weil4.py` | exact Walsh coefficients for Lemma 2.7 |
| `exhaust_F1.py` | exhaustive `|F| = 1`, `n = 4, 6` |
| `dualpair.py` | Prop 2.6 values |
| `milp_tern.py` | complete ternary MILP, `q <= 43` |
| `min4.c`, `run_min4.sh`, `summarize_min4.py` | exhaustive 4-term minima, regime primes, comparison with the bounds |
| `c0table.py` | rigorous `LB(q)` and `c_0(q)` |
| `regime_spot.py` | `q = 3251`: Lemma 1.1 vs direct pattern sums (12 random `lambda`, 0 mismatches); `sigma >= LB` on 60 random non-pair `lambda` |
| `anneal.py` | (evidence only) annealing search; it misses the dual pairs, illustrating that the intermediate range is structured, not random |
