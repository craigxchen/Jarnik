# Referee report on `round8/weakest.md` (rank reduction, simultaneous Lang–Waldschmidt, barriers)

Referee: adversarial, round 8. Scripts: `round8/referee_weakest_checks/`. Each runs from that directory with
`python3 <script> [args]`; its output is stored next to it as `*_out.txt`. The repository, `weakest.md` and
`weakest_checks/` were not modified.

## 0. Verdict

**The headline result does not survive. Every conditional theorem of §2 rests on a hypothesis that is false.**

1. **(fatal) `SLW(k, g, kappa)` is false for every `k`.** Proposition R1 (§1) refutes it for every `k`, every
   `kappa >= 0`, every constant `c` and every `g < 1`. At `k = 5` it already fails for every `g < 1/4` with
   `{0, +-1}` coefficients. The counterexamples are one-parameter polynomial templates
   `beta_i(t) = prod_l (t/b_l + i)^(sigma_il)` with `sum_l sigma_il b_l^(2j-1) = 0`. They have fixed coefficients,
   independent rows and union height `~ t^(2D)`, while `max |Lambda_i| ~ t^(-(2s+1))`. They are not circle
   configurations, and that is why the write-up's search for rank-5 *circle* families missed them.
   * Theorem 2.3 ("`SLW(5, g)` gives `M <= 458`") is vacuous.
   * Proposition 2.7's window `[1/4, 4/15)` contains no true instance (`SLW(4, g)` fails for `g < 3/10`).
   * The claims "five forms is the first `k` that can hold" and "the Dirichlet heuristic predicts `SLW(k, g)` for
     `g > 1/k`" are false. Exact computation gives finite violations of `SLW(5, 1/4)` with `c = 1/2`, `kappa = 1`
     (§1.3).
2. **(fatal) `SLW^log(theta, g, kappa)` is false.** Proposition R2 (§2) refutes it for every `theta > 0` and every
   `g < 1/2`, using templates with `k = s >= theta log(n+2)` forms. So Theorem 2.4, the claimed growth improvement
   `M = O((log R)^(1-eta))`, has a false hypothesis. Its support ("implied by every `SLW(k_0, g)`") is void.
3. **Against the fake circles** the logic is fine: the SLW inputs are violated by the round-5 fakes. The failure is
   different: the inputs are also violated by genuine algebraic numbers.
4. **What can be repaired (§3).** Weighting the heights by the coefficients
   (`prod_j Mah(alpha_j)^(-g B_j)`) keeps both proofs.
   * The weighted hypothesis is still false at `k <= 5` for `g < 1/4`.
   * At `k >= 6` it is open. It is new, and the Dirichlet heuristic, which is beaten at `k = 4, 5`, does not
     support it.
   * The circle-restricted version is a restatement of the conclusion.

   No conditional growth theorem from a standard, or even a supported, hypothesis remains.
5. **What survives, all re-derived (§4):**
   * the unconditional rank reduction (Theorem 1.3, checked exactly on 73,728 windows of actual circles);
   * the support lemma (Lemma 3.2), whose stated minimum `(q+7)/2` is wrong: the true minimum is `(q+3)/2`;
   * the LW blindness of the skew double core (Theorem 3.3, Lemma 3.4). It covers Lang's product form, an
     improvement on `lw.md`, within the scope noted in §5;
   * the abc, Subspace, Vojta and Bombieri–Lang propositions, with minor slips.

The verdict of §8.2.1 stands: no standard named conjecture weaker than Vojta's has been shown to suffice. Item 2 of
§8.2 ("the weakest hypothesis found that gives growth is `SLW^log`") must be withdrawn.

## 1. Fatal: the hypothesis `SLW(k, g, kappa)` is false for every `k`

### 1.1 Where the write-up goes wrong

`SLW(k, g, kappa)` (§2.1 of the write-up) is a statement about **arbitrary** systems of `k` forms in logarithms of
unit-circle Gaussian rationals. Its height factor `prod_(j in U) A(alpha_j)^(-g)` counts each number once,
whatever the size of its coefficients; the coefficients are charged only polynomially, by
`prod_j (1 + B_j)^(-kappa)`. The write-up supports the hypothesis by the Dirichlet heuristic (Remark 2.1.3) and
by the absence of rank-5 *circle* families (Remark 2.6, §7, §8.2.3: "`k = 5` is consistent with all known
data"). Circle families are a very special subclass of systems. One-parameter polynomial families with **fixed**
integer coefficients and growing heights are not covered by the heuristic, and they beat it for every `k`.

### 1.2 The counterexample

**Proposition R1.** Let `k >= 1`, `s >= 0`, let `b_1 < ... < b_D` be positive integers, and let
`sigma_1, ..., sigma_k in Z^D` be linearly independent with

```text
sum_(l=1..D) sigma_il b_l^(2j-1) = 0        for j = 1..s and i = 1..k.                          (R.1)
```

Then for every `c in (0, 1]`, every `kappa >= 0` and every `g < (2s+1)/(2D)`, the inequality (SLW) fails for some
system of `k` forms. In particular:

* (a) `b = (1, ..., 6)`, `s = 1` and the five `{0, +-1}` vectors
  `(1,1,-1,0,0,0), (1,0,1,-1,0,0), (1,0,0,1,-1,0), (1,0,0,0,1,-1), (0,1,1,0,-1,0)` show that
  **`SLW(5, g, kappa)` is false for every `g < 1/4`**: the hypothesis of Theorem 2.3 at `k = 5` is false;
* (b) `b = (1, ..., 5)`, `s = 1`, `(1,1,-1,0,0), (1,0,1,-1,0), (1,0,0,1,-1), (0,1,1,0,-1)`: **`SLW(4, g)` is false for
  every `g < 3/10`**, so the window `[1/4, 4/15)` of Proposition 2.7 contains no true instance;
* (c) for every `k` and `s`, `D = k + s` and `b = (1, ..., D)`: the `s x D` system (R.1) has a kernel of dimension
  `>= k`, so **`SLW(k, g, kappa)` is false for every `g < (2s+1)/(2(k+s))`, hence for every `g < 1`**.

So no instance of `SLW(k, g, kappa)` used anywhere in the write-up (all have `g < 1/4`) is true, for any `k`,
`kappa` or `c`.

*Proof.*
1. *Choice of `t`.* Put `B = 2 b_D` and `Q = lcm(b_1, ..., b_D) * prod_(p <= B) p`. Let `t` be a positive multiple
   of `Q`, `m_l = t/b_l` and `u_l = m_l + i`. Every prime `p <= B` divides `m_l`, so `N(u_l) = m_l^2 + 1 = 1 (mod p)`.
   Hence `N(u_l) > 1` and every prime factor of `N(u_l)` exceeds `B`; such a prime is `= 1 (mod 4)` (it divides
   `m_l^2 + 1`) and splits.
2. *Private primes.* Let `pi` be a Gaussian prime dividing `u_l`, `p = N(pi) > B`. Suppose `pi` also divides
   `u_l'` with `l' != l`, or `conj(u_l')` for some `l'`. Since `b_l u_l = t + i b_l`, `pi` divides
   `i(b_l - b_l')`, `i(b_l + b_l')` or `2 i b_l`, a nonzero element `i x` with `|x| <= B`. Then `p | x^2`, so
   `p <= B`, a contradiction. So the sets `P_l` of Gaussian primes of `u_l` are nonempty, pairwise disjoint, and
   `P_l` meets no `conj(P_l')` (including `l' = l`).
3. *The forms.* For a split `p` with a prime above it in `P_l`, write `u_l = eta prod pi_p^(e_p)` or with `conj(pi_p)`,
   and put `epsilon_p = +1` or `-1` accordingly. Then `u_l / conj(u_l) = (unit) prod_(p in P_l) alpha_p^(epsilon_p e_p)`,
   `alpha_p = pi_p/conj(pi_p)`. Put `b_ip = sigma_(i,l(p)) epsilon_p e_p`. The rows `(b_ip)_p` are linearly independent,
   because `sigma -> (sigma_(l(p)) epsilon_p e_p)_p` is injective (each `P_l` is nonempty). The `alpha_p` are distinct
   split-prime quotients, hence multiplicatively independent modulo roots of unity. Since
   `exp(sum_p b_ip Log alpha_p) = (unit) exp(2 i phi_i)` with
   `phi_i := sum_l sigma_il arctan(b_l/t)` (`arg u_l = arctan(1/m_l)`), there is `b_i0 in Z` with
   `Lambda_i = b_i0 Log i + sum_p b_ip Log alpha_p = 2 i phi_i`, and `|b_i0| <= 1 + 2 sum_p |b_ip|`.
4. *Smallness.* For `t >= 2 b_D`, expand `arctan x = sum_(j>=0) (-1)^j x^(2j+1)/(2j+1)` with `x = b_l/t <= 1/2`. By
   (R.1) the terms `j < s` cancel, so `|Lambda_i| = 2|phi_i| <= (8/3) ||sigma_i||_1 (b_D/t)^(2s+1)`.
5. *The right side of (SLW).* Let `n` be the number of `p`. Then `n <= sum_l omega(m_l^2 + 1) = o(log t)`.
   `sum_p log(1 + B_p) <= n log(1 + ||sigma||_inf) + sum_l log d(m_l^2 + 1) = o(log t)` (divisor bound), and
   `log(1 + B_0) = O(log log t)`. Since the primes are private, `sum_p log A(alpha_p) <= sum_p log(2p)
   <= sum_l log(m_l^2 + 1) + n log 2 <= 2D log t + o(log t)`. So the right side of (SLW) is
   `>= t^(-2Dg - o(1))`, while the left side is `<= K t^(-(2s+1))` with `K` independent of `t`. If `2Dg < 2s + 1`,
   (SLW) fails for all large `t` in `Q Z`.
6. *(a), (b), (c).* The listed vectors satisfy (R.1) with `s = 1` and are independent (checked exactly). In (c) the
   kernel of an `s x D` matrix has dimension `>= D - s = k`, and `(2s+1)/(2(k+s)) -> 1` as `s -> infinity`. QED

### 1.3 Exact computation (`slw_template_counterexample.py`)

Exact factorisation of the `N(u_l)`, exact ranks of the coefficient rows, angles at 300 digits,
`gamma := log(1/max|Lambda_i|) / sum_(p in U) log H(alpha_p)` with `H(alpha_p) = max(p, 2|Re pi_p^2|)` (the write-up's
own normalisation), and `gamma_eff`, the exponent that remains after charging every factor of (SLW) with
`c = 1/2`, `kappa = 1` (including `B_0`). (SLW) with constants `(c, kappa, g)` fails exactly when `gamma_eff > g`.

| system | `k` | coefficients | `log t` | `n` | rank | `gamma` | limit | `gamma_eff(1/2, 1)` |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| (a) `D = 6`, `s = 1` | 5 | `{0,+-1}` | 32.56 | 24 | 5 | 0.2435 | `1/4` | 0.146 |
| (b) `D = 5`, `s = 1` | 4 | `{0,+-1}` | 30.16 | 18 | 4 | 0.2955 | `3/10` | 0.195 |
| `D = 7`, `s = 2` | 5 | integers, entries `<= 3` | 37.07 | 21 | 5 | 0.3507 | `5/14` | 0.270 |
| `D = 9`, `s = 4` | 5 | integers, entries `<= 19` | 41.70 | 30 | 5 | 0.4898 | `1/2` | **0.351** |

The last two rows are *finite* counterexamples to `SLW(5, 1/4, 1)` with `c = 1/2`: at the stated `t` the five forms
have independent rows and violate (SLW) outright. Row (a) approaches `1/4` from below at the same rate as the
golden family that the write-up itself uses to refute `SLW(4, g < 1/4)` (Prop. 2.5(a), Table 7.1: `0.2440` at
`log N' = 135.6`).

### 1.4 Consequences for the write-up

* Theorems 2.3, 2.4 and Proposition 2.7 are correct implications with false hypotheses. They prove nothing.
  In particular "`SLW(5, g)` gives `M <= 458`" and "`SLW^log` gives `M = O((log R)^(1-eta))`" (§0, §2, §8) have no
  content.
* Proposition 2.5 is true but far from sharp: the falsity thresholds are `1` (not `1/4`, `2/7`, `1/2`) for every
  `k`. The claims "five forms is the first `k` that can hold" (§0, Remark 2.1.4, Theorem 2.3, §8.2.3), "the golden
  family attains, but does not beat, the Dirichlet exponent" (Remark 2.1.3) and "`SLW(5, g)` for `g` in `(1/5, 1/4)` is
  consistent with everything known" (Remark 2.6) are false.
* The Dirichlet heuristic of Remark 2.1.3 is not a valid heuristic for this hypothesis: one-parameter polynomial
  families have density zero among systems of bounded height but carry the systematic cancellation (R.1).
* Remarks 2.4b, 2.4c, 3.6.3 ("`SLW(2, g)` would see the profile for every `g < 2/3`"), the §8.1 ladder, §8.2.2–5 and
  the `k = 2, 3` item of §8.3 rest on `SLW` and must be withdrawn.
* Remark 2.4a explains the slip. It realises a system `Log(gamma_i / conj gamma_i)` as a circle configuration of
  norm `prod_i N(gamma_i)`, and it needs "squarefree normalised norm" to identify that norm with the union height.
  In the templates the five `gamma_i` share the primes of the `u_l` with mixed signs, so the circle norm is about
  `t^(2 sum_i ||sigma_i||_1) = t^30` while the union height is `t^12`. `SLW` charges only the union height.

## 2. Fatal: `SLW^log(theta, g, kappa)` is false

The write-up's support for `SLW^log` is that it is implied by every `SLW(k_0, g)` (Lemma 2.2). By §1 those are
all false, so this support is void. `SLW^log` asks for (SLW) only when `k >= max(k_0, theta log(n+2))`, so a
fixed-`k` template with `t -> infinity` (and `n -> infinity`) is not automatically covered. It is covered once `k`
grows slowly.

**Proposition R2.** For every `theta > 0`, `g < 1/2`, `kappa >= 0`, `k_0 >= 1` and `c in (0, 1]` there is a system with
`k >= max(k_0, theta log(n + 2))` forms (in `n` numbers, independent rows) that violates (SLW). So
`SLW^log(theta, g, kappa)` is false, in particular in the range `theta < 1/log 5.5`, `g < 1/4` of Theorem 2.4.

*Proof.* Use Proposition R1's construction with `D = 2s`, `b = (1, ..., 2s)`, `k = s` forms, and `s` independent
kernel vectors of (R.1) taken among the Cramer vectors (signed `s x s` minors). Their entries satisfy
`||sigma||_inf <= S_s := s! (2s)^(s(2s-1))`, so `log S_s = O(s^2 log s)`. Here `Q_s = lcm(1..2s) prod_(p <= 4s) p = e^(O(s))`.
Let `T` be large, and choose `s = s(T) := max(k_0, ceil((theta + 1) log log log T))`.

1. *Few primes for most `t`.* Consider `t in (T, 2T]` with `Q_s | t`; there are `T/Q_s + O(1)` of them, and
   `Q_s <= T^(1/4)`. Fix `l`. A prime `p <= 4s` never divides `m_l^2 + 1` (step 1 of R1). A prime
   `4s < p <= sqrt T` divides `m_l^2 + 1 = (t/l)^2 + 1` for at most `2(T/(Q_s p) + 1)` of these `t` (two classes of
   `t/Q_s` mod `p`). At most 4 prime factors of `m_l^2 + 1 <= 5 T^2` exceed `sqrt T`. Summing,
   `sum_t omega(m_l^2 + 1) <= 2 (T/Q_s)(log log T + O(1)) + 2 pi(sqrt T) + 4 #t`, so the mean of
   `sum_l omega(m_l^2 + 1)` over these `t` is `<= 2s(2 log log T + O(1))`. By Markov's inequality at least half of the
   `t` have `n <= n_T := 4s(2 log log T + O(1))`. Fix such a `t`.
2. *Coverage.* `theta log(n + 2) <= theta (log log log T + log(8s) + O(1)) <= s` for large `T`, by the choice of `s`.
   Also `s >= k_0`. So `SLW^log` applies to this system.
3. *Violation.* As in R1 (steps 4–5), with explicit dependence on `s`:
   * `log(1/max|Lambda_i|) >= (2s+1) log T - (2s+1) log(2s) - log(8 s S_s) = (2s+1) log T - O(s^2 log s)`;
   * `(n+1) log(1/c) + kappa sum_j log(1 + B_j) <= (n+1) log(1/c) + kappa (n log(1 + S_s) + n log(1 + log_5(5T^2)) + O(log(s S_s log T)))`,
     which is `O_(c,kappa)(s loglog T (s^2 log s + loglog T)) = o(log T)`;
   * `g sum_p log A(alpha_p) <= g (4s log(2T) + n log 2)`.

   (SLW) would force `(2s + 1 - 4sg) log T <= o(log T) + O(s^2 log s)`. Since `g < 1/2`, `2s + 1 - 4sg >= 1`, and
   `s^2 log s = o(log T)`, this fails for large `T`. QED

*Illustration* (`slw_log_illustration.py`): for `s = 3, 4, 5, 6` (`k = s` forms, `D = 2s`) at `log t` between 16 and 39,
`n <= 40`, the systems satisfy `k >= theta log(n+2)` at the largest admissible `theta = 1/log 5.5` (the threshold is
`<= 2.2`), and `gamma` is between `0.527` and `0.560` (limit `(2s+1)/(4s)`). Even after charging `c = 1/2`, `kappa = 1`,
`gamma_eff` lies between `0.297` and `0.376`. So every instance of `SLW^log(theta, g, kappa)` with
`theta <= 1/log 5.5`, `k_0 <= 6`, `c >= 1/2`, `kappa <= 1` and `g < 0.29` is already violated by one explicit system
(the `s = 6` row at `log t = 38.66`).

**Against the fakes.** The growth claim of Theorem 2.4 does use an input that the round-5 fake circles violate
(rank-`k` sub-configurations of a fake have `max |Lambda_i| <= C e^(-W/4)` with union height `<= e^(W + n log 2)`).
That check passes. The defect is the other one: the input is false for genuine algebraic numbers too.

## 3. What can be salvaged from Section 2

The proofs of Theorems 2.3 and 2.4 apply the hypothesis only to systems coming from one normalised circle:
rows `a_(x_i) - a_(x_0)` with `|b_ij| <= e_j`. Two repairs keep those proofs and drop the template systems.

### 3.1 Weighted heights

Let `Mah(alpha)` be the Mahler measure (`Mah(alpha_p) = p` for `alpha_p = pi_p/conj(pi_p)`). Define `SLW_w(k, g, kappa)` as
`SLW(k, g, kappa)` with `prod_(j in U) A(alpha_j)^(-g)` replaced by `prod_j Mah(alpha_j)^(-g B_j)`, and `SLW_w^log`
likewise.

* *Liouville baseline.* For the split-prime numbers `alpha_p`, the write-up's own Lemma 1.1 gives
  `|Lambda| >= sqrt2 prod_p p^(-|b_p|/2)` for every nonzero single form, so the weighted inequality holds at
  `g = 1/2` for every `k`. `SLW_w(k, g)` with `g < 1/4` is therefore a *simultaneous improvement of Liouville* from
  exponent `1/2` to below `1/4`.
* *The proofs go through.* In the circle application `B_j <= e_j`, so `sum_j B_j log Mah(alpha_j) <= W` (even without the
  `r log 2` term). Theorem 2.3 holds with `SLW_w(k, g, kappa)`, `g < 1/4`, in place of `SLW`, and Theorem 2.4 with
  `SLW_w^log` in place of `SLW^log`. I re-derived both proofs (the displayed inequality in the proof of 2.3 and the
  choice of `P` are correct).
* *Still false at `k <= 5`.* The `{0, +-1}` templates of R1(a), (b) have `B_p = e_p`, so their weighted height is
  `sum_l log N(u_l) <= 2D log t` and `gamma_w -> (2s+1)/(2D)`: `1/4` for `k = 5` and `3/10` for `k = 4`
  (`slw_template_counterexample.py`: `gamma_w = 0.2479` at `log t = 32.6`, and `0.2958`). So `SLW_w(5, g)` is false
  for every `g < 1/4`, and `SLW_w(4, g)` for every `g < 3/10`. A second rank-5 `{0, +-1}` family, with `s = 2`, is
  `b = (1,2,4,5,7,9,10,11,13,14)` (`weighted_pm1_search.py`), also with `gamma_w -> 1/4`. The theorem therefore needs
  `k >= 6`, and the bound becomes `M <= (5.5^5 + 1)/2 = 2516`.
* *Open at `k >= 6`.* Nothing I tried refutes `SLW_w(6, g < 1/4)`:
  * `s = 1` templates have `gamma_w <= 3/(2(k+1)) <= 3/14`;
  * integer-coefficient templates pay for their coefficients (`gamma_w = 0.176` and `0.043` for the `D = 7` and `D = 9`
    rows of §1.3);
  * `{0, +-1}` templates with `s = 2` have rank `<= 5` for every `b` in `binom({1..20}, 8)`, `binom({1..18}, 9)` and
    `binom({1..20}, 10)` (`weighted_pm1_search.py`; ranks 2, 3, 5 respectively), and a weighted exponent `>= 1/4`
    needs `D <= 10` at `s = 2`.

  A rank-6 circle family on arcs `C sqrt R` would refute it at the boundary; none is known.
* *Support.* None beyond the absence of counterexamples. The Dirichlet heuristic of Remark 2.1.3 predicts the
  threshold `1/k`, and it is beaten at `k = 4` (`3/10 > 1/4`) and `k = 5` (`1/4 > 1/5`) by one-parameter families. So it
  cannot be quoted for `k = 6`. The hypothesis is new and unnamed.
* *Comparison with Vojta.* `SLW_w(6)` concerns 7-point configurations, as `V(7, 1/4)` does. The "one point earlier"
  advantage claimed in Remark 2.4c, §0 and §8.2.4 disappears.

### 3.2 Circle-restricted form

`SLW_circ(k, g)`: there is `c > 0` such that for every normalised configuration on `x^2 + y^2 = N'` and every `k + 1` of
its points whose exponent vectors have affine rank `k`, `max_i |theta_(x_i) - theta_(x_0)| >= c N'^(-g)`. Theorem 2.3
holds verbatim with `SLW_circ(k, g)`, `g < 1/4`. By Remark 2.4a this says that rank-`k` sub-configurations never lie
on arcs of length `c R^(1-2g)`. It is the conclusion of the uniform theorem, restricted to rank `k` and taken with a
power margin; Theorem 1.3 is the only content added. It is consistent with the random model for `g > 1/5`
(`k = 5`). Vertical circle templates do not refute it: fibres of `(sum sigma b, sum sigma b^3)` over `{+-1}^D` have at
most 3, 3, 4 points for `D = 8, 9, 10` and `b` in `{1..18}` (`vertical_fibres.py`), which agrees with §7.3.

### 3.3 Net effect

No conditional growth theorem survives from a hypothesis that is not either false or a restatement of the
problem. The corrected conditional statements are:

* `SLW_w(k, g, kappa)` for some `k >= 6` and `g < 1/4` implies the uniform theorem with
  `M <= (5.5^(k-1) + 1)/2`;
* `SLW_w^log(theta, g, kappa)` with `k_0 >= 6`, `theta < 1/log 5.5` and `g < 1/4` implies
  `M = O((log R)^(1-eta))`.

Both hypotheses are new, open and unsupported. Both use an input that the fakes violate, so the barrier map does
not exclude them. Neither is "the weakest standard hypothesis".

## 4. The unconditional results, re-derived

| claim | what I did | verdict |
| --- | --- | --- |
| Lemma 1.1 (gap) | re-derived: `beta = prod pi^(v+) conj(pi)^(v-)` is primitive of odd norm `e^(w(v))`, at distance `>= 1/sqrt(2 N)` (relative) from the axes and diagonals | correct |
| Lemma 1.2 (difference forms; bound `5 + 2 sum e_j` on the `Log i` coefficient) | re-derived the `2 pi i m` bookkeeping | correct |
| Theorem 1.3 (rank reduction) | re-derived steps 1–5; exact check on actual circles | correct |
| Corollary 1.4 | both directions; `(5.5)^(theta log(2+W'/log 5))` arithmetic | correct (a reduction, not a growth result) |
| Lemma 2.2, Theorems 2.3–2.4 as implications | re-derived the log inequality and the choice of `P` | correct, but vacuous (§1–2) |
| Prop. 2.5 | re-derived (b): Pell identity, gcd, Størmer, rank 3, radical `7/8` | correct, not sharp (§1.4) |
| Lemma 3.1 | skew-Hadamard block computation, cut sizes, pair separation, character lattice | correct |
| Lemma 3.2 (support `>= (q+1)/2`) | re-derived both cases | correct; the stated minimum is wrong (minor 5.1) |
| Theorem 3.3, Lemma 3.4, Cor. 3.5 | re-derived every step | correct within scope; minor 5.2, scope 5.3 |
| Props 4.1–4.2, Cor. 4.3 | re-derived the triangle computation and the radical count | correct as scoped |
| Lemma 5.1, Theorem 5.2, Prop. 5.3 | re-derived the heights, the hyperplane count and the pigeonhole and line counts | correct |
| Props 6.1–6.4 | re-derived the projection argument, the Cayley table and the denominators | correct; minor 5.4–5.5 |

### 4.1 Theorem 1.3

* *Step 1* needs `Delta < pi/2`. That is exactly why `W >= W_0(C)` is assumed: distinct exponent vectors fail at
  `Delta = pi/2`, for example on `N = 65` with `(4 - 7i)/(-7 - 4i) = i`, which my checker hit before I added the
  filter.
* *Step 2.* `|theta_x - theta_y - theta_x' + theta_y'|/2 <= Delta`. Lemma 1.1 then gives `w(v) >= rho`.
* *Step 3.* `w` restricted to `V = Lambda (x) R` is a norm. Every `u in D` has `w(u) <= W` because
  `|a_xj - a_yj| <= e_j`. The volume comparison gives `|D| <= (1 + 2W/rho)^k`.
* *Step 4.* `0` and `+-(a_x - a_1)` are `2M - 1` distinct elements of `D` (lexicographic order).
* *Step 5.* `W_0 = 18(2 log^+ C + log 2)`, so `rho >= 4W/9`.

*Exact check* (`rank_packing_actual.py 20000`). Every `N <= 20000` with all prime factors `= 1 (mod 4)` and
`r_2(N) >= 16`; every window of 2 to 8 consecutive points with `Delta < pi/2`; normalised; `C` the exact arc constant
of the window.
* 73,728 windows and 5,820,936 pairs `(u, u')`. The smallest slack is `w(u - u') - rho = 0.087`, and it is never
  negative.
* The packing bound was applicable (`rho > 0`) in 34,512 windows and was never violated.
* Ranks seen: at most 3.

The theorem is a clean and correct reduction: the uniform theorem is equivalent to bounded rank, and growth
follows from rank `<= theta log W`. Since the round-5 fakes have full rank `M - 1`, any rank bound must use
information the fakes violate. The write-up says so (§8.3); it is consistent with the barrier map.

### 4.2 Lemmas 3.1 and 3.2

Both cases of Lemma 3.2 are right.
* If `P != 0 (mod q)`, then `|Z| <= d`.
* If `P = 0`, then the moments `sum_x mu_x x^m` vanish for `m <= d` (using `sum_x x^m = 0` for `m <= q - 2`). A zero
  set `|Z| >= d + 1` leaves `|supp mu| <= d`, so the Vandermonde forces `mu = 0 (mod q)` and then `q | lambda`.

`skew_support_minimum.py` confirms:
* the lattice of rational characters `Z^M_0 + Z lambda*`: integrality of the cut values is equivalent to membership
  on 9,000 random vectors, `q = 7, 11, 19`;
* the bound: no set of `2q - (q+1)/2` cuts carries a nonzero zero-sum vector. All 1,001 sets (`q = 7`) and all
  74,613 sets (`q = 11`) have full rank modulo `2^61 - 1`, and full rank modulo a prime implies full rank over `Q`.

### 4.3 Theorem 3.3, Lemma 3.4, Corollary 3.5

* *Theorem 3.3.* The case split `B = |b_0|` versus `B = ||v||_inf` is right, and so is
  `(1+eps)(q+1)/2 - q(1+eps/2)/2 = (1 + eps + q eps/2)/2 >= eps q/4`.
* *Lemma 3.4.* Steps 1–4 are right:
  * step 1: the measure bound `(4l/pi + 2) rho`, and `||lambda||_inf >= 2/M` from the lattice;
  * step 2: both block maps are injective; I re-derived `lambda(S'_a) = (K lambda~)_a/2` and `(K lambda~)_inf = -sum lambda = 0`;
  * step 3: the second large coordinate, from `K^T K = M I`;
  * step 4: `||v||_inf^(-eps) |v_(T**)|^(-(1+eps)) <= (2 sqrt q)^(1+eps) |v_(T*)|^(-(1+eps))`.

  The union bound converges, because injectivity reduces the sum to one block, giving
  `sum <= 2 (2 sqrt q)^(1+eps) (1 + 2 zeta(1+eps))^q`.
* *What this improves.* This covers Lang's *product form*, with `B` in the numerator, every `b_0` (including the
  homogeneous case `b_0 = 0`) and every constant `C_*`. The round-8 LW referee found that `lw.md`'s union bound covers
  only the max form. On this point `weakest.md` improves on `lw.md`, within the skew-double-core class.
* *Corollary 3.5.* `P(success and LW) >= P(success) - P(LW fails)` holds whatever the correlation between success and
  `delta`, so the corollary is right as stated, given its inputs.

## 5. Minor issues

1. **Lemma 3.2: the minimum is `(q+3)/2`, not `(q+7)/2`.**
   * The write-up says the minimum found is `(q+7)/2`, "attained by `K(e_a - e_b)`, so the lemma is sharp up to `+3`".
     Its own output reports 6, not 7, at `q = 7`.
   * The single columns `lambda = K e_a` are zero-sum and have exactly `(q+3)/2` nonzero cuts for every `q` tested
     (`7 <= q <= 67`). The values are `lambda(S_b) = -(M/2) delta_ab` and `lambda(S'_b) = 1 + chi(a - b)` for `b != a`.
   * For `q = 7, 11` the exact minimum is `(q+3)/2` (`skew_support_minimum.py`, certified as in §4.2). This agrees with
     the certified `M/2 + 1` of the round-8 LW referee.
   * So the lemma is sharp up to `+1`. The limiting support fraction `1/4` and every use of the lemma are unchanged.
2. **Lemma 3.4, `tau_0`.** Step 5 needs `log(6M) + (1+eps) log(2 sqrt q) <= 2q`. This fails for large `eps` at
   small `q` (e.g. `eps = 10`, `q = 7`). Restrict to `eps <= 1`, or add `(8/(eps q)) (1+eps) log(2 sqrt q)` to `tau_0`.
   The conclusions are unaffected.
3. **Scope of the LW barrier (Corollary 3.5).** The statement "LW plus local data cannot prove growth" holds for this
   input class only, under four restrictions:
   * it is conditional on `exact.md` Theorem 5.1, an unrefereed draft (stated);
   * it covers cluster-determined character forms only (Remark 3.6.1, stated);
   * it covers forms written in the basis `alpha_T`. LW applied to other algebraic numbers, such as composite
     quotients `z_x/z_y` or Kummer roots, has a different height product. The valuation lemma that controls this over
     `Q(i)` is `lw.md` Lemma 2.2, which is neither cited nor proved here;
   * the angle-model class does not model additive cluster data. For example `Log((z_x + z_y)/(2 z_u)) = O(Delta)`
     is cluster-controlled but is not a character form (the point the LW referee raised for `lw.md`). The class also
     does not model the prime-level parity law of `round6_sharp/sharp.md`.
4. **Vojta threshold (§6.1, §8.1).**
   * The non-density threshold for `V(m, eps)` is `eps < (m-6)/2`. From `h_K >= ((m-6)/2) log R' - O(1)`, this gives
     `eps < 1/2` at `m = 7`, which is the referee-corrected `conditional.md`.
   * The write-up's `(m-6)/4` is the particular choice `eps_alpha` of `conditional.md` Theorem 5.2. It is sufficient
     but is not the threshold.
5. **Proposition 6.1(c)** ("there is no Vojta-type hypothesis that gives growth without giving the uniform bound") is a
   meta-claim and is not proved. What is proved:
   * every instance `V(m, eps)` with `eps < (m-6)/2` gives the uniform bound;
   * non-density for one `m` is equivalent to the uniform theorem.

   Statements with `m` growing with `R`, for example Vojta with explicit dependence of the exceptional set on `m`,
   are not excluded.
6. **§3.1, normalisations of LW.** "Every standard normalisation with `C(eps)^n` or `C(eps)` in front implies
   `LW(eps, C_*)` with `C_* = max(1, C(eps))`" is backwards when `C(eps) < 1`. It then implies `LW(eps, C(eps))`, not
   `LW(eps, 1)`. This is harmless, because Theorem 3.3 and Lemma 3.4 hold for every `C_* > 0`.
7. **Evidence misapplied.** §7.2 (rank-5 windows of actual circles, median `gamma = 0.105`), §7.3 and the "rank `<= 4`
   for every known family" observation are evidence about *circle* configurations. They are presented as support for
   `SLW(5)`, which is a statement about *arbitrary* systems and is false (§1). They remain valid as evidence for
   `SLW_circ`.
8. **"The Dirichlet exponent for that many forms is `1/(theta log n) -> 0`, far below the `1/4` used"** (after
   Theorem 2.4). The template systems of R2 have `k = s >= theta log(n+2)` forms and exponent `> 1/2`. The heuristic
   does not describe the truth.
9. Remark 2.1.1 and §5.2 quote the literature from memory (flagged in the write-up). In particular the attribution
   of simultaneous lower bounds to Philippon–Waldschmidt, and Evertse's `(2^(60 m^2) delta^(-7m))^s`, could not be
   checked here either. Neither affects a proof.

## 6. Corrected statement

> **Status: no growth result.** The round's question remains unanswered: no standard named conjecture weaker than
> Vojta's is known to give `M = o(log R/loglog R)`.
>
> **Conditional (corrected).**
> 1. The hypotheses `SLW(k, g, kappa)` of §2.1 are false for every `k >= 1`, `kappa >= 0`, `c in (0,1]` and
>    `g < 1`. For `k = 5` they fail already with `{0, +-1}` coefficients for every `g < 1/4`, and for `k = 4` for every
>    `g < 3/10`. `SLW^log(theta, g, kappa)` is false for every `theta > 0` and `g < 1/2` (Props R1, R2: polynomial
>    templates `prod_l (t/b_l + i)^(sigma_il)` with `sum_l sigma_il b_l^(2j-1) = 0`). Theorems 2.3, 2.4 and
>    Proposition 2.7 are therefore vacuous, and Proposition 2.5 is not sharp.
> 2. With heights weighted by the coefficients, `prod_j Mah(alpha_j)^(-g B_j)`, the proofs of Theorems 2.3–2.4 go
>    through:
>    * `SLW_w(k, g, kappa)` for some `k >= 6` and `g < 1/4` implies the uniform theorem, with `M <= (5.5^(k-1) + 1)/2`
>      for `R >= R_0(C)`;
>    * `SLW_w^log(theta, g, kappa)` with `k_0 >= 6`, `theta < 1/log 5.5` and `g < 1/4` implies `M = O((log R)^(1-eta))`.
>
>    `SLW_w` holds at `g = 1/2` (Liouville). It is false at `k <= 5` for `g < 1/4`. At `k >= 6` it is open, new, and
>    supported by no heuristic, since the Dirichlet heuristic is beaten at `k = 4, 5`. The circle-restricted variant
>    is a restatement of the conclusion with a power margin.
>
> **Unconditional (correct).**
> 3. *Rank reduction* (Theorem 1.3). A normalised configuration on an arc of length `<= C sqrt R` with exponent lattice
>    of rank `k` and `log N' >= 18(2 log^+ C + log 2)` has `M <= ((1 + 2W/rho)^k + 1)/2 <= (5.5^k + 1)/2`.
>    * The uniform theorem is equivalent to bounded rank.
>    * Rank `<= theta log W'` with `theta < 1/log 5.5` would give `M = O((log R)^(1-eta))`.
>
>    No rank bound is known, and the round-5 fakes have full rank.
> 4. *Support lemma* (Lemma 3.2). On the skew double core every nonzero rational character is nonzero on at least
>    `(q+1)/2` cuts. The true minimum is `(q+3)/2` for `q = 7, 11`, and column characters attain `(q+3)/2` for every
>    `q <= 67` tested.
> 5. *LW blindness* (Theorem 3.3, Lemma 3.4). Lang's product-form LW (`eps <= 1`, any constant) is satisfied, with
>    probability `>= 1 - e^(-eps q tau/8)`, by every cluster-determined character form of a random-position angle model
>    on the skew double core. With `exact.md` Theorem 5.1 (unrefereed), single-form LW applied to character forms in
>    the basis `alpha_T`, together with the point-level inputs of `A_1(2M)`, cannot prove growth.
> 6. *Other routes.*
>    * abc is vacuous on squarefree profiles for pair, triangle, Ptolemy and `n`-term relations (Prop. 4.2).
>    * Bombieri–Lang is vacuous on these rational varieties (Prop. 6.2).
>    * `QS(m, delta, f = o(s))` would give growth (Theorem 5.2), and is false for `m <= 3` (Prop. 5.3).
>    * One instance `V(m_0, eps)` with `eps < (m_0 - 6)/2` gives the uniform theorem (Prop. 6.1(a, b)).
>    * The integrality set `S` of the chord data contains every inert prime `<= M - 2` (Props 6.3–6.4).

## 7. Files (`round8/referee_weakest_checks/`)

| script | checks | type | output |
| --- | --- | --- | --- |
| `slw_template_counterexample.py` | Prop. R1: templates for `k = 5` (`{0,+-1}`, `D = 6`; integer, `D = 7, 9`) and `k = 4` (`D = 5`). Exact factorisation, private primes asserted, exact row ranks, angles at 300 digits; `gamma`, `gamma_w`, `gamma_eff(c = 1/2, kappa = 1)` | exact plus high-precision angles | `slw_template_counterexample_out.txt` |
| `slw_log_illustration.py` | Prop. R2: `D = 2s`, `k = s = 3..6`; coverage `k >= theta log(n+2)` at `theta = 1/log 5.5`; `gamma`, `gamma_eff` | same | `slw_log_illustration_out.txt` |
| `weighted_pm1_search.py D B` | §3.1: rank of the `{0,+-1}` solutions of `sum sigma b = sum sigma b^3 = 0` over all `b` in `binom({1..B}, D)` | exact search (float rank on 0/±1 matrices of rank `<= 6`) | `weighted_pm1_search_out.txt` |
| `vertical_fibres.py D B` | §3.2: fibres of vertical *circle* templates (bounded `C`), sizes and affine ranks | exact | `vertical_fibres_out.txt` |
| `skew_support_minimum.py` | Lemma 3.1 lattice; Lemma 3.2 bound (certified for `q = 7, 11` by rank mod `2^61 - 1`); exact minimum `(q+3)/2`; column characters for `q <= 67` | exact | `skew_support_minimum_out.txt` |
| `rank_packing_actual.py X` | Theorem 1.3 steps 1–3 on every window of 2–8 consecutive points, `N <= X` (`X = 20000`) | exact, angles at 50 digits | `rank_packing_actual_out.txt` |

`slw_log_illustration.py` imports `slw_template_counterexample.py` from the same directory.
