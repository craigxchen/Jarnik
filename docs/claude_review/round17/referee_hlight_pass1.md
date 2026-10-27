# Referee report on round17/hlight.md (H_light for sdc_q), adversarial

Scripts are mine and import nothing from `hlight_checks/`. They live in `round17/referee_hlight_checks/`, and each `x.py` / `x.c` writes `x_out.txt`.
The pattern family is built from the record's definition (round5/exact.md Thm 4.1: Paley-I skew-Hadamard `K`, normalised
`DKD`, `S_a = {x : K(x,a) = -1}`, `S'_a = S_a xor {inf, a}`), **not** from the author's `P_a = a - D`, `P'_a = a + D`.

## 0. Verdict

**The result survives.** Theorem 2.8 (H_light for `sdc_q`) is correct, with all quantifiers.
* Every lemma and proposition was re-derived by hand. The only correction is a q = 7 slip in the displayed constant (§2, issue M1).
* My independent exact computations reproduce every number in the write-up: the pattern formula, pair margins, AGL(1,q) invariance, the 4-term minima for all 93 primes `7..1091` and for `3251, 3259, 3271, 5479`, and the dual pairs.
* Weil's bound is used correctly. Its hypotheses are a squarefree monic quartic and a quadratic character.
* The consequence for barrier Lemma 3.4 / Thm 3.5 is stated correctly: H_light is removed, and CIRCLE(q) itself stays undecided.

There are no fatal or major issues. There are four minor ones (§2), and none of them changes the conclusion.

## 1. Re-derivation (all steps checked)

* **Pattern family.** For the Paley-I `K` (`K(x,y) = delta_xy + chi(y-x)`, `K(inf,.) = 1`, `K(.,inf) = -1`, then `DKD`),
  `S'_a = {a} ∪ {x : chi(x-a) = -1}` equals the author's `P_a`, and `complement(S_a)` equals the author's `P'_a`. With the other sign convention the two are swapped.
  For `sum lambda = 0` the value `|lambda(S)|` does not change under complement, so sigma is identical. This was checked
  computationally against the definition (`check1`, `q = 7..83`, 0 mismatches over 720 random integral `lambda` with
  entries up to 3, including `inf`).
* **Lemma 1.1.** `lambda(P_a) = f(a) + (F - f(a) + c(a))/2` follows from `1_QR = (1 - delta_0 + chi)/2`. Correct.
* **Lemma 1.2.**
  * (1) gives `chi*chi = 1 - q delta_0` from `sum_{u != 0} chi(u-1) = -chi(-1) = 1`.
  * (2) follows because the inner sum is `chi(-1)(chi*chi)(t'-t)`.
  * (3) and (4) are immediate.
  * All are correct.
* **Lemma 2.1.**
  * (1): `N >= sum|h| >= (q+1)|F|`.
  * (2): `q f(b) = F - sum_{a != b} c(a) chi(b-a)` gives `||c||_1 >= 2q - |F| + |c(b)| >= 2q - 1`.
  * Both are correct, and the remaining case `|F| <= 1`, `||f||_inf <= 1` is exactly the ternary case.
* **Lemma 2.2.** `|t| = (2/n)t^2 + (2/n)|t|(n/2-|t|)` holds identically, and `sum t_S^2 = (q+1)n/2`. For an integer
  `1 <= |t| <= n/2 - 1` we have `|t|(n/2-|t|) >= n/2 - 1`. Correct.
* **Lemma 2.3.** Euler's criterion holds at `u = 0` because `m >= 1`. The argument has three steps:
  * `(X-t)^m` is expanded with `binom(m,j) != 0 mod q` since `m < q`.
  * The Vandermonde step on the support uses the hypothesis "some `f(t)` is not divisible by `q`".
  * Hence `deg P = m - j_0 >= 0` and `P != 0`.

  Correct. Tightness at pairs is not only exact: it follows from `sum_a chi((a-x)(a-y)) = -1`, giving `(q-3)/2` zeros.
  My tests found 0 violations (`check3`):
  * all `f` with values in `[-3,3]` at `q = 7`;
  * all ternary `f` at `q = 11`;
  * 200k random ternary `f` at `q = 19`;
  * 18k random `f` with values up to `3q` (some divisible by `q`) at `q = 23..103`.
* **Prop 2.4.**
  * The case split on `U + V = h(a) in {1,2}` is right.
  * `phi(n) = q + 2 - n/2 - 2(q+1)/n` is concave, with `phi(6) = 2(q-2)/3` and `phi(q+1) = (q-1)/2`.
  * At `n = 4`, `nu >= q-3`.

  Correct. Exhaustive at `q = 19, 23` (`check5`): the `n = 4` minimum is `(q-1)/2` (attained), and the `n = 6` minimum is `13, 17`, above `(q+1)/2`.
* **Prop 2.5.**
  * On `T`, `nu >= n`. Off `T`, `c` is even and `|c| <= n`.
  * Zeros satisfy `N_0 <= (q-3)/2` (Lemma 2.3, `|T| <= (q+1)/2`).
  * `4M + n^2E <= qn - n`. If the bracket in the bound for `M` is negative, the conclusion is still implied by (1).
  * `g(n) = (q(n-2) + 4n - n^2)/(n+2)`, with `g(4) = q/3` and `g(6) = (q-3)/2`.
  * `(n+2)(g - (q-3)/2) = (1/2)[q(n-6) + 11n - 2n^2 + 6] >= 0` for `q >= 2n+1`.

  All were re-derived and are correct. Exhaustive `n = 6` at `q = 19, 23` gives minima `11, 13`, above `(q-3)/2 = 8, 10`. Sampled `n = 4..22` at `q = 43` shows no violation.
* **Prop 2.6 (dual pairs).** `F = 0`, `||lambda||_1 = (q+1)/2`, `sigma = (q-1)/2`. The first check is from the definition at
  `q = 7..83`, and the second uses the exact formula at `q = 3251, 3259, 5479` (`check4`). Correct.
* **Lemma 2.7.**
  * Walsh expansion: with `x_i = f_i eps_i`, `|sum x_i| = 3/2 + (1/2)e_2(x) - (1/2)x_1x_2x_3x_4`. I solved this from the three orbit values `4, 2, 0`, which gives `|coef| = 1/2` as claimed.
  * `|S_I| <= 1 + 2 = 3` for `|I| = 2`.
  * `|S_T| <= 3 sqrt q`, by Weil for the squarefree monic quartic `prod(X - t_i)`, which is not a constant times a square.
  * Constants: `-6 - 9 + 4 = -11`.

  Correct.
* **Thm 2.8 / Cor 2.9.** The case assembly is correct. Multiples of pairs give `|l|(q+1) - q`, checked exhaustively at `q = 7`
  with entries up to 3: minimum 9 at `|l| = 2`.

## 2. Issues

**M1 (minor; statement slip at q = 7).** The displayed `LB(q) = min(max(q/3, q/2 - 11 - 1.5 sqrt q), (q-3)/2)` equals
`2 < 7/3` at `q = 7`. So "in particular H_light(1/6) holds for every q" does **not** follow from the displayed formula at
`q = 7`, and the §3 table entry "proved LB/r = .1667" at `q = 7` contradicts that formula.
* The claim itself is true. At `q = 7` an `F = 0` ternary character with `n >= 6` has `n > (q+1)/2`, so Prop 2.5(1)
  gives `sigma >= n - 1 >= 5`.
* The author's `c0table.py` uses the true case split, not the displayed min.
* Exhaustively, `min sigma = 3` at `q = 7` (`check2`, entries up to 3).
* Fix: state `LB(q)` for `q >= 11`, and `LB(7) = 3`. Alternatively, replace `(q-3)/2` by "`(q-3)/2`, or `n - 1 >= (q+1)/2` when `n > (q+1)/2`". Either way `(q-3)/2 >= q/3` iff `q >= 9`.

**M2 (minor; overstatement).** §0 says "the base 4-terms `e_x + e_y - e_z - e_inf` have `sigma = (q-1)/2` exactly". This is
false as a statement about all base 4-terms: they take exactly two values, `(q-1)/2` and `(q+3)/2` (`check2`, all
`q = 7..47`). Prop 2.4's "attained by `e_0 + e_1 - e_z - e_inf` for a suitable `z`" is the correct statement.

**M3 (minor; citation status).** Weil's bound is quoted with its true hypothesis: `g` is not a constant times a square, and `d` counts distinct roots.
* I could not locate the theorem number by search. My searches returned only neighbouring citations, e.g. [LN97, Thm 5.38] for a different Weil sum, so the label should stay [M]/[L-snippet].
* [M]: the standard reference is Lidl–Niederreiter, *Finite Fields*, Thm 5.41.
* [M] remark: for a squarefree quartic, Hasse's genus-1 bound gives `|S_T| <= 2 sqrt q + 1`. That would improve Lemma 2.7 to `q/2 - 10.5 - sqrt q`, which is `0.2396 r` at `q = 3251`. This is optional.
* Independent safeguard: at `q = 3251, 3259, 3271, 5479`, my exhaustive enumeration verifies `sigma >= q/2 - 11 - 1.5 sqrt q` for **every** 4-term `F = 0` character directly (0 violations, integer test). So `c_0 >= 0.2349` at these regime primes holds even without Weil. For general regime `q` it rests on Weil, a classical theorem.

**M4 (minor; incomplete threshold list).** §4 item 2 lists `tau >= max(tau_D, tau_2, 4)`. Barrier Thm 2.3 also requires
`sigma*(tau) <= r + 1` and `||lambda||_1 Delta <= 1/2` for light `lambda`. §0 says "thresholds of barrier Thm 2.3", which
covers them; item 2 should say the same. All of these are `O(log q)`, so the barrier scale is unaffected.

**Points checked and found correct:**
* **Regime.** `q >= 3251` is the first sdc regime prime (`d = q >= 316 ln 8q`: `316 ln 26008 = 3212.4`). `w = sqrt(8 q ln 8q)`. `m' = 2057 < q`.
* **Normalisation.** The `m < N` relaxation is used, not `2m' < q`. The cap threshold `eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})` and `W = 2q tau` are carried over unchanged.
* **`tau_D`.** It is Lemma 3.4's formula with `r = 2q`, `n = q+1`, `2/c_0 = 12`. For `c_0 = 0.2349` the factor is `2/c_0 = 8.514 <= 8.52`.
* **`kappa_H` / `kappa_V`.** The write-up consistently uses `kappa_H` for the cap test and touches neither side of `kappa_H <= kappa_V`.
* **Relevance.** All relevant `lambda` are integral for `k = e_1`, so Thm 2.8 covers exactly H_light's domain, and more (light and heavy).
* **Prop 2.5.** It needs `o > 1/(2 c_0)` together with `K >= K_0(c_0, o)`. That gives `o >= 4` with `c_0 = 1/6`, and `o >= 3` with `c_0 = 0.2349 > 1/6`. The referee's `tau_L` fix is retained.
* **Evidence labels.** Evidence is labelled as evidence (`anneal.py`). The exact minima are correctly stated as unused in any proof.

## 3. Independent computations

| check | result | file |
| --- | --- | --- |
| sdc_q from exact.md Thm 4.1 (Paley-I): AGL(1,q) invariance, pair `sigma` | invariant; `sigma_xy = 1` for all pairs incl. `inf`; `q = 7..83` | `check1_out.txt` |
| Lemma 1.1 formula vs definition | 0 mismatches / 720 random `lambda` (`|lambda_x| <= 3`, incl. `inf`, 6-term ternary) | `check1_out.txt` |
| exhaustive global min, `q = 7` (`|lambda_x| <= 3`), `q = 11` (ternary) | 3, 5 (`= (q-1)/2`); multiples of pairs min 9 | `check2_out.txt` |
| base 4-terms `e_x+e_y-e_z-e_inf` | values `{(q-1)/2, (q+3)/2}`, `q = 7..47` (M2) | `check2_out.txt` |
| Lemma 2.3 | 0 violations (exhaustive `q = 7, 11`; random `q = 19..103`, values up to `3q`) | `check3_out.txt` |
| 4-term `F = 0` minima, all 93 primes `7..1091` + `3251, 3259, 3271, 5479` (C, AGL normal form `d_0 + d_1 - d_c - d_e`) | identical to the author's table, e.g. `1571 (0.2416 r)` at `3251`, `2667` at `5479`; `< (q-3)/2` for all `q >= 19`; `q = 7, 11` not below; 0 violations of the Lemma 2.7 bound | `r4_out.txt`, `r4_all_out.txt`, `r4_all_summary.txt` |
| AGL normal form vs brute force (no symmetry) | equal: 7, 9, 13 at `q = 19, 23, 31` | `check4_out.txt` |
| dual pairs | `sigma = (q-1)/2`, `||lambda||_1 = (q+1)/2`, `q = 7..83` (definition), 3251, 3259, 5479 (formula) | `check4_out.txt` |
| `|F| = 1` (`n = 4, 6`), `F = 0` (`n = 6`) exhaustive at `q = 19, 23`; `F = 0`, `n = 4..22` sampled at `q = 43` | all above the proved bounds | `check5_out.txt` |

## 4. Corrected statement

**Theorem (H_light for sdc_q; proved).** Let `q = 3 mod 4` be prime with `q >= 7`, and `r = 2q`. Let `lambda` be integral with `sum lambda = 0` and not an
integer multiple of a pair. Then `sigma_lambda >= ceil(q/3)`, so H_light(1/6) holds for every such `q`; at `q = 7`
this uses Prop 2.5(1) for `n >= 6`. For `q >= 11`, moreover,
`sigma_lambda >= LB(q) = min(max(q/3, q/2 - 11 - 1.5 sqrt q), (q-3)/2)`.
* The bound `q/2 - 11 - 1.5 sqrt q` uses Weil's bound for a squarefree quartic ([M]: Lidl–Niederreiter Thm 5.41). At `q = 3251, 3259, 3271, 5479` it is also verified directly by exhaustive exact enumeration.
* Hence `c_0 >= 0.2349` at every regime `q >= 3251`, and `c_0(q) >= 1/4 - O(q^{-1/2})`.
* Dual pairs give `sigma = (q-1)/2` for every `q`, so `c_0 <= (q-1)/(4q)`.
* Exact 4-term minima give a sharper upper bound, e.g. `0.2416` at `q = 3251`.
* Some base 4-terms have `sigma = (q-1)/2`; others have `(q+3)/2`.

**Consequences.**
* Barrier Lemma 3.4 holds unconditionally, with `tau_D = max(12(2 ln(8q+5) + ln(1000/C)), 4 ln(8q+5) + 8, 2 ln(1000 q (q+1)^2/C))`.
* Hence, with the referee's corrections (`m < N`, AGL(1,q)), barrier Thm 3.5 is the following unconditional cap test. Suppose that at one regime prime `q >= 3251` and `W = 2q tau`, `tau = O(log q)` above `tau_D` and **all** thresholds of barrier Thm 2.3 (`tau_2`, `4`, `sigma* <= r+1`, `||lambda||_1 Delta <= 1/2`), we have `inf_theta kappa_H(hat u_theta) <= eps_0/(2 sqrt(m'+1) e^{pi^2/(4q)})`. Then route 5 is capped.
* CIRCLE(q) is **not** decided. Growth and uniformity are unchanged.
