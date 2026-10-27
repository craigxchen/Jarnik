# Referee report on `round17/survive.md` (CIRCLE(q), survival side)

Referee scripts: `round17/referee_survive_checks/`. Each `x.py` writes `x_out.txt`, except `toy.py`, which writes
`toy_<q>_<L>_out.txt`. All code is independent of `survive_checks/`. The record (`docs/claude_review/`) was used
read-only.

## 0. Verdict

**The result survives, with corrections. None of them is fatal.**

* The main new theorems are correct as stated: Lemma 2.1, the inequality in (2.2), Theorem 3.1, Theorem 3.4,
  Corollaries 3.2 and 3.3, and Propositions 4.1 and 4.2. I re-derived every step and checked them by independent
  exact computation.
* CIRCLE(q) is **not decided**. The write-up says so correctly, and no growth claim is made.

Corrections:

1. **Major (Remark 4.4; the number is false, the conclusion holds after repair).**
   * The claim "`min_p G = +1823`" is false under the author's own formula. Their grid `p = 10^{-e/20}` never goes
     above `p = 0.891`, and with the stated `lambda = p/2` one gets `G(p) -> -infinity` as `p -> 1`. For example,
     `G(0.999) = -1487` and `G(1 - 10^-6) = -6213`. So the computation as specified does **not** show that the chain
     fails for every `p`.
   * The formula is also wrong. The smallest atom of the `p`-biased law is `min(p/2, 1 - p)`, not `p/2`.
   * With the correct atom, an interval-partition bound gives `G(p) >= 2848 > -ln 2` on all of `(0,1)`. So the
     conclusion holds (`g44_cert_out.txt`).
2. **Minor (the budget depends on `tau_D`).**
   * "Budget `142.5 q`" rests on `tau_D(c_0 = 0.184)`. That value of `c_0` is evidence only (H_light is unproved),
     and `C = 1` is assumed.
   * A robust form exists. The third term of `tau_D` does not depend on `c_0`, so for `C = 1` and **every** `c_0`
     one has `tau_D >= 62.34` and budget `>= 25.6 q`.
   * The budget is an artefact of the crude union bound in Lemma 3.4. If a sharper Diophantine lemma allowed the
     circle cap at `tau <~ 13.3`, the budget would drop below `q ln 3`. At `tau <= 11.1` the survival criterion
     (3.2) becomes empty.
3. **Minor (equation (2.2)).** It writes `Cheb^C_A(M_i) = sup{|<c,M_i>|/|c|_1 : c real}`. Only `>=` holds in general.
   Counterexample: three points, `f` = the cube roots of unity, constants. Then `Cheb^C = 1`, while the supremum over
   real `c` is `sqrt 3/2`. The text then says "`<=` is not claimed", so the display should be an inequality.
   Nothing downstream uses the equality.
4. **Minor (the calibration table, §6.3).**
   * The reported "factor 7 to 33" compares the LP value with the bound **before** using `max_s |sigma_s alpha| >= 1`
     and `binom >= 1`. Corollary 3.2 uses the bound **after** those steps, and on the toys that bound is up to about
     `5 * 10^4` below the LP value. For example, `(7,2)` at `m = 7` gives LP `62.76` against `(3.1)` bound `0.0013`.
   * The factors also depend on which witness is chosen. With shorter witnesses I get factors of 21 to 47 (§7).
   * The row `(7,2)` lists `m = 1..6` as the degrees with `dim V > 0`, but `m = 7, 8` also have `dim V > 0`
     (44 and 1).
5. **Minor (wording).**
   * "Counting cannot produce it" (Remark 4.3) is a statement about one counting argument (a Bombieri-Vaaler basis
     plus pigeonhole on the codimension of `V'`). It is not a theorem about all counting methods.
   * "In practice at about `4j`" (§7) rests on sampled `h'` data, but the text calls it "exact cost data".
   * "A witness of degree between `8jq/(q+1)` and `8j+8`": a witness has order exactly `m'`. The bracket only locates
     `m'`.
   * Corollary 3.3(1), "at every barrier `tau`", should read "at every `tau >= tau_D` meeting the thresholds of
     barrier Thm 2.3".

## 1. Setup checks (regime, patterns, normalisations)

* **The patterns are the record's patterns** (exact; `pats.py`).
  * Construction: I rebuilt `sdc_q` from the record's definition. That is the Paley-I skew-Hadamard matrix,
    normalised on column 0, with cuts `S_a`, `S'_a = S_a xor {0,a}` and relative pattern = the side avoiding the
    base row (round12 `flib.py` conventions).
  * Result: this set equals survive's `{P_a, P'_a}` for all 39 primes `q = 3 mod 4` with `7 <= q <= 383`. There are
    `2q` patterns of size `(q+1)/2` each.
  * `AGL(1,q)` maps the pattern set onto itself for every such `q <= 83`, checked exactly over all `q(q-1)` maps.
* **Regime arithmetic** (50 digits; `regime.py`). `q = 3251` is prime, `q = 3 mod 4`, and it is the first regime
  prime: `316 ln 4r = 3212.51 <= d = q`. Also `w = 514.1998`, `j = 257`, `m' = 2057 < q`,
  `8j = 2056 <= m' <= 8j + 8 = 2064`.
* **Thresholds.** These match barrier.md §1, §3.3 and README-16 §1: `d = q`, `r = 2q`, `W = 2q tau`,
  `eps_0 = C q e^{-q tau/2}/(m'+1)` and `ln H' = m' q tau/20`.
  * `kappa_c = |Phi_c|/(|c|_1 d^m)` matches the record. So does `V_H` (integral `c` in `V(A,m')` with
    `|c|_inf <= H`).
  * `kappa(hat u) = kappa(u)/|u|_inf^m >= kappa(u)` holds because `Phi_c` is homogeneous of degree `m` and
    `|u_theta|_inf <= 1`.
* **`tau_D` and budget.** Lemma 3.4's formula with `c_0 = 0.184` and `C = 1` gives `tau_D = 296.09` and budget
  `463236.8 = 142.49 q`. The author's numbers are reproduced. Further values:

  | `c_0` | `tau_D` | budget/`q` |
  | --- | --- | --- |
  | 0.25 | 217.9 | 103.4 |
  | 0.5 | 109.0 | 48.9 |
  | `>= 0.87` | 62.34 | 25.6 |

  At the record's `W_min = 3n ln n - 14n` (`tau = 5.13`) the budget would be `-9718`. That is irrelevant only
  because Theorem 3.5 needs `tau >= tau_D`.
* **`kappa_V` versus `kappa_H`.** This is handled correctly. Theorem 3.1 bounds `kappa_H` for every
  `H >= |c|_inf`, since every certificate used is an `AGL`-image of one integral `c` with the same `|.|_inf`.
  Corollary 3.2 needs `ln |c|_inf <= ln|c|_1 < q tau/2 < ln H' = m' q tau/20`, which holds because `m' >= 10`.
  §5.2 correctly declines to claim `kappa_H >~ kappa_V` on `C_1`.

## 2. Lemma 2.1 and the isotypic reduction

I re-derived each identity:

* `s P_a = P_{sa}` or `P'_{sa}` according to `chi(s)`;
* `(tau k)_x = k_{x-1}` gives `l_s(tau k) = zeta^s l_s(k)`;
* `l_1(mu_s k) = l_s(k)` and `l_{-1}(mu_s k) = l_{-s}(k)`;
* `<g c, f> = sum c_k f(g k)`.

All correct. The exponents `2i - m` with `0 <= i <= m < q` are distinct mod `q`, since `q | 2(i-l)` forces `i = l`.
So the projection (2.1) is correct.

In the proof of (2.2), the identity `l_{u_theta} o tau^{-a} = l_{u_{theta - 2 pi a/q}}` holds (re-derived), and the
`theta`-twisted average gives `2^{-m} binom(m,i) M_i` exactly. So
`Cheb_A(u_theta) >= 2^{-m} binom(m,i) Cheb^C_A(M_i)` holds. The equality sign in the display is not justified
(item 3 of §0).

## 3. Theorem 3.1 (equivariant witness bound): correct

* **Item 1.** `alpha` lies in `Z[zeta]` and is nonzero, so `|N(alpha)| >= 1`. The Galois group is all of
  `F_q^*`, so the largest of the `q-1` conjugates has modulus `>= 1`.
* **Item 2.** `<mu_s c, M_i> = sigma_s(alpha)`, because `c` is integral and `M_i(mu_s k) = sigma_s(M_i(k))`.
  Write `f_a(theta) = Phi_{tau^a mu_s c}(u_theta)`. The twisted mean of the `f_a` equals
  `2^{-m} binom(m,i) e^{i theta(2i-m)} sigma_s(alpha)`. A mean of `q` numbers has modulus at most their largest
  modulus, so some `|f_a(theta)|` reaches the bound. This works for every `theta` and needs no Diophantine input.
* **Item 3.** `c_*` lies in `V_H` and has the same `|.|_1`.

So the referee obstruction `e^{-Theta(q^2)}` is genuinely removed **for a given witness**. The loss is now
`|c|_1^{-1}` instead of roughly `(|c|_1 q^m)^{-(q-3)/2}`. It is not removed for exact non-membership with no height
control: the minimal witness height is then unbounded except by Cramer/Hadamard-type bounds (Remark 4.3's scale).
The write-up says this correctly in §5.2 and in Corollary 3.3(3).

Independent exact checks (`toy.py`):
* integral `c` from an exact `fmpz` nullspace plus LLL;
* `alpha_i` exact in `Z[x]/(x^q - 1)`, with `alpha = 0` iff all coefficients are equal;
* conjugate moduli at 40 digits;
* LP Chebyshev values from HiGHS (evidence).

| `(q, L)` | `m` | `dim V` | witness `|c|_1` | `max_s |sigma_s alpha|` | LP Cheb (min over 4 `theta`) | (3.1) bound | item-2 check, `min_theta` ratio |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (7,1) | 1, 2, 3 | 56, 49, 21 | 2, 4, 768 | 1, 3.05, 1398.6 | 1.90, 1.81, 1.71 | 0.25, 0.125, 4.9e-4 | 1.95, 2.16, 2.51 |
| (7,1) | 4 | 0 | — | — | — | — | — |
| (7,2) | 2, 4, 5, 6, 7 | 749, 644, 483, 232, 44 | 4, 12, 24, 72, 126 | 3.80, 34.2, 170.1, 2186.6, 9391.3 | 9.69, 46.2, 81.3, 111.6, 62.8 | 0.125, 0.031, 0.013, 0.0043, 0.0013 | 1.93, 2.61, 3.38, 3.09, 3.74 |
| (7,2) | 8 | 1 | 1008 (`|c|_inf = 42`) | 28224 (rational) | 7.65625 (all `theta`) | — | 1.000000 |
| (11,1) | 1, 2, 3 | 132, 121, 55 | 2, 4, 6 | 1, 3.60, 11.56 | 1.92, 1.88, 1.88 | 0.25, 0.125, 0.0625 | 1.98, 2.06, 2.60 |
| (11,1) | 4 | 0 | — | — | — | — | — |

* **Item-2 check.** For every case and each of 61 values of `theta`,
  `max_a |Phi_{tau^a mu_s c}(u_theta)| / (2^{-m} binom |sigma_s alpha|) >= 1`. This confirms the uniformity mechanism.
* **Bound against LP.** The LP value always exceeds the bound.
* **The sdc7 datum** (`L = 2`, `m = 8`) is reproduced exactly: `dim V = 1`, `|c|_1 = 1008`, `|c|_inf = 42`,
  `<c, M_4> = +-28224 in Z`, and `|Phi|/|c|_1 = 7.65625` (`m = 8 > q` here; only `i = 4` pairs).
* **Witness choice matters.** My witnesses at `(7,2)`, `m = 2, 4, 5, 6`, are shorter than the author's
  (`4, 12, 24, 72` against `6, 16, 32, 128`). So the author's "shortest" witnesses are not the shortest. With mine,
  the "LP/bound-before-`max >= 1`" ratio is 21, 47, 41, 13, not the reported 5.7, 5.6, 4.9, 4.9. This is evidence
  only; see item 4 of §0.

## 4. Theorem 3.4 and Corollaries 3.2 and 3.3: correct

* **Theorem 3.4.** `V_H` is spanned by integral vectors of height `<= H`. If `Phi_c` restricted to `C_1` is nonzero
  for some `c` in `V_H`, then it is nonzero for some spanning vector `c_0`, which has `|c_0|_1 <= |A| H`. Then
  apply Theorem 3.1. The two cases are mutually exclusive, since the bound is positive.
* **Corollary 3.2.** Rearranging (3.2) gives `2^{-m'} q^{-m'}/lambda_q > C q e^{-q tau/2}/(m'+1) = eps_0`, and
  `eps_0 > T_cap`. The claim `lambda_q < infinity iff C_1 is not in TC_V` holds because `V(A,m)` is defined over
  `Q`.
* **The window.** At `q = 3251` and `tau = tau_D(0.184)`, the window `[ln H_*, ln H'] = [~4.60e5, 9.90e7]` is
  reproduced. The `-ln|A|` term costs `<= q ln 3 = 3572`.
* **Scope.** Corollary 3.2 shows that the **circle test** does not cap route 5. It does not show that route 5
  survives, since other barrier points of `G` are not excluded. The text says "circle test is passed",
  consistently.

## 5. Propositions 4.1 and 4.2: correct

**Lemma 4.0.**
* `P'_a = a + ({0} u QR)` because `chi(-1) = -1`.
* `|hat 1_{D'}(b)| = |1 + chi(b) i sqrt q|/2 = sqrt(q+1)/2`, using the Gauss sum `g = i sqrt q`. This is valid for
  `q = 3 mod 4` prime, which is the hypothesis in force.
* Parseval then gives `(q+1)(K^2 + |v|_2^2)/4`.
* Exact check (`parseval.py`): 0 failures of the identity and of `h' >= (q+1)|v|_2^2/(2|v|_1)` over 2000 random
  vectors each at `q = 7, 11, 19, 23, 43`, and 300 at `q = 3251`. Also `h'(e_x) = q+1`.

**Proposition 4.1.**
* The extreme-point additivity is valid: a generic `lambda` has unique maximisers, and the coefficient at the
  maximiser of the convolution is `prod_l c_l(k_l^+)`.
* On `A`, `range_P <= 2j`. Each coordinate lies in exactly `q+1` patterns.
* Type (b) satisfies `m_l <= sum_x (|T_x| - 1) <= sum_x range(T_x)`.
* Type (a) satisfies `sum_P |v(P)| >= (q+1)/2`.
* Hence `(q+1) M/2 <= 4jq`, and `8jq/(q+1) < 8j <= 4w <= m'`.
* Exact orders multiply to a nonzero leading form, so in fact no product-type element of `V(A,m')` exists at all.
  The text claims somewhat less than this.

**Proposition 4.2.**
* `|Y| <= 2j` follows from choosing `a` not in `Y`, using `P_a u P'_a = F_q` and `P_a n P'_a = {a}`.
* The monomial-basis argument (minimal support, characters on `{+-1}^Z`) is correct.

**Remark 4.3.**
* `ln r_0 >= ln sum_{s <= j} binom(q,s) 2^s = 1073.27`. This is reproduced, and the surjectivity it rests on follows
  from Proposition 4.2.
* "Cannot" refers to one argument only (item 5 of §0).
* The Fukshansky citation exists: Monatsh. Math. 147 (2006) 25-41, arXiv math/0409374. Its abstract describes
  small-height integral points outside a hypersurface, with bounds in degree and number of variables. Abstract-level
  use is labelled.

**§5.3 (Sylvester transfer).** The conclusion also follows without the quoted normaliser fact: maps `x -> x^s` on
`F_{2^k}^*` are `F_2`-linear only for `s` a power of 2. So only Frobenius conjugates are realised, as stated.

## 6. Remark 4.4 (the error)

`biased_l2.py` (author) uses the grid `p = 10^{-e/20}`, `e = 1..399`, that is `p <= 0.891`, with
`lambda = p/2`. My recomputation (`g44.py`, `g44_cert.py`):

* **Author's formula.** `G(0.891) = 1825`, `G(0.99) = 96`, `G(0.999) = -1487`, `G(1 - 10^-6) = -6213`, and
  `G -> -infinity` as `p -> 1`. The reason: the tail lower bound `(1-p)^{|P|-j-1}` goes to 0, while the
  hypercontractive term stays bounded. So "`min_p G = +1823`, chain fails for every `p`" does not follow from the
  stated computation.
* **Corrected formula.** The minimal atom is `min(p/2, 1-p)`, and it is `1 - p` for `p > 2/3`.
  * A certified interval-partition lower bound on `[1e-30, 1 - 1e-30]` is `G >= 2848.26`, with the minimum near
    `p = 2/3`.
  * The two ends are handled analytically. The exponent margins are `m' - j - 1 > 0` at `p -> 0` and
    `m' - (|P| - j - 1) = 689 > 0` at `p -> 1`.
  * So `G > -ln 2` for every `p` in `(0,1)`, and the remark's conclusion is restored.
* **Scope.** Even corrected, the remark only says that this chain, with worst-case `(2,4)` hypercontractivity and
  this tail lower bound, cannot close. It is not a statement about `kappa`. The author's label says the same.

## 7. Evidence versus proof

* Labels are mostly accurate. The toys are labelled evidence, and the LP values floating point.
* Two exceptions:
  * "in practice at about `4j`" is presented as following from "exact cost data", but those data are sampled for
    support 4 to 8;
  * the calibration factor quoted in the summary ("bound within 7-33 of LP") does not measure the loss of the bound
    used in Corollary 3.2 (item 4 of §0).
* No growth or uniformity claim is made. The growth rate stays `M <= (2/3+eps) log R/loglog R`.

## 8. Corrected statement

The following is proved for `sdc_q`, with `q = 3 mod 4` prime, `A = A_q(L)` and `1 <= m < q`.

* **(Theorem 3.1).** Let `c` be in `V(A,m) n Z^A` with `alpha = <c, hat k(1)^i hat k(-1)^{m-i}> != 0`. Then for
  every `theta` and every `H >= |c|_inf`:
  `kappa_H(hat u_theta) >= 2^{-m} binom(m,i) max_s |sigma_s alpha| / (|c|_1 q^m) >= 2^{-m} binom(m,i)/(|c|_1 q^m)`.
* **(Theorem 3.4).** On the Fourier circles, `kappa_H` is either identically 0 or at least `2^{-m}/(|A| H q^m)`.
* **(Corollary 3.2) at regime parameters.** These are `L = j`, `m = m'`, `W = 2q tau`. At `q = 3251`: `j = 257`,
  `m' = 2057 < q`. If an integral witness has `ln|c|_1 < q tau/2 - m' ln 2q - ln(Cq/(m'+1))`, then
  `inf_theta kappa_{H'}(hat u_theta) > eps_0 > T_cap` on every `C_b`. So the circle test does not cap route 5 at
  that `tau`.
* **The budget.** At `tau = tau_D` (Lemma 3.4, `C = 1`) the budget is `>= 25.6 q` for every `c_0`, and `142.5 q` if
  `c_0 = 0.184` (evidence value). It falls below `q ln 3` only if the cap could be run at `tau < 13.3`, which needs
  a sharper Lemma 3.4.
* **(Propositions 4.1 and 4.2).** No product-type measure (zonotope, sub-cube tensor, or convolution of these) lies
  in `V(A_q(j), m')` at all. Nor does any measure supported on a union of coordinate sub-cubes inside `A`.
* **(Remark 4.4), corrected.** The minimal atom is `min(p/2, 1-p)`, which gives `G >= 2848` on `(0,1)`. The
  `p`-biased `L^2` / hypercontractive chain cannot close for any `p`. The author's "+1823" is a grid artefact of a
  formula whose infimum is `-infinity`.

**Not proved:**
* the existence of any regime witness (W_q);
* the cap `C_1 in TC_{H'}`.

CIRCLE(q) is undecided. Per the round rule, route 5 stops unless the cap side succeeds.

Sources: [arXiv math/0409374](https://arxiv.org/pdf/math/0409374),
[Claremont: Integral Points of Small Height Outside of a Hypersurface](https://scholarship.claremont.edu/cmc_fac_pub/58).
