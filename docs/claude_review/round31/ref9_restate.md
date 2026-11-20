# Referee report, Part IX, block "restate": IX.1 and the last two paragraphs of IX.6

Refereed text: `astra_free_IX/answer.md`, lines 1961-1978 (IX.1), 2239-2257 (the rational-direction fake) and 2263
(what the Gram matrix still permits). Sources checked against: `round31/README_VIII.md` §§1.1, 1.3, 1.7, 1.8;
`round5/construct.md` §§1.1-1.2, 2, 5 (Theorem 2 and its proof); `astra_free_IX/reply9.txt` ((G1)-(G4)); IX.3 of
the answer (lines 2042-2053) for the pair-surplus estimate. Scripts (my own, exact integer arithmetic unless
stated) are in `round31/ref9_restate/`.

## 0. Verdicts in brief

| Item | Verdict |
|---|---|
| IX.1, restated VIII.1 with the `(1/2) log 2` term | **Proved** (matches README_VIII §1.1 and §1.8(a) word for word in substance) |
| IX.1, new: `M <= 2 + s` for `C' < 2^{3/4}` (block mod 2) | **Proved.** Classes, collision gain and threshold all correct. It is not new: README_VIII §1.1 already states it. It is also not the best this method gives (§1.2: `C' < 2` at affine cost 3, and a 2-adic tower up to `(5/4) log 2`) |
| IX.1, radical corollary `M <= (2/kappa + o_C(1)) log R/loglog R`, uniform in `kappa` | **Proved** |
| IX.1, (VIII.8) with quantifiers | **Proved**, consistent with README_VIII §1.3 (the judge's corrected wording) |
| IX.6, growth-scale `A_1(2m)` fake in the direction window of `b = 2+i` | **Proved, with corrections**: `K` and `a` depend on `b` (absolute only because `b` is fixed); all extra monomials have `n` even; the scope must be stated precisely (§2.4) |
| IX.6, "what the Gram matrix still permits" | **Proved for IX.3.** Novelty is minimal: the claim about the mean follows from `W = O(m log m)` and the all-ones vector alone. The max claim needs the two-sided pair estimate. That estimate is checked for IX.3; for IX.5 it is inherited from IX.5 and not rechecked here. The remaining content is correctly identified as the Gaussian-integer cocycle and integrality |

No step is false. The most serious point is one of interpretation. The rational-direction fake shows only that the
integer projection inequality (IX.1) is not a consequence of the `A_1(2m)` axioms. It does not show that (IX.1) is
independent of the stronger P3-char axioms of IX.3. It also says nothing about genuine clusters outside rational
direction windows of small height (§2.4).

## 1. IX.1

### 1.1 Restated VIII.1 (Proved)

**Statement (as in README_VIII §1.1, which Astra reproduces).** Let `M >= 2` *distinct* normalised marks
`z_i'` lie on `|z|^2 = N'`. Let `C' = max_{i != j} |z_i' - z_j'| / N'^{1/4} <= C`, `s = omega(N')`, and `S` a finite
set of primes `l = 3 mod 4` with `(1/2) log 2 + sum_{l in S} log l/(l+1) > log C'`. Append to the exponent rows the
one-hot indicators of the class of `z_i' mod l` in the norm fibre of size `l+1`. Then the enlarged rows are
affinely independent, and `M <= 1 + s + K` with `K = sum_{l in S} l`.

Each step of Astra's sketch is the corresponding step of the judge's check (README_VIII §1.8(a), items 1-9). The
facts used are: `E|X-Y| <= 2u(1-u) <= 1/2`; the collision probability is at least `1/(l+1)` by Cauchy-Schwarz,
because an affine dependence forces equal residue laws; colliding inert primes divide `b_ij` and not `g_ij`; `1+i`
divides every `b_ij` because `N'` is odd and `N(g_ij) | N'`; a block of `l+1` classes has affine cost `<= l`.
Distinctness of the marks is needed and is assumed ("distinct marks").

### 1.2 The block modulo 2, and the 2-adic tower (Proved; Computed)

**The classes.** `Z[i]/(2)` has four elements, `0, 1, i, 1+i`. An odd-norm Gaussian integer `x + iy` has `x + y`
odd, so it is `= 1 mod 2` (`x` odd, `y` even) or `= i mod 2` (`x` even, `y` odd). So there are exactly two classes,
as Astra says. Modulo `(1+i)^3` there are four unit classes `+-1, +-i`.

**The gain.** Equal residue laws for the two draws give collision probability `P(z_i' = z_j' mod 2) >= 1/2` by
Cauchy-Schwarz on two classes. A collision means `2 = -i(1+i)^2` divides `z_i' - z_j'`. Since `g_ij` has odd norm,
`v_{1+i}(b_ij) >= 2`. So `E v_{1+i}(b_ij) >= 1 + 1/2`. Since `log|1+i| = (1/2) log 2`, the gain is
`(1/2)(1/2) log 2 = (log 2)/4`. With `S` empty, the hypothesis becomes `(3/4) log 2 > log C'`, that is
`C' < 2^{3/4} = 1.68179...`, exactly as Astra says. The block has two classes, so its affine cost is one, and
`M <= 1 + s + 1 = 2 + s`. **Proved.**

This refinement is not new: README_VIII §1.1 already states "a block `z' mod 2` (classes `1, i`, affine cost 1)
gives `M <= omega(N') + 2` for `C' < 2^{3/4}`." Astra's sentence restates it.

**Referee addition: the tower (Proved).** The same argument works modulo `(1+i)^k` for every `k`. Let `f_k` be the
number of unit classes mod `(1+i)^k` that contain a Gaussian integer of norm `n`, for a fixed odd `n = 1 mod 4`.
The classes of all marks lie in such a fibre, since `N(z)` mod a suitable power of 2 is determined by the class.
Then `E v_{1+i}(b_ij) >= sum_{j <= k} 1/f_j`, and the enlarged rows (exponents plus a one-hot block of `f_k`
classes mod `(1+i)^k` plus the blocks at `S`) are affinely independent whenever

`(1/2) log 2 * sum_{j<=k} 1/f_j + sum_{l in S} log l/(l+1) > log C'`, at extra affine cost `f_k - 1`.

**Computed** (`two_adic.py`, brute force over residues): `f_1, ..., f_8 = 1, 2, 4, 4, 8, 8, 16, 16`. Hence:

| `k` | `sum_{j<=k} 1/f_j` | threshold for `C'` (S empty) | bound |
|---|---|---|---|
| 1 | 1 | `sqrt 2` | `M <= 1 + s` |
| 2 | 3/2 | `2^{3/4} = 1.682` | `M <= 2 + s` (Astra) |
| 4 | 2 | `2` | `M <= 4 + s` |
| 6 | 9/4 | `2^{9/8} = 2.181` | `M <= 8 + s` |
| `infinity` | 5/2 | `2^{5/4} = 2.378` | cost unbounded |

The mod-4 block (`k = 4`) gives `M <= 4 + s` for every `C' < 2`, at the same affine cost as `k = 3`. This matches the
observed first bare-row deficiency at chord constant exactly 2 (README_VIII §1.7). The 2-adic part of `log|b_ij|`
can contribute at most `(5/4) log 2` in this argument.

**Computed** (`mod2block2.py`, exact ranks over Q): every circle `N <= 10^6` with at least 6 lattice points, and
every window of `>= 3` consecutive points of arc constant `<= 4`, was normalised and tested. That is 108,320 windows.
For each `k <= 6` and every window with `C'` below the `k`-threshold, the rows enlarged by the mod-`(1+i)^k` block
were tested for affine independence, and `M <= s + f_k` was checked. **Failures: 0.** The numbers tested at
`k = 1, ..., 6` were 7,944, 17,936, 26,984, 35,464, 39,808 and 43,560. The bare exponent rows are dependent in 16
windows with `2 <= C' < 2.181`, for example `N = 25, 50, 100, 125`. In all 16 the 2-adic blocks at `k = 5, 6` restore
independence. So the tower does nontrivial work in the data.

### 1.3 Radical corollary (Proved)

**Statement.** `M <= 1 + K_C + omega_1(rad N)`. With `R = sqrt N` and `kappa = log N / log rad N`,
`M <= (2/kappa + o_C(1)) log R/loglog R` as `N -> infinity`, with the `o_C(1)` uniform in `kappa`.

*Check.* `s = omega(N') <= omega_1(N)`, because `N' | N` and `N'` is split-only. The `k`-th prime is `>= k+1`, so `s`
distinct primes have product `>= (s+1)!`. Put `U = log rad N = log N / kappa`, so `log (s+1)! <= U`.
* If `U >= log N/(loglog N)^2`, then `U -> infinity` and `log U = loglog N + O(logloglog N)`. Both hold uniformly in
  `kappa`, because they depend only on `N`. Stirling gives `s <= (1 + o(1)) U/log U = (1 + o(1)) log N/(kappa loglog N)`.
  Since `log N = 2 log R`, `loglog N ~ loglog R` and `kappa >= 1`, this is `(2/kappa + o(1)) log R/loglog R`.
* If `U < log N/(loglog N)^2`, then `s <= U/log 5 = o(log R/loglog R)`, which is again within the bound.

`1 + K_C` is `o_C(log R/loglog R)`. The corollary beats `2/3` only when `kappa > 3` (README_VIII §1.7).

### 1.4 (VIII.8) with explicit quantifiers (Proved; consistent with README_VIII §1.3)

Astra's quantifiers are as follows. For each fixed `rho` in `(0.99, 1)` and all `q >= q_0(rho)` she gives
`m = (rho/2 + o_rho(1)) q`, `P_+ <= 2 q^{1.99} <= m^{1.995}` and `W <= (7.96/rho + o_rho(1)) m log m`. Along a
diagonal `rho_j -> 1`, with `q_j` chosen for each `j`, the coefficient is `7.96 + o(1)`. No uniform starting
threshold in `rho` is claimed. This is exactly the judge's corrected wording.

Arithmetic check: `W <= 2(q-1) log(2q^{1.99}) = (3.98 + o(1)) q log q`, and `m log m = (rho/2 + o(1)) q log q`, so
the coefficient is `7.96/rho`. `2 q^{1.99} <= m^{1.995}` needs `q^{0.005} >~ 2 (2/rho)^{1.995}`, that is
`q >~ 10^{184}`. This agrees with the judge's remark that the bound needs `q > 10^{180}`.

"Arc constant one" is correct: the model's arc constant is `pi/20 <= 1` (README_VIII §1.5). The statement about
(VIII.10) (`23/(3 rho)` at fixed density, `23/3` on the diagonal) is also consistent.

## 2. IX.6: the growth-scale fake in a prescribed rational-direction window

### 2.1 Corrected statement (Proved)

Fix `b = 2+i`, `beta = arg b`. Let `K >= 1` and `a > 0` be constants with
`dist(k beta, (pi/4) Z) >= a (2+|k|)^{-K}` for all integers `k != 0`. These exist by the Baker-Wustholz bound for
the two logarithms `log((2+i)/(2-i))` and `log(-1)`. **They depend only on `b`** (degree 2, height `(1/2) log 5`). They
are "absolute" only because `b` is fixed. Fix `A > 8(K+1)`. Then for all sufficiently large primes `q = 3 mod 4`,
with `m = q+1`, there is an angle model in `A_1(2m)` (round 5, Definition 1.2: all Lemma 0.1 monomial gaps for every
`v in Z^s \ {0}`, axiom (R) for all four units at all odd prime powers `<= 2m`, zero windings) with:
`s = 9m - 8` primes `= 1 mod 4` in a window `[P, P e^{1/(4s)}]`, `m^A <= P <= 2 m^A`;
`W <= s(log P + 1/(4s)) = (9A + o(1)) m log m`; arc constant 1; and all `m` angles in
`[beta + m^{5/4} Delta, beta + (m^{5/4} + 1) Delta]`, `Delta = e^{-W/4}`. Consequently the real numbers
`Re(conj(b) z_x)`, with `z_x = R e^{i theta_x}` the model points, lie in an interval of length
`sqrt5 (m^{5/4} + 1/2) + o(1) = O(m^{5/4})`.

The `O(m log m)` constant is `9A` and depends on `b` through `K`.

### 2.2 Step-by-step check

1. **Profile (Computed and Proved).** I rebuilt the fully flipped Paley-I profile `P(H_q, id, 8, a*)` of
   construct.md §2 (`fakecheck.py`, exact, `q = 7, 11, 19, 23`). It has `M = q + 1` rows and
   `r = 8q + (q+1) = 9M - 8` columns. Its rank is `M`, and for every row `i` some pair of columns differs by
   exactly `+-2 e_i` (the flipped copy of row `i` and an unflipped copy of `a(i)`). So `s = 9m - 8` and full row rank
   hold.
2. **Integrality of the extra monomials (Proved, with an addition).** If `v in Z^s` lies in the real row span,
   `v = S^T c'`, then each twin difference gives `+-2 c'_i in Z`. So `c = 2c'` is integral and unique by full row
   rank, as Astra says. *Addition:* every entry of `S` is `+-1`, so `(S^T c)_j = sum c = n (mod 2)` for all `j`
   (checked on random `c`). Hence `v = S^T c/2` is integral **iff** `c in Z^m` and `n` is even. The phase is then
   `<v, phi> = (n/2) theta_0 + c.epsilon/2`, with `n/2` an integer. This only simplifies the two-logarithm step.
   It also shows that the zero-sum case `n = 0` is exactly the round-5 character set `L`.
3. **Round-5 hypotheses (Proved).** `P_* = max(160 M^2 L/C, 21 r^2)` with `L ~ 2M log M` (for `Q_0 = 2M`) is
   `~ 320 m^3 log m < m^A` once `A > 3`. So `A > 4` suffices, as stated. Pigeonholing `[m^A, 2m^A]` into
   `ceil(4 s log 2)` windows gives the prime window, exactly as in the round-5 window lemma.
4. **Independence of `theta_0` (Proved).** For zero-sum `c`, `<v(c), phi> = (1/2) c.epsilon` (construct.md §1.1).
   `psi_u(c)` depends only on `u` and `k = 0`. The Lemma 4.1 residues `s omega_x / conj(omega_x)` do not involve the
   angles. So Steps 1-4 of the proof of Theorem 2 (probability `> 1/2` that (5.1) holds for every nonzero zero-sum
   `c`) concern `epsilon` only. Prescribing `theta_0` is legitimate. *Remark:* in round 5 the monomials of the row
   span with `n != 0` were handled implicitly by the Haar completion, in which `theta_0` is free. Astra correctly
   sees that prescribing `theta_0` makes them deterministic functions of `epsilon`, so they need separate protection.
   That is the content of the next two steps.
5. **Small extra monomials, `t = ||c||_1 <= e^{kappa W}`, `kappa = 1/(16(K+1))` (Proved).** `w(v) >= B log P`
   with `B = ||c||_inf`: for a coordinate with `|c_i| = B`, two entries of `v` differ by `B`, so their absolute
   values sum to at least `B`. The displacement of `<v, phi>` from `n beta/2` is at most
   `(|n| m^{5/4} + t) Delta/2 <= (1 + m^{5/4}) t Delta/2`. This is
   `<= m^{5/4} 3^{K+1} e^{(K+1) kappa W - W/4} (2+t)^{-K} = o((2+t)^{-K})`, uniformly, because `W >= s log P` grows
   faster than `log m`. With `t <= mB`, `2 + mB <= (3m)^B` and `P >= m^A`, `A > 2K` already gives
   `2 P^{-B/2} <= (a/4)(2+mB)^{-K}` for all `B >= 1` once `m` is large; `A > 8(K+1)` is more than enough. Both
   Lemma 0.1 inequalities follow from `dist(<v, phi>, (pi/4) Z) >= arcsin(e^{-w(v)/2})`. **Correct.**
6. **Tail `t > e^{kappa W}` (Proved).** Then `B >= e^{kappa W}/m`, and at most `(2B+1)^m` vectors have sup norm
   `B`. Condition on all offsets but one coordinate with `|c_i| = B`. The phase is then uniform on an interval of
   length `B Delta/2`. That interval meets at most `2B Delta/pi + 2` forbidden arcs, each of length
   `<= 2.1 P^{-B/2}`. So the bad probability is `O((1 + 1/(B Delta)) P^{-B/2}) = O((1 + Delta^{-1}) P^{-B/2})`.
   With `Delta^{-1} = e^{W/4}`, Astra's inequality `W/4 + m log(2B+1) <= B log P/4` holds on the whole tail: at
   `B_0 = e^{kappa W}/m`, `m log(2B_0 + 1)/B_0 -> 0` and `W/B_0 -> 0`, and `log(2B+1)/B` decreases. So the tail sum
   is `<= sum_{B >= B_0} O(P^{-B/4}) = o(1)`. **Correct.** The estimate is unconditional under the product law of
   `epsilon`, so it may be intersected with the round-5 event of probability `> 1/2`.
7. **Haar completion (Proved; Computed).** The realisations with fixed `theta_0, epsilon` form a coset of the image
   `T` of `ker S` in `(R/2 pi Z)^s`. `T` is a closed connected subtorus because `S` is integral. An integral `v` outside the
   real row span restricts to a nontrivial character of `T`, so its value is Haar-uniform on the circle. The bad
   measure for one `v` is `<= (8/pi) * 1.01 e^{-w(v)/2}`, and
   `sum_{v != 0} e^{-w(v)/2} = prod_p (1 + 2/(sqrt p - 1)) - 1`. Computed (`fakecheck.py`, the `s` smallest
   primes `= 1 mod 4` above `m^A`):

   | `A` | `m = 20` | `m = 100` | `m = 1000` |
   |---|---|---|---|
   | 4 | 1.35 | 0.195 | 0.0181 |
   | 5 | 0.212 | 0.0180 | 0.000569 |
   | 8 | 0.00215 | 1.78e-5 | — |

   The product minus one is `~ 2s/sqrt P ~ 18 m^{1 - A/2}`. This is `o(1)` for every `A > 2`, so certainly for
   `A > 8(K+1)`. **Correct.**
8. **Points and projection (Proved; Computed).** The angles are `theta_0 + epsilon_x`, in an arc of length
   `R Delta = e^{W/4} = sqrt R`, so the arc constant is 1. `Re(conj(b) z) = sqrt5 R cos(theta - beta)` with
   `theta - beta in [x_0, x_0 + Delta]`, `x_0 = m^{5/4} Delta`. The variation is
   `<= sqrt5 R (x_0 + Delta) Delta = sqrt5 (m^{5/4} + 1)`, because `R Delta^2 = 1`. Computed with mpmath at
   `W = 200`: length / `m^{5/4}` = 2.2989, 2.2396, 2.23608, 2.236068 at `m = 10, 10^2, 10^4, 10^6`, tending to
   `sqrt 5`.
9. **Two-logarithm input (Proved, with source caveat).** `e^{i beta} = (2+i)/sqrt5`, so `e^{2 i beta} = (2+i)/(2-i)`
   and `e^{4 i n beta} = ((2+i)/(2-i))^{2n}`. This base has modulus one and is not a root of unity. The bound
   `|alpha^k - 1| >= c k^{-K'}` is a two-logarithm Baker-Wustholz bound for `k log alpha - 2 pi i j`. It is linear in
   `log k`, which is what gives a polynomial `K`. An elementary bound gives only `5^{-k/2}`
   (`alpha^k - 1 = ((2+i)^k - (2-i)^k)/(2-i)^k`). That is too weak here, because `5^{-mB/2}` cannot be beaten by
   `P^{-B/2} = m^{-AB/2}`. So Baker's theory is genuinely needed. I did not reread the Baker-Wustholz paper. The
   record's refereed `round31/ref_baker.md` is the source, as for VIII.2. Numerically (`fakecheck.py`, `n <= 2*10^5`),
   `n * dist(n beta, (pi/4) Z)` stays between 0.09 and 0.53 at the record minima. This is consistent with an
   irrationality exponent near 2, far better than what the proof needs.

### 2.3 Points needing correction (wording only)

* "absolute `a > 0, K >= 1`": they depend on `b` (here fixed), hence so do `A`, `kappa` and the constant in
  `W = O(m log m)`.
* `n` is always even for integral `v` (§2.2 item 2). Harmless.
* "Taylor's theorem ... interval of length `O(m^{5/4})`": the constant is `sqrt5 + o(1)` times `m^{5/4}`.

### 2.4 What the fake shows, and what it does not (the fairness question)

In the fake the projections `Re(conj(b) z_x)` are **real numbers**, not integers. The larger sieve behind (IX.1)
counts integers in an interval by their residues modulo small primes. A model supplies residues only of the abstract
points (`rho_x^{(l)}`, `l <= 2m`), with no integer whose residues they are. So no contradiction with (IX.1) arises,
and none is claimed. Precisely:

* **Shown.** (IX.1) is not a consequence of the `A_1(2m)` axioms. These are the monomial gaps for all integral
  monomials, axiom (R) at all odd prime powers `<= 2m` with all four units, and the arc constraint. Read for models,
  (IX.1) says: "for every `0 != b in Z[i]`, the real projections `Re(conj(b) z_x)` have diameter `>> M^{3/2}`". This
  fails in `A_1(2m)` at growth scale (`M = m ~ s/9`, `W = O(m log m)`), already for the single direction `2+i` and
  with an explicitly moving common phase (no drifting phase is replaced by zero). Equivalently, the restricted rate
  theorem of IX.6 cannot be proved from `A_1(2m)`. Its hypothesis `|b|(1 + sqrt(R') |sin(theta - arg b)|) = o(s^{3/2})`
  holds in the fake (the left side is `~ sqrt5 m^{5/4}`), yet `M ~ s/9` is not `o(s)`. So any proof of IX.2/IX.6 must
  use integrality of the projected coordinates. As Astra says, it does.
* **Not shown.**
  1. Nothing about genuine clusters. The window around `2+i` is excluded for genuine clusters by IX.2 itself, not
     by the fake.
  2. Not that (IX.1) is independent of the stronger axioms of IX.3 (P3-char at every cutoff, prime-level factored
     residues). This fake has only round-5 point residues up to `2m`, and its primes are at height `m^A` with
     `A > 8(K+1)`, far outside the deep window. Whether the IX.3 or IX.5 fakes can be placed in a rational-direction
     window was not addressed.
  3. Nothing about directions of large height. A genuine cluster need not have any `b` with
     `|b|(1 + sqrt(R') |sin|) = o(s^{3/2})`. By Dirichlet one can find `|b| <= H` with `|sin(theta - arg b)| <~ 1/(H|b|)`.
     Optimising `H` gives a product of order `R'^{1/4}`, and for a typical (badly approximable) `theta` nothing much
     smaller is possible. `R'^{1/4}` is exponentially large in `s`, so the criterion bites only on exceptionally
     well-approximable orientations. Astra says this herself ("a restricted theorem on the scale of orientations, without a
     parent-only consequence").

Labels: the fake is **Proved**. "For any genuine integer marks the larger sieve would give `m << m^{5/6}`" is a
correct arithmetic consequence of (IX.1) (`m^{3/2} << m^{5/4}`), conditional on IX.2, which this block did not
referee.

## 3. "What the Gram matrix still permits"

**Setting.** reply9 (G3): `G_ij = W - 2 d_ij = 4 log l_ij - 4 log|b_ij|`, with `l_ij = |z_i' - z_j'|/N'^{1/4}`. For
a model, define `B_ij := l_ij exp((d_ij - W/2)/2)`. Then `G = W I + 4(L - log B)` off the diagonal holds by
definition. `G = sum_p log p sum_tau sigma sigma^T` is positive semidefinite, because it is built from the sign
matrices of the model's own exponent rows. The identity `log|b_ij| = log l_ij + (d_ij - W/2)/2` is the genuine one
(`log|g_ij| = (W - d_ij)/2`), so real magnitudes satisfying (G3) exist in every angle model.

**Check of the pair surplus in IX.3 (Proved).** In IX.3, `a = 1 + K` and `e = 2` on the `q - 1` labels, and the
marks lie in the class `K_{i0} = 1`, so column 0 contributes nothing. Then `d_ij - W/2 = sum_j c_j log p_j` with
`c_j = |K_ij - K_i'j| - 1 in {-1, 0, 1}` and `sum c_j = (q+1) - (q-1) = 2` (row distance `q+1`). Hence
`d_ij - W/2 = 2 log P_0 + sum_j c_j log(p_j/P_0)`, with `log(p_j/P_0) in [0, sigma]` and
`sigma = O_C(q log q/P) -> 0`. Under the random permutation the mean of the error is in `[0, 2 sigma]`. Its
deviation is `O(sigma sqrt(q log q)) = o(1)` simultaneously over the `O(q^2)` pairs, by bounded differences. So the
estimate is two-sided: `d_ij - W/2 = 2 log P_0 + o(1)`. Then
`log B_ij <= log C' + log P_0 + o(1) = O_C(log m)`, since `P_0 = q^{d + o(1)}`.

**IX.5 (not rechecked).** The `g`-column statement `2 g log P_0 + o(1)` needs all `g(q-1)` primes in one narrow
window `[P_0, P_0 e^sigma]` with `sigma sqrt(g q log q) -> 0`. IX.5's "details of the replicated profile" assert
owner norms between `P_0^g` and `(P_0 e^sigma)^g`, which suggests this. I did not verify that the conditional
angular-sum construction keeps such a window. That belongs to the IX.5 block. If the primes only lie in `[P, 2P]`,
the max claim could fail (the error could reach `O(sqrt(q log q))`), but the mean claim below still holds.

**What the conclusion amounts to (Proved; novelty minimal).** The mean statement needs none of this.
Positivity of `G` on the all-ones vector gives, in **every** angle model, `mean_{i != j} log B_ij <=
W/(4(M-1)) + log C'` (reply9 (G3)). So every class containing fakes with `W = O(M log M)` already permits
`mean log B_ij = O(log M)`. That includes `A_1(2m)` (round 5), the VIII.3 window fakes, and IX.3 with P3-char at all
cutoffs. So "the Gram identity with real magnitudes plus the residue axioms of IX.3 cannot force
`mean log|b_ij| / log m -> infinity`" is correct. It is a corollary of the existence of the IX.3 fakes at
`W <= (22/3 + o(1)) m log m`, together with the observation that the mean of `log|b|` and `W/(M log M)` are tied by
the all-ones vector. The two-sided estimate adds the stronger fact that the **maximum**, not just the mean, is
`O(log m)`. So no other eigenvector or minor of `G` can force a single large real cofactor either.

**What remains outside the statement.** The following genuine constraints are not used by the real-magnitude model
and are not imitated by these fakes, as far as is shown:
1. integrality: `|b_ij|^2` is a rational integer, and `b_ij in Z[i]` with `(1+i) | b_ij`;
2. the **cocycle** `g_ij b_ij + g_jk b_jk + g_ki b_ki = 0` with `g_ij` the actual Gaussian gcds determined by the rows
   and the prime factors `pi_p`. This is a complex, not merely real, relation;
3. the coupling of the `b_ij` across all pairs through one set of Gaussian primes. In the fake, `arg b_ij` is not
   tied to the angles `phi_p`. The genuine relation `z_i' = eps_i prod pi_p^{a_ip} conj(pi_p)^{e_p - a_ip}` fixes
   both `g_ij` and `b_ij`.

Astra names (2), as "(G4)" in her text. The cocycle is listed in reply9's (G4) paragraph, and the label is a
cross-reference to that paragraph, not to a separate axiom. Items (1) and (3) are implicit in "Gaussian".

## 4. Files

* `ref9_restate/two_adic.py`: 2-adic unit classes and fibre sizes `f_k`.
* `ref9_restate/mod2block.py`, `mod2block2.py`; outputs `mod2block_1e6.out`, `mod2block2_1e6.out`: exact window
  tests of the 2-adic blocks, `N <= 10^6`, arc `<= 4`.
* `ref9_restate/fakecheck.py`: Paley profile (rank, twins, parity), Haar union product, projection length,
  irrationality data for `arctan(1/2)`.
