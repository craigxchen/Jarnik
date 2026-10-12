# Round 8: which Diophantine hypothesis gives a growth improvement? Rank reduction, simultaneous Lang–Waldschmidt, and barriers for LW, abc, Subspace and Bombieri–Lang

Round 8. The target was the weakest standard Diophantine hypothesis that implies a **growth
improvement** `M = o(log R/loglog R)` for lattice points on arcs of length `C sqrt R` of `x^2 + y^2 = R^2`
(a uniform bound is a bonus). For each candidate, I either prove that it implies growth or prove that it
does not suffice. Context: `README.md` files of `growth/`, `threshold/`, `round4/`, `round5/`; the
barrier summary of the task. Scripts are in `round8/weakest_checks/` (§9).

Conventions. *Theorem*, *Proposition* and *Lemma* mean proved here with all quantifiers, or cited
exactly (refereed notes of this review, or classical results). *Hypothesis* means an assumption,
stated precisely. *Exact check* means integer or `Fraction` arithmetic, with angles at 60 digits.
*Evidence* means a computation that no proof uses. The literature could not be consulted in this session
(network access to arXiv and authors' pages was blocked), so statements about published results are
from memory and are marked where it matters.

## 0. Verdict

**No standard named conjecture weaker than Vojta's was found to suffice. The obstruction was located
exactly: it is the rank of the exponent lattice. A simultaneous Lang–Waldschmidt bound for five linear
forms in logarithms, with a power to spare relative to the Dirichlet exponent, implies the uniform
theorem. A version for logarithmically many forms implies growth with a power saving. The standard
single-form Lang–Waldschmidt conjecture, abc, Bombieri–Lang and the quantitative Subspace Theorem
(with any known `S`-dependence) do not suffice.** For LW the blindness is proved unconditionally, and
the existence of LW-consistent fake circles at `W = O(M log M)` uses the draft Theorem 5.1 of `exact.md`.
The simultaneous hypothesis is introduced here; I do not claim it is a named conjecture.

| # | result | type |
| --- | --- | --- |
| Thm 1.3 | **Rank reduction.** A normalised cluster whose exponent lattice has rank `k` has `M <= (5.5^k + 1)/2` points (for `W >= W_0(C)`). So the uniform theorem is equivalent to bounded rank (Cor. 1.4), and growth follows from rank `<= theta log W` with `theta < 1/log 5.5` | unconditional, new |
| Thm 2.3 | **`SLW(k, g, kappa)` with `g < 1/4` implies the uniform theorem**, with `M <= (5.5^(k-1)+1)/2` and rank `<= k - 1`. For `k = 5`: `M <= 458` | conditional, new |
| Thm 2.4 | **`SLW^log(theta, g, kappa)` (at least `max(k_0, theta log(n+2))` forms in `n` logarithms), with `theta < 1/log 5.5` and `g < 1/4`, implies `M = O((log R)^(1-eta))`**. This is a growth improvement with a power saving. The hypothesis is implied by every `SLW(k_0, g)` | conditional, new |
| Prop 2.5 | **`SLW(4, g)` is false for `g < 1/4`** (golden eight-point family, rank 4, attains the Dirichlet exponent `1/4`). **`SLW(3, g)` is false for `g < 2/7`** (Cilleruelo–Granville family at Pell parameters, `(a+1)^2+1 = 5 (x^2+y^2)^2`). **`SLW(1, g)` is false for `g < 1/2`**. So five forms is the first case in which the hypothesis can hold | unconditional, new |
| Prop 2.7 | `SLW(4, g)` in the narrow window `1/4 <= g < 4/15` also implies the uniform theorem (via 5-point agreement) | conditional, new |
| Lemma 3.2 | **Support lemma for the skew double core.** Every nonzero rational character is nonzero on at least `(q+1)/2 = M/2` of the `2q` cuts (prime-field uncertainty for the pair `K, K^T`). The bound is sharp up to `+3` | unconditional, new |
| Thm 3.3, Lemma 3.4, Cor 3.5 | **Single-form Lang–Waldschmidt (any `eps`, any constant) is blind on the skew double core.** Every cluster-determined form has LW bound `<= e^(-eps q tau/4) e^(-W/4)`. Random-position fake-circle constructions satisfy all these LW inequalities with probability `1 - e^(-Omega(W))`. With the draft fake circles of `exact.md` Thm 5.1 (`W = O(M log M)`), **LW plus local data cannot prove growth** | unconditional (conditional on the draft for the final clause), new |
| Props 4.1–4.2, Cor 4.3 | **abc (over `Q` or `Q(i)`) is vacuous on squarefree profiles**, and all fake circles are squarefree, so abc does not suffice. abc gives only the restricted theorem for `rad(N) <= N^(1/4-eta)` and, for two-level profiles, `rad(N) >= N^(1/2-o(1))` | conditional / barrier |
| Thm 5.2, Prop 5.3 | **Subspace Theorem with `S = {inf} u primes(N)`.** If the exceptional subspaces of the cluster system `Sigma_m` number `f(s) = o(s)`, growth follows (bounded `f`: uniform). The known bounds are `exp(O(s))`, which is useless. For `m <= 3` every admissible `f` is `>= s^(2-delta-o(1))` (explicit pigeonhole solutions), so the needed hypothesis is false there. For `m >= 4`, `delta > 1` it is open and non-standard | conditional + barrier, new |
| Prop 6.1 | **Vojta:** one instance `V(m_0, eps)` for one fixed `m_0 >= 7` suffices for the uniform theorem. Dimension dependence is irrelevant, larger `m_0` is weaker, and no Vojta-type statement gives growth without the uniform bound | unconditional reduction |
| Props 6.2–6.4 | **Bombieri–Lang** is vacuous (all auxiliary varieties, including the chord-data variety, are rational). The bounded reduced chords `beta_xy` are an `S`-integrality condition (Cayley coordinates), with `S` containing every inert prime `< M - 1`. Lang–Vojta would need `S`-uniformity, which is non-standard | unconditional, barrier |

**Comparison with the known conditional theorem** (`growth/conditional.md`: `V(7, 1/4)` implies the uniform
theorem). `SLW(5, g)` and `V(7, 1/4)` are not known to be comparable. Both are non-sharp Dirichlet-type
statements:
* Vojta's canonical threshold `2/5` against the needed `1/2` (7 points);
* SLW's Dirichlet exponent `1/5` against the needed `1/4` (5 forms, 6 points).

SLW reaches the problem one point earlier because it sees the rank. Vojta's conjecture is exactly
borderline at 6 points, where the Cilleruelo–Granville and axis-Pell families live; they have rank
`<= 3`. SLW's conclusion is effective in `M` (`458`). The Dirichlet-sharp SLW for all `k`, like Vojta,
implies the full Cilleruelo–Granville `R^(1-eps)` conjecture (Remark 2.4b). `SLW` is not a named
conjecture that I could verify in the literature. Its `k = 1` case is Lang–Waldschmidt type at exponent
`g` instead of `1 + eps`.

**Evidence (§7).**
* The golden family attains `gamma_4 -> 1/4` (rank 4).
* The Cilleruelo–Granville–Pell family attains `gamma_3 -> 2/7`.
* Among all circles with `N <= 10^8`, rank-5 windows have SLW exponents with median `0.105` and maximum
  `0.234`.
* Clusters of rank 5 with `C <= 6` exist at small `N` (17 circles up to `10^8`) and grow more slowly than
  lower ranks; the data cannot test the asymptotic statement.
* Every infinite bounded-`C` family whose rank I determined has rank `<= 4`. A polynomial template of rank 5
  would need degree `>= 8` (§7.3), and none is known.

## 1. Setting, characters, and the rank reduction

### 1.1 Conventions

These follow `two_thirds/proof.md` §0–1 and `round4/outside.md` §1.

* `N >= 1`, `R = sqrt(N)`, `W = log N`. A *configuration* is a set of `M >= 2` distinct Gaussian
  integers `z_x` with `|z_x|^2 = N` on one arc of length `<= C sqrt(R)`. Its angular width is
  `Delta <= C R^(-1/2) = C e^(-W/4)`. Real lifts `theta_x` of `arg z_x` are chosen in one interval of
  length `Delta`.
* *Normalisation* (`two_thirds` Lemma 1.1). Dividing by the Gaussian gcd keeps `M`, keeps the arc
  constant `<= C` and does not increase `N`. A normalised configuration has `N` odd, every prime
  factor of `N` split, and
  `z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`, with `0 <= a_xj <= e_j`, `N = prod_j p_j^(e_j)`
  (`j = 1..r`), a fixed Gaussian prime `pi_j` above `p_j`, and a unit `eps_x`. If `g` is the gcd and
  `R' = R/|g|`, then `R' >= R^(1/2)/C` (`growth/conditional.md` Lemma 4.1(a): `|g| <= |z_x - z_y| <= C sqrt R`).
  So `W' := log N' >= log R - 2 log C` tends to infinity with `R`.
* `phi_j = arg pi_j`, `w(v) = sum_j |v_j| log p_j` for `v in Z^r`, `alpha_j = pi_j / conj(pi_j)`.
* The *exponent lattice* of the configuration is `Lambda = span_Z {a_x - a_y}` in `Z^r`. Its rank `k` is
  the **rank of the configuration**. It does not depend on the choice of the `pi_j`.

**Lemma 1.1 (gap; `outside.md` Lemma 0.1, `walsh-extraction.md` Lemma 1.1).** For every
`v in Z^r \ {0}`, `dist(<v, phi>, (pi/4) Z) >= arcsin(e^(-w(v)/2) / sqrt 2)`.

**Lemma 1.2 (difference forms).** For all `x, y`, `<a_x - a_y, phi> = (theta_x - theta_y)/2 (mod (pi/4) Z)`.
Moreover `Log(z_x/z_y) = b_0 Log i + sum_j (a_xj - a_yj) Log alpha_j` for an integer `b_0` with
`|b_0| <= 5 + 2 sum_j e_j`, where `Log` is the principal logarithm.

*Proof.* `arg z_x = arg eps_x + <2 a_x - e, phi> (mod 2 pi)`, and `arg eps_x in (pi/2) Z`. Subtract and
halve. For the second claim, `z_x/z_y = (eps_x/eps_y) prod_j alpha_j^(a_xj - a_yj)`. So `Log(z_x/z_y)`
differs from `Log(eps_x/eps_y) + sum_j (a_xj - a_yj) Log alpha_j` by `2 pi i m` with
`2 pi |m| <= pi/2 + pi + pi sum_j e_j` (all principal logarithms have modulus `<= pi`, and
`|Log(z_x/z_y)| <= Delta <= pi/2`). Put `b_0 = Log(eps_x/eps_y)/(i pi/2) + 4m`. QED

### 1.2 The rank reduction

**Theorem 1.3 (points are exponential in the rank).** Let a normalised configuration of rank `k`
lie on an arc of length `<= C sqrt R`, and assume `W >= W_0(C) := 9 (4 log^+ C + 2 log 2)`. Then the `M`
exponent vectors `a_x` are distinct, and

```text
M  <=  ( (1 + 2W/rho)^k + 1 ) / 2  <=  ( 5.5^k + 1 ) / 2,        rho := W/2 - 2 log C - log 2.
```

*Proof.*
1. *Distinct vectors.* If `a_x = a_y` with `x != y`, then `z_x/z_y` is a unit other than 1, so
   `|theta_x - theta_y| >= pi/2`. But `Delta <= C e^(-W/4) < pi/2` once `W >= W_0(C)`.
2. *Separation.* Let `D = {a_x - a_y}` be the difference set. Take `u = a_x - a_y` and
   `u' = a_x' - a_y'` in `D` with `u != u'`, and put `v = u - u' != 0`. By Lemma 1.2,
   `<v, phi> = (theta_x - theta_y - theta_x' + theta_y')/2 (mod (pi/4) Z)`, and the right side has
   modulus `<= Delta`. Lemma 1.1 gives `e^(-w(v)/2)/sqrt 2 <= Delta <= C e^(-W/4)`, that is
   `w(v) >= rho`. Taking `u' = 0` shows `w(u) >= rho` for `u != 0` as well.
3. *Packing.* `w` is a norm on `R^r` (a weighted `l^1` norm), hence on `V = Lambda (x) R`, which has
   dimension `k`. Since `|a_xj - a_yj| <= e_j`, every `u in D` has `w(u) <= W`. The `w`-balls of radius
   `rho/2` around the points of `D` are pairwise disjoint (step 2) and lie in the ball of radius
   `W + rho/2`. Comparing volumes in `V`: `|D| (rho/2)^k <= (W + rho/2)^k`.
4. *Counting.* Order `R^r` lexicographically and label the (distinct) exponent vectors
   `a_1 < ... < a_M`. Then `0` and the `2(M-1)` vectors `+-(a_x - a_1)` (`x >= 2`) are distinct
   elements of `D`. So `|D| >= 2M - 1`.
5. *Constants.* Since `log C <= log^+ C`, `rho >= W/2 - 2 log^+ C - log 2`. If `W >= W_0(C)` then
   `2 log^+ C + log 2 <= W/18`, so `rho >= 4W/9` and `1 + 2W/rho <= 11/2`. QED

**Corollary 1.4.**
(a) *Uniform theorem iff bounded rank.* Fix `C`. Arcs of length `C sqrt R` carry boundedly many
lattice points for `R >= R_0` if and only if normalised configurations on such arcs have bounded rank
for `R >= R_0`. (`M >= k + 1` gives one direction; Theorem 1.3 and `W' -> infinity` give the other.)
(b) *Growth criterion.* If, for some `theta < 1/log 5.5` and all large `W'`, every normalised
configuration on an arc of length `C sqrt R` has rank `<= max(4, theta log(2 + W'/log 5))`, then
`M <= max(5.5^4, (2 + W'/log 5)^(theta log 5.5)) = O((log R)^(1 - eta))` with `eta = 1 - theta log 5.5 > 0`.
This is far below `log R/loglog R`.

**Remark 1.5 (where the difficulty sits).** For squarefree `N` the trivial bound is `M <= 2^k`.
Theorem 1.3 extends an exponential bound to all exponent profiles, using only the gap Lemma 1.1 on
four-term characters. Every infinite family with bounded `C` whose rank I determined has rank at most 4
(§2.4, §7):
* the golden eight-point family has rank 4;
* the Cilleruelo–Granville six points have rank 3;
* symmetric clusters `{z, conj z}`, such as the axis Pell family, have rank `<= M/2 = 3`, because
  conjugation sends `2a - e` to `-(2a - e)`;
* five-point families (the cyclotomic family of the notes) have rank `<= 4` trivially.

I did not compute the rank of the translated-Pell eight-point family of the notes.

The near-extremal profiles of the fake circles (Hadamard and skew double cores) have the full rank
`M - 1`. So any growth or uniform theorem has to exclude high-rank configurations, and Corollary 1.4 says
this is all it has to do. The rank is invisible to every single linear form (§3). It is visible to
simultaneous forms (§2).

## 2. Simultaneous Lang–Waldschmidt: five forms suffice; four only in a narrow window

### 2.1 The hypothesis

For `alpha in Q(i)^*`, `H(alpha)` is the naive height: the largest absolute coefficient of the primitive
minimal polynomial of `alpha` over `Z`. Put `A(alpha) = max(H(alpha), e)`. For `alpha_j = pi_j/conj(pi_j)`
the minimal polynomial is `p_j X^2 - 2 Re(pi_j^2) X + p_j`, which is primitive, so
`p_j <= H(alpha_j) <= 2 p_j` (`round4/flipped.md` §0).

**Hypothesis `SLW(k, g, kappa)`** (simultaneous Lang–Waldschmidt, `k` forms, height exponent `g`).
There is `c in (0, 1]` with the following property. Let `n >= 1`, and let `alpha_1, ..., alpha_n in Q(i)`
have modulus 1 and be multiplicatively independent modulo roots of unity. Let `(b_ij)`
(`1 <= i <= k`, `0 <= j <= n`) be an integer matrix whose rows `(b_i1, ..., b_in)` are linearly
independent over `Q`. Put `Lambda_i = b_i0 Log i + sum_(j=1..n) b_ij Log alpha_j` (principal logarithms).
Then

```text
max_i |Lambda_i|  >=  c^(n+1)  prod_(j=0..n) (1 + B_j)^(-kappa)  prod_(j in U) A(alpha_j)^(-g),       (SLW)
B_j = max_i |b_ij|,   U = { 1 <= j <= n : B_j > 0 }.
```

(Each `Lambda_i` is nonzero: a vanishing `Lambda_i` would be a multiplicative relation between the
`alpha_j` modulo roots of unity.)

**Hypothesis `SLW^log(theta, g, kappa)`.** There are an integer `k_0 >= 1` and one constant `c` such that
(SLW) holds for every system whose number of forms satisfies `k >= max(k_0, theta log(n + 2))`.

**Remarks 2.1.**
1. *Status.* For `k = 1` this is a Lang–Waldschmidt-type lower bound (Lang, *Elliptic Curves:
   Diophantine Analysis*, Grundlehren 231, 1978, pp. 212–217; Waldschmidt, *Diophantine Approximation
   on Linear Algebraic Groups*, Springer 2000, Ch. 1). There the exponent is `1 + eps` on the naive
   heights, and the factor in front is `C(eps)^n B` or `C(eps) B`, where `B` is the largest coefficient
   (the normalisations differ between sources; `flipped.md` §4 explains why the naive height is the
   right one here). Simultaneous linear forms in logarithms are studied by Philippon and Waldschmidt.
   Their proved bounds are of Baker type (products of heights), and these are useless here for the
   reason given in `flipped.md` Prop. B. **I could not consult these texts in this session.** Network
   access to arXiv and the authors' pages was blocked. The precise form `SLW(k, g, kappa)` is therefore
   introduced here as a hypothesis. I do not claim it is a named conjecture in the literature. What
   supports it is items 3–4 below and §7.
2. *Monotonicity* (Lemma 2.2). `SLW(k, g, kappa)` implies `SLW(k', g', kappa')` whenever `k' >= k`,
   `g' >= g` and `kappa' >= kappa`. In particular `SLW(k_0, g, kappa)` implies `SLW^log(theta, g, kappa)` (with
   the same `k_0`) for every `theta`. Larger `k` gives a weaker hypothesis.
3. *Dirichlet heuristic.* Let the `alpha_j` range over unit-circle elements with `A(alpha_j) <= A_j`.
   There are about `prod A_j` choices. Treat the `k` values `Lambda_i` as independent and uniform near
   0. Then the chance that all `k` are `<= delta` is about `delta^k`. The expected number of systems with
   `max |Lambda_i| <= (prod A_j)^(-g)` is then about `(prod A_j)^(1 - gk)`. So (SLW) should hold for every
   `g > 1/k` and fail below. Every theorem below uses `g < 1/4` with `k >= 5`, which leaves room
   `1/4 - 1/k >= 1/20` in the exponent. The golden family (Prop. 2.5(a)) shows the heuristic is attained, but
   not beaten, at `k = 4`.
4. *Sharpness and falsity for small `k`.* Proposition 2.5 shows that `SLW(4, g)` is false for every
   `g < 1/4` (golden family, which attains the Dirichlet exponent `1/4` with `k = 4`). It also shows that
   `SLW(3, g)` is false for `g < 2/7`, and `SLW(1, g)` is false for `g < 1/2`. So `k = 5` is the first
   value at which the hypotheses used below can hold.

**Lemma 2.2 (monotonicity).** If (SLW) holds for `(k, g, kappa)` with constant `c`, it holds for
`(k', g', kappa')` with the same `c` whenever `k' >= k`, `g' >= g` and `kappa' >= kappa`.

*Proof.* Take `k'` forms with independent rows. Choose `k` of them with independent rows. Their
maximum is at most the maximum over all `k'` forms. Their union support `U_k` is contained in `U_(k')`,
and their `B_j` are `<=` those of the `k'` forms. Apply (SLW) to the chosen forms in the
`n_k = |U_k| <= n` numbers they involve. Then
`c^(n_k+1) >= c^(n+1)`, `prod_(U_k) A^(-g) >= prod_(U_(k')) A^(-g')` (as `A >= e`), and the factors
`(1 + B_j)^(-kappa)` only grow. QED

### 2.2 Five simultaneous forms imply the uniform theorem

**Theorem 2.3.** Assume `SLW(k, g, kappa)` for some `k >= 1`, `g < 1/4` and `kappa >= 0`. Then for every
`C > 0` there is `R_0(C)` such that for all `R >= R_0(C)` every normalised configuration on an arc of
length `<= C sqrt R` of `x^2 + y^2 = R^2` (in any position) has rank `<= k - 1`. Consequently every such
arc contains at most `(5.5^(k-1) + 1)/2` lattice points.

By Proposition 2.5 the hypothesis is false for `k <= 4`, so the theorem has content only for `k >= 5`.
`SLW(5, g, kappa)` with any `g < 1/4` gives at most `floor((5.5^4 + 1)/2) = 458` points.

*Proof.* Normalise. By §1.1 it suffices to treat normalised configurations with `W = log N` as large as
we like. Suppose the rank is `>= k`. Choose points `x_0, ..., x_k` with `a_(x_i) - a_(x_0)` (`i = 1..k`)
linearly independent; this is possible because `Lambda` is spanned by the `a_x - a_(x_0)`. Put
`Lambda_i = Log(z_(x_i)/z_(x_0))`. By Lemma 1.2 these are `k` forms as in (SLW). Their rows are the
`a_(x_i) - a_(x_0)`, so they are independent. The numbers are `alpha_j = pi_j/conj(pi_j)` for `j` in the
union support `U`, and these are multiplicatively independent modulo roots of unity by unique
factorisation in `Z[i]`. Also `|Lambda_i| = |theta_(x_i) - theta_(x_0)| <= Delta <= C e^(-W/4)`.

Bound the right side of (SLW) with `n = |U| <= r`:
* `A(alpha_j) <= 2 p_j`, so `sum_(j in U) log A(alpha_j) <= sum_j log p_j + r log 2 <= W + r log 2`;
* `B_j <= e_j` for `j >= 1`, and `B_0 <= 5 + 2 sum_j e_j <= 5 + 2W/log 5` (Lemma 1.2).

Fix `P >= 5` and split the primes at `P`. There are at most `P` split primes `<= P`, and those above `P`
number at most `W/log P`, so `r <= P + W/log P`. For `p_j > P`,
`log(1 + e_j) <= e_j log 2 <= e_j log p_j (log 2/log P)`. For `p_j <= P`, `log(1 + e_j) <= log(1 + W)`.
Hence

```text
sum_(j>=1) log(1 + B_j)  <=  W log 2 / log P  +  P log(1 + W).
```

Taking logarithms in (SLW) and using `|Lambda_i| <= C e^(-W/4)`:

```text
W/4 - log C  <=  (P + W/log P + 1)(log(1/c) + g log 2)  +  kappa ( W log 2/log P + P log(1+W) + log(6 + 2W) )  +  g W.
```

Choose `P = P(c, g, kappa)` with `(log(1/c) + log 2 + kappa log 2)/log P <= (1/4 - g)/2`. Then

```text
(1/4 - g) W / 2  <=  log^+ C + O_(c, g, kappa)( log(2 + W) ),
```

which fails for `W >= W_1(C, c, g, kappa)`. So for `W >= max(W_0(C), W_1)` the rank is `<= k - 1`, and
Theorem 1.3 bounds `M`. Since `W' >= log R - 2 log C`, both thresholds hold for `R >= R_0(C)`. QED

**Theorem 2.4 (growth from logarithmically many forms).** Assume `SLW^log(theta, g, kappa)` with
`theta < 1/log 5.5` and `g < 1/4`. Then for every `C` and `R >= R_0(C)`, every arc of length
`<= C sqrt R` contains at most

```text
max( 5.5^4 ,  (2 + 2 log R / log 5)^(theta log 5.5) )  =  O( (log R)^(1 - eta) ),      eta = 1 - theta log 5.5 > 0,
```

lattice points, with `5.5^4` replaced by `5.5^(k_0 - 1)` if the hypothesis holds with `k_0 > 5`. This is a
growth improvement with a power saving.

*Proof.* As above, with `k_1 = max(k_0, ceil(theta log(r + 2)))`. If the rank is `>= k_1`, choose `k_1 + 1`
points of rank `k_1`. They involve `n <= r` numbers, so `k_1 >= max(k_0, theta log(n+2))` and
`SLW^log` applies. The displayed inequality in the proof of Theorem 2.3 does not involve `k`, so the same
`W_1` works. Hence the rank is `<= k_1 - 1`. Then Theorem 1.3, `r <= W/log 5` and `W <= 2 log R` give the
bound. QED

`SLW^log(theta, g, kappa)` is implied by `SLW(k_0, g, kappa)` for every `k_0` (Lemma 2.2), so it is implied by the
hypothesis of Theorem 2.3. It is weaker in form: for systems in `n` numbers it asks for a lower bound only
when there are at least `theta log(n+2)` forms. The Dirichlet exponent for that many forms is
`1/(theta log n) -> 0`, far below the `1/4` used.

### 2.3 What the hypothesis means for circles

**Remark 2.4a (lattice-point reformulation).** Let `beta_i = gamma_i/conj(gamma_i)` (`i = 1..k`) with
`gamma_i` conjugate-primitive Gaussian integers. Put `gamma_0 = 1`, and define
`z_i := gamma_i prod_(l != i) conj(gamma_l)` for `i = 0..k`. All `z_i` have norm `prod_l N(gamma_l)`, and
`z_i/z_0 = beta_i`. So every system `Lambda_i = Log beta_i` is the system of a `(k+1)`-point
configuration on one circle, of rank `k` if the `beta_i` are independent. For exponent vectors with
entries in `{0, +-1}` and squarefree normalised norm (so `prod_U A = N'^(1+o(1))`), (SLW) with `g < 1/4`
therefore says: *`k+1` lattice points of rank `k` do not lie on an arc of length `R^(1 - 2g - o(1))`
(for large `R`).* For
`k = 5` and `g` in `(1/5, 1/4)`, that arc is longer than `C sqrt R` by the power `R^(1/2 - 2g)`. This is
the precise sense in which `SLW(5, g)` is an "a priori" form of the uniform theorem: it predicts that
the rank of clusters, not their size, is bounded, and it predicts this on longer arcs.

**Remark 2.4b (Dirichlet-sharp SLW implies the `R^(1-eps)` conjecture).** Suppose `SLW(k, (1+eps')/k)`
holds for every `k`. Run the proof of Theorem 2.3 on an arc of length `C R^alpha` (`alpha < 1`). The
angular width is `C e^(-(1 - alpha) W/2)`. Theorem 1.3 holds with `rho = (1 - alpha) W - O(log C)` and
the factor `1 + 2/(1 - alpha) + o(1)` in place of `5.5`. The contradiction appears as soon as
`(1 + eps')/k < (1 - alpha)/2`. So the rank is `<= 2(1+eps')/(1-alpha)`, and
`M <= (1 + 3/(1-alpha))^(2(1+eps')/(1-alpha))`. This is the Cilleruelo–Granville conjecture for
arcs `R^(1-eps)`.

**Remark 2.4c (comparison with Vojta's threshold).** In the units of (SLW), a configuration on an
arc `C sqrt R` has `log(1/|Lambda|) / log prod A ~ (1/2 log R)/(2 log R) = 1/4`. Vojta's conjecture on
`X_m^#` (`growth/conditional.md` §5) makes `m`-point tuples degenerate once `lambda_Y/h > 2/(m-2)`, which
is exponent `1/(m-2)` in these units. `m` points give `k = m - 1` forms, so Vojta needs
`1/(k-1) < 1/4`, that is `k >= 6` (`m = 7`), while `SLW(k, g)` needs `k >= 5` (`m = 6`). At `m = 6` Vojta's
conjecture is exactly borderline and says nothing. Indeed the Cilleruelo–Granville and axis-Pell
six-point families exist. SLW decides the six-point case because it sees the rank: those six-point
families have rank `<= 3`. The two hypotheses are not known to be comparable.

### 2.4 Fewer than five forms

**Proposition 2.5.** For every `c in (0,1]` and `kappa >= 0`:
(a) `SLW(4, g, kappa)` is false for every `g < 1/4`; hence (Lemma 2.2) `SLW(k, g, kappa)` is false for
all `k <= 4` and `g < 1/4`;
(b) `SLW(3, g, kappa)` is false for every `g < 2/7`;
(c) `SLW(1, g, kappa)` is false for every `g < 1/2`.

*Proof.* In each case an infinite family of configurations is exhibited, of rank `>= k`, with
`max |Lambda_i| <= C' N^(-1/4)` (resp. `<= 2 N^(-1/2)` in (c)). Each family has
`log prod_(j in U) A(alpha_j) <= log rad(N) + omega(N) log 2`, `n = omega(N) = o(log N)` and
`sum_(j>=1) log(1 + B_j) <= log d(N) = o(log N)`, `log(1 + B_0) = O(log log N)`. So (SLW) would force
`log(1/max |Lambda_i|) <= (g + o(1)) log rad(N)`. The facts used are
`omega(N) = O(log N/loglog N)` and `log d(N) = O(log N/loglog N)` (Hardy–Wright, Thm 317 and §22.10).

(a) *Golden family* (`growth/families.md` Prop. G, refereed). For `n = 1 mod 3` the eight odd words
`w in {0,1}^4` give points `Z_w(n)` of one norm. They are divisible exactly by `(1+i)^2`, the
normalised norm is `N' = phi^(8n+16)/20 (1 + o(1))`, and the angular spread is
`2 |q|^n (2 sqrt5 - 3)/phi (1+o(1)) = O(N'^(-1/4))`.
*Rank 4.* The norm of `J_r = f_(r+1) + i f_r` is `f_(2r+1)`, and `gcd(f_(r+1), f_r) = 1`. So no rational
prime divides `J_r`, and for an odd prime `p | f_(2r+1)` exactly one of `pi_p, conj(pi_p)` divides `J_r`.
For the four odd indices `2n+1, ..., 2n+7`, `gcd(f_a, f_b) = f_(gcd(a,b))` with `gcd(a,b) in {1,3}`
(it divides `2, 4` or `6` and is odd), and `f_1 = 1`, `f_3 = 2`. So the odd parts of the four norms
are pairwise coprime. For `m >= 13` odd, `f_m` has a primitive prime divisor (Carmichael 1913). This
prime is odd and is not 5, because 5 is primitive for `m = 5`. Hence each `J_(n+s)` has a private odd prime
`p_s != 5`. At `p_s` the exponent of `Z_w` is `e w_s` or `e (1 - w_s)`. The differences of the
odd words span `Q^4` (for example `1000, 0100, 0010, 0001, 1110` have independent differences). So
these five points have rank 4 in prime coordinates.
*Conclusion.* `log rad N' <= log N'`, so (SLW) would give `(1/4) log N' - O(1) <= (g + o(1)) log N'`,
which fails for `g < 1/4` and large `n`. The family attains the Dirichlet exponent `1/4` for `k = 4`
exactly (§7, Table 7.1).

(b) *Cilleruelo–Granville family at Pell parameters.* Let `X^2 - 5Y^2 = 1` with `X < 0`, put
`x = X - 2Y`, `y = Y` and `a + 1 = 2(x^2 - y^2) - 2xy`. Then `a + 1 + i = (2 + i)(x + iy)^2`, because the
imaginary part of the right side is `(x^2 - y^2) + 4xy = (x + 2y)^2 - 5y^2 = 1`. Hence
`(a+1)^2 + 1 = 5 (x^2 + y^2)^2`. The orbit `(X, Y) -> (9X - 20Y, -4X + 9Y)` from `(-9, 4)` gives
`a + 1 = 682, 219602, 70711162, ...` (infinitely many, increasing).
* *Points.* `z_sigma = prod_(j=1..4) (a + j + i sigma_j)`, `sigma in {+-1}^4`, `sum sigma = 0`: six points
  of norm `N = prod_j F_j`, `F_j = (a+j)^2 + 1`, with normalised constant `C -> 4 sqrt 2`
  (`families.md`; exact values in Table 7.1).
* *Gcd.* A Gaussian prime `w != 1+i` dividing all six points divides some `a + j + i` and some
  `a + j' - i`. Two factors with the same sign of `i` are coprime except at `1+i` (their difference is a
  rational integer of modulus `<= 3`, and 3 is inert, so it cannot divide `a + j + i`). Hence `w` divides at
  most one `+` factor and at most one `-` factor, and the pattern taking the other signs there avoids `w`
  (patterns have two `+` and two `-`). So the gcd is a power of `1+i`. Since `F_j = 2 (mod 4)` for odd
  `a + j`, the gcd is `(1+i)^2` and `N' = N/4`.
* *Rank 3.* Primes shared by two `F_j` divide the resultants `d^2(d^2+4)` (`d = 1, 2, 3`), so they lie in
  `{2, 5, 13}`. By Størmer's theorem, `x^2 + 1 = 2^a 5^b 13^c` has only finitely many solutions. So for
  all but finitely many members every `F_j` has a private prime `> 13`, and the exponent map is
  injective on sign patterns. Four of the six patterns with affinely independent sign vectors (for
  example `++--, +-+-, +--+, -++-`) then have rank 3.
* *Radical.* `rad(F_1) <= 5(x^2+y^2) <= sqrt(5 F_1)`, so `rad(N') <= sqrt5 (a+5)^7` while
  `N' >= a^8/4`. Hence `log rad N' <= (7/8) log N' + O(1)`, and (SLW) would give
  `(1/4) log N' - O(1) <= (7g/8 + o(1)) log N'`. This fails for `g < 2/7`.

(c) *Liouville-extremal pairs* (`round4/hardness.md` Prop. 3.2(3)): `z = u`, `w = i conj(u)` with
`u = A + (A-1)i`, `n = A^2 + (A-1)^2`, `|Log(z/w)| = 2 arcsin((2n)^(-1/2)) <= 2 n^(-1/2)`, rank 1, and
`log prod_U A <= log n + omega(n) log 2`. QED

**Remark 2.6.** In (b) the exponent `2/7` comes from one square factor. Any rank-4 family with a free
parameter, in which congruence conditions could be imposed, would refute `SLW(4, g)` above `1/4` in
the same way. The golden family has no free parameter. For `k = 5` no rank-5 family is known at all
(§7). So `SLW(5, g)` for `g` in `(1/5, 1/4)` is consistent with everything known, and the Dirichlet
heuristic predicts it with room to spare.

**Proposition 2.7 (four forms in the window `1/4 <= g < 4/15`).** Assume `SLW(4, g, kappa)` with `g < 4/15`.
Then the uniform theorem holds: there is an absolute `M_0` such that arcs of length `C sqrt R` carry at
most `M_0` lattice points for `R >= R_0(C)`. By Proposition 2.5(a) such a hypothesis can hold only for
`g >= 1/4`. It is consistent with the golden family if and only if that family has
`log rad(N') >= (1/(4g) - o(1)) log N'`. This is expected (Fibonacci numbers are believed to be
squarefree up to bounded factors) but not proved.

*Proof.*
1. *Sub-configurations of rank `<= 3` are small.* In the proof of Theorem 1.3, steps 2–4 apply verbatim to
   any subset `A'` of the exponent vectors of a normalised configuration, with the `W` of the whole
   configuration. So a subset spanning an affine space of dimension `<= 3` has at most
   `L := (5.5^3 + 1)/2 < 84` points, provided `W >= W_0(C)`.
2. *Few rank-deficient 5-subsets.* A 5-subset of affine rank `d <= 3` is determined by `d + 1` affinely
   independent points together with `4 - d` further points of their affine span, and that span contains
   `< 84` points. So there are at most `sum_(d <= 3) M^(d+1) 84^(4-d) <= 4 * 84^4 M^4` such subsets, a
   fraction `f_M <= 4 * 84^4 * 120 * 2^5 / M` of all `binom(M, 5)` (for `M >= 8`).
3. *Averaging the support.* Let `X` be a uniformly random 5-subset and `U(X)` the set of primes at which
   the points of `X` do not all have the same exponent.
   * For `e_j = 1`, convexity gives `P(j notin U(X)) = (binom(m,5) + binom(M-m,5))/binom(M,5)`, which is
     at least `(1/16)(1 - 10/(M-4))`.
   * For `e_j >= 2`, `log p_j <= e_j log p_j / 2`.

   Hence `E sum_(j in U(X)) log p_j <= (15/16 + 1/(M-4)) W`. Some `X` of affine rank 4 therefore has
   `sum_(j in U(X)) log p_j <= (15/16 + 1/(M-4)) W/(1 - f_M)`.
4. *Contradiction.* Apply `SLW(4, g)` to the four forms `Log(z_(x_i)/z_(x_0))` of that `X`, and bound
   exactly as in Theorem 2.3, with `sum_(U) log A <= sum_(U) log p_j + r log 2`. For `M >= M_1` (so that
   `g (15/16 + 1/(M-4))/(1 - f_M) < 1/4 - (1/4 - 15g/16)/2`) and `W >= W_1(C, c, g, kappa)` this is
   impossible.
5. *Conclusion.* So for `W' >= max(W_0, W_1)` either the rank is `<= 3` (and `M < 84`) or `M < M_1`. QED

The window `[1/4, 4/15)` is narrow, and a rank-4 family with a free parameter would close it (Remark 2.6).
`SLW(5, g)` with `g < 1/4` needs no averaging, and its exponent can be taken anywhere in
`(1/5, 1/4)`. That makes it the more robust hypothesis.

## 3. Single-form Lang–Waldschmidt does not suffice: the skew double core is LW-blind

### 3.1 The hypothesis

**Hypothesis `LW(eps, C_*)`** (Lang–Waldschmidt over `Q(i)`, naive heights, as in `flipped.md` (H_kappa)
but with `n` numbers). Let `n >= 1`, `alpha_1, ..., alpha_n in Q(i)^*` and `alpha_0 = i`. Let
`b_0 in Z` and `b_1, ..., b_n in Z \ {0}` be such that `Lambda = sum_(j=0..n) b_j Log alpha_j != 0`. Then

```text
|Lambda|  >=  C_*^(n+1) B prod_(j=0..n) ( max(1, |b_j|) A_j )^(-(1+eps)),     B = max_j |b_j|,   A_j = max(H(alpha_j), e^|Log alpha_j|, e).
```

Every standard normalisation with constant `C(eps)^n` or `C(eps)` in front implies `LW(eps, C_*)` with
`C_* = max(1, C(eps))`. So consistency with `LW(eps, C_*)` for all `C_*` covers all of them. (The `n` numbers
are those with a nonzero coefficient; with `C_* > 1` dummy numbers must not be allowed.)

Single forms of pair type are already Liouville-optimal (`hardness.md` Thm 3). LW gains only through
large coefficients on low-height numbers. That is the mechanism of `flipped.md` Theorem L, where a label
character of a single Hadamard core has support weight `W/M + F`.

### 3.2 The skew double core

Let `q = 3 (mod 4)` be prime and `M = q + 1`, with rows indexed by `{inf} u F_q`. Let `chi` be the Legendre symbol mod `q`.
Define `S(inf, a) = -1`, `S(a, inf) = 1`, `S(a, b) = chi(b - a)` (`a, b in F_q`), `S(inf, inf) = 0`, and
`K = I + S`. The profile (`round5/exact.md` §4, a draft; the facts used are re-proved or checked here)
has `2q` cuts, each carrying one split prime `p_T` with exponent 1:

```text
S_a = { x : K(x, a) = -1 },      S'_a = S_a xor {inf, a}      (a in F_q).
```

Point `x` has exponent `[x in T]` at `p_T`. For zero-sum `lambda in Q^M` the character is `v_T = lambda(T) = sum_(x in T) lambda_x`.

**Lemma 3.1 (structure).** `K` is a skew-Hadamard matrix (`K K^T = M I`, `K + K^T = 2I`). Every cut has `M/2`
elements, and every pair of rows is separated by exactly `M` of the `2q` cuts. The rational characters
(zero-sum `lambda` with `v_lambda in Z^(2q)`) are `Z^M_0 + Z lambda*`, `lambda* = (2/M) 1 - 2 e_inf`. In
particular every nonzero rational character has `||lambda||_inf >= 2/M`.

*Proof.* `S` is skew because `chi(-1) = -1`. The Jacobsthal block `Q(a,b) = chi(b-a)` satisfies
`Q Q^T = qI - J` and `Q 1 = 0`. A block computation then gives `S^2 = -qI`, so
`K K^T = (I+S)(I-S) = I - S^2 = MI`. Column `a in F_q` of `K` is orthogonal to column `inf` (all ones),
so it is balanced, and `|S_a| = M/2`. `S'_a` has the same size because `inf in S_a` and `a notin S_a`.
`S_a` is the cut of columns `a, inf` of `K`, and `S'_a` is the cut of columns `a, inf` of `K^T`
(`exact.md` §4 checks this row by row). Two orthogonal rows of a Hadamard matrix with a constant
column disagree, relative to that column, in exactly `M/2` columns. This happens once for `K` and once for `K^T`.
*Characters.* `lambda(S'_a) - lambda(S_a) = lambda_a - lambda_inf in Z` gives `lambda = lambda_inf 1 + c`
with `c in Z^M`. Then `lambda(S_a) = lambda_inf M/2 + c(S_a) in Z` gives `lambda_inf in (2/M) Z`; the
converse is clear. Exact check: `check_skew_support.py`, all assertions for `q = 7, 11, 19, 23, 31, 43`. QED

With all `log p_T` in `[tau, (1 + eta) tau]`, `W <= 2q(1 + eta) tau`, and every pair character has weight
`>= M tau`. With equal weights, `W = 2q tau` and the pair margin is `w(a_x - a_y) - W/2 = tau`.

**Lemma 3.2 (support lemma; new).** Every nonzero zero-sum `lambda in Q^M` has
`#{T : lambda(T) != 0} >= (q+1)/2 = M/2`.

*Proof.* Scaling does not change which `lambda(T)` vanish, so we may assume `lambda in Z^M` is primitive.
For `a in F_q`, `lambda(S_a) = -(K^T lambda)_a / 2` (since `sum lambda = 0`) and
`lambda(S'_a) = lambda(S_a) - lambda_inf + lambda_a`. Let `Z = {a in F_q : lambda(S_a) = lambda(S'_a) = 0}`.
For `a in Z`: `lambda_a = lambda_inf`, and
`0 = (K^T lambda)_a = -lambda_inf + lambda_a + sum_(x in F_q, x != a) chi(a - x) lambda_x = sum_(x in F_q) chi(a - x) lambda_x`.
Let `d = (q-1)/2`. Since `chi(y) = y^d (mod q)` for all `y in F_q`, the polynomial
`P(X) = sum_(x in F_q) lambda_x (X - x)^d in F_q[X]` vanishes on `Z`. Its coefficient of `X^(d-m)` is
`(-1)^m binom(d, m) sum_x lambda_x x^m`, and `binom(d, m) != 0 (mod q)` for `0 <= m <= d`.

*Case `P != 0`.* Then `|Z| <= d`, so at least `q - d = (q+1)/2` values `a` have `lambda(S_a) != 0` or
`lambda(S'_a) != 0`. Hence at least `(q+1)/2` cuts carry a nonzero value.

*Case `P = 0`.* Then `sum_x lambda_x x^m = 0 (mod q)` for `m = 0..d`. Put `mu_x = lambda_x - lambda_inf`
(`x in F_q`). Since `sum_(x in F_q) x^m = 0 (mod q)` for `0 <= m <= q - 2` (with `0^0 = 1`) and
`d <= q - 2`, also `sum_x mu_x x^m = 0 (mod q)` for `m = 0..d`. Now `mu` vanishes on `Z`. If
`|Z| >= d + 1`, then `|supp mu| <= q - d - 1 = d`, and the invertible Vandermonde system on `supp mu`
(using `m = 0..|supp mu| - 1`) gives `mu = 0 (mod q)`. So `lambda_x = lambda_inf (mod q)` for all `x`.
Then `0 = sum lambda = (q+1) lambda_inf = lambda_inf (mod q)`, so every entry is divisible by `q`,
contradicting primitivity. Hence `|Z| <= d` again. QED

*Check* (`check_skew_support.py`). The bound holds on 452,480 vectors for `q = 7, ..., 43`
(exhaustive over small supports, structured, random). The minimum found is `(q+7)/2`, attained by
`lambda = K(e_a - e_b)`, so the lemma is sharp up to `+3`. As `q -> infinity` the support fraction
`(q+7)/(4q)` tends to `1/4`: the skew double core sits exactly on the LW threshold `W/4`, from above.

### 3.3 LW is blind on this profile

**Theorem 3.3.** Let `eps > 0`, `C_* > 0`, `C_1 = max(1, C_*)`. Suppose that
`tau <= log p_T <= (1 + eps/2) tau` for all `T`. For a nonzero rational character `lambda` with
`v = v_lambda`, put

```text
rho_lambda := C_1^(2q+1) e^(-(1+eps)(q+1) tau/2) ||v||_inf  prod_(T in supp v) |v_T|^(-(1+eps)).
```

Then for every `b_0 in Z`, the right side of `LW(eps, C_*)` for the form
`Lambda = b_0 Log i + sum_T v_T Log alpha_T` is at most `rho_lambda`, and

```text
rho_lambda  <=  C_1^(2q+1) e^(-eps q tau/4) e^(-W/4).
```

So on this profile, for every cluster-determined form, LW only excludes values within
`C_1^(2q+1) e^(-eps q tau/4)` times the arc's angular scale `e^(-W/4)`. That factor is `e^(-Omega(W))` once
`tau >= tau_0(eps, C_*)`. The upper bound that the arc condition gives for these forms is
`||lambda||_1 Delta` (`outside.md` (1.2)). So no single form yields a contradiction.

*Proof.* The `n` numbers with nonzero coefficient are the `alpha_T`, `T in supp v`, so `n <= 2q`.
`A_0 >= 1` and `A_T >= H(alpha_T) >= p_T >= e^tau`. If `B = |b_0|`, then
`B max(1,|b_0|)^(-(1+eps)) <= 1`. Otherwise `B = ||v||_inf` and `B max(1,|b_0|)^(-(1+eps)) <= ||v||_inf`.
In both cases the bound is at most
`C_1^(n+1) ||v||_inf prod_(supp v) |v_T|^(-(1+eps)) e^(-(1+eps) tau |supp v|)`, and
`|supp v| >= (q+1)/2` by Lemma 3.2. Since `prod_(supp) |v_T| >= ||v||_inf`, we get
`rho_lambda <= C_1^(2q+1) e^(-(1+eps)(q+1)tau/2)`. Finally `W <= 2q(1 + eps/2) tau` and
`(1+eps)(q+1)/2 - q(1+eps/2)/2 >= eps q/4`. QED

**Lemma 3.4 (random positions).** Keep the assumptions of Theorem 3.3. Let `Delta in (0, 1/2]`,
`k in Z^M`, and let `delta` be uniform on `[0, Delta]^M`. For a rational character `lambda` put
`chi_lambda(delta) = (1/2) sum_x lambda_x (delta_x + (pi/2) k_x)`. In an angle model this is the exact value of
`<v_lambda, phi>` (`outside.md` (1.2)). Then

```text
P( exists lambda != 0 : dist(chi_lambda(delta), (pi/4) Z) < rho_lambda / 2 )
   <=  (6M / Delta) (2 sqrt q)^(1+eps) (1 + 2 zeta(1+eps))^q  C_1^(2q+1)  e^(-(1+eps)(q+1) tau/2).
```

If moreover `Delta >= C e^(-W/4)` (arc constant `C`) and `tau >= tau_0(eps, C_*) + (8/(eps q)) log^+(1/C)`,
where `tau_0(eps, C_*) := (8/eps)(log((1 + 2 zeta(1+eps)) C_1^3) + 2)`, this is at most `e^(-eps q tau/8)`.

*Proof.*
1. *One character.* Pick `x` with `|lambda_x| = ||lambda||_inf >= 2/M` (Lemma 3.1). Conditionally on the
   other coordinates, `chi_lambda` is uniform on an interval of length `l = ||lambda||_inf Delta/2`. The set of
   points within `rho/2` of `(pi/4) Z` meets such an interval in measure `<= (4l/pi + 2) rho`. So the
   probability is `<= (2 + 4/(||lambda||_inf Delta)) rho_lambda <= (3M/Delta) rho_lambda` (using `M >= 8`,
   `Delta <= 1/2`).
2. *Injectivity.* By Lemma 3.1 the block values `(lambda(S_a))_a` determine `lambda`, since
   `(K^T lambda)_inf = sum lambda = 0` and `K` is invertible. The same holds for `(lambda(S'_a))_a`, since
   `lambda(S'_a) = (K lambda~)_a / 2` with `lambda~ = lambda - 2 lambda_inf e_inf`, `(K lambda~)_inf = 0`
   (`exact.md` §4.1(iii)).
3. *A second large coordinate.* `sum_a lambda(S_a)^2 = ||K^T lambda||^2/4 = M ||lambda||_2^2/4`, and likewise
   for the `S'_a`. Let `T*` carry `||v||_inf <= ||lambda||_1 <= sqrt M ||lambda||_2`. The other block then has
   squared norm `>= ||v||_inf^2/4`, hence a coordinate `T**` with `|v_(T**)| >= ||v||_inf/(2 sqrt q)`.
4. *Summation.* `||v||_inf prod_(supp) |v_T|^(-(1+eps)) <= ||v||_inf^(-eps) |v_(T**)|^(-(1+eps)) prod_(T in B(T*), T != T*) |v_T|^(-(1+eps))`,
   where `B(T*)` is the block of `T*` and the other factors were bounded by 1. This is
   `<= (2 sqrt q)^(1+eps) prod_(T in B(T*), v_T != 0) |v_T|^(-(1+eps))`. Sum over `lambda` by summing over the
   values on the block of `T*` (injective, step 2) and over the two choices of block:
   `sum_lambda (...) <= 2 (2 sqrt q)^(1+eps) (1 + 2 zeta(1+eps))^q`.
5. *Constants.* Combine steps 1 and 4 with `1/Delta <= C^(-1) e^(W/4)`, `W/4 <= q(1 + eps/2) tau/2` and
   `(1+eps)(q+1)/2 - q(1+eps/2)/2 >= eps q/4`. QED

**Corollary 3.5 (LW-consistent fake circles).** Consider any random construction of angle models on
the skew double core in which:
* the positions `delta` are uniform on `[0, Delta]^M`;
* the primes satisfy `tau <= log p_T <= (1 + eps/2) tau` with `tau >= tau_0`;
* every other ingredient (windings, residue data, the free part of `phi`) is independent of `delta`;
* the construction succeeds with probability `>= p_0`.

Then with probability `>= p_0 - e^(-eps q tau/8)` its output also satisfies `LW(eps, C_*)` for **every** form
whose value is determined by the cluster: all rational characters and all `b_0`.

The union-bound constructions of `construct.md` (refereed) and `exact.md` §5 (draft) draw `delta` uniformly
in exactly this way. `exact.md` Theorem 5.1 claims such angle models in the class `A_1(2M)` on this
profile with `W <= (12 + o(1)) M log M`. Its primes lie in `[P, P e^(1/(4r))]` with
`tau = log P = 4 log M + 2 log L + O(log(1/C)) -> infinity`, so the spread condition holds for every fixed
`eps`. It succeeds with probability `>= 1/2` (its union bound loses `<= 0.414 + 10^(-100)`; the free part
of `phi` is chosen off a Haar-null set, independently of `delta`). **If that draft theorem holds, then no argument that combines single-form
Lang–Waldschmidt bounds for cluster-determined forms with the local inputs of `A_1(2M)` can prove
`M = o(log R/loglog R)`.** The only unproved input is the draft theorem itself. Its key combinatorial
claims (Lemma 3.1, `sigma_c >= ||c||_1/2`) are re-checked exactly in `check_skew_support.py` and
`exact.md`'s own checkers.

**Remarks 3.6.**
1. *Forms not determined by the cluster.* A form whose coefficient vector is not a rational character
   takes a value that depends on the free part of `phi` (`outside.md` (1.2) and Theorem A.1). For an
   actual circle, LW applied to it says something about the prime arguments alone and does not involve
   the cluster. Making the free part of a fake circle satisfy all these inequalities as well is a
   genericity statement of the type of `outside.md` Theorem A.1. I do not prove it here. The barrier is
   therefore stated for arguments that apply LW to cluster-determined forms.
2. *Single core versus double core.* A single Hadamard core (with copies and flips) always has a
   character supported on one label plus the flipped primes. LW sees that character (Theorem L). The
   skew double core uses the two orthogonal bases `K` and `K^T`. Lemma 3.2 is a prime-field uncertainty
   principle for this pair: no rational combination of rows is sparse in both. This confirms, for an
   explicit and provable profile, the lead's suspicion about profiles built from two or more
   independent bases (there stated for random permutations).
3. *The threshold is exactly LW's.* The support fraction tends to `1/4`, and LW's exponent `1 + eps` on
   the support weight barely misses `W/4`. A single-form bound with exponent `< 1` on the support would
   see the profile. Such a bound is false in general by Liouville-extremal pairs (Prop. 2.5(c)). It is
   simultaneity that sees it: with `lambda_1 = K(e_a - e_b)` and `lambda_2 = K(e_a - e_c)`, the union support
   is `3q/4 + O(1)` columns, that is `(3/8 + o(1)) W` (`check_skew_support.py`, second part: minimum over
   `b, c` equal to 36, 54, 81, 117 for `q = 43, 67, 103, 151`). So `SLW(2, g)` would see the profile for
   every `g < 2/3`.

## 4. abc over Q(i): vacuous on squarefree profiles

**Hypothesis abc.** *Over `Q`* (Masser–Oesterlé): for every `delta > 0` there is `K_delta` such that coprime
positive integers with `a + b = c` satisfy `c <= K_delta rad(abc)^(1+delta)`. *Over `Q(i)`* (the
number-field form; Elkies, Vojta): coprime `a, b, c in Z[i]` with `a + b + c = 0` satisfy
`max(|a|,|b|,|c|)^2 <= K_delta N(rad(abc))^(1+delta)`, where `N(rad(abc))` is the product of the norms of
the distinct Gaussian primes dividing `abc`.

**Proposition 4.1 (what abc says about clusters).** Assume abc over `Q`, and let `delta > 0`.
(a) *Pairs* (`growth/conditional.md` Prop. 8.1, refereed). In a normalised configuration on an arc of
length `<= C sqrt R`, every pair satisfies
`sum_(j : a_xj != a_yj) log p_j = log rad(n_xy) >= (1/4 - delta) W - log C - O_delta(1)`, where
`n_xy = prod_j p_j^(|a_xj - a_yj|)` is the pair norm.
(b) *Restricted theorem* (`conditional.md` Cor. 8.2). If `rad(N) <= N^(1/4 - eta)`, every arc of length
`C sqrt R` has at most one lattice point once `R >= R_0(C, eta)`.
(c) *Two-level profiles.* Suppose that at every prime the points occupy only the levels `0` and `e_j`
(as in the squarefree case). Then at most `M^2/4` pairs are separated at each prime. Summing (a) over
pairs gives `log rad(N) >= (2(M-1)/M)((1/4 - delta) W - log C - O_delta(1))`. So such profiles carrying many
points have `rad(N) >= N^(1/2 - o(1))`.

**Proposition 4.2 (vacuity on squarefree profiles).** Let the normalised `N` be squarefree and
`W >= 8(2 log^+ C + log 2)`.
(a) The pair consequence 4.1(a) follows from the gap Lemma 1.1, because `rad(n_xy) = n_xy` and
`log n_xy = w(a_x - a_y) >= W/2 - 2 log C - log 2 >= W/4`.
(b) Let `G = gcd(z_x, z_y, z_w)`. Divide the triangle relation `(z_x - z_y) + (z_y - z_w) + (z_w - z_x) = 0`
by `G`. It becomes `A_xy c_xy + A_yw c_yw + A_wx c_wx = 0`, where `A_xy` is the product of the `pi_j^(+-)`
over the primes at which `x, y` agree and `w` differs, and `c_xy = (z_x - z_y)/g_xy`. The `A`'s are
squarefree and pairwise coprime, and `N(A_xy A_yw A_wx) = N/N(G)`. After removing any common factor `h`
of the three terms, the abc inequality holds for every `delta >= 0` and every `K_delta >= 1` as soon as
`C^2 <= N^(1/2)`:

```text
max |term|^2 / N(h)  <=  C^2 N^(1/2) / (N(G) N(h))  <=  N / (N(G) N(rad h))  <=  N(rad(product of the three terms))
```

(`N(h) >= N(rad h)`, and the Gaussian primes of the `A`'s that do not divide `h` still divide the
quotient). The same argument covers Ptolemy relations and the `n`-term abc conjecture of
Browkin–Brzeziński (exponent `2n - 5 + delta`); see `conditional.md` Prop. 8.3 for the general
statement (pairwise coprime squarefree blocks with bounded cofactors).

**Corollary 4.3 (abc does not suffice).** The fake circles of `construct.md` (§5: `e = 1`), `lower.md`
(§3: `e = 1`) and `exact.md` (§5: `e = 1`) all have squarefree `N`, and each satisfies the Liouville pair
inequalities. By Proposition 4.2, every abc consequence obtained through
pair, triangle, Ptolemy or `n`-term relations holds automatically for any squarefree configuration that
satisfies those pair inequalities. So no argument using abc through these relations, together with the
local inputs of the fake-circle classes, can prove `M = o(log R/loglog R)`. abc gives the restricted
theorem 4.1(b) for very powerful `N`, and the structural statement 4.1(c). Unconditionally, by
contrast, two-level prime powers behave like squarefree columns of weight `e_j log p_j` in the cut identity,
the chain lemma of `two_thirds` (which charges nothing for them) and the pair gaps. So these local inputs
do not force `rad(N)` to be close to `N`. The "close to squarefree" heuristic is therefore an abc phenomenon, and it
is beside the point: the obstruction is already squarefree.

*Precise scope.* abc is a statement about additive relations among integers. A proof that "abc
cannot help" must fix which relations are used. Proposition 4.2 covers every relation whose terms are
products of chords `z_x - z_y` with pairwise coprime squarefree prime blocks. It does not cover
relations with high-degree cancellation, which would have to come from algebraic identities among the
points. Those are the auxiliary-polynomial methods, barred by `threshold/upper.md`.

## 5. The quantitative Subspace Theorem with `S = {inf} u primes(N)`

### 5.1 The system and the reduction

Work over `K = Q(i)` with normalised absolute values: `|x|_inf = |x|^2` at the complex place and
`|x|_w = N(w)^(-v_w(x))` at a Gaussian prime `w`. Then `prod_v |x|_v = 1`, and for a primitive
`x in Z[i]^m` the height is `H(x) = prod_v max_i |x_i|_v = max_i |x_i|^2`. For a finite set `P` of split
rational primes, let `S_P` consist of `inf` and the places above `P`, with `s = |P|`.

**The system `Sigma_m(P)`.** At `inf` use the forms `X_1` and `X_i - X_1` (`i = 2..m`). At each finite
`w in S_P` use `X_1, ..., X_m`. All determinants are 1. For `delta > 0` consider

```text
prod_(v in S_P) prod_(i=1..m) |L_(i,v)(x)|_v / ||x||_v  <=  H(x)^(-m-delta),      x in K^m \ {0}.                    (5.1)
```

**Lemma 5.1.** Let `z_1, ..., z_m` be distinct lattice points of an arc of length `<= C sqrt R` on
`|z|^2 = N`, let `P` contain every prime factor of `N`, and let `x = (z_1, ..., z_m)/g` with
`g = gcd(z_i)`. Then `x` is primitive, `H(x) = N/|g|^2 >= R/C^2`, and the left side of (5.1) is
`<= C^(2(m-1)) H(x)^(-m-(m-1)/2)`. So `x` solves (5.1) for every `delta < (m-1)/2` once
`H(x) >= C^(4(m-1)/((m-1) - 2 delta))`.

*Proof.* The `x_i = z_i/g` have norm `N_x = N/|g|^2`, and all their prime factors lie above `P`. So at the
finite places of `S_P`, `||x||_w = 1` and `prod_w |x_i|_w = N(x_i)^(-1) = N_x^(-1)`. The finite part is
therefore `N_x^(-m)`. At `inf` the product is
`prod_(i>=2) |z_i - z_1|^2 / |z_1|^2 <= (C^2 R / N)^(m-1) = (C^2 N^(-1/2))^(m-1) <= C^(2(m-1)) N_x^(-(m-1)/2)`.
Finally `|g| <= |z_1 - z_2| <= C sqrt R` gives `N_x >= R/C^2`. QED

The excess `(m-1)/2` over the Subspace exponent `m` is large. But `S_P` depends on `N`, and the classical
theorem only gives finitely many subspaces for each fixed `S`.

**Hypothesis `QS(m, delta, f, H_0)`.** For every finite set `P` of split primes, with `s = |P|`, the
solutions of (5.1) with `H(x) >= H_0` lie in the union of at most `f(s)` proper linear subspaces of `K^m`.

**Theorem 5.2.** Assume `QS(m, delta, f, H_0)` for some `m >= 2` and `0 < delta < (m-1)/2`. Then every
arc of length `<= C sqrt R`, with `R/C^2 >= max(H_0, C^(4(m-1)/((m-1)-2delta)))`, carries
`M <= f(omega(N')) + m - 1` lattice points, where `N'` is the normalised norm. Consequently:
* if `f` is bounded, the uniform theorem holds;
* if `f(s) = o(s)`, then `M = o(log R/loglog R)` (a growth improvement);
* if `f(s) <= s^(1-eta)`, then `M <= (log R)^(1-eta+o(1))`.

*Proof.* All ordered `m`-tuples of distinct points of the arc solve (5.1) with `P = primes(N)` (Lemma 5.1).
There are `M(M-1)...(M-m+1)` of them, and by `QS` they lie in `<= f(s)` proper subspaces. A proper subspace
lies in a hyperplane `sum c_i X_i = 0` with some `c_(i0) != 0`. In that hyperplane `z_(i0)` is determined by
the other coordinates (the relation is homogeneous, so the scaling by `g` does not matter). So the hyperplane
contains at most `M(M-1)...(M-m+2)` of the tuples. Hence `M - m + 1 <= f(s)`. For the growth statement,
`s = omega(N')` and `W' >= sum_(k <= s) log(4k+1) >= s log(4s/e)` (the `k`-th split prime is `>= 4k+1`), so
`s = O(W'/log W')`. QED

### 5.2 How much `S`-dependence is unavoidable

The known bounds are exponential in `s`. Evertse (*Compositio Math.* 101 (1996) 225–311) bounds the
number of subspaces containing the solutions of height `>= max(m^(4m/delta), H(L))` by
`(2^(60 m^2) delta^(-7m))^s`. Evertse–Schlickewei (*J. reine angew. Math.* 548 (2002)) and Evertse–Ferretti
(*Ann. of Math.* 177 (2013)) improve the dependence on `m` and `delta` in the absolute and parametric
formulations. These statements are as I recall them; the texts were not accessible in this session. In
any case `QS` with `f(s) = exp(O(s))` gives only `M <= exp(O(s))`, which is weaker than the trivial
`M <= 4 d(N)`. The next proposition shows that for small `m` any admissible `f` grows faster than `s`.

**Proposition 5.3.** Let `P_s` be the first `s` split primes.
(a) For `m = 2` and `0 < delta < 1`, the system `Sigma_2(P_s)` has at least `c_delta s^(2-delta) (log s)^(-delta)`
solutions of height `>= H_0` that are pairwise distinct points of `P^1(K)`, for `s >= s_0(delta, H_0)`. So
`QS(2, delta, f, H_0)` forces `f(s) >= c_delta s^(2-delta-o(1))`.
(b) For `m = 3` and `0 < delta < 1`, the system `Sigma_3(P_s)` has at least `c_delta s^(4-2delta-o(1))`
solutions in pairwise distinct points of `P^2(K)`. Moreover any proper subspace of `P^2(K)` (a line)
contains at most `O(s^(2-delta))` of them. So `QS(3, delta, f, H_0)` forces `f(s) >= c s^(2-delta-o(1))`.
For `m = 2, 3` the range `delta < (m-1)/2` needed in Theorem 5.2 lies inside `delta < 1`. So **`QS(m, delta, o(s))`
is false for `m <= 3`.**

*Proof.* (a) Let `phi_j = arg pi_j`, reduce mod `pi/4`, and put `L = ceil((pi/2) P^delta)` with `P` the
largest prime of `P_s`. Partition `[0, pi/4)` into `L` intervals. By convexity at least
`s(s - L)/(2L)` unordered pairs `{i, j}` fall into a common interval. For such a pair put
`u = pi_i conj(pi_j)`. Then `dist(arg u, (pi/4) Z) <= pi/(4L) <= P^(-delta)/2`. Choose the unit `eta` with
`|eta conj(u)/u - 1| <= 2 dist(arg u, (pi/4) Z)` and put `x = (u, eta conj u)`. This is primitive with
`H(x) = n := p_i p_j <= P^2`. The left side of (5.1) is `<= 4 dist^2 n^(-2) <= P^(-2 delta) n^(-2) <= n^(-2-delta)`.
The ratio `x_2/x_1` determines `{i, j}`, so distinct pairs give distinct points. Pairs with
`p_i p_j < H_0` number `O_(H_0)(1)`. Since `P = O(s log s)` (prime number theorem in progressions) and
`L <= s/2` for `s >= s_0(delta)`, the count is `>= s^2/(4L) - O(1)`.
(b) Take two pairs from (a) with disjoint prime sets, `u = pi_i conj(pi_j)` and `u' = pi_k conj(pi_l)`, and
put `x = (u u', (eta conj u) u', u (eta' conj u'))`. All coordinates have norm `n n'`. The point `x` is primitive
(a common prime would divide `u'` and `conj(u')`, or `u` and `conj(u)`), and `H(x) = n n' <= P^4`. Its archimedean
factor is `|eta conj(u)/u - 1|^2 |eta' conj(u')/u' - 1|^2 <= 16 dist^2 dist'^2 <= P^(-4 delta) <= (n n')^(-delta)`, and
its finite factor is `(n n')^(-3)`. So (5.1) holds.
*Counting points.* In `P^2(K)` the point is `(1 : rho : rho')` with `rho = eta conj(u)/u` and `rho' = eta' conj(u')/u'`,
and `rho` (resp. `rho'`) determines the first (resp. second) pair. So distinct ordered choices give
distinct points. Pairs of pairs sharing a prime number at most `4 B p_2`, where `p_2` is the number of
pairs and `B <= s` is the largest bin. So there are `>= p_2 (p_2 - 4s) >= p_2^2/2` points once `p_2 >= 8s`,
which holds for large `s` by (a), since `delta < 1`.
*Counting per line.* Take a line `c_1 X_1 + c_2 X_2 + c_3 X_3 = 0`. If `c_3 != 0`, then `rho'` is determined by
`rho`. If `c_3 = 0`, then `rho` is fixed. Either way the line contains at most `p_2` of the points. So at
least `p_2/2 >= c s^(2-delta-o(1))` lines are needed. QED

*Check* (`check_subspace_count.py`, `delta = 0.45`). Every pair solution is verified at 60 digits. There
are 10, 33, 109, 333, 903 and 2635 of them for `s = 25, ..., 800` (against the predicted `s^2/(2L)`), and
3.45 million triples at `s = 800`. At `delta = 0.9` the pigeonhole gives none for `s <= 800`: here the
Liouville floor `dist >= (2n)^(-1/2)` is still of the order of the bin width. This is a finite-size effect;
asymptotically (a) needs only `delta < 1`.

**Remark 5.4 (what would suffice).** Product ("cube corner") configurations of `m - 1` pigeonhole pairs
solve (5.1) for every `delta < 1`. By Liouville on each pair (`|rho_i - 1| >= sqrt2 n_i^(-1/2)`), they do not
solve it for `delta > 1 + o(1)`. So for `m >= 4` and `delta in (1, (m-1)/2)` the constructions above give
nothing, and `QS(m, delta, f)` with `f(s) = o(s)` is neither refuted nor supported by them. By Theorem 5.2 it
would give growth; with `f` bounded it is an `S`-independent Subspace Theorem for large `delta`, which
would give the uniform theorem. Neither is a standard conjecture. The Subspace route therefore does not
give a "weakest standard hypothesis". It does give a precise target: sublinear `S`-dependence of the number
of exceptional subspaces of `Sigma_m` for some `m >= 4` and `delta in (1, (m-1)/2)`. For `m <= 3` the
corresponding bound is false (Prop. 5.3).

## 6. Vojta with growing dimension; Bombieri–Lang for chord data

### 6.1 Dimension dependence in Vojta's conjecture is irrelevant

Notation of `growth/conditional.md`: `X_m = {|z_1|^2 = ... = |z_m|^2}` in `P^(2m-1)`, `Y` the diagonal line,
`T^(m)_C` the short-arc `m`-tuples. `V(m, eps)` is Vojta's main conjecture (`D = 0`, `A = H`) for a smooth
model of `Bl_Y X_m`. It implies that `T^(m)_C` is not Zariski dense as soon as `eps < (m-6)/4`
(`conditional.md` Thm 5.2 at `alpha = 1/2`).

**Proposition 6.1.**
(a) If `T^(m')_C` is not Zariski dense in `X_(m')`, then `T^(m)_C` is not Zariski dense in `X_m` for every
`m >= m'`.
(b) Hence a single instance `V(m_0, eps)` with `m_0 >= 7` and `eps < (m_0 - 6)/4` gives the uniform theorem
(`conditional.md` Prop. 3.2(i)). No uniformity in `m`, and no explicit dependence on the dimension, is
needed for growth or for the uniform bound.
(c) Conversely the uniform theorem is equivalent to "`T^(m)_C` is not dense for some `m`"
(`conditional.md` Prop. 3.2(ii)). So within this family the weakest useful statement is the non-density
itself. Increasing `m` only weakens the required instance. There is no Vojta-type hypothesis that gives
growth without giving the uniform bound.

*Proof of (a).* The coordinate projection `pi : X_m -> X_(m')` (forget `z_(m'+1), ..., z_m`) is defined and
dominant on the torus, which contains every point of `T^(m)_C` (`conditional.md` Lemma 2.1(4)), and
`pi(T^(m)_C)` lies in `T^(m')_C`. If `Z'` is a proper closed set containing `T^(m')_C`, the closure of
`pi^(-1)(Z')` is a proper closed subset of the irreducible `X_m` containing `T^(m)_C`. QED

The Vojta threshold in SLW units is `1/(m-2)` (Remark 2.4c). It reaches the needed `1/4` at `m = 7`, one
point later than `SLW(5)`, because the canonical class does not see the rank.

### 6.2 Bombieri–Lang, and what the chord data really give

**Proposition 6.2 (Bombieri–Lang is vacuous for the natural auxiliary varieties).** `X_m`, its blow-up
along `Y`, the configuration spaces `M_(0,m+2)` (which carry the Ptolemy data, `growth/ptolemy.md`), and the
chord variety `V_m` below are all rational over `Q(i)`. Their canonical classes are not big, and
Bombieri–Lang says nothing about their rational points.

*Proof.* `X_m` is birational to `(P^1)^m` (`conditional.md` §7). `M_(0,n)` is an open subset of
`(P^1)^(n-3)`. `V_m` is isomorphic to an open subset of `A^(m-1)` (below). Blow-ups and open subsets of
rational varieties are rational. QED

**Proposition 6.3 (chord data are integrality data).** Let `z_1, ..., z_M` be a normalised
configuration and `rho_x = z_x/z_1`. Put `t_x = i(rho_x + 1)/(rho_x - 1)`, which lies in `Q` because
`|rho_x| = 1`. Then `rho_x = (t_x + i)/(t_x - i)`, and for `x != y`

```text
rho_x / rho_y = (T_xy + i)/(T_xy - i),        T_xy = (t_x t_y + 1)/(t_y - t_x).
```

Moreover the denominator of `T_xy` (in lowest terms) is at most `|beta_xy|/sqrt 2`, and every odd rational
prime dividing `beta_xy` divides it, where `beta_xy = (z_x - z_y)/g_xy` is the reduced chord. So
`t = (t_x)` is a `(1/D)`-integral point of the affine variety
`V_M = { (t_x, T_xy)_(2 <= x < y <= M) : T_xy (t_y - t_x) = t_x t_y + 1 }`, with `D = lcm_(x,y) denominator(T_xy)`.
Here `t_x = T_x1`. The rational points of `V_M` are those of `A^(M-1)` minus the diagonals, because
`t_x = t_y` forces `t_x^2 = -1`.

*Proof.* The Cayley formula is direct: `(t_x + i)(t_y - i) = (t_y - t_x)(T_xy + i)` and
`(t_x - i)(t_y + i) = (t_y - t_x)(T_xy - i)`. By the pair identity (`hardness.md` Lemma 3.1),
`z_x/z_y = u/(eps' conj u)` with `u = z_x/g_xy = A + Bi`, `gcd(A, B) = 1`, `A + B` odd, and `eps'` a unit,
and `beta_xy = u - eps' conj(u)`. Solving `(T+i)/(T-i) = u/(eps' conj u)` case by case gives:

| `eps'` | `T` | `beta_xy` |
| --- | --- | --- |
| `1` | `A/B` | `2Bi` |
| `-1` | `-B/A` | `2A` |
| `-i` | `(A-B)/(A+B)` | `(A+B)(1+i)` |
| `i` | `-(A+B)/(A-B)` | `(A-B)(1-i)` |

In each case the denominator is `|B|, |A|, |A+B|` or `|A-B|`, which is `|beta_xy|/2` or `|beta_xy|/sqrt 2`
(`gcd(A+B, A-B) = 1` since `A + B` is odd). QED

**Proposition 6.4 (the integrality set must grow).** In a configuration with `M` points, every inert
prime `l <= M - 2` divides the denominator `D` of Proposition 6.3. So the set `S` of places at which `t` fails
to be integral contains all inert primes below `M - 1`, and `|S| >= pi(M-2; 4, 3) ~ M/(2 log M)`.

*Proof.* The points lie in `l + 1` classes modulo `l` (`two_thirds` Lemma 3.3(I)). For `M > l + 1` two of
them coincide, so `l | beta_xy` (Lemma 3.3(I) again). By the table, `l` divides the denominator of `T_xy`
(it cannot divide `1 +- i`). QED

**Consequence.** Bombieri–Lang (rational points) does not apply. The cluster data are an
integrality condition on a rational variety. The relevant conjecture is Lang–Vojta for `S`-integral
points of a log-general-type pair: a partial compactification of `V_M`, essentially the `M_(0,M+2)`
picture of `growth/ptolemy.md`. For a fixed `S` that conjecture gives degeneracy, with an exceptional set
(and finitely many further points) depending on `S`. Proposition 6.4 shows that `S` must grow with `M`.
In near-extremal clusters `|beta_xy| <= C e^(mu_xy/2)` with pair margins `mu_xy = w_xy - W/2`. The cut
bound gives `sum_(x<y) mu_xy <= M W/4`, so on average `|beta_xy| <= C M^(K/4)` when `W = K M log M`. A bound
`M <= F(S)` from Lang–Vojta would need control of the exceptional sets as `S` grows. That is an
`S`-uniformity hypothesis of the same non-standard kind as in §5. So the chord-data route reduces to §5
and does not supply a standard hypothesis.

## 7. Evidence

Nothing in this section is used in a proof.

### 7.1 Known families (exact; `check_families_slw.py`)

`gamma_k := log(1/max_i |Lambda_i|) / sum_(j in U) log H(alpha_j)`, maximised over the `(k+1)`-subsets of rank
`k` and their base points. (SLW) with constant `c` predicts `gamma_k <= g + o(1)`.

| family | parameter | `log N'` | `log N'/log rad N'` | rank | `C` | max `gamma_k` |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| golden (`k = 4`) | `n = 4` | 20.10 | 1.087 | 4 | 5.897 | 0.2222 |
| | `n = 16` | 66.30 | 1.025 | 4 | 5.898 | 0.2403 |
| | `n = 34` | 135.59 | 1.024 | 4 | 5.898 | 0.2440 |
| CG generic (`k = 3`) | `a = 10^6 + 1` | 109.14 | 1.015 | 3 | 5.657 | 0.2428 |
| CG–Pell (`k = 3`) | `a = 681` | 50.83 | 1.213 | 3 | 5.657 | 0.2815 |
| | `a = 70711161` | 143.21 | 1.152 | 3 | 5.657 | 0.2808 |
| | `a = 22768774561` | 189.40 | 1.150 | 3 | 5.657 | 0.2776 |

The golden values approach `1/4` from below. The deficit `~ omega(N') log 2/log N'` comes from
`H(alpha_j) <= 2 p_j`. Asymptotically `gamma_4 >= (1/4)(log N'/log rad N') - o(1) >= 1/4 - o(1)`
(Prop. 2.5(a)). The CG–Pell values approach `(1/4)(8/7) = 2/7 = 0.2857`, and `log N'/log rad N'`
approaches `8/7` (Prop. 2.5(b)).

### 7.2 All circles up to `10^8` (`check_slw_data.py 1e8`, `check_rank_counts.py 1e8`)

There are 351,061 circles `N <= 10^8` (all prime factors `= 1 mod 4`) with at least 48 lattice points.
Every window of `k+1` angularly consecutive points of rank `k` was examined:

| `k` | windows of rank `k` | median `gamma_k` | 99.9% quantile | max `gamma_k` | Dirichlet `1/k` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 16,740,216 | 0.098 | 0.286 | 0.467 | 0.333 |
| 4 | 6,068,288 | 0.085 | 0.226 | 0.331 | 0.250 |
| 5 | 293,576 | 0.105 | 0.190 | 0.234 | 0.200 |

At `N <= 10^8` the constants (`c^(n+1)`, the arc constant, prime powers lowering `rad N`) dominate. The
maxima for `k = 3, 4` exceed the Dirichlet exponents for that reason, so these numbers cannot test an
asymptotic exponent. What they show is that rank-5 windows are much rarer and have smaller `gamma`.

Largest clusters with normalised `C <= 6`, by rank (circles with at least 24 points; counts are of
circles having such a cluster of rank `>= k`):

| decade of `N` | circles | rank `>= 3` | rank `>= 4` | rank `>= 5` | rank `>= 6` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `[10^4, 10^5)` | 1,127 | 252 | 10 | 0 | 0 |
| `[10^5, 10^6)` | 12,956 | 1,092 | 72 | 0 | 0 |
| `[10^6, 10^7)` | 142,080 | 4,526 | 288 | 6 | 0 |
| `[10^7, 10^8)` | 1,512,833 | 22,863 | 980 | 14 | 0 |

Rank-5 clusters with `C <= 6` exist at small `N`: 17 circles up to `10^8`, all with six points; the
largest normalised `log N'` among them is `18.3`. The counts grow more slowly with each extra rank
(ratios over the last decade: `5.1, 3.4, 2.3`). The random model predicts a factor `10^(-1/4)` per extra
rank per decade, times polylogarithmic factors. It also predicts that the total number of rank-5
clusters on arcs `C sqrt R` is finite (`sum_N r(N)^6 N^(-5/4) < infinity`) but grows up to astronomically
large `N` (the polylogarithmic factors win until `log_10 N ~ 10^2`). So the data are consistent with
`SLW(5)` and with the random model, and cannot test them.

### 7.3 Polynomial templates and rank

Write `z_sigma(u) = prod_(k=1..D) (u + c_k^(sigma_k))` with `c^(+1) = c`, `c^(-1) = conj c`. By `families.md`
Lemma P, bounded `C` forces `sum_k sigma_k Im(c_k^r) = const` for `1 <= r < D/2`. For large `u` the
factors have private primes, so the rank of the family is the affine rank of its sign vectors. That
rank is at most `D - rho`, where `rho` is the rank of the system `(Im c_k^r)_(k, r)`.
* `D = 4` gives rank `<= 3` (Cilleruelo–Granville).
* `D = 6` gives rank `<= 4`, unless all roots lie on one vertical line (`Re c_k` constant). Then only odd
  `r` count and the condition is `sum sigma_k b_k = const`. By Stanley's theorem on the Erdős–Moser
  problem, at most 5 subsets of 6 distinct positive integers have equal sums, so there are at most 5 points.
  (`vertical_templates.py`: no 6-set in `{1..20}` has a fibre of size `>= 6`.)
* A rank-5 six-point template therefore needs `D >= 8`. For vertical templates of degree 8 the conditions
  are equal sums and equal cube sums. `vertical_templates8.py` finds no fibre of size `>= 4` for 8-sets in
  `{1..18}`. The general degree-8 search of `families.md` (box `X <= 7`, `Y <= 2`) has largest bounded-`C`
  fibre 4.

So no rank-5 polynomial template is known. Whether one exists at larger degree is a
Prouhet–Tarry–Escott-type question that I have not settled. A positive answer would refute `SLW(5, g)` for
`g < 1/4`. Theorem 2.3 with `k = 6` would then still apply, with `SLW(6, g)` as the hypothesis.

## 8. Comparison and what remains

### 8.1 The ladder of hypotheses

| hypothesis | standard? | conclusion for actual circles | where |
| --- | --- | --- | --- |
| Vojta `V(7, 1/4)` on `Bl_Y X_7` | yes (Vojta 1987) | uniform, `M_0` ineffective | `conditional.md` Thm 5.1 |
| Vojta, one instance `V(m_0, eps < (m_0-6)/4)`, any `m_0 >= 7` | yes | uniform; larger `m_0` is weaker; no dimension uniformity needed | Prop. 6.1 |
| `SLW(5, g, kappa)`, any `g < 1/4` | no: simultaneous LW, Dirichlet exponent `1/5` | uniform, `M <= 458`, rank `<= 4` | Thm 2.3 |
| `SLW(k, g, kappa)`, `k >= 5`, `g < 1/4` | no (weaker as `k` grows) | uniform, `M <= (5.5^(k-1)+1)/2` | Thm 2.3 |
| `SLW^log(theta < 1/log 5.5, g < 1/4)` | no (implied by every `SLW(k_0, g)`) | growth, `M <= (log R)^(1-eta)` | Thm 2.4 |
| `SLW(4, g)`, `1/4 <= g < 4/15` | no | uniform | Prop. 2.7 |
| `SLW(4, g < 1/4)`, `SLW(3, g < 2/7)`, `SLW(1, g < 1/2)` | — | **false** (golden, CG–Pell, Liouville pairs) | Prop. 2.5 |
| Lang–Waldschmidt, single forms, any `eps`, `C_*` | yes (Lang 1978) | **no information** on skew double cores; does not suffice (given `exact.md` Thm 5.1) | Thm 3.3, Lemma 3.4, Cor. 3.5 |
| abc over `Q` or `Q(i)` | yes | only `rad(N) <= N^(1/4-eta)` circles; vacuous on squarefree `N`; does not suffice | Props 4.1–4.2, Cor. 4.3 |
| Subspace Theorem with known `S`-dependence | theorem | nothing (`exp(O(s))` subspaces) | §5.2 |
| `QS(m, delta, f = o(s))` | no | growth (`f` bounded: uniform) | Thm 5.2 |
| `QS(m, delta, f = o(s))` for `m <= 3` | — | **false** | Prop. 5.3 |
| Bombieri–Lang | yes | nothing (rational varieties) | Prop. 6.2 |
| Lang–Vojta for integral points, fixed `S` | yes | nothing (`S` must grow with `M`) | Props 6.3–6.4 |

### 8.2 Answer to the question of the round

1. **No standard named conjecture weaker than Vojta's has been shown to suffice.** Lang–Waldschmidt
   (single forms), abc, Bombieri–Lang, Lang–Vojta with fixed `S`, and the Subspace Theorem with any
   known `S`-dependence all fail. Each failure is proved:
   * for LW, through the LW-blind skew double core (conditional only on the draft fake-circle theorem of
     `exact.md` for the final "does not suffice" statement);
   * for abc, by vacuity on squarefree profiles;
   * for Bombieri–Lang and Lang–Vojta, by rationality and the growth of `S`;
   * for the Subspace Theorem, by Proposition 5.3 and the exponential known bounds.

   Within Vojta's conjecture, one instance on one fixed variety suffices. The required instance gets
   weaker as the number of points grows, and no dimension uniformity is needed. Vojta-type hypotheses
   cannot give growth without the uniform bound (Prop. 6.1(c)).
2. **The weakest hypothesis found that gives growth** is `SLW^log`: lower bounds for systems of at least
   `max(k_0, theta log(n+2))` simultaneous linear forms in `n` logarithms of Gaussian rationals on the unit
   circle, at height exponent `< 1/4`. It gives `M <= (log R)^(1-eta)`, and it is implied by each fixed-`k`
   hypothesis `SLW(k_0, g)`. **Among fixed-`k` hypotheses**, `SLW(k, g)` with `g < 1/4` gets weaker as `k`
   grows, and each one gives the uniform theorem. `k = 5` is the smallest `k` for which it can hold
   (Prop. 2.5). There it gives `M <= 458`, together with the structural conclusion that clusters have rank at
   most 4, which is exactly what every known infinite family has. Both are simultaneous strengthenings of Lang–Waldschmidt with a power to spare relative
   to the Dirichlet exponent. I do not claim they are named conjectures in the literature (Remark 2.1.1).
3. **The threshold "five forms" is sharp in both directions.**
   * Single forms (`k = 1`, the standard LW) are blind on the skew double core.
   * `k <= 4` forms at exponent `< 1/4` are false (golden family). `k = 4` survives only in the window
     `[1/4, 4/15)`.
   * `k = 5` is consistent with all known data. For `m = 6` points Vojta's conjecture is borderline, and
     the six-point families exist with rank `<= 3`.
4. **Strength relative to the known conditional theorem.** `V(7, 1/4)` and `SLW(5, g)` are not known to be
   comparable. They play parallel roles: each is a non-sharp Dirichlet-type statement with room to spare
   (Vojta: canonical threshold `2/5` against the needed `1/2` for 7 points; SLW: Dirichlet `1/5` against the
   needed `1/4` for 5 forms). SLW reaches the problem one point earlier (6 points instead of 7) because it
   sees the rank. Its conclusion is effective in `M` (`458`) though not in `R_0`. Vojta's conclusion is
   ineffective in both. The Dirichlet-sharp SLW for all `k` implies the full Cilleruelo–Granville
   `R^(1-eps)` conjecture (Remark 2.4b), as Vojta does.
5. **Honest assessment of `SLW(5, g)`.** By Remark 2.4a, the special case of `SLW(5, g)` that the proof uses
   says: six lattice points of rank 5 never lie on an arc of length `R^(1-2g+o(1))`. For `g` close to
   `1/4` this is only slightly stronger than the rank-5 case of the uniform theorem itself. The hypothesis
   has content beyond restating the problem for three reasons:
   * it is a special case of a uniform, Dirichlet-consistent transcendence statement about logarithms of
     arbitrary unit-circle Gaussian rationals;
   * it is needed only for rank-5 configurations, the smallest rank that no known family reaches, and
     not for large clusters;
   * Theorem 1.3 turns it into a bound for every cluster.
   The same caveat applies to `V(7, 1/4)`. The special case it uses is the non-density of `T^(7)_C`, a
   statement of the same nature as the uniform theorem (`conditional.md` Prop. 3.2).

### 8.3 Open points

* *`exact.md` Theorem 5.1* (draft) is the only unproved input of Corollary 3.5. If it fails, the LW barrier
  still holds at the level of Theorem 3.3 (no single cluster-determined form ever contradicts), but
  the existence of LW-consistent fake circles at `O(M log M)` would need another construction on an
  LW-blind profile.
* *`k = 2, 3`.* `SLW(2, g)` sees the skew double core for `g < 2/3` (Remark 3.6.3). Whether `SLW(2)` or
  `SLW(3)` at admissible exponents (`>= 2/7` by Prop. 2.5(b)) implies growth for all profiles is open. A
  fake circle on a profile in which every two (resp. three) independent characters have union support
  `>= W/(4g)` (resp. `>= 3W/(4g)`) would show it does not. Profiles built from `t >= 4` generic Hadamard bases
  are the natural candidates (the lead's construction). Their fake-circle property is not proved.
* *Rank-5 families.* A single infinite family of six-point configurations of rank 5 on arcs `C sqrt R`
  would refute `SLW(5, g)` for every `g < 1/4` (and move the threshold to `k = 6`). None is known (§7). Polynomial
  templates need degree `>= 8` for this (§7.3), which is a Prouhet–Tarry–Escott-type question. The template theorem of the notes
  bounds quadratic-unit templates. For golden templates with `b` blocks `J_n, ..., J_(n+b-1)`, two distinct
  words differ in angle by `|q|^n |D| (1+o(1))/phi` with `D = 2 sum_s (w_s - w'_s) q^s != 0`. `D` is nonzero
  because a `{0, +-1}` polynomial cannot vanish at `q = -phi^(-2)`: its conjugate root `-phi^2` has modulus
  `> 2`. Since `R ~ phi^(bn)`, the arc constant is of order `phi^(n(b/2 - 2))`, so `b <= 4`, which is rank `<= 4`.
  The exhaustive polynomial-template search of `families.md` §4 found only the six-point
  Cilleruelo–Granville templates, which have rank 3.
* *Non-character forms in Corollary 3.5* (Remark 3.6.1): a genericity extension of `outside.md`
  Theorem A.1 to LW inequalities.
* *Unconditional content.* Theorem 1.3 and Corollary 1.4 reduce the uniform theorem to bounding the
  rank, and growth to rank `<= theta log W`. No unconditional handle on the rank is known. The fake
  circles have full rank, so any such handle must use information beyond the local data.

## 9. Files (`round8/weakest_checks/`)

Every script runs from that directory with `python3 <script> [args]`. Its output is stored next to it.

| script | what it checks | type | output |
| --- | --- | --- | --- |
| `g8common.py` | exact Gaussian arithmetic, factorisation into split primes, exponent vectors, exact ranks, SLW exponents | library | — |
| `check_families_slw.py` | Prop. 2.5(a)(b): golden family (rank 4, `gamma_4`), generic CG (rank 3), CG–Pell (rank 3, `gamma_3 -> 2/7`); Pell identity asserted | exact arithmetic, 60-digit angles | `check_families_slw_out.txt` |
| `check_skew_support.py` | Lemma 3.1 (skew-Hadamard, balance, pair separation `M`) and Lemma 3.2 (support `>= (q+1)/2`; minimum `(q+7)/2`), `q = 7..43`; second part: union supports of two label-difference characters, `q = 43..151`; third part: the block identities of Lemma 3.4 steps 2–3 | exact | `check_skew_support_out.txt`, `check_skew_support_part2_out.txt`, `check_skew_support_part3_out.txt` |
| `vertical_templates.py B D`, `vertical_templates8.py B` | §7.3: equal-sum fibres of vertical templates (degree 6: sums; degree 8: sums and cube sums) | exact | `vertical_templates_out.txt`, `vertical_templates8_out.txt` |
| `check_subspace_count.py` | Prop. 5.3: pair and triple solutions of `Sigma_2, Sigma_3` for the first `s` split primes | exact construction, 60-digit verification | `check_subspace_count_out.txt` |
| `check_chord_cayley.py` | Prop. 6.3: Cayley coordinates, denominators against reduced chords, addition formula | exact (`Fraction`) | `check_chord_cayley_out.txt` |
| `check_slw_data.py X` | §7: SLW exponents of all windows of rank `k = 3, 4, 5` and ranks of the largest clusters, all `N <= X` | evidence (floating point) | `check_slw_data_out_X1e8.txt` |
| `check_rank_counts.py X` | §7: clusters of rank `>= 3, 4, 5, 6` with `C <= 6` per decade of `N` | evidence | `check_rank_counts_out.txt` |

The repository `/home/user/Jarnik` and the earlier round folders were not modified.
