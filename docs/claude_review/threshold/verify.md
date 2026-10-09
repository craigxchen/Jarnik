# Verification of the auxiliary-polynomial reduction in `tau_problem.md`

Referee: independent audit, 2026-10-09. Scripts and logs are in `scratchpad/tau/`. My files are
`taulib.py`, `e1*.py`, `e2_exact.py`, `e3_curves.py`, `e4_liouville.py` and `e5_alt_ball.py`,
with their `*_output.txt` and `e1_*.txt`, `e5_*.txt` logs. The other files in that directory
(`upper.md`, `upper_checks/`, `hilb.py`, `peel.py`, ...) belong to a concurrent agent.

Conventions. *Theorem*, *Lemma* and *Proposition* mean proved here. *Data* means exact finite
computation. Data is never used as proof of an asymptotic statement.

## 0. Verdict

| item in `tau_problem.md` | verdict |
|---|---|
| **Reduction**: some `F` with `F|X_M != 0` and `ord_Y(F|X_M) = m > 2 deg F` implies the uniform theorem | **Correct** (Theorem 2.5). §2 gives a shorter, fully elementary proof with explicit constants. The conclusion is stronger than stated: an *absolute* bound `2D` for `R >= R_0(C)`, and bounded counts on arcs `C R^alpha` for every `alpha < 1 - D/m`, where `1 - D/m > 1/2`. |
| Liouville/Taylor step, uniformly on compact `X_M(C)`, including the torus-chart boundary and the real locus | Correct. All real points of `X_M` lie in the open torus, so the trigonometric chart covers them and the isotropic points of `Y` play no role (§2.6). The algebraic version also holds on all of `X_M(C)` (Remark 2.7). |
| Integrality of `F(z)` | Correct. Rational coefficients suffice: the conditions are `Q`-linear. |
| Homogeneity: `F` is nonzero on every product of circles `C_N^M` | Correct, and trivial in the trigonometric form: `F|C_N^M = N^{D/2} f(theta)`. |
| Grid lemma "with distinctness" | Correct but unnecessary. Liouville applies to tuples with repeated points too, so `F` itself vanishes on `S^M`. This gives `|S| <= 2D`, better than `2 deg F'` (Lemma 2.4). |
| Small-`R` step | Correct. |
| Decoupling by `s = sum k` | Correct (Lemma 3.1). |
| **Decoupling by the l1-shell `t`** ("the N-weights separate shells") | **False.** On `X_M` every shell carries the same weight `N^{D/2}`. Explicit counterexamples at `M = 2` and `M = 3` are in §3.2, and the per-slice failures are tabulated there. The correct invariant uses the full slices `A_{s,D} = {sum k = s, |k|_1 <= D, |k|_1 = D mod 2}`. The shell quantity `tau_sh(M)` defined in `tau_problem.md` is only a lower bound for the true threshold `tau(M)`. The implication "`tau_sh(M) > 2` implies the theorem" survives, because `tau_sh <= tau`. |
| `tau` = pseudo-effective threshold of `pi^*H - tE` on `Bl_Y X_M` | Correct for the true `tau(M)` (Proposition 4.1). For `tau_sh` only `<=` is proved. |
| `beta_M <= tau(M)` (Ru–Vojta) | Correct (Proposition 4.2). It is far from equality: at `M = 3`, `beta ~ 0.77` while `tau = 3/2`. |
| "Dimension count (generic) gives ratio `< 2 M^{-1/(M-1)}`" | Wrong direction. The generic count on full slices gives `tau(M) >= (binom(2n,n)/2^n)^{1/n}` (`n = M-1`). This lies in `[2M^{-1/(M-1)}, 2)` and tends to 2 (Proposition 4.4). |
| "nothing proves `tau <= 2` yet" | **New barrier (Theorem 4.5):** `tau(M) <= M/2` for every `M >= 3`. Hence `tau(3) = 3/2` exactly, `3/2 <= tau(4) <= 2`, and `5/3 <= tau(5) <= 5/2`. No auxiliary polynomial on `X_3` or `X_4` can prove the theorem. The first open case is `M = 5`. |

**Concurrent claim.** `tau/upper.md`, written by another agent this round and not refereed by
me, claims that `tau(M) = 2 - 1/ceil(M/2) < 2` for every `M`. That would make the reduction
vacuous. I did not check its proof. Its explicit slice bound
`phi_M(p,n) = (M-1) min_j (p/(M-j) + n/j)` (with `p = (D+s)/2`, `n = (D-s)/2`) does agree with
all my independent full-slice data:

* `reg(A_{s,D}) = floor(phi_M(p,n))` holds in 202 of the 203 computed slices (`M = 2..5`);
* the one exception (`M = 5`, `D = 5`, `s = 1`: 7 < 8) is strictly below the bound, which is
  consistent with `phi_M` being an upper bound;
* `phi` and my Theorem 4.5 agree at `M = 3`.

**Overall.** The reduction is a valid implication. It is the Cilleruelo–Córdoba (Ramana
determinant) Liouville argument with the Vandermonde replaced by an arbitrary form, plus a
Bezout grid lemma. Its combinatorial reformulation must use full slices, not shells. Its
hypothesis is provably unsatisfiable for `M <= 4`, and it is open, or false by the concurrent
claim, for `M >= 5`.

## 1. Setting

* `X_M = {u_1 v_1 = ... = u_M v_M}` in `P^{2M-1}`, where `u_j = x_j + i y_j` and
  `v_j = x_j - i y_j`.
* `Y` is the diagonal line, and `t = u_1 v_1`.
* For a form `F` of degree `D` write `F|X_M` in the monomial basis of the coordinate ring (a
  complete intersection, so it is projectively normal):

  ```text
  F  =  sum_k c_k  t^{(D-|k|_1)/2} u^{k+} v^{k-}   (mod I(X_M)),     |k|_1 <= D,  |k|_1 = D mod 2.   (1.1)
  ```

* On the real circles `z_j = R e^{i theta_j}`, with `N = R^2` and `u_j = z_j`,
  `v_j = conj(z_j)`, this reads

  ```text
  F(z)  =  sum_k c_k N^{(D-|k|_1)/2} prod_j z_j^{k_j^+} conj(z_j)^{k_j^-}  =  R^D f(theta),
  f(theta) = sum_k c_k e^{i k.theta}.                                                  (1.2)
  ```

  This is the lead's formula, with `z^k := prod z^{k+} conj(z)^{k-}`. As a *Laurent* polynomial
  the same term is `N^{(D-s)/2} z^k`, where `s = sum k`. The weight depends on `s`, not on `t`.
* Conversely, every finitely supported `c` with that support defines a form by (1.1). Its
  coefficients lie in `Z[i]` if the `c_k` do, and then `F(z)` is in `Z[i]` at Gaussian integers.

**Lemma 1.1 (three equivalent vanishing conditions).** For an integer `m >= 0` the following
are equivalent.

* (a) `ord_Y(F|X_M) >= m`, meaning the order at the generic point of `Y`.
* (b) For every `s`, and every polynomial `P` of degree `< m` on the hyperplane `{sum k = s}`,
  `sum_{sum k = s} c_k P(k) = 0`.
* (c) All partial derivatives of `f` of order `< m` vanish on the real diagonal
  `Delta = {theta_1 = ... = theta_M}`.

*Proof.*

* (a)⇔(b). Near the generic point of `Y`, use torus coordinates `u_j = u w_j`,
  `v_j = v/w_j`, `w_1 = 1`. Then (1.1) becomes `sum_s u^{(D+s)/2} v^{(D-s)/2} P_s(w)` with
  `P_s(w) = sum_{sum k = s} c_k w^k`. For different `s`, the monomials `u^a v^b` are
  independent over the completed local ring `C(v/u)[[w-1]]`, so `ord_Y = min_s ord_{w=1} P_s`.
  The order of a Laurent polynomial at `w = 1` is the least `|alpha|` with
  `((w d/dw)^alpha P_s)(1) = sum c_k k^alpha != 0`.
* (b)⇔(c). At `theta^0 = (phi, ..., phi)`,
  `d^alpha f(theta^0) = i^{|alpha|} sum_s e^{i s phi} sum_{sum k = s} c_k k^alpha`.
  This vanishes for all `phi` if and only if every inner sum vanishes (finite Fourier series
  in `phi`). `[]`

Also: `F|X_M != 0` if and only if `c != 0`, if and only if `f != 0`. Indeed the real torus
`(S^1)^M` is Zariski dense in `(C^*)^M`, and homogeneity moves any slice `t = N` to every other
slice.

## 2. The reduction

**Lemma 2.1 (Taylor bound on the real torus).** Assume (c) of Lemma 1.1. Then for every `theta`
in `R^M`,

```text
|f(theta)|  <=  K delta(theta)^m,      K = ||c||_1 D^m / m!,      delta(theta) = max_j |theta_j - theta_1|.
```

*Proof.* Put `theta^0 = (theta_1, ..., theta_1)` in `Delta`, `h = theta - theta^0` and
`g(tau) = f(theta^0 + tau h)`.

* `g^{(r)}(0) = sum_k c_k (i k.h)^r e^{i k.theta^0}`. This is a combination of the moments in
  (b) with `P(k) = (k.h)^r`, so it vanishes for `r < m`.
* By the integral form of Taylor's theorem, `|g(1)| <= sup_tau |g^{(m)}(tau)|/m!`.
* `|g^{(m)}| <= sum_k |c_k| |k.h|^m <= ||c||_1 (D |h|_inf)^m`, since `|k|_1 <= D`. `[]`

No chart, compactness argument or ideal-membership statement is needed: the trigonometric chart
covers all of `X_M(R)` (§2.6).

**Lemma 2.2 (Liouville).** Let `c` have entries in `Z[i]` (so `F` has `Z[i]` coefficients),
`c != 0`, and assume (b) holds with `m >= 2D + 1`. Let `R^2 = N` be an integer and let
`z_1, ..., z_M` be lattice points, **not necessarily distinct**, on one arc of length
`<= C sqrt R` of `|z| = R`. If `R > R_0 := (K C^m)^2`, then `F(z) = 0`.

*Proof.*

* Choose arguments `theta_j` in an interval of length `<= C sqrt R / R`. Then
  `delta <= C R^{-1/2}`.
* `F(z)` is in `Z[i]`, and by Lemma 2.1, `|F(z)| = R^D |f(theta)| <= K C^m R^{D - m/2}`.
* Since `m >= 2D + 1`, this is at most `K C^m R^{-1/2}`, which is `< 1`. `[]`

**Lemma 2.3 (non-vanishing on every circle).** `F|C_N^M = N^{D/2} f`. So `F` is not
identically zero on `C_N^M` for any `N > 0` if and only if `c != 0`. (The lead's homogeneity
argument is correct; (1.2) makes it a one-liner.)

**Lemma 2.4 (grid lemma, no distinctness needed).** Let `Theta` be a set of `n >= 2D + 1` angles
that are distinct mod `2 pi`. If `f` vanishes on `Theta^M`, then `c = 0`.

*Proof.* Induct on `M`.

* For fixed `(theta_2, ..., theta_M)` in `Theta^{M-1}`, the function
  `e^{i D theta_1} f` is a polynomial of degree `<= 2D` in `e^{i theta_1}`.
* It has `n >= 2D + 1` zeros, so each coefficient `a_{k_1}(theta_2, ..., theta_M)` vanishes on
  `Theta^{M-1}`.
* Each `a_{k_1}` is a trigonometric polynomial with `|k'|_1 <= D`, so the induction hypothesis
  applies to it. `[]`

The lead's version multiplies `F` by `prod_{i<j} |z_i - z_j|^2` to restrict to distinct tuples.
It is also valid: the Bezout degree in each pair is `D + 2(M-1)`, which gives
`|S| <= 2(D + 2M - 2) <= 2 deg F'`. It is unnecessary, because Lemma 2.2 already covers repeated
points.

**Theorem 2.5 (the reduction, verified).** Suppose there are `M`, `D`, `m >= 2D + 1` and
`c != 0` with entries in `Z[i]` (equivalently in `Q(i)` or `C`, see Remark 2.8) satisfying (b).
Then for every `C > 0` and every arc of length `<= C sqrt R` on `x^2 + y^2 = R^2` (`R^2` in `Z`,
any position):

```text
#(lattice points on the arc)  <=  2D                     if R > R_0(C) = (K C^m)^2,
#(lattice points on the arc)  <=  K C^{m+1} + 1          otherwise.
```

So `M(C) <= max(2D, K C^{m+1} + 1)`, with `K = ||c||_1 D^m/m!`.

*Proof.*

* Let `S` be the set of lattice points on the arc, with angle set `Theta`.
* If `R > R_0`, Lemma 2.2 gives `F = 0`, hence `f = 0`, on `Theta^M`. Lemma 2.4 gives
  `|S| <= 2D`.
* If `R <= R_0`, distinct lattice points are at distance `>= 1`. So an arc of length `L` holds
  at most `L + 1` of them, and `L + 1 <= C sqrt(R_0) + 1 = K C^{m+1} + 1`. `[]`

**2.6 Audit of the lead's step 1 and the "real-locus subtlety".**

* Every real point of `X_M` has `u_j v_j = |z_j|^2 = t > 0`. So `X_M(R)` is the compact torus
  `(S^1)^M/±1`, inside the open torus of `X_M`.
* The points of `Y` outside the torus are the two isotropic points `(x:y) = (1:±i)`. They are
  non-real, and real tuples never approach them.
* The Taylor bound of step 1 is therefore exactly Lemma 2.1, and it is uniform in the arc
  position because the constant `K` does not depend on `theta`.
* The hypothesis is the complex algebraic order. It equals the real trigonometric order
  (Lemma 1.1), and real vanishing to a higher order is irrelevant because only an upper bound for
  `|F|` is used.

**Remark 2.7 (algebraic form).** The lead's bound also holds on all of `X_M(C)`, isotropic points
included:

```text
|F(p)| <= K' ||p||^{D-m} max_j |l_j(p)|^m,      l_j in {x_j - x_1, y_j - y_1}.
```

* `Y` is smooth and lies in the smooth locus of `X_M` (`conditional.md`, Lemma 2.1(4)).
* In a regular local ring `A` with prime `P` such that `A/P` is regular, `gr_P A` is a
  polynomial ring over `A/P`. So `P^m` is `P`-primary, and generic order `>= m` gives
  `F in I_Y^m` at *every* point of `Y`.
* Compactness of `Y(C)` and of the unit sphere of the cone then gives `K'`.

This step is correct, but Lemma 2.1 makes it unnecessary.

**Remark 2.8 (coefficients).** Condition (b) is a system of `Q`-linear equations in `c` (the
moments `sum c_k k^alpha`). If it has a nonzero complex solution, it has a nonzero rational one,
and then a `Z[i]` one after clearing denominators.

**Remark 2.9 (strength: a sanity flag).** Theorem 2.5 gives more than the target.

* Absolute bound: `2D` points for `R >= R_0(C)`, whatever `C` is.
* The same proof with arcs of length `C R^alpha` gives `|F(z)| <= K C^m R^{D - m(1-alpha)}`. So
  arcs of length `C R^alpha` carry at most `2D` points for large `R` whenever
  `alpha < 1 - D/m`, and `1 - D/m > 1/2`.

So `tau(M) > 2` would imply bounded counts on arcs `R^{1/2+delta}`. Vojta predicts this
(`conditional.md`, Theorem 5.2), but it is far beyond current knowledge. This is consistent
with the barrier results of §4.

* The Pell family (8 points on arcs of fixed `C`) forces `D >= 4` for any such `F`.
* The Cilleruelo–Granville 10-point family on arcs `~R^{3/5}` forces `2D >= 10` if `m/D > 5/2`.

Neither is a contradiction.

**Data (E4, `e4_liouville.py`).**

* On 156 sub-tuples of the Cilleruelo–Granville six-point clusters (`a` up to 3162,
  normalized arc length `C ~ 8`), `F(z)` is a Gaussian integer.
* The Ramana identity `Norm(F) N^{s'^2} = prod_{i<j} Norm(z_i - z_j)` holds exactly for `F_3`
  and `F_5`.
* The bound of Lemma 2.1 holds with margins `e^{-3.47}` (`M = 3`) and `e^{-15.8}` (`M = 5`).
* `F` vanishes at tuples with a repeated point.

## 3. The combinatorial reformulation

**Lemma 3.1 (decoupling by `s`, correct).**
`ord_Y(F|X_M) = min over s with c|_s != 0 of ord(c|_s)`, by the proof of Lemma 1.1. Hence

```text
tau(M) := sup_F ord_Y F / deg F  =  sup_{D,s} reg(A_{s,D}) / D,
A_{s,D} = {k in Z^M : sum k = s, |k|_1 <= D, |k|_1 = D mod 2},
```

where `reg(A) = min{d : polynomials of degree <= d interpolate A}`. This equals the largest
order of a nonzero measure on `A`.

* The quantity `r(D) = max_s reg(A_{s,D})` is superadditive: multiplying forms convolves
  measures, and orders and degrees both add.
* So `tau(M) = lim r(D)/D` (Fekete). Finite-`D` data give only lower bounds for `tau`, and
  mod-p data give upper bounds for `r(D)` at that `D`.

**3.2 Decoupling by shells is false.** The shell `t` enters (1.1) only through the factor
`t^{(D-t)/2} = (u_1 v_1)^{(D-t)/2}`, which is `(uv)^{(D-t)/2}` near `Y`. It is a unit with no
`w`-dependence and cannot separate the shells.

* **`M = 2`, `(s, D) = (0, 2)`.** The measure `delta_{(1,-1)} - 2 delta_0 + delta_{(-1,1)}`
  on `A_{0,2}` corresponds to `F = u_1 v_2 + u_2 v_1 - 2 u_1 v_1 = -|z_1 - z_2|^2`. It has
  order 2. The shells at `s = 0` are `S_{0,0} = {0}` (order 0) and
  `S_{0,2} = {±(1,-1)}` (order `<= 1`).
* **`M = 3`, `(s, D) = (0, 4)`.** `F_3^2`, where `F_3` is the alternating orbit sum of
  `(1,0,-1)`, has a measure with 19 points on the shells `t = 0, 2, 4`. It has exact order 6.
  Each of its three shell components separately has order 0. Every measure on the single shell
  `S_{0,4}` has order `<= 4`: `reg_p = 4` at two primes, and `rank_p <= rank_Q` makes this a
  rigorous upper bound. So `ord` does not decompose over shells.
* **Systematic data (E1).** The number of slices `(s, D)` where `reg(A_{s,D})` exceeds every
  `reg(S_{s,t})` with `t <= D`:

  | `M` | slices with strict gain | range |
  |---|---|---|
  | 2 | 36 | `D <= 12` |
  | 3 | 72 | `D <= 18` |
  | 4 | 12 | `D <= 9` |
  | 5 | 3 | `D <= 7` |

  Examples: `M = 4`, `(D,s) = (4,0)`: ball 6 vs shells 5. `M = 5`, `(7,1)`: ball 11 vs shells 10.

**Consequence.** The quantity `tau_sh(M)` defined in `tau_problem.md` (sup over single shells
`S_{s,t}` of `ord/t`) satisfies `tau_sh(M) <= tau(M)`. A shell measure on `S_{s,t}` *is* a form
of degree `t`. Whether `tau_sh = tau` is not proved.

* In the computed range the overall maxima agree:

  | `M` | `D` range | max ratio, shells | max ratio, slices |
  |---|---|---|---|
  | 2 | `<= 12` | 1 | 1 |
  | 3 | `<= 18` | 3/2 | 3/2 |
  | 4 | shells `<= 9`, slices `<= 10` | 3/2 | 3/2 |
  | 5 | shells `<= 7`, slices `<= 8` | 5/3 | 5/3 |

* Antisymmetric measures (E5) agree as well, up to `t_min + 10` (M=5), `+8` (M=6), `+6` (M=7):
  maxima `5/3`, `5/3`, `7/4`, attained only by the Vandermonde measure.
* My recomputation of the lead's shell numbers for `M = 2..4` (`t <= 9`) and `M = 5`
  (`t <= 7`) agrees with `tau_problem.md`.
* For `M = 3`, `reg(A_{s,D}) = (3D - s)/2` exactly in all 99 slices with `D <= 18`. The upper
  bound comes from mod-p rank. The lower bound comes from `F_3^{(D-s)/2} (u_1 - u_2)^s`.

**Correction needed in `tau_problem.md`.** Replace the shells `S_{s,t}` by the slices
`A_{s,D}`, and the denominator `t` by `D`. The sentence "the N-weights separate shells" should be
deleted. Question (C)'s "several shells with N-dependent coefficients" is already inside the
reduction when the coefficients are constant. With `N`-dependent coefficients the form is not
fixed and Lemma 2.1 loses uniformity, unless the dependence is polynomial in `N`. That case is
again a homogeneous form, of higher degree.

## 4. Geometry: pseudo-effective threshold, beta, and a barrier for M <= 4

**Proposition 4.1.** `tau(M) = sup{t : pi^*H - tE is pseudo-effective on Bl_Y X_M}`.

*Proof.*

* `H^0(Bl, D pi^*H - mE) = {F in H^0(X_M, O(D)) : ord_Y F >= m}`. This uses two facts:
  `pi_* O(-mE) = I_Y^m` (Y is smooth in the smooth locus), and
  `H^0(X_M, O(D)) =` forms of degree `D` mod `I(X_M)` (projective normality).
* So `tau` is the sup over `Q`-effective classes.
* If `pi^*H - t_0 E` is pseudo-effective and `t < t_0`, then
  `pi^*H - tE = (1 - t/t_0) pi^*H + (t/t_0)(pi^*H - t_0 E)` is big + pseudo-effective, hence
  big, hence `Q`-effective. The two sups therefore coincide. `[]`

For the shell quantity, only `tau_sh <= tau` is available.

**Proposition 4.2.** `beta_M = beta(H, E) <= tau(M)`.

*Proof.* `h^0(NH - mE) <= h^0(NH)`, and `h^0(NH - mE) = 0` for `m > tau N`. So
`sum_{m >= 1} h^0(NH - mE) <= tau N h^0(NH)`. `[]`

The inequality is far from equality: `conditional.md` measured `beta_3 ~ 0.77` (at `N = 24`),
while `tau(3) = 3/2` (Theorem 4.5). Both `beta` and `tau` measure approximation to the rational
line `Y` at a single place. For a *rational* target, Liouville with the pseudo-effective
threshold is stronger than the Subspace/Ru–Vojta route. This is already the case on `P^1`
with a rational point: `tau = 1` gives Liouville's exponent, while `beta = 1/2` gives only Roth's.

**Proposition 4.3 (monotonicity).** `tau(M+1) >= tau(M)`. A form in the first `M` points has
the same order along `Y`, by Lemma 1.1(b), since moments in fewer coordinates are unchanged. It
is nonzero on `X_{M+1}`, because `X_{M+1} -> X_M` is dominant. With the Vandermonde measures,
`sup_M tau(M) >= lim 2M/(M+1) = 2`.

**Proposition 4.4 (generic count, corrected).**
`tau(M) >= gamma_M := (binom(2n,n)/2^n)^{1/n}`, with `n = M - 1`.

*Proof.*

* Counting by the positive part gives the exact formula
  `#{k : sum k = 0, |k|_1 = 2a} = sum_{p,q >= 1} M!/(p! q! (M-p-q)!) binom(a-1,p-1) binom(a-1,q-1)`.
* Hence `n! |A_{0,D}|/D^n -> 2^{-n} sum_p binom(M,p) binom(M-2,p-1) = binom(2n,n)/2^n`, by
  Vandermonde's identity.
* A nonzero measure of order `m` exists once `binom(m-1+n, n) < |A_{0,D}|`. `[]`

Values for `M = 3..12`: 1.225, 1.357, 1.446, 1.511, 1.560, 1.600, 1.632, 1.659, 1.681, 1.701.

* `gamma_M >= 2M^{-1/(M-1)}`. This is the averaged-slice bound; the lead's sentence has the
  inequality reversed.
* `gamma_M < 2`, since `binom(2n,n) < 4^n`, and `gamma_M -> 2`.
* `gamma_M` is always below the Vandermonde value `2M/(M+1)`.
* E2 checks the exact count formula against enumeration (`M = 3, 4, 5`, `D <= 8`) and the
  limit at `D = 4000`.

**Theorem 4.5 (new barrier).** For every `M >= 3`, `tau(M) <= M/2`. Consequently:

* `tau(3) = 3/2`;
* `3/2 <= tau(4) <= 2`;
* `5/3 <= tau(5) <= 5/2`.

In particular, no form on `X_3` or `X_4` satisfies the hypothesis `m > 2D` of the reduction.

*Proof.* Let `r` be in `(C^*)^M` with all `r_j != 1`. Fix `g != 0` and define `b_j` by
`(1 + b_j + ig)/(1 + b_j - ig) = r_j`. Consider the polynomial curve

```text
u_j(e) = (1 + (b_j - ig) e)  prod_{l != j} (1 + (b_l + ig) e),
v_j(e) = (1 + (b_j + ig) e)  prod_{l != j} (1 + (b_l - ig) e),          j = 1..M.
```

(i) **The curve lies in the cone of `X_M`.** Each product `u_j v_j` equals
`prod_l (1 + (b_l + ig) e)(1 + (b_l - ig) e)`, independent of `j`.

(ii) **Contact of order at least 2 with `Y` at `e = 0`.**

* `phi(0) = (1, ..., 1)` lies in `Y` and inside the torus, so it is not a base point.
* The coefficient of `e` in `u_j` is `sum_l b_l + ig(M-2)`, independent of `j`. The same holds
  for `v_j`.
* So `w_j(e) := u_j(e)/u_1(e) = 1 + O(e^2)`, and likewise for the `v`'s.

(iii) **The curve passes through a general point at `e = 1`.**

* Put `q_l = 1 + b_l - ig`, so that `1 + b_l + ig = r_l q_l`. Then
  `(u_j(1), v_j(1)) = (prod q) (prod_{l != j} r_l, r_j)`.
* Given a torus point `x = (u_j, v_j)` with `u_j v_j = T`, choose `mu` with
  `mu^{M-2} = T/prod v` and put `r_j = mu v_j`.
* Then `prod_{l != j} r_l = mu^{M-1} prod_{l != j} v_l = mu u_j`, so `phi(1)` is proportional
  to `x`.
* This covers a nonempty Zariski-open set of points (`r_j != 1`). The log-Jacobian of
  `r -> (u_j/v_j)_j` is `J - 2I`, with determinant `(M-2)(-2)^{M-1} != 0` for `M >= 3`.

(iv) **The bound.** Let `F` have degree `D`, `F|X_M != 0`, and `ord_Y F >= m`. Pick such an `x`
with `F(x) != 0`.

* `F(phi(e))` is a polynomial in `e` of degree `<= MD`, and it is nonzero at `e = 1`.
* Near `e = 0`, (1.1) and the torus chart give
  `F(phi(e)) = sum_s t(e)^{(D-s)/2} u_1(e)^s P_s(w(e))`.
* Each `P_s` vanishes to order `>= m` at `w = 1` and `w - 1 = O(e^2)`, so
  `ord_{e=0} F(phi(e)) >= 2m`.
* Hence `2m <= MD`. `[]`

For `M = 3` the upper bound is attained by `F_3` (order 3, degree 2). Remarkably,
`F_3(phi(e)) = c e^6` exactly, which was verified on an actual lattice triple.

In geometric terms: the curves form a covering family with `pi^*H . Gamma <= M` and
`E . Gamma >= 2`, so `pi^*H - tE` is not pseudo-effective for `t > M/2`.

*Exact checks (E3, `e3_curves.py`, Gaussian-rational arithmetic).*

* For `M = 3..7` with random `r` and `g`: the curve lies on the cone, its contact order is 2,
  its degree is `M`, and `phi(1)` is proportional to the prescribed point.
* `F_3` (and `F_5` for `M >= 5`) composed with `phi` has order exactly `2m` and degree `<= MD`.
* Through an actual lattice triple on `N = 77068225`, `F_3 o phi = c e^6`.

**4.6 Beyond `M = 4`.** The same construction is limited for larger `M`.

* *Contact 2.* Sign-vector (Cilleruelo–Granville type) curves
  `prod_j (L_j + i sigma_j K_j)` need `rank Sigma = M` to pass through a general point. So
  `d >= M`, and contact 2 gives only `d/2 >= M/2`.
* *Higher contact.* Contact `kappa` adds `(kappa - 1)(M - 1)` conditions on `2d` parameters.
  The resulting dimension count is heuristic and fails at `M = 3`, where it would contradict
  `tau(3) = 3/2`. So it cannot be trusted to decide `M >= 5`.

The question "`tau(M) > 2` for some `M >= 5`?" is not settled here. The concurrent `upper.md`
claims "no", for all `M`.

## 5. Novelty

**Project notes** (`research/docs`, 1170 files; grep for auxiliary, determinant, Ramana,
Vandermonde, Seshadri, threshold, "order of vanishing", pseudo-effective, Liouville,
diagonal + vanish, blow-up).

* No note defines `sup ord_Y F / deg F` on `X_M`. None uses the pure Liouville bound with an
  arbitrary form, mentions Seshadri or pseudo-effective thresholds, or proves `tau <= 2`.
* *Special case in the notes.* The Ramana determinant
  (`codex_uniform_bound_research.md` §5, `bounded_split_support_reduction.md`,
  `inverse_valuation_profile.md`) is exactly the Vandermonde instance of the reduction: ratio
  `(2s'+1)/(s'+1)`, with `Norm(J_s) N^{s^2} = prod Norm(z_j - z_i)`, checked in E4.
* *Different mechanisms.* Several notes treat vanishing-order budgets, all for different
  mechanisms:
  - `hermitian_polynomial_barrier.md` and `plucker_compatibility_route.md`: averaged orders
    along isotropic or small-diagonal subspaces in half-angle or bracket coordinates, combined
    with *divisibility*;
  - `polarized_coefficient_support_multiplicity.md` and `gale_torus_newton_cut_orders.md`:
    Gale pullbacks at eight points;
  - `six_point_descent_osculating_hyperplane.md`: a one-point Liouville inequality for a
    specific hyperplane.
* *This round's write-ups.* `growth/auxiliary.md` Prop. A covers only interpolation
  determinants (`V x Phi`). `growth/conditional.md` has `X_M`, `Y`, the Bezout lemma (with
  distinctness) and `beta`, but uses Ru–Vojta rather than Liouville. The reduction is therefore
  new to the project as a formulation, but it is the direct generalization of a mechanism the
  notes already use.

**Literature.** I used web search only; arXiv and university servers were unreachable by
direct fetch.

* Cilleruelo–Córdoba, "Trigonometric polynomials and lattice points", Proc. AMS 115 (1992).
  Citing papers describe its key tool as the product-of-distances inequality for lattice points
  on a circle. That is `|F_V(z)| >= 1` for the Vandermonde form, i.e. the `F = F_V` case of
  Lemma 2.2, with the classical exponent `1/2 - 1/(4 floor(k/2) + 2)`. I could not read the
  paper itself.
* The general principle — a single section vanishing to high order along a rational subvariety
  gives a Liouville-type approximation bound off its base locus, governed by an effective or
  pseudo-effective threshold — is standard. See McKinnon–Roth's analogue of Liouville's theorem
  with the asymptotic base locus (Eur. J. Math. 2016, arXiv:1306.2977) and work on
  approximation constants of closed subschemes.
* Bombieri–Pila and Heath-Brown determinant methods build point-dependent interpolation
  determinants. `growth/auxiliary.md` Prop. A shows these reduce to Vandermonde times a
  symmetric factor here. The fixed-form `tau` question is different from them.
* I found no published computation of the pseudo-effective threshold of `Bl_Y X_M`, and no
  statement equivalent to Theorem 4.5.

**Assessment.** The reduction is correct but standard in kind: Cilleruelo–Córdoba/Liouville
plus Bezout. Its value is that it isolates one clean invariant, `tau(M)`. The new rigorous
content of this audit is:

* the shell correction;
* the identification `tau = mu_eff >= beta`;
* the generic count `gamma_M`;
* Theorem 4.5: `tau <= M/2`, `tau(3) = 3/2`, and no use of `X_3` or `X_4`.

## 6. Reproducibility

All scripts run with `python3` (numpy 2.x) from `scratchpad/tau/`.

| script | what | output |
|---|---|---|
| `taulib.py` | slices and shells, modular rank (two 31-bit primes), exact kernel certification | – |
| `e1_shell_ball.py`, `e1b.py`, `e1c.py` | `reg` of shells vs full slices, two primes. `reg_p >= reg_Q`, so these are upper bounds. | `e1_M{2..5}_*.txt` (no prime mismatches) |
| `e2_exact.py` | §3.2 counterexamples, Vandermonde orders (`M <= 6`), exact count formula, monotonicity | `e2_exact_output.txt` |
| `e3_curves.py` | Theorem 4.5 curve family (exact `Q(i)`) | `e3_curves_output.txt` |
| `e4_liouville.py` | Liouville sanity checks and the Ramana identity on CG clusters | `e4_liouville_output.txt` |
| `e5_alt_ball.py` | antisymmetric measures: shells vs full slices, `M = 5, 6, 7` | `e5_M{5,6,7}.txt` |

Lower bounds on `reg` that are claimed as exact come from explicit measures (Vandermonde,
products), verified in exact integer arithmetic. Upper bounds come from full rank modulo a prime,
which implies full rank over `Q`.
