# Referee report on `tau/verify.md` (the audit of the auxiliary-form reduction in `tau_problem.md`)

Referee: adversarial check, 2026-10-09. My scripts and logs are in
`scratchpad/tau/referee_verify_checks/`. The library `rlib.py` is written from scratch: it does
not import `taulib.py`. It uses different primes (`2^30 - 35`, `2^30 - 41` and twelve primes
below `2^29`) and a different coordinate projection (it drops `k_M`, not `k_1`).

Conventions. *Proved* means that I re-derived the argument line by line. *Certified* means an
exact computation over `Q`: an upper bound on `reg` from full rank modulo a prime, which implies
full rank over `Q`, together with a lower bound from an explicit integer measure checked in
Python integers. *Evidence* means anything else.

## 0. Verdict

**The result survives.** Every statement that `verify.md` labels Theorem, Lemma or Proposition
is correct; I re-derived each one. Every numerical claim I re-ran agrees, and I have now
certified each of them exactly over `Q`. `verify.md` itself had only modular upper bounds for
some of them.

One status claim must be corrected:

> "No auxiliary form on `X_3` or `X_4` can prove the theorem; `M >= 5` stays open."

The concurrent note `tau/upper.md` (Theorem 3.1) proves `reg A_{s,D} <= phi_M((D+s)/2, (D-s)/2)`.
I re-derived that proof line by line (§6) and found no gap. It gives

```text
tau(M) = tau_sh(M) = 2 - 1/ceil(M/2) < 2      for every M,
```

so **no auxiliary form on any `X_M` satisfies `m > 2D`.** The reduction is a correct but
vacuous implication. Theorem 4.5 of `verify.md` (`tau <= M/2`) is correct but superseded for
`M >= 4`; at `M = 3` the two bounds coincide. `verify.md` said explicitly that it had not
refereed `upper.md`, so this is a change of status, not an error in its mathematics.

Summary of findings:

| # | claim in `verify.md` | finding | severity |
|---|---|---|---|
| 1 | Reduction (Theorem 2.5), Taylor Lemma 2.1, Liouville Lemma 2.2, grid Lemma 2.4, small-`R` step | **Correct** (§1). Constants re-derived. | – |
| 2 | Shells do not decouple; counterexamples at `M = 2, 3` | **Correct.** Recomputed with exact rational rank (§2). | – |
| 3 | Strict-gain counts 36/72/12/3; `phi_M` matches 202 of 203 slices | **Correct, now certified.** `verify.md`'s `e1` logs are mod-`p` only, so the ball-side lower bounds for `M = 4, 5` and the "exact match" were evidence. I certified them exactly (§3). | minor (rigor) |
| 4 | `tau` = pseudo-effective threshold, `beta <= tau`, monotonicity, Fekete, `gamma_M` | **Correct** (§4). The "inequality reversed" remark about the lead's dimension count is interpretive. | minor |
| 5 | Theorem 4.5: `tau(M) <= M/2` | **Correct** (§5). Also tested adversarially on non-Vandermonde forms that mix several `s`. | – |
| 6 | "`M >= 5` stays open"; "equality of the shell and slice sups is not proved" | **Superseded** by `upper.md` Theorem 3.1, which I re-derived: `tau = tau_sh = 2 - 1/ceil(M/2)` for every `M` (§6). | major (status) |
| 7 | Novelty: "no project note formulates the fixed-form threshold"; the identification `tau = mu_eff` listed as new | **Overstated** (§7). The full-slice formula and the pseudo-effective threshold were already in `growth/conditional.md` Prop. 6.2 and `growth/referee_checks/tau_threshold.py`. A bound of the same strength, `M/2` from covering lines, is in `upper.md` Remark 3.3, written concurrently. | minor |

## 1. The reduction, re-derived

Notation as in `verify.md` §1. Write `F mod I(X_M) = sum_k c_k t^{(D-|k|_1)/2} u^{k+} v^{k-}`.
The coefficients `c_k` are unique, because these monomials form a basis of the degree-`D` part
of the coordinate ring. To check linear independence I used the torus chart
`u_j = u w_j`, `v_j = v/w_j`: the monomial `(e, k)` goes to `u^{e+|k+|} v^{e+|k-|} w^k`,
distinct pairs give distinct monomials, and the torus is dense because `X_M` is irreducible.
The `t = 0` part of the cone has dimension `M`, against `M + 1` for the cone itself.

1. **Real points and the trigonometric form.** At a real point, `u_j v_j = x_j^2 + y_j^2 = t`.
   If `t = 0`, all coordinates vanish, so `t > 0` and every real point lies in the open torus.
   On `|z_j| = R`, each basis term equals `R^{D-|k|} R^{|k|} e^{i k.theta}`, so
   `F(z) = R^D f(theta)`. Correct.
2. **Lemma 2.1 (Taylor bound).**
   * Put `g(tau) = f(theta^0 + tau h)` with `h_1 = 0`. Then
     `g^{(r)}(0) = i^r sum_s e^{i s theta_1} sum_{sum k = s} c_k (k.h)^r`.
   * `(k.h)^r` is a polynomial of degree `r < m` in `k`, so the inner sum vanishes by (b).
   * The integral remainder gives `|g(1)| <= sup|g^{(m)}|/m!`, which is valid for
     complex-valued `g`.
   * `|k.h| <= |k|_1 |h|_inf <= D delta`.
   * Hence `K = ||c||_1 D^m/m!`. Correct.
3. **Lemma 2.2 (Liouville).**
   * An arc of length `<= C sqrt R` has angular width `<= C R^{-1/2}`. Wrap-around is harmless,
     because Lemma 2.1 holds for all real `theta`.
   * `F(z)` lies in `Z[i]` when `F` has `Z[i]` coefficients in `x, y`, or in `u, v`.
   * The reduced coefficients `c_k` may lie in `Z[i][1/2]`. This does not matter: only
     `||c||_1` enters `K`, and integrality is read off the original `F`.
   * With `m >= 2D + 1` and `R >= 1`, `|F(z)| <= K C^m R^{-1/2} < 1` for `R > (K C^m)^2`.
     Correct.
4. **Lemma 2.4 (grid lemma).**
   * `e^{iD theta_1} f` is a polynomial of degree `<= 2D` in `e^{i theta_1}`.
   * Its coefficients `a_{k_1}` have `|k'|_1 <= D - |k_1| <= D`, so the induction closes with
     the same `n >= 2D + 1`. The base case is `M = 1`.
   * Liouville applies to tuples with repeated points, so `f = 0` on all of `Theta^M`.
   * Hence `|S| <= 2D`. Correct, and better than the lead's `2 deg F'`.
5. **Small `R`.** Consecutive lattice points on the arc are at chord distance `>= 1`. So the
   count is `<= floor(L) + 1 <= C sqrt(R_0) + 1 = K C^{m+1} + 1`. Correct.
6. **Arcs `C R^alpha`.** The bound becomes `|F| <= K C^m R^{D - m(1-alpha)}`, which is `< 1`
   eventually iff `alpha < 1 - D/m`. Since `m > 2D`, `1 - D/m > 1/2`. Correct.
7. **Remark 2.7 and Remark 2.8.** `Y` lies in the smooth locus: the singular locus of `X_M` is
   where at least two pairs `(u_j, v_j)` vanish, and that misses `Y`. So `P^m = P^{(m)}`
   locally, and the moment conditions are `Q`-linear. Correct.

*Evidence (`chkE2_taylor_hp.py`).* I checked Lemma 2.1 in 80-digit decimal arithmetic on
generic, non-Vandermonde kernel measures: `M = 4` (order 6), `M = 3` with `s = 2` (order 5),
and `M = 5` with `s = 1` (order 11). The worst value of `|f|/(K delta^m)` was `0.042`.

A caution for future checks: my first attempt used double precision (`chkE_gamma_taylor.py`)
and gave ratios around `1e12`. That came entirely from floating-point cancellation at
`delta < 1e-2`, not from a flaw in the lemma.

## 2. Shell decoupling: exact recomputation (`chkA_counterexamples.py`)

Ranks here are computed by exact `Fraction` elimination.

* **`M = 2`.** The measure `(1, -2, 1)` on `A_{0,2}` has exact order 2. The identity
  `F = -|z_1 - z_2|^2` on `X_2` was confirmed numerically. Exact values:
  `reg_Q(S_{0,2}) = 1`, `reg_Q(S_{0,0}) = 0`, `reg_Q(A_{0,2}) = 2`.
* **`M = 3`.**
  * `F_3` has order 3.
  * `F_3^2` has support 19, which is all of `A_{0,4}`, on the shells `t = 0, 2, 4`, and has
    exact order 6.
  * Its three shell parts have total masses `-6, 12, -6`, so each has order 0.
  * Exact values: `reg_Q(S_{0,4}) = 4`, `reg_Q(S_{0,2}) = 3`, `reg_Q(A_{0,4}) = 6`.

So `ord_Y` is not the minimum over shells, and at `(s, D) = (0, 4)` the slice beats every shell
in order. The underlying reason is right: in the torus chart, the exponents of `u` and `v` are
`(D+s)/2` and `(D-s)/2` whatever `|k|_1` is.

## 3. Slice data: recomputed and certified exactly (`chkB_gain.py`, `chkB_certify.py`, `chkD_extend.py`)

**Strict-gain slices** (`reg A_{s,D} > max_{t <= D} reg S_{s,t}`). My counts, from my own code
modulo two primes:

| `M` | range | count |
|---|---|---|
| 2 | `D <= 12` | 36 |
| 3 | `D <= 18` | 72 |
| 4 | `D <= 9` | 12 |
| 5 | `D <= 7` | 3 |

These are identical to `verify.md`, slice by slice. The listed examples also agree:
`M = 4`, `(4, 0)`: 6 vs 5; `M = 5`, `(7, 1)`: 11 vs 10.

**Rigor gap, now closed.** Modular rank bounds `reg_Q` only from above, so a count of strict
gains needs exact lower bounds as well. `verify.md`'s `e1` logs contain only mod-`p` values.
Exact lower bounds were justified there only for `M = 2` (collinear points) and `M = 3` (the
products `F_3^a (u_1 - u_2)^s`). I certified every slice exactly:

* **Gain slices.** An explicit integer measure on `A_{s,D}` of order `reg_p(A)`, verified in
  exact integers, so `reg_Q(A) = reg_p(A)`. For the shells, `reg_Q <= reg_p` suffices.
* **Non-gain slices.** An explicit exact measure, of order `reg_p(A)`, on the best shell.

All `48 + 99 + 29 + 19 = 195` slices are certified, and so are the six `M = 4`, `D = 10`
slices (`chkD_extend.txt`):

* `M = 2` is exact by hand. Slices and shells are collinear point sets, so `reg = |A| - 1`.
  This gives `reg A_{s,D} = D`, while shells have `reg` equal to `s` (if `t = s`) or 1. There
  is a strict gain iff `2 <= D` and `s < D`, which gives `sum_{D=2}^{12} floor(D/2) = 36`.
* `M = 3, 4, 5` were certified by `chkB_certify.py` (`chkB_cert_M{3,4,5}.txt`). The counts 36/72/12/3 are therefore exact
theorems for these ranges. In every certified slice `reg_Q(A_{s,D}) = floor(phi_M(p, n))`,
except `M = 5`, `D = 5`, `s = 1`, where `reg_Q = 7 < 8`. That is `verify.md`'s single
exception, now exact. No slice exceeds `phi_M`.

**Extension (`chkD_extend.py`, `chkD_extend.txt`).** Every slice below has an exact integer
measure for the lower bound and full rank modulo two primes for the upper bound, so `reg_Q`
is exact.

| slices | exact `reg_Q` | vs `floor(phi_M)` |
|---|---|---|
| `M = 5`, `D = 8`, `s = 0, 2, 4, 6, 8` | 13, 12, 12, 11, 8 | equal in all five |
| `M = 6`, `D = 4, 5, 6`, all `s` (10 slices) | – | equal in eight; strictly below at `(D, s) = (5, 1)`, 7 < 8, and at `(6, 2)`, 9 < 10 |

The `M = 5`, `D = 8`, `s = 0, 2` slices are the last two of `verify.md`'s 203. So **all 203
slices used in `verify.md` are now certified exactly over `Q`**: 202 equal `floor(phi_M)` and
one lies strictly below. In all, 216 slices are certified, and none exceeds `phi_M`. The best
`M = 5` ratio in this range is `13/8 < 5/3`.

## 4. Geometry (`verify.md` §4.1–4.4)

* **Prop. 4.1 (`tau` = pseudo-effective threshold).** Correct.
  * `X_M` is an irreducible complete intersection.
  * Its singular locus, where two or more pairs vanish, has codimension 3. So `X_M` is normal
    and projectively normal.
  * `Y` lies in the smooth locus, so `pi_* O(-mE) = I_Y^m`.
  * "Big + pseudo-effective is big" closes the sup argument.
* **Prop. 4.2 (`beta <= tau`).** Correct: `h^0(NH - mE) = 0` for `m > tau N`.
* **Monotonicity and Fekete.** Correct.
  * Products are nonzero because `X_M` is integral.
  * `ord_Y` is additive because it is a valuation: `Y` is smooth and lies in the smooth locus.
* **Prop. 4.4 (`gamma_M`).** Correct.
  * The count formula is right (multinomial times compositions), and so is the leading
    coefficient `binom(2n,n)/2^n`, by Vandermonde's identity.
  * `gamma_M >= 2M^{-1/(M-1)}` follows from `binom(2n,n) >= 4^n/(2 sqrt n) >= 4^n/(n+1)`.
  * Exact rational checks for `M = 2..400` (`chkE_gamma_taylor.py`) confirm
    `2M^{-1/(M-1)} <= gamma_M < 2`, `gamma_M < 2M/(M+1)` and
    `gamma_M <= 2 - 1/ceil(M/2)`; the last is consistent with §6.
  * The criticism that the lead's sentence "has the inequality reversed" is interpretive. On
    full slices, the generic count gives `gamma_M >= 2M^{-1/(M-1)}`, as `verify.md` says. In
    the lead's own shell framework, the generic count gives ratios tending to 0 as `t` grows.
    `verify.md`'s value is the correct one for the right invariant.

## 5. Theorem 4.5 (`tau(M) <= M/2`), re-derived and attacked (`chkC_curves.py`)

**(i) The curve lies on the cone.** `u_j v_j = prod_l (1 + (b_l + ig)e)(1 + (b_l - ig)e)`,
which does not depend on `j`. Correct.

**(ii) Contact of order 2 with `Y`.** The coefficient of `e` in `u_j` is
`sum_l b_l + ig(M-2)`, independent of `j`. Correct.

**(iii) The curve passes through a general point.**

* Solving `(1 + b + ig) = r (1 + b - ig)` gives `1 + b = -ig (1+r)/(1-r)`, and then
  `q = 1 + b - ig = -2ig/(1-r)`, which is nonzero.
* `phi(1)` is proportional to `(prod_{l != j} r_l, r_j)`.
* `r_j = mu v_j` with `mu^{M-2} = T/prod v` works for `M >= 3`.
* The admissible points (`r_j != 1`) form a dense open set.

**(iv) The bound.** Near `e = 0` the curve stays in the torus, and
`F(phi(e)) = sum_s u_1^{(D+s)/2} v_1^{(D-s)/2} P_s(w(e))` with `w(e) - 1 = O(e^2)`. Each `P_s`
is a Laurent polynomial that is analytic at `w = 1` and vanishes there to order `>= m`. So the
order at `e = 0` is `>= 2m`, while the degree is `<= MD` and the value at `e = 1` is nonzero.
Hence `2m <= MD`. Correct.

**Adversarial test.** I used exact generic kernel measures of maximal order, not Vandermonde
forms, including forms that mix two weights `s`: `M = 3, 4, 5`, `D` up to 7, orders up to 11.
They were composed with random exact Gaussian-rational curves from my own parametrisation. In
every case the curve has contact 2, `ord_{e=0} F(phi) >= 2m`, `deg <= MD` and `F(phi(1)) != 0`;
the order bound is attained with equality in all 14 trials. No counterexample.

## 6. Status of `M >= 5`: `upper.md` Theorem 3.1, re-derived

`verify.md` reports this concurrent claim without refereeing it. Because it decides the
bottom line, I checked its proof (`upper.md` §2–3) line by line.

* **Lemma 2.1 (peeling).** Given `g` on `A`:
  * interpolate `g` on `A_1` by `P_1`;
  * interpolate `(g - P_1)/L_1` on `A' = A \ Z(L_1)` by `P'`;
  * then `P = P_1 + L_1 P'` interpolates `g` on `A`, and
    `deg P <= max(reg A_1, reg A' + 1)`.

  Induction on `r` finishes the proof. Correct.
* **Lemma 2.2.** `Q^M(p,n) cap {k_1 = c}` is `Q^{M-1}(p-c, n)` for `0 <= c <= p`, and
  `Q^{M-1}(p, n-|c|)` for `-n <= c < 0`. Correct; the parity condition is automatic. `reg` is
  invariant under the affine identification.
* **Lemma 2.3 (sorted peeling).** `N(r_(i)) >= i`. Correct.
* **Theorem 3.1, induction on `M`.**
  * Base case: `Q^2(p,n)` is `p + n + 1` collinear points, so `reg = p + n = l^2_1`.
  * Case `j = 1`, with `K = p + mn`: the bounds `r^+_c <= V - c` and `r^-_c <= V` give
    `N(v) <= K - v + 1` in both sub-cases. Correct.
  * Case `j = m`: the case `j = 1` applied under `k -> -k`. Correct.
  * Case `2 <= j <= m-1`:
    - `1/alpha + 1/beta = 1`.
    - Identity (3.1): `V+/alpha + V-/beta = mp/(m+1-j) + mn/j = l^M_j(p,n)`; I verified it
      by hand.
    - Every slice value is `<= min(V+, V-)`, by monotonicity of `l_{j-1}` and `l_j`. This
      makes `x+, x- >= 0`.
    - Counting gives `N(v) <= floor(x+) + 1 + floor(x-) <= K - v + 1`.

    Correct.
* **Corollary 3.2.**
  * Even `M = 2k`: `l_k = (2 - 1/k)(p + n)`.
  * Odd `M = 2k + 1`: the mean of `l_k` and `l_{k+1}` is `(2 - 1/(k+1))(p + n)`.

  Correct.

I found no gap. `construct.md` reports an independent re-check. Numerically, none of my
216 certified slices (§3) exceeds `phi_M`.

**Lower bound.** The minimal Vandermonde measure lies on a **single shell** and has exact
order `C(M,2)`. I checked this exactly for `M = 2..6` and in closed form for `M <= 40`
(`chkF_vandermonde_shell.py`). Its ratio is `2 - 1/ceil(M/2)`.

**Consequences for `verify.md`.**

* `tau_sh(M) = tau(M) = 2 - 1/ceil(M/2)` for every `M`. The "equality of the overall sups ...
  is not proved" in `verify.md` §3 is now proved. The shell error in `tau_problem.md` is
  genuine slice by slice, but it does not change the value of the threshold.
* `tau(4) = 3/2` and `tau(5) = 5/3`, sharpening `3/2 <= tau(4) <= 2` and
  `5/3 <= tau(5) <= 5/2`.
* The hypothesis `m > 2D` of the reduction is unsatisfiable for every `M`. "`M >= 5` stays
  open" is false as a statement about the problem, though true of `verify.md`'s own methods.
* `beta_M <= tau(M) < 2` for every `M`. This answers the open sub-question of `conditional.md`
  §6 negatively; `upper.md` itself notes this.
* `upper.md` Remark 3.3 gets `ord <= min(p + (M-1)n, (M-1)p + n) <= MD/2` from lines through
  `e` in the Cremona-graph model. That is the same strength as Theorem 4.5, from a simpler
  covering family, written concurrently.

## 7. Novelty

* **Research notes.** Grepping `research/docs/` (1170 files) for pseudo-effective, Seshadri,
  covering family, movable curve, BDPP, contact order and `ord_Y` finds nothing that
  anticipates the fixed-form threshold or any bound `tau <= c`. `verify.md`'s claim holds
  for `research/docs`.
* **Growth files.** The claim is overstated for this project's `growth/` files:
  * `growth/conditional.md` Prop. 6.2, steps 2–3, already computes `ord_Y` as a minimum over
    `s` of orders on the full slices `A_s`. That is the correct form of the "shell
    correction". `verify.md` cites this in its summary, but §5 still lists "the shell
    correction" as new content.
  * `growth/referee_checks/tau_threshold.py`, dated before `verify.md`, has the docstring
    "Pseudo-effective threshold `tau(H,E)` on `Bl_Y X_M` = `max_s dstar(A_s)`". Its log gives
    `tau_N(3) = 1.5` at `N = 6`. So the identification `tau = mu_eff` on full slices was
    already formulated, as a numerical definition without proof. `verify.md` §5 says "No note
    defines `sup ord_Y F / deg F` on `X_M`".
* **What is new in `verify.md`.**
  * The written proof of Prop. 4.1.
  * The elementary Liouville proof without a distinctness factor.
  * `gamma_M`.
  * Theorem 4.5 with its contact-2 curve family. It was written concurrently with
    `upper.md`, which proves a stronger statement.

## 8. Issues

1. **Major (status).** "`M >= 5` stays open" and "the reduction's hypothesis is open for
   `M >= 5`" are superseded by `upper.md` Theorem 3.1, re-derived above. The correct
   statement is that `tau(M) = 2 - 1/ceil(M/2) < 2` for all `M`, so the reduction never
   applies.
2. **Minor (rigor of data).** The strict-gain counts and the "202 of 203 exact matches" with
   `phi_M` rested on mod-`p` ranks, which bound `reg_Q` only from above. They are now
   certified exactly (§3), and all agree.
3. **Minor (novelty).** The full-slice formula and the pseudo-effective-threshold formulation
   predate `verify.md` in `growth/conditional.md` and `growth/referee_checks/tau_threshold.py`.
   The `M/2`-type curve bound is also in `upper.md` Remark 3.3, concurrently.
4. **Minor (interpretation).** The "dimension count inequality reversed" remark depends on
   reading the lead's sentence as being about full slices. The value of `gamma_M` is correct.

## 9. Corrected statement

The reduction is a correct implication. Suppose `F`, a form of degree `D` with `Z[i]`
coefficients, is nonzero on `X_M` and has `m = ord_Y(F|X_M) > 2D`. Then:

* arcs of length `C sqrt R` carry at most `2D` lattice points once `R > (K C^m)^2`, where
  `K = ||c||_1 D^m/m!`;
* `M(C) <= max(2D, K C^{m+1} + 1)`;
* the same holds on arcs of length `R^alpha` for every `alpha < 1 - D/m`.

The proof is elementary (trigonometric Taylor bound plus a coefficientwise grid lemma).

The order along `Y` decouples by `s` but not by `l1`-shell, as the exact counterexamples at
`M = 2` and `M = 3` show. The right invariant is `tau(M) = sup reg(A_{s,D})/D`. It equals the
pseudo-effective threshold of `pi^*H - tE` on `Bl_Y X_M`, and `tau(M) >= beta_M`.

Covering curves of degree `M` with contact 2 give `tau(M) <= M/2`, hence `tau(3) = 3/2`.
However, by the concurrent `upper.md` Theorem 3.1, re-derived here, `reg A_{s,D} <= phi_M` and

```text
tau(M) = tau_sh(M) = 2 - 1/ceil(M/2) < 2      for every M,
```

so `tau(4) = 3/2` and `tau(5) = 5/3`. No auxiliary form on any `X_M` satisfies the hypothesis.
The reduction is valid but vacuous, and `M >= 5` is not open.

## 10. Reproducibility (`scratchpad/tau/referee_verify_checks/`)

| script | content | log |
|---|---|---|
| `rlib.py` | independent library: brute-force slices and shells, modular rank (primes `~2^30`), exact integer kernel via CRT and rational reconstruction, verified in Python integers, `phi_M` | – |
| `chkA_counterexamples.py` | §2: exact `Fraction` ranks | stdout reproduced in §2 |
| `chkB_gain.py`, `chkB_certify.py` | §3: strict-gain recount and exact certification | `chkB_M{4,5}.txt`, `chkB_cert_M{3,4,5}.txt`, `chkB_M*_D*.json` |
| `chkC_curves.py` | §5: curve family on generic kernel forms | stdout reproduced in §5 |
| `chkD_extend.py` | §3: remaining slices `M = 4`, `D = 10`; `M = 5`, `D = 8`; `M = 6`, `D <= 6` | `chkD_extend.txt` |
| `chkE_gamma_taylor.py`, `chkE2_taylor_hp.py` | §1, §4: `gamma_M` inequalities (exact); Taylor bound (80 digits) | stdout |
| `chkF_vandermonde_shell.py` | §6: Vandermonde measures sit on one shell, ratio `2 - 1/ceil(M/2)` | stdout |
