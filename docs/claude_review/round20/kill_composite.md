# Kill report: `round20/composite.md` (composite-supported analytic inputs)

*Adversarial killer, round 20. Scripts and outputs: `kill_composite_checks/` (absolute path
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_composite_checks/`). All my
scripts are independent of the proposal's `composite_checks/common.py`. Labels: **proved** (re-derived here),
**exact**, **evidence**, **heuristic**.*

## 0. Verdict

**KILLED as stated. A narrower statement survives.**

*What survives (re-derived line by line; no mathematical error found):* Lemmas 1.2-1.4, Prop. 1.5, Theorem M, Theorem P
(P1, P2, P2', P3, P4, P5) and Prop. 3.1. Corollary 2.2 also survives, with one constant corrected (F1). So
**`(c, gamma)`-robust composite inputs with `gamma` fixed before `c_0`** are capped at the rate `log R/loglog R`.
This covers inputs used without cross multiplicativity (minimal graft) and, with multiplicativity, inputs in the
forms P1-P5. P1 is a correct answer to round-10 referee Slip 3. This part is a cap; it gives no growth.

*What is killed:*
1. **The headline "composite support adds no separating power"** and **"the only non-robust leftover is the sharp large
   sieve on primitive support"** are false as stated (§5). Two non-robust composite inputs are fake-violated and are
   not shown to collapse:
   - (F3) smoothed-radius sector counts with badly approximable edges. They are true with error `O(1)` to `O(log X)`
     (evidence to `X = 10^7`) and are violated by the minimal graft by a factor of about 16 (evidence).
   - (F4) smooth all-ideal Hecke sums used **without** multiplicativity. Prop. 3.1 needs complete additivity, so for
     them nothing is proved either way.
2. **The classification of the leftover is wrong under multiplicativity** (F2, major). The "transfer of the round-10
   referee's argument" gives the max-scale cross-gap only for the **minimal** graft. With cross multiplicativity,
   primitive support contains every mixed product `d*beta`, and the sharp large sieve then demands
   `|sin(psi_d + arg beta - arg alpha')| >= c/max(Nd*Nbeta, N alpha')` for all `beta, alpha'`. This is a product-scale
   (LCG-type) gap, not the max-scale gap. Exact spacing computation: a badly approximable fake angle passes the
   max-scale form (minimal graft, normalised spacing `>= 0.23`), yet fails the multiplicative form (`5.6e-5` at
   `X = 10^6`; genuine `1.003`).
3. **The "needed statement" (§5) is mis-specified.** Uniformity in sector position is not necessary, and the
   sharp/smooth dichotomy is not exhaustive: F3 is neither sharp nor smooth, is position-dependent, and has
   precision `X^(o(1))`.
4. Minor (F1): Corollary 2.2 repeats the constant slip that the round-10 referee already corrected (Slip 1). The
   `lw.md` 6.1 fakes need `A > 4`, not `A >= 4`.

None of the killed items yields growth or a uniform bound. Each is a **classification** error in the barrier map:
the composite class is capped only in its robust part, and **at least three** non-robust composite items stay
unclassified, not one.

## 1. What the proposal claims

* **Two grafts** of a fake `(S, phi)` onto the ideals of `Z[i]`. The *minimal* graft changes only the S-smooth ideals.
  The *multiplicative* graft rotates each S-fibre `F_d = {d beta : beta coprime to S}` by `rho(d)`.
* **Prop. 1.5.** The optimal all-ideal large-sieve constant satisfies `Lambda* >= max(#I(X), (2K+1) floor(X^(1/2)))`, so
  every true all-ideal Hecke large sieve is robust.
* **Thm M.** Robust composite inputs (C1)-(C3) survive the minimal graft up to `1 + O(t^(1/2))` and `1 + eta`.
* **Thm P.** P1-P5 survive the multiplicative graft.
* **Cor. 2.2.** Grafts of the `lw.md` 6.1 fakes and of the P3-char fakes give the cap.
* **Prop. 3.1.** Smooth Hecke sums together with multiplicativity force `phi_p = +-theta_p (mod pi/2)` (collapse).
* **§3.3.** The only non-robust leftover is the sharp primitive large sieve, equivalent to the max-scale cross-gap.
* **§5.** A separating input needs precision `X^(gamma(c_0))`, `gamma(c_0) < 2/c_0`, uniformly in position. Sharp
  counts fail this; smooth sums collapse.

## 2. Re-derivation of Theorem M and the lemmas

**Lemma 1.2.**
* (1) There are `2s` prime ideals above `S`, each of norm `>= P_0`. With `t <= 1/4`, `-log(1-t) <= 2t`, so
  `log prod <= 4sP_0^(-sigma) = u <= 1`, and `e^u - 1 <= (e-1)u <= 2u`. Correct. The weighted form is correct for
  `t <= 1/2`.
* (2) The count `binom(2s+L, L)` is correct.
  * For `L <= 2s`: `binom(2s+L, L) <= (e(2s+L)/L)^L <= (4es)^L`.
  * For `L > 2s`: `binom(2s+L, L) <= 2^(2s+L) <= 4^L`.
  * `d(d) <= 2^Omega(d)` because `v + 1 <= 2^v`. Correct.

**Lemma 1.3.** Re-derived (convex sector, parallel bodies `K_(+-r)`, `r = 2^(-1/2)`). The sharp error is in fact
`rP + 1`; the stated `rP + pi/2` is conservative. `rP <= 0.7071 * 3.5708 y^(1/2) = 2.525 y^(1/2)`. For `y_2 < 1` the
annulus case is also fine, since the main term `<= (pi/4) y_2 <= 8.3 y_2^(1/2)`. Correct.

**Lemma 1.4 (Gallagher).** Steps checked:
* `int |h'| <= 2 ||T||_2 ||T'||_2 <= 4 pi K ||b||^2`.
* With `delta = 1/K`, an arc of length `1/K` in `t = 2 vartheta/pi` is an arc of length `pi/(2K)` in `vartheta`.
* So `Q <= pi Y/(4K) + 2.53 Y^(1/2) + 1.58`, and the total is
  `(1 + 4 pi)(pi Y/4 + 4.11 K Y^(1/2)) = 10.66 Y + 55.8 K Y^(1/2) <= 56(Y + K Y^(1/2))`.
* At `K = 1` the closed arc of length `pi/2` double-counts the edge ideals. This is harmless for an upper bound.

Correct. Independent numerical check: `c1_out.txt` (§7).

**Theorem M.**
* (1) `|D_k| <= 2 sum |a_alpha|` over the changed (S-smooth, `!= 1`) ideals, which is `<= 2 Psi_S(X)^(1/2) ||a||`.
  * Summed over `2K + 1 <= 3K` frequencies: `<= 12 K Psi_S ||a||^2`. Minkowski then gives `(1 + t_X^(1/2))^2`.
  * For `X >= P_0`: `Psi_S X^(-1/2) <= X^(kappa - 1/2) <= P_0^(kappa - 1/2) = 4es P_0^(-1/2)`, since `kappa < 1/2`.
  * So `t <= 48 e s P_0^(-1/2)/c = 130.5 s P_0^(-1/2)/c`. Correct.
* (2) The change is `<= 2 C_c 2^(BL) Psi_S(X_1) <= 2 C_c X_1^(kappa_B)`. Divided by `c X_1^gamma` it is
  `<= (2 C_c/c) P_0^(kappa_B - gamma)` once `kappa_B < gamma`. Correct.

Theorem M is a Hamming-distance argument. It is valid for **any** functional that is Lipschitz in the changed
ideals (Remark 2.1), but **only at a precision coarser than the footprint `Psi_S(X) <= X^((1+o(1))/A)`**. That is
exactly why it says nothing about the inputs of §5 (F3, F4).

## 3. Re-derivation of Theorem P and Prop. 3.1

**P1.**
* The duality step is correct: the operator and its adjoint have the same norm, for any support `T`.
* Fibre `d = 1`: the sum over `beta in T` coprime to `S` is a sub-sum of the genuine dual sum over `T`, so it is
  `<= Lambda ||b||^2`.
* Fibre `d != 1`: Lemma 1.4 at height `X/Nd`.
* Summing over `d != 1` gives `56(X tau_1 + K X^(1/2) tau_(1/2)) <= 56 tau_(1/2)(X + K X^(1/2))`.

Correct. This genuinely answers `round10/referee_analytic.md` Slip 3.

**P2, P2'.** `lambda^x_k = lambda_k u_k` with `u_k` completely multiplicative and unimodular, so a free-coefficient
class closed under unimodular twists absorbs the twist. This is trivial and correct. The `D^2` change is
`<= sum 2/N frak p = 4s/P_0`.

**P3.**
* `alpha = d beta`, with `d` S-smooth and coprime to `q`, and `beta` coprime to `qS`.
* Inclusion-exclusion over the squarefree `e | rad(qS)`, with `vartheta(e)` shifts.
* The main terms depend on the arc only through `|I'|`.
* The error is `16.6 (x/Nd)^(1/2) sum_e Ne^(-1/2)`, where `sum_e Ne^(-1/2) = Pi_q (1 + tau'_S)`.

Correct. *Scope remark:* `Pi_q` can be as large as `q^(o(1))`. So P3 is uniform in a **fixed** `q`. Congruence
families with `q` up to `x^theta` (Bombieri-Vinogradov type) carry an extra factor `x^(o(1))`, which can exceed
`M^(A/2 - 1)` at large `x`. This is minor, and the proposal states P3 for a fixed `q` only.

**P4.**
* Hyperbola split at `y^(1/2)`: the three pieces `A + B - C` each have a rotation-invariant main term.
* The error is `<= 3 * 8.3 * 8.4 (1 + tau'_S) y^(3/4) = 209.2 (1 + tau'_S) y^(3/4)`.
* `d(d beta) = d(d) d(beta)` for coprime arguments.

Correct.

**P5.**
* Box version: `F_d(I')` expands via `w 1_coprime = w * h`, and both rotations share the main term. Correct.
* Twisted version: the fibre `d = 1` is the same on both sides. The difference is
  `<= 2 C_E X_1^gamma tau^(w)_gamma (1 + tau^(h)_gamma)`. Correct.

**Prop. 3.1 (entire-L collapse).**
* Mellin identity: `int G(X) X^(-s-1) dX = f~(s) L'(s)`. The `O(X^(-delta_f))` decay makes `f~ L'` analytic on
  `Re s > -delta_f`.
* Poles of `E_p` at `s_m = i(arg b_p + 2 pi m)/log p`. A numerator of another `E_(p')` cancels at most one `m` per
  `(p', sign)`, by unique factorisation (`p'^(m_1 - m_2) = p^(m'_1 - m'_2)`).
* `L(it, lambda_k) != 0`: `Gamma(it + 2|k|)` is finite and nonzero for `k != 0`, and `L(1 - it, conj lambda_k) != 0`.
* Test function `f = f_0 u^(-i t_m)` gives `f~(s_m) = int f_0/u > 0`.
* Conclusion `{b_p, conj b_p} = {a_p, conj a_p}`. At `k = +-1` this gives an actual configuration of lattice points of
  norm `N` with the same profile up to flips, on the same arc. So the fake is an actual cluster.

Correct. **Hypothesis to keep in view:** `vartheta'` is completely additive and agrees with `vartheta` off `S`. Without
complete additivity the proof does not start (§5, F4).

## 4. Corollary 2.2 and Prop. 1.5

**Prop. 1.5.**
* The lower bound `(2K+1) floor(X^(1/2))` comes from the rational integers: `lambda_k(n) = 1`, `N(n) = n^2`.
* `#I(X) >= X/2` for `X >= 100`.
* "Max >= mean" gives `Lambda* >= (X + K X^(1/2))/4`.

Correct, and confirmed by exact Gram eigenvalues (§7, `c1`). The statement "every true all-ideal large sieve is
automatically robust" holds for supports containing `I(X)` or the integers up to `X^(1/2)`. As the proposal itself
says, it does not hold for primitive or other "one ideal per angle" supports.

**Corollary 2.2: constant slip (F1, minor).**
* The proposal uses "`A > max(2, 1/gamma)`, `A >= 4`" for the `lw.md` 6.1 fakes. This repeats
  `round10/analytic.md` Thm B verbatim, including the slip that `round10/referee_analytic.md` Slip 1 corrected:
  `lw.md` Thm 6.1 needs `P >= (40M/C_0)^4 > M^4`, so `P = M^A` is admissible only for `A > 4`.
* Corrected:
  * `lw.md` 6.1 grafts: `A > max(4, 1/gamma)`, `W <= (2A + o(1)) M log M`;
  * P3-char grafts: `A > max(2, 1/gamma)`, `W <= (bA + o(1)) M log M`, `b >= 33`.
* The rate is unchanged.
* Apart from this, the asymptotics `kappa, kappa_B = (1 + o(1))/A`, `t = O(M^(1 - A/2))` and
  `eta = O(2^B M^(1 - A gamma))` are correct for `s < 2M` (resp. `b(M - 1)`), `P_0 >= M^A`.

**Barrier logic.** It is valid **for inputs with `(c, gamma, B, C_c)` fixed before `c_0`**, exactly as Thm B. The
proposal's Remark 2 notes that, unlike the prime case, **true** composite inputs exist at every precision. So for the
composite class the R1 caveat is not a technicality: the inputs the cap leaves out exist and are true. That is the
point of §5.

## 5. The leftover and the needed statement: three classification errors

### 5.1 F2 (major): with multiplicativity the sharp primitive large sieve is NOT the max-scale cross-gap

**The equivalence it rests on.** The round-10 referee's equivalence is: sharp large sieve at height `X` iff all
points of the support are `c/X`-separated. One direction is the two-point argument; the other is Montgomery-Vaughan.
I accept this equivalence; it is a statement about point sets.

**Which point set.** The point set is the angle multiset of the **support** under the grafted angle function. On
primitive support:

* **Minimal graft.** Only the S-smooth primitive `d` move (to `psi_d = <v, phi>`). The conditions are fake-fake
  spacing, which follows from Lemma 5.1 since `e^(-w(v - v')/2) >= 1/max(Nd, Nd')`, and the fake-genuine condition
  `|sin(psi_d - arg alpha)| >= c/max(Nd, N alpha)`. This is the proposal's composite max-scale cross-gap. **Correct
  for this graft.**
* **Multiplicative graft.** Every primitive `alpha = d beta` with `d != 1` S-smooth moves, to `psi_d + vartheta(beta)`.
  The fake-genuine spacing now reads
  `(**)  |sin(psi_d + arg beta - arg alpha')| >= c / max(Nd * Nbeta, N alpha')`
  for all `beta` coprime to `S` with `d beta` primitive, and all primitive `alpha'` coprime to `S`.
  * Put `gamma = alpha' conj(beta)`. Then (**) says `|sin(psi_d - arg gamma)| >= c/max(Nd Nbeta, N gamma/Nbeta)` for
    every factorisation of `gamma`.
  * For `gamma` with a divisor of norm about `(N gamma/Nd)^(1/2)`, this is a **product-scale** gap
    `>= c/(Nd N gamma)^(1/2)`. That is the LCG shape of `round10/analytic.md` Thm D, not the max-scale shape.

**Exact check** (`c2_sharp_primitive.py`, `c2_out.txt`).
* Setup: one fake prime `p = 101`, fake angle `psi` with `tan psi = sqrt7 - 2`. This is a quadratic irrational, so
  every `j psi` is badly approximable, and the max-scale gap is plausible.
* Quantity: `X * min spacing` over **all primitive ideals** of norm `<= X`.

| `X` | genuine | minimal graft | multiplicative graft |
|---|---|---|---|
| `2*10^3` | 1.041 | 0.833 | `1.16e-1` |
| `5*10^4` | 1.010 | 0.513 | `1.14e-2` |
| `10^6` | 1.003 | 0.228 | `5.6e-5` |
| `4*10^6` | 1.001 | 0.456 | `2.3e-4` |

* The genuine value is always `>= 1`, the lattice fact.
* The minimal graft stays bounded below. This is consistent with the max-scale (badly approximable) condition.
* The multiplicative graft decays roughly like `p/X`, as expected for about `X^2/p` "random" cross pairs.
* So the same fake angle **passes** the max-scale form and **fails** the multiplicative sharp large sieve by four
  orders of magnitude.

**Status of the multiplicative form (heuristic).** Let `gamma` run over the lattice approximants (convergents) of
the direction `psi`. Taking `beta = 1` in (**) forces bounded partial quotients. Taking `beta` a divisor of `gamma`
forces the S-coprime part of every large convergent to have no divisor with norm in `(1/(ca), ca N gamma/Nd)`
(`a` the next partial quotient). So every large convergent must be "small times prime". With bounded partial
quotients, the branching tree of continued fractions has about `K` children per level, against a prime probability of
about `1/(n log K)` at depth `n`, so it dies out. Hence heuristically **only lattice directions pass**: this form
behaves like LCG/Thm D (collapse), not like an unclassified badly-approximable condition.

**Consequence.** The proposal's sentence "the only non-robust leftover is the sharp large sieve on primitive support
... equivalent to a composite max-scale cross-gap" is wrong whenever the argument also uses multiplicativity. The
correct split is:
* (a) without multiplicativity: max-scale cross-gap, open, unclassified (inherits round 10);
* (b) with multiplicativity: product-scale gap (**), heuristically collapsing, not proved.

Neither gives growth.

### 5.2 F3 (major): a non-robust, position-dependent composite input with precision `X^(o(1))`

**The input.** Let `f` be a fixed `C^infinity` bump on `(1/2, 2)` and `I = [psi_1, psi_2]` an arc whose edge slopes
`tan psi_i` are quadratic irrationals (badly approximable). Put
`S_I(X) = sum_alpha f(N alpha/X) 1[vartheta(alpha) in I]`.
* The radial weight is smooth, so the circle-problem error disappears.
* The angular cutoff is sharp. By 2D Poisson summation, the remaining error comes from dual lattice points near the
  edge normals. It is controlled by the continued fraction of the slopes, so it is `O(1)` to `O(log X)` for badly
  approximable edges (Ostrowski-type; **heuristic as a theorem here, I did not write the proof**).
* For rational edges the error is `asymp X^(1/2)` (lattice points on the edge).

**Evidence** (`c3_ba_sector.py`, `c3_out.txt`). Genuine error over `10^4 <= X <= 10^7`:
* BA edges `(sqrt2-1, (sqrt5-1)/2)`: `|err| <= 0.131`.
* BA edges `(sqrt7-2, sqrt3-1)`: `|err| <= 0.118`.
* Rational edges `(1/3, 1/2)`: `+2.25, +7.11, +22.5, +45.0, +71.1`, i.e. `0.0225 X^(1/2)`.

**Fake-violated** (`c3b_ba_graft.py`, `c3b_out.txt`; evidence).
* Setup: 30 random sectors with quadratic-irrational edges; a fake with `s = 40` primes in `[1009, 1601]` and
  Haar-random `phi`; `X = 1.5*10^6`, about `P_0^2`.
* Genuine `|err|`: max 0.448, mean 0.141.
* Minimal graft `|err|`: max 8.06, mean 2.30, i.e. 16 times larger on average. (At `X = 2.5*10^6` the window weight
  of the S-smooth ideals is smaller, and the ratio is 4.)

**Why neither theorem covers it.**
* Theorem M needs the input's error to be `>= c X^gamma` with `gamma > kappa_B ~ 1/A`. Here the error is
  `O(log X)`.
* Theorems P3 and P5 need the genuine statement **uniformly in position**, because the fibre `d` sees the sector
  rotated by `rho(d)`, which destroys the badly approximable edges.
* The input holds only for badly approximable positions, with position-dependent constants. The proposal's open
  point 1 notes the gap only "with cross multiplicativity". F3 is uncovered even **without** it.
* Whether F3 collapses fakes (with or without multiplicativity) is **not decided** here. The proposal does not list it.

### 5.3 F4 (major): smooth sums without multiplicativity are not "irrelevant"

* Prop. 3.1 assumes `vartheta'` is completely additive. An argument that uses smooth all-ideal Hecke sums **without**
  using multiplicativity only needs **some** extension of the fake that satisfies the sums.
* Both standard grafts violate them (the proposal's own `smooth_check`: `0.03` to `3.2` against `1e-12`). But other,
  compensating extensions are not excluded.
  * For the asymptotic form ("for each `f`, `O_f(X^(-delta_f))`"), a graft that changes only the finitely many ideals
    the argument uses (the fake's points and divisors) is invisible. So the asymptotic form alone **cannot** separate:
    it is trivially fake-satisfied.
  * For effective forms (explicit `C(f, k, B)` at height `X ~ N`), the question is whether the about
    `prod(e_p + 1) - M` free ideals of norm `N` can compensate the `K` power sums that the precision makes visible.
    This is a counting question that neither the proposal nor I decide.
* The proposal lists this as open point 3 and calls it "irrelevant for growth". That judgment is unsupported. It is
  an unclassified non-robust composite input class in its own right.

### 5.4 The "needed statement" (§5 of the proposal) is mis-specified

The proposal states:

> a separating composite input needs error `<= X^(gamma(c_0))`, `gamma(c_0) < 2/c_0`, **uniformly in sector position**
> at heights `>= M^(c_0/2)`; sharp counts fail it, smooth sums collapse.

Three problems:
* **Uniformity in position is not necessary.** The fake's information sits on the S-smooth ideals in **every**
  direction, so an argument may use sectors at positions it chooses (for example badly approximable edges), not at
  the cluster's arc. F3 meets the precision requirement non-uniformly.
* **The dichotomy "sharp counts / smooth sums" is not exhaustive.** F3 is smooth in norm and sharp in angle.
* **"Smooth sums collapse" holds only with multiplicativity** (F4).

The correct reading is that the cap holds for inputs whose precision is fixed before `c_0` **and** that are either
fibre-blind (Thm M) or used through P1-P5. Everything finer is unclassified. None of it is known to give growth.

## 6. Barrier checklist, item by item

The proposal is a cap, so the checklist is applied to (i) its cap, (ii) its positive candidates, and (iii) the extra
candidates F2(b), F3 and F4 found here.

1. *Fake-satisfied inputs (pairs, gcds, cut, residues, characters, GAP-Fourier, energy, `L^p`).*
   * (i) The cap adds robust composite ensemble inputs to this list. Valid.
   * (iii) F2(b), F3 and F4 are not of these kinds: they are pointwise or precision-`X^(o(1))` ensemble statements.
     F3 is violated by the minimal graft (evidence). F2(b) is violated by the multiplicative graft (exact spacing).
2. *Robust or averaged statements (Thm B, Prop C, Thm D).*
   * (i) Thm M and Thm P extend Thm B correctly.
   * F2(b) is a pointwise statement, heuristically of LCG/Thm-D type, so a collapse.
   * F3 is not averaged over position. It is a fixed sector with badly approximable edges, so Prop C (which concerns
     certifying emptiness by absolute-value Fourier counting) does not apply to it as stated.
   * F4 is pointwise in the smooth weight.
3. *Auxiliary polynomials on `X_M` / `Mbar_(0,k)`.* Not used.
4. *Single linear form; SLW.* Not used by any of the candidates.
5. *S-integer thresholds through `|S|`.*
   * The cap's bounds depend on `S` only through `(s, P_0)`. They are footprints, not thresholds. Fine.
   * Any use of F3 at a cluster would have constants depending on the edge continued fractions. Those are independent
     of `S`, so no `|S|`-threshold arises.
6. *Route 5 box statement / CIRCLE(q).* Not used.
7. *`(L_5)`/`(L_6)` mechanisms at `k = 8`; computational `(L_8)`.* Not used.
8. *abc unit-pencil lane.* Not used.
9. *Exact Dirichlet exponent (`round4/hardness.md`).*
   * F2(b) sits at the exact exponent: it is a product-scale gap with constant `c`, like LCG.
   * F3 sits at precision `O(log X)`, an absolute-constant-type statement.
   * Neither is reachable by Roth/Schmidt-type tools in the form needed. Neither is shown to give growth.

**Project constraints.**
* No fixed-S Subspace step appears.
* No gcd or denominator cost is discarded: the grafts keep the profile.
* No drifting phase is set to zero: `rho(d)` is kept.
* The cap's threshold `M_4` depends only on the input's constants, as in Thm B.
* The proposal's §5 "needed statement" does **implicitly** treat position-uniformity as necessary. That is not a
  violation of the constraints, but it is an unjustified restriction (§5.4).

## 7. Computations

Every script is in `kill_composite_checks/`, run under `prlimit --as=3500000000`, with output streamed to `*_out.txt`.

| script | label | result |
|---|---|---|
| `c1_ls_constant.py` | evidence (eigenvalues, double precision) | `Lambda*(X, K)` as the top eigenvalue of the Toeplitz Gram matrix, `X = 10^2 ... 10^5`, `K = 0, 1, 5, 20, 100`. Prop. 1.5 lower bound and Lemma 1.4 upper bound hold in all 20 cases. `Lambda*/(X + K X^(1/2))` lies in `[0.70, 1.83]`. |
| `c2_sharp_primitive.py` | exact integer data, floating-point angles | Normalised min spacing on primitive support, `p = 101`, `tan psi = sqrt7 - 2`: genuine `>= 1.001`; minimal graft `0.23-0.83`; multiplicative graft `1.2e-1` falling to `5.6e-5` (F2). |
| `c3_ba_sector.py` | evidence | Smoothed sector counts. Badly approximable edges: `|err| <= 0.131` up to `X = 10^7`. Rational edges: `0.0225 X^(1/2)` (F3). |
| `c3b_ba_graft.py` | evidence | 30 badly approximable sectors, `s = 40` fake primes. Genuine mean `|err|` 0.141 against minimal-graft mean 2.30 at `X = 1.5*10^6` (F3 is fake-violated). |

The proposal's own outputs (`composite_checks/*_out.txt`) are consistent with my re-derivations. In particular,
`Psi_S = binom(2s+L, L)` exactly, and the graft changes to the large-sieve constant are small.

## 8. What survives (exact statement)

> **Surviving theorem (cap for robust composite inputs).** Fix `c > 0`, `gamma in (0, 1]`, `B >= 0`, `C_c >= 1`
> before the target constant.
> * Let `A > max(4, 1/gamma)` for the `round8/lw.md` 6.1 fakes, or `A > max(2, 1/gamma)` for the P3-char fakes.
> * Then for all large Paley orders `M`, the minimal graft satisfies every `(c, gamma)`-robust composite input
>   (C1)-(C3), up to `1 + O(M^(1 - A/2) + 2^B M^(1 - A gamma))`.
> * The multiplicative graft satisfies P1 (any support with `Lambda >= c(X + K X^(1/2))`), P2, P2', P3
>   (`gamma >= 1/2`, fixed `q`), P4 (`gamma >= 3/4`) and P5 (position-uniform families).
> * Hence such inputs, combined with the profile, the monomial gaps (or P3-char residues at moduli `<= 3M^A`) and
>   Haar-generic properties, cannot prove `M = o(log R/loglog R)`.
> * The R1 caveat (`gamma`, `A` chosen after `c_0`) applies, and here it is not vacuous: true composite inputs exist
>   at every precision.
>
> Also surviving: Prop. 1.5 (for supports containing `I(X)` or `{n <= X^(1/2)}`) and Prop. 3.1 (smooth Hecke-sum decay
> plus complete additivity forces the fake to be an actual cluster).

**Not surviving.**
* "Composite support adds no separating power."
* "The only non-robust leftover is the sharp primitive large sieve, equivalent to the max-scale cross-gap."
* The §5 "exact strength requirement" with position-uniformity.

**Corrected barrier-map entry.**
* §5.4 gets a new row: robust composite inputs, kind F, the witness above, with the R1 caveat.
* §5.5 keeps "composite-supported analytic inputs, non-robust part", with at least these unclassified items:
  * (a) sharp primitive large sieve without multiplicativity, i.e. the max-scale cross-gap;
  * (b) the same with multiplicativity, i.e. the product-scale gap (**), heuristically collapsing;
  * (c) smoothed sector counts with badly approximable edges (precision `O(log X)`, position-dependent);
  * (d) effective smooth Hecke sums used without complete additivity.
* None of them is known to give growth or a uniform bound.
