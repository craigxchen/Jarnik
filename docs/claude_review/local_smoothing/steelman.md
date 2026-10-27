# Steelman: the best case that the local smoothing / Kakeya / restriction news helps us

Role: argue as hard as honesty allows that the claimed theorems (openai/math #079, "Local smoothing for the Euclidean
wave equation", dated 2026-10-25), or the auxiliary-PDE method behind them, could give a growth improvement over
`log R/loglog R` or the uniform bound, or could decisively unblock route 5 (MIN^pair), P3-all\* or H^gen/L7.
Written 2026-10-27. It builds on `facts.md`, `chains.md` and `lifts.md` in this directory and does not repeat them.

Labels:
* **[V]** read in the primary source (the cloned repository, `dl/math`);
* **[P]** proved here, elementary, with all quantifiers;
* **[C]** checked numerically here; the script is in `steel_checks/`;
* **[L]** literature located this session, from search snippets only (arXiv, EMS and Cambridge hosts were unreachable);
* **[M]** from memory, not re-checked;
* **[S]** speculation.

## 0. Verdict

**Not useful, even steelmanned.** I found no mechanism that plausibly gives `M = o(log R/loglog R)` or unblocks an
open route. Every candidate hits at least one of four walls:
* **W1, single shell.** On a single circle the real-time inequality is vacuous (lifts F1).
* **W2, resolution.** L^p sees an arc cluster only for `p > log R/log M`. All L^p and energy chains are capped at
  `Theta(log R)` (chains Prop 2, lifts F2).
* **W3, Euclidean.** Fakes and the transcendental configuration `Z` of chains Prop 1 satisfy every Euclidean input.
* **W4, regime.** The new results are fixed-dimension, continuous upper bounds. Route 5 needs a discrete statement
  in growing dimension, at exponential precision, that sees arithmetic.

What the steelman adds beyond the three assessments, in order of strength:

1. **S3 (the best case, but it fails).** Route 5's residual inequality is an *anti-concentration* statement, the
   same kind of statement the new method proves.
   * [P] MIN^pair (in its H = infinity form) follows if no degree-`m'` polynomial with leading form `<k',u>^(m')`, for
     `u` in the gap region, can be `d^(m') eps_0`-small on half of the ternary cube.
   * [P] For a Gaussian (continuous) proxy of the cube, Carbery-Wright proves this at once, with room
     `exp(-Omega(W/m'))`, for **every** `u != 0`.
   * So the whole content of MIN^pair is discrete. The known discrete-to-continuous transfers (invariance principle,
     polynomial Littlewood-Offord) are vacuous at degree `m' ~ 4 sqrt(8 d ln 4r) >> log d`.
   * The paper supplies only continuous anti-concentration, so it does not touch that gap. By fake Thm 2.1, closing
     the gap is the barrier question itself.
2. **S4 (a real conceptual match, but no mechanism).** The paper's two tools are a Monge-Ampere potential and an
   Aubin-Mabuchi energy. In pluripotential theory these are the classical tools for *both* halves of the determinant
   method:
   * the homogeneous Bernstein-Walsh inequality (Robin function) is a continuous MIN;
   * Berman-Boucksom show that the Aubin-Mabuchi energy governs transfinite diameters, i.e. Vandermonde sizes.

   The transfer to the finite set `A` breaks at degree 3: `k_x^3 - k_x` vanishes on `A`, so `A` is never norming.
   The "pluripotential theory of `A` at degree `m'`" is exactly the record's leading-form ideal `J(A)`.
3. **S1 (a correction to the scope of lifts F1).** F1 kills the real-time inequality on a single circle, but not the
   proof's auxiliary objects.
   * [P][C] The imaginary light-ray translation, restricted to `|xi| = R`, is a Gaussian arc filter of angular
     width `(l R)^(-1/2)`. At displacement `l ~ 1` this is exactly our arc scale `C/sqrt R`.
   * Its Bessel expansion is a sum of Hecke-character values. I checked the identity to `10^(-9)` on four norms.
   * But the paper only carries *lower* bounds through that scale. Its quantities are L^4 in the plane (W2), and
     the Fourier side is the capped "character" class.
4. **S2 (a precedent, with a reason it does not transfer).** Euclidean incidence and decoupling methods have already
   beaten divisor-bound losses `exp(c log N/loglog N)`, which have our `log/loglog` shape:
   * Herr-Kwak 2024 [L]: the sharp `L^4` Strichartz estimate on `T^2`, via Szemeredi-Trotter;
   * Guth-Maldague-Wang 2024 [L]: `(log R)^c` decoupling, with an appendix on sixth-order correlations of
     `x^2 + y^2 = m`.

   In both cases the divisor bound sat inside an *average*. Our `M` is the *supremum* itself, so W2 and W3 apply.
5. S5-S10 (Kakeya/incidence for clusters, multi-shell lifts, square functions for route 5's forms, P3-all\* plates,
   H^gen, formal infrastructure): each fails, for the explicit reasons given in §2.

**What would change the verdict:** a *discrete* version of the paper's positivity method: an anti-concentration or
inverse-Littlewood-Offord theorem for polynomials of degree about `sqrt(d log r)` on the ternary cube cut by slabs,
whose exceptional directions are collision strata. Nothing of the kind is claimed or in sight (§3).

## 1. The walls, with the exact numbers used below

* **W1** (lifts F1, [P]). If `f` has frequencies on one sphere `|xi| = rho`, then `e^(-iT|D|) f = e^(-iT rho) f` and
  `J^s f = (1 + rho^2)^(s/2) f`. Every local smoothing, Bochner-Riesz and square-function inequality is then a scalar
  identity.
* **W2** (chains Prop 2, lifts F2, [P]). Let `F = sum over cluster e(<z,x>)` on `T^2`. Then:
  * `|F| >= M/2` on a rectangle of area `>= c C^(-3) R^(-1/2)`, so `||F||_p^p >= c (M/2)^p C^(-3) R^(-1/2)`;
  * any `||F||_p <= K_p ||F||_2` therefore gives only `M <~ R^(1/p)`;
  * even the strongest Lambda(p) or energy hypotheses give at best `M = O(log R)`.

  The obstruction is the volume factor `R^(-1/2)`, not a loss. **An eps-free, log-free estimate changes nothing.**
* **W3** (chains Prop 1, record fake theorems). Two kinds of configuration satisfy every perturbation-stable
  Euclidean inequality:
  * `Z = {R e^(ik/R)}`, which has `floor(C sqrt R)` points;
  * the record's fakes, which have genuine primes, exponent profile and gap inequalities, but real angles.
* **W4** (lifts F3/F4). The new results are L^p upper bounds in fixed dimension `n`, with curvature and with
  `C_(n,p,eps)` of unstated growth in `n`. Route 5 lives in dimension `d = n - 1 -> infinity`, on flat ternary
  frequency sets, at precision `e^(-W/4)`.

A useful input must get past at least one wall without running into another.

## 2. Candidates

### S1. The light-ray filter is non-vacuous on a single circle, at exactly our scale

**Observation [P].** The paper's mechanism [V] (introduction; `current.tex`) translates spacetime by
`-i l (1, omega)` with `|omega| = 1`. This multiplies the mode `e^(i(xi.x - T|xi|))` by
`exp(-l(|xi| - xi.omega))`. On the circle `|xi| = R`, with `theta_xi = arg xi`, this is

```text
w_l(xi) = exp(-l R (1 - cos(theta_xi - omega)))  ~  exp(-l R (theta_xi - omega)^2 / 2),
```

a Gaussian *arc filter* of angular width `(l R)^(-1/2)`. In the paper the auxiliary radius `s = |tau|` runs from
`2a = 2Kh` to `3`, with `h = 1/N`. I read `s` as playing the role of `l`: the barrier potential `C s^-2 e^(-cs/h)`
lives on the scale `s ~ h`, which fits. I did not verify the exact identification in the current construction.

On that reading:
* `|tau| ~ 2a` is a fixed-aperture filter, of width `K^(-1/2)`;
* `|tau| ~ 1` is a filter of width `N^(-1/2) = R^(-1/2)`. That is the canonical-sector scale, i.e. an arc of length
  `~ sqrt R`.

So, unlike the real-time inequality, the proof's auxiliary objects are not vacuous on one circle. Evaluated at
`x = 0` with unit coefficients, the filtered wave is the smoothed arc count

```text
S_omega(l) = sum_{|xi|^2 = N} exp(-l R (1 - cos(theta_xi - omega))),
S_omega(l) >= e^(-l C^2/8) * #(points within angle C/(2 sqrt R) of omega),
```

and `max_omega S_omega(1)` is comparable to `M` up to constants that depend on `C`.

**Fourier side [P].** `e^(kappa cos phi) = sum_k I_k(kappa) e^(ik phi)`. So

```text
S_omega(l) = e^(-lR) sum_k I_k(lR) e^(-ik omega) sum_xi e^(ik theta_xi).
```

Only `k = 0 mod 4` survives, and the inner sums are Hecke Grossencharacter values at `N`. Since
`I_k(kappa) ~ e^kappa e^(-k^2/2kappa)/sqrt(2 pi kappa)`, frequencies `|k| <~ sqrt(lR)` contribute.

**Numerical check [C].** `steel_checks/arcfilter.py` uses `l = 1` and `C = 1`:

| `N` | points | `M` | `max S` | Bessel series |
| --- | --- | --- | --- | --- |
| 243061325 | 384 | 2 | 2.97764164 | 2.97764165 |
| 157163452745 | 1024 | 3 | 2.96771074 | 2.96771074 |
| 303460625 | 200 | 2 | 2.15906618 | 2.15906617 |
| 1349232625 | 256 | 2 | 1.97349735 | 1.97349735 |

**Needed.** An *upper* bound on the arc-scale concentration `max_omega S_omega(1)`, or on the arc-scale energy
`F(T, |tau| ~ 1)`, of size `o(log R/loglog R)` or `O(1)`.

**Available [V].** The paper's Proposition "fixed-barrier" (`barrier.tex`) gives three inequalities:
* `F(0, |tau| = 4a) <= C a^2 ||f||_(p_c)^(p_c)`;
* `F(T, |tau| = 2a) >= c a^2 ||e^(-iT|D|) f||_(p_c)^(p_c)`;
* `Delta F <= C s^-2 e^(-cs/h) F` on `2a <= |tau| <= 3`.

The two circle comparisons carry **lower** bounds on `F` from `|tau| = 2a` out to `|tau| ~ 1` and back. Each
comparison loses a log. The arc-scale circle is a waypoint that receives lower bounds. No upper bound for the energy
there is asserted.

**Gap.**
1. *Direction.* The paper concludes that concentration at coarse aperture forces energy at fine aperture. A
   cluster *makes* the fine-aperture energy large, which is consistent with every inequality in the paper.
2. *W2.* `F` detects `int mu^(n/(n-1))`, which is an L^4 quantity in the plane, since `p_c = 4`. Take the
   whole-circle function `f = sum_{|xi|^2=N} e(xi.x)`. Zygmund gives `||f||_4^4 <= 3 r_2(N)^2`. The arc-filtered
   piece has `||.||_4^4 >= c M^4 C^-3 R^(-1/2)`. So the best conceivable comparison gives only
   `M <~ C^(3/4) R^(1/8) r_2(N)^(1/2)`.
3. *Capped class.* The Fourier side is a Gaussian-weighted Hecke-character sum. That is the record's "character /
   Fourier" class, which is capped.

**Status: dead.** It is still worth recording: it is the one place where the new proof touches our exact scale,
and it shows that the touch is through L^4 lower-bound bookkeeping.

### S2. Log-loss and eps-free forms, extremisers, and the divisor-bound precedent

**Needed.** Either:
* single-circle L^p bounds at `p ~ log R/loglog R` with constants `e^(o(p))` (W2); or
* an inverse theorem: near-extremisers of an L^p inequality for single-circle sums must be arc clusters, with
  quantitative stability.

**Available [V].**
* At `p = p_c`, a loss `(log N)^(2/p_c)`.
* A lower bound `(log N)^(1/p_c)` from *radial focusing*. Radial focusing uses many shells and all directions, so it
  says nothing about one circle (W1).
* Above `p_c`, log-Sobolev losses. Nothing for `p` growing with `N`.

**Gap.** W2 is a volume obstruction. Even a loss-free estimate gives `M <~ R^(1/p)` at fixed `p`. No inverse theorem
is claimed.

**The precedent a steelman must address [L].**
* *Bourgain's bounds.* Bourgain's `L^4` Strichartz bound on `T^2` and his `L^6` bound on `T` had losses
  `exp(c log N/loglog N)`. These came from divisor bounds, which have the same `log/loglog` shape as our bound.
* *Herr-Kwak* (Forum Math. Pi 2024, arXiv 2309.14275) proved the sharp `L^4` estimate with `(log N)^(1/4)` loss.
  They used Szemeredi-Trotter, which is Euclidean incidence geometry.
* *Guth-Maldague-Wang* (JEMS 26 (2024) 875-917, arXiv 2009.07953) proved `(l^2, L^6)` decoupling for the parabola
  with `(log R)^c` constant. Their appendix applies it to the sixth-order correlation of `x^2 + y^2 = m`.

So Euclidean methods *have* replaced `log/loglog`-type arithmetic losses by log powers.

**Why it does not transfer [P].** In both cases the divisor bound was a crude bound for one circle *inside an
average* over circles or boxes. The Euclidean method improves the average and never the supremum.

Our `M` is the supremum over all `R` and arcs. W2 shows that any route through averaged L^p or energy quantities is
capped at `Theta(log R)`. W3 shows that perturbation-stable inputs are satisfied by `Z` and by fakes.

The precedent confirms the direction found in `chains.md`: arc counts are an input to these estimates, never an
output.

**Status: dead.**

### S3. Route 5's MIN^pair is an anti-concentration statement (strongest candidate)

**Setting** (record, `round12/README.md` §1, §2.3). Use the notation `d`, `A`, `m'`, `Phi_c`, `kappa`, `eps_0`, `G`
of the record:
* `A` is a slab-cut subset of the ternary box `T = {k : |k_x| <= 1}`, with `3^d` points and `|A| >= 3^d/2`;
* `m' = ceil(4w)`, `w = sqrt(8 d ln 4r)`;
* `eps_0 = C d e^(-W/4)/(m'+1)`;
* `G` is the gap region `|u_x - u_y| >= e^(-eps_xy/2)`.

**Reformulation [P].** Use the nearstrata duality: `d^(m') kappa(u)` is the best uniform approximation error of
`<k',u>^(m')` on `A` by polynomials of degree `< m'`, where `kappa` is taken over all of `V(A, m')`. Then:

```text
kappa(u) > eps_0   <=>   every polynomial Q of degree m' with leading form <k',u>^(m') has max_A |Q| > d^(m') eps_0.
```

Since `|A| >= |T|/2`, this follows from the anti-concentration statement

```text
AC(u):  Pr_{k uniform on T} ( |Q(k)| <= d^(m') eps_0 ) < 1/2   for every such Q.
```

* *Caveat.* MIN^pair uses `kappa_H`, which ranges over the height-`<= H` span `V_H`. `AC` controls `kappa`, i.e. the
  `H = infinity` form. Whether `V_H = V(A, m')` is a separate question.
* *Why this answers part of W4.* `AC(u)` is an upper bound on a sublevel-set measure in the variable `k`, holding
  **for each** `u in G`. So it is not a measure-in-`u` statement, and the `u`-measure caps (round 10 Cor. B.1) do not
  apply to it directly.

**The continuous proxy is easy [P, using Carbery-Wright [M]].** Replace the uniform measure on `T` by a standard
Gaussian on `R^d`:
1. Every lower-degree `P` is orthogonal to the top Hermite chaos. So `Var Q >= ||(l_u^(m'))^(=m')||^2 = m'! |u|_2^(2m')`,
   which is `>= m'!` because `|u|_inf = 1`.
2. Carbery-Wright (2001) gives
   `Pr(|Q| <= eps sigma) <= C m' eps^(1/m')`.
3. With `eps sigma = d^(m') eps_0`, we have `eps <= d^(m') eps_0` and `ln(1/eps_0) = W/4 - O(log d)`. So the
   probability is at most `C m' d e^(-W/(4m') + O(log d/m'))`. In the regime, `W >~ 3 n ln n` and
   `m' = O(sqrt(n log n))`, so this is `exp(-Omega(sqrt(n log n)))`.

So in the continuous proxy, `AC(u)` holds with enormous room for **every** `u != 0`: there is no tangent cone at all.

**Hence the entire content of MIN^pair is discrete [P].** On the actual box it fails without the gap condition. For
a column direction `beta_p` the record exhibits `Q = <k', beta_p>^(m'-h-1) prod_t (<k', beta_p> - t)`, which vanishes
on `A` at degree `m'` (cheap directions, `h = h_A(beta_p) <= 2j < m'`).

So what is needed is an **inverse Littlewood-Offord theorem of degree `m' ~ 4 sqrt(8 d ln 4r)`**, at scale
`e^(-W/4)`, on the ternary cube cut by slabs. Its conclusion must be: if `Q` concentrates at `0` on half the cube,
then `u` lies near a collision stratum, with per-pair radii.

**What is available.**
* *Invariance principle* (Mossel-O'Donnell-Oleszkiewicz) [M]. The transfer error is
  `O(m' tau^(1/(8m')))`, where `tau >= 1/d` is the maximal influence. Since `m' >> log d`, this error is `>= 1`, so
  the principle is vacuous.
* *Polynomial Littlewood-Offord* (Meka-Nguyen-Vu, Theory of Computing 12 (2016) [L]). The bound is roughly
  `(log r)^(O(m log m))/sqrt(r_rank)`, and the rank is at most `d/m'`. Vacuous at `m' ~ sqrt(d log r)`.
* *Inverse Littlewood-Offord.* Known for linear forms (Tao-Vu, Nguyen-Vu) and quadratic forms [M]. Nothing at growing
  degree.
* *The new paper.* Continuous anti-concentration: density concentration forces a deep Monge-Ampere potential
  (Alexandrov), measured by an energy. It works in fixed dimension, for wave densities. Nothing discrete.

**Gap.** The discrete transfer at degree `~ sqrt(d log r)` is the whole problem. It is exactly the record's algebra of
`J(A)` and its fat-point base locus (rounds 12-13).

**The fake test applies [P from record].** fake Thm 2.1: a nondegenerate real tangent-cone point `u* in G` of a
(Mar) fake pattern box gives a `Q` vanishing on `A`, so `AC(u*)` is false there. Any proof of `AC` must therefore
decide the barrier question. No robust or "generic" anti-concentration argument can do so.

No finite test is available either:
* the small models sdc7 and sdc11 have certified nondegenerate tangent-cone points, so `AC` fails there exactly;
* the smallest LP-supercritical Sylvester model has `n >= 32`, beyond exact computation (record).

**Status: not unblocked.** This is the strongest honest link: the new method's *kind* of statement is the kind route
5 needs. But the new method's *setting* (continuous, fixed dimension) is exactly the part that is easy.

### S4. Monge-Ampere potential plus Aubin-Mabuchi energy are the classical tools of the determinant method

**Match [P, given standard pluripotential theory [M] and Berman-Boucksom [L]].** Let `K` be a non-pluripolar compact
set in `C^d` with Robin function `rho_K`.

1. *Homogeneous Bernstein-Walsh.* For `deg p = m`, `log|p(lambda z)| - m log|lambda| -> log|p_hat(z)|`, and
   `V_K(lambda z) >= (1/m) log(|p(lambda z)|/||p||_K)`. Together these give

   ```text
   |p_hat(z)| <= ||p||_K exp(m rho_K(z)).
   ```

   For `Q` with leading form `l_u^m` this gives `||Q||_K >= sigma_K(u)^m`, where
   `sigma_K(u) = sup_z |<u,z>| e^(-rho_K(z))`. **That is a continuous MIN:** a lower bound for the directional
   Chebyshev error, in the sense of Zaharjuta's directional Chebyshev constants [M].
2. *Berman-Boucksom* (Invent. Math. 181 (2010) 337-394, arXiv 0803.1950 [L]). The Monge-Ampere (Aubin-Mabuchi)
   energy of the extremal weight governs the volume of sup-norm balls of sections, i.e. the transfinite diameter and
   the asymptotics of Vandermonde determinants. That is the archimedean half of route 5 (Step 0, the determinant size).

So the paper's two named tools are the classical tools for both halves of the determinant method. This is the closest
conceptual match I found.

**Gap [P].** Transferring 1 from `K` to `A` needs a norming inequality `||Q||_K <= L ||Q||_A` at degree `m'`.
* For `A` inside the ternary box, `L = infinity` for every `m >= 3`, because `k_x^3 - k_x` vanishes on `A` and not on
  any `K` with interior.
* `A` is finite, hence pluripolar. Its degree-`m'` "extremal function" is finite-dimensional linear algebra,
  `J(A)_(m')^perp`, which is precisely what the record studies.
* Consistency check: inequality 1 would give MIN for **all** `u`, but cheap directions exist. So the transfer must
  fail, and it fails exactly through the ternary relations.

The paper's own contribution is a *quaternionic* Monge-Ampere equation on spatial fibres of a wave current. It adds
nothing for finite sets.

**Status: analogy only.**

### S5. Kakeya and incidence inputs for clusters

**Needed.** An existence-level (L^infinity) statement: no circle of radius `R`, centred at a lattice point, carries
`M` lattice points on an arc of length `C sqrt R`. It must use integrality in a way that `Z` and the fakes do not
satisfy.

**Available.** All of these are L^p (`p <= n`) overlap bounds, or counts of *how many* rich curves exist:
* Kakeya and Nikodym maximal bounds with `(log 1/delta)^(4/p_c)` [V];
* Wolff-type circular Kakeya and tangency bounds (2-D local smoothing implies circular maximal bounds) [M];
* Szemeredi-Trotter (as in Herr-Kwak) [L].

**Gap.**
* One rich arc is not excluded by any bound on the number of rich curves.
* The `O(M^2)` chord directions lie within angle `C/sqrt R` of each other. That is the "pair chords" class, which is
  capped.
* `Z` and the fakes satisfy every such bound (W3).

**Status: dead.**

### S6. Multi-shell lifts that escape W1

The idea: build from one cluster a frequency set spread over many shells, where local smoothing is not an identity.

**[P].** For `z_i` on `|z| = R`,

```text
|z_1 + ... + z_k|^2 = k^2 N - (1/2) sum_{i,j} |z_i - z_j|^2 .
```

So the `k`-fold sums of an arc cluster lie in the shell `kR - k C^2/4 - o(1) <= |xi| <= kR`, of thickness at most
`k C^2/4`. For `k < 4/C^2`, the window `T in [1, 2]` sees about one shell: `e^(-iT|D|)` acts as `e^(-iTkR)` times a
multiplier within `O(k C^2)` of 1. So W1 holds approximately.

Larger `k`, or products `z_i w` with small Gaussian integers `w` (shells `R|w|`, effective `C' = C |w|^(1/2)`), do
give many shells. But:
* the resulting inequalities are Euclidean, and the analogous sets built from `Z` or from fakes satisfy them (W3);
* the only arithmetic they could carry is additive coincidences, i.e. `E_k`, which is capped (W2);
* long-time versions (`T ~ R` on `T^2`) would see arithmetic but are not claimed, and are L^p anyway.

**Status: dead.**

### S7. Square-function or decoupling input for route 5's forms

**Needed [P].** Write `Phi_c(u) = i^(-m') (d/dt)^(m') g_c(tu)|_(t=0)`, where `g_c(delta) = sum_k c_k e^(i<k,delta>)`
has frequencies in `A ⊂ {-1,0,1}^d` and vanishes to order `m'` at `0`. MIN^pair asks: for each `u in G`, some
integer `c` of height `<= H` has `|d^(m')/dt^(m') g_c(tu)|_(t=0)| > |c|_1 d^(m') eps_0`. That is a pointwise reverse
Bernstein inequality at a single point.

**Available.** Square-function and decoupling estimates are L^p equivalences for frequency pieces on *curved* sets,
at scale tending to infinity, in fixed dimension.

**Gap.**
* The frequency set is flat: a ternary cube of radius `<= sqrt d` with separation 1. For flat sets, decoupling
  constants are the trivial `(#pieces)^(1/2 - 1/p)`.
* The dimension `d` grows.
* The needed statement is pointwise and a lower bound.

**Status: dead.**

### S8. P3-all\* / Rig^SL: one angle vector in many slabs

**Needed.** By round 14 Lemma J0, there is one `phi in T^r` per monomial and unit inside an integer-normal slab.
Rigidity means that the slabs have **no common point**; construction means that they have one.

**Available.** Kakeya and Nikodym maximal bounds control `||sum 1_plates||_(L^p)` for `p <= n`. They never control
the `L^infinity` (maximal-overlap) norm that decides emptiness.

**Gap.**
* The plates' normals are integer vectors set by labels, i.e. arithmetic. They are not a direction-separated family.
* Feasibility is a Kronecker / geometry-of-numbers question.

**Status: dead.**

### S9. H^gen / L7

The paper contains no primes. Prime-supported analytic inputs are capped by round 10 Thm B whenever they are robust.

The record's one uncapped analytic window is `P_+(N) < M^(2 - delta)` (round 10, `analytic.md` §6). It still needs
the statistics of selected Gaussian-prime residues, which Euclidean geometry cannot see.

**Status: dead.**

### S10. Non-mathematical uses

* The Lean library (`lean/OAI/Analysis/Fourier/...`, about 190k lines in the import closure, not built) may contain
  Fourier-analytic lemmas reusable in our own Lean files (e.g. `GaussianChain/ObtuseBessel.lean`). That is engineering
  only.
* Same-repository entries outside this lens, headlines only, not assessed:
  * 017, the irrationality exponent of pi: one specific linear form in logs, i.e. the capped single-form class;
  * 003, quasi-RH for Dirichlet and Q(sqrt -3) Hecke L-functions: robust, capped by Thm B, wrong field;
  * 028, the Gaussian moat: a sieve-type statement about prime gaps, not angles.

  None bears on the uniform bound as far as headlines show.

## 3. Ranking

| candidate | wall it gets past | wall it hits | outcome |
| --- | --- | --- | --- |
| S3 anti-concentration form of MIN^pair | W4 "direction" (it is an upper bound in `k`, pointwise in `u`); Cor. B.1 | W4 discreteness and degree `~ sqrt(d log r)`; fake Thm 2.1 | best case; the paper supplies only the easy continuous half |
| S4 pluripotential match (Robin function, Aubin-Mabuchi, Vandermonde) | none in substance | the norming constant is infinite at degree `>= 3` on ternary sets | analogy |
| S1 arc filter at scale `R^(-1/2)` | W1, for the proof's auxiliary objects | W2 (L^4), direction (only lower bounds), capped character class | dead, but informative |
| S2 log-loss and divisor precedent | none | W2, W3 (average vs supremum) | dead |
| S5-S9 | none | W2/W3/W4 | dead |

**Probability that the claimed results lead to a growth improvement or unblock a route** [S]: well under 5%. That
residue sits almost entirely on S3, and S3 would need a discrete analogue of the method that the paper does not
contain.

**What would change the verdict:**
* a discrete anti-concentration theorem in the spirit of the paper's positivity method, for polynomials of degree
  `~ sqrt d` on slab-cut ternary cubes, with collision strata as the only exceptional directions;
* a periodic or long-time version on `T^2` would not suffice, because of W2.

Recommendation: record as **not useful**, with S3 noted as the reason not to dismiss positivity methods for route 5
in general. Nothing in the announcement is worth a round.

## 4. Files

* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/steelman.md`
  (this note)
* `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/steel_checks/arcfilter.py`
  (S1 check; `python3 -I arcfilter.py`)

## Sources

* Primary [V], cloned repository `dl/math` (openai/math, commit 7ef5c0c):
  * `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/build/sections/introduction.tex`,
    `barrier.tex` (Prop. fixed-barrier), `pencil.tex`, `current.tex`;
  * `CONTENTS.md` (headlines).
* Herr-Kwak, *Strichartz estimates and global well-posedness of the cubic NLS on T^2*, Forum Math. Pi (2024),
  https://arxiv.org/abs/2309.14275 [L, snippet]
* Guth-Maldague-Wang, *Improved decoupling for the parabola*, JEMS 26 (2024) 875-917,
  https://arxiv.org/abs/2009.07953 [L, snippet]
* Meka-Nguyen-Vu, *Anti-concentration for polynomials of independent random variables*, Theory of Computing 12
  (2016), https://arxiv.org/abs/1507.00829 [L, snippet]
* Berman-Boucksom, *Growth of balls of holomorphic sections and energy at equilibrium*, Invent. Math. 181 (2010),
  https://arxiv.org/abs/0803.1950 [L, snippet]
* Carbery-Wright 2001 (anti-concentration for polynomials, Gaussian/log-concave); Mossel-O'Donnell-Oleszkiewicz 2010
  (invariance principle); Zaharjuta 1975 (directional Chebyshev constants); Klimek, *Pluripotential Theory* (Robin
  function) [M]
* Record: `docs/claude_review/README.md`; `round10/README.md` (Thm B, the open window); `round12/README.md`
  (§1 notation, §2.3 nearstrata duality and MIN^pair, §2.4 fake Thm 2.1, §4.2 J6); `round13/README.md` §0;
  `round14/README.md` §0; `round15/README.md`
* Companion assessments: `facts.md`, `chains.md` (Props 1-2), `lifts.md` (F1-F4) in this directory
