# Adversarial skeptic: does any proposed mechanism survive?

Role: take every *positive* mechanism proposed by the steelman (`steelman.md`, S1-S10) and by the three assessors
(`facts.md`, `chains.md`, `lifts.md`), and try to refute it. Tests: direction of inequalities, losses, dimension, whether
fake circles (or the transcendental configuration `Z`) satisfy the input, whether the 2-D version was already known, and
whether the mechanism actually reaches arcs of length `C sqrt R` with `O(1)` or `o(log R/loglog R)` points. Default:
refuted unless a concrete argument bridges the gap. Written 2026-10-27.

Labels:
* **[V]** read in the primary source (`dl/math`, openai/math commit 7ef5c0c);
* **[P]** proved here;
* **[C]** computed here (`skeptic_checks/`);
* **[L]** literature, search snippet only;
* **[M]** from memory;
* **[S]** speculation.

## 0. Verdict

**Nothing survives. Every positive mechanism is refuted as a way for the announcement to help us.** The steelman's
verdict "not useful" stands. Four of the arguments behind the negative case are wrong or overstated, and are corrected
in §3. None of the corrections helps the positive case.

| proposed mechanism | proposed by | refuted by | §  |
| --- | --- | --- | --- |
| S1 light-ray arc filter at scale `R^(-1/2)` | steelman | W2 volume factor (not "direction"; see §3.1) | 2.1 |
| S2 log-loss / eps-free forms; divisor-bound precedent | steelman | volume obstruction; average vs supremum | 2.2 |
| S3 MIN^pair as anti-concentration AC(u) | steelman (best case) | paper contributes nothing; classical CW gives the easy half; AC is a strict strengthening | 2.3 |
| S4 Monge-Ampere / Aubin-Mabuchi = determinant-method tools | steelman | analogy; Step 0 already proved; `A` pluripolar | 2.4 |
| S5 Kakeya / incidences for clusters | steelman | counts rich curves, not one rich arc; `Z`/fakes | 2.5 |
| S6 multi-shell `k`-fold lifts | steelman | one sector; `L^(kp)` in disguise; `Z` | 2.6 |
| S7 decoupling for route 5's forms | steelman | flat, growing dimension, pointwise lower bound | 2.7 |
| S8 Kakeya plates for P3-all\* | steelman | `L^inf` emptiness, arithmetic normals | 2.7 |
| S9 H^gen / L7 | steelman, chains (entry 003) | no primes; robust inputs capped (Thm B) | 2.7 |
| S10 Lean reuse; entries 017/003/028 | steelman | engineering only; capped classes | 2.7 |
| method as a positivity template | facts, chains | outputs `L^p` upper bounds, capped | 2.8 |
| "revisit if periodic / long-time on `T^2`" | facts | W1 holds at every time length | 2.8 |
| `n >= 3` results via `X_M` lifts | facts | no curvature/Fourier structure in `X_M` | 2.8 |
| single-circle `Lambda(p)` at `p ~ log R/loglog R` ("least remote link") | lifts F2 | chains Prop 2: `>= 4 log R` even with optimal constants | 2.8 |
| tensor lift to `S^3`, exponent-lattice lift, Turan-Nazarov/Remez | lifts | as stated in lifts.md; nothing to add | 2.8 |
| reverse direction (our bound feeds toral eigenfunction bounds) | chains | helps others, not us | 2.8 |

**Probability that the announcement itself leads to a growth improvement or unblocks an open route** [S]: below 1%.
The steelman puts the residue on S3. But the announcement contributes nothing to S3 (§2.3): the continuous half was
already known (Carbery-Wright 2001), and the discrete half is the record's open problem, unchanged.

## 1. What I verified

1. **[V] The light-ray filter and its scale.** `current.tex` §"The symbol determines the Gaussian": translation by
   `-i rho(1, omega)` multiplies a mode by `exp(-rho |xi| |xi_hat - omega|^2/2)`, with angular width `(rho N)^(-1/2)`.
   `barrier.tex` line 72: the energy's Gaussian is `Psi = (1/2) kappa y |W - v|^2` with `y = |Y| ~ s = |tau|`
   (`energy.tex` line 40). So at `|tau| ~ 1` the aperture is `N^(-1/2)`, the canonical-sector scale, which equals our
   arc scale. **The steelman's unverified identification in S1 is correct.**
2. **[V] What the paper proves at general `s = |tau|`** (Prop. fixed-barrier, `barrier.tex` lines 153 and 229), with
   `sigma = s/h`:
   ```text
   F(T,tau) <= C a^2 sigma^m ||g_T||_pc^pc          (upper, every s; lossy by N^m at s ~ 1)
   F(T,tau) >= c a^2 sigma^(-m) e^(-C3 gamma^2 sigma) ||g_T||_pc^pc   (lower, every s; e^(-cN) at s ~ 1)
   ```
   Both are useless at `s ~ 1`. Only the fixed scaled radii `2K` and `4K` are used. The circle passages are Lemma
   `lem:cylinder` (`propagation.tex`): `H(0,r) >= c(M) log(r/b)/D^2 * int_1^2 H(T,b) dT`, a lower-bound propagation.
3. **[V] Hypotheses of the core estimate:** frequencies in a fixed annulus `c_1 <= h|xi| <= c_2` and a fixed small
   aperture `|v| <= gamma`, on `R^n/(100 Z^n)` with `T` in `[1,2]`. A single circle's lattice points in one
   `gamma`-cap qualify. A single circle is the case where the real-time estimate is a scalar identity (W1).
4. **[C] Carbery-Wright room in S3** (`skeptic_checks/cw_room.py`, §3.2): at the regime threshold the
   Gaussian-proxy bound is **> 1/2**, so the "enormous room" claim holds only asymptotically.
5. **[L] Herr-Kwak.** Forum Math. Pi 2024, arXiv 2309.14275: sharp `L^4` Strichartz on `T^2` via Szemeredi-Trotter.
   It improves Bourgain's divisor-type loss; the extremal example shows `(log N)^(1/4)` growth. This confirms the
   steelman's S2 citation (abstract-level only). Guth-Maldague-Wang was not re-checked.
6. **[L] Expert reaction.** A fresh search found no expert assessment of, and no reported gap in, the local smoothing
   claim. The results were about other OpenAI claims: unit distance, Navier-Stokes, the 722 Erdos papers.

## 2. Refutations

### 2.1 S1: the arc filter at `|tau| ~ 1`

Claim: the proof's auxiliary energy, unlike the real-time inequality, is non-vacuous on one circle, and acts at exactly
our scale.

* **Granted:** the identification (§1.1), and that `max_omega S_omega(1)` is comparable to `M`.
* **The steelman's "direction" objection is wrong** (§3.1): for single-shell data, an *upper* bound at the arc scale
  does follow from the paper's inequalities.
* **Refutation (W2, volume):**
  * The energy at `|tau| = s` is an integral over the whole fiber: real `X` and complex direction `W`. So the arc
    filter enters only *averaged over directions*.
  * **[S/P]** Grant the most favourable lower bound the construction could give at `s ~ 1` (not stated in the paper):
    `H(0,1) >~ a^2 avg_W int |g_W|^4`, where `g_W` is the filtered wave at aperture `R^(-1/2)` (`p_c = 4` in the plane).
  * The cluster contributes only on `W` within angle `~ R^(-1/2)` of its centre, times `x` in a dual box of area
    `~ R^(-1/2)`. That gives `avg_W int |g_W|^4 >= c r'^(-2) M^4 R^(-1)` for normalised data on `r'` cap points.
  * Combined with `H(0,1) <~ log R * a^2 ||f||_4^4` and `||f||_4^4 <= 2` (additive energy of points on a circle),
    this gives `M <~ R^(1/4) r'^(1/2) (log R)^(1/4) = R^(1/4 + o(1))`.
  * The steelman's `R^(1/8) r^(1/2)` omits the direction average. Either way the bound is useless, and the obstruction
    is the measure of the set where the filtered wave is large, not a loss.
* **Capped class:** the filter's Bessel expansion is a weighted Hecke-character sum, so it is in the record's capped
  character class.
* **Fakes / `Z`:** `S_omega(l)` is defined for any real point set on the circle, and the inequalities hold for all data.

**Refuted.**

### 2.2 S2: eps-free or log-loss forms, and the divisor-bound precedent

* The best form the paper gives is `(log N)^(2/p_c)` at `p_c`. The paper's own Prop. 8.1-8.2 (per `chains.md`, `facts.md`)
  show that *some* log loss is necessary there.
* On a single circle the real-time inequality is an identity (W1), so the loss is irrelevant: there is nothing to lose.
* For any `L^p` chain on `T^2`, chains Prop 2 caps the output at `Theta(log R)`, which is worse than the proved
  `(2/3) log R/loglog R`. The obstruction is the box volume `V ~ C^3 sqrt R`, so it is loss-free.
* **Precedent.** Herr-Kwak [L] beat the divisor bound `exp(c log N/loglog N)` in an `L^4` Strichartz estimate,
  which is an average over the resonance (circles weighted by time integration). Our `M` is a supremum over a single
  arc, and arc counts are an *input* to such estimates. The precedent shows these methods can improve averages of
  divisor-type quantities. It gives no instance of a supremum of a single-shell count being improved. The one known
  route to single-circle structure (Cilleruelo-Cordoba) is arithmetic.

**Refuted.**

### 2.3 S3: MIN^pair as anti-concentration (the steelman's best case)

Claim: route 5's residual inequality is an anti-concentration statement AC(u), "the same kind of statement" the new
method proves.

1. **The paper contributes nothing to either half.**
   * *Continuous half.* Carbery-Wright (2001) [M] already gives the Gaussian- or log-concave-proxy statement. The paper
     is neither needed nor used.
   * *Discrete half.* This is the record's open problem: the algebra of `J(A)`, an inverse Littlewood-Offord theorem at
     degree `~ sqrt(d log r)`, decided on fake boxes by fake Thm 2.1. It is untouched.
2. **"Same kind of statement" conflates two different anti-concentrations.**
   * The paper's positivity argument (Alexandrov: density `mu >= m` on a ball forces Monge-Ampere depth) bounds the
     *large values* of a physical-space density. Its output is an `L^(n/(n-1))` bound on `mu = |wave|^2`, i.e. an
     `L^(p_c)` bound on the wave, in fixed dimension, for Lebesgue measure.
   * AC(u) is a *small-ball* bound at one scale `delta = d^(m') eps_0` for the pushforward of counting measure on
     `{-1,0,1}^d` under a polynomial of growing degree. It is an `L^infinity`-at-a-point statement about a discrete
     measure.
   * The norms, the measures and the dimension regime all differ. There is no reduction either way.
3. **AC is a strict strengthening of what route 5 needs, and is proved for the wrong functional** [V, record round 12
   nearstrata §1-2].
   * AC quantifies over every half `S` of the cube, not just over the actual slab box `A`.
   * AC controls `kappa_V = kappa` (all of `V(A, m')`). Thm 4.1 uses `kappa_H`, and `kappa_H <= kappa_V`, with
     equality only if `V_H = V(A, m')`. Step 0 holds only on `V_H` (integrality of `c` is essential), so the
     `H = infinity` form does not feed Thm 4.1 without the unproved equality.

   The steelman notes the second point as a caveat. It is a genuine gap, not a technicality.
4. **Even the "easy" continuous half is only asymptotic** (§3.2). At `W = 3 n ln n` the Carbery-Wright bound is
   `> 1/2` for `n` below about 6,000 (`r = n`) or about 12,000 (`r = n^2`). Route 5's machinery starts at
   `n >= 2900 / 5809` (`m_max >= 4w`).
5. **Fake test.** Any proof of AC(u*) at a nondegenerate real `TC_H` point `u*` of a (Mar) fake box would contradict
   fake Thm 2.1. So a proof must decide the barrier question. No robust or continuous method, including this paper's,
   can do that.

**Refuted as a contribution of the announcement.** AC itself is a reformulation of an open problem that predates the
announcement, and is stronger than needed.

### 2.4 S4: pluripotential tools are the determinant method's tools

* The match itself is correct [M/L]: homogeneous Bernstein-Walsh gives a continuous MIN, and Berman-Boucksom relate the
  energy at equilibrium to transfinite diameter and Vandermonde asymptotics. All of this is classical, dating from
  1975-2010, and needs nothing from the paper.
* The archimedean "determinant size" half (Step 0) is **already proved** in the record (round 12 t1rigour Thm 4.2).
  Supplying it again gains nothing; the bottleneck is MIN^pair.
* Transfer to `A` fails: `A` is finite (pluripolar), and `k_x^3 - k_x` kills every norming inequality at degree
  `>= 3`. The steelman proves this.
* The paper's quaternionic, fibrewise Monge-Ampere equation for a wave density has no finite-set analogue in it.

**Refuted (analogy only).**

### 2.5 S5: Kakeya, Nikodym, circular maximal, Szemeredi-Trotter

* These are `L^p` overlap bounds for `p <= n`, or counts of how many curves are `k`-rich. One `M`-rich arc is
  consistent with every such bound.
* The chord directions of a cluster form the capped pair-chord class.
* `Z = {R e^(ik/R)}` has `C sqrt R` points on such an arc and satisfies every incidence bound that does not use
  integrality.
* 2-D Kakeya, with log loss, is classical (Cordoba 1977), and Cordoba's `log^(1/2)` beats the paper's 2-D log power.

**Refuted. Not new in 2-D.**

### 2.6 S6: multi-shell `k`-fold lifts

* [P, rechecked] `|z_1 + ... + z_k|^2 = k^2 N - (1/2) sum_(i,j) |z_i - z_j|^2`, so the shell thickness is
  `<= (k-1) C^2/4`.
* All `k`-fold sums lie in one *angular* sector of width `C/sqrt R`, i.e. `~ C sqrt k` canonical sectors at radius
  `kR`.
* Local-smoothing gains come from many sectors. With `K' = C sqrt k` sectors, the gain is at most `K'^O(1)`.
* Since `||F^k||_p = ||F||_(kp)^k`, a gain `G` on `F^k` is a gain `G^(1/k) -> 1` on `F`, in `L^(kp)`. That is W2 again,
  capped at `Theta(log R)`.
* `Z`'s `k`-fold sums are all distinct (transcendence of `e^(i/R)`), so `Z` saturates nothing.

**Refuted.**

### 2.7 S7-S10

* **S7.** Route 5's forms have frequencies in a flat ternary cube, with unit separation, in growing dimension `d`,
  and the requirement is a pointwise *lower* bound at one point. Decoupling for a flat set is trivial, and the paper's
  constants `C_(n,p,eps)` have no stated dependence on `n`. **Refuted.**
* **S8.** P3-all\* rigidity is emptiness of an intersection of slabs with integer, label-determined normals: an
  `L^infinity` / Kronecker question. Kakeya bounds `||sum 1_plates||_p` for `p <= n` only. **Refuted.**
* **S9.** The paper contains no primes. Entry 003 (quasi-RH for Dirichlet and finite-order Hecke L-functions over
  `Q(sqrt -3)`) is the wrong field and the wrong characters: we need the infinite-order angle characters of `Q(i)`. It
  is also a robust ensemble input, capped by round 10 Thm B. **Refuted.**
* **S10.** Lean reuse is engineering. Entry 017 is a single linear form in logs (capped class). Entry 028 is about
  prime gaps, not angles. **Refuted.**

### 2.8 Assessor-proposed residual angles

* **Positivity template** (facts, chains). The method converts density concentration into an `L^(p_c)` upper bound. On
  one circle the real-time output is an identity (W1). On many circles it is an average, capped by chains Prop 2. It
  produces no auxiliary polynomial and no statement about `J(A)`. **Refuted.**
* **"Revisit if the method is extended to periodic or long-time settings on `T^2`"** (facts). This trigger is wrong.
  * For data on one circle, `e^(-iT|D|)` (and `e^(iT Delta)`) acts by a scalar phase **for every `T`**, so a long-time
    estimate is still an identity on single-shell data.
  * Arithmetic enters only through many shells, i.e. Strichartz-type averages, where arc counts are inputs and chains
    Prop 2 caps the output at `Theta(log R)`.

  A periodic, long-time extension would therefore not change the verdict either. **Refuted.**
* **`n >= 3` results via equal-norm varieties `X_M`** (facts). `X_M` questions are about rational points and
  Diophantine structure; they have no curved-hypersurface Fourier structure. No implication is known or sketched.
  **Refuted.**
* **Single-circle `Lambda(p)` at `p ~ log R/loglog R`, the "least remote link"** (lifts F2). Even a hypothetical optimal
  bound cannot help.
  * The box lower bound plus `||F_cluster||_p <= K_p sqrt M` gives `M <= 4 K_p^2 V^(2/p)`.
  * At `p = log R/loglog R`, `V^(2/p) ~ R^(1/p) = log R`, so `M <= 4 K_p^2 log R >= 4 log R`, since `K_p >= 1`.
  * Minimising over `p` (with `K_p >= c sqrt p` beyond `p ~ 1`) gives `Theta(log R)` (chains Prop 2).

  So F2's "least remote link" is dead even granting an optimal theorem. This should be recorded as a cap, not a link.
  **Refuted.**
* **Tensor lift to `S^3`, exponent-lattice lift, Turan-Nazarov / Remez** (lifts). The refutations in `lifts.md` stand.
  The exponent-lattice lift is a flat slab with only Diophantine content (capped). Turan-Nazarov / Remez lose
  `exp(3^d)` or give only sublevel measure (capped by Cor. B.1). **Refuted.**
* **Reverse direction** (chains). Our arc bound may feed toral eigenfunction or nodal estimates. That helps other
  problems, not ours, so it fails the usefulness criterion by definition.

## 3. Corrections to the negative case (accuracy; none helps the positive case)

### 3.1 S1 "direction" is not a valid refutation

The steelman says the paper only carries lower bounds through `|tau| ~ 1`, so no upper bound is available there. For
single-shell data this is wrong.

* **[P, modulo one inference]** If `f` has frequencies on one circle, then `g_T = e^(-iT rho) f`. Every mode carries
  the same unimodular factor, also after the complex shifts. So the fiber density `mu` and the energy `F(T,tau)` are
  independent of `T`.
  * *The inference:* `T` enters `F` only through `g_T`. I read this from `current.tex` and `energy.tex` but did not
    check every line.
* Then `Delta_tau F <= V F` on `2a <= |tau| <= 3`, and the circular mean `H(s)` is a radial supersolution.
* Compare on `[2a, 1]` with the solution of `phi'' + phi'/s = V phi`, `phi(2a) = 0`, `phi(1) = H(1)`. Using `F >= 0`,
  the minimum principle, and the moment bound of `lem:cylinder`:
  `H(4a) >= c log 2 / log(1/2a) * H(1)`. Hence

  ```text
  H(1) <~ log N * H(4a) <~ log N * a^2 ||f||_pc^pc       (by barrier-input at |tau| = 4a).
  ```

So there *is* an upper bound at the arc scale. S1 still dies, but by W2 (§2.1), not by direction.

### 3.2 S3's "enormous room for every `u`" is asymptotic only

[C] `skeptic_checks/cw_room.py` uses the record's `eps_0 = C d e^(-W/4)/(m'+1)` with `C = 1`, a Carbery-Wright
constant of 1 (optimistic), `sigma >= sqrt(m'!)`, and `m' = ceil(4 sqrt(8 d ln 4r))`.

| `n` | `r` | `W` | `m'` | `W/(4m')` | log10 of the CW bound (AC needs `< -0.30`) |
| --- | --- | --- | --- | --- | --- |
| 2900 | `n` | `3 n ln n` | 1864 | 9.30 | **+1.27** |
| 5809 | `n^2` | `3 n ln n` | 3731 | 10.12 | **+1.37** |
| 10^4 | `n^2` | `3 n ln n` | 5035 | 13.72 | **+0.11** |
| 10^4 | `n` | `3 n ln n` | 3683 | 18.76 | -2.15 |
| 10^5 | `n` | `3 n ln n` | 12850 | 67.20 | -21.91 |
| 2900 | `n` | `10 n ln n` | 1864 | 31.01 | -8.15 |

The crossover at `W = 3 n ln n` is `n ~ 6,000` (`r = n`) and `n ~ 12,000` (`r = n^2`). It is larger with a realistic
Carbery-Wright constant. So the continuous proxy is easy only for `n` beyond the record's thresholds, or for
`W >> n ln n`. This does not affect the verdict.

### 3.3 lifts F2 should be recorded as a cap, not as the "least remote link"

See §2.8: chains Prop 2 already shows that the hypothetical theorem F2 asks for yields at best `>= 4 log R`.

### 3.4 W2's scope

Chains Prop 2 is proved for the box-lower-bound plus `Lambda(p)` / `k`-fold-energy chain applied to the cluster's own
exponential sum, or a sum containing it. "All `L^p` / energy chains are capped at `Theta(log R)`" is broader than what
is proved: weighted or multi-shell `L^p` arguments that use integrality in other ways are not formally covered.

This does not matter here. The paper's inequalities hold for all data, so in any such argument the arithmetic content
comes from the other inputs, and W1 / W3 apply.

## 4. What would change the verdict

Nothing in harmonic analysis of this type, including endpoint, loss-free, periodic or long-time versions (§2.2, §2.8).

Only results of a different kind would matter:
* a *discrete* small-ball (inverse Littlewood-Offord) theorem for polynomials of degree `~ sqrt(d log r)` on slab-cut
  ternary cubes, with collision strata as the only exceptional directions, **for `kappa_H`** (or a proof that
  `V_H = V(A, m')`);
* a statement about Gaussian-prime residue statistics (H^gen / L7) that is not robust;
* a construction or rigidity theorem for P3-all\*.

The announcement contains none of these, and its method has no visible discrete analogue.

Recommendation: record the announcement as **not useful**. No round is warranted. Record the S1 identification
(verified) and the §3 corrections for accuracy.

## 5. Files and sources

* This note: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/skeptic.md`
* Check:
  `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/local_smoothing/skeptic_checks/cw_room.py`
  (`python3 -I cw_room.py`)
* Primary source [V], openai/math commit 7ef5c0c, local clone `dl/math`,
  `preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/build/sections/`:
  * `introduction.tex` (Thm 1.1, Thm 1.2, §"A conserved current and two circles");
  * `current.tex` (Gaussian width `(rho N)^(-1/2)`, cap hypotheses);
  * `barrier.tex` (Prop. fixed-barrier; (barrier-upper-all), (barrier-lower-all) with `sigma = s/h`);
  * `energy.tex` (`y ~ s`);
  * `propagation.tex` (Lemma `lem:cylinder`).
* Herr-Kwak, *Strichartz estimates and global well-posedness of the cubic NLS on T^2*, Forum Math. Pi (2024),
  https://arxiv.org/abs/2309.14275 [L]
* Search for expert reaction (none found on local smoothing): eweek.com, scientificamerican.com, letsdatascience.com
  results about other OpenAI claims [L]
* Carbery-Wright (2001); Berman-Boucksom (Invent. Math. 2010); Cordoba (1977) [M]
* Record: `docs/claude_review/README.md`; `round12/README.md` §1, §2.3 (`kappa_H <= kappa_V`, equality iff
  `V_H = V(A,m)`), §2.4 (fake Thm 2.1); `round12/nearstrata.md` lines 79-166
* Companion notes: `steelman.md`, `facts.md`, `chains.md` (Props 1-2), `lifts.md` (F1-F4)
