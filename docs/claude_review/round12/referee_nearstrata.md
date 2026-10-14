# Referee report on `round12/nearstrata.md` (near-strata Lojasiewicz bound, arithmetic separation, assembly)

Referee scripts: `round12/referee_nearstrata_checks/` (each `x.py` writes `x_out.txt`; none imports the author's code).

## 0. Verdict

**The proofs survive; several framing claims do not.** I re-derived every numbered statement. Lemma 1.1, Props 2.1 and 2.6-2.7,
Lemmas 2.2-2.3 and 2.5, Theorems 2.4-2.5, Lemmas 3.1-3.3, Theorem 3.4, Cors 3.5-3.6, Theorems 4.1 and 4.3, Cor 4.2,
Lemma 5.0 and Props 5.1-5.3 are correct as stated. I found no wrong quantifier in any theorem, and no cheap subtorus that
Theorem 3.4 misses. Every scan threshold is reproduced exactly (§3).

Nothing here is a growth improvement or a uniform bound, and the write-up says so. The defects are in what the results
are said to mean:

* **(Major) MIN\* is not the minimal statement relative to the pair data (§4.3).** The pair excesses
  `eps_xy = d_xy - W/2 + 2 ln C'` are fixed by the profile, its primes and `W`. Lemma 3.1 gives each pair its own gap,
  `|hat u_x - hat u_y| >= e^{-eps_xy/2}`. MIN\* keeps only their average (`Gamma <= abar`). The statement that is minimal
  relative to the pair data is the profile-specific **MIN^pair** (§4). It is implied by (MIN) and strictly weaker.
  A concrete example with genuine primes (Sylvester `k = 16`, `r = 65535` primes `p = 1 mod 4` in
  `[e^25, e^25 + 6*10^6]`) makes this exact:
  - (MIN) is false at every order in the Prop 5.2 window, since `Lambda = 0` at `g >= 3` and `abar = 12.5`.
  - The individual gaps alone admit no `u` at all: `(n-1) e^{-max eps/2} = 104.9 > 2`.

  So Prop 5.2 ("(MIN) is false at large orders") and Prop 5.1 obstruct the *averaged* route only. They say nothing about
  the pair-data route, and "the open content of route 5 at `T = 1` is MIN\*" is overstated.
* **(Moderate) "(MIN) forces `P(L) >= o_L n/2`" (summary, §0 item 7, §7.2) is wrong as worded.** Prop 5.1(ii) proves this
  for MIN\* at *all* large levels `g`. The constant `Gamma_L` depends on the stratum and is not controlled, so at the single
  level `g = abar` (Theorem 4.1), or at the levels `g <= (1/c + o(1)) ln n` that suffice for the growth question (§4.3),
  no such necessary condition follows.
* **(Moderate, new information) MIN\* is provably false for a regime profile at `m' = ceil(4w)`.** Take the Sylvester
  profile with `n = 4096`, `r = 4095`. It satisfies (R1) (`ceil(4w) = 2256 <= m_max = 2631`). Its 3-column span is cheap,
  with cost `6j = 1686`, so its proved transversal order is `>= 570`, and `P = 1,046,528`. Then
  `o N_2/P >= 4567.8 > d = 4095`, and Prop 5.1(ii) with explicit constants shows that MIN\* fails for every
  `g > 328.4`. This does not contradict any theorem: Cor 4.2 quantifies over `n >= n_1`, and the author's
  `n_B = 5185` already signals it. But it forces `n_1 > 4096`. It also shows that the hypothesis of Cor 4.2, "(MIN) for
  every `W > W_1`", fails through clusters with absurdly large `W` (`W > 2.7*10^6` at `n = 4096`), not through
  near-extremal ones.
* **Minor** (details in §2):
  - The sufficient condition "`theta < d` and `ln(1/gamma) <= (3/8 - o(1)) n ln n`" (§0 item 5, Remark 2.8 reading, §4.3,
    and the round summary) is false for `theta in (d/2, d)`. The correct budget is `(1 - theta/d)(3/4 - o(1)) n ln n`;
    Theorem 4.3(a) itself is right.
  - In Prop 5.3 the asymptotic constants are attached the wrong way round. The order-free condition is `n > (64+o(1))w`
    and the proved-order condition is `n > (8+o(1))w`.
  - Cor 3.7: `rho_F^min` is a lower bound for the smallest admissible radius, not that radius.
  - The thresholds are not monotone in `r`. At `n = 5185` the proved-order condition holds for `r = n` but fails for
    `r in [5138, 5177]`, so "`n >= 5185` (`r = n`)" cannot be read as "`r <= n`".
  - Prop 5.2 asserts `Gamma <= 3` for all `k` without a proof of monotonicity. This does not matter, since `abar > 3.4`.

## 1. What was checked and how

| item | method | result |
| --- | --- | --- |
| Lemma 1.1 (Hoeffding, Siegel count, rigidity, Step 0) | by hand against [D] Thm 2.1(c) and `theory.md` Lemma 2.4 (sharp remainder `|x|^{m+1}/(m+1)!`, so Step 0 holds exactly) | correct |
| Prop 2.1 (duality), Lemma 2.2, Lemma 2.3 | by hand | correct |
| Thm 2.4 (uncut cube) | by hand | correct |
| Lemma 2.5 (strata in `TC`) | by hand, including the Newton-interpolation degree `s + L` on `A2`; exact nullspace check of `V(A_{S,L}, m)` on 7 small boxes (`onecut_tc_exact.py`): every sampled point of types I/II/III is an exact common zero; 42 random points off `Sigma` are not | correct |
| Thm 2.5(a) gadgets | by hand (hexpair's 8 points, all four cases); **exhaustive** check over all class structures with `s <= 9`, all `1 <= L < s` (`gadget_bruteforce.py`, 1795 instances): the maximum gadget degree at cost `<= L` equals `tau = min(2 sigma, 2L+2t, s+L)` exactly (so `nu_eta` is also optimal for this gadget set) | correct |
| Thm 2.5(b), (2.1) | by hand | correct |
| Props 2.6, 2.7 | by hand (BCR 2.6.7 hypotheses: `kappa_H^{-1}(0) in h^{-1}(0)` is ZE-TC^braid) | correct; the "Reading" after 2.7 is an interpretation, correctly argued |
| Lemmas 3.1-3.3, Thm 3.4, Cors 3.5, 3.6 | by hand (gcd norm `N e^{-d_xy}`; convexity of the Plotkin sum; Markov; window/outlier count; `(m'-1)/(2j) < 5` for `w >= 10`) | correct |
| Cor 3.7 thresholds | own scan (`scans.py`) | reproduced: 5185 / 6425 / 10280 |
| Prop 5.3 thresholds | own scans (`scans.py`, `nA_scan.py`, exact rational recheck) | reproduced: `n_A` = 473,716 / 564,084 / 947,350; `n_B` = 5185 / 6398 / 10280 |
| (R1) thresholds | exact `m_max` by big-integer binomials (`r1_check.py`) | reproduced: 2900 / 3623 / 5809 |
| Prop 5.2 numbers | exact `m_max` (`params.py`); `Gamma` by closed form and by direct summation (`sylv_gamma.py`) | reproduced: `2jk = 40896 < m_max = 42065` at `k = 16`; `26340 > 21034` at `k = 15`; `Gamma = 2.99962` |
| Thms 4.1, 4.3, Cor 4.2 | by hand | correct implications |
| Lemma 5.0, Prop 5.1 | by hand (constants `4^m`, `Gamma_L`, choice of `eta`) | correct |

## 2. Section-by-section findings

### 2.1 §1 and §2 (Lojasiewicz side): correct

* **Lemma 1.1.**
  - Hoeffding: the ranges are `<= 4`, so the bound is `2 exp(-w^2/(8d)) = 1/(2r)`.
  - The rigidity inequality `(m'/2 - w) W <= 2 ln|c|_1 + 2m' ln(dC') - 2 ln m'!` follows from `|F_c|^2 = Norm F_c`.
  - With `m' >= 4w` this gives `W <= W_1`.
  - The Siegel count matches [D] Thm 2.1(c) at `B = 1`; the bound `>= U/2` needs `ln U >= 3.3`, which holds.
  - Step 0 holds exactly, because `theory.md` Lemma 2.4 uses the sharp remainder.
* **Thms 2.4-2.5.**
  - I checked the case analysis line by line.
  - The exhaustive check confirms the combinatorial core: no class structure with `s <= 9` defeats the gadget count, and
    the count is attained exactly.
  - Part (b) and the deduction of `TC_H = Sigma(S,L,m)` are correct.
  - These are statements about zonotopal model boxes. The author says so (Remark 2.10) and does not present them as
    evidence about regime boxes.
* **Prop 2.6(ii)** is the BCR Lojasiewicz inequality for continuous semialgebraic functions on a compact set, applied
  correctly.
* **Prop 2.7.** The test point satisfies every claimed bound. The pair `{3, n}` is the only short pair, and
  `Gamma <= g_0 + (2/N_2) ln(1/tau)`.
* **Remark 2.8: minor error.** The display is right, but the gloss "`theta < d` and
  `ln(1/gamma) <= (3/8 - o(1)) n ln n`" (repeated in §0 item 5, §4.3 and the round summary) is right only for
  `theta <= d/2`. With `W >= 3n ln n - 14n` (T3, `C <= 1`), Thm 4.3(a) needs
  `ln(1/gamma) < (1 - theta/d)(3/4) n ln n - O(n)`. For `theta` close to `d` the budget tends to `0`.

### 2.2 §3 (arithmetic side): correct; Corollary 3.7's wording overstates

* **Lemma 3.1.**
  - `Norm gcd(z_x, z_y) = N e^{-d_xy}`: the split primes give `p^{e_p - |a_x - a_y|}`, and the inert primes and `1 + i`
    appear to their full common power.
  - Hence `e^{-d_xy/2} <= |delta_x - delta_y| <= C e^{-W/4}`.
* **Lemma 3.2.** The Plotkin sum is convex in each `a_{p,x}`, so `sum_{x<y} d_xy <= (n^2/4) sum_p e_p log p <= n^2 W/4`.
  Then `sum eps_xy <= nW/4 + 2N_2 ln C' = N_2 abar`. No hypothesis on `N` is used.
* **Lemma 3.3, Thm 3.4, Cor 3.5.** Correct.
  - Theorem 3.4 holds for every integer direction, so rescaled profile columns (for example `1_S` in place of a sign
    column `beta = 2*1_S - 1`) are covered.
  - **Twin and duplicate effects.**
    - Duplicate columns change only `r`.
    - Twin points cannot occur, because equal valuations force `z_x = z_y` on a short arc.
    - A near-twin direction `e_x - e_y` has width `<= 4`. It is an outlier direction (`ell = 0`, `e = 2`) and is counted
      in `E`.
    - Near-duplicate columns `v_1, v_2` span the same subspace as `v_1, v_1 - v_2`, and the latter is an outlier
      direction.

    So no width-cheap span escapes the class count `q_F <= 16`. Width-cheap is not the same as projection-cheap:
    spans whose projection `pi(A)` has interpolation degree below `m'` for reasons other than widths are not covered.
    The author lists this as open (§7.2(b)).
* **Cor 3.6.** Correct. The type-II count uses two groups of total size `n - |X|`.
* **Cor 3.7.** The scan is correct and reproduced. Two wording problems:
  - **(i) What `rho_F^min` is.** Lemma 5.0 gives `kappa_H <= 4^{m'} eta^{o_F}`, so a near-strata radius `rho_F` must
    satisfy `rho_F >= (eps_0 4^{-m'})^{1/o_F}`. Because `o_F >= m' - d_F` and `eps_0 4^{-m'} < 1`, this is
    `>= rho_F^min`. So `rho_F^min` is a lower bound for the smallest admissible radius, not that radius, and
    `sigma_F > rho_F^min` is a necessary compatibility condition only. The order-free version is the `n_A` condition of
    Prop 5.3, which needs `n >= 473,716` (`r = n`).
  - **(ii) Which `r` the thresholds cover.** The thresholds are for `r in {n, n ln n, n^2}` only, and they are not
    monotone in `r` (`r_monotone.py`). At `n = 5185` the proved-order condition (the `W`-coefficient condition B of
    Cor 3.7) holds at `r = n`, with margin `-33`, but fails for every `r in [5138, 5177]`. A statement for every cluster
    needs a scan over all `r <= W/ln 5`. One more detail: at `r = n ln n` and `n = 6424` the `W`-coefficient condition
    holds with exact equality (`o Q = d`, exact rational), so the stated failure "just below 6425" is an artefact of
    requiring strict inequality.

### 2.3 §4 (assembly): implications correct; §4.3's minimality claim is false

* **Theorem 4.1, Cor 4.2 and Thm 4.3(a)-(c)** are correct implications. I re-derived the inequalities in (a) and (c).
  Cor 4.2 also needs `n >= n_2(C)` for (R2), and the author notes this.
* **§4.3, "Why this is the minimal statement relative to the arithmetic used". This is false.**
  - **What the pair data actually say.** Lemma 3.1 is a family of `N_2` individual inequalities
    `|hat u_x - hat u_y| >= e^{-eps_xy/2}`. Each `eps_xy = d_xy - W/2 + 2 ln C'` is computable from the profile
    (`a_{p,x}`, `e_p`, `log p`) and `W`, which are the same data that define the slab box `A`.
  - **What MIN\* keeps.** MIN\* keeps only the averaged consequence `Gamma <= abar`. The author's argument treats the
    `eps*_xy` as free, so it shows minimality only for *profile-blind* pair data.
  - **What a proof may use.** A proof is entitled to the profile-specific gaps, and Theorem 4.1 holds verbatim with them
    (§4, MIN^pair).
  - **The example.** In the Sylvester `k = 16` example of §0, (MIN) fails because the averaged sublevel set contains
    the collision-free binary-fraction point. The individual gaps admit no point at all: `eps_xy in [12.4997, 12.8739]`,
    so every pair needs `|u_x - u_y| >= 1.6*10^{-3}`, and 65536 points do not fit in `[-1,1]` with that spacing
    (`sylv_pairdata.py`).
  - **Consequences.**
    - The "minimal open statement" is MIN^pair, a statement about the profile *with its weights*, not "about the
      profile's slab box alone".
    - Props 5.1-5.2 are obstructions to MIN/MIN\* and not to MIN^pair. A collision-poor stratum `L` defeats MIN^pair
      only if some point of `L`, after perturbation, meets every individual gap. The ordering of points along `L`
      must then be compatible with which pairs have large `eps_xy`.
  - **Scope of the example.** It proves only that MIN\* is not minimal. It does not show that MIN^pair escapes the
    obstructions for near-extremal regime profiles. For those, the individual `eps_xy` must be very non-uniform
    (sorted neighbours need `eps >~ 2 ln(n/2)`, against the average `abar ~ 1.5 ln n`). Whether that non-uniformity is
    compatible with collision-poor strata is a new, well-posed question.
* **"Qualitative core" and "open content of route 5 at `T = 1`".** These should be restated relative to MIN^pair. With
  the individual gaps the qualitative core is not ZE-TC^braid on all of `{Gamma < infinity}`. It is the absence of
  zeros of the `Phi_c` in the closed gap region `{|u_x - u_y| >= e^{-eps_xy/2}}`, which for the actual profile is a
  compact set disjoint from `Br`. On that set, `Lambda > 0` is exactly the statement "no `Phi_c`-common zero meets every
  individual gap", which is weaker than ZE-TC^braid.

### 2.4 §5 (necessary conditions): correct propositions, overstated summaries

* **Lemma 5.0, Prop 5.1.** Correct.
  - The choice `eta = (3n/2) exp(-(N_2/(2P))(g - a_y))` meets `eta <= delta_y/4` exactly when `g >= Gamma_L`.
  - **But the necessary condition `o N_2/P(L) <= d` is derived only for MIN\* as `g -> infinity`.**
    - `Gamma_L = a_y + 2 ln(6n/delta_y)` is not controlled; it depends on the least gap of a generic point of `L`.
    - So "(MIN) forces `P(L) >= o_L n/2`" (round summary, §0 item 7, §7.2) holds only for MIN\* at all levels, which is
      the hypothesis of Cor 4.2. It does not hold for (MIN) at one level `abar`.
    - Nor does it hold for the bounded-level version (`g <= (1/c + o(1)) ln n`) that §4.3 says suffices for the growth
      question. For the growth question the necessary condition is not established.
  - In the example of §2.5 the failure starts only at `g > 328`.
* **Prop 5.2.** Correct. My values: `m_max(2^16) = 42065`, `2jk = 40896`, `m_max(2^15) = 21034`, `2jk = 26340`, and
  `Gamma = 2.999621` at `k = 16`. `Gamma(n)` increases towards `3` (2.0665, 2.8263, 2.9838, 2.99953, 2.99962, 2.99980 at
  `n = 8, 64, 1024, 4096 (2^12), 2^16, 2^17`). Monotonicity is not proved, but `abar >= 3.4` makes it irrelevant.
  "(MIN) fails for every cluster with this profile" is true, and vacuous for weightings like those of §0.
* **Prop 5.3.** The thresholds are reproduced. As noted in §0, the asymptotic constants are attached the wrong way
  round. The checks are:
  - `n_A = 473,716`: `w = 7401` and `64w = 473,700`.
  - `n_B = 5185`: `w = 642` and `8w = 5136`.

  So the order-free condition is `n > (64 + o(1)) w`, and the proved-order condition is `n > (8 + o(1)) w`.

### 2.5 A regime profile on which MIN\* is provably false at `m' = ceil(4w)` (new)

Take the Sylvester profile with `k = 12`: `n = 4096`, `r = 4095`, columns `1[a.x = 0]`, `e_p = 1`
(`sylv12_minstar.py`, `sylv4col_out.txt`).

1. **Regime.** `w = 563.82`, `j = 281`, `m' = 2256` and `m_max(4096) = 2631` (exact). So (R1) holds, and `w >= 10`.
2. **The cheap span.** The three columns `1[x_i = 0]`, `i = 1, 2, 3`, have widths `<= 2j` each. They form a cheap
   family of cost `6j = 1686 <= m' - 1`.
3. **Order and collisions.** Lemma 5.0 gives transversal order `o >= 570`. The span has 8 classes of 512 points, so
   `P = 1,046,528` and `N_2/P = 8.0137`.
4. **Conclusion.** `o N_2/P >= 4567.8 > d = 4095`. Prop 5.1(ii), evaluated at `y` = the 3-bit fraction point
   (`delta_y = 1/7`, `Gamma_L = 30.20`), gives `ln(1/Lambda_A(g)) > (d/2) g + ln((m'+1)/d)` for every `g > g* = 328.4`.

So **MIN\*(A) is false** for this regime box. It fails only at levels `g > 328`, that is `W > 2.7*10^6`, far above
`W_min(4096) = 4.5*10^4`. This is consistent with the author's `n_B = 5185`, but it is stronger: it is an actual failure,
not a failure of a sufficient inequality. Two consequences follow:

* Cor 4.2 needs `n_1 > 4096`.
* The large-`g` part of its hypothesis carries obstructions that the near-extremal regime never sees.

For `k >= 13` the proved orders no longer violate the condition (table in `sylv4col_out.txt`). For `k <= 18` the
4-column spans pass only if their *true* orders are `<= 2^{k-4} - 1` (for example `<= 4095` at `k = 16`). Their proved
lower bounds are `<= 8`, and their upper bounds are `m'`. So between `n_B` and `n_A` the true transversal order of
4-column Sylvester spans decides MIN\*; this is §7.2(c) made concrete.

## 3. Independent computations

All of these are in `referee_nearstrata_checks/`, and none imports the author's code. *Exact* means integer or rational
arithmetic. *Float* means double precision, used only where the margin is large and stated.

| script | content | status | result |
| --- | --- | --- | --- |
| `params.py` | `w, j, ceil(4w), 2jk, m_max` for Sylvester `k = 14..17` | exact `m_max` (big-integer binomials), `w` at 50 digits | window `(2jk, m_max]` empty at `k = 15` (26340 > 21034), nonempty at `k = 16` (40896 < 42065) and at `k = 17` |
| `r1_check.py` | (R1) thresholds | exact `m_max` | last failure 2899 / 3622 / 5808 for `r = n, n ln n, n^2` |
| `scans.py` | Cor 3.7 (both conditions) and Prop 5.3 `n_B` | integer parameters, float ratios | 5185 / 6425 / 10280 (Cor 3.7); 5185 / 6398 / 10280 (`n_B`) |
| (inline) | exact rational recheck of `n_B` near 6398-6425 | exact | `max(oQ) - d = 11.01` at 6397, `-30.03` at 6398, `0` exactly at 6424 |
| `nA_scan.py` | Prop 5.3 `n_A` | float scan to `4*10^6` (`r = n`), `2*10^6` (others); exact rational at the threshold | 473,716 / 564,084 / 947,350; exact margins `-7.6e-5`, `-6.4e-5`, `-6.0e-5` |
| `r_monotone.py` | threshold monotonicity in `r` | float | at `n = 5185`, 40 values `r in [5138, 5177]` fail; none fail at `n >= 5500` (for `r <= n`) |
| `gadget_bruteforce.py` | Thm 2.5(a): exhaustive maximum gadget degree vs `tau` | exact enumeration, `s <= 9` | 1795 instances; 0 violations; equality in all |
| `onecut_tc_exact.py` | Lemma 2.5 / Thm 2.5: exact nullspace of the moment matrix of `A_{S,L}`, membership of strata points in `TC_V` | exact (FLINT `fmpz_mat`, `Fraction`) | 7 boxes (`d = 5, 6, 7`); every sampled stratum point is in `TC_V`; 42 random points are not |
| `sylv_gamma.py` | `Gamma` of the binary-fraction point | 40-digit closed form and direct sum | 2.999621 at `k = 16`; increasing towards 3 |
| `sylv_pairdata.py` | genuine-prime Sylvester `k = 16` profile: Walsh-Hadamard pair distances, `eps_xy`, `abar`, feasibility of the individual gaps | sieve exact; sums in float64 (error `< 10^{-8}`, margins `> 7`) | 65535 primes `= 1 mod 4` in `[e^25, e^25 + 6*10^6]`; `W = 1.638*10^6` (T3-consistent); `abar = 12.50`; `eps in [12.4997, 12.8739]`; `(n-1) e^{-max eps/2} = 104.9 > 2`, so no `u` satisfies the individual gaps |
| `sylv4col.py` | 3- and 4-column Sylvester spans at `m' = ceil(4w)`, `k = 12..20`: cost, proved order, critical order `dP/N_2` | integer | 3-column span violates Prop 5.1 at `k = 12` (`570 > 511`); passes for `k >= 13` |
| `sylv12_minstar.py` | explicit failure level of MIN\* at `k = 12` | float, large margin | MIN\* false for `g > 328.4` |

## 4. Corrected statement

**Setting.**
- **Fixed throughout:** `T = 1`; `w = sqrt(8 d ln 4r)`; `j = floor(w/2)`; `H = U^2`; any profile
  `beta_p in [-1,1]^n`. Each statement says at which order `m'' in [ceil(4w), m_max]` it holds.
- **Pair excesses:** `eps_xy = d_xy - W/2 + 2 ln C'`, and `abar = W/(2(n-1)) + 2 ln C'`.

**Proved (as in the write-up):**

1. **Duality (Prop 2.1).**
2. **Model bounds (Thms 2.4-2.5).** The effective bound
   `kappa_H >= (m!/(2d)^m)(dist_inf(u, TC_H)/(n-1))^m` holds on the uncut ternary cube and on every one-cut box
   (`H >= 2^m`), with the tangent cone an explicit union of rational subspaces. These boxes are zonotopal, and no actual
   regime box is (Remark 2.10).
3. **Qualitative vs quantitative (Props 2.6-2.7).** `Lambda_A(g) > 0` for all `g` iff ZE-TC^braid. That gives CL only
   ineffectively. No bound that depends only on dimension, degree, height and "zero set in `Br`" can give exponent `d/2`.
4. **Separation (Lemmas 3.1-3.3, Theorem 3.4, Cors 3.5-3.6).** For any cluster on any norm:
   - each pair satisfies `|hat u_x - hat u_y| >= e^{-eps_xy/2}`, and `sum eps_xy <= N_2 abar`, so `Gamma(hat u) <= abar`;
   - `dist(hat u, L) >= (1/2) e^{-N_2 abar/(2P(L))}` for every subspace `L` with `P(L)` collision pairs;
   - every width-cheap span at `m' = ceil(4w)` (`w >= 10`) has at most 16 classes outside `<= (m'-1)/2` exceptional
     points.
5. **Assembly (Thm 4.1, Cor 4.2, Thm 4.3).** Proved implications; they hold verbatim at every order in
   `[ceil(4w), m_max]`.

**Strengthened assembly (referee).** Theorem 4.1 holds with (MIN) replaced by the weaker, profile-specific

```text
(MIN^pair)   min{ kappa_H(u) : |u|_inf = 1, |u_x - u_y| >= e^{-eps_xy/2} for all x < y in [n] }  >  C d e^{-W/4}/(m''+1),
```

with the minimum taken as `+infinity` if the set is empty. The proof is the same, using Lemma 3.1 in place of Lemma 3.2.
The data are the profile, the primes and `W`. (MIN) implies (MIN^pair), and the converse fails: the genuine-prime
Sylvester profile at `k = 16` violates (MIN) at every order in `(2jk, m_max]` and satisfies (MIN^pair) vacuously. So:

- MIN^pair, not MIN\*, is the minimal statement relative to the pair inequalities.
- Its qualitative core is "no common real zero of the `Phi_c` meets every individual gap", which is weaker than
  ZE-TC^braid.

**Necessary conditions (corrected wording).**
- If MIN\*(A) holds at *all* levels `g >= 1.5 ln n - 7`, then every linear subspace `L` of `TC_H` with transversal order
  `>= o` has `o N_2/P(L) <= d` (Prop 5.1).
- No such condition is proved for (MIN) at a single level, or for the bounded levels `g <= (1/c + o(1)) ln n` that suffice
  for the growth question.
- MIN\* is false on the regime Sylvester box `n = 4096` (`m' = 2256`) for `g > 328.4`. So any `n_1` in Cor 4.2 exceeds
  4096.
- (MIN) is false on Sylvester boxes at orders in `(2jk, m_max]`, `k >= 16` (Prop 5.2). That is an obstruction to the
  averaged route only.
- At `m' = ceil(4w)`, every width-cheap span meets the MIN\* necessary condition:
  - whatever its true transversal order, for `n >= 473,716` (`r = n`; 564,084 for `r = n ln n`; 947,350 for `r = n^2`);
  - with the proved orders, for `n >= 5,185 / 6,398 / 10,280`.

  These are scans up to `4*10^6` and `2*10^4` respectively, for those three `r` only; the thresholds are not monotone
  in `r`.
- Cor 3.7: for those `r`, `sigma_F > rho_F^min` for `n` in the scanned ranges. This is a necessary compatibility
  condition, since `rho_F^min` is a lower bound for the true minimal radius.

**Sufficient condition (corrected).** CL with `theta < d` implies (MIN) when
`ln(1/gamma) + theta ln C' + ln(Cd/(m''+1)) < (1 - theta/d) W/4`. At `W >= 3n ln n - 14n` this allows
`ln(1/gamma) <= (1 - theta/d)(3/4 - o(1)) n ln n`, which is `(3/8 - o(1)) n ln n` only when `theta <= d/2`.

**Not proved:** any growth improvement or uniform bound, MIN\*, MIN^pair, or ZE-TC^braid for any regime box.
