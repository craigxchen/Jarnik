# Round 21 judge: the low-height window P_+(N') <= M^(2-delta)

Two write-ups this round. An adversarial referee attacked each one, and each report sits next to its write-up:

| key | write-up | referee report | side of the window |
| --- | --- | --- | --- |
| prove | `prove.md` + `prove_checks/` (notes `notes.md`) | `referee_prove.md` + `referee_prove_checks/` | actual clusters: counting, the large sieve, the graft barrier inside the window |
| fakes | `fakes.md` + `fakes_checks/` | `referee_fakes.md` + `referee_fakes_checks/` | fake circles: Haar scarcity, rigidity of complete prime segments, transplants, large sieve |

Paths are relative to this directory (`round21/`, scratchpad). `overview` means
`docs/claude_review/overview/overview.md` in the repository, which was only read; the edited copy is `round21/overview.md`.
The labels are those of the overview: proved / conditional / model / exact / evidence / heuristic / unrefereed.
Statements in §1 are the **referee-corrected** ones. Remarks added by this README are marked *judge* and are
**unrefereed**; the judge's computations are in `judge_checks/` (§4).

**Notation.** `N'` is the normalised norm, `M` the number of points, `W = log N'`, `Y = P_+(N')`, and "the window"
means `Y <= M^(2-delta)` for a fixed `delta > 0`. `s` (prove) and `r` (fakes) both denote the number of primes. The
free dimension is `m = s - rank B` (prove's `B` has rows `2(a_x - a_1)`; fakes' `L` is the row lattice of the same
matrix, so `m = r - dim L_R`). A profile is *pinned* or *rigid* when `m = 0`.

## 0. Bottom line

**No growth improvement and no uniform bound were proved this round,** either unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` is still the best bound for actual circles, inside the window as
everywhere else.

**Status of the low-height window `P_+(N') <= M^(2-delta)`: OPEN.** There is no restricted growth theorem inside it,
and it is not capped by fakes. Both write-ups survive refereeing as partial results. Every formal statement is
correct; five major overstatements were corrected (§1).

What the round established (refereed unless marked):
1. **Restricted growth below the window** (prove Cor 2.2; proved; minimal novelty). If `P_+(N') = o(log R)`, then
   `M = o(log R/loglog R)`, uniformly in the arc's position. Near-extremal clusters have `P_+(N') >= c_eps M log M`, so
   this misses them. **Restricted progress only.**
2. **The large sieve in the window reduces to the cross-gap item.**
   * Lemma 3.2: the large sieve is equivalent to arc counts up to an absolute factor.
   * Cor 4.3: the realisation condition alone never forces concentration for unpinned profiles.
   * Theorem LS: every angle model satisfies the cluster-only sieve at the genuine constant and the union sieve at
     twice it.
   * The union sieve at the genuine constant is **uncapped** in the window (row 7 below `2D`).
3. **The record's fake constructions cannot be run in the window.**
   * Theorem A: the Haar success probability on any realisation coset is `<= exp(-2.17 sum_J p_j^(-1/2))`, which is
     `<= exp(-c M^(delta/2))` when `m >= cM`.
   * Cor A.1: the `r^2` law is intrinsic to Haar completion.
   * J1 (*judge, unrefereed, exact check*): equal-weight windows need `P >= 16 r(r-1) e^(-1/(4r))` by counting alone.
   * None of this proves that window fakes do not exist.
4. **Rigidity and transplants.**
   * Theorem C: angle models are genuine on every complete initial segment of split primes.
   * Lemma T: arithmetic models are actual clusters.
   * The window's `A_0` question becomes (Q): do non-arithmetic gap-admissible points survive on realisation cosets
     as `Lambda -> infinity`? Neither answer gives growth (J2, *judge*). "Yes" caps only `A_0` and the row-3 profile
     inputs.
5. **New uncapped items inside the window** (found by the referees):
   * square-root-quality (GRH) prime sums and box counts;
   * the union large sieve at the genuine constant;
   * robust inputs on primes in `[c M log M, M^(1+eta)]`;
   * the pinned case `m = 0`, in which all four record clusters have `m <= 2` (J3).

**Scores** (rubric of `round10/README.md` §1; §6): fakes 4, prove 3.

**Recommendation** (§7): **stop.** The window is the general problem with smaller primes, not a cheaper sub-problem.
The only tractable decisive item left there is to extend Proposition A' to realisation cosets, which would cap the
window for `A_0`. That has map value only (about 10-15% in one round). A restricted growth theorem inside the window
in one more round: about 2%. Unconditional growth: < 0.1%.

## 1. Referee-corrected statements
### 1.1 prove (`prove.md`; `referee_prove.md`): the partial result survives, with its scope corrected

The referee re-derived every proof and recomputed every number with independent code. No formal statement is false.
Three narrative claims were overstated (two of them major), and one remark was labelled proved when it is a heuristic.

**Proved (referee-corrected).**
1. **Thm 2.1, Cor 2.2, Cor 2.3: restricted growth below the window.** These are corollaries of `round8/stability.md`
   Thm 3.3(3) combined with `omega(N') <= pi(Y; 4, 1)`. Every cluster with `M >= 5` satisfies
   `M < 4 omega(N')(1 + 2M log^+ C/W') <= 4 pi(P_+(N'); 4, 1)(1 + 2M log^+ C/W')`. Hence, for fixed `C` and
   `M >= M_1(C)`, `M < 12.56 max((log R)^(1/2), P_+(N'))/loglog R`. So **clusters with `P_+(N') = o(log R)` satisfy
   `M = o(log R/loglog R)`**, uniformly in the position of the arc, with thresholds depending only on `C`.
   Near-extremal clusters (`M >= (eps/12.56) log R/loglog R`) therefore have `P_+(N') > (eps/12.56) log R >= c_eps M log M`.
   Novelty is minimal: for `C <= 1` the record's Prop 3.4 already gives `M <= (2 + o(1)) pi(Y; 4, 1)`.
   Checked on all 724 windows of the record circles (referee) and all 267 windows of the seven record norms (author).
2. **Lemma 3.2 (large sieve = arc counts, up to constants).** For any finite set of angles, the optimal dual
   large-sieve constant at integer frequencies `|k| <= K` lies between `(2K + 1) m''/sqrt 2` and `32 K m''`, where
   `m''` is the maximal number of points in an arc of length `1/(8K)` (unit circle normalised to length 1). For a
   cluster's prime angles the arc counts are at most the lattice count `eta Y/2 + 2.53 Y^(1/2) + 1.58` (Lemma 3.1).
   In the window this bound is non-trivial by a factor `>> M^(delta/2)` (Prop 3.3). Nothing is said about
   sharp-constant large sieves.
3. **Lemma 4.1, Prop 4.2, Cor 4.3.** `Real(Delta)` is invariant under `H = ker(B mod pi/2)`. If all coordinates are
   pairwise independent in `Cok = Z^s/B^T Z^(M-1)` (which needs `rank B <= s - 2`), every non-empty `Real(Delta)`
   contains points whose large-sieve constants and arc counts are within a factor 3 of the genuine prime-angle
   bounds. So a large sieve with constants relaxed by a fixed factor, **together with the realisation condition
   alone**, cannot exclude such profiles. **Nothing is proved for pinned or partly pinned profiles, or for the large
   sieve combined with the monomial gaps.**
4. **Props 5.2-5.3 and Cor 5.4 (conditional).** Grafts whose angles are spread, and whose primes' genuine angles are
   also spread, satisfy every `(c, gamma)`-robust prime-supported large sieve up to a factor
   `1 + O(M^(-eta/4)(log M)^(1/4))`, even when `s >> X^(1/2)`. They satisfy robust prime sums and box counts **only for
   `gamma > 1/(1 + eta)`**, which is `> 1/2` inside the window. Hence: **if** spread gap-admissible angle models exist
   with `W = O(M log M)` and all primes in `[M^(1+eta), M^(2-delta)]`, **and** the genuine angles of those primes are
   spread (an unverified upper-bound-sieve assumption, Caveat 3), then `A_0` plus these restricted robust inputs cannot
   prove growth for such clusters. **Square-root-quality inputs (GRH-quality Hecke prime sums, box counts with error
   `x^(1/2)`) are not covered anywhere in the window, and neither are primes in `[c M log M, M^(1+eta)]`.**
5. **Prop 6.1 (first moment).** Haar-random angles violate on average `>= 1.08 (exp(2 sum_j (p_j^(1/2) + 1)^(-1)) - 1)`
   monomial gaps, and `Prob(all single-prime gaps hold) <= exp(-2.17 sum_j p_j^(-1/2))`. So first-moment (union-bound)
   constructions on the free coordinates need `sum_j p_j^(-1/2) = O(1)`. All refereed fakes have
   `sum_j p_j^(-1/2) <= 0.22` (corrected from "`<= 1/32`"; `round5/construct.md` Thm 2 reaches `21^(-1/2)`). A window
   object with `s >= cM` primes has sum `>> M^(delta/2)`. **`s >= cM` is proved for actual clusters; for `A_0` fakes it
   needs a transfer of `stability.md` Thm 3.3 that is not proved** (probably repairable from the pair chords).
6. **Prop 6.2** (= `round10/analytic.md` Thm D): lattice-angle realisations are shadows of actual window clusters.
   **Prop 8.1:** for `A_0` alone, growth in the window is equivalent to `W_min^win(M)/(M log M) -> infinity`. This
   restates the definition of a fake. With analytic inputs, only the conditional half (item 4) is available.

**Heuristic, not proved:** Remark 2.4 ("counting provably cannot pass exponent 1"; it was labelled proved); Remark 4.4
("in the pinned case the gaps already give more than the large sieve": a comparison of magnitudes, and below a
profile-dependent threshold `c(B)`); "the large sieve plays no part / a window theorem cannot come from the large
sieve" (§8.1, §8.4); Heuristic 6.3 ("window fakes plausibly exist").

**Evidence (reproduced by the referee).**
* The exhaustive scan of the 40,715,511 split-only 200-smooth norms `<= 10^18`.
* The best `C_k` per decade beyond `10^9` trends upward, **but not monotonically** (e.g. `k = 7`: 7.12 then 6.61), and
  never beats the records for `k = 7, 9, 10`.
* **No record cluster lies in the window for any `delta > 0`:** every one has `P_+ > k^2`, with ratios 1.13 to 3.64.
  Restricting to `P_+ <= k^(7/4)` worsens the best constants by factors up to about 12.
* The record data: Smith invariants, `D_j = (24, 8, 6, 4, 24)` for the 7- and 10-point records, and rank deficits.
* Gap Monte Carlo on the **full torus** only: measure about `exp(-(1.1 to 1.3) E)` for `E <= 10`.
### 1.2 fakes (`fakes.md`; `referee_fakes.md`): the partial result survives; three stated consequences overreach

The referee re-derived every numbered proof and recomputed every number with code written from scratch. All proofs
are correct. Three stated consequences are overstatements (all major), and the evidence scope is narrower than claimed.

**Proved (referee-corrected).**
1. **Theorem A (Haar scarcity on realisation cosets).** Let `J` be any set of coordinates whose unit vectors are
   linearly independent modulo `L_R`; one can take `|J| = m`. On every coset of `H_0` (every non-empty `Real(delta)`
   is a finite union of such cosets), the Haar proportion of gap-admissible vectors is
   `<= prod_{j in J} mu_1(G_1(p_j)) <= exp(-2.1735 sum_{j in J} p_j^(-1/2))`, with `kappa = (4/pi)(1 + 2^(-1/2))`.
   This is the coset form of prove's Prop 6.1(3) and has the same constant.
2. **Corollary A.1.** Haar completion with success probability `>= 1/2` needs `P_+ > 9.83 m^2`. For the skew double
   core (`m = M - 1`) this is `9.83 (M - 1)^2`.
3. **Corollary A.2.** In the window, **if `m >= cM`**, the Haar success probability is `<= exp(-2.17 c M^(delta/2))`. So
   the record's union-bound existence proofs cannot be run there. **This is not a proof that no Haar point succeeds**
   (corrected from "the record's constructions provably fail"). The capped inputs force `r > M/4`, not `m >= cM`, so
   `m >= cM` is a hypothesis.
4. **Proposition A'.** On the **full torus** (no profile, no coset), `mu(G) > 0` up to `Lambda ~ (1/4) log P`. This is
   a `log^2` gain over the union bound, still far from the window.
5. **Proposition B, Theorem C, Corollary C (rigidity of complete segments).** If every split prime below `p` is among
   the model's primes, with genuine angles, the gap inequalities force `phi_p = ± arg pi_p (mod pi/2)`. Hence an angle
   model is genuine on every complete initial segment of split primes. The Step-1 covering is branch-independent, and
   the referee checked exactly that it leaves only `± arg pi_p` for all 180 split `p <= 2477` (the author's 63-branch
   certificate to 53 is redundant).
6. **Lemma T, Corollary T (transplants).** Every gap-admissible model with `phi_j = arg beta_j`,
   `Norm beta_j <= p_j`, that realises a profile is an actual cluster on a circle of norm `<= N`, with points distinct
   up to units. Coprimality is not needed. So barrier-scale such models for unbounded `M` would disprove
   Cilleruelo-Granville.
7. **Theorem LS (large sieve).** Every angle model satisfies the cluster-only Hecke large sieve over split primes at
   `D = 2K + pi X/2` (`K` an integer). The union with the external split-prime angles satisfies it at `2D`.
   **The union inequality at the genuine constant `D`, or at any constant below `2D`, is not shown for window fakes.**
   At `K <= X^(1/2)` it is a row-7 robust (R1) input. Theorem B caps such inputs only for grafted primes `>= M^A` with
   `A > 2`, so they **remain uncovered in the window**. The window's large-sieve question is reduced to the cross-gap
   item, not eliminated (corrected from "the large sieve adds nothing").

**Corrected consequence for (Q).** "(Q) yes => the window is capped" holds only for **class `A_0` and the row-3
profile inputs**. Rows 2 and 4 are not capped: their residue designs need all primes `> X = M^A` with `A >= 2`, which
is impossible in the window. Row 7 is not capped either (item 7), nor are Haar-a.e. inputs on `T^r`, which may fail
on a whole null coset. "Arithmetic" has to mean `phi_j = arg beta_j` with `Norm beta_j <= p_j`.

**Heuristic, not proved:** the grouping heuristic ("the `r^2` law comes from Haar completion, not from the
equal-weight window"; stated flatly in §0 and the summary); the Poisson law and its extrapolation; "window fakes
cannot come from a positive-probability random choice". The remark that transplants satisfy the cross-gap is
**false** whenever some `beta_j` is associate to an external prime of norm `<= X`; for those, `sin = 0`.

**Evidence (reproduced by the referee).** All of it is on the **full torus** `T^r`, with `r <= 12` primes forming
complete segments minus 5 and `Lambda <= 1.8`. None of it is on realisation cosets or sparse window sets.
* Positive measures, consistent with a Poisson-type decay `exp(-c(Pi - 1))` with `c` between 1.1 and 1.8 (9 points,
  `Lambda <= 1.03`).
* Non-arithmetic gap-admissible vectors, verified against all monomials of weight `<= 34` at `Pi = 17.4` (360,050)
  and of weight `<= 33` at `Pi = 26.4` (1,043,472). The tails beyond the truncation are 0.022 and 0.147. The sentences
  "points ... **exist**" (§5.5) and "rigidity is **false**" (§5.4) are evidence phrased as fact.

## 2. Status of the low-height window

**Status: OPEN.** No restricted growth theorem has been proved inside the window, and the window is not capped by
fakes. Precisely:

| part of the problem | status after round 21 | where |
| --- | --- | --- |
| clusters with `P_+(N') = o(log R)` (**below** the window) | **restricted growth proved**: `M = o(log R/loglog R)`, uniformly in the arc's position (minimal novelty) | prove Thm 2.1, Cor 2.2 |
| near-extremal clusters | must have `P_+(N') > (eps/12.56) log R >= c_eps M log M`; in units of `M` the window is `[c M log M, M^(2-delta)]` | prove Cor 2.3 |
| window, inputs of class `A_0` (pairs, monomial gaps) | **open** = question (Q); union-bound fake constructions cannot run there when `m >= cM` (Thm A, Cor A.2); single equal-weight windows are impossible (J1); arithmetic models are actual clusters (Lemma T); complete initial segments are rigid (Thm C) | fakes §3, J1 |
| window, robust large sieve with relaxed constants | not capped. Conditionally capped for spread models with primes `>= M^(1+eta)` (Cor 5.4). Cluster-only LS at `D` and union LS at `2D` follow from the gaps (Thm LS) | prove §5, fakes §4 |
| window, union large sieve at the genuine constant (row 7 below `2D`) | **uncapped** (new) | referee_fakes §5 |
| window, square-root-quality prime sums and box counts (`gamma <= 1/(1+eta)`, GRH quality) | **uncapped** (new) | referee_prove §4 |
| window, primes in `[c M log M, M^(1+eta)]` for robust analytic inputs | **uncapped** (new) | referee_prove §4 |
| window, pinned/rigid profiles (`m = 0`) | **untouched** by every measure statement; profile-dependent concentration forced | prove Rem 4.4, fakes open 5, J3 |
| window, residue inputs (P3-all, rows 2 and 4) | **uncapped**: the residue designs of rows 2 and 4 need primes `> M^A`, impossible in the window | referee_fakes §0 item 3 |

**Is the window populated? (evidence)** At the sizes reachable (`M <= 11`), no record cluster has `P_+ <= k^2`, so none
lies in the window for any `delta > 0`. Small clusters with `P_+ <= k^2` exist, but their constants are worse. This
says nothing about the asymptotic window.

**Has the window's status changed?** Its *label* has not: it was "no fakes known; open", and it is "no fakes known;
open". Its *content* has changed in four ways:
1. The record's explanation ("every known fake needs `P >~ r^2 log r`") is now proved twice over. Equal-weight windows
   need `P >= 16 r(r-1) e^(-1/(4r))` by counting (J1). Haar/union-bound completion needs `P_+ > 9.83 m^2` (Cor A.1).
2. The obstruction has been localised to (Q) for `A_0`, plus the uncapped items in the table.
3. The large sieve has been reduced to the cross-gap item. Inside the window, it is no longer a separate unclassified
   input.
4. Three robust analytic sub-items that Theorem B caps outside the window are now known to be uncapped inside it.

## 3. Additions to the barrier map

These are proposals for overview §5.4/§5.5. The overview copy in this directory implements the §5.5 change (§8).
1. **§5.5 "low-height window" item** (replaces the round-10 wording; still uncapped):
   * Equal-weight window designs are impossible there (J1, *judge, unrefereed*), and union-bound/Haar completion has
     success probability `<= exp(-2.17 c M^(delta/2))` when `m >= cM` (fakes Thm A; proved, refereed).
   * The cluster-only large sieve at `D` and the union large sieve at `2D` follow from the monomial gaps (fakes Thm LS;
     proved, refereed). The union large sieve at the genuine constant is the cross-gap item and is uncapped.
   * Spread grafts with primes `>= M^(1+eta)` satisfy every robust large sieve, and robust prime sums and box counts
     with `gamma > 1/(1+eta)`. This is conditional on the existence of spread window models (prove Cor 5.4). It does
     not cover GRH-quality prime sums and box counts, or primes below `M^(1+eta)`.
   * Arithmetic angle models (`arg beta_j`, `Norm beta_j <= p_j`) are actual clusters (Lemma T). Complete initial
     segments of split primes are rigid (Thm C).
   * **Open: (Q)**, non-arithmetic gap-admissible points on realisation cosets as `Lambda -> infinity`. A yes caps
     only `A_0` and the profile inputs. A no makes `A_0` lossless in the window and does not give growth (J2).
   * The pinned case `m = 0` is untouched.
2. **Restricted growth below the window** (a restricted result, not progress on the rate):
   `P_+(N') = o(log R) => M = o(log R/loglog R)` (prove Cor 2.2; proved, refereed; a corollary of `stability.md` Thm 3.3).
3. **Refuted/corrected claims** for §8.2:
   * "The large sieve adds nothing in the window" (fakes);
   * "The short-arc condition never forces concentration" (prove (F2): unpinned profiles, realisation alone);
   * "Counting provably cannot pass exponent 1" (prove Rem 2.4: heuristic);
   * "All refereed fakes have `sum p^(-1/2) <= 1/32`" (the correct bound is `0.22`);
   * "(Q) yes caps the window" (only `A_0` and row 3).

## 4. Judge's own checks (unrefereed)

### J1. A single equal-weight window is impossible in the low-height window (proof + exact check)

*Statement (judge; unrefereed).* Let `r >= 2` and `a > 0`. The interval `I = [a, a e^(1/(4r))]` contains at most
`a e^(1/(4r))/(16 r) + 1` integers `= 1 mod 4`. Hence, if `r` split primes all lie in one window `[a, a e^(1/(4r))]`,
then `a >= 16 r (r - 1) e^(-1/(4r))`. In particular, a profile whose `r >= M/4` primes lie in one such window has
`P_+ >= a >= M^2 (1 - 4/M) e^(-1/M)`, which exceeds `M^(2-delta)` once `M >= M_0(delta)`.

*Proof.* `|I| = a(e^t - 1)` with `t = 1/(4r)`, and `e^t - 1 <= t e^t`, so `|I| <= a e^t/(4r)`. A closed interval of
length `y` contains at most `y/4 + 1` integers in one residue class mod 4. If `I` holds `r` such integers then
`r <= a e^t/(16 r) + 1`. For `r >= M/4` the bound reads `16 r(r-1) >= M(M - 4)`. QED

This turns referee_prove's counting remark ("a window holds `~ P/(8 r log P)` split primes") into an elementary
inequality: the equal-weight window design of every refereed fake needs `P >~ r^2` for trivial counting reasons.
Brun-Titchmarsh in the Montgomery-Vaughan form adds the factor `log`. It does **not** decide whether other designs
work in the window; that is (Q).

*Exact check* (`judge_checks/j1_window.py` -> `j1_out.txt`, re-verified with 50-digit endpoints by
`j1_verify.py` -> `j1_verify_out.txt`). For each `r`, the table gives the least split prime `q` such that
`[q, q e^(1/(4r))]` holds `>= r` split primes. Since the window can be shifted to start at its smallest prime, `q`
bounds below the smallest prime of any `r` split primes in one window.

| `r` | 2 | 5 | 10 | 20 | 50 | 100 | 200 | 400 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| least `q` | 37 | 977 | 5,189 | 23,537 | 175,753 | 906,289 | 4,368,113 | 19,950,661 |
| `16r(r-1)e^(-1/(4r))` | 28.2 | 304 | 1,404 | 6,005 | 39,005 | 158,005 | 636,005 | 2,552,005 |
| `q/(r^2 log r)` | 13.3 | 24.3 | 22.5 | 19.6 | 18.0 | 19.7 | 20.6 | 20.8 |

The lower bound holds in every case, and `q ≈ 20 r^2 log r` throughout (consistent with the PNT heuristic `16 r^2 log r`).

### J2. Neither branch of (Q) gives a growth theorem in the window (logical remark)

*Judge; unrefereed.* fakes §0 and §6.3 say that under scenario R ("(Q) no"), "no window fakes" would be a restricted
growth theorem. That is true as an equivalence, but it does not follow from R. Under R every barrier-scale window angle
model is arithmetic, hence an actual cluster (Lemma T). So for `A_0`, "no window fakes at `W = O(M log M)`" becomes
**the same statement as** "no actual window clusters at `W = O(M log M)`", which is the restricted growth statement
itself. Put differently, R makes the `A_0` relaxation lossless in the window; it does not bound `W_min^win`.
Under scenario P the window is capped for `A_0` (§1.2). Hence:
* deciding (Q) cannot by itself produce a window growth theorem;
* "yes" closes the window for the `A_0` method; "no" shows that `A_0` loses nothing there, so whatever proves growth
  must use the actual arithmetic of `N'` beyond the pair and monomial gaps, or the uncapped analytic inputs of §3.

prove's Prop 8.1 already says this for `A_0` (growth in the window is equivalent to `W_min^win/(M log M) -> infinity`);
J2 only records that (Q) does not decide `W_min^win`.

### J3. The pinned/rigid case is the common blind spot, and the records sit in it (read from refereed data)

*Judge; unrefereed reading of exact data that the referee reproduced* (`referee_prove_checks/records_out.txt`). The
free dimensions `m = s - rank B` of the four record clusters are:
* 10-point and 7-point records: `rank B = s = 5`, so `m = 0` (pinned);
* golden 8-point: `rank B = s - 1`, so `m = 1`;
* 5-point: `rank B = s - 2`, so `m = 2`.

So on every record cluster:
* Theorem A's bound `exp(-2.1736 sum_{j in J} p_j^(-1/2))` has `|J| <= m <= 2`. Its value is exactly 1 for the two
  pinned records, at least `exp(-2.1736 * 5^(-1/2)) = 0.378` for the golden 8-point cluster, and at least
  `exp(-2.1736 (5^(-1/2) + 13^(-1/2))) = 0.207` for the 5-point cluster (computed inline).
* Cor 4.3 applies only to the 5-point cluster.

Both write-ups list the pinned/rigid case as open: prove open point 3 (with the profile-dependent threshold `c(B)`)
and fakes open point 5 (`m = 0`, "a discrete problem about characters of `Z^r/L`"). Every measure-theoretic statement
of the round is vacuous there. The pinned case is also where prove's Remark 4.4 says concentration onto
`(pi/2)Z/D_j` **is** forced, below a threshold that depends on the profile. That threshold may not be treated as
uniform in `R` (owner constraint).

## 5. What this means for growth and the uniform bound

**Growth: unchanged.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound valid for all
`R`, including inside the window. No rate `o(log R/loglog R)` follows from anything proved this round or from any
standard conjecture.

**Restricted progress (clearly restricted).** One restricted growth theorem was proved: clusters with
`P_+(N') = o(log R)` have `M = o(log R/loglog R)` (prove Cor 2.2). It applies **below** the window, to a class that
near-extremal clusters provably avoid (Cor 2.3), and it is a two-line corollary of `stability.md` Thm 3.3. So it does
not bear on the rate. No restricted growth theorem was proved **inside** the window.

**Uniform bound: unchanged.** It is still known only conditionally: on `V(7, 1/4)`, `V_8`/`(L_8)`, or `V_E(5)`/`(L_5^circ)`.
Nothing this round touches the conditional lane.

**What the round did establish.** The low-height window is not a soft spot of the barrier map that an input could
exploit cheaply:
* The large sieve in the window reduces to the cross-gap item (fakes Thm LS; prove Lemma 3.2).
* Counting alone gives restricted growth only up to `P_+ = o(log R)`. The window starts at `c M log M`.
* The `A_0` part reduces to (Q), and by J2 neither answer to (Q) gives growth.
* What remains uncapped in the window is the same as everywhere else: residue inputs (P3-all), sharp-constant and
  genuine-constant large sieves, and square-root-quality analytic inputs. To these the window adds the pinned case.

The one genuinely new structural fact is the rigidity of complete initial segments (Thm C): **models are forced to be
genuine on any complete initial segment of split primes they contain.** It is proved and exactly checked. Its
relevance to growth depends on a sparse version (fakes open point 1). Near the 2/3 bound, actual clusters contain
almost complete initial segments (`stability.md` Thm 5.2, as the referee notes). But a growth improvement must handle
general ratios `A`, where clusters are sparse.

**Rubric (b) did not move.** No write-up produced an input that is true for actual circles, provably false for every
fake of size `O(M log M)`, and weaker than the problem. Theorem C's hypothesis (complete segments) is satisfied by
actual clusters only near the extremal ratio, and fakes need not meet it.

## 6. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing; (b) fake-violated input; (c) movement of growth or
uniformity. Each sub-score is out of 10. The overall score is a judgement that weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **fakes** | 7 | 2 | 2 | **4** | The stronger write-up. It proves real results: Theorem A on realisation cosets with `|J| = m`; Cor A.1, which shows the `r^2` law is intrinsic to Haar completion; Proposition A'; Theorem C, a clean rigidity theorem whose Step 1 the referee checked exactly to `p = 2477`; Lemma T in the strong form; and Theorem LS. All of them were re-derived without error. The evidence is careful and reproducible. Deductions: three major overstatements ("provably fail", "the large sieve adds nothing", "(Q) yes caps the window"); evidence on the full torus and on complete segments presented as window-like; evidence phrased as fact; a false cross-gap remark. Also, as J2 notes, scenario R is presented as giving a growth theorem when it does not. (b) = 2: Theorem C is an input that actual clusters satisfy and generic models violate, but only on complete segments, which fakes need not contain. (c) = 2: a proved correction of the barrier map (the large sieve reduced to the cross-gap item, Haar scarcity) and a sharp formulation (Q). |
| **prove** | 6 | 1 | 2 | **3** | A careful, honest attempt that correctly claims no window growth. Every formal statement is true, and the referee reproduced every computation exactly, including the scan of `4.07 * 10^7` norms. Lemma 3.2 (the large sieve equals arc counts up to an absolute factor) and Props 5.2-5.3 (Theorem B extended into the window, conditionally) are useful map results. Deductions: two major narrative overstatements, (F2) and (F3); a heuristic labelled proved (Rem 2.4); an unproved transfer to `A_0` (Prop 6.1(4)); a false side claim (`1/32`); misreadings of its own table. Thm 2.1 has minimal novelty. (b) = 1: no fake-violated input. (c) = 2: a restricted growth theorem below the window, plus the referee-extracted identification of three uncapped analytic sub-items inside the window. |

*Referee reports.* Both were decisive and careful. Each re-derived every proof and recomputed every number with
independent code. Between them they found every overstatement, and they extracted the round's most useful map
content: the uncapped GRH-quality and genuine-constant items inside the window, and the exact check of Theorem C's
covering to 2477, well beyond the author's certificate.

## 7. Single most promising next step, and whether to stop

**Recommendation: stop internal rounds on growth and uniformity.** Rounds 19 and 20 already recommended this. Round 21
tested the one remaining item that looked like a soft spot, the low-height window: the only regime with no known fake,
and the only one where the large sieve is non-trivial for the cluster. The result:
* the large sieve there is equivalent to arc counts and reduces to the cross-gap item;
* counting gives restricted growth only below the window;
* the `A_0` part reduces to (Q), and neither answer to (Q) yields growth (J2);
* what remains open in the window (residues, genuine-constant and square-root-quality analytic inputs, the pinned
  case) is what remains open everywhere.

The window is therefore not a cheaper sub-problem. It is the general problem with smaller primes.

**The single most promising next step (map-only, and only if the window is to be settled for its own sake).** Prove
fakes open point 2, **positivity of `nu(G ∩ coset)` on realisation cosets by integrating one coordinate of `H_0` at a
time**, for barrier-scale window profiles with `m >= cM`, as `Lambda -> infinity`. This is Proposition A' carried from
the full torus to `H_0`. As in Proposition A', integrate the free coordinates in order and bound, at each step, the
measure of the bad set of the monomials that first involve that coordinate. The new difficulty is that the remaining
`r - m` coordinates are affine functions of the free ones on the coset, so monomials in them are not independent
events. (*Judge's sketch; unrefereed; not attempted.*)
* *Success* proves scenario P. Window fakes then exist for `A_0` and the row-3 profile inputs, and the overview's
  window item becomes "capped for `A_0`". This is not growth, and it leaves rows 2, 4 and 7 open.
* *A proof of the opposite* (scenario R, via a sparse Theorem C) would make `A_0` lossless in the window. By J2 that
  is not growth either.

So the step has map value only. It is the most tractable decisive item left in the window: Proposition A' already
reaches `Lambda ~ (1/4) log P`, and the only new difficulty is the deterministic coordinates of the coset.

**Honest probabilities.**
* Unconditional growth `M = o(log R/loglog R)` from any further internal round: **< 0.1%** (unchanged).
* An unconditional uniform bound: **essentially 0**.
* A restricted growth theorem **inside** the window (for all window clusters with `P_+(N') <= M^(2-delta)`) from one
  more internal round: **about 2%**. It would need an input that uses the actual arithmetic of `N'` beyond pairs and
  gaps (J2); none is in sight.
* Proposition A' extended to cosets with `Lambda -> infinity` (window capped for `A_0`): **about 10-15%** in one
  focused round. Value: map only.
* A sparse Theorem C strong enough to force scenario R: **about 3%**. Value: map only (J2).

**Decision: stop.** The record's lasting value is unchanged:
* the `(2/3 + eps)` theorem;
* the conditional theorems (`V(7, 1/4)`, `(L_8)`, `V_E(5)`/`(L_5^circ)`);
* the refereed barrier map, now with a referee-corrected low-height-window entry.

The external package of `round19/README.md` §6 and `round20/README.md` §7 needs no new item from this round. If
anything is added, it should be one sentence: "the low-height window is not a soft spot: the large sieve there reduces
to the cross-gap item, and the `A_0` relaxation reduces to the existence question (Q)".

## 8. Files

* Write-ups: `prove.md` (`prove_checks/`, notes `notes.md`), `fakes.md` (`fakes_checks/`).
* Referee reports: `referee_prove.md` (`referee_prove_checks/`), `referee_fakes.md` (`referee_fakes_checks/`).
* Judge: this README; `judge_notes.md`; `judge_checks/j1_window.py` -> `j1_out.txt`; `judge_checks/j1_verify.py`
  -> `j1_verify_out.txt`. J3 reads `referee_prove_checks/records_out.txt`.
* Overview copy: `overview.md`. The round-21 placeholder in the abstract has been replaced by a referee-corrected
  summary. The §5.5 and §7.5 "low-height window" items have been updated, and a round-21 paragraph has been added to
  §8.2. The repository is unchanged.
