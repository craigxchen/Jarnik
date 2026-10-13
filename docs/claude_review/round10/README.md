# Round 10: L4(K) exact test and the fresh-mechanism panel, after refereeing

Four write-ups this round, each attacked by an adversarial referee whose report sits next to it
(`referee_*.md`):
* `L4K.md`: the scale-dependent inverse step L4(K) that `round9/README.md` §6 recommended;
* `analytic.md`: Hecke characters and Gaussian primes in sectors;
* `additive.md`: the relation lattice and additive gaps;
* `determinant.md`: auxiliary polynomials with forced divisibility in the growth regime.

**§§2-5 give the referee-corrected statements.** Remarks added by this README are marked *README remark*. Nobody
but the judge has checked them, so they are unrefereed. Conventions and labels are those of `round7/README.md` and
`round9/README.md`:
* `N = R^2`, `W = log N`, and `n = M` is the number of points;
* near-extremal means `W <= c_0 M log M`;
* regime (1.4) means `p_j <= M^(h_0)` with `h_0 = 2c_0 - 3 + o(1)` and `X = 3M^A`, with `A` chosen after `c_0`;
* the labels are proved / conditional / exact / evidence / heuristic / unrefereed.

## 0. Bottom line

**No growth improvement for actual circles and no uniform bound was proved this round**, either unconditionally or
from a standard conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound.
**No new cap on P3-all\* was proved.** Its bracket `3 M log M - O(M) <= W_min^all <= (66A + o(1)) M log^2 M/loglog M`
is unchanged.

What changed (all refereed):
1. **Route 3 as formulated is refuted** (`L4K.md` Thm 2). L4(K) of `round9/README.md` §6 is false for every
   `K >= 2`. Bounded norms `<= K` force collisions only up to `M^(floor(K/2)+o(1))` (Thm 1), which sharpens round-9
   Prop. H. The surviving question L4\*(K) is a pure label lemma. The referee proved that every fake circle whose
   primes lie in the regime, given the genuine labels of those primes, satisfies its hypotheses *and* conclusion. So
   route 3's fake-separating power, if any, must come entirely from the label-angle step (L7 / slab conditions).
2. **Analytic number theory, used robustly, is capped** (`analytic.md` Thm B, the graft barrier). This covers
   prime-supported inputs with power-saving error terms or large-sieve constants of size `c(X + K X^(1/2))`, with
   `gamma`, `c` and `A` fixed before `c_0`. Used pointwise, it is the lattice itself (Thm D: LCG holds iff each angle
   is the argument of a Gaussian integer of norm `<= p_j`). Two inputs remain unclassified: arithmetic fibre inputs,
   which are P3-all, and the sharp large sieve.
3. **Additive data, used in full, are the problem itself** (`additive.md` §2). The additive gap `AG(M^A)` is violated
   by almost every fake, at the corrected height `M^(1 + c_0/4 + o(1))`. But its universal form, Lemma (AG), is
   strictly stronger than growth. The referee rejected the headline for this reason.
4. **Forced divisibility is not "exactly critical" in the growth regime** (`determinant.md` Thm 2.1). This corrects
   `round8/README.md` §2.4 and the summary of `threshold/`. Once the auxiliary form is chosen after the profile, the
   LP is supercritical by a factor `>= 0.12 sqrt((n-1)/ln 4r)`. Actual clusters in the regime then satisfy an exact
   algebraic system (Thm 4.1). Almost every fake violates it. Under a zero-estimate hypothesis ZE, which is
   non-standard and unproved, clusters on squarefree norms satisfy `M = O(loglog R)` (Thm 5.7).

**Scores** (§6): determinant 5, L4K 3, analytic 3, additive 2.

**Recommended next step** (§8): decide ZE for the slab-cut forced-divisibility linear system. The question is
whether its base locus reaches nondegenerate points near the identity of the torus, i.e. whether it has components
there other than subtori through the identity (non-toric pieces, or translates `gH`). A finite first test is
specified.

## 1. Scoring rubric

Each mechanism gets three sub-scores and an overall score, all out of 10. The overall score is a judgement, not an
average. It weighs (b) and (c) above (a), because a correct result that leaves the barrier untouched does not move
the problem.
* **(a) Correctness after refereeing.** 10 means every claim survived as stated. Deductions for major corrections,
  for claims labelled proved that were not, and for headlines that did not survive.
* **(b) Fake-violated input.** This asks whether the mechanism uses an input that is true for actual circles and
  provably false for the fakes.
  * 0: no such input. The input holds for every fake.
  * 3: violated only by almost every fake, or violated by every fake only because it is equivalent to the original
    problem.
  * 5: violated by almost every fake, with a residual statement that is no longer arithmetic.
  * 10: provably violated by every fake of size `O(M log M)`, by an input weaker than the problem.
* **(c) Movement of growth or uniformity.** 0 means none. Partial credit goes to a conditional growth theorem, a
  proved correction to the barrier map, or a refutation that redirects a route. 10 means unconditional
  `M = o(log R/loglog R)`.

## 2. L4K: the inverse step L4(K) (refereed: refutation survives, corrected question redone)

Write-up: `L4K.md`. Report: `referee_L4K.md`. Scripts: `L4K_checks/` and `referee_L4K_checks/`. The referee's code
is independent and imports nothing from the author's. All of this concerns label systems (designs) in P3-all\*.
Actual circles give the genuine system, which satisfies every statement here. **Nothing here bounds actual
circles.**

**Setting.** Labels `ell_j(q)` take values in `Z/m_q`. The parity is forced: `ell_j(q) = [(p_j/q) = -1] mod 2`.
* `n_q = m_q/4`.
* *Good up to norm `K`*: every `v` with `1 <= ||v||_1 <= K` and every unit `s` satisfy the per-unit budget
  `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`.
* *System-level coherent* at `q` with `zeta`: `ell(q) = u_q L^zeta(q) + (m_q/4) w`, with one unit `u_q` shared by
  all columns.
* `S_h = {v : ||v||_1 <= h}`.

1. **Theorem 1, the forced/free boundary (proved; independently checked: ALL OK).** Let `h = floor(K/2)`.
   * If `n_q < |S_h|` (about `(2r)^h/h!`), every label vector at `q` has a level-1 collision of norm `<= 2h`.
     For odd `K` the bound improves by a factor of about 2.
   * If `n_q > A_K = K(2(K+1)^2 B + 2K + 1)`, the explicit parity-respecting vectors (2.1),
     `a_j = 2(Lambda + (K+1) beta_(sigma(j))) + nu_j` with `beta` a `B_h` set, have no collision of norm `<= K` at
     any level, for every permutation `sigma`.
     * `B = r` for `K <= 3`.
     * For larger `K`, `B < (2r)^h` by Bose-Chowla (cited from memory), or `B <= h(r+h)^(2h-1) + 1` greedily.

   **Bounded norms `<= K` force nothing above `M^(floor(K/2) + o(1))`.** This sharpens round-9 Prop. H, which had
   `M^(K + o(1))`.
   * *Exact:* the forced/free thresholds are `q* = 31-89` for `K = 2` and `83-199` for `K = 3` (`r = 4..8`, primes
     `<= 200`), against `Q_K = 1297-117913`. Between 97% and 99.9% of the `Q_7` log-mass below `Q_K` is free.

2. **Theorem 2, refutation of L4(K) (proved, with the referee's repairs m1, m3, m5, m9).**
   * *Hypotheses.* Fix `K >= 2`, `h >= 0`, `delta in (0,1)`. Take `r <= rho M` with `r >= C(K,h,delta) log M`, and
     `Y >= Q_K` with `log Y = O(log M)`. Let `p_j = 1 mod 4` be any primes.
   * *Construction.* Genuine data at `q <= q_T = 4A_K + 5` and at `q | 2N`. At every other `q`, the vector (2.1)
     with independent uniform `sigma_q`.
   * *Conclusion.* For `M >= M_0` some realisation is parity-respecting and good for every norm `<= K` at every
     level and every `X`. It is not system-level coherent with any `zeta` of norm `<= M^h` on any set of `Q_7`
     moduli carrying more than `2 delta + o(1)` of `theta_7([1, Y])`. Moduli dividing `Norm zeta` may be counted
     either way.
   * **Hence L4(K) as posed in `round9/README.md` §6 is false for every `K >= 2`.**
   * *Not fakes.* These systems violate the norm-`(2h+2)` budgets at **every** level on a positive proportion of
     the free moduli: 63% (`r = 4`) and 100% (`r = 8`) of `q <= 20000` at `K = 2`. So they are not in P3-all\*.
     The input they falsify is goodness of the monomials of norm `> K`, which actual circles satisfy.
   * *Correction m5.* "Near-extremal clusters have `r = Theta(M)`" is unsupported (only `r = O(M)` is proved), and
     it is not needed.

3. **Lemma P (proved; level 1).** Equality of the full collision lattices gives `c' = u c`. Short generation is
   needed only to pass from the `V_K`-collision sets to the lattices.

4. **The threshold (major M1; corrected).** "`K >= 2L(h_0 + h)`, twice the README estimate", labelled proved, is
   withdrawn.
   * Linnik's `C_L k^L` is an upper bound for the least prime, so it cannot show necessity.
   * The Hamming conclusion of L4 feeds round-9 **Theorem B** (Heath-Brown), not Theorem A. Theorem B needs only
     `Y >= p_j Norm zeta`. So `floor(K/2) >= (h_0 + h)(1 + o(1))`, i.e. `K ~ 2(h_0 + h)`, suffices for genuine height
     (`v_(p_j)(Norm zeta_j)` odd) on all but `o(M)` columns.
   * The genuine-class (rigid-core) form needs neither Theorem A nor Theorem B. Its required length is set by L7,
     which nobody has quantified.

5. **The corrected question (major M2).** The write-up's displayed L4^forced(K) is withdrawn.
   * Its set `F_K = {n_q <= A_K}` consists mostly of free moduli, and there the coherent share is only `0.13-0.25`.
   * It also drops the budgets of norm `> K`, which actual circles and P3-all\* satisfy.

   The referee's replacement:

   > **L4\*(K) (open).** Regime (1.4), `r <= rho M`, `h >= h_0`. Suppose a parity-respecting system is good for
   > every monomial of norm `<= K'`, at all moduli of `F_(K')^ex`, for some `K' >= K + 1` (or all norms, as in
   > P3-all\*). Here `F_(K')^ex` is the exact forced set of the largest norm used. Must the system be in the genuine
   > class (one shared unit, conjugation per column) on an initial segment `q <= q_core` of `F_K^ex`? In Hamming
   > form: must it be coherent with some `zeta` of norm `<= M^h` on `(1 + delta)/2` of `theta_7(F_K^ex)`?

   *Exact data and evidence* (`r = 4`, primes `{5, 13, 17, 29}`):
   * The good level-1 systems on the `K = 2` forced primes number `5,066,549,580,791,808`, of which `2^37` carry the
     genuine signature everywhere (exact, reproduced).
   * *Rigid core.* At `K = 4`, every good system on `{3..47}` is in the genuine class, whether the rest of the
     forced range is genuine or carries a non-genuine top (exact: 1,966,080 systems, reproduced). Rigidity at
     `q = 3, 5, 7` is vacuous (m4).
   * *Norm-4 flexible top.* At norm `<= 4`, good systems exist that are outside the genuine class at all 18 moduli
     in `[107, 199]` (exact, reproduced). Their system-level coherent share of the forced `Q_7` mass lies between
     0.5 and about 0.71 (m10). So at Theorem B's threshold this is not a counterexample.
   * *With norm 5 added* (evidence):
     * one modulus at a time, the moduli 107, 113, 127 and 151 become rigid;
     * 167 and 199 keep 96 and 40 non-genuine survivors, against 90,048 and 500,160 at norm 4;
     * no randomized joint search finds a non-genuine top;
     * 211, the first norm-5-forced modulus above the norm-4 range, is flexible (832 survivors).

     This is consistent with a flexible top that moves up with the norm.
   * Not known: whether the rigid segment grows polynomially in `M`, or whether it survives primes with
     `log p_j/log r` bounded away from 1.

6. **Fake test (referee §7; proved).** Every fake circle has genuine primes. Give it the genuine labels
   `log(pi_j/conj pi_j)`. The result is parity-respecting, good for every norm at every level (because
   `m_(v,s) <= |D_s(Z_v)| <= nu_s^(-1) e^(w(v)/2)`), and coherent with `zeta = pi`. So the hypotheses **and** the
   conclusions of L4(K), L4^forced and L4\* hold for every fake whose primes are in the regime.
   * **L4\* uses no fake-violated input.**
   * Route 3 = L4\* + Thm B + L7 can beat the barrier only through L7, the step that links the labels to the
     angles (the slab conditions of P3-G), with `A` chosen after `c_0`.

**Status.** Route 3 as formulated in round 9 is refuted. L4\*(K) is open. It is a structural lemma about labels,
not a barrier-beating input.

## 3. analytic: Hecke characters and Gaussian primes in sectors (refereed: barrier survives, classification incomplete)

Write-up: `analytic.md`. Report: `referee_analytic.md`. Scripts: `analytic_checks/` and `referee_analytic_checks/`.
**Status: barrier. No growth improvement.**

1. **Proposition A (proved; scope corrected).** Some fibre statements are derived from identities and inequalities
   valid for every real angle vector: sector counts in the norm fibre, sizes of `S_k(N)`, Riesz products. Every angle
   model satisfies these, so they are capped by `outside.md` A.3 and `lw.md` 6.1.
   * **Not covered (referee, major):** *arithmetic* fibre inputs.
   * Example: `N^(2k) S_k(N) = sum_(|z|^2 = N) z^(4k)` is a rational integer. For `e = 1` it equals
     `4 prod_j 2 Re(pi_j^(4k))`. Its residues mod `q` are determined by the `pi_j`.
   * These use only the primes of `N`, hold for actual circles, and fail for the fakes. On four explicit fakes the
     distance of `N^2 S_1(N)/4` to `Z` is 0.095, 0.236, 0.196 and 0.417 (`construct.md` with `q = 11`; skew double
     core with `q = 3, 7, 11`).
   * Their per-prime form `p_j^(2k) e^(4ik phi_j) in Z[i]` pins `phi_j` to a finite set, so it is the original
     problem again. Their residue form is the open class P3-all.
   * So "Hecke characters on the fibre add nothing" and "separation needs external primes" are false as stated.
2. **Theorem B, the graft barrier (proved; ranges corrected).** Fix `c > 0` and `gamma in (0, 1]`. Every
   `(c, gamma)`-robust, *prime-supported* ensemble input is satisfied up to a factor `1 + o(1)` by the genuine
   Gaussian-prime angle function with a fake's `O(M)` angles grafted in. The inputs covered:
   * the Hecke large sieve restricted to primes, with constant `>= c(X + K X^(1/2))`;
   * Hecke prime sums with error `>= c X^gamma` (zero-free region, zero density, GRH);
   * Kubilius/Ricci/Harman-Lewis polar-box counts, and sieve bounds in sectors and short intervals;
   * the congruence-class versions of all of these.

   The grafted fakes come from two families:
   * the `lw.md` 6.1 fakes, with `A > max(4, 1/gamma)` (not `A >= 4`) and `W <= (2A + o(1)) M log M`;
   * the `sharp.md` P3-char fakes, with `A > max(2, 1/gamma)`, `W <= (bA + o(1)) M log M` and `b >= 33`.

   Limits:
   * The cap bars an argument only if `gamma`, `c` and the residue exponent `A` are fixed before the target
     constant. This is the round-9 R1 caveat. Published inputs have absolute `gamma`.
   * Composite-supported inputs are not covered by the proof: the all-ideal Hecke large sieve and Type I/II sums.
     They need a multiplicative graft, which is plausible but unproved.
3. **Proposition C (proved).** Counting prime tuples with a Fourier majorant and absolute-value bounds on Hecke sums
   gives a bound that is invariant under column rotations. So it is `>= 1` and can never certify that no
   configuration exists, also with distinct-tuple inclusion-exclusion.
   * The amplification statements and "Kubilius-type averaging gives exceptional sets, never empty ones" are
     heuristic.
4. **Theorem D (proved; sharpened by the referee).** `LCG_n(phi; n)` says
   `|sin(phi - arg alpha)| >= (n Norm alpha)^(-1/2)` for every `alpha` with `Norm alpha <= n` that is not parallel to
   `e^(i phi)`. It holds **iff** `phi = arg beta` with `beta` primitive and `Norm beta <= n`.
   * So angle models satisfying LCG at every prime are exactly the shadows of Gaussian-integer tuples. That class
     has exactly the actual `W_min`.
   * *Exact:* LCG holds at every actual prime, and fails at 142/142 primes of explicit realisations of the
     `construct.md` fake and of the skew double core (`q = 3, 7, 11`), with witnesses of norm `<= p_j`.
   * LCG is fake-killing only because it *is* the original problem.
5. **Proposition D.1 and Lemma D.2 (proved, conditional on Conjecture GPA).**
   * GPA: every line through 0 that is not a lattice direction of norm `<= h` passes within `h^(-1/2)` of a Gaussian
     prime of norm `!= h`. The region is a strip at all scales, not a sector at one scale.
   * GPA is *not known* to follow from GRH; the explicit formula falls short by `h^(1/2)(log X)^3`. At lattice
     directions it is Bunyakovsky-type.
   * Under GPA, the prime cross-gap PCG collapses to Theorem D.
   * D.2 (the genuine angle is isolated under PCG): its hypothesis holds for all 8977 split `p <= 2*10^5` (exact).
6. **The sharp large sieve (referee, major: misclassified).** Up to constants, the sharp `(X + K)` large sieve is
   equivalent to the max-scale cross-gap `|sin(phi_p - arg varpi)| >= c/max(p, Norm varpi)`, not to PCG.
   * *Exact:* 12,546 composite lattice directions `beta` with `p < Norm beta <= 4p` (`p = 101, 1009, 10009`) satisfy
     it against every Gaussian prime of norm `<= 10^5`. All of them violate LCG, and 9,360 violate PCG.
   * So Proposition D.1 does not apply, and its status is open.
   * At best it reconstitutes the actual problem up to `O(M)` in `W`.

**Dichotomy, as corrected.** Narrow-sector distribution of Gaussian primes, used through robust prime-supported
statements, cannot see the `O(M)` selected primes (Thm B, Prop. C). Used pointwise at LCG or PCG strength, it is the
lattice (Thm D, D.1). The dichotomy is **not exhaustive**: arithmetic fibre inputs (P3-all) and the sharp large sieve
remain unclassified. Neither is a new growth mechanism.

**Open** (the write-up's §6): the regime `P_+(N) < M^(2 - delta)`. There Theorem B does not apply, and no fakes are
known; every known fake needs `P >~ r^2 log r`.

## 4. additive: relation lattice and additive gap (refereed: theorems correct, headline does not survive)

Write-up: `additive.md`. Report: `referee_additive.md`. Scripts: `additive_checks/` and `referee_additive_checks/`.
**Status: partial. No growth improvement.** The referee re-derived every proof and ran 4,602 exact runs (1,790 with
`M >= 4`) without a counterexample. The headline "fake-violated input" does not survive.

1. **Collapse theorems (proved).**
   * *Thm 2.1.* The relation lattice `K = {c : sum c = 0, sum c_x z_x = 0}` has rank `M - 3` iff the points lie in a
     lattice coset.
   * *Thm 2.2 (rigidity).* If `M >= 5` points of a circle satisfy all relations of an actual cluster, they are that
     cluster up to similarity.
   * *Thm 2.3.* A lattice coset of `>= 5` points on a circle is a norm-form problem for some imaginary quadratic
     order, and nothing bounds the discriminant.
   * *Thm 2.6.* Every `c` with `sum c z = 0 mod N` and `||c||_1 < N/L` is an exact relation. So residues mod `N` plus
     small lifts are the actual problem.
   * *Prop 2.7.* The invariants of `K` are local data: the index is even, its odd primes are those where the cluster
     has `<= 2` residues, and projections have index `A_123/I`. Item (2) also holds for odd `q | N`.
   * *New exact identity (referee).* `covol(K)^2 I^2 = sum` over triples of `(2 Area)^2`. It holds on 1,790/1,790
     clusters.
   * *Consequence.* The geometric-mean height of the `M - 3` relations is `<= M^(c_0/4 + 1/2 + o(1))` in `l_2`. The
     write-up's `c_0/2 + 1/2` loses a curvature factor.
2. **The additive gap `AG(H)`, with corrected heights.** `AG(H)` says that every zero-sum `c` with `||c||_1 <= H` has
   `sum c z = 0` or `|sum c z| >= 1`. Actual circles satisfy it at every height. Take `M` points on an arc
   `C N^(1/4)` with `W <= c_0 M log M`.
   * **Thm R2 (referee; proved).** Every relation-free configuration violates `AG` at `l_1` height
     `M^(1 + c_0/4 + o(1))`.
     * The write-up's Thm 2.4(1), as written, gives `M^(1 + c_0/2 + o(1))`. Its `M^(c_0/2)` is the box height, not
       the `l_1` height.
     * The improvement from `c_0/2` to `c_0/4` uses curvature: the values lie in a rectangle, not a disc.
     * LLL on random arcs confirms the exponent: slopes 0.332 (`M = 6`) and 0.198 (`M = 8`), against `1/(M-3)`.
   * `AG(M^A)` forces difference-module rank `rho <= 2 + c_0 M/(4(A-1)) + O(M/log M)`. AG at height `e^(W/4+O(M))`
     forces a lattice coset.
   * The claim "violated at `M^(c_0/2+o(1))`" is proved only for `c_0 >= 4`. Heuristically it is false for
     `3 <= c_0 < 4`.
   * The ladder ceilings should read `1 + c_0/4` (relation-free) and `1 + c_0/2` (imposed relations, heuristic).
3. **Where AG's content lies.**
   * **Prop R1 (referee; proved).** For every angle model and every `c` with `v(c) != 0`,
     `|sum c w| >= sqrt2 e^((W - w(v(c)))/2) - o(1)`.
     * So the monomial gaps (class `A_0`, which the fakes satisfy) already give AG on every `c` with
       `w(v(c)) <= W + log 2 - o(1)`.
     * Every exact relation has `w(v(c)) >= 3W/2 - O(log(nC^3))`; on actual data the minimum is `1.414 W`.
   * **Prop R3 (referee; proved).** For any rank-`M` profile, the gap-admissible positions have empty interior and
     positive measure: a fat Cantor set.
4. **Fake test (referee, major: the headline does not survive).**
   * *Generic only.* `AG(M^A)` with `A > 1 + c_0/4` is violated by **almost every** member of each barrier family:
     * Cor 2.5 holds in `phi`-space. All barrier families have rank `M` and uniform `delta`.
     * Exact LLL on the skew double core finds violations at `l_1` height about `e^14` (`M = 8`, `A_0`) and `e^20`
       (`M = 8`, `A_1(2M)`).
     * 53-bit positions create spurious violations, a precision trap.
   * *Universal form.* "Every fake of size `O(M log M)` violates AG" is **Lemma (AG)**, which is open.
     * It is **strictly stronger than growth**: actual clusters are angle models satisfying `AG(infinity)`, so
       restricting Lemma (AG) to them *is* the growth statement.
     * Its `rho = 2` case contains the Gaussian-profile circle problem for every imaginary quadratic order.
     * Hybrids (an `O`-lattice cluster carrying a gap-admissible Gaussian profile) are not excluded.
   * So AG is generic-fake-violated, like any exact identity of actual circles.
5. **Pairs and lattice lines.**
   * *Prop 3.1 (proved; class `A_0`, capped).* Every pair outside one matching has `Norm g_xy >= N^(1/4)/(sqrt2 C)`.
     This answers round-9 L1/L2 open point 1 with `delta = 1/4`.
   * *Thm 3.2 and Cor 3.3 (proved; actual circles only; 3.2(c) needs `N^(1/4) >= 2C/pi`).* In low position
     (`h* <= M^b`) they give a restricted (R1-arch). They are nearly sharp on data (slacks `-0.010` and `-0.025`).
     They are not fake-violated in their conclusions.
   * *Lemmas 4.1-4.2, Prop 4.3, Remark 4.4 (proved).* These are Newton-polygon structure; on four-point relations it
     is the triangle identity.
   * *§5 (proved).* Sums, Sidon sets, the `sigma`-cocycle, the `F_q` polynomial method, determinants and Bombieri-Pila
     all reduce to capped classes.

## 5. determinant: forced divisibility in the growth regime (refereed: survives as stated)

Write-up: `determinant.md`. Report: `referee_determinant.md`. Scripts: `determinant_checks/` and
`referee_determinant_checks/`. **Status: partial. No growth improvement.** The referee found no error affecting a
theorem. The corrections are interpretive, plus wording and constant slips.

Notation:
* `n` points; `r` split primes; profile columns `beta_p in [-1, 1]^n`, weights `omega_p` with `sum omega_p <= 1`.
* `delta`: relative angles. `f_delta(k) = e^(i<k,delta>)` for `k` in the zero-sum lattice `Z^n_0`.
* `F_c(z) = N^(D/2) sum_k c_k e^(i<k,delta>)`, with `D = 2n(n-1)`. It is a Gaussian integer for actual clusters.
* `V(A, m)`: integer vectors `c` supported on `A` that annihilate all polynomials of degree `< m` on `A`.
* `V_H`: the span of its integer members with `|c|_inf <= H`.

1. **Theorem 2.1, LP supercriticality (proved; thresholds reproduced with exact big integers).** Let `n >= 64`, take
   any `r` columns and weights, and set `T = n` and `w = n sqrt(8(n-1) ln 4r)`. The slab-box is
   `A = {k in Z^n_0 : |k_x| <= T for x < n, |<k, beta_p>| <= w for all p}`. Let `m` be the largest integer with
   `binom(m+n-2, n-1) <= (2T+1)^(n-1)/8`, and let `h_A` be the width. Then:
   * `|A| >= (2n+1)^(n-1)/2`;
   * `m >= 0.7 n(n-1)` and `h_A(beta_p) <= 2w`;
   * `dim V_H >= |A|/2`;
   * `m/(2w) >= 0.12 sqrt((n-1)/ln 4r)`.

   `m >= 4w` holds exactly from `n = 2137` (`r = n`), `2679` (`r = n ln n`) and `4278` (`r = n^2`), and always when
   `n - 1 >= 300 ln 4r`. Further:
   * *Thm 2.3.* The Hadamard profile has ratio `<= sqrt n`.
   * *Cor 2.4.* An LP-critical profile needs `r >= e^(0.0144(n-1))/4`.
   * **Consequence (correction to the growth map).** "Forced divisibility is exactly critical"
     (`threshold/README.md`, `round8/README.md` §2.4) holds only for a form fixed before the cluster (cube measure,
     `2^n` columns). With the support chosen from the cluster's own profile, the LP is supercritical in the whole
     growth regime. The obstruction is nonvanishing, not criticality.
2. **Theorem 4.1, rigidity (proved; checked end to end on 240 kernel vectors over actual clusters).** Suppose an
   actual cluster has `n >= 64`, `n - 1 >= 300 ln 4r` and `W >= 16 ln n + 8 log^+ C + 14`. Then `F_c(z) = 0` for
   every integer `c in V(A, m)` with `|c|_inf <= H`, hence on `V_H`. Equivalently, `f_delta|_A in V_H^perp`, a space
   of dimension `<= |A|/2`.
   * *Proof idea.* A Taylor/Liouville upper bound meets a lower bound from the integrality and forced
     `pi`-divisibility of `F_c(z)`.
   * The theorem needs `>= 2137` points. No such cluster is known, so it may be vacuous, and it gives no bound by
     itself.
3. **Proposition 4.2 and Corollary 4.3 (proved; "almost every").** `Z_P = {delta : (4.2)}` is a closed real-analytic
   Lebesgue-null set. So relation (4.2) fails for almost every `delta` of every barrier family at the same
   `(n, r, W)`.
   * The barrier families are the round-5 fakes (`M >= 3527`) and the skew double core (`M >= 3030`), at Paley
     orders.
   * *Small-scale illustration:* 400/400 random fake positions are certified off `Z_P`.
   * **Interpretive caveat (referee, major).** The separation is generic: the property "points are Gaussian
     integers" separates in exactly the same way.
     * "*Every* fake violates (4.2)" is a special case of ZE. A fake steered onto a nondegenerate point of `Z_P`
       would refute ZE and carry the barrier over to the nonvanishing step.
     * The separating input is integrality plus forced divisibility of multi-term forms `F_c(z)`, which
       `outside.md` §6 item 1 already places outside the barrier.
     * **What is new is that this input is LP-supercritical in the growth regime**, so it yields an exact
       relation.
4. **Rank (proved).**
   * *Lemma 3.1.* Squarefree clusters have exponent rank `>= n/4 - 2`.
   * *Lemma 3.2.* Every normalised cluster has a nondegenerate full-rank subcluster. The proof uses the `(pi/4)Z` gap
     of `outside.md` Lemma 0.1.
   * *Prop 5.2.* Low-rank profiles (rank `<~ sqrt(n/log r)`) lie in `Z_P` automatically. With `kappa = 0`, as the
     referee justifies by Cramer's rule, they are invisible to the method.
5. **Theorem 5.7 (proved, conditional on ZE).**
   * **Hypothesis ZE(n).** For every profile with `n - 1 >= 300 ln 4r`, no nondegenerate `delta`
     (`sum v_x delta_x != 0` for all `0 != v in Z^n_0`) with `max |delta_x| <= 1/(4D)` lies in `Z_P`.
   * **Theorem.** If ZE(`n'`) holds for all `n' >= n_1`, normalised clusters on squarefree norms satisfy
     `n < 4 max(n_1, 300 ln 4r + 1) + 4`, i.e. `n = O(log W) = O(loglog R)`.
   * *Status of ZE.* Non-standard and unproved. It is an unlikely-intersection statement: a toric variety against a
     linear space of codimension `>= |A|/2`.
   * *The weakest form that suffices* (write-up open point 1). Restrict ZE to algebraic points
     `w = e^(i delta) in (Q(i) cap S^1)^(n-1)` that are multiplicatively independent, of height `<= W`, and within
     `e^(-W/4)` of `1`.
   * *ZE must also exclude translates `gH` of subtori near `1`.* Bombieri-Masser-Zannier does not do this
     effectively.
   * **Prop R1 (referee; proved).** If the slabs do not cut (`A` = full box, `2^(2T(n-1)) <= H`), then
     `V_H = V(B_T, m)` and `Z_P` is exactly the sparse strata, which are degenerate. So ZE holds there, and **all of
     the difficulty is in the slab cuts**.
   * *Evidence.* At the write-up's sizes (`n = 4, 4, 5`), `V_H = V(A, m)` and the base locus is only the identity and
     the sparse strata. These sizes are far below the regime, and the slabs do not cut there.
6. **Negative results of §5 (proved, with corrected constants).**
   * Parallelepiped certificates reduce to single characters (Prop 5.3; (iii) should read
     `X_v >= W/2 - 2 log(C|v|_1/(2 sqrt 2))`).
   * Valuation certificates are impossible at `e_p <= 3`, and residue collisions favour vanishing (Lemma 5.4,
     Prop 5.5).
   * Grid-lemma arguments are critical on balanced profiles (Prop 5.6). The expectation identity is proved; the
     "at best critical" reading is heuristic.
   * Prop 5.1 holds for points coinciding *with the base point `n`*.

## 6. Scores

Sub-scores follow the rubric of §1: (a) correctness, (b) fake-violated input, (c) movement of growth or uniformity.

| mechanism | (a) | (b) | (c) | **overall** | reason |
|---|---|---|---|---|---|
| **determinant** | 8 | 5 | 3 | **5** | All theorems survive. Thm 2.1 is a proved correction to the growth map: the form chosen after the profile makes forced divisibility supercritical in the growth regime. This is the only mechanism of the round whose output is an exact constraint that fakes do not satisfy, with a residual hypothesis (ZE) that no longer involves arithmetic. The conditional conclusion (`M = O(loglog R)` on squarefree norms, where every known fake lives) is far beyond the target. Deductions: (b) the separation is only "almost every", and "every fake" is part of ZE; (c) ZE is non-standard and unproved; Thm 4.1 needs `>= 2137` points; non-squarefree norms need a separate rank bound. |
| **L4K** | 7 | 0 | 2 | **3** | Two clean proved theorems: Thm 1 sharpens Prop. H to `M^(floor(K/2))`, and Thm 2 refutes L4(K) for every `K`. The exact data reproduce to the last digit. This decisively redirects route 3. Deductions: two major corrections (the threshold claim labelled proved; a corrected question that fails its first test). (b) is 0 because the referee proved that every fake with genuine labels satisfies L4\*, so no input of the surviving lemma is fake-violated. |
| **analytic** | 7 | 3 | 1 | **3** | Thm B is a correct and useful barrier: robust prime-supported analytic inputs cannot see `O(M)` grafted angles. Thm D characterises LCG exactly. Deductions: two major misclassifications (arithmetic fibre inputs; the sharp large sieve), so the advertised dichotomy is not exhaustive. Its fake-violated input (LCG) is the original problem. |
| **additive** | 5 | 3 | 1 | **2** | The collapse theorems, Prop 3.1 and Thm 3.2 are correct and nearly sharp on data. The headline does not survive: the AG heights are wrong (`M^(c_0/2)` should be `M^(1 + c_0/4)`), the violation is generic only, and Lemma (AG) is strictly stronger than growth. Prop 3.1 settles a round-9 open point but is capped. |

**Reading.** No mechanism reaches 6. That would need a growth theorem from a hypothesis that is standard or at least
plausibly provable, or a fake-violated input weaker than the problem with a credible route to growth. The
determinant mechanism scores highest. It is not close to a proof. It is the only route of this round, and the first
since round 8, where the arithmetic input has been *spent* to produce a geometric statement (ZE) that can be tested
and attacked without reference to fakes.

## 7. Updated status

**Unchanged.**
* No growth improvement and no uniform bound: `M <= (2/3 + eps) log R/loglog R`.
* The only conditional route known to give the uniform bound is still Vojta's conjecture for one fixed variety.
* P3-all\* is open, with the bracket of §0.
* The round-9 R1 caveat stands: the residue caps are proved only at fixed `A`. The caps of Thm B (§3) carry the same
  caveat in `gamma` and `A`.

### 7.1 Additions to the barrier map (refereed this round)

| input or method | status | source |
|---|---|---|
| robust prime-supported analytic ensemble inputs (Hecke large sieve over primes with constant `>= c(X + K X^(1/2))`; Hecke prime sums with error `>= c X^gamma`; Kubilius/Ricci/Harman-Lewis sector counts; sector and short-interval sieve bounds; congruence versions), with `c`, `gamma`, `A` fixed before `c_0` | **capped** by grafts of `lw.md` 6.1 (`A > max(4, 1/gamma)`, `(2A + o(1)) M log M`) or of `sharp.md` P3-char fakes (`A > max(2, 1/gamma)`, `(bA + o(1)) M log M`) | `analytic.md` Thm B |
| absolute-value Fourier counting of prime tuples | **cannot certify emptiness** (rotation invariant, `>= 1`) | `analytic.md` Prop C |
| real-analytic fibre statements (sector counts, `S_k(N)` sizes, Riesz products) | **capped** (they are functions of the angle vector) | `analytic.md` Prop A |
| lattice cross-gap LCG at every prime | **equivalent to the problem** (LCG iff `phi_j = arg beta_j`, `Norm beta_j <= p_j`) | `analytic.md` Thm D (sharpened) |
| prime cross-gap PCG | equivalent to the problem **under GPA** (not known from GRH) | `analytic.md` Prop D.1 |
| full relation lattice; residues mod `N` plus small lifts | **equivalent to the problem** (norm-form problem for an arbitrary imaginary quadratic order) | `additive.md` Thms 2.1-2.3, 2.6 |
| additive gap `AG` on combinations of weight `w(v(c)) <= W` | **implied by the monomial gaps**, hence capped | referee Prop R1 (`referee_additive.md`) |
| `AG(M^A)`, `A > 1 + c_0/4` | violated by almost every fake; the universal form is Lemma (AG), **strictly stronger than growth** | `additive.md`, referee Thm R2 |
| pair gcd bound `Norm g_xy >= N^(1/4)/(sqrt2 C)` off one matching | **capped** (class `A_0`) | `additive.md` Prop 3.1 |
| norm-`<= K` label budgets | force collisions only for `n_q <~ (2r)^(floor(K/2))`; nothing above `M^(floor(K/2)+o(1))` | `L4K.md` Thm 1 |
| label lemmas of L4 type (parity plus budgets, any norms) | **no fake-violated input**: every fake with genuine labels satisfies hypotheses and conclusion | `referee_L4K.md` §7 |
| auxiliary forms with forced divisibility, support chosen from the profile | **not capped; LP-supercritical** in the growth regime; outputs the exact relation (4.2); growth reduces to the nonvanishing hypothesis ZE (squarefree norms) | `determinant.md` Thms 2.1, 4.1, 5.7 |

### 7.2 Unclassified after this round

These are true for actual circles, and are neither capped nor shown to give growth:
* arithmetic fibre inputs (integrality and residues of `N^(2k) S_k(N)`): the P3-all class;
* the sharp `(X + K)` large sieve, which is equivalent to the max-scale cross-gap;
* composite-supported analytic inputs (all-ideal large sieve, Type I/II sums);
* the regime `P_+(N) < M^(2 - delta)`, where no fakes are known;
* ZE and its algebraic-point form.

### 7.3 Routes

* **Route 3** (labels, i.e. L4\* + Thm B + L7) is reformulated. L4(K) is refuted, and L4\*(K) is open. Its analytic
  threshold is `K ~ 2(h_0 + h)` via Theorem B. Since L4\* uses no fake-violated input, the step that must beat the
  barrier is L7 (label-angle). No handle on L7 has appeared in rounds 7-10.
* **Routes 1, 2 and L5:** unchanged from `round9/README.md`. Prop 3.1 of `additive.md` gives route 1 the
  `delta = 1/4` pair bound, which is capped.
* **Route 5 (new; determinant).** Thm 4.1 (proved) together with ZE (open) gives `M = O(loglog R)` on squarefree
  norms (Thm 5.7). Non-squarefree norms additionally need a linear lower bound on the rank (`determinant.md` open
  point 3).

### 7.4 Corrections to earlier READMEs

The repository is read-only. These are to be applied when round 10 is merged.
1. **`round9/README.md` §6.** L4(K) as displayed is false for every `K >= 2` (`L4K.md` Thm 2). The scale estimate
   "`K >~ L(h_0 + h)`" should read `K ~ 2(h_0 + h)`, because the Hamming conclusion feeds Theorem B.
2. **`round9/README.md` §5 Prop. H.** Sharpen from `M^(K + o(1))` to `M^(floor(K/2) + o(1))` (`L4K.md` Thm 1).
3. **`round8/README.md` §2.4 and the `threshold/` summary in `docs/claude_review/README.md` §3.**
   * "Forced divisibility is exactly critical" holds only for an auxiliary form fixed before the cluster.
   * With the support chosen from the profile, the LP is supercritical by a factor `>= 0.12 sqrt((n-1)/ln 4r)`
     (`determinant.md` Thm 2.1), and the obstruction is nonvanishing.
   * The single-place `tau* < 2` statement is unaffected.
4. **`round8/README.md` §2.4,** row "counting under GRH for Hecke L-functions (heuristic)". There is now a proved
   cap for robust prime-supported inputs, with fixed `gamma` and `A` (`analytic.md` Thm B).
5. **`round7/README.md` §2.1 row 9** (additive data, open). In full, additive data are the problem
   (`additive.md` §2). On weights `<= W` they are implied by the monomial gaps (referee Prop R1). In between, the
   class is open.

## 8. The single most promising next step

**Recommended: decide Hypothesis ZE for the slab-cut forced-divisibility system (route 5).**

### 8.1 Statement to decide

Fix a profile on `n` points with `r` columns, and let `(A, m)` be its slab-box and order from Theorem 2.1. Let
`L(A, m)` be the linear system of Laurent polynomials `f = sum_(k in A) c_k w^k` on the torus `T = (C^*)^(n-1)` with
`ord_e f >= m`, where `e` is the identity; its coefficient space is `V(A, m)`. Its base locus is
`Bs(A, m) = Z(V(A, m))`. Since `V_H ⊂ V(A, m)`, `Z_P ⊇ Bs(A, m) ∩ (real torus)`.

> **ZE-real(A, m).** Near `e`, every point of `Bs(A, m)` on the real torus `(S^1)^(n-1)` is degenerate, i.e. lies on
> a proper rational subtorus *through `e`*. So neither translates `gH` (`g notin H`) of subtori nor non-toric
> components may reach nondegenerate points near `e`. The height-restricted version uses `V_H` in place of
> `V(A, m)`.
>
> **ZE-alg** (the weakest form Theorem 5.7 needs; `determinant.md` open point 1). Consider only points
> `w in (Q(i) cap S^1)^(n-1)` with multiplicatively independent coordinates, height `<= W`, and within `e^(-W/4)`
> of `e`. Then none lies on `Z_P`.

### 8.2 Why this step

1. **The arithmetic has already been spent.** Theorem 4.1 is unconditional and refereed. It converts integrality and
   forced divisibility, which are true for actual circles and false for fakes off `Z_P`, into the exact relation
   (4.2).
   * What remains, ZE-real, mentions no fakes, labels, residues or primes. It is a statement about an explicit
     linear system on a toric variety.
   * Every other open route still has its missing lemma *inside* the barrier problem: L7 for route 3, Lemma (AG),
     and P3-all itself.
   * This is the classical shape of a transcendence proof: auxiliary polynomial, then vanishing, then zero
     estimate.
2. **The payoff is far beyond the target.** ZE-alg gives `M = O(loglog R)` on squarefree norms (Thm 5.7). Every known
   fake is squarefree (`round8/README.md` §2.4), so squarefree norms are exactly where the barrier lives.
3. **The question is already localised.**
   * Without slab cuts ZE holds (referee Prop R1).
   * Low-rank strata are in the base locus for a known reason (Prop 5.2), and they are degenerate.
   * So the open content is the effect of the cuts on components through nondegenerate points.
4. **It is finite-testable now**, at model sizes (§8.4).

### 8.3 What can make ZE-real fail (*README remark*, unrefereed)

This is standard intersection theory on the projective toric variety `X_A` of `conv(A)` with its line bundle `L_A`.
Note that `A` is the full lattice-point set of the polytope `B_T ∩ slabs`.
* **Criterion.** If an irreducible curve `Y ∋ e` satisfies `L_A · Y < m mult_e(Y)`, then `Y ⊂ Bs(A, m)`. Restricting
  a section that vanishes to order `m` at `e` to `Y` would otherwise give a degree below its number of zeros.
* For a curve in the torus, `L_A · Y = sum_i h_A^+(u_i)`. Here the `u_i` are the rays of its tropicalisation,
  oriented so that a one-parameter subgroup in direction `u` gets the width `h_A(u)`, with multiplicity, and
  `h_A^+(u) = max_(k in A) <k, u>`.
* *Toric curves.* Every rational direction `u` with `h_A(u) < m` gives a subgroup in the base locus. Examples are
  `u = beta_p` (`h_A <= 2w < m`) and the coordinate directions `e_x`, `x < n` (`h_A <= 2T`). They lie in the
  degenerate strata of Props. 5.1-5.2.
* *Elementary fact.* An irreducible curve `Y` in the torus lies in a coset of the subtorus spanned by its rays: a
  character that vanishes on all the rays has no zeros or poles on `Ybar`, so it is constant on `Y`. If `e in Y`, the
  coset is the subtorus itself. So a curve through `e` with nondegenerate points needs balanced rays spanning all of
  `R^(n-1)`, hence **at least `n` rays**.
* *Non-toric curves.* A curve through `e` with `L_A · Y < m mult_e(Y)` and a nondegenerate real branch would refute
  ZE-real.
  * On blowups of toric surfaces at the identity point, non-toric negative curves through the identity do exist.
    These are the examples of González Anaya, González and Karu, and the earlier Goto-Nishida-Watanabe and
    Castravet-Tevelev examples, all quoted from memory.
  * So this failure mode is real, not hypothetical.
* *Where the margin is (heuristic).* By the elementary fact, a nondegenerate curve with `mult_e = 1` has
  `L_A · Y >= min sum_i h_A^+(u_i)` over balanced, spanning ray multisets. This is a convex-geometry quantity of
  the slab polytope.
  * Spanning forces several rays in directions where `A` is not thin. For example, the coordinate fan
    `e_1, ..., e_n` (a line-type curve) gives `(n-1)T + max_A k_n >= (n-1)T`, which is about `n^2` at `T = n`, against
    `m ~ 0.74 n^2`.
  * *Exact* (`judge_checks/m_ratio.py`): `m/(n(n-1)) = 0.7445, 0.7366, 0.7362` and `(n-1)T/m = 1.343, 1.358, 1.358`
    at `n = 64, 2137, 4278`. Also `m/(2w) = 2.000` at `n = 2137`, `r = n`.
  * Rays in the `beta_p` directions cost about `w` each, and spanning needs `>= n` of them, so about `n w >> m`.
  * Raising `mult_e` costs about `n - 1` conditions per unit, against one parameter per extra ray.
  * So degree-forced nondegenerate curves are not expected, but the margins are constant factors, not powers of
    `n`. This ignores curves of positive genus, higher-dimensional components, translates `gH`, and base-locus
    components that are not forced by degree. It is a heuristic only.

### 8.4 First test (exact linear algebra, plus certified numerics; evidence)

1. **Model sizes.** Take `n - 1 = 2, 3, 4`, squarefree profiles (random `+-1` columns, `r = n` to `3n`), and also
   the skew-double-core profile at `q = 3`. Choose `(T, w, m)` so that:
   * (i) the slabs cut, i.e. `A` is a proper subset of `B_T`;
   * (ii) `2w < m`, so the model is supercritical in every `beta_p` direction;
   * (iii) `|A| > binom(m+n-2, n-1)`, so `V(A, m) != 0`.

   Compute an exact integer basis of `V(A, m)` and of `V_H`. (The referee's `ref_VH_small.py` and the author's
   `baselocus*.py` are the starting harness.)
2. **Tangent cone.** Find the real common zeros on the unit sphere of `Phi_c(u) = sum_k c_k <k, u>^m`,
   `c in V(A, m)`. Are they all on rational subspaces?
3. **Base locus.** Compute the components of `Bs(A, m)` through or near `e`, by Gröbner bases for `n - 1 = 2, 3` and
   numerical irreducible decomposition beyond. Classify each as a subtorus through `e` (degenerate), a translate `gH`
   with `g notin H`, or non-toric. For `n - 1 = 2`,
   check every low-degree curve through `e` against the criterion of §8.3.
4. **Freedom.** Intersect over several admissible `(T, w)` and over sub-clusters (`determinant.md` open point 2).
   Does any non-toric component survive all of them?
5. **Fakes.** Draw gap-admissible positions for the model profiles and certify them off `Z_P` (as in
   `ref_fakes.py`). Then search *adversarially* for gap-admissible nondegenerate `delta` on `Z_P`. ZE implies there
   are none.

### 8.5 Reading of the outcomes

* **Near `e`, only subtori through `e` in every model.** This is evidence for ZE. The next step is a proof of ZE-real through
  Seshadri-type lower bounds for `L_A` at `e`, of the form `L_A · Y >= m mult_e(Y)` for every non-subtorus
  `Y ∋ e`, with the analogue for higher-dimensional `Y`. Tropical intersection numbers make this a statement about
  the slab polytope and balanced fans.
* **A non-toric component, or a translate `gH`, with nondegenerate real points near `e` that survives the intersection
  of step 4.** Then ZE-real is false, and only ZE-alg (the arithmetic of `w`) is left.
  * If, in addition, a gap-admissible fake lies on it, the barrier transfers to the nonvanishing step
    (`referee_determinant.md` §4). Route 5 is then capped as formulated.
  * Either way the outcome is decisive for route 5.

**Fake test for this step.** The input is the integrality and forced divisibility of the multi-term forms `F_c(z)`,
`c in V_H`. It is true for every actual cluster (Thm 4.1). For a fake, the dichotomy "`F_c = 0` or
`Norm F_c >= N^D prod p^(-h)`" fails at every `delta` off `Z_P`, which by Prop 4.2 is almost every `delta` of every
barrier family. Whether it fails for *every* fake is exactly what test 5 probes.

### 8.6 Also to do (cheaper, or needed for the record)

1. **Non-squarefree norms for route 5.** Prove a linear lower bound on the exponent rank of near-extremal
   non-squarefree clusters, the analogue of Lemma 3.1 (`determinant.md` open point 3).
2. **Route 3.** L4\*(K) at `K = 5, 6`, `r = 4..6`, by ILP over signatures. Does the rigid core grow with `K`? This is
   worth doing only alongside some handle on L7, since L4\* alone cannot beat the barrier (§2 item 6).
3. **Analytic leftovers.** Classify the sharp large sieve (capped, or lattice-only). Decide the regime
   `P_+(N) < M^(2 - delta)`: build a fake there, or prove a restricted theorem.
4. **Round-9 record items.** Referee sketch R2: `A`-uniform caps for `A_X` and P3-char. Prop. H is now superseded by
   `L4K.md` Thm 1.

**Not recommended.**
* Lemma (AG): it is strictly stronger than growth.
* Further robust analytic ensemble inputs: they are capped by Thm B at fixed `gamma`, `A`.
* L4(K) and L4^forced in their displayed forms: refuted, or failing the first test.

## 9. Files

All paths are relative to this folder (`round10/`). Each script runs as `python3 <script>` from its own directory.
Outputs sit next to the scripts.

| write-up | report | scripts (author / referee) |
|---|---|---|
| `L4K.md` (notes: `L4K_notes.md`) | `referee_L4K.md` | `L4K_checks/` / `referee_L4K_checks/` (notes in `notes.md` there) |
| `analytic.md` (notes: `analytic_notes.md`) | `referee_analytic.md` (notes: `referee_analytic_notes.md`) | `analytic_checks/` / `referee_analytic_checks/` |
| `additive.md` (notes: `additive_notes.md`) | `referee_additive.md` (notes: `referee_additive_notes.md`) | `additive_checks/` / `referee_additive_checks/` |
| `determinant.md` (notes: `determinant_notes.md`) | `referee_determinant.md` | `determinant_checks/` / `referee_determinant_checks/` (notes in `notes.md` there) |

The judge's working notes are in `judge_notes.md`. `judge_checks/m_ratio.py` holds the exact computation of the
order `m` quoted in §8.3.

**Caveats.**
* Several analytic inputs are quoted from memory, and only the shape of their error terms is used: Linnik's
  exponent, Heath-Brown's quadratic large sieve, Harman-Lewis, Kubilius, Ricci, the Hecke large sieve, Bose-Chowla.
* The toric-blowup literature in §8.3 (González Anaya-González-Karu, Goto-Nishida-Watanabe, Castravet-Tevelev) is
  also quoted from memory. No proof here uses it.
* Two of the author's `L4K` runs were still in progress at write-up time: `deform_4_4_last.txt` (complete to
  `q_0 = 151`) and `joint_4_4_first.txt`. No refereed statement depends on them.
* The `r = 8` coherence sampling of `L4K.md` was abandoned for run time.

---------------------------------------------------------------------------------------------

**Repository note.** `L4K_checks/forced_sig_2_4_first.pkl` (2.1 MB) is not committed; `forced_sig.py`
regenerates it.
