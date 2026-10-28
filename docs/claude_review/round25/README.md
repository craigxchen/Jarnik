# Round 25 — judge report (intra-cluster multi-point Liouville; Proposition A' on realisation cosets)

| key | write-up | referee report | question |
| --- | --- | --- | --- |
| intra | `intra.md` (`intra_checks/`, `intra_notes.md`) | `referee_intra.md` (`referee_intra_checks/`, `referee_intra_notes.md`) | Does a Thue/Siegel/Dyson multi-point argument on the `binom(M,5)` gcd-reduced 5-subtuples of one cluster (points of `T_5` near `e` that **share** primes) beat the one-point caps? This is round 24's open point 1. |
| coset | `coset.md` (`coset_checks/`, `coset_notes.md`) | `referee_coset.md` (`referee_coset_checks/`) | Round 21's recommended item: extend Proposition A' to window realisation cosets with `m >= cM`, which would decide scenario P versus R for the `A_0` question (Q). |

All paths are relative to `scratchpad/round25/`, and the repository was only read. Labels follow the overview:
proved / conditional / exact / evidence / heuristic / unrefereed. §1 gives the **referee-corrected** statements.
Remarks marked *judge* are mine and are **unrefereed**. My computations are in `judge_checks/` (§5).

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the
best bound for every `R`. Neither write-up proves `M = o(log R/loglog R)`, `(L_5^circ)`, `A_5(delta)` for any
`delta < 1/30`, or a restricted growth theorem in the low-height window. Both survive refereeing in their mathematics.
Both headlines needed scope corrections, and one of them was **reversed**.

* **intra: reduces to route 5 (confirmed).** The intra-cluster method gains nothing for growth or uniformity.
  * All subtuple points are projections `pr_S(Z)` of one point `Z` of the cluster torus `T_M`. Every intra-cluster
    `n`-point comparison is one-point toric Liouville on `T_M` at `Z` (Theorem 3.1).
  * Conversely, a star partition `T_M ≅ T_5^{(M-1)/4}` makes route 5's forced-divisibility forms intra-cluster forms
    (Theorem 4.1). The two coincide as classes of auxiliary polynomials with their non-archimedean bounds.
  * Round 24's open point 1 (shared primes) is answered **negatively** for intra-cluster configurations, by the
    referee and against the author's headline. Shared primes **do** beat the one-point `15/4` cap at the data level
    (referee Prop R1, extended-BCH cut designs). They do it through route 5's LP supercriticality, so the output is an
    exact algebraic relation and not a contradiction. Turning that relation into a bound needs route 5's open zero
    estimate (ZE-col\*; MIN^pair / (R_H)).
  * The cap `15/4` survives only on the "cube" class. There it gives a proved barrier against the **uniform** bound
    for Taylor/Dyson-type archimedean bounds. It is not a growth barrier.
* **coset: partial (confirmed).** Proposition A' was extended to realisation cosets (Theorem 4.1). Then the author
  proved that certificates of that form **cannot** reach the window when `m >= cM` (Theorem 5.1 / Corollary 5.2). So
  round 21's recommended route is closed negatively, for its certificate form. Three other results hold:
  * on the full torus, `mu(G(p)) > 0` iff `5 ∉ p`, at every `Lambda` (Theorem 2.2);
  * Haar positivity on a coset is **equivalent** to the existence of a point with valid characters and uniformly slack
    non-characters (Theorem 3.1 with the referee's Prop R1);
  * arithmetic seeds are actual clusters, so seeding with them is circular (Prop 6.1).

  Scenario P is proved only for small free dimension and conditionally (Corollary 7.3). The window stays **open** for
  `A_0`, and by round-21 J2 neither answer to (Q) gives growth.

**Scores** (§6): intra 3, coset 3. **Decision** (§8): stop internal rounds, as rounds 19, 20 and 22–24 already
recommended. Fold §7 into the expert note. The intra and route-5 lines now end at one missing lemma, route 5's zero
estimate at the actual point. It has been open since round 11 and is the record's one internal decisive item. Expected
yield from another internal round on it is very low (§8).

## 1. Referee-corrected statements

### 1.1 intra (`intra.md`; `referee_intra.md`): survives; status `reduces_to_route5` confirmed; headline item 4 reversed

Setting. A normalised cluster `z_1..z_M` on an arc of length `C sqrt R` is one point `Z = [z] in T_M(Q)`, where
`T_M = (R^1_{Q(i)/Q} G_m)^M / diag`. For a 5-subset `S`, `y_S = pr_S(Z) in T_5(Q)` is the gcd-reduced subtuple. `L_intra`
is the class of comparisons "archimedean upper bound for `|F(y_{S_1},...,y_{S_q})|` `>=` product-formula lower bound",
with toric forced divisibility at all places, shared primes, integrality and conjugates/group words. `F` may be chosen
after the profile.

> **Round 25 intra (referee-corrected).** Status: reduces to route 5. No growth bound and no uniform bound.
> 1. *(proved)* **One point.** The heights, gcds, valuations and shared primes of all `y_S` are linear functions of
>    `v_pi(Z)` (Lemma 2.1). The gcd reduction changes no inequality (Lemma 2.2, a multi-point Lemma K1); the common
>    factor is tracked exactly.
> 2. *(proved)* **Reduction and converse.** Every `L_intra` comparison is implied by toric Liouville `(TL_M)` at `Z`
>    for `G = Phi^*F`, with the same value `G(Z) = F(y)`, `RHS_G <= RHS_F` and `ord_e G >= idx F` (Theorem 3.1).
>    Conversely, a star partition gives `T_M ≅ T_5^{(M-1)/4}` over `Q`, and route 5's slab-box forms and its forced
>    divisibility (I1) are intra-cluster forms (Theorem 4.1; block degree `4n`, since route 5's box has `T = n`). The two
>    coincide **as classes of auxiliary polynomials with their non-archimedean bounds**. `L_intra`'s anisotropic Dyson
>    weights are an archimedean refinement that route 5 does not use.
> 3. *(proved)* **Cube barrier (uniform bound only).** Suppose `Z` is simple and `eta`-balanced on **all** `2^{M-1}-1`
>    cuts of `[M]`, with `eta < 1/(2^M-1)`. Then no comparison with an isotropic Taylor upper bound excludes `Z`, for
>    any `G`, uniformly in `G` once `w >= 4(log(kappa/C) + 1.39 + o(1))/(1-(2^M-1)eta)`; the uniformity comes from
>    referee Lemma R2, `C_T^-(G) <= 1.39 ord_e G`. The Dyson-form analogue holds only for families with
>    `C^-(F) = o(idx F · w)`. Zero estimates do not change this. The subtuples sit at `mu >= 4 - 2^{3-M} - O(eta)`. So
>    within Taylor/Dyson-type archimedean bounds the class proves neither the uniform bound nor anything below
>    `M ≈ log_2 log R`. Direction-sensitive bounds (route 5's leading form `Phi_c(u)`) are **not** covered.
> 4. *(proved from the refereed route-5 Theorems 2.1(d) and 4.1, plus exact design facts; referee Prop R1)* **Round 24
>    open point 1, intra-cluster case: there is NO data-level cap in general.** Take a 5-wise uniform linear cut design.
>    For example, use extended double-error-correcting BCH with `t >= 13`, `M ≡ 1 (4)`, `M >= 8189` and
>    `r = 2^{2t+1} - 2` primes. Then every 5-subtuple is exactly balanced and lies strictly inside `mu < 4(1-eta)`.
>    Yet the intra-cluster `n`-point comparison on the star subtuples fails unless `F(y) = 0`, by route 5's LP
>    supercriticality. So shared primes beat the one-point cap `15/4`, but only by forcing **exact relations**. Using
>    those relations needs route 5's zero estimate, which is open. The `15/4` cap holds on the cube class (item 3).
>    The disjoint-primes hypothesis of round 24 Prop 4.3 and of the expert note's anisotropy paragraph is
>    **essential**. "Maximal prime sharing" has no meaning here: what separates the two cases is the number of
>    distinct cuts (`log r` against `M`).
> 5. *(proved for lacunary degrees, assuming `R'_1 >= rho d_1`; heuristic otherwise)* With lacunary block degrees, the
>    `n`-point comparison is the leading-block one-point comparison up to `O(kappa omega/rho)`, where
>    `kappa <= 2 + O(log C/W)`. The claim that comparable degrees leave only route 5's zero estimate as a zero step is
>    heuristic. Esnault–Viehweg/Nakamaye-type Dyson lemmas are weak there, not inapplicable.
> 6. *(proved)* **Growth regime.** In route 5's regime the method is LP-supercritical, and its output is exact
>    relations `F(y_{S_1},...,y_{S_q}) = 0` among star subtuples (Cor 7.1). Growth needs route 5's zero estimate at
>    the actual point (ZE-col\*; MIN^pair or (R_H) would suffice), which is open. `T_5` counting cannot supply it,
>    because the base locus depends on the cluster's own profile (template-dependent).
> 7. *(evidence)* At `M = 10`, within the scanned slab family, nothing beyond the monomial and pair relations is
>    forced, and the monomial relations are fake-blind.

*Corrections adopted.*
* **M1 (major; headline reversed).** The author's §0 item 4 and summary said Theorem 5.1 "settles round 24's open
  point 1 for intra-cluster configurations, with cap `15/4` even with maximal prime sharing". Theorem 5.1 needs `Z`
  balanced on every cut, hence `r >= 2^{M-1}-1` primes. Subtuple balance needs only a 5-wise uniform cut family, which
  exists with `r ≈ 2M^2` (exact: 32/63, 64/255, 64/127 and 256/32767 cuts carried). In route 5's regime such data
  violate the intra comparison unless `F(y) = 0` (item 4). Judge check J1 extends the exact design facts to `t = 5`.
* **M2 (major, repaired for Taylor).** The threshold in Cor 5.2(2)/5.3 contained `|C(F)|/idx F`, which grows with
  `log H(F)`. Cor 5.3 then asserted one family of data defeats every inference: "for each `F` there is `w`" was
  presented as "there is `w` for every `F`", which is the owner's forbidden step. Only `C^-` matters. Lemma R2 bounds
  it uniformly for the isotropic Taylor bound; the Dyson-form version is unproved.
* **M3 (major, scope).** Def 3.0 admits **any** archimedean bound. Theorem 5.1 and Cor 5.3 cover only Taylor/Dyson
  bounds. Cor 4.2 ("the classes coincide") holds only for the polynomial class and the non-archimedean side.
* **Minor.** Prop 6.1's zero-lemma dichotomy is heuristic. `kappa <= 2` and `h_B >= W/2` need `C <= 1`. "Strict in
  47/47" is built into the construction (overlapping blocks). The sharing gain is sample-dependent (author 3.5–64,
  referee 0–134), and comparing it with `W` means nothing. "Nothing supercritical at `M = 10`" holds within the scanned
  family only. `Phi` is a closed immersion after padding, not a surjection. Cor 5.3 should read `(2^{M-1}-1)(1-eta)w`.
* *Unaffected.* Lemmas 2.1, 2.2, 3.2, 3.3, 4.0; Theorems 3.1, 4.1(1)–(3) and 5.1 as stated; Cor 7.1. The referee
  reproduced the record 10-point cluster (`C = 8.2386`, `lambda(Z) = 3.1127`). `F(y) = G(Z)` holds exactly in `Q(i)`
  on three clusters, Lemma 3.3 integrality holds, and (5.1) holds on random and structured supports (`reg(R_5) = 3`,
  `reg(2R_5) = 6`). No Gaussian gcd or denominator cost is dropped, no common phase is fixed, no fixed-`S` Subspace
  exceptional set is used.

### 1.2 coset (`coset.md`; `referee_coset.md`): survives; status `partial` confirmed

Setting (`round21/fakes.md` §1.1). Take a profile `(p, e, (a_x))`, `u_x = 2a_x - e`, `lambda_x = u_x - u_1`. `G(p) ⊂ T^r`
(`T = R/(pi/2)Z`) is the closed set of angle vectors satisfying every monomial gap with weights `p`. A **realisation
coset** is `C = g + H_0`, of free dimension `m`, cut out by the arc conditions. A **character** is `v in L^sat`, which
is constant on `C`. **Scenario P**: `nu(G(p) ∩ C) > 0` on window realisation cosets with `m >= cM`. **Scenario R**:
`G(p) ∩ C` is empty or consists of near-arithmetic points.

> **Round 25 coset (referee-corrected).** Partial. No growth and no uniformity.
> 1. *(proved)* **Lemma 2.1.** A point satisfying all non-character gaps for weights `p/tau^2` (`tau > 1`) and all
>    character gaps for `p` is a Lebesgue density point of `G(p) ∩ C`, on any coset.
> 2. *(proved)* **Theorem 2.2 (full torus).** `mu(G(p)) > 0` iff `5 ∉ p`, for every `Lambda`. Transplant vectors
>    (`arg beta_j`, `N(beta_j)` the preceding split prime) are slack points. This settles the main question of fakes open
>    point 4. It does **not** settle the genuine-background fibre version, since a genuine background is tight. It
>    supersedes Proposition A' as a positivity statement only; quantitatively A' is far stronger where it applies.
> 3. *(proved)* **Corollary 2.3.** `log mu(G(p)) >= -(2.54 + o(1)) r^2 log P` at `tau = 2`. This refutes the
>    **full-torus** Poisson extrapolation of fakes §5.1 for large `Lambda`, with crossover near `Lambda ~ log r`.
>    Fakes open point 3 (on realisation cosets) is **untouched**. *(Numerical correction: for primes 13..73,
>    `tau^2 = 41/37`, `n_0 = 1859`, `log mu >= -3.2e4`.)*
> 4. *(proved; Theorem 3.1 sharpened by referee Prop R1, a Borel–Cantelli argument)* **Equivalence.**
>    `nu(G(p) ∩ C) > 0` iff `C` contains a point whose characters hold for `p` and whose non-characters hold for
>    `p/tau^2`, for some `tau > 1`. Moreover `nu`-a.e. point of `G(p) ∩ C` is such a point. Haar positivity is
>    therefore existence with slack; the measure form of (Q) adds nothing.
> 5. *(proved)* **Theorem 4.1 / Corollary 4.2.** Proposition A' holds on cosets, with exact step sums in graph
>    coordinates.
> 6. *(proved)* **Theorem 5.1 / Corollary 5.2 (certificate barrier).** No certificate of the Theorem 4.1 form proves
>    positivity when `m >= cM` and `M >= M_0(c, delta)`. That form means full-circle Fubini steps with per-step
>    fibre union bounds, in any lattice or isogeny basis and in any order. Examples: `M_0 = 4.5e3, 4.6e5, 7.0e8,
>    4.9e17, 4.3e20` at `(c, delta) = (1/4,1), (1,1/2), (1/4,1/2), (1/2,1/4), (1/4,1/4)`; `M_0` is uniform in the
>    profile. Not covered: exact or non-uniform evaluation of the step measures, and localized (Lemma 2.1-type) steps.
> 7. *(proved)* **Proposition 6.1.** Transplant points with `N(beta_j) <= p_j` on `Real(delta)` are actual clusters
>    with `log N_beta <= W` and `C_beta <= C`. Seeding with such points is circular (round-21 J2). "Arithmetic" must
>    mean `N(beta_j) <= p_j`, because arguments of Gaussian integers of unbounded norm are dense.
> 8. *(proved)* **Proposition 7.1, Lemma 7.2, Corollary 7.3 (conditional P, small `m`).** Positivity holds on
>    `Delta`-translates of a Gaussian-seeded `(pi/4)`-base coset whenever `(1/2)(4 eps_*/pi)^m > H_s(rho)`, with
>    `rho = C kappa_inv e^{-W/4}` **kept explicit**. Here `kappa_inv` is the norm of a right inverse of `(lambda_x)`,
>    a determinant denominator. Such seeds force `arcsin(y_v/sqrt2) <= |v|_1 rho` for every character. The range
>    `m <~ M^{(2-delta)/5}` is a sketch and also needs `log kappa_inv + log b* = o(W)`. It fails by a power of `M` for
>    `m >= cM`.
> 9. *(exact)* `G(p) ∩ line = ∅` for 61 explicit lines in `T^8` (primes 13..73), with certified overlaps `>= 4.4e-6`.
>    10 of them pass truncated character tests; the `eta`-line passes to `T = 40`. These lines are not realisation
>    cosets and not in the window regime, so they are evidence for heuristic A only.
> 10. *(open)* **(Q_seed)**, with the corrected heavy condition `arcsin(y_v/sqrt2) <= |v|_1 C kappa_inv e^{-W/4}` and
>    heavy mass `H_s(rho)`. It is a **sufficient** route to scenario P, not a reduction of the `A_0` window question.
>    By round-21 J2 neither answer gives growth.

*Corrections adopted.*
* **Major (owner's denominator constraint).** (Q_seed), §0 item 6 and the reading of Cor 7.3 wrote the heavy mass as
  `e^{-W/4+O(r)}` and the character condition as `w(v) >= W/2 - O(log W)`. Both drop `log(C kappa_inv)`, which
  Hadamard bounds only by `O(r log(r e_max))`, the same order as `W` at barrier scale. The referee's example is
  unit-triangular differences with `-1` above the diagonal (`e_j = 2`). Their naive right inverse has norm
  `2^{Theta(r)}`; judge check J2(c) confirms that the exact inverse has entries `2^{j-i-1}` and row `l^1` norm
  `2^{r-1}`. The referee did not check whether the `mod pi/2` freedom in `eta` lowers the norm; nothing bounds it.
* **Minor.** "Reduces to (Q_seed)" is one-directional. "Fails for the coset `T^r`" refutes only the full-torus
  extrapolation, because `T^r` is a realisation coset only when `M = 1`. Theorem 2.2 settles only the main part of
  open point 4. "Also necessary" in §6 needed Prop R1. The slack for 13..73 was wrong. "Arithmetic" needed
  `N(beta_j) <= p_j`. The count of lines passing the character tests is 10, not 1. The barrier is scoped to Theorem
  4.1-form certificates.
* *Unaffected.* Every theorem as stated. Independent code reproduced the 61 covers, the Cor 4.2 formula (relative
  error `<= 1.4e-7`, truncation), the `M_0` table, the transplant density table, and Lemma 7.2 in exact Gaussian
  arithmetic over all 256 shifts (minimum margin 0.026). The common phase `theta_0` is kept exactly. No fixed-`S`
  Subspace input is used.

## 2. Does the intra-cluster multi-point argument give a growth-relevant gain?

**Answer: no. It reduces to route 5 in the growth regime and to one-point caps in the lacunary or cube regimes. The
only place where sharing beats a one-point cap is route 5's own LP supercriticality.**

| regime of the auxiliary data | what the intra-cluster method is | label |
| --- | --- | --- |
| any `F`, non-archimedean side | `(TL_M)` at the one point `Z` for `G = Phi^*F`; same value, weaker width and order | proved (Thm 3.1) |
| star blocks, slab boxes | route 5 exactly (forms and (I1)) | proved (Thm 4.1) |
| lacunary block degrees (what Roth-type zero lemmas need) | the leading-block **one-point** comparison, up to `O(kappa omega/rho)` | proved under `R'_1 >= rho d_1` (Prop 6.1) |
| `Z` balanced on all cuts (cube class) | capped at the one-point `15/4`, uniformly in `G`, for Taylor/Dyson bounds | proved (Thm 5.1 + Lemma R2) |
| 5-wise uniform designs in route 5's regime | **beats** the one-point cap, but its output is `F(y) = 0` (exact relations) | proved (referee Prop R1) |
| growth regime `W = O(M log M)` | LP-supercritical; exact relations among star subtuples; needs route 5's ZE | proved (Cor 7.1); ZE open |

*Judge remarks (unrefereed).*
* **J-a (exact; `judge_checks/j2_scale.out`). Neither the cube data nor the design data are in the growth regime.**
  * Cube data have `r >= 2^{M-1}-1` primes, so `W >= (2^{M-1}-1) log 5`. That is `822` at `M = 10` against
    `M log M = 23`.
  * The extended-BCH data of Prop R1 have `r = 2^{2t+1}-2 >= 2M^2 - 2`. At `t = 13`, `M = 8189`,
    `w_0 >= 57 log 2`, the bound is `W >= 5.30e9 >= 7.1e4 · M log M`.

  Theorem 1.1 is far from tight on both families. So item 3 of §1.1 is a barrier against the **uniform** bound only,
  and item 4 is a statement about the method's power, not about growth. The intra method touches growth only through
  item 6, which is route 5.
* **J-b (logical). The anisotropic Dyson weights, the one ingredient route 5 lacks, cannot remove the zero-estimate
  obstruction.**
  * In route 5's regime the comparison is already supercritical: for every `F` in the forced space `V`, it concludes
    `F(y) = 0`.
  * A stronger archimedean bound can only enlarge `V`.
  * A growth bound needs some `F in V` with `F(y) != 0` at the **actual** point. That is a zero estimate at a point
    whose base locus depends on the cluster's profile, the template-dependent step the owner forbids treating as
    uniform.

  Whether a larger `V` makes that zero estimate easier is part of route 5's open ZE question. The intra write-up does
  not address it.
* **J-c (exact; `judge_checks/j1_design.out`).** The design facts behind Prop R1 extend beyond the referee's `t = 3, 4`
  to the extended BCH design at `t = 5` (`M = 32`, `m = 11`):
  * all `201376` 5-subsets of columns are independent over `F_2`;
  * every one of 202 sampled 5-subtuples is exactly balanced;
  * `Z` carries 1024 of the `2^31 - 1` cut classes.

  The route-5 regime condition `M - 1 >= 300 ln 4r` first holds at `t = 13` (`8188 >= 6030.4`). This confirms the
  referee's threshold.

**Verdict on the lens.** Rounds 24 and 25 together close the multi-point Thue/Siegel/Dyson family on `T_5`:
* disjoint primes: one-point caps (round 24 Prop 4.3);
* shared primes, intra-cluster: route 5 (this round);
* lacunary degrees: one-point caps (Prop 6.1).

The only part not closed is round-25 intra open point 2: shared-prime tuples that are not images of one balanced point,
for example subtuples of different circles. Those are not data of a single cluster. Nobody has proposed an argument
that would import them, and none is in sight.

## 3. Status of the low-height window after the coset result

**Status: OPEN, unchanged in truth value.** No restricted growth theorem holds inside `P_+(N') <= M^{2-delta}`, and
the window is not capped for `A_0`. What changed is the map around (Q).

| item (round 21 numbering) | before round 25 | after round 25 |
| --- | --- | --- |
| `A_0` question (Q) on realisation cosets, `m >= cM` | open; the recommended route was Proposition A' on cosets (10–15% map value) | **still open.** A' on cosets is proved (Thm 4.1). Certificates of that form provably **cannot** reach the window for `M >= M_0(c, delta)` (Cor 5.2). The round-21 route is closed for its certificate form |
| Haar positivity versus existence | Theorem A: Haar probability `<= exp(-c M^{delta/2})`; positivity undecided | **equivalent** to a slack point with valid characters (Thm 3.1 + Prop R1). Theorem A's smallness never decides positivity. The measure form of P is the existence form |
| full-torus picture | Prop A' up to `Lambda ~ (1/4) log P`; Poisson extrapolation `exp(-c(Pi-1))` proposed | positivity for **every** `Lambda` when `5 ∉ p` (Thm 2.2). The full-torus Poisson extrapolation is false at large `Lambda` (Cor 2.3) |
| fakes open point 4 (rigid prime sets) | open | main question answered: no prime set without 5 is rigid. The genuine-background fibre version is open |
| fakes open point 3 (doubly exponential upper bound on realisation cosets) | open | **open** (only the full-torus extrapolation is refuted) |
| seeding by arithmetic points | Lemma T: arithmetic models are actual clusters | Prop 6.1: arithmetic slack seeds on realisation cosets are actual clusters (`N(beta_j) <= p_j`), so circular |
| scenario P | unproved | proved conditionally, for small `m` only (Cor 7.3; range `m <~ M^{(2-delta)/5}` a sketch needing `log kappa_inv = o(W)`). Fails by a power of `M` at `m >= cM` |
| scenario R | unproved | unproved. Exact emptiness on 61 lines in `T^8` (not realisation cosets, not window-scale) is evidence for heuristic A only |
| row 7, square-root-quality inputs, primes in `[cM log M, M^{1+eta}]`, pinned `m = 0` | uncapped / untouched | unchanged |

*Judge reading.* By Prop R1, the `A_0` window question on a realisation coset is now exactly this:

> **(Q')** Does `C` contain a point satisfying the character gaps for `p` and the non-character gaps for `p/tau^2`
> (some `tau > 1`) that is **not** within the Lemma-2.1 neighbourhood of a transplant with `N(beta_j) <= p_j`?

A "yes" on every window coset with `m >= cM` caps `A_0` and the row-3 profile inputs there (scenario P). A "no" makes
`A_0` lossless there (scenario R). By round-21 J2 neither gives growth: P caps one input class, and R still needs a
bound on `W_min^win` from inputs not yet in hand. Round 25 therefore resolves the round-21 recommendation as a
**negative** for the proposed tool, with map value. The window is the general problem with smaller primes, as round
21 concluded, and the coset analysis gives no reason to revise that.

## 4. Consequences for growth and for the uniform bound

Only growth (`M = o(log R/loglog R)` for every `R`) or uniform progress (`M(C) < infinity`) counts. **Neither moved.**

* **Growth.** No new route, and no route shortened.
  * The intra line merges into route 5. The missing step is still route 5's zero estimate at the actual point
    (ZE-col\*; MIN^pair or the sharper (R_H) of round 16), which has been open since round 11.
  * The coset line shows that the window's `A_0` question (Q') does not decide growth (round-21 J2), and the
    proposed tool for deciding it has a proved certificate barrier.
  * The `(2/3 + eps)` bound stands.
* **Uniform bound.** Question B (`(L_5^circ)`, or `A_5(delta)` with `delta < 1/30`) is untouched. Two pieces of map
  value:
  1. A further proved barrier. Intra-cluster Liouville with Taylor/Dyson archimedean bounds cannot prove the uniform
     bound, nor anything below `M ≈ log_2 log R`, because cube-balanced clusters are invisible to it (§1.1 item 3).
     This is in line with round 24: pointwise Liouville-type arguments on `T_5` are capped at `15/4`.
  2. A correction to the record. Round 24's suggestion that shared primes give "the same cap" is false for
     intra-cluster data (§1.1 item 4). The disjoint-primes hypothesis in the expert note's anisotropy paragraph is
     therefore load-bearing and must stay explicit. The failure is through exact relations, which leads back to route
     5's zero estimate. It does not open a route to Question B.
* **Owner constraints.** Two violations were found and corrected by the referees:
  * a template-dependent threshold presented as uniform (intra M2), repaired for Taylor bounds by Lemma R2;
  * a dropped determinant denominator `log kappa_inv` (coset major), now explicit.

  No fixed-`S` Subspace exceptional set is used, no Gaussian gcd cost is discarded (multi-point Lemma K1 tracks it),
  and no drifting common phase is set to zero (`theta_0` is kept in `Real(delta)`, and the intra evaluations are of
  degree-0 characters).

## 5. Judge's own checks (unrefereed)

All runs use `prlimit --as=3500000000`, one process at a time. The scripts are written from scratch and import nothing
from the authors' or referees' code.

| id | script → output | claim checked | result |
| --- | --- | --- | --- |
| J1 | `judge_checks/j1_design.py` → `j1_design.out` (exact, `F_2` and `GF(2^t)` arithmetic) | extended double-error-correcting BCH columns `(1, x, x^3)`, `x in GF(2^t)`: any 5 independent, all 5-subtuples exactly balanced (15 cuts × `2^{m-5}` primes, gcd class `2^{m-5}-1`), cut classes carried | `t = 3`: 56/56 independent, 0 balance failures, 64 of 127 cuts. `t = 4`: 4368/4368, 0 failures, 256 of 32767. `t = 5` (**new**): 201376/201376 independent, 0 failures on 202 sampled subtuples, 1024 of `2^31 - 1` cuts |
| J2(a) | `judge_checks/j2_scale.py` → `j2_scale.out` (exact rationals) | referee Prop R1(d): `(1-2^{-m})(1+eta)/(1-eta) < 1-eta` for `m = 2t+1`, `eta = 2^{-(2t+4)}` | holds for every `t = 3..30` |
| J2(a') | same (floats; integer comparison with a margin `>= 2000`) | route-5 regime `M - 1 >= 300 ln 4r` with `r = 2^{2t+1} - 2`, `M = 2^t - 3` | fails for `t <= 12`; holds from `t = 13` (`8188` vs `6030.4`) |
| J2(b) | same | scale of the barrier data (§2, J-a) | R1 data: `W >= 5.30e9 = 7.2e4 · M log M`. Cube data: `W >= (2^{M-1}-1) log 5`, e.g. `822` at `M = 10`, `8.4e5` at `M = 20` |
| J2(c) | same (exact rationals) | referee coset finding 1: `U = I - (strict upper ones)` has inverse entries `2^{j-i-1}` | true for `r = 5, 10, 20`; max row `l^1` norm `2^{r-1}` |
| J2(d) | same | Lemma R2's constant: `1 + log m <= (1.39 - (log 4 - 1)) m` | true for `m = 1..10^6`, and for all `m >= 1` by monotonicity of `1.0037 m - 1 - log m` |

The logical remarks J-a, J-b (§2) and the reading (Q') (§3) are mine and unrefereed.

## 6. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing, (b) fake-violated input, (c) movement of growth
or uniformity. The overall score is a judgement that weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | overall | reasons |
| --- | --- | --- | --- | --- | --- |
| **intra** | 6 | 3 | 2 | **3** | **Correct core.** Thm 3.1, Thm 4.1 and Thm 5.1 as stated were re-derived and checked exactly on three actual clusters; the bookkeeping is exact. **Deductions.** The headline item 4 was reversed by the referee (Prop R1). The cube-barrier threshold was template-dependent (repaired only for Taylor). Its "every inference / any zero estimate" scope was too wide. The zero-lemma dichotomy is heuristic. **(b) = 3.** The input is route 5's forced-divisibility system: almost every fake violates its exact relations, but the residual step (ZE) is the open problem. Nothing is added beyond route 5. **(c) = 2.** Proved map corrections: the intra/route-5 identification closes the round-24 multi-point lens, and the cube barrier is a uniform-bound barrier; via the referee, round 24's open point 1 is answered negatively. No growth movement. |
| **coset** | 7 | 1 | 2 | **3** | **Correct core.** Every theorem survives as stated, and the referee added the equivalence Prop R1. Exact covers and constants reproduced. **Deductions.** One major point: the determinant denominator was dropped in the headline open question (owner constraint). Several overstatements: "reduces to (Q_seed)", "fails for `T^r`", "supersedes A'", "one line passes every character test", and a wrong slack constant. **(b) = 1.** The work is on the fake-existence side. The only fake-relevant necessary condition (Gaussian seeds force heavy characters) concerns one seeding route, not an input that fakes violate. **(c) = 2.** It closes round 21's recommended route negatively for its certificate form. It settles the main question of fakes open point 4, refutes the full-torus Poisson extrapolation, and turns Haar positivity into existence with slack. No growth movement, by J2. |

Both are honest pieces of work whose headlines did not fully survive. Neither moves growth or uniformity, which caps the
overall score at 3.

## 7. Proposed additions to the expert note (text only)

These are text proposals only; `expert_note/expert_note.md` was not edited. Each uses the note's labels.

**7.1 Append to the "Anisotropy" paragraph of §3, after "...with disjoint primes nothing in $\mathcal L$ beats $15/4$
($\delta=1/30$)...":**

> The disjoint-primes hypothesis is essential (round 25 `intra.md`, `referee_intra.md`, **proved**). The
> $\binom M5$ gcd-reduced 5-subtuples $y_S$ of one cluster are the projections $\mathrm{pr}_S(Z)$ of a single point
> $Z\in T_M(\mathbb Q)$, $T_M=(R^1_{\mathbb Q(i)/\mathbb Q}\mathbb G_m)^M/\mathrm{diag}$. Every multi-point
> Liouville/Dyson comparison among them, with toric forced divisibility across shared primes, is one-point toric
> Liouville on $T_M$ at $Z$. Conversely, a star partition $T_M\cong T_5^{(M-1)/4}$ identifies route 5's
> forced-divisibility forms with intra-cluster forms. Two results follow.
> * If $Z$ is balanced on all $2^{M-1}-1$ cuts of $[M]$, nothing in this class with a Taylor-type archimedean bound
>   beats $15/4$, uniformly in the auxiliary form.
> * For 5-wise uniform cut designs (e.g. extended BCH, $M\ge8189$) every 5-subtuple is exactly balanced with
>   $\mu<4(1-\eta)$, yet route 5's LP supercriticality forces $F(y)=0$.
>
> So shared primes do beat the one-point cap, but only by producing exact algebraic relations. Using those relations
> requires route 5's open zero estimate at the actual point. This gives nothing towards Question B.

**7.2 Replace the clause "and the low-height window $P_+(N')\le M^{2-\delta}$ remain uncapped; ... (`round21/README.md`,
status **open**)" in §1 by:**

> and the low-height window $P_+(N')\le M^{2-\delta}$ remains uncapped (`round21/README.md`, status **open**). Round 25
> (`coset.md`, `referee_coset.md`, **proved**) reduced the window's question for pair and monomial-gap inputs to pure
> existence.
> * On a realisation coset, the gap-admissible set has positive Haar measure iff the coset contains a point that
>   satisfies the character gaps and satisfies the non-character gaps with uniform slack.
> * Slack arithmetic points (transplants with $N(\beta_j)\le p_j$) are actual clusters.
> * Sequential union-bound certificates (Proposition A' on cosets) provably cannot reach the window when the free
>   dimension is $\ge cM$.
> * On the full torus the gap set has positive measure for every prime set without 5, at every height.
>
> Neither answer to the existence question gives growth.

**7.3 Add to §5, item 3 ("Auxiliary polynomials and single-place Diophantine approximation"), a fifth bullet:**

> * Multi-point (Thue/Siegel/Dyson) arguments on the 5-subtuples of one cluster are one-point toric Liouville on
>   $T_M$ at the cluster point, which is the determinant method with forced divisibility ("route 5"). In the growth
>   regime they force exact algebraic relations and stop at a zero estimate at the actual point. With lacunary degrees
>   they reduce to one-point caps (round 25 `intra.md`, **proved**; the comparable-degree dichotomy is **heuristic**).

**7.4 Add to the pointers table of §6:**

> | round 25: intra-cluster multi-point Liouville; Prop A' on cosets, Haar positivity = existence with slack | `round25/README.md`, `round25/intra.md`, `round25/coset.md` (with `referee_*.md`) |

## 8. Single most promising next step

**Decision: stop internal rounds.** This confirms rounds 19, 20, 22, 23 and 24.
* Round 25's two lenses were the last open items flagged by rounds 21 and 24 (round-21 §7; round-24 open point 1).
* Both are now decided at the level of the map. The intra line reduces to route 5, and the coset route is closed for
  its certificate form.
* Neither produced a new input that actual circles satisfy and fakes provably violate.

**The single most promising next step is to fold §7 into the expert note and put Questions A and B to experts.**
Question B (`(L_5^circ)`, or `A_5(delta)` with `delta < 1/30`) is the sharpest statement in the record that would
give the uniform bound. Both round-25 lines, and every multi-point Liouville variant, are now proved to stay below it.

*If one more internal round is run anyway*, the only internal item on which several routes converge is **route 5's
zero estimate at the actual point**: (R_H) of round 16, equivalently ZE-col\*. The intra line (Cor 7.1) and route 5
(rounds 10–17) both stop exactly there. A useful formulation would avoid the owner's forbidden template-dependence:
show that, for the cluster's own regime box and gap sphere, the forced space `V` contains a form that does not vanish
at `Z`. Round 16 proved that degree/variance/rank-blind anti-concentration tools are vacuous here, so any attempt must
use `u in G` (the gap sphere).

*Honest probabilities (judge's estimates; not derived).*

| outcome of one more internal round on (R_H) / ZE-col\* | probability |
| --- | --- |
| a proved map correction or barrier around (R_H) | about 15% |
| (R_H) or ZE-col\* proved at the actual point (by round 16, (R_H) gives MIN^pair and `M = O(loglog R)` on all norms) | about 1% |
| unconditional growth `M = o(log R/loglog R)` for every `R` | < 0.5% |
| the uniform bound, `(L_5^circ)` or `A_5(delta < 1/30)` | < 0.1% |
| a restricted growth theorem inside the window `P_+(N') <= M^{2-delta}` | about 1% |

No probability is offered for expert circulation. Question B is a McKinnon–Roth/Vojta-type statement, and no
unconditional theorem of that kind is known (expert note §3), so the value of circulating it lies in an expert's
judgement, not in a forecast.

## 9. Files

* Write-ups: `intra.md`, `coset.md`; notes `intra_notes.md`, `coset_notes.md`; checks `intra_checks/`, `coset_checks/`.
* Referee reports: `referee_intra.md` (notes `referee_intra_notes.md`, checks `referee_intra_checks/`),
  `referee_coset.md` (checks `referee_coset_checks/`).
* Judge checks (unrefereed): `judge_checks/j1_design.py` → `j1_design.out`; `judge_checks/j2_scale.py` →
  `j2_scale.out`.
* This report: `README.md`.
