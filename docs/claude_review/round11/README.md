# Round 11: Hypothesis ZE for the determinant route, after refereeing

## 0. Bottom line

**No growth improvement for actual circles and no uniform bound was proved this round**, either
unconditionally or from a standard conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still
the best bound.

**Hypothesis ZE is open: neither proved nor refuted in the regime.**
* The evidence for it is weak: exact theorems on two zonotope-free toy models that are LP-subcritical.
* There is no evidence against it in the regime.
* It fails exactly where it was expected to fail: for underdetermined systems and thin families.

What changed (refereed, except the items marked J, which are this README's unrefereed findings):
1. **Every stronger or delta-free reformulation of ZE proposed this round is false.** Only the local ZE-tor and
   the restricted ZE-col\* (both open) survive.
   * *ZE-braid* (compute.md) fails at Sylvester orders `n = 2^k >= 2^16`, because of `log_2 n`-column rational
     subtori in the base locus (referee).
   * The delta-free *SL(n)* (theory.md) fails because of twin and duplicate columns of the barrier profiles, and for
     every odd `m` (referee).
   * Its repair *SL\** fails at Sylvester orders `n = 2^k >= 2^13`, because each column is two other columns plus a
     quarter-size indicator (**J1**, §5.1).
   * *ZE-col* fails for thin families and must be restricted to ZE-col\* (referee).
   * Common cause: the cheap rational subtori of the base locus, which are rich for Hadamard-type profiles. They
     contain all profile-level data of the actual point, such as its `p`-adic valuation vectors, but not the point
     itself, which is nondegenerate.
2. **The round-10 finite test cannot see the regime.**
   * At `n - 1 = 2, 3, 4` every slab-box examined is settled by a zonotope (product, i.e. character) certificate.
     Such certificates never exist for actual clusters in the regime (compute.md Props 2.1-2.2).
   * The first zonotope-free models (`H_8`, `T = 1`) have toric base loci (Theorem 5.2, exact). But they are
     LP-subcritical and lack multi-column strata.
   * A zonotope-free model with `dim V = 1` (skew double core, `q = 7`) has nondegenerate zeros (Prop 6.2). **ZE is
     an overdetermination statement.**
3. **Conditional results.**
   * SL and SL\* imply `n = O(loglog R)`, but both hypotheses are false, so these implications are vacuous.
   * New: ZE-col\* implies `M = O_C((loglog R)^2)` **on all norms** (subcluster.md Thm 5.4, corrected).
   * Unconditional new tools: the dichotomy on general supports with `V_H` of relative codimension `O(1/n)`
     (theory.md Lemmas 1.1-1.2); the signed mean-width bound `reg(S) <= E_eps h_S(<eps,.>)` (subcluster.md
     Thm 4.4); and "blind to the universal family implies `M < (300 ln 4r + 2)^2`" (Cor 4.2).
4. **Box size (J2, J3).**
   * With round 10's box `T = n`, the leading-order (tangent-cone) criterion can never certify a near-extremal
     cluster. A Chebyshev bound shows it needs `W >~ n^2` (§5.2).
   * Route 5 works equally well at `T = 1`. Theorem 4.1 holds from `n ~ 2900` (`r = n`), and the order can be taken
     as small as `m = ceil(4w) = Theta(sqrt(n log n))`.
   * There the actual position (`e^(-Theta(n log n))`) lies deep inside the Taylor regime. ZE then reduces to a
     **real tangent-cone statement** (ZE-TC) plus a near-strata separation that the pair inequalities supply for
     large `n` (heuristic).
   * In the two `T = 1` models the tangent-cone statement holds (J4).

**Scores** (rubric of `round10/README.md` §1):

| write-up | (a) correctness | (b) fake-violated input | (c) movement | **overall** | reason |
|---|---|---|---|---|---|
| **subcluster** | 7 | 3 | 3 | **4** | Most results survive. New conditional theorem on all norms (`(loglog R)^2`, ZE-col\*), unconditional Cor 4.2, new Thm 4.4. Deductions: the hypothesis had to be restricted (ZE-col false for thin families), ZE-tor had to be made local, two claims were withdrawn, and the "empty below 2137" computation was wrong. |
| **compute** | 8 | 3 | 2 | **4** | All exact results survive independent recomputation. Shows that the round-10 test sizes are vacuous (zonotopes), gives the first zonotope-free models with an exact toric base locus, and shows that ZE needs overdetermination (Prop 6.2). Deductions: the recommended target ZE-braid is false in the regime, and the `n = 8` models are LP-subcritical, so they were overstated as evidence. |
| **theory** | 6 | 1 | 1 | **3** | Every proof is a correct implication. Lemmas 1.1-1.2 are useful, and the `p`-adic mechanism (Thm 3.1) is new and correct. But the headline hypothesis SL(n) is false, its repair SL\* is false too (J1), Thm 3.1 is empty for odd `m` and for Sylvester-Hadamard profiles, the evidence models had the wrong level structure, and (2.1) is vacuous at `T = n` (J2). |

**Recommended next step** (§8): run route 5 at `T = 1` with `m = ceil(4w)`, prove the reduction of ZE to the
quantitative real tangent cone, and decide the qualitative tangent-cone statement ZE-TC for slab-cut ternary
cubes. ZE-TC says that common real zeros of the `Phi_c` lie in spans of cheap families. Start with few cuts and
with the Sylvester-Hadamard family.

## 1. Setting and notation (from round 10)

Round 11 had three write-ups, each attacked by an adversarial referee whose report sits next to it:
* `compute.md`: the finite test of ZE proposed in `round10/README.md` §8.4;
* `theory.md`: attempts to prove ZE, and a delta-free substitute;
* `subcluster.md`: what intersecting the rigidity relations of many sub-clusters can force.

Conventions and labels are those of `round10/README.md`. *Proved* means proved with all quantifiers in a
refereed write-up. *Exact* means integer or exact-rational computation. *Evidence* means floating point or a
finite experiment that no proof uses. *README remark* marks material added by this judge. Nobody but the judge
has checked it, so it is **unrefereed**; the judge's scripts are in `judge_checks/`.

**The objects** (`round10/determinant.md`, cited as **[D]**).
* *Cluster.* `n` points on `|z|^2 = N`, `W = log N`, arc constant `C`.
* *Profile.* Relative angles `delta_x`, exponent columns `a_p`, centred columns `beta_p`, `r` split primes.
* *Slab-box* of [D] Thm 2.1. With box size `T` (round 10: `T = n`) and `w = T sqrt(8(n-1) ln 4r)`:
  `A = {k in Z^n_0 : |k_x| <= T (x < n), |<k,beta_p>| <= w for all p}`. The order `m` is the largest integer
  with `8 binom(m+n-2, n-1) <= (2T+1)^(n-1)`, and `D = 2(n-1)T`.
* `V(A,m)` is the space of measures on `A` that annihilate all polynomials of degree `< m`. `V_H` is the span of
  its integer members with `|c|_inf <= H`.
* *Base locus.* `Bs(A,m)` is the common zero set of `P_c(w) = sum_k c_k w^k` on `(C^*)^(n-1)`, and
  `Z_P = {delta : sum_k c_k e^(i<k,delta>) = 0 for all c in V_H}`.
* *Tangent-cone forms.* `Phi_c(u) = sum_k c_k <k,u>^m`.
* *Cheap family.* Integer directions `v_j` with `sum_j h_A(v_j) < m`. Their real span always lies in
  `Bs`/`Z_P` and consists of degenerate points (theory.md Prop 2.2).

**What round 10 established.**
* **Thm 4.1 (proved, unconditional).** Every actual cluster with `n - 1 >= 300 ln 4r` and
  `W >= 16 ln n + 8 log^+ C + 14` has `delta in Z_P`.
* **Hypothesis ZE(n).** No nondegenerate `delta` with `max|delta_x| <= 1/(4D)` lies in `Z_P`. Here
  *nondegenerate* means `sum v_x delta_x != 0` for every `0 != v in Z^n_0`.
* **Thm 5.7 (proved, conditional).** ZE(`n'`) for all `n' >= n_1` implies `n = O(log W) = O(loglog R)` on
  squarefree norms.
* The round-10 README recommended deciding ZE, starting with a finite test (§8.4) at `n - 1 = 2, 3, 4`.

## 2. Referee-corrected statements: compute (finite ZE tests)

Write-up: `compute.md`. Report: `referee_compute.md`. Scripts: `compute_checks/` and
`referee_compute_checks/`. The referee used an independent kernel, factorisation and Groebner engine (flint
instead of sympy). **Status: evidence. The exact results survive; two interpretive claims do not. No growth
bound, no proof of ZE.**

1. **Zonotope certificates (Prop 2.1, proved).** Suppose `A` contains a lattice zonotope `k_0 + sum_i [0,L_i] v_i`
   with `sum L_i >= m`. Then an explicit `c in V_H` (a product of finite differences, `|c|_1 <= 2^m`) satisfies
   `P_c(w) = w^(k_0) prod_i (w^(v_i) - 1)^(alpha_i)`. So `Bs(A,m)` lies in finitely many subtori `{w^(v_i) = 1}`,
   and every point of it is degenerate. Write `zeta(A)` for the largest such `sum L_i`. Then
   `zeta(A) <= rho(A) := max{j : V(A,j) != 0}`.
2. **Zonotopes never exist in the regime (Props 2.2 and 6.1, proved).** Let an actual cluster with nondegenerate
   positions satisfy [D] (4.1). Then its own slab-box has `zeta(A) < m`, because a zonotope would violate the
   (corrected) character inequality. The same holds for every support `A' ⊂ B_T` with
   `sum_p omega_p h_(A')(beta_p) <= m'/2`. The hypothesis `A' ⊂ B_T` was added by the referee.
   * So `zeta(A) < m <= rho(A)` ("faithful") is a **necessary** condition for a model to resemble a regime
     support. It is **not sufficient** (items 4 and 6).
   * The remark "every model with `m <= 2T` is unfaithful" needs `L >= T`.
3. **The round-10 test sizes are settled by zonotopes.**
   * `n - 1 = 2`: every slab-box is a lattice zonotope (a rectangle, parallelogram or hexagon) with
     `zeta = rho`. This is proved; the referee's proof needs no condition (ii).
   * `n - 1 = 3`: all 355 distinct admissible slab-boxes from full-rank sign profiles (one common `L`, `T <= 4`,
     `L <= 3`, condition (ii')) have `zeta = rho`. Exact and reproduced: 980 / 339 admissible triples under
     (ii') / (ii).
   * `n - 1 = 4`: every sampled set has `zeta = rho` (the author's 74, the referee's 75; exact, sample only).
   * So in every case computed, the README §8.4 test is decided by character-type certificates and says nothing
     about the regime.
4. **The first zonotope-free models (exact).** Take Sylvester `H_8`, base point 8, `T = 1`. Prop 3.1 (proved)
   gives `zeta <= ((n-1)/n) max_p h_A(beta_p)` for Hadamard profiles.
   * `M8 = (L = 2, m = 8)`: `|A| = 1037`, `dim V = 42`.
   * `M8' = (L = 3, m = 11)`: `|A| = 1877`, `dim V = 7`.
   * Both satisfy `zeta < m = rho`, `V_H = V(A,m)` and (ii'). Both **fail** the round-10 condition (ii)
     `h_A(beta_p) < m` (`8 = 8`, `12 > 11`), and miss the regime margin `h_A <= m/2` by a factor `>= 2`.
   * For `H_8`, condition (ii) cannot be met at `T = 1` (any `L`) nor at `T = 2`, `L <= 2` (referee).
5. **Theorem 5.2 (exact; independently confirmed).** In `(C^*)^7`:
   * every point of `Bs(A_3, 11)` (model `M8'`) has `w_x = 1` or `w_x = w_y`;
   * every point of `Bs(A_2, 8)` (model `M8`) lies on `{w_x = 1}`, on `{w_x = w_y}`, or on one of finitely many
     explicit subtori `{w^a = w^(a')}` through `e`. The list of these depends on the basis: six in the author's,
     seven in the referee's.
   * So no nondegenerate `delta` lies on `Z_P` for these two `(A,m)`.
   * Method: each basis polynomial is a monomial times two-term factors times irreducible factors that are
     **forms in `x = w - 1`**. The ideal of those forms is saturated by the collision hyperplanes (Bayer's
     criterion) and becomes the unit ideal.
   * This is a model statement with ad hoc `(T, L, m)`. It is not ZE(8), which is vacuous, and not ZE at
     [D] Thm 2.1's parameters, where the slabs do not cut at `n = 8`.
   * The remark that each single `P_c` has nondegenerate *real* zeros is unproved. Floating point is consistent
     with it.
6. **Prop 6.2 (proved from exact data; independently confirmed).** Take the skew double core at `q = 7`, `T = 1`,
   `L = 2`, `m = 8`. Here `zeta <= 7 < 8 = rho` and `dim V = 1`; `c` is even, and `Phi = 8! g` with `g`
   irreducible over `Q` and changing sign. Then `Z_P` contains nondegenerate real points arbitrarily close to
   `e`. **So `zeta < m` alone does not imply ZE-real; ZE depends on the linear system being large.**
7. **ZE-braid is false in the regime (referee, proved).** ZE-braid (§5.3 of the write-up) was the recommended next
   target: realisable symbols have no common zero off the braid arrangement. The referee's counterexample:
   * Take the Sylvester-Hadamard profile of order `n = 2^k >= 2^16` with [D] Thm 2.1's slab-box.
   * The `k` coordinate-hyperplane columns cost `k * 2 floor(w/2) < m` (for example `2.682e9 < 3160214892` at
     `n = 2^16`).
   * So by [D] Prop 5.2 (with `kappa = 0`), a `k`-dimensional rational subtorus lies in `Bs(A,m)`, and its
     generic points have no collisions.
   * These points are degenerate, so ZE survives. The braid pattern at `n = 8` is an artefact of `m/w ~ 2`, where
     at most one column subgroup is affordable.

*README remark (checked).* My `judge_checks/sl_star_hadamard_out.txt` recomputes `m(2^16) = 3160214892` (60-digit
log-gamma) and the same inequality.

**Reading.** The two computable zonotope-free models have toric base loci, while single forms have non-toric zero
sets. This is consistent with ZE, but it is weak evidence: both models are LP-subcritical, and they lack the
multi-column cheap subtori that the regime has.

## 3. Referee-corrected statements: theory (delta-free SL route)

Write-up: `theory.md`. Report: `referee_theory.md`. Scripts: `theory_checks/` and `referee_theory_checks/`.
**Status: the headline does not survive.** The proposed delta-free replacement SL(n) is false, so Theorem 3.3 is
vacuous. Every proof is correct as an implication. ZE is not proved.

1. **Lemma 1.1 (proved).** `codim(V_H(S) in V(S,m)) <= (U ln 2 + E_S ln(UB)) / ln(2H+1)`.
2. **Lemma 1.2 (proved).** The dichotomy of [D] Thm 4.1 holds on *any* support `S ⊂ B_T` whose cost is at most
   `3mW/8`, with heights up to `H = e^(mW/40)`, when `W >= 27 log^+ C + 40`.
   * With that `H`, `V_H` has relative codimension `O(1/n)` in `V(A,m)`. This answers the round-10 worry that
     `Z_P` might be cut out by only half of `V(A,m)`, but only for ZE.
   * At the height `H_0 = (|S| n^(m-1))^2` used in SL, the codimension bound is only about `E/2`.
3. **Props 2.1-2.5 (proved; wording fixed).**
   * The tangent cone of `Z_P` at `e` lies in the common real zeros of the `Phi_c` (the degree-`m` least space).
   * The spans of cheap families lie in the base locus and are degenerate.
   * Cube slices reduce to products, i.e. to characters.
   * The leading-order criterion (2.1) of Lemma 2.4 provably fails at nondegenerate points within relative
     distance `~ |u|/(10 n^(3/2))` of a toric stratum. There every Taylor order is `O(eta)`; "inconclusive" refers
     to the sufficient condition (2.1).
   * README §5.2 shows that with round 10's `T = n`, (2.1) in fact fails *everywhere* for near-extremal clusters.
4. **Theorem 3.1 and Cor 3.2 (proved, unconditional, delta-free).** Widen the slab of one light split prime `p`
   (`log p <= W/8`) to `|<k,beta_p>| <= m/e_p`; call the result `A_p`. Then for clusters in the regime no vector
   of `V_0(A_p) = V_(H_0)(A_p)` has a unique `p`-adic minimiser (`UM_p`). The proof uses the exact `pi`-adic
   valuation of `F_c(z)`. The referee's qualifications:
   * **Lemma R2 (referee, proved; exact checks on 60 models).** On `A_p` the integer `l_p = <k,a_p>` takes at
     most `2 floor(m/2) + 1` values.
     * If `m` is odd, no `c` has `UM_p`, for **every** profile. Theorem 3.1 is then empty.
     * `m(n)` is odd for 4076 of the 8137 values `n in [64, 8200]`, including `2137, 3068, 3528, 4096, 8192`.
     * If `m` is even, `UM_p` can only occur on the bottom level `l_p = -m/2`.
   * So "the level obstruction disappears on `A_p`" is false.
5. **Hypothesis SL(n) is false (referee Prop R1, proved; exact on `q = 7..23`).**
   * Duplicate or twin columns (`a_p - a_q` supported on at most two non-base points) keep `l_p` within an
     interval of length `<= w + 4T < m`. The level-sum lemma then excludes `UM_p`.
   * The round-5 Paley fake profile (`M >= 3528`) and the skew double core (`M >= 3068`) are admissible. They
     fail SL at all but at most 1 (resp. 2) columns, at every Paley order in the regime. SL(n) needs 8.
   * A Hadamard profile with doubled columns fails at all columns.
   * Theorem 3.3 ("SL(n) implies `n = O(loglog R)`") is therefore a correct implication with a hypothesis that
     never holds.
6. **Repair (referee): Theorem R3 (proved) and Theorem R4 (proved, conditional).** Widen a light group `Q` of
   columns by factors `lambda_q`, with `sum_(q in Q) lambda_q log q <= W_P/8`. Then actual clusters admit no
   `UM_p` vector in `V_0(A_(Q,lambda))` for any `p in Q` (Theorem R3). Hypothesis **SL\*(n)** says every admissible
   profile has some light `(Q, lambda)` and some `p in Q` with such a vector. Theorem R4: SL\*(`n'`) for all
   `n' >= n_1` implies `n = O(loglog R)` on squarefree norms.
   * The referee left SL\* open, noting that a profile whose columns are short rational combinations of other
     columns in many disjoint ways would defeat it.
   * **README §5.1 shows that the Sylvester-Hadamard profiles are such profiles, so SL\* is false too**, and
     Theorem R4 is vacuous.
7. **Props 3.4, 3.5 and 4.1 (proved).**
   * Prop 3.4(ii)-(iii) characterise SL-full (the full space `V(S,m)`), not SL.
   * Prop 3.5 (obstructions) is the tool used in README §5.1. If `a'_p` lies in the rational span of a cheap
     family for `(S, m)`, no `c in V(S,m)` has `UM_p`.
   * Prop 4.1: the actual point `w` has `pi`-adic valuation vector `+-a'_p` at every split place.
8. **Evidence.** All the evidence models have at least `2m + 1` levels (the Boolean and Fourier sets, and the
   `had_sl` box), so they do not model `A_p`. On the faithful box the exact behaviour is that of Lemma R2.

## 4. Referee-corrected statements: subcluster (combined systems)

Write-up: `subcluster.md`. Report: `referee_subcluster.md`. Scripts: `subcluster_checks/` and
`referee_subcluster_checks/`. **Status: partial. Most results survive. The conditional theorem needs a restricted
hypothesis, and two interpretive claims and one computational claim are withdrawn.** The setting is normalised
clusters (`N` odd, all primes split, `W = log N`, `C' = max(C,1)`).

1. **Lemmas 1.0-1.2 (proved).**
   * `delta` lies in the column space `L = span(1, columns of (a_x, kappa_x))` and is nondegenerate relative to
     `L`.
   * Every integer `c` on `Z^M_0` whose order *relative to `L`* passes the LP test (1.2) gives an exact relation
     `sum c_k e^(i<k,delta>) = 0`. This is the *universal family*. Sub-cluster slab-box systems are instances of it,
     and they are genuinely new relations.
2. **Props 2.1-2.2 (proved).** The combined sub-cluster locus has dimension `<= n - 2` for every `M`. Each fibre
   over `n - 1` coordinates is either the whole circle or has `<= 2T` points. This says nothing about one actual
   point.
3. **Theorem 3.1 and Prop 3.2 (proved).**
   * On `L` the relations become trigonometric polynomials whose coefficients are the push-forwards `pi_* c`, and
     `ord_L(c) = ord(pi_* c)`.
   * For an affinely independent `J` the push-forward is injective. So Theorem 3.1 isolates only what
     *degenerate* sub-clusters add. In the full-rank case, the intersection over nondegenerate `J` is the open
     question itself; "nondegenerate sub-clusters add nothing" is withdrawn.
   * Product push-forward certificates work only in an `O(log)` window around the character threshold. The
     reading "no stronger than characters" is a leading-order statement; inside that window products can still
     matter for the constant.
4. **Theorem 4.1 and Corollary 4.2 (proved, unconditional).** Assume `W >= 16 ln M + 8 log^+ C + 14` and
   `W > 8 ln M + 8 ln(2W) + 12 ln C'`.
   * If `rho_hat >= 300 ln 4r`, the cluster's column torus is visible to the combined locus.
   * If `M >= (rho_hat+1)(rho_hat+2)`, an admissible *linear form in the points*, `sum_x c_x z_x = 0` with
     `ord_L >= 3`, sees the cluster.
   * Hence a cluster blind to the universal family has
     `M < (rho_hat+1)(rho_hat+2) < (300 ln(4W/ln 5) + 2)^2`.
   * Nothing shows that actual clusters are blind.
5. **Theorem 4.4 (proved, new; exact tests agree).** `reg(S) <= E_eps h_S(<eps,.>)` for every finite
   `S ⊂ Z^rho`, with equality on boxes and `l^1`-balls. It strengthens `threshold/construct.md` Theorem B and
   depends on coordinates.
   * Corollary 4.5: cross-polytope column geometries with weight ratio `gamma > 1/2` are blind to all slab-box
     sub-cluster systems.
   * The known examples are too small to have admissible sub-clusters. In the example with primes in
     `[P, P^1.9]`, pair separation can fail.
6. **ZE-tor must be local (referee, major).** ZE-tor asks that `Z_F`, inside the lifted ball of radius `1/(4D)`,
   is a finite union of rational subspaces containing `1`. At the fixed torus radius `1/(4D)`, subtori of large
   height have parallel sheets, and intersections have torsion cosets. With the local definition, Lemma 5.3(i)-(ii)
   and Prop 5.1 hold.
   * Prop 5.1: under local ZE-tor for every single sub-cluster system, every actual cluster in range is blind to
     the combined locus.
   * The rank target "`rho_hat <= theta ln W`" is then equivalent to the open combinatorial statement B(theta),
     with a sufficient pure-linear-algebra form B'(theta).
7. **ZE-col is false for thin families (referee Prop R, exact).** A single integer measure of order 2 on `Z^3` has
   a curve of zeros with nondegenerate points arbitrarily close to `0`. The proof uses the transversal common
   zero of its quadratic and cubic leading forms at a direction `(1, x, y)` with `[Q(x):Q] = 4`.
   * So Theorem 5.4 as stated, which assumes ZE-col(univ) for all profiles and all points, rests on a hypothesis
     that fails for thin families.
   * **Corrected Theorem 5.4 (proved, conditional).** Assume **ZE-col\***: for every actual cluster whose profile is
     robustly visible, the actual `delta` is not in `Z_univ`. Robustly visible means (V1) `rho_hat >= 300 ln 4r`,
     or (V2) `M >= (rho_hat+1)(rho_hat+2)` with an admissible basis of the degree-2 kernel.
   * Then `M <= max(K(C), (300 ln(4W/ln 5) + 2)^2) = O_C((loglog R)^2)` **on all norms**.
8. **Computations.**
   * Withdrawn: "below `n ~ 2137` the honest system is empty". At `r = n` the [D] slab-box is LP-supercritical for
     every `n in [440, 2400)` and Lemma-1.1-admissible for every `n in [468, 2400)`. Below `300 ln 4r` the
     combined locus is trivial by definition.
   * The §6.5 tangent cones are exactly `{0}` in all 560 systems, but only because `m <= 2T`, so this says nothing
     about ZE.
   * "Fakes certified off the loci" is automatic, since the loci are null sets.
   * Theorem 4.4 passed 950 random and 50 adversarial exact tests.

*README remark.* Remark 5.5 of the write-up says more concretely what (V2) asks. With `W >> E ln(MB)`, every
low-rank cluster satisfies `z_x = Q(a_x - a_b)` for one complex quadratic `Q` on the exponent lattice. Its
product-type consequences are only character inequalities. For example, a four-term progression in
`{a_x - a_b}` would give a pair with `X_v <= W/3`, so it is excluded by Prop 3.2. What (V2) asks beyond
characters is the open part of ZE-col\*.

## 5. Judge's findings (README remarks: unrefereed, with exact checks)

Four findings. J1 is a short proof from refereed ingredients. J2 and J3 are elementary but change the reading of
the round. J4 is evidence.

### 5.1 J1: Hypothesis SL\* is false (Sylvester-Hadamard profiles)

**Proposition J1.** Let `n = 2^k` with `k >= 13`. Build the Sylvester-Hadamard profile:
* points `x in F_2^k`;
* columns `a_p = 1[a.x = 0]` for `a in F_2^k \ {0}`, so `r = n - 1` and every `e_p = 1`;
* primes with all `log p` in a window `[tau, tau(1 + 2/n)]`.

Then the profile is admissible in the sense of SL\*: exponent rank `n - 1`, `n - 1 >= 300 ln 4r`, and pair
separation `>= 1/2 - 1/(2n)`. Moreover, for every base point, every light `(Q, lambda)` and every `p in Q`, **no
`c in V(A_(Q,lambda), m)` has `UM_p`**. Hence SL\*(n) fails at every such `n`, and Theorem R4 of
`referee_theory.md` is vacuous, like Theorem 3.3 of `theory.md`.

*Proof.*
1. **Column identity.** Fix distinct nonzero `a, b` and one of the four classes `c` of `(a.x, b.x)`. Comparing
   the four classes gives integers `alpha, beta, gamma` and `delta = +-2` with
   `1[a.x=0] = alpha + beta 1[b.x=0] + gamma 1[(a+b).x=0] + delta 1[(a.x,b.x) = c]`. For example, for
   `c = (1,1)`: `1[a.x=0] = 1 - 1[b.x=0] + 1[(a+b).x=0] - 2 * 1[a.x = 1, b.x = 1]`.
   * Exact check: `judge_checks/sl_star_hadamard.py`, all `(a, b, c)` for `k = 3..7`, 0 failures. The identity
     only involves the four classes, so it holds for every `k`.
   * Choose `c` different from the class of the base point `x_0`. In base-point coordinates this gives
     `a'_p = beta a'_b + gamma a'_(a+b) + delta sum_(x in V_c) e_x`, with `|V_c| = n/4` and `x_0` not in `V_c`.
2. **Choosing the pair.** Lightness gives `|Q| tau <= W_P/8 <= (n-1) tau (1 + 2/n)/8`, so
   `|Q| < (n-2)/2 - 1`. The `n - 2` nonzero vectors other than `a` form `(n-2)/2` pairs `{b, a+b}`. Pick a pair
   disjoint from `Q`; both of its columns keep their narrow slabs in `A_(Q,lambda)`.
3. **The family is cheap.** On `S = A_(Q,lambda)`:
   * `h_S(e_x) <= 2T = 2n`, because `S` lies in `B_T`;
   * `h_S(a'_b), h_S(a'_(a+b)) <= 2 floor(w/2)`, because `|<k,beta_q>| <= w` and `e_q = 1`.

   So the cost is at most `4 floor(w/2) + (n/4) 2n`.
   * By the proved bound `m >= 0.7 n(n-1)` ([D] Thm 2.1(b)) this is `< m` at `n = 2^14`:
     `1.737e8 < 1.879e8`.
   * The ratio `(n^2/2 + 2w)/(0.7 n(n-1))` decreases in `n`, so the bound holds for every `k >= 14`.
   * At `n = 2^13` it needs the exact `m = 49388999 > 4.708e7` (60-digit log-gamma, equal to the referee's
     big-integer value).
   * Table: `judge_checks/sl_star_hadamard_out.txt`.
4. **Conclusion.** `theory.md` Prop 3.5(ii) (refereed as correct) now says that no `c in V(S,m)` has `UM_p`.
   That proposition rests on the `p`-adic point `p^(a'_p)` lying on the cheap subtorus, hence in the base locus.
   Since `V_0(S)` is contained in `V(S,m)`, SL\* fails.
5. **Admissibility.**
   * The functions `1[a.x=0]` (`a != 0`) together with `1` span all characters of `F_2^k`, so the rows are affinely
     independent.
   * Each pair of points is separated by exactly `n/2` of the `n - 1` columns. With the prime window this gives a
     separation `>= (n/(2(n-1)))/(1 + 2/n) >= 1/2 - 1/(2n)`.
   * `n - 1 >= 300 ln(4(n-1))` for `n >= 2^12`. QED

**Reading.**
* This is the `p`-adic counterpart of the referee's refutation of ZE-braid. The valuation vector `a'_p` of the
  actual point lies in the span of a cheap family, so a single-place valuation certificate cannot see the cluster.
* Theorem 3.1 / R3 therefore give no information on Sylvester-Hadamard profiles (`n >= 2^13`), at any widening.
* A residue-level variant (the face residue sums, `theory.md` SL-res) is not touched by this argument. Round-10
  Prop 5.5(iii) shows, however, that residue collisions push towards vanishing.
* The referee's open question ("do admissible profiles with columns that are short combinations of other columns
  in many disjoint ways exist?") has a positive answer, so SL\* is false.

### 5.2 J2: at T = n the leading-order (tangent-cone) criterion is vacuous for near-extremal clusters

**Proposition J2.** Let `c in V(S,m)` and `|u|_inf = 1`, and suppose `|<k',u>| <= R` on `supp c`. Then
`|Phi_c(u)| <= 2 |c|_1 (R/2)^m`. Hence `kappa(u) <= 2^(1-m)` in the normalisation of `theory.md` (2.1), where
`R <= D/2`. For an actual cluster `|u|_inf >= e^(-W/2)`, because distinct lattice points are at distance `>= 1`.
So criterion (2.1) can hold only if

```text
W  >  W_TC := 2(m-1) ln 2 + 2 ln((D/2)/(m+1)).
```

*Proof.* The push-forward `mu` of `c` under `k -> <k',u>` is a measure on `[-R, R]` whose moments of order `< m`
vanish. So `int x^m dmu = int (x^m - p(x)) dmu` for every `p` of degree `< m`. Take `x^m - p` to be the monic
Chebyshev polynomial `2^(1-m) R^m T_m(x/R)`, whose sup-norm on `[-R, R]` is `2(R/2)^m`. QED

**Exact data** (`judge_checks/scale_T.py`, output `scale_T_Tn_out.txt`). At `T = n`, `W_TC ~ 1.02 n^2`:
* `W_TC = 4.66e6` at `n = 2137`, against `W_min = 3n ln n - 14n = 19236`, a factor 242;
* the factor is 306 at `n = 3000`.

So with round 10's `T = n`, the sufficient condition (2.1) of `theory.md` Lemma 2.4 **never** holds for a
cluster with `W < ~n^2`. That includes every near-extremal cluster (`W = O(n log n)`).
* The write-up's "(2.1) is a very weak requirement" is wrong.
* The obstruction of Prop 2.5 is not confined to neighbourhoods of the strata; it is everywhere.
* Only the sufficient criterion with the uniform remainder fails. A sharper comparison with the `(m+1)`-st moment
  of the same `c` is not excluded.
* *Heuristic reading:* at `T = n` the actual position (`|delta| ~ n^(-3n/4)`) is far outside the scale (`~2^(-m)`,
  `m ~ 0.74 n^2`) on which order-`m` Taylor data control `Z_P`. So ZE at `T = n` is a semi-global statement, not a
  question about the germ at `e`.

### 5.3 J3: route 5 works at box size T = 1, where the scales separate

**Proposition J3 (constants checked exactly in the range shown; a referee should recheck them).** Run [D]
Theorem 2.1 and Theorem 4.1 with `T = 1` in place of `T = n`. So `B_1 = {k in Z^n_0 : |k_x| <= 1 (x < n)}`,
`w = sqrt(8(n-1) ln 4r)`, `D = 2(n-1)`, `m` is the largest integer with `8 binom(m+n-2, n-1) <= 3^(n-1)`, and
`H = U^2`.
* (a) `|A| >= 3^(n-1)/2`, by the same Hoeffding argument, which does not depend on `T`.
* (b) `m/n` lies between 0.6421 and 0.6440 for `n` in `[1000, 8000]` (exact big-integer binomials). The entropy
  bound `binom(N,K) <= e^(N h(K/N))` gives `m >= c_0 n - O(log n)` with `c_0 ~ 0.6418`, from
  `(1+c) ln(1+c) - c ln c = ln 3`.
* (c) `dim V_H >= |A|/2`, by the same Siegel count with `B = T^(m-1) = 1`. It needs `ln U > 3.3`, which holds
  for `n >= 5`.
* (d) `m >= 4w` for every `n` in `[2900, 9000]` (`r = n`), `[3623, 9000]` (`r = n ln n`) and `[5809, 9000]`
  (`r = n^2`), exactly. Compare round 10's 2137, 2679 and 4278 at `T = n`.
* **Theorem 4.1 at `T = 1`.** The comparison `(m/2 - w)W <= 2 ln|c|_1 + 2m ln(DC/2) - 2 ln m!`, with
  `|c|_1 <= U^3` and `U <= 3^(n-1)`, gives a contradiction once `W > W_1(n, C)`, where
  `W_1 = 51.75 + O(log^+ C)` at `n = 3000` (`C = 1`, maximal `m`), and this `W_1` decreases in `n`. So every actual cluster with
  `m >= 4w` and `W > W_1` satisfies (4.2) for the `T = 1` slab-box. The conditional Theorem 5.7 carries over
  with these thresholds.
* **Scale separation.** At `T = 1`, `W_TC ~ 1.39 m ~ 0.89 n`, which is far below `W_min` (2672 against 30057 at
  `n = 3000`). At `W = W_min` the Chebyshev ceiling `2^(1-m)` exceeds the requirement of (2.1) at
  `|u| <= e^(-W/4)` by a factor `e^6178` (`n = 3000`), `e^9101` (`n = 4000`) and `e^22363` (`n = 8000`).

So at `T = 1`, nondegenerate actual positions (`|delta| <= C e^(-W/4) = e^(-Theta(n log n))`) lie deep inside the
region where the order-`m` Taylor data decide (`m ~ 0.64 n`). Away from the cheap strata, ZE then becomes a
statement about the **real tangent cone** `TC_H(A) = {u : Phi_c(u) = 0 for all c in V_H}`, together with a
lower bound for `sup_c |Phi_c(u)|/(|c|_1 (D/2)^m)` of size `>= e^(-W/4) poly(n, C)`.
* That bound is generous: the Chebyshev ceiling of the quantity is `2^(1-m) = e^(-0.45 n)`, against
  `e^(-W/4) <= n^(-3n/4 + O(1))`. At the smaller order `m = ceil(4w)` the ceiling is `e^(-O(sqrt(n log n)))`;
  see "Choice of the order" below.
* The compute.md models `M8` and `M8'` are `T = 1` models.

**Near the strata (heuristic, partly proved).**
* *Proved.* Let `F` be a cheap family of cost `d_F`, and `L_F` its span. Then every `Phi_c` vanishes to order
  `>= m - d_F` transversally along `L_F`.
  * Indeed `<k,u_L>^(m-i)` agrees on `A` with a polynomial in the cheap coordinates of degree
    `<= min(m - i, d_F)`, so the coefficient of `eta^i` in `Phi_c(u_L + eta v)` vanishes for `i < m - d_F`.
* *Heuristic.* At relative distance `eta` from `L_F`, (2.1) therefore needs roughly
  `eta^(m - d_F) >~ e^(-W/4) 2^(O(m))`.
  * At `T = 1` and `d_F <= m/2` this is `eta >~ n^(-2.3)`, a **polynomial** separation, since
    `W/(2m) ~ 2.34 ln n` at `W = 3 n ln n`.
  * At `T = n` it is `eta ~ 1`, which is impossible.
* *Arithmetic supply.* The pair inequality gives `eta >= e^(-(d_xy - W/2)/2)/(C sqrt n)` for any pair `(x, y)`
  orthogonal to `L_F`: same pattern on the cheap columns, neither point in the coordinate set.
  * For near-extremal clusters the average of `d_xy - W/2` over pairs is `<= W/(2(n-1)) ~ 1.5 ln n`.
  * So at least half of all pairs give a polynomial `eta`.
  * Whether one of them is orthogonal to a given cheap stratum is a combinatorial question about the profile.
* So at `T = 1` the obstruction of `theory.md` Prop 2.5 becomes polynomial, where at `T = n` it was exponential.
  Whether the quantitative (Lojasiewicz-type) lower bound on the twisted leading forms holds with the needed
  constants is open; this is the main risk of the programme in §8.

**Choice of the order: take `m` as small as the LP allows, `m' = ceil(4w)`.**
[D] Thm 4.1 needs only `m >= 4w`. Any `m'` in `[4w, m]` works:
* `V(A, m')` contains `V(A, m)`;
* the Siegel count needs only `binom(m'+n-2, n-1) <= |A|/4`;
* the comparison gives `W_1(m') = 53` at `n = 3000` and `631` at `n = 10^6`, against `W_min = 3.0e4` and `2.7e7`.

The data are in `judge_checks/scale_m4w_out.txt` (`r = n`, `C = 1`; floating point except the log-gamma `m`).
Three things improve.
1. **Few long directions are cheap together.** *Sketch (unrefereed).*
   * Call an integer direction `v` long if `h_A(v) >= w - O(1)`; a balanced `e_p = 1` column attains its slab and
     is long.
   * Since `A` contains every `sum_P e_x - sum_N e_y` with `|P| = |N| <= w/2`, a direction with fewer than `~w/2`
     non-modal entries has width at least about that number.
   * So a cheap family (total width `< ceil(4w)`) has at most four long directions, plus short directions of total
     support `O(w)`.
   * If the long directions are 0/1-valued off `O(w)` points, as columns are, a cheap span separates the points
     into at most `2^4` classes plus `O(w)` exceptional points. So for `n >> w` it lies in the braid arrangement.
   * Long directions taking many values need a separate argument; this is part of the sketch.
   * The referee's ZE-braid counterexample needs about `log_2 n` cheap columns, so it does not apply at this order.
2. **The near-strata radius shrinks to `e^(-W/(4m')) = e^(-Omega(sqrt(n log n)))`.** This beats the pair
   separation `eta ~ e^(-W/(4(n-1)))/sqrt n` of a Plotkin-tight pair for every `n >= 10^4` in the table (at
   `n = 10^5`: `ln rho_0 = -40.0` against `ln eta = -10.9`). It does not do so at `n = 3000` or `5000`.
3. **The slack of the quantitative tangent-cone bound becomes huge.** The Chebyshev ceiling of
   `sup_c kappa_c` is `2^(1-m') = e^(-O(sqrt(n log n)))`, against the required `e^(-W/4) = n^(-3n/4 + O(1))`. A
   lower bound that loses up to `e^((3/4) n ln n - O(sqrt(n log n)))` against the ceiling would still suffice.

### 5.4 J4: the tangent-cone form of ZE holds in the two T = 1 models

**Corollary J4 (follows from the refereed data of Theorem 5.2).** For `M8` and `M8'`, every real or complex common
zero `u` of the forms `Phi_c`, `c in V(A,m)`, lies on one of finitely many rational hyperplanes through 0:
* the collision hyperplanes `u_x = 0` and `u_x = u_y`;
* in `M8`, also the tangent hyperplanes `<a - a', u> = 0` of the two-term factors.

*Proof.* Take a basis polynomial of "form type",
`P_c(w) = const * w^a * prod (w^b - w^(b'))^(e) * prod_f H_f(w - 1)^(e_f)`, with each `H_f` a nonzero form.
Substitute `w = e^z`. The lowest-order part of `P_c(e^z)` is `const * prod <b - b', z>^e * prod H_f(z)^(e_f)`.
On the other hand, the degree-`m` part of `P_c(e^z) = sum_k c_k e^(<k,z>)` is `Phi_c(z)/m!`. So, for these `c`,
`Phi_c` equals `m! const` times the product of the linear tangent forms of the two-term factors and the `G_j`.
If `u` avoids all those hyperplanes, then `G_j(u) = 0` for every `j in J`. The saturation of Theorem 5.2
(the ideal `(G_j : j in J)` saturated by the collision forms is the unit ideal) then puts `u` on a collision
hyperplane. QED

So the tangent-cone statement of §8 holds in both `T = 1` models. This is evidence of the same weak kind as
Theorem 5.2: the models are LP-subcritical. An independent check, a direct Groebner saturation of the 42 octics of
`M8` (`judge_checks/tc_m8.py`), was abandoned after about 25 minutes: the initial Groebner basis did not finish
with flint's naive Buchberger. The corollary does not depend on it.

## 6. Status of ZE

**ZE is neither proved nor refuted in the regime.** The evidence for it is weak and structurally limited. There is
no evidence against it in the regime. Every *stronger* or *delta-free* reformulation tried this round is false.

| statement | status after round 11 | source |
| --- | --- | --- |
| **ZE(n)** ([D] §5.5): no nondegenerate `delta`, `max|delta_x| <= 1/(4D)`, in `Z_P` | **open** | — |
| ZE-alg ([D] open point 1): only algebraic `w in Q(i)`, independent, height `<= W`, within `e^(-W/4)` of `e` | **open** (weakest form that suffices) | — |
| ZE-tor, local (lifted ball; `Z_P` = finite union of rational subspaces containing `1`) | **open**; implies ZE | referee_subcluster §7.1 |
| ZE-braid: realisable symbols vanish only on the braid arrangement | **false** in the regime (Sylvester `n = 2^k >= 2^16`) | referee_compute §2.8 |
| SL(n): delta-free, one widened slab | **false** (twins and duplicates; odd `m`) | referee_theory Prop R1, Lemma R2 |
| SL\*(n): delta-free, grouped widening | **false** (Sylvester `n = 2^k >= 2^13`) | README §5.1 (J1, unrefereed) |
| ZE-col(univ) for all profiles and all points of `L` | **false** for thin families | referee_subcluster Prop R |
| ZE-col\*(univ): actual `delta`, robustly visible profiles | **open**; implies `M = O((loglog R)^2)` on all norms | referee_subcluster §7.2 |
| ZE-real when `V` is small (`dim V = 1`, zonotope-free) | **fails** (exact) | compute Prop 6.2 |
| ZE on `M8`, `M8'` (Sylvester `H_8`, `T = 1`, zonotope-free, LP-subcritical) | **holds** on the whole complex torus (exact) | compute Thm 5.2 |
| tangent-cone form on `M8`, `M8'` | **holds** | README §5.4 (J4) |
| ZE at the README §8.4 sizes `n - 1 = 2, 3, 4` | holds, but only through zonotope (character) certificates, which never exist in the regime | compute Props 2.1-2.4 |

**What the round established about ZE.**
1. **It is genuinely an overdetermination statement.**
   * Single forms always have non-toric zero sets through `e` (compute §5.2; the real version is evidence).
   * Systems with `dim V = 1` (compute Prop 6.2), and thin families (Prop R), have nondegenerate zeros arbitrarily
     close to `e`.
   * So any proof must use many sections at once. The regime does have this: `dim V_H >= |A|/2`, and by
     theory.md Lemma 1.1 `V_H` has relative codimension `O(1/n)` in `V(A,m)` at height `e^(mW/40)`.
2. **The toric part of the base locus is large and structured.**
   * `Bs(A,m)` contains the span of every cheap family. These combine up to `~0.37n` coordinates (`T = n`) with up
     to `~m/w ~ 0.25 sqrt(n/ln 4r)` columns.
   * For Sylvester-Hadamard profiles there are many such families, some with collision-free generic points
     (the ZE-braid counterexample).
   * Their `p`-adic shadows contain the valuation vectors `a'_p` of actual points (J1).
   * Any correct reformulation must allow all rational subtori; only the archimedean position of the actual
     `delta` (its nondegeneracy) distinguishes it from them.
3. **Small models cannot test the regime.**
   * Zonotope-free (faithful) models first exist at `n = 8`.
   * The regime's LP margin `h_A <= m/2`, together with multi-column strata, needs `m >= 4w`. That is
     `n >= 2137` at `T = n` and `n >= 2900` at `T = 1` (`r = n`).
   * For `H_8` even the weaker condition (ii) is out of reach at `T = 1`, and at `T = 2` with `L <= 2`.
4. **At round 10's box size `T = n`, ZE is not a tangent-cone question** (J2): the leading-order criterion can
   never certify a near-extremal cluster. **At `T = 1` it is** (J3): the actual position lies deep inside the
   Taylor regime, and Prop 2.5's near-strata obstruction drops from exponential to polynomial size (heuristic).

**Evidence for ZE:**
* compute Theorem 5.2, with the referee's caveats: two zonotope-free models where only the intersection of many
  sections is toric;
* J4: the tangent-cone version in the same models;
* round 10's Prop R1: no cuts implies ZE.

**Evidence against ZE in the regime:** none. ZE fails only where it was expected to fail: `dim V = 1` and thin
families. Prop 6.2 also shows that, for an underdetermined model, a nondegenerate position *on* `Z_P` exists for
the fake profile (skew double core, `q = 7`). That is the barrier-transfer scenario of `round10/README.md` §8.5,
realised only outside the regime.

## 7. What this means for the determinant route and the growth rate

**Growth rate: unchanged.** `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for
actual circles. Nothing this round improves it unconditionally or from a standard conjecture. The uniform bound
is not touched.

**Conditional theorems on record after round 11.**

| hypothesis | conclusion | status of the implication | status of the hypothesis |
| --- | --- | --- | --- |
| ZE(`n'`) for all `n' >= n_1` | `n = O(loglog R)`, squarefree norms | proved ([D] Thm 5.7) | open |
| same, with the `T = 1` slab-box | same | README §5.3 (J3, constants exact in range; unrefereed) | open; the tangent-cone reduction of §8 applies here |
| ZE-col\*(univ) (actual `delta`, robustly visible) | `M = O_C((loglog R)^2)`, **all norms** | proved (subcluster Thm 5.4 as corrected) | open |
| local ZE-tor for single sub-cluster systems | actual clusters are blind to `Z_comb`; rank target `<=>` B(theta) | proved (Prop 5.1 with local ZE-tor) | open |
| SL(n) | `n = O(loglog R)` | proved (theory Thm 3.3) | **false** (vacuous) |
| SL\*(n) | `n = O(loglog R)` | proved (referee Thm R4) | **false** (J1; vacuous) |

**Unconditional gains** (all proved and refereed):
* Lemmas 1.1-1.2 of theory.md: `V_H` is almost all of `V(A,m)`, and the dichotomy holds on any support of cost
  `<= 3mW/8`.
* Theorem 3.1 / R3: a delta-free necessary condition on profiles. It is empty for odd `m` and for
  Sylvester-Hadamard profiles (J1).
* Prop 2.1 and Prop 2.2 of compute.md: zonotope certificates, and their absence in the regime.
* Theorem 4.4 of subcluster.md (signed mean-width bound).
* Corollary 4.2 of subcluster.md: a cluster blind to the universal family has
  `M < (300 ln(4W/ln 5) + 2)^2`.
* Lemma R2 (levels on widened boxes).

None of these bounds an actual cluster.

**Meaning for route 5.**
1. **The arithmetic is spent, and only the position is left.**
   * Thm 4.1, Lemma 1.2 and Thm R3 extract everything that integrality, forced divisibility and single-place
     valuations give.
   * Every attempt to finish without the archimedean position of `delta` has failed: SL and SL\* (profile only),
     and ZE-braid (algebraic, braid only). Each failed through the same object, the cheap subtori.
   * The only input that separates the actual point from those subtori is its nondegeneracy (Lemma 3.2 of [D],
     Lemma 1.0 of subcluster.md), together with its quantitative form (the pair inequalities).
   * So the missing step must be ZE-like: a statement about `Bs(A,m)` near nondegenerate positions, possibly
     restricted to arithmetic points (ZE-alg).
2. **The round-10 test plan is retired.**
   * The finite test of `round10/README.md` §8.4 is vacuous at its sizes (zonotopes), and only toy-faithful from
     `n = 8`.
   * The Seshadri-type degree bounds of §8.5 are only necessary conditions (compute §8 item 5).
   * The tangent-cone test (§8.4 step 2) is meaningless at `T = n` (J2), but meaningful at `T = 1` (J3).
3. **Non-squarefree norms.** subcluster.md gives a second conditional route: ZE-col\*, `(loglog R)^2`, no rank lemma
   needed. The low-rank part of ZE-col\* (V2) is the statement "`e^(i<y,s>)` is quadratic on the exponent set
   only for degenerate `s`" (Remark 5.5). It is open, and finite to state.
4. **Barrier logic is unchanged.**
   * The input that fakes violate is still "almost every fake" (Cor 4.3 of [D]). "No fake lies on `Z_P`" is part
     of ZE.
   * Prop 6.2 shows that, for an underdetermined model at `n = 8`, nondegenerate positions on `Z_P` exist for a
     fake profile. Whether a gap-admissible fake can be steered onto `Z_P` in the regime remains open. If one can,
     ZE fails and route 5 is capped as formulated.
5. **Route 5 is still the best-scoring open route, but its residual hypothesis is no closer to proof.**
   * The round sharpened what ZE must be: overdetermined, allowing rational subtori, using the archimedean
     position.
   * It also found the box size at which ZE becomes a tangent-cone question (J3).

## 8. The single most promising next step

**Recommended: run route 5 at box size `T = 1` and the smallest admissible order `m = ceil(4w)`. There ZE
reduces to the real tangent cone. Then decide the tangent-cone hypothesis ZE-TC for slab-cut ternary cubes.**

### 8.1 Why this step

* **It avoids every failure mode found this round.**
  * Delta-free versions (SL, SL\*) fail on Hadamard and twin profiles.
  * Braid-only versions fail on the multi-column subtori.
  * Thin and underdetermined systems have nondegenerate zeros.
  * ZE-TC keeps every rational subspace (the cheap spans), keeps the full overdetermined system `V_H`, and uses
    the archimedean position. Those are the three features the round showed to be necessary.
* **At `T = 1` the scales separate** (J3).
  * The actual position has `|delta| <= C e^(-W/4) = n^(-3n/4 + O(1))`.
  * The order-`m` Taylor data live at scale `2^(-m)`: `e^(-0.45 n)` at the maximal `m`, and
    `e^(-O(sqrt(n log n)))` at `m = ceil(4w)`.
  * So ZE becomes a question about the real tangent cone, plus a near-strata lemma whose required separation is
    polynomial at the maximal `m` and `e^(-Omega(sqrt(n log n)))` at `m = ceil(4w)`.
  * At `m = ceil(4w)` at most four long directions (for example balanced columns) are cheap together (sketch,
    §5.3). So the Hadamard multi-column strata that refute ZE-braid at `T = n` are not present.
  * At `T = n` neither is true (J2), and ZE is semi-global.
* **It is a real-algebraic statement in which no exponentials appear**, and it is decidable for each `(A, m)`.
  * On the uncut ternary cube it holds exactly: there the leading-form ideal is `(k_x^3 : x < n)`, and the tangent
    cone is precisely the union of the sparse coordinate spans (`|supp u| < m/2`).
  * In the two zonotope-free `T = 1` models it holds (J4).

### 8.2 The statements

Fix an admissible profile and its `T = 1` slab-box `(A, m, V_H)` (README §5.3), with `D = 2(n-1)`. Take
`m = ceil(4w) = Theta(sqrt(n log r))`; any `m` in `[4w, m_max]` is admissible, but the small one is better. For
`|u|_inf = 1` put `kappa_c(u) = |Phi_c(u)|/(|c|_1 (D/2)^m)`.

* **Step 0 (proved; theory.md Lemma 2.4, read in reverse).** If an actual cluster has `delta in Z_P`, then
  `u_hat = u/|u|_inf` satisfies `kappa_c(u_hat) <= eps` for every `c in V_H`, with
  `eps = (D/2)|u|_inf/(m+1) <= C n e^(-W/4)`.
* **ZE-TC(A, m), qualitative (open).** Every real `u != 0` with `Phi_c(u) = 0` for all `c in V_H` lies in the real
  span of a cheap family. Equivalently, by compute.md R2, `<k,u>^m` lies in the degree-`m` part of the
  leading-form ideal of `I(A)` only for such `u`.
* **ZE-TC_eps (open; the form needed).** Every `u_hat` with `kappa_c(u_hat) <= eps` for all `c` lies within
  relative distance `rho_F` of the span `L_F` of some cheap family `F` of cost `d_F`. Here `rho_F` is of order
  `(eps 2^m)^(1/(m - d_F))` up to polynomial factors. This is the size forced by the transversal vanishing order
  `m - d_F`, which is proved in §5.3.
* **Arithmetic separation (open, combinatorial; heuristically easy for large `n`).** The actual `u_hat` has
  relative distance `> rho_F` from every cheap span.
  * The pair inequality supplies this whenever some pair `(x, y)` orthogonal to `L_F` (same pattern on the cheap
    columns, both outside the short supports) has `d_xy - W/2 < W/(2(m - d_F)) - ln n - 2 ln C`.
  * At `m = ceil(4w)` the right side is `Omega(sqrt(n log n))`, while the average of `d_xy - W/2` over pairs is
    `<= W/(2(n-1)) = O(log n)`. So all but a vanishing fraction of pairs qualify.
  * At most `2^4` column classes occur (sketch, §5.3), so orthogonal qualifying pairs exist once `n >> w`.
  * The Hadamard `k`-column strata, which separate all points, are not cheap at this order.

Step 0, ZE-TC_eps and arithmetic separation together contradict [D] Thm 4.1 at `T = 1`. With [D] Lemmas 3.1-3.2
this gives `n = O(log W) = O(loglog R)` on squarefree norms. In the heuristic table of §5.3 the scale comparison
favours the arithmetic only from `n ~ 10^4`, which is harmless for an asymptotic statement. Turning the sketches
of §5.3 into proofs is part of the step.

### 8.3 First concrete tasks

1. **Referee J1-J3** (this README, §5).
   * J3 needs a closed-form proof that the maximal `m` at `T = 1` is at least `4w` for all `n >= n_0`. The entropy
     bound `binom(N,K) <= e^(N h(K/N))` gives `m >= 0.64 n - O(log n)`.
   * J3 also needs the `T = 1` thresholds and the constant `W_1`.
2. **Qualitative ZE-TC for the ternary cube with few cuts.**
   * The ideal is `(k_x^3) +` the leading forms of the slab polynomials.
   * Decide ZE-TC with one cut and with two cuts, for arbitrary sign columns. Then decide it for the
     Sylvester-Hadamard family, using the `AGL(k,2)` symmetry.
   * Look for an induction on the number of cuts. This is the tangent-cone analogue of round-10 Prop R1, which
     handles zero cuts.
   * The compute.md method (exact factorisation into two-term factors and forms in `x = w - 1`, then saturation by
     rational hyperplanes) is the computational tool for small instances.
3. **The near-strata lemma at `T = 1`.**
   * Prove the transversal-order statement of §5.3 for the twisted forms.
   * Count near-tight pairs orthogonal to each cheap span on admissible profiles. This needs a better count than
     Markov's inequality on the pair excesses.
4. **Fake test for the step.** The skew double core at `q = 7` (`T = 1`, `dim V = 1`) already has nondegenerate
   real zeros of its single tangent form (compute Prop 6.2). So ZE-TC, like ZE, can hold only for overdetermined
   systems. In the regime, a gap-admissible fake whose direction `u_hat` lies in `TC_eps` would refute ZE-TC_eps
   and cap route 5. Fakes satisfy the pair inequalities, so such a fake is automatically far from the cheap
   spans. Search for one first in the Paley and skew families, using the symmetry.

### 8.4 Reading of the outcomes

* **ZE-TC fails for a regime-type `T = 1` box** (a nondegenerate real direction killed by every `Phi_c`). Then look
  along that direction at the higher-order terms.
  * A nondegenerate real branch of `Z_P` near `e` refutes ZE-real.
  * Only ZE-alg (the arithmetic of `w`) is then left, and route 5 is probably capped.
* **ZE-TC holds with polynomial Lojasiewicz control.** Then route 5 reduces to the combinatorial separation lemma
  about pair excesses and cheap spans, a statement about profiles of actual clusters, not about positions.

**Not recommended:**
* further delta-free single-place valuation hypotheses (the SL family; J1);
* algebraic reformulations that exclude rational subtori (ZE-braid);
* finite tests at `T = n` or at `n - 1 <= 4` (zonotopal);
* the tangent-cone criterion at `T = n` (J2).

## 9. Files

All paths are relative to this folder (`round11/`). Each script runs as `python3 <script>` from its own directory,
and outputs sit next to the scripts.

| write-up | report | scripts (author / referee) |
| --- | --- | --- |
| `compute.md` (notes: `notes.md`) | `referee_compute.md` | `compute_checks/` / `referee_compute_checks/` (notes in `notes.md` there) |
| `theory.md` (notes: `theory_notes.md`) | `referee_theory.md` | `theory_checks/` / `referee_theory_checks/` |
| `subcluster.md` (notes: `subcluster_notes.md`) | `referee_subcluster.md` | `subcluster_checks/` / `referee_subcluster_checks/` |

**Judge's checks** (`judge_checks/`; working notes in `judge_notes.md`):

| script | what | status |
| --- | --- | --- |
| `sl_star_hadamard.py` | J1: the column identity on `F_2^k` (all `(a, b, class)`, `k = 3..7`); the cheapness inequality `4 floor(w/2) + n^2/2 < m` for `n = 2^k`, `k = 11..20`, with the proved bound `m >= 0.7n(n-1)` and with 60-digit log-gamma `m` (`k <= 16`); regime conditions | exact identity; integer and proved-bound inequality |
| `scale_T.py` | J2 (`T = n`: `W_TC` against `W_min`) and J3 (`T = 1`: exact big-integer `m`, thresholds for `m >= 4w` over `n in [1000, 9000]`, Theorem 4.1 constant `W_1`, room for (2.1)). Outputs: `scale_T_Tn_out.txt`, `scale_T_T1_out.txt` | exact `m`; elementary floating point elsewhere |
| `scale_m4w.py` | J3b: `T = 1` with `m' = ceil(4w)`, giving `W_1(m')`, the number of cheap columns, and the heuristic scale comparison `rho_0` against `eta` | log-gamma `m`; heuristic comparison |
| `tc_m8.py` | J4 cross-check: direct Groebner saturation of the 42 octic tangent-cone forms of `M8` by rational hyperplanes | **abandoned** (initial Groebner basis unfinished after ~25 min); not needed for Corollary J4 |

**Caveats.**
* J1-J4 and the programme of §8 are unrefereed. J1 rests on refereed ingredients (theory.md Prop 3.5(ii) and
  [D] Thm 2.1(b)) plus the exact checks above.
* J3's thresholds are exact only in the scanned range. All-`n` statements need the entropy bound, which is
  sketched but not written out.
* The near-strata analysis of §5.3 (transversal order proved; the rest heuristic) is a scale count, not a proof.
* The "at most four cheap columns, hence cheap spans lie in the braid arrangement" claim at `m = ceil(4w)` is a
  sketch.
* The repository is read-only. The corrections below should be applied when round 11 is merged:
  * `round10/README.md` §8.4-§8.5: the finite test is vacuous at `n - 1 <= 4`; the Seshadri route gives only
    necessary conditions; the tangent-cone test is meaningless at `T = n`;
  * `round10/determinant.md` §5.5: the tangent-cone first test should be run at `T = 1`;
  * `theory.md` Lemma 2.4: the sentence "(2.1) is a very weak requirement" is false at `T = n` (J2).
