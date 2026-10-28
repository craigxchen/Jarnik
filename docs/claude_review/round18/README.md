# Round 18 judge: the arithmetic Ptolemy lemma (L_8)

Three write-ups this round, each attacked by an adversarial referee whose report sits next to it:
* `census.md` (round-17 task 1: a census of small-height balanced anchored configurations), with `referee_census.md`;
* `lift.md` (round-17 task 2: do the `(L_6)`-refuting families `E_t` lift to `k = 7, 8`?), with `referee_lift.md`;
* `theory.md` (proof lane and refutation lane for `(L_8)`), with `referee_theory.md`.

§2 gives the referee-corrected statements. Remarks added by this README are marked *README remark*. Nobody but the
judge has checked them, so they are **unrefereed**. Labels used: proved / conditional / exact data / evidence /
heuristic / unrefereed. Paths are relative to this directory unless they start with `docs/` (repository, branch
`main`) or `branch:` (`origin/agent/uniform-gcd-formalization`).

## 0. Verdict

**No growth improvement and no uniform bound were proved this round**, unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`docs/claude_review/two_thirds/proof.md`) is still the best bound for
actual circles.

* **`(L_8)`: open.** It is neither proved nor refuted. No restricted version holds on a Zariski-dense set with uniform
  constants. **No concrete exceptional set `Z` was found.** All three referees confirm that the object studied is
  exactly `(L_8)` of `ptolemy.md` §7 as Steps 0-5 of §6.2 use it (§1).
* **What the round established** (all refereed; none of it moves growth):
  1. *Reduction to contact profiles* (census Prop 1.3). This is a sandwich, not an equivalence:
     `(L_8^c)^all => (L_8) => (L_8^f) => (L_8^c)^cat` (theory, as corrected). The chain hypothesis equals the
     caterpillar condition on the block system (theory Lemma 1.1, checked exactly on all 39208 trees with 8
     leaves). On pairwise-coprime profiles, which every proof must cover, chains, anchoring and residues carry no
     information (theory Prop 2.1).
  2. *Emptiness below an explicit scale* (census Thm E). A balanced 8-datum needs at least 24 primes. At
     `eta = 0.1` the class is empty for `e^w < 59.2`, and every minor `det(P_i,P_j)` with `i, j >= 1` is at least
     `10^49.4`. **No feasible search can test `(L_8)`.** Block-sampling lattice searches return `H_Sigma ≍ w` even at
     `k = 5, 6`, where bounded-`H` points are known to be abundant.
  3. *Loci on which `(L_8)` holds, with locus-dependent constants* (census Prop 7.1, 7.4, 7.5). These are:
     every subvariety missing some `D_T` (all symmetric loci, all forgetful fibres, the boundary), and each fixed
     curve of genus `<= 1` with `K.C <= 0`. Covering genus-`<= 1` families have class imbalance `>= 5/13`.
     **This does not exclude a dense supply with finitely many points per member** (referee M1: the `Q`-lines of
     Kapranov's `P^5` pass through every rational point).
  4. *The lifting test, curve level* (lift Thm 4.2, unconditional, exact over `Q` plus openness). For `t` in a
     nonempty Zariski-open `U_0`, every class-`beta_7` curve over `E_t` is a pencil lift of genus 9 (or 5). **The
     `(L_6)` mechanism does not lift through class `beta_7`.** Classes `c beta_7` with `c >= 2`, and all
     point-level lifts, are open.
  5. *A genuine gap in the barrier map* (theory Prop 2.6, conditional on abc). A unit pencil `lambda A/C` on
     `Mbar_(0,8)` with `gain = 10 - K^log_rel.beta > 0` would give `(L_8^pc)`. The refereed abc barrier does not
     cover such pencils. **No pencil with positive gain is known.** Morphisms have gain `<= -32`, and so do all
     tested pencils:
     * every unit pencil with `lambda = ±1` and `D <= 4` (277,940 pencils);
     * every rational `lambda` with `D <= 3` (referee).

     Even a positive pencil would give only `(L_8^pc)`, a strict consequence of `(L_8)`. That is not enough for
     Theorem D without a further, unproved extension (§3.3).
  6. `(L_k) => (L_(k+1))` (census Lemma 7.7). Expected dimensions of balanced curves are `-10d + 2` (genus 0) and
     `-10d` (genus 1), against `7d + 1` at `k = 7`.

**Scores** (rubric of `docs/claude_review/round10/README.md` §1, each out of 10; details in §6):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| census | 7 | 1 | 1 | **5** |
| theory | 7 | 1 | 2 | **5** |
| lift | 6 | 1 | 1 | **4** |

**Recommendation (§7). Stop internal rounds on growth and uniformity.** The round-17 target `(L_8)` has been
examined from three directions. The outcome:
* no method in either direction;
* no concrete exceptional set;
* a proof that no feasible computation can test it.

The consolidated expert document `docs/claude_review/overview/overview.md` already exists. The single most useful
next step is external: add a short round-18 addendum to it (§7.3) and put `(L_8)` before experts in Diophantine
approximation, in its lattice form `(L_8^pc)` together with the abc-gain question. The only internal item with a
non-negligible chance of a definite answer is the abc-gain search at degree 5 (§7.2). It cannot give unconditional
growth.

## 1. Object and normalisation (checked against ptolemy.md §7 and §6.2)

All three write-ups and all three referees use the same object. The judge re-checked it item by item.

**`(L_8)`** (`docs/claude_review/growth/ptolemy.md` §7). There are `c > 0`, `eta_0 > 0` and one proper Zariski-closed
`Z ⊂ Mbar_(0,8)`, fixed before the data, such that the following holds. Let `Phi = [P_0 : ... : P_7] ∈ M_(0,8)(Q) \ Z`
with:
* anchored primitive integer vectors, `P_0 = (1,0)`, `Y_j = det(P_0, P_j)`;
* `det(P_i, P_j) = t_ij prod_(T ∋ i,j) n_T`;
* at each prime, the blocks `n_T` with `p | n_T` form two nested chains on disjoint row sets (item 532);
* `|log n_T - w| <= eta_0 w`.

Then `H_Sigma := sum_j log|Y_j| + sum_(i<j) log|t_ij| >= c w - O(1)`. Here `H_Sigma` has `7 + 21` terms.

**Normalisations, as fed by Steps 0-5:**
* *Scale.* `w = W/2^m = W/128` (Step 0).
* *Balance.* Step 1 delivers balance for every nonempty `T ⊆ [7]`, with `eta(M) -> 0`, so one fixed `eta_0`
  suffices. Imposing balance only on the 119 visible blocks (`2 <= |T| <= 6`) is equivalent:
  - singletons enter no minor;
  - the top block can be inflated by a fresh prime through `z -> nz`, which fixes `Phi` and adds a one-block chain.

  census §1.1 gives this argument; both the census and theory referees accepted it.
* *Exceptional set.* `Z` is quantified before the data. Step 2 removes it at cost `c(Z) M^7`, uniformly in `R`. No
  fixed-`S` exceptional set enters anywhere.
* *Residues and Gaussian costs.* `(L_8)` is form-blind. The Gaussian gcd and denominator costs are already paid in
  Step 0 and in the item-532 extraction behind Theorem D, and nothing this round discards them.
* *Anchor.* Anchor 0. The chain condition does not depend on the anchor (theory Lemma 1.1). The exceptional set does
  depend on the labelling, so a statement admitting an arbitrary anchor must symmetrise `Z` (census m1).

**Typo in the record (census finding, confirmed by its referee and by J1).** `ptolemy.md` §1 gives the number of
visible blocks as `2^m - m - 1`. The correct count is `2^m - m - 2`, which is 119 at `m = 7`, equal to the number of
boundary divisors of `Mbar_(0,8)`. Nothing depends on it.

## 2. Referee-corrected statements

Each statement below is the referee's corrected version, condensed. The full versions are in §8 of
`referee_census.md`, §8 of `referee_lift.md` and §5 of `referee_theory.md`.

### 2.1 census (`census.md`; report `referee_census.md`; outcome (c), survives)

Every formal result was re-derived. Every computation was reproduced exactly by independent code, including:
* a C census at `k = 5..8` (for example 3,646,069 configurations at `k = 5` with best `eps = 0.1333`, and the single
  `k = 8` configuration with 29/119 coverage);
* the symmetric-loci minima 3/6/25/48;
* the `k = 6` family table for `N = 2..14`;
* the lattice bounds, with the index checked exactly.

1. **Sandwich reduction (Prop 1.3; proved).**
   * `(L_8^c)(eta_1, Z)` for all `Phi` implies `(L_8)` with `eta_0 = eta_1/2` and `c = eta_1/4`.
   * `(L_8)` implies `(L_8^c)(eta_0, Z)` for `Phi` admissible from anchor 0 (any anchor with symmetrised `Z`).

   Here `(L_8^c)` says that balanced contact profiles (all 119 contact integers `d_T` with
   `|log d_T - w| <= eta w`) do not occur for large `w` off `Z`. **This is not an equivalence** (m2). Small-`H`
   data may have non-caterpillar trees at primes that divide residues.
2. **Theorem E (proved).**
   * Each prime divides at most `k - 3` visible blocks, so a balanced datum needs at least 5 / 9 / 14 / 24 primes at
     `k = 5..8`.
   * At `k = 8`, `eta = 0.1`: the class is empty for `e^w < 59.2`, and the 21 minors with `i, j >= 1` are
     `>= 10^49.4`. At `eta = 1/2`: `e^w >= 19.9`, minors `>= 10^20.1`.
3. **`(L_8)` on special loci, with locus-dependent constants (Prop 7.1, Cor 7.2, 7.3, Prop 7.4, 7.5(4); proved).**
   * Every subvariety missing some `D_T` carries no balanced points. This covers all symmetric loci `Y_sigma`, with
     at least 48 identically trivial contacts (exact over `S_8`), all forgetful fibres (48/80/100/112 bounded
     contacts), and the boundary.
   * On each fixed curve of genus `<= 1` with `K.C <= 0`, `(L_8)` holds for `eta_0 < 5/141`, with constants
     depending on `C`.
4. **Families (Prop 7.5(1)-(3), Cor 7.6; proved, with the M1 correction).**
   * Covering families of genus-`<= 1` curves satisfy `7 K.gamma >= a(70 - 182 eta(gamma))` and `K.gamma <= 0`,
     hence class imbalance `eta >= 5/13`.
   * So no bounded-degree covering family can have members that *each carry infinitely many* balanced small-`H`
     points. Criterion (b) of `round4/L7.md` Thm 2 never applies at `k = 8`, and criteria (a) and (c) need
     unbounded degree.
   * **Correction M1 (major, interpretive).** This does **not** exclude a Zariski-dense balanced supply with finitely
     many points on each member of a bounded-degree covering family. Every rational point lies on `Q`-lines of
     Kapranov's `P^5`. The sentences of census §6.3, §7.2 (end) and §8(a) that say otherwise are withdrawn. They
     treated a curve-dependent threshold as uniform, against the owner's constraint.
   * Further corrections: "rigid" should read "non-covering" (m5). The statement "`Z` contains every class-`d beta`
     curve with infinitely many points" cites the wrong direction (m3). What is true is that `Z` must contain every
     curve carrying infinitely many balanced small-`H` points, which is `round4/L7.md` Thm 2(a) (genus 0, under its
     hypotheses (i)-(iv)) or Thm 2(c) (genus 1, sketched).
5. **Lemma 7.7 (proved).** `(L_k)` with `(c, eta_0, Z_k)` implies `(L_(k+1))` with `(2c, eta_0, pi^(-1) Z_k)`.
6. **Exact data and heuristics.**
   * The exhaustive boxes reach at most 18/25, 23/56 and 29/119 boundary divisors at `k = 6, 7, 8`.
   * On the `k = 6` family `eps -> 0` with `H_Sigma = O(1)` (`eps` 0.39 -> 0.078, not monotone).
   * CRT/lattice upper bounds `H <= ~8w, 42w, 150w, 465w` at `k = 5..8`.
   * The heuristic success rates `e^(-5w), e^(-17w), e^(-49w), e^(-129w)` explain why block sampling is blind to the
     `k`-threshold.
   * Minor corrections: m4 (the cost of moving `Sigma`-parts is `<= 6·79`, still bounded), m5 (21% is full coverage,
     not balance), m6 (a translation is needed in Lemma 1.1).

### 2.2 lift (`lift.md`; report `referee_lift.md`; partial, survives)

1. **Lemma 1.1 (proved).** Restricting to a sub-configuration on `K ∋ 0` keeps the residues and multiplies the
   blocks in groups of `2^(k-k')`, with the same `eta`. So an 8-lift at scale `w` gives 7-lifts at scale `2w`, and
   forces `Phi_6` itself to be `eta`-balanced.
2. **Theorem 4.2 (proved; exact over `Q` at 9 parameter points plus openness).** Let `t` lie in a nonempty
   Zariski-open `U_0 ⊆ U ∩ {e_3^2 - e_1 e_2 e_3 + e_1^2 e_4 ≠ 0}`. `U_0` contains `t_0 = (2,3,-5,1)`, item 109's
   `(1/3, 3, 2, 5/3)` and 7 other points.
   * (a) For `L ∈ Pic^8(E_t ⊗ Qbar)`, `L ≠ O(2H)`, no nonzero section of `L ⊠ O(2)` vanishes at the 25 points
     `(B_A, v_A)`. For `L = O(2H)` these sections are exactly the pencil `h^2 P_lam - a^2 Q_lam`.
   * (b) Every irreducible class-`beta_7` curve meeting `M_(0,7)` and mapping onto `E_t` is a pencil lift
     `C'_lam`, of geometric genus 9 (or 5 for at most two `lam`). **No class-`beta_7` lift of genus `<= 4` exists.**

   The referee re-proved (a) by a shorter route: the 2-dimensional `W = Pi * Pencil`, and symbolic coplanarity of
   `aux_1, aux_2, aux_3, B_bc`.
3. **Corollary 4.3 (proved).** A class-`beta_8` lift of genus `<= 4` must have both 7-point projections birational,
   of class `2 beta_7` and degree 4 over `E_t`. Pencil `beta_8` lifts exist, of genus 33.
4. **Correction (major).** The headline "the `E_t` mechanism does not lift, already at the first step" holds only
   for class `beta_7`. Genus-one positive-rank lifts of class `c beta_7`, `c >= 2`, are the same mechanism and are
   **open**. So are all point-level lifts.
5. **Counts (heuristic; arithmetic correct).** `-9` per `E_t` for one added point, `-18` for class `2 beta_7`, `-42`
   for `beta_8`; integer analogues `e^(-9w)` and `e^(-42w)`. The bookkeeping reproduces `-K.beta_8 = -10`.
6. **Data (evidence; nearly vacuous).** For `N <= 7` the configurations `Phi_N(t_0)` are not themselves 1/5-balanced
   (`eta_min = 1.00, 0.46, 0.44, 0.22, 0.26`). By Lemma 1.1 that alone excludes 1/5-balanced lifts. At `eta = 1/5`
   the skew test needs negative skew, so it excludes nothing. Only `N = 8` (`eta_min = 0.157`, 9 unsplittable
   contacts, assuming `Sigma`-free lift blocks) carries information.

### 2.3 theory (`theory.md`; report `referee_theory.md`; partial, survives)

1. **Lemma 1.1 (proved; exact on 26/236/2752/39208 trees, four criteria, 0 mismatches).**
   * For split systems, "two nested chains on disjoint row sets" for one anchor, for every anchor, and the
     caterpillar condition are equivalent.
   * Applied to an `(L_8)`-datum, this is the caterpillar condition on the *block system* at each prime. It equals the
     `p`-adic tree of `Phi` only at primes dividing no residue (item H).
2. **Prop 2.1 (proved as a reduction).** Every proof of `(L_8)` must prove `(L_8^pc)`, the pairwise-coprime case.
   There the chain, anchoring and residue hypotheses hold automatically. The gloss "no structure other than
   `K.beta = 10` survives" is not a theorem (item G).
3. **Rate statement (proved, as corrected; item A).** For any `f -> infinity`:
   `(L_8^c)^all => (L_8) => (L_8^f)` (for `f <= c w`) `=> (L_8^c)^cat`. So no slowly-growing variant avoids
   excluding dense caterpillar balanced points. It is **not** shown to be equivalent to `(L_8)`.
4. **Lattice form (Lemma 2.2, Prop 2.3; proved; index checked by HNF on 80/80 configurations).**
   * `[Z^7 : Lambda(x)] = prod_(|T|>=3) d_T^(|T|-1)`.
   * `kappa = sum (|T|-2) ell_T - 6 log|x|` is the **finite part** of `h_K` (item D).
   * A balanced point is a vector of its own congruence lattice at `e^(31w)`, against the Minkowski scale
     `e^(49w)` and the trivial bound `e^(30w)`.
   * At positive lattice density `kappa <= -sum log e_T`, so density, box and sieve mechanisms (the `(L_5)`
     mechanism) cannot produce balanced points at `k = 8`. This holds in the Kapranov chart, and by the referee's
     hand count in the integral chart too.
   * Heuristic: `98 - 108 = -10 = -K.beta`.
5. **Lemma 2.4 and Theorem 2.5 (proved).**
   * Lemma 2.4: matching coefficients are pairwise coprime.
   * Theorem 2.5 (Ridout): for a fixed prime set `S_0`, `H_Sigma >= 8(1-eta)(1-eps) w - C(S_0, eta, eps)` when all
     blocks are `S_0`-units.

   Because `C` depends on `S_0` itself, Theorem 2.5 is **useless for the circle problem**, as stated. It respects
   the owner's first constraint.
6. **abc lane (conditional on abc).**
   * Prop 2.6 (proved): a unit pencil `lambda A/C`, `lambda ∈ Q^*`, with
     `gain = L.beta - F_red.beta = 10 - K^log_rel.beta > 0` gives `(L_8^pc)` off `{B = 0}`.
   * Equivalently, `gain = sum_Y ((m_Y - 3)/3) Y.beta` over the components of the three fibres.
   * Morphisms have gain `<= -32` (Bruno–Mella, **read from the abstract only**).
   * Data, exact; a rigorous upper bound for the pencils enumerated. With `lambda = ±1` and `D <= 4`: 277,940 pencils,
     maximum gain `-32`. All special rational `lambda` with `D <= 3` (17 values, referee): maximum `-32`.
     `D = 4` with `lambda ≠ ±1`, and every `D >= 5`, are untested (item F).
   * The "abc + auxiliary pencil" gap in the barrier map is genuine. The refereed abc barrier
     (`conditional.md` Prop 8.3, `round8/weakest.md` Prop 4.2) covers only relations among products of differences.
7. **Refutation lane.**
   * Expected dimensions `-10d + 2` (genus 0) and `-10d` (genus 1); the first was already in `L8.md` §7.
   * No genus-`<= 1` curve carrying infinitely many balanced points belongs to a covering family.
   * The computed rigid curves with `K > 0` (F-curves of types 2222 and 3221) lie in the boundary.
   * **Correction B (major, interpretive).** A refutation of `(L_8)` refutes `V_8` at balanced points, but **not
     conversely**. Refuting `(L_8)` needs a Zariski-dense supply of balanced points that are caterpillar at every
     prime up to `o(w)` of residue. Dense non-caterpillar balanced points would refute only `V_8`. The "iff" of
     theory §3.1 and the "exactly" of §3.4 and §5(1) are withdrawn.
   * Item C: genus-`>= 2` covering families are not "excluded". Item E: the closing split into curves with infinitely
     many points and sporadic points is a description, not a theorem.

## 3. Status of (L_8)

### 3.1 Verdict: open

| possible status | outcome this round |
| --- | --- |
| proved | **no.** No proof of `(L_8)`, of `(L_8^pc)`, or of any `f(w) -> infinity` version. |
| refuted | **no.** No dense balanced supply is known, at `k = 8` or even at `k = 7`. |
| restricted theorem | **only on special loci, with locus-dependent constants** (§2.1 item 3), or for a fixed prime set `S_0` with `C = C(S_0)` (theory Thm 2.5). Neither is a statement on a Zariski-dense set with uniform constants, so neither is a piece of Step 4. |
| conditional | `V_8 => (L_8)` (ptolemy Thm D, refereed earlier). New: abc plus a unit pencil of positive gain `=> (L_8^pc)`, and no such pencil is known. |
| **status** | **open** |

### 3.2 The exceptional set `Z`: no concrete candidate

The round-17 census was meant to locate `Z`. It did not, and Theorem E explains why a search could not: balanced
8-configurations begin at minors of size about `10^49`. What is now known about `Z`, all refereed:

| locus | must `Z` contain it? | source |
| --- | --- | --- |
| subvarieties missing some `D_T`: symmetric loci, forgetful fibres, boundary strata | **no**. They carry no balanced points, so they are harmless whether or not they are in `Z`. | census Prop 7.1, Cor 7.2, 7.3 |
| a fixed curve of genus `<= 1` with `K.C <= 0` | **no**, if `eta_0 < 5/141` (only finitely many balanced small-`H` points on it) | census Prop 7.4, 7.5(4) |
| pull-backs of the `(L_6)` curves `E_t` through class-`beta_7` lifts | **no**: such lifts have genus 9 or 5, so finitely many rational points | lift Thm 4.2, Faltings |
| curves carrying infinitely many balanced small-`H` points | **yes**. Such a curve has class `∝ beta` (census Prop 7.4). It lies in no covering family (census Prop 7.5). None is known; expected dimension `-10d + 2` or `-10d`. | census Prop 7.4, 7.5; `round4/L7.md` Thm 2(a), (c) |
| genus-one lifts of class `c beta_7` (`c >= 2`) over `E_t`, and their `beta_8` combinations | unknown (expected dimensions `-18` and `-42`) | lift Cor 4.3, §6.1 |

So the only shape `Z` is forced to have is "non-covering curves of genus `<= 1` and class `d beta` with infinitely
many rational points". Whether any exists is open, and the expected dimension is negative. **The first concrete
candidate `Z` asked for in round 17 does not exist in the record**, and the record contains no method that would
produce one.

### 3.3 The logical map around `(L_8)` after round 18

```text
(L_8^c)^all  =>  (L_8)  =>  (L_8^f) (f <= c w)  =>  (L_8^c)^cat  =>  (L_8^pc)              [census 1.3, theory A, 2.1]
V_8          =>  (L_8)                                                                    [ptolemy Thm D]
(L_8)        =>  uniform theorem                                                           [ptolemy Steps 0-5]
(L_7)        =>  (L_8)                                                                     [census 7.7]
abc + unit pencil with gain > 0  =>  (L_8^pc)                                              [theory 2.6]
not (L_8)    =>  dense balanced caterpillar points (up to o(w) residue)  =>  not V_8      [theory B]
```

*README remark (structural, unrefereed; it uses only the refereed implications above).* The abc lane ends at
`(L_8^pc)`, the bottom of the first chain. Theorem D needs `(L_8)` itself. A positive-gain pencil would therefore
give the uniform theorem conditionally on abc only after a further step: the gcd and radical bookkeeping of
Prop 2.6 would have to be extended from pairwise-coprime profiles to nested chains with residues. theory open
point (b) names this step, and nobody has done it. For nested blocks, `gcd(A(x), C(x))` can exceed
`prod d_T^(mu_T)` by a positive multiple of `w`, which lowers `h(f(Phi))`. So the extension is a genuine gap, not
bookkeeping.

### 3.4 Why the round's evidence is weak, and cannot become strong

* **No feasible computation can test `(L_8)` (census Thm E, proved).** At `eta = 0.1` a balanced datum has all 21
  minors `>= 10^49.4`. Exhaustive boxes reach only 29/119 boundary divisors.
* **Block-sampling searches return `H ≍ w` at every `k`, including `k = 5, 6`**, where bounded-`H` points are known
  to be abundant. So a search failure at `k = 8` says nothing.
* **The structural tests passed** (density mechanisms, covering families, symmetric loci, forgetful fibres,
  class-`beta_7` lifts). Each excludes one *mechanism* of refutation. Together they confirm that `(L_8)` is not
  mis-stated in a way that the known `(L_5)` and `(L_6)` mechanisms exploit. They do not bear on whether `(L_8)` is
  true. Its truth remains equivalent in strength to an unconditional `D = 0` case of Vojta's conjecture, on a
  5-dimensional rational variety, for unboundedly many primes.

## 4. Consequences for the uniform bound and the growth rate

**Unchanged.** `M <= (2/3 + eps) log R/loglog R` remains the best bound for actual circles. No result of this round
touches actual clusters or the rate `log R/loglog R`. The uniform bound is still known only conditionally, on
Vojta's conjecture for one fixed variety: `Bl_Y X_7`, or `Mbar_(0,8)` via `(L_8)`.

**Conditional theorems: changed or added rows only.**

| hypothesis | conclusion | implication | hypothesis status |
| --- | --- | --- | --- |
| abc + a unit pencil `lambda A/C` on `Mbar_(0,8)` with `gain > 0` | `(L_8^pc)` | proved (theory Prop 2.6) | **open**; morphisms have gain `<= -32`, all tested pencils (`D <= 4` with `lambda = ±1`; `D <= 3` with any `lambda`) have gain `<= -32` |
| `(L_8^pc)` | uniform theorem | **not proved**: Theorem D needs full `(L_8)` (§3.3) | — |
| `(L_7)` | `(L_8)`, hence uniform | proved (census Lemma 7.7) | open; predicted false (`K.beta_7 = -7`) |

**Rubric (b), fake-violated input: unchanged.** `(L_8)` is a statement about integer `2 x 8` matrices. Fake circles
do not test it, and no input this round is shown to be true for circles and false for fakes. The work is about the
one lemma that is predicted true by a standard conjecture and blind to fakes. It is not about the fake barrier.

**What the record says about growth after round 18** (R4-R18, lit18). The round-17 §4 list stands. Two additions:
* **Ptolemy route, unconditional:** open with no method in either direction. No feasible test exists (Thm E). No
  explicit `Z` exists in the record. The `(L_5)`/`(L_6)` refutation mechanisms provably fail at `k = 8` (density:
  theory Prop 2.3; covering families: census Prop 7.5; class-`beta_7` lifts: lift Thm 4.2).
* **Barrier map, one genuine gap:** abc applied to auxiliary pencils with high-multiplicity fibres. This is not
  covered by the abc barrier (products of differences only) or by `L8.md` Thm A (Liouville only). It is closed for
  morphisms and for low-degree pencils (data). Even if it opens, it yields `(L_8^pc)`, not `(L_8)`.

## 5. Judge's own checks (unrefereed)

Scripts are in `judge_checks/`, and each writes `*_out.txt`. All arithmetic is exact (`Fraction`, integers);
logarithms are used only to report sizes. **None of this has been refereed.**

* **R18-J1 (`j1_combinatorics.py` -> `j1_out.txt`): exact reproduction of the combinatorial constants the round
  relies on.** All agree with the write-ups and reports.
  * Visible blocks: `2^m - m - 2 = 10, 25, 56, 119` at `k = 5..8`, equal to `#D_T = 2^(k-1) - k - 1`. This confirms
    the `ptolemy.md` §1 typo.
  * `K.beta` from `K = sum_s (s(n-s)/(n-1) - 2) D_s`:
    * `n = 8`: coefficients `-2/7, 1/7, 2/7` on 28, 56, 35 divisors give `K.beta = 10`;
    * `n = 7`: coefficients `-1/3, 0` on 21, 35 give `K.beta = -7`;
    * Kapranov form: `-6·31 + 196 = 10`, with `H.beta = 31`.
  * Theorem E: `r_k = 5, 9, 14, 24` and `p_(r_k) = 11, 23, 43, 89`. At `eta = 0.1`, emptiness for
    `e^w < 8.8, 17.3, 30.5, 59.2` and minors `>= 10^2.6, 10^7.8, 10^20.0, 10^49.4`. At `eta = 1/2`, `k = 8`:
    `e^w >= 19.9` and minors `>= 10^20.1`.
  * Prop 7.5 threshold: `70 - 182·(5/13) = 0`, and `-56 + 126 = 70 = 7 K.beta`.
  * Forgetful fibres: 48, 80, 100, 112 vertical (bounded) contacts over `Mbar_(0,q)`, `q = 4..7`, by enumeration,
    equal to `(2^(q-1) - q - 1) 2^(8-q)`.
  * Cross-ratio pencil: three fibres of 16 boundary divisors each, so `L.beta = 16`, `F_red.beta = 48`,
    `gain = -32`. This is the value at which all of theory's enumerations stall.
* **R18-J2 (`j2_gcd_excess.py` -> `j2_out.txt`): the gap in §3.3 is real.**
  * Setup: unit pencil `A = (x_1 - x_2)(x_4 - x_5)`, `C = (x_1 - x_3)(x_2 - x_3)`, with one prime `p = 10007`
    carrying the nested chain `{1,2} ⊂ {1,2,3}`: `p^2 | x_1 - x_2` and `p || x_1 - x_3`.
  * Result: `v_p(gcd(A(x), C(x))) = 2`, whereas `sum_T mu_T v_p(d_T) = 1`. The gcd exceeds Prop 2.6's value by a full
    block `d_12`.
  * Reading: on caterpillar profiles, `h(f(Phi))` can fall below `(L.beta) w` by a positive multiple of `w`. So
    extending Prop 2.6 from `(L_8^pc)` to `(L_8)` is not routine. (The radical also counts `p` once, which might
    compensate. That is not analysed.)
* **No other new mathematical claims.** Everything else in this README is a referee-corrected statement or a reading
  of the record.

## 6. Scores

Rubric of `docs/claude_review/round10/README.md` §1, each sub-score out of 10:
* (a) correctness after refereeing;
* (b) fake-violated input;
* (c) movement of growth or uniformity.

The overall score is a judgement, and it weighs (b) and (c) above (a). Nothing this round moves growth or
uniformity, so no (c) exceeds 2. Partial (c) credit is given only for a proved correction to the barrier map.
(b) is 1 throughout: `(L_8)` is fake-blind by design, and no input shown true for circles and false for fakes was
produced.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **census** | 7 | 1 | 1 | **5** | It did its assigned task as far as the task can be done. It proves (Thm E) that the literal census is infeasible, then replaces it with exact statements: the Prop 1.3 sandwich, the loci theorems Prop 7.1/Cor 7.2/7.3, the families threshold `eta >= 5/13` (Prop 7.5), and Lemma 7.7. Its calibration at `k = 5, 6` shows that block sampling is blind to the threshold. Every number was reproduced exactly by independent C and Python code. Deductions: the major interpretive overclaim M1 ("no bounded-degree family produces balanced points densely"), which treats a curve-dependent threshold as uniform, plus six minor items. It did not locate `Z`, and Thm E explains why no search could. |
| **theory** | 7 | 1 | 2 | **5** | The round's one new mechanism-level finding is a genuine gap in the barrier map (abc on auxiliary pencils). It comes with an exact criterion (`gain = 10 - K^log_rel.beta`), a morphism bound, and an enumeration of 277,940 pencils, upheld entry by entry. The caterpillar lemma, the lattice form and the reduction to `(L_8^pc)` sharpen what any proof must do. Deductions: the major "iff / exactly a refutation of `V_8`" overclaim (B), the "all unit pencils" overclaim (F, only `lambda = ±1`), and six wording overstatements. The (c) credit is for the barrier-map correction. The lane ends at `(L_8^pc)`, not at the uniform theorem (§3.3). |
| **lift** | 6 | 1 | 1 | **4** | Theorem 4.2 is a clean, unconditional, exact classification. The referee re-proved it by a shorter route. It decides the cheapest curve-level version of round-17 task 2. Deductions: the headline ("does not lift") holds only for class `beta_7` (major); the §5 integer data are nearly vacuous (the `Phi_N` are not themselves balanced for `N <= 7`, and the skew test is void at `eta = 1/5`); `U_0` is not explicit. The point-level question, which is what `(L_8)` is about, is untouched, and the write-up says so honestly. |

## 7. Single most promising next step (whole record)

### 7.1 Where the record stands

The barrier constraints of `round17/README.md` §7.1 (items 1-6) all stand. Round 18 adds:

7. **No refutation mechanism of the `(L_5)`/`(L_6)` type exists at `k = 8`.** Density mechanisms fail (theory
   Prop 2.3). Covering genus-`<= 1` families with infinitely many points per member fail (census Prop 7.5).
   Class-`beta_7` lifts of the `(L_6)` curves fail (lift Thm 4.2). What survives as a possible refutation is a
   dense supply of caterpillar balanced points with finitely many on each curve of every family, or infinitely many
   non-covering genus-`<= 1` curves of class `∝ beta`. Either would also refute `V_8`.
8. **No computation can test `(L_8)`** (Thm E: minors `>= 10^49` at `eta = 0.1`). Sampled searches return
   `H ≍ w` at every `k`.
9. **The only uncovered conditional mechanism, abc on auxiliary pencils, stops at `(L_8^pc)`.** Extending it to
   `(L_8)` meets a genuine gcd excess on nested chains (J2).

### 7.2 The candidates, re-rated

| candidate | would give | method known? | judge's estimate (one further internal round) |
| --- | --- | --- | --- |
| proof of `(L_8)` | the uniform bound (unconditional) | no; every unconditional technique in the record is excluded by a refereed barrier | **`< 0.1%`** (lowered from `< 1%`: round 18 removed the hope that `(L_8)` reduces to something weaker; see the sandwich, Prop 2.1 and the rate statement) |
| refutation of `(L_8)` | would show the Ptolemy route is mis-specified; it would also refute `V_8` | no; every family mechanism fails (7.1 item 7) | `~0` |
| explicit `Z` (non-covering genus-`<= 1` curves of class `d beta` with infinitely many points) | the first piece of the exceptional set; no growth | no; expected dimension `-10d + 2`; class-`2 beta_7` lifts need `3^25` incidence patterns | `< 2%`, and worthless for growth |
| **abc-gain search**: unit pencils with `D = 5` (all rational `lambda`), plus `D = 4` with `lambda ≠ ±1` | if `gain > 0` is found: abc `=> (L_8^pc)`, and a uniform theorem from abc only after the unproved nested-chain extension (§3.3, J2) | yes, finite and exact (theory/referee code; `D = 4` needed 66 processes, so `D = 5` is about 10x that, under the memory cap) | finding `gain > 0`: **`~3-5%`** (data stuck at the morphism value `-32` for all `D <= 4`). An abc-conditional uniform theorem in the same round: `< 1%`. Unconditional growth: `0` |
| `(L_7)` refutation (round 7 Q1-Q3) | nothing for `(L_8)`: by Lemma 7.7 the implication goes `(L_7) => (L_8)` only | finite computations | refutation plausible; value `~0` for growth |
| route 3, route 5, `QS^prod`, composite-supported analytic inputs | as in `round17/README.md` §7.2 | no | unchanged (`< 0.1%` to `~0`) |

### 7.3 Recommendation

**Stop internal rounds on growth and uniformity.** Eighteen rounds have produced a single growth-relevant proved
theorem, the `(2/3 + eps)` constant. They have produced a refereed barrier map that excludes every unconditional
mechanism the program knows. The last target with all three good properties (it implies the uniform theorem, it is
predicted true by a standard conjecture, and fakes do not test it) has now been examined by census, lifting and
theory. The outcome:
* no method in either direction;
* no explicit exceptional set;
* a proof that it cannot be tested computationally.

Further internal reformulation rounds are very unlikely to change this.

**The single most promising next step is external.** It costs one short writing task:
1. Append a round-18 addendum to `docs/claude_review/overview/overview.md` §7.1 and §6.6. The "smallest useful
   experiments" paragraph there describes exactly what round 18 did, so replace it with the results:
   * Thm E (no feasible census);
   * the sandwich `(L_8^c)^all => (L_8) => (L_8^c)^cat => (L_8^pc)`;
   * the lattice form of `(L_8^pc)` (theory §5(2) and open point (a): short vectors at `e^(31w)` in a lattice of
     Minkowski scale `e^(49w)`, trivial bound `e^(30w)`);
   * the loci on which `(L_8)` holds (census Prop 7.1, 7.4, 7.5, corrected per M1);
   * lift Thm 4.2;
   * the abc-gain question: does any pencil on `Mbar_(0,8)` have `K^log_rel.beta < 10`?
   * Also fix the `2^m - m - 1` typo in `ptolemy.md` §1.

   The addendum must carry referee corrections M1, B and F verbatim.
2. Put `(L_8)` before outside experts in Diophantine approximation. The relevant expertise is in gcd inequalities
   and Vojta's conjecture for blow-ups of projective space, and in abc over function fields and number fields. The
   record has nothing more to contribute without a new idea, and lit18 found no applicable tool.

**If the owner wants one more internal round anyway,** the only item with a finite, exact, decisive outcome is the
abc-gain search of §7.2, with this stop rule. If every unit pencil with `D = 5` (all `lambda`) and `D = 4`
(`lambda ≠ ±1`) has gain `<= -32`, close the abc lane and stop. If one has positive gain, the next task is the
nested-chain extension of Prop 2.6, which J2 shows is not routine. Run it under the memory rules: one process per
`A`-orbit and `C`-chunk, `prlimit --as=3500000000`, results streamed to disk.

**Honest assessment.**
* Unconditional growth `M = o(log R/loglog R)` from any further internal round: essentially `0`.
* The uniform bound unconditionally: essentially `0`.
* A conditional uniform theorem from abc instead of Vojta, via the pencil lane: under 1% in one round. Even then it
  is conditional, so it would not count as growth progress.
* The record's lasting value is the `(2/3 + eps)` theorem, the conditional theorems and the barrier map. All are
  already consolidated in `overview/overview.md`.

**Not recommended:**
* any further `(L_8)` census or search, which Thm E makes futile;
* lifting `E_t` at class `c beta_7`, `c >= 2`, which is infeasible exactly and could at best supply a piece of `Z`;
* proof attempts on `(L_8)` by Liouville, Subspace or density arguments, all excluded;
* `(L_7)` work;
* everything in the round-17 "not recommended" list.

## 8. Records to correct; files

**Records to correct.**
* `growth/ptolemy.md` §1: the number of visible blocks is `2^m - m - 2`, not `2^m - m - 1` (census; J1).
* `overview/overview.md` §7.1: replace "smallest useful experiments" with the round-18 results (§7.3 above).
  §6.6 / route table: chance of a proof of `(L_8)` `< 0.1%`, and add the abc-pencil gap to §5.5 ("not capped").
* `round17/README.md` §7.3 called the census "likely to produce usable evidence: an explicit candidate `Z`". It
  produced none, and census Thm E shows a census cannot. Record this as a finding, not an error.
* census.md §6.3, §7.2 (end), §8(a); theory.md §3.1, §3.4, §5(1), §2.4(2); lift.md §0 and §6.1 headline: read them
  only with referee corrections M1, B, F and the lift "class `beta_7` only" correction.

**Files (this directory, `round18/`).**
* Write-ups:
  * `census.md` + `census_checks/`, with notes in `census_notes.md`;
  * `lift.md` + `lift_checks/`;
  * `theory.md` + `theory_checks/`, with notes in `theory_notes.md`.
* Referee reports: `referee_census.md`, `referee_lift.md`, `referee_theory.md`, with their `referee_*_checks/`.
* Judge: this `README.md`; `judge_checks/j1_combinatorics.py` -> `j1_out.txt`; `judge_checks/j2_gcd_excess.py` ->
  `j2_out.txt`; `judge_checks/notes.md`.
