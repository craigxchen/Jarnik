# Round 20 judge: fresh-mechanism panel II (four lenses, four kill reports)

Four write-ups this round. An adversarial killer attacked each one, and each report sits next to its write-up:

| key | write-up | kill report | lens |
| --- | --- | --- | --- |
| gauss | `gauss.md` + `gauss_checks/` | `kill_gauss.md` + `kill_gauss_checks/` | full Gaussian model at small `k`: gcd reduction, `(L_k^circ)`, `A_k(delta)` |
| arakelov | `arakelov.md` + `arakelov_checks/` | `kill_arakelov.md` + `kill_arakelov_checks/` | arithmetic intersection theory on the model of `uv = N t^2` |
| composite | `composite.md` + `composite_checks/` | `kill_composite.md` + `kill_composite_checks/` | composite-supported analytic inputs (the last unclassified row of overview §5.5) |
| wild | `wild.md` + `wild_checks/` | `kill_wild.md` + `kill_wild_checks/` | unlikely intersections in tori (BMZ, Habegger, Rémond, anomalous Zilber-Pink) |

Paths are relative to this directory (`round20/`, scratchpad). `overview` means
`docs/claude_review/overview/overview.md` in the repository, which was read only. Labels: proved / conditional / exact /
evidence / heuristic / unrefereed, as in the overview. **Nothing written this round has been refereed by anyone other
than the killers.** Statements in §1 are the killer-corrected ones; remarks added by this README are marked *judge*
and are **unrefereed**. The judge's computations are in `judge_checks/` (§4).

## 0. Bottom line

**No growth improvement and no uniform bound were proved this round,** either unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for actual circles. The
uniform theorem is still known only conditionally on Vojta. **All four lenses were killed; no live lead survives.**

What the round did change, all at the level of the barrier map or the conditional lane:
1. **gauss** (killed as a route; the most substantial mathematics of the round). Proved: `A_k(delta) => uniform`
   (Thms 3.4, 3.4'), a refutation of the 4-point (and, by the killer, 3-point) version (Thm 4.1), and a new non-vacuous
   Vojta instance on a 4-fold, `V_E(5)` (`D = E`, `S = {inf}`), that implies the uniform theorem (Thm 4.2,
   conditional). Killed because the needed statement `(L_5^circ)` is exactly Liouville-critical *by invariance*
   (killer Lemma K1: dividing by the tuple gcd changes nothing any method in the record can see), the target point
   is rational (so Roth/Schmidt do not apply), fakes violate it, and only Vojta predicts it. Same standing as `(L_8)`.
2. **arakelov** (killed; nothing survives as a lead). The dictionary (Thm A) is exact. Every inequality the lens
   produces is row 2 (pairs, axiom (R)), row 8 (`tau*`) or Lim (rank / `|S|`) of overview §5.4. The residual needed
   statement is `(L_8)` / Vojta, which is at least as strong as the target.
3. **composite** (the cap survives, the classification claims are killed). **Robust composite-supported inputs are now
   capped** (Thms M, P; kind F; R1 caveat). The non-robust part stays unclassified, as at least four items. The judge
   narrows two of them (§4, J1: unrefereed): smoothed-norm, sharp-angle sector counts used **without**
   multiplicativity are capped by a "swap extension", and so are smooth Hecke sums without multiplicativity at
   precisions coarser than `X^{(1+eps)/A - 1/4}`.
4. **wild** (killed). Unlikely intersections give nothing. The fixed carriers are vacuous (`V^oa` is empty), the
   moving-line comparisons are false, and the subgroup-membership ceiling is `omega(N')` (= lit18 Prop B). Its one
   open residue `(S_min)` is **false for every `n`**, already on genuine short-arc 3-clusters (killer Prop K.1,
   re-checked exactly by the judge, J2).

**Scores** (rubric of `round10/README.md` §1; §6): gauss 5, composite 4, wild 3, arakelov 2.

**Recommendation** (§7): stop internal rounds on growth and uniformity, as round 19 already recommended. This round
tested the four remaining "fresh lens" families, and each reduced to an existing cap or to the Vojta lane. The one
useful addition to the external package is gauss's 5-point Vojta instance `V_E(5)` with its `k = 4` refutation:
it is the smallest-dimensional variety on record whose Vojta instance gives the uniform theorem.

## 1. Killer-corrected outcome of each lens

### 1.1 gauss (`gauss.md`; `kill_gauss.md`): killed as a route; four lemmas survive

**Setting.** A primitive `k`-tuple `zeta_1..zeta_k` of lattice points on a circle of radius `R`, up to similarity, is a
point of the norm-one torus `T_k` inside the twisted Losev-Manin space (the twisted `Mbar_{0,k+2}` with the two
isotropic marks of `threshold/L8.md` §6.1). The Gaussian block norms are its toric heights. A short arc is
archimedean contact with the identity `e` (the exceptional divisor `E` of `Bl_e`).

**Survives (killer-corrected).**
1. **Thm 3.4 / 3.4' (proved, elementary).** Let `k >= 3` and `0 < delta < 1/(2^k - 2)`. Define `A_k(delta)`: there
   are a proper Zariski-closed `Z ⊂ T_k` and an `R_0` such that every primitive `k`-tuple of distinct lattice points
   of a circle of radius `R >= R_0`, with all distances `<= R^{1/2 - delta}`, has class `[zeta] ∈ Z`. Then `A_k(delta)`
   implies `M(C) < inf` for every `C`. For `delta <= 2^{-k-1}` the proof uses only the convexity averaging of
   Lemma 3.2. Otherwise it uses item 532's balanced selection, **within one Gaussian unit class** (killer minor 1:
   across unit classes `t_ab` need not be an integer; 408 of 1866 tested pairs have `|t_ab| ∈ sqrt(odd/2)`).
   Also `A_k(delta) => (L_k^circ) => uniform`. `(L_k^circ)` must be stated within one unit class, or with
   `t_ab := |c_ab|/2` real and the bound `H_Sigma >= c w - O(1)`. The converse `(L_5^circ) => A_5` is **not** proved
   ("equivalently" in the proposal is withdrawn). The mechanism is that the tuple gcd shortens the arc: in the
   reduced radius `R''` the tuple lies on an arc `R''^{1/2 - 1/(2^k - 2)}`, i.e. `R''^{7/15}` at `k = 5`.
2. **Thm 4.1 (proved; exact rank reproduced at 6 rational points).** `A_4(delta)` is false for every
   `delta < 1/14`, and `(L_4^circ)` is false. An explicit rational 2-parameter family of Gaussian linear blocks gives
   balanced, bounded-residue 4-tuples at exponent `3/7`, dense in `T_4`. *Killer corollary (proved by projection,
   `k3_projection.py`):* `A_3(delta)`, `delta < 1/6`, and `(L_3^circ)` are false. So `k = 5` is the smallest `k` at which
   the statement can hold.
3. **Thm 4.2 (conditional).** Let `V_E(k, eps)` be Vojta's Main Conjecture for `(Bl_e Lbar_k^tw, D = E, S = {inf})`.
   It implies `A_k(delta)` whenever `(k-1)(1 + 2 delta)/4 > 1`. At `k = 5` this holds for every `delta > 0`, so
   **`V_E(5)` implies the uniform theorem**. The margin is `kappa = 2 delta` (`(K + D_PQ).gamma_hat = 1/16` at `k = 5`,
   recomputed exactly for `k = 3..9`). The `D = 0` version needs `k = 6`.
4. **Thm 4.3(ii) (proved).** For `k >= 5` there is no bounded-degree covering family of curves with infinitely many
   A-points on its members. Thm 4.3(i) is narrowed to: genus 0, approached point `O` real algebraic, contact order
   bounded below by Liouville (`O` rational) or Roth (otherwise).

**Killed.**
* **The headline "subcritical exponent".** Killer Lemma K1 (proved): for every bihomogeneous `F` of bidegree
  `(a, b)`, `F(z) = g^a gbar^b F(zeta)`, and the tropical forced divisor scales in the same way. Angles, pair
  weights, residues and block angles are unchanged. So every forced-divisibility, Liouville, pair, cut, character and
  single-form inequality holds for `z` iff it holds for `zeta = z/g`. On the `k` chosen rows the profile of a
  balanced cluster is the cube profile, so `construct.md` Cor C3 caps every auxiliary form at ratio exactly 1, and
  single-form Lang-Waldschmidt is exactly tight (`|Lambda| >= e^{-(1+eps)4w}` against a needed `e^{-4w}`). The move
  from exponent `1/2` to `7/15` is a change of units, not of difficulty. Only the heuristic class count (`e^{-w}` per
  scale) and Vojta's `D = E` term see the reduction, and neither is a proof method.
* **The needed statement as a route.** `(L_5^circ)` / `A_5(delta < 1/30)` asks that the Liouville exponent at the
  **rational** point `e` of the 4-fold `T_5` not be attained by a Zariski-dense set of balanced points. This fails at
  `k = 4` (Thm 4.1) and on `P^1`. Roth, Schmidt and Ru-Vojta beat Liouville only at irrational or Galois-split
  targets. Fakes violate the statement, so a proof would need genuine Gaussian angles at all primes of `N''` jointly,
  which no method in the record supplies. Its standing is that of `(L_8)`, with a thinner heuristic margin
  (`kappa = 2 delta ≈ 1/15` against `K.beta = 10`).
* Minor: the density mechanism is excluded at `k = 5` only heuristically (not by Thm 4.3); the phi-test of §2.1
  carries no weight (a form-blind point has a 2-parameter family of circle models); the ratio range for Thm 4.1 is
  `[8.31, 142.28]`, not `[8.3, 43.4]`; attribution: the model and the "Gaussian analogue with 6 points" pointer are in
  `threshold/L8.md` §6.1. The new content is the reduction written out, `D = E` (6 points to 5), Thm 4.1, and the
  circle analogue of round-18 item 7.

**Outcome.** Killed as a route. Four lemmas survive (one of them conditional). `(L_5^circ)`, `(L_6^circ)` and
`(L_7^circ)` are **open, Liouville-tight by invariance, fake-violated and Vojta-predicted, with no method in either
direction**.
### 1.2 arakelov (`arakelov.md`; `kill_arakelov.md`): killed; nothing survives as a lead

**Survives as correct mathematics (reformulation; no growth input).**
1. **Thm A (exact; re-verified independently on all 26 configurations).** On the regular model `X_N` of
   `uv = N t^2` over `Z[i]` (with `A_{e-1}` chains at each `pi | N`), the exponent levels are the components that the
   sections meet, `(P_i . P_j)_fin = log Norm c_ij` (reduced chords), the vertical correction of
   `D_0 = sum P_j - (M/2)(O_u + O_v)` is the cut slack `Q_p`, and the cut identity of the `(2/3 + eps)` proof is the
   genus-0 Hodge-index equality written place by place.
2. **Thm B (proved, as a reformulation).** The `(2/3 + eps)` finite inequality is the adelic `M`-point Fekete
   inequality **with the crude `(1+i)` bound** `binom(M,2) log 2`. The genuine 2-adic minimum is larger by a factor of
   about 1.9-2.4 in contacts (85, 2809, 48336 against 45, 1225, 19900 at `M = 10, 50, 200`). Constants only.
3. **Prop B.2 (proved, repaired).** The adelic-capacity class ADC is capped. It holds for `X >= 2M`, and the `(1+i)`
   place needs the `round5/lower.md` witnesses at `(48 + o(1)) M log M`: the cited `construct.md` witness has residues
   at odd prime powers only. ADC is row 2 of overview §5.4 restated in Arakelov language.
4. **Prop 3.3(b) (proved).** Zhang / Northcott / arithmetic bigness on a fixed model would need `tau*(M) > 2`. This is
   row 8 restated.
5. **Thm C (proved; near-tautological).** Heights on `Mbar_{0,M}` built from local data plus product formulas reduce
   to pair product formulas. With genuine places only it gives the per-pair inequality, i.e. row 2 again.
6. **Evidence.** The §4.2 slack table is reproduced. The §4.3 totals are 29 pairs (expected 25.67), not 28 / 25.5, and
   1 triple (expected 0.023). Cutting at `X = 2M` (not `M`, since an `A_X` fake designs every prime up to `X`):
   `Gamma^{>2M}_T = 0` for all `|T| >= 3` on all 15 arc clusters. Non-separation also follows structurally from Thm C.

**Killed.**
* The residual needed statement. Its precise form ("equivalently `(L_8)`, or Vojta `D = 0` on `Bl_Y X_7` /
  `Mbar_{0,8}`") **implies the uniform theorem** (overview §4.2, Thms 4.1, 4.3), so it is at least as strong as the
  target. Its vague form ("global coherence, weaker than lattice gaps at every prime, not fixed-`S`") names no
  inequality, and both of its boundary conditions are already in the record (round 10 Thm D; lit18 Props A, B).
* Every inequality the lens can produce lands in row 2, row 8 or Lim (Mordell-Weil packing gives counts
  `e^{O(rank)}`; bounded rank is equivalent to the uniform theorem, overview Thm 3.12).

**Outcome.** Killed. The dictionary is a clean and exact re-labelling of the `(2/3 + eps)` proof, worth one paragraph in
an expert account. It is not a route.
### 1.3 composite (`composite.md`; `kill_composite.md`): the cap survives; the classification claims are killed

**Survives (re-derived line by line; no mathematical error found).**
1. **Lemmas 1.2-1.4.** `S`-smooth masses `tau_sigma <= 8 s P_0^{-sigma}`; `Psi_S(X) <= binom(2s+L, L) <=
   X^{log(4es)/log P_0}`; sector lattice counts with error `<= 2.53 sqrt(y) + 1.58`; an all-ideal Gallagher large sieve
   with constant `<= 56(Y + K sqrt Y)`.
2. **Prop 1.5.** The optimal all-ideal large-sieve constant is `>= max(#I(X), (2K+1) floor(sqrt X))`, so every true
   all-ideal Hecke large sieve is robust. This holds for supports containing `I(X)` or the integers up to `X^{1/2}`.
3. **Thm M (minimal graft)** and **Thm P (multiplicative graft: P1, P2, P2', P3 with fixed `q`, P4, P5).** P1 answers
   `round10/referee_analytic.md` Slip 3.
4. **The cap (corrected constant, killer F1).** Fix `c, gamma, B, C_c` before the target constant. Take
   `A > max(4, 1/gamma)` for the `round8/lw.md` 6.1 fakes, with `W <= (2A + o(1)) M log M`, or `A > max(2, 1/gamma)` for
   the P3-char fakes, with `W <= (bA + o(1)) M log M`, `b >= 33`. Then for all large Paley orders `M`, the minimal graft
   satisfies every `(c, gamma)`-robust composite input (C1)-(C3), up to `1 + O(M^{1-A/2} + 2^B M^{1 - A gamma})`. The
   multiplicative graft satisfies P1-P5. So robust composite inputs, combined with the profile, the monomial gaps (or
   P3-char residues) and Haar-generic properties, cannot prove `M = o(log R/loglog R)`. **R1 caveat:** precision fixed
   before `c_0`. Here the caveat is not vacuous, because true composite inputs exist at every precision.
5. **Prop 3.1 (entire-L collapse; proved).** Smooth all-ideal Hecke sums with decay `O_f(X^{-delta_f})`, together with
   **complete additivity** of the angle function, force `phi_p = ± theta_p (mod pi/2)`: the fake is an actual cluster.
   This is the composite analogue of round 10 Thm D, i.e. the problem itself.

**Killed.**
* "Composite support adds no separating power" and "the only non-robust leftover is the sharp primitive large sieve,
  equivalent to the max-scale cross-gap".
  * **F2 (major).** With multiplicativity, primitive support contains every mixed product `d beta`, and the sharp large
    sieve then demands the **product-scale** gap
    `|sin(psi_d + arg beta - arg alpha')| >= c/max(Nd Nbeta, N alpha')`. This is LCG-shaped, not max-scale. Exact
    spacing at `p = 101`, `tan psi = sqrt7 - 2`, `X = 10^6`: `X · minspacing` is `1.003` genuine, `0.228` minimal
    graft, `5.6e-5` multiplicative graft. Heuristically only lattice directions pass, i.e. a collapse.
  * **F3 (major).** Smoothed-norm, sharp-angle sector counts with badly approximable (quadratic-irrational) edges are
    true with error `O(1)` to `O(log X)` (evidence to `10^7`; rational edges `0.0225 X^{1/2}`). The minimal graft
    violates them (mean `|err|` 2.30 against 0.141 genuine at `X = 1.5·10^6`). Neither Thm M nor Thm P covers them.
  * **F4 (major).** Smooth Hecke sums used **without** complete additivity are not shown to collapse: Prop 3.1 needs
    complete additivity.
* The §5 "needed statement" (position-uniformity required; a sharp/smooth dichotomy) is mis-specified. Position
  uniformity is not necessary, and F3 is neither sharp nor smooth.

**Judge narrowing (§4, J1, unrefereed).** F3 used **without** multiplicativity is capped by a swap extension, at every
precision `>= X^{2/A - 1/2 + o(1)}` (proof below, plus an exact-data check that the swap restores the genuine error
to within 0.0015). F4 used without multiplicativity is capped at precisions `>= X^{(1+eps)/A - 1/4}` (proof sketch). F4
at super-polynomial precision is **not** capped, and the evidence shows the swap does not reach it (`|F_k|` 0.016-0.039
after the swap, against `~1e-12` genuine).

**Outcome.** A proved correction to the barrier map: robust composite inputs move to the capped table (kind F). The
non-robust part stays in §5.5 as four sub-items (§3 below), none of which is known to give growth.
### 1.4 wild (`wild.md`; `kill_wild.md`): killed

**Survives as barrier facts (every step re-derived by the killer).**
1. **Thm 6.2.** The fixed carriers have empty open-anomalous locus: `X_M`, the chord space `L_M` (`M >= 4`), and
   `M_{0,k}` in its Plücker torus (`k >= 5`). Habegger / anomalous Zilber-Pink say nothing on them. The cross-ratio
   characters have rank `k(k-3)/2` (exact, `k = 4..9`).
2. **Lemma 6.3, Prop 6.4.** The moving carrier `ell_u = span_Q(Re u, Im u)` is never in a coset. For a cluster,
   `h(ell_z) <= (3/2) log R + log C`. **This is an upper bound only.** Genuine short-arc 3-clusters
   `((b-i)^2, b^2+1, (b+i)^2)` have Plücker gcd `2b` and `h(ell_z) = log(2(b^2 - 1)) = log R + O(1)`.
3. **Prop 6.6.** (UI-ell) is false for every `kappa < 1`, and (UI-ell_Sigma) is false for every `kappa_Sigma(n) < n/2`.
4. **Lemma 6.7, Prop 6.8, Cor 6.9.** Inputs acting through subgroup membership see only rank-deficient clusters, so
   their ceiling is `M <= omega(N') + O_C(1)`, the lit18 Prop B ceiling. **The Rémond-cone case is heuristic only.**

**Killed.** The one open residue, `(S_min)`, is **false for every `n`**.
* `n = 2`: on the proposal's own `(1 : g : g^2)` with `g = A + i`, `A -> inf`, the ratio is 1 > 2/3.
* `n = 2`, on actual circles (Prop K.1, proved): `z = ((b-i)^2, b^2+1, (b+i)^2)`, `b` even, has chord diameter
  `-> 4 sqrt R`, rank 1, and `Sigma_min / h(ell_z) -> 1`. *Judge re-check J2 (exact, independent code):* equal norms
  `R^2`, Gaussian gcd 1, `(diam/sqrt R)^2 = 15.84, 15.998, ... -> 16`, Plücker gcd exactly `2b`,
  `H(ell) = 2(b^2-1)`, and `exp(2 Sigma_min) = R^2` exactly, for `b = 10, ..., 10^20`.
* `n >= 3` (Prop K.2, proved): `u = (1, X_1 + i, ..., X_{n-1} + i, (X_1 + i)(X_2 + i))` with `X_x = T + x` has one
  multiplicative relation and ratio `-> n/2 > n/3`. The rescues fail too (exact evidence): rank `<= sqrt(2n)` gives
  ratio/`n` of 0.64-0.81, and a single full-support circuit gives ratio/`n` `-> 1`.

**Outcome.** Killed. The lens contributes a barrier-map row (§3). The structural reason it fails: the equal-norm
constraint is what makes clusters special, and restoring it returns the problem to Vojta / `(L_k)`.

## 2. Surviving live leads

**None.** A live lead here means a precisely stated input that is (i) not capped by any checklist item, (ii) not at
least as strong as the target or as Vojta, and (iii) has a plausible method. No item from this round meets all three.

| candidate left standing | (i) uncapped? | (ii) weaker than target / Vojta? | (iii) method? | verdict |
| --- | --- | --- | --- | --- |
| `(L_5^circ)` / `A_5(delta < 1/30)` (gauss) | **no**: item 3 caps it exactly (Lemma K1 + `construct.md` C3 on `k` rows, ratio 1), item 4 caps it exactly (single-form LW tight at `eps = 0`), and the `hardness.md` item applies (Liouville at a rational point) | it implies the target; predicted only by `V_E(5)` | none in either direction; `k <= 4` refuted | **not a lead**: same standing as `(L_8)` |
| "global coherence" / `(L_8)` (arakelov) | Vojta lane | **no**: implies the target | none | not a lead |
| non-robust composite items (a)-(d) (composite, §3 below) | unclassified | (b), (c) with multiplicativity, and (d) at super-polynomial precision are heuristically or provably "the problem itself" once used in full (Prop 3.1, Thm D shape) | none known to give growth | **map items, not leads** |
| `(S_min)` (wild) | n/a | n/a | n/a | **false** for every `n` |

So the last section of the round-17 candidate table ("composite-supported analytic inputs: ~0 for growth; a
classification would be map-only") is confirmed: the classification was map-only, and it is now mostly done.

## 3. Additions to the barrier map

Proposed edits to overview §5.4 / §5.5 / §7.6. All of them are **unrefereed** except where they restate a killer's
re-derivation; "proved" below means proved in the write-up and re-derived by its killer.

**New capped rows (§5.4).**

| # | input class | witness | `W` of witness | why capped | kind; status | source |
|---|---|---|---|---|---|---|
| 17 | **Robust composite-supported analytic inputs** (all-ideal Hecke large sieve over any support with constant `>= c(X + K X^{1/2})`; Type I/II, multilinear and Halász bounds with free coefficients; divisor-weighted sums; Selberg/GPY-weighted sums; polar-box counts with error `>= c X^gamma`), with or without multiplicativity in the forms P1-P5, precision fixed before `c_0` | minimal and multiplicative grafts of the `lw.md` 6.1 fakes (`A > max(4, 1/gamma)`) or of the P3-char fakes (`A > max(2, 1/gamma)`) | `(2A + o(1)) M log M`, resp. `(bA + o(1)) M log M` | grafts change only `Psi_S(X) <= X^{(1+o(1))/A}` ideals (minimal), or rotate each `S`-fibre by a constant (multiplicative). Every true all-ideal large sieve is robust (Prop 1.5) | F; proved (composite Thms M, P, Cor 2.2 as corrected by kill F1); R1 caveat, not vacuous here | `round20/composite.md`, `kill_composite.md` §8 |
| 17' | *(judge J1)* smoothed-norm, sharp-angle sector counts over **all ideals**, for any fixed finite family of sectors (rational, badly approximable or any edges), used **without** multiplicativity, at any precision `>= X^{(2+eps)/A - 1/2}`; and Lipschitz test functions of (norm, angle), including smooth Hecke sums, at precision `>= X^{(1+eps)/A - 1/4}` | swap extension of the row-17 witnesses (§4, J1) | as row 17, with `A > 4 + 2eps` (resp. `A > 4 + 4eps`) | the swap moves each fake ideal's genuine angle to a non-`S`-smooth ideal of nearby norm, so the multiset of (norm, angle) changes only by small norm shifts | F; **unrefereed** (judge; proof §4, exact-data check) | this README §4 |
| 18 | **Arithmetic intersection theory on genus-0 models** (arithmetic Hodge index, Faltings-Hriljac, adelic energy/capacity/Fekete, Bost slopes on `P^1`, equidistribution energies) | `round5/lower.md` witnesses (2-adic collisions included at `X = 2M`) | `(48 + o(1)) M log M` | in genus 0 every global relation is a pair relation and every finite configuration is ultrametric, so all of these are row 2 (ADC); Zhang/Northcott/bigness on a fixed model is row 8 (`tau* > 2` needed); Mordell-Weil packing is Lim (rank) | Red to rows 2, 8, Lim; proved (arakelov Prop B.2 repaired, Prop 3.3(b), Thm C) | `round20/arakelov.md`, `kill_arakelov.md` |
| 19 | **Unlikely intersections in tori** (BMZ, Habegger, Rémond cones, anomalous Zilber-Pink) | — (fixed carriers); full-rank fake profiles (subgroup case) | rows 1-3 | `V^oa` is empty for `X_M`, `L_M`, `M_{0,k}`; moving-line height comparisons are false at every useful constant, already on actual 3-clusters; subgroup membership sees only rank-deficient clusters, ceiling `omega(N')` | V + Red to lit18 Prop B; proved (wild Thm 6.2, Props 6.6, 6.8; kill_wild K.0-K.2); Rémond-cone ceiling heuristic | `round20/wild.md`, `kill_wild.md` |
| 20 | **gcd-reduced `k`-point circle statements** `(L_k^circ)`, `A_k(delta)`, `k >= 5`, attacked by forced divisibility, auxiliary forms, pair/cut/character data or single-form LW | — | — | killer Lemma K1: dividing by the tuple gcd leaves every such inequality literally unchanged; `construct.md` C3 on `k` rows gives ratio exactly 1; the target `e` is rational, so Roth/Schmidt do not apply | V (exact criticality); proved (kill_gauss §6.1, re-derived from the proposal's own §4.4) | `round20/gauss.md`, `kill_gauss.md` |

**Refuted (closed, not caps).** `A_4(delta)` for `delta < 1/14` and `(L_4^circ)`; `A_3(delta)` for `delta < 1/6` and
`(L_3^circ)` (gauss Thm 4.1 + killer corollary; proved). `(UI-ell)` for `kappa < 1`, `(UI-ell_Sigma)` for
`kappa_Sigma(n) < n/2`, and `(S_min)` with all its natural rescues (wild Prop 6.6; kill_wild K.0-K.3).

**§5.5 (not capped), replacing the line "Composite-supported analytic inputs: unclassified".** Only the non-robust
part stays, as four sub-items. None is known to give growth.
* (a) Sharp large sieve on primitive support, **without** multiplicativity = the composite max-scale cross-gap
  `|sin(<v,phi> - arg alpha)| >= c/max(Nd, N alpha)`. Open; it joins the existing sharp-large-sieve item of round 10.
* (b) The same **with** multiplicativity = the product-scale gap
  `|sin(psi_d + arg beta - arg alpha')| >= c/max(Nd Nbeta, N alpha')`. LCG-shaped; heuristically collapsing (only
  lattice directions pass); exact spacing evidence `5.6e-5` against `1.003` genuine. Not proved.
* (c) Smoothed sector counts with badly approximable edges, used **with** complete multiplicativity, or used
  uniformly over a whole bounded-type class of edges. *Judge J1 caps the finite-family, non-multiplicative case.*
  The uniform-in-type case is heuristically capped as well (§4, J1 remark).
* (d) Effective smooth Hecke sums **without** complete additivity, at precision finer than `X^{(1+eps)/A - 1/4}`.
  (Coarser precisions are capped by J1. With complete additivity, Prop 3.1 shows they are the problem itself.)

**§7.6 / round-17 §7.1 list: two new items.**
9. *No gcd/normalisation trick.* Dividing a sub-tuple by its Gaussian gcd changes no forced-divisibility, Liouville,
   pair, cut, character or single-form inequality (Lemma K1). A "subcritical exponent in the reduced radius" is not a
   gain unless the method sees the absolute size of `N''` (Vojta's `D = E` term does; nothing proved does).
10. *No genus-0 arithmetic intersection theory and no unlikely-intersection theorem on fixed carriers.* These reduce to
    pairs (row 2), `tau*` (row 8), or rank (`omega(N')` ceiling).

**Constraint audit.** No write-up assumed a uniform exceptional set from a fixed-`S` Subspace Theorem, discarded a
Gaussian gcd or denominator cost, set a drifting phase to zero, or treated a template-dependent threshold as uniform in
`R` (each killer checked this item by item). J1 likewise keeps every fake angle and every gcd; its threshold
`M_0(A, E, eps)` depends only on the fixed sector family and on `A`, which are fixed before the target constant.

## 4. Judge's own checks (unrefereed)

Scripts and outputs are in `judge_checks/`, run under `prlimit --as=3500000000`. The notes file is
`judge_checks/notes.md`. **None of this has been refereed.**

### J1. The swap extension: sector counts without multiplicativity are capped (proof + exact-data check)

*Setting.* Ideals of `Z[i]` are represented by Gaussian integers in the first quadrant, with genuine angle
`theta(alpha) ∈ [0, pi/2)`. A fake `(S, phi)` has `s <= b M` primes, all `>= P_0 = M^A`. `D` is the set of `S`-smooth
ideals `!= (1)`, i.e. the ideals all of whose prime factors lie above `S`; it contains the fake's points, divisors
and primes. The minimal graft gives `d ∈ D` the fake angle `psi_d = <v_d, phi>` (mod `pi/2`). For a finite family `F`
of closed arcs `I ⊂ [0, pi/2)` with edge set `E`, and `f ∈ C^1_c((1/2, 2))` with `0 <= f <= 1`, put
`S_I(X; theta') = sum_alpha f(N alpha / X) 1[theta'(alpha) ∈ I]`.

> **Proposition J1 (judge; unrefereed).** Fix `F`, `f`, `eps > 0` and `A > 4 + 2 eps`. Let `phi` avoid the
> measure-zero set where some `psi_d ∈ E`. Then for `M >= M_0(A, E, eps)` there is an injection
> `sigma : D -> Id \ D` (`Id` = the set of nonzero ideals) such that the angle function
> `theta_sw = psi_d on d ∈ D`, `theta_sw(sigma(d)) = theta(d)`, `theta_sw = theta` elsewhere
> satisfies, for every `I ∈ F` and every `X > 0`,
> `|S_I(X; theta_sw) - S_I(X; theta)| <= ||f'||_inf · sum_{d ∈ D, Nd ∈ [X/4, 4X]} |N sigma(d) - Nd| / X
> <= C_f X^{(2+eps)/A - 1/2}`,
> and the left side is `0` for `X < P_0/4`.

*Proof.*
1. *Pairing.* Suppose `theta(sigma(d))` and `psi_d` lie in the same cell of `[0, pi/2) \ E`. Then
   `1[psi_d ∈ I] = 1[theta(sigma d) ∈ I]` for every `I ∈ F`. The pair `{d, sigma(d)}` therefore contributes
   `(f(Nd/X) - f(N sigma d/X))(1[psi_d ∈ I] - 1[theta(d) ∈ I])` to the difference. Its absolute value is at most
   `||f'|| · |Nd - N sigma d| / X`. It vanishes unless `Nd` or `N sigma d` lies in `(X/2, 2X)`, and
   `|N sigma d - Nd| <= Nd/2` then forces `Nd ∈ [X/4, 4X]`. Ideals outside `D ∪ sigma(D)` do not change.
2. *Partners.* Process the `d ∈ D` in increasing norm. Let `Y = Nd >= P_0`, `H = Y^{1/2 + 1/A}`,
   `eta = Y^{-1/2 + 2/A}`. Let `J_d` be an arc of length `eta/2` with one endpoint at `psi_d`, lying in the cell of
   `psi_d`. This is possible once `eta` is below half the minimal gap of `E ∪ {0, pi/2}`, i.e. for `Y >= Y_0(E)`, which
   holds for `M` large because `Y >= M^A`. The polar box `B_d = {alpha : |N alpha - Y| <= H, theta(alpha) ∈ J_d}`
   has area `eta H / 2`. It is a Jordan domain with perimeter `L <= 2(H Y^{-1/2} + eta Y^{1/2})` (for `H <= Y/2`).
   Counting the unit squares that meet its boundary, it contains `>= eta H/2 - 4(L + 1) = Y^{3/A}/2 - O(Y^{2/A})`
   lattice points.
3. *Room.* The excluded points of `B_d` are the elements of `D` and the partners already used. Both have norm
   `<= 2Y`, so there are at most `2 Psi_S(2Y)`. By composite Lemma 1.2, `Psi_S(2Y) <= (2Y)^{log(4es)/log P_0}`, and
   `log(4es)/log P_0 <= (1 + eps)/A` for `M >= M_0`. Since `Y^{1/A} >= M`, we get
   `2 Psi_S(2Y) < Y^{3/A}/8` for `M >= M_0`. So a free partner `sigma(d) ∈ B_d` exists.
4. *Sum.* The number of `d` with `Nd ∈ [X/4, 4X]` is `<= Psi_S(4X)`, and each has `|N sigma d - Nd| <= H(4X)`. So the
   bound is `||f'|| Psi_S(4X)(4X)^{1/2 + 1/A}/X <= C_f X^{(2+eps)/A - 1/2}`. There are no `d` with `Nd < P_0`. `□`

*Consequences (unrefereed).*
* `theta_sw` agrees with the fake on every `S`-smooth ideal. So all local data (rows 1-4) are those of the fake, and
  Thm M's Hamming argument still applies, with `2 Psi_S` in place of `Psi_S`.
* Hence the composite class of row 17, together with any F3-type input for a fixed finite family of sectors (any
  edges), used **without** multiplicativity, is capped by the row-17 witnesses with `A > max(4 + 2eps, ...)`. This
  holds for every allowance at least the genuine error plus `X^{-delta}` for a fixed `delta > 0` (the usual
  "constants fixed before `c_0`" proviso).
* The killer's F3 example (precision `O(log X)`) is covered.
* *Lipschitz variant (sketch, same proof).* Take a test function `g(N alpha/X, theta(alpha))` with `|d_u g| <= L_1` and
  `|d_theta g| <= L_2`, for example `f(u) e^{4ik theta}`. Choose `eta = Y^{-1/4}`, `H = Y^{3/4}`, with no cell
  condition. Each pair changes the sum by `<= 2(L_1 H/X + L_2 eta)`, so the total change is
  `<= C(L_1 + L_2) X^{(1+eps)/A - 1/4}`. **So smooth Hecke sums used without complete additivity are capped at every
  precision `>= X^{(1+eps)/A - 1/4}`** (`A > 4 + 4 eps`). Item (d) of §3 is what remains.
* *Remark (heuristic).* For an input uniform over **all** edges of partial-quotient type `<= K_0`, the cell condition
  must hold for the closed nowhere-dense set `E_{K_0}`. Its Hausdorff dimension is `1 - O(1/K_0)`. A Borel-Cantelli
  count over the `Psi_S(X) ≈ X^{1/A}` values `psi_d` per scale, at distance `X^{-1/2+2/A}` from `E_{K_0}`, converges once
  `A >> K_0`. That suggests the uniform-in-type case is capped too, for `K_0` fixed before `A`. This is not proved: the
  `psi_d` are lattice combinations of `phi`, and the marginals of `<v, phi>` on the fake's realising set were not
  checked.

*Exact-data check* (`j1_swap.py` -> `j1_out.txt`). Same fake (seed 7, `s = 40` primes in `[1009, 1601]`), same 30
badly approximable sectors, same `f` as the killer's `c3b`. The genuine and minimal-graft numbers reproduce the killer's
exactly.

| `X` | # `S`-smooth in window | genuine mean / max `\|err\|` | minimal graft | swap extension | max `\|swap - genuine\|` | J1 bound (`\|\|f'\|\| sum \|dn\|/X`) |
|---|---|---|---|---|---|---|
| `1.5·10^6` | 3240 | 0.141 / 0.448 | 2.298 / 8.056 | 0.141 / 0.449 | 0.0014 | 0.144 |
| `2.5·10^6` | 2835 | 0.144 / 0.560 | 0.587 / 1.163 | 0.144 / 0.559 | 0.0015 | 0.089 |

There were no partner failures, with `eta = X^{-0.4}`. *Smooth Hecke sums* (`j1b_hecke.py` -> `j1b_out.txt`,
`X = 1.5·10^6`): `|F_k|` is `8e-12, 5e-12, 2e-12` genuine, `0.19, 1.08, 6.22` for the minimal graft, and
`0.016, 0.029, 0.039` for the swap extension (`k = 1, 2, 3`). The swap removes 1-2 orders of magnitude, as the
Lipschitz variant predicts, but it does not reach super-polynomial precision. This is consistent with item (d) staying
open.

### J2. Re-check of kill_wild Prop K.1 (exact integers; `j2_k1.py` -> `j2_out.txt`)

For `z = ((b-i)^2, b^2+1, (b+i)^2)` and `b = 10, 100, 10^4, 10^8, 10^20`:
* all three norms equal `R^2`, with `R = b^2 + 1`;
* the Gaussian gcd is 1;
* `(chord diameter / sqrt R)^2 = 15.8416, 15.9984, ... -> 16`;
* the Plücker gcd is exactly `2b`, and `H(ell_z) = 2(b^2 - 1)`;
* `exp(2 Sigma_min) = R^2` exactly.

So `Sigma_min / h(ell_z) -> 1`. Confirmed. (The killer's "arc length `-> 4 sqrt R` from above" concerns the arc; the
chord tends to `4 sqrt R` from below. Both are consistent.)

No other new mathematical claim is made in this README. Everything else is a killer-corrected statement or a reading of
the record.

## 5. What this means for growth and the uniform bound

**Growth: unchanged.** `M <= (2/3 + eps) log R/loglog R` is still the best bound valid for all `R`. No intermediate rate
`o(log R/loglog R)` follows from anything proved or from any standard conjecture. Every unconditional statement this
round is a cap, a refutation or a reformulation:
* robust composite inputs: capped;
* genus-0 arithmetic intersection theory: reduces to rows 2, 8 and Lim;
* unlikely intersections: vacuous, false, or capped at `omega(N')`;
* gcd-reduced `k`-point statements: exactly Liouville-critical.

**Uniform bound: unchanged, with one new conditional form.** The uniform theorem is still conditional only. The round
adds a fourth conditional route, `V_E(5)`: Vojta for `(Bl_e Lbar_5^tw, D = E, S = {inf})` on a **4-fold**. It is
non-vacuous (`(K + D_PQ).gamma_hat = 1/16 > 0`) and proved to imply the uniform theorem (gauss Thm 4.2 with Thm 3.4';
re-derived by the killer). Compared with the record's other instances:
* `V(7, 1/4)` on `Bl_Y X_7`: a 7-fold, `D = 0`;
* `V_8` / `(L_8)` on `Mbar_{0,8}`: a 5-fold, `D = 0`, `K.beta = 10`;
* `V_E(5)`: a 4-fold, but with `D = E` and a much thinner margin, `kappa = 2 delta ≈ 1/15`. The `D = 0` form needs
  `k = 6`.

This is a different instance, not a weaker hypothesis in any proved sense. It does not change the probability of an
unconditional result. It is, however, the smallest variety on record for an expert to look at. Its 4-point analogue
**fails** (Thm 4.1), which shows that the instance is sharp in `k`.

**What counts.** Only unconditional growth or a uniform bound counts as progress. Neither moved. Rubric (b) moved
nowhere: no write-up found an input that is true for actual circles, provably false for every fake of size
`O(M log M)`, and weaker than the problem. The record's diagnosis stands: such an input must be global, multi-place
and tied to the actual arithmetic of `N`. The only known statements of that kind are Vojta-type, and their sharpest
arithmetic forms are `(L_8)` and now `(L_5^circ)`, both open with no method in either direction.

## 6. Scores

Rubric of `round10/README.md` §1: (a) correctness after refereeing (here: after the kill report); (b) fake-violated
input; (c) movement of growth or uniformity. Each sub-score is out of 10. The overall score is a judgement that weighs
(b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **gauss** | 7 | 2 | 3 | **5** | The most substantial new mathematics of the round. It proves an elementary reduction `A_k => uniform`, an exact refutation at `k = 4` (extended to `k = 3` by the killer), a new non-vacuous Vojta instance on a 4-fold that gives the uniform theorem, and the circle analogue of round-18 item 7. All of it re-derived and recomputed. Deductions: the headline ("strictly below the critical exponent") is a change of units, as killer Lemma K1 shows, and the write-up did not see that the tightness is forced by invariance; `t_ab` integrality across unit classes; Thm 4.3(i) overclaimed; the density mechanism excluded only heuristically; "equivalently". (b) = 2: the needed statement is fake-violated but not known to be true and is Vojta-strength, not a weaker input. (c) = 3: a new conditional form of the uniform theorem plus a proved refutation and a proved cap (Lemma K1, with the killer). |
| **composite** | 6 | 2 | 3 | **4** | Correct, careful analysis (Lemmas 1.2-1.4, Prop 1.5, Thms M and P, Prop 3.1, all re-derived without error). It answers round-10 Slip 3 and moves robust composite inputs into the capped table. That is a proved correction to the barrier map, and it closes most of the last unclassified analytic row. Deductions: three major classification errors (F2: the multiplicative sharp sieve is product-scale; F3: a position-dependent input was missed; F4: "irrelevant" without proof), the mis-specified §5 requirement, and a repeat of a constant slip already corrected in round 10. (b) = 2: Prop 3.1 is fake-violated only because it is the problem itself. |
| **wild** | 6 | 0 | 2 | **3** | It chose a genuinely unused, `S`-uniform theory and settled it cleanly: fixed carriers vacuous, line comparisons false, subgroup ceiling `omega(N')`. All of this re-derived. Deductions: its one open residue `(S_min)` is false for every `n`, and the write-up wrongly said its own examples do not refute it; Prop 6.4(iii) was read as sharp, but it is only an upper bound on genuine clusters; Cor 6.9 overclaimed for cones. (c) = 2 for a proved barrier-map row. |
| **arakelov** | 6 | 0 | 1 | **2** | Theorem A is exact and pleasant, and the reduction of genus-0 arithmetic intersection theory to rows 2, 8 and Lim is correct. Deductions: its "needed statement" is at least as strong as the target ("equivalently `(L_8)`"); the proof of Prop B.2 omits the 2-adic place (repaired with a different witness and constant); scope gap `M <= X < 2M`; arithmetic slips; and the multi-point evidence adds nothing beyond Thm C. (c) = 1: it re-expresses existing rows. |

*Killer reports.* All four were decisive and careful. Each re-derived the proofs, recomputed numbers with independent
code, and found the real weaknesses: Lemma K1 (gauss), the 2-adic gap (arakelov), F2-F4 (composite), and `(S_min)`
false for every `n` (wild). Lemma K1 and kill_wild K.0-K.2 are proved new content and enter the barrier map (§3).

## 7. Single most promising next step, and whether to stop

**Stop internal rounds on growth and uniformity.** Round 19 recommended stopping. This round was an extra panel of four
untried lenses, and it confirms that recommendation from four independent directions:
* the Gaussian structure at small `k` reduces to a Liouville-tight, Vojta-strength statement;
* genus-0 arithmetic intersection theory reduces to pairs, `tau*` and rank;
* composite analytic inputs are capped in their robust part, and their non-robust part is the problem itself where it
  bites;
* unlikely intersections are vacuous or capped at `omega(N')`.

No lens produced an input that is fake-violated, weaker than the problem, and provable.

**The single most promising next step (external, not another internal round).** Add one short item to the expert package
of `round19/README.md` §6: **`(L_5^circ)` and `V_E(5)`**, stated as follows.

> *Question.* Let `T_5` be the norm-one torus of similarity classes of 5-tuples on circles. Are there `c > 0`,
> `eta_0 > 0` and a proper Zariski-closed `Z ⊂ T_5` such that no five primitive lattice points of a circle of radius
> `R`, with `eta_0`-balanced Gaussian profile over all 15 cuts and class outside `Z` (within one unit class), lie within
> `R^{7/15 + c}` of each other once `R` is large?

* A yes gives the uniform theorem (gauss Thm 3.4'; proved). Vojta for `(Bl_e Lbar_5^tw, D = E, S = {inf})` predicts
  yes (Thm 4.2).
* The analogue at 4 points is false (Thm 4.1): there is an explicit dense rational family at exponent `3/7`.
* Every method in the record is exactly critical on it (Lemma K1). This is a statement about approximation to a
  **rational** point of a 4-fold beyond Liouville, which is the McKinnon-Roth regime.

It sits next to `(L_8)` as the second self-contained Diophantine formulation. It is smaller (4-fold, 5 points), it
needs an archimedean divisor, and its margin is thinner. The expertise needed is the same: Vojta for blow-ups of toric
varieties at a rational point, and McKinnon's conjecture on best approximations.

*First concrete task, if anyone continues internally (not recommended):* a refutation attempt only. Search for a
2-parameter rational family of balanced 5-tuples at exponent `>= 7/15` in the shape of Thm 4.1's Gaussian linear
blocks, with block degrees chosen freely. The killer's expected-dimension count
`D(2 - (k-1)(1-e)) + k - 3 >= k - 2` already says that one-variable polynomial families cannot be dense at `k = 5`
for any profile, and Thm 4.3(ii) excludes bounded-degree covering families. So such a search is predicted to fail, and
success would refute `V_E(5)`. The `N <= 10^12` scan cannot test it, because balanced 5-tuples need at least 4 split
primes and the predicted count is `e^{-w}` per scale.

**Map-only follow-ups (optional; no growth either way).**
* Referee J1 (§4), and if it holds, record row 17' in overview §5.4.
* The two composite items that are plausibly closable: (b), a proof that the product-scale gap collapses (only lattice
  directions pass); and the uniform-in-type case of (c), via the Borel-Cantelli remark.
* Update overview §5.4/§5.5/§7.6 with §3 above.

**Honest probabilities.**
* Unconditional growth `M = o(log R/loglog R)` from any further internal round: **< 0.1%**.
* An unconditional uniform bound: **essentially 0**.
* A proof of `(L_5^circ)` or `(L_8)` in a further internal round: **< 0.5%** (each is Vojta-strength at a rational point,
  and Liouville-tight for every known method).
* A refutation of `(L_5^circ)` (which would refute `V_E(5)`): **~1%**. Such a refutation would not count as progress.
* That an external expert recognises `(L_8)` or `(L_5^circ)` as a known or tractable case of Vojta/McKinnon for blow-ups
  of toric varieties: unknown, and the only path with non-negligible value.

**Decision: stop.** The record's lasting value is unchanged: the `(2/3 + eps)` theorem, the conditional theorems (now
four forms, the smallest on a 4-fold), and the refereed barrier map. The map now has rows for composite analytic
inputs, genus-0 Arakelov theory and unlikely intersections, and a proved invariance cap on gcd normalisation.

## 8. Files

* Write-ups: `gauss.md` (`gauss_checks/`, notes `notes.md`), `arakelov.md` (`arakelov_checks/`), `composite.md`
  (`composite_checks/`), `wild.md` (`wild_checks/`, `wild_notes.md`).
* Kill reports: `kill_gauss.md` (`kill_gauss_checks/`, `kill_gauss_notes.md`), `kill_arakelov.md`
  (`kill_arakelov_checks/`, `kill_arakelov_notes.md`), `kill_composite.md` (`kill_composite_checks/`,
  `kill_composite_notes.md`), `kill_wild.md` (`kill_wild_checks/`).
* Judge: this `README.md`; `judge_checks/j1_swap.py` -> `j1_out.txt`; `judge_checks/j1b_hecke.py` -> `j1b_out.txt`;
  `judge_checks/j2_k1.py` -> `j2_out.txt`; `judge_checks/notes.md`.
* The repository was read only; the overview was not edited this round.
