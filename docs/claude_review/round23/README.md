# Round 23 — judge report (split calibration; tau_e width theorem)

| key | write-up | referee report | lens |
| --- | --- | --- | --- |
| split | `split.md` (`split_checks/`, `split_notes.md`) | `referee_split.md` (`referee_split_checks/`, `referee_split_notes.md`) | the split Losev-Manin 4-fold `LM_5` at `e`: is split Question B known, false, easier, or a dead end, and does anything transfer to `T_5`? |
| width | `width.md` (`width_checks/`, `width_notes.md`) | `referee_width.md` (`referee_width_checks/`, `referee_width_notes.md`) | the weighted mean-width inequality (W) and `tau_e(-K_{LM_5})` |

All paths are relative to `scratchpad/round23/`. The repository was read only. Labels are those of the overview:
proved / conditional / exact / evidence / heuristic / snippet / unrefereed. §1 gives the **referee-corrected**
statements. Remarks marked *judge* are mine and are **unrefereed**; my computations are in `judge_checks/` (§5).

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best bound for every
`R`. Neither write-up reaches `A_5(delta)` for any `delta < 1/30`, `alpha^ess_e(-K) > 15/4`, or `(L_5^circ)`, and
neither claims to. Both survive refereeing. Only growth or uniform progress counts, and there is none.

What changed (all refereed unless marked *judge*):
1. **`tau_e(-K_{LM_5}) = 10/3` is a theorem** (width Thm 2.5): `t_max(m) = floor(10m/3)` for every `m`, for the split and
   twisted forms. With toric Thm 3.3, toric Liouville (any finite family of auxiliary Laurent polynomials, arbitrary
   profiles) gives exactly Cilleruelo-Córdoba, `delta = 1/10`, and never better. This settles round 22's open point 1
   (map value). It **lowers** a section cap; it opens no lane. General (W) is still open.
2. **The split case is a calibration, not a source of mechanism** (split Thms 1, 2 and the parity lemma).
   * Split `A_5(delta < 1/30)` implies Erdős #887 in its absolute-`K` form, uniformly in the position on the
     hyperbola, and Ruzsa's #886 for `eps > (15/64)(1+2 delta)`. Both problems are open. So the split statement is
     **not known**, and it is at least as hard as an open Erdős problem.
   * The July 2026 five-divisor construction for #887 is a cut-balanced rational curve at exactly the Liouville constant
     `15/4`. It lies in a 2-parameter family that sweeps an irreducible hypersurface `W`, so every admissible split `Z`
     contains `W`.
   * The twisted form cannot carry such curves: Galois forces `E_A.C = E_{A^c}.C`.
   * Split Question B is **undecided** (not known, not refuted). It is not easier. No unconditional split mechanism
     exists to transfer.
3. *Judge (unrefereed, §5 J4).* The random-5-subset proof of split Thm 1 runs verbatim in the Gaussian setting:
   **`A_5(delta)` for one `delta < 1/30` implies an absolute Cilleruelo-Granville bound**, i.e. one `K` with
   `M_inf(C) <= K` for every `C`. Contrapositive: any construction with `M_inf(C) -> inf` as `C -> inf` would refute
   Question B. No such construction is known; Oganesyan's `R^alpha` disproof is withdrawn.

**Scores** (§6): split 4, width 4. **Decision** (§8): stop internal rounds. Fold §7 into the expert note and put
Question B, now with an Erdős-problem calibration, to experts.

## 1. Referee-corrected statements

### 1.1 split (`split.md`; `referee_split.md`): survives, status `partial`; one major error outside the main results

Setting: `T = G_m^5/G_m` over `Q`, `LM_5` its permutohedral compactification, `B = -K`, `e = [1:...:1]`. A class has
a unique primitive positive representative `u` near `e`. `H_{-K}([u]) = (max u/min u)·lcm(u)` (proved; round 22 J2
checked it exactly) and `dist([u],e) ≍ (max u - min u)/min u`.

* **Split Question B (precise).** Is `alpha^ess_e(-K_{LM_5}) > 15/4` on `T(Q)`? Equivalently, is split `A_5(delta)`
  true for some `delta < 1/30`? *Elementary form:* are there a proper Zariski-closed `Z ⊂ T` and `gamma > 15/4` such
  that only finitely many primitive 5-tuples `u` of distinct positive integers with `[u] ∉ Z` satisfy
  `((max u - min u)/min u)^gamma · (max u/min u) · lcm(u) <= 1`? Up to factors absorbed by the open condition: five
  factorisations of `N = lcm` within hyperbolic angle `N^{-(1+2 delta)/4}` on `xy = N`. The window `(15/4, 4)`, the
  section cap and the Vojta prediction `>= 4` are the same as in the twisted case.
* **Theorem 1 (proved; every step re-derived by the referee).** Let split `A_5(delta)` hold for one
  `delta ∈ (0, 1/30)`, with `Z`. Then there is `K = K(delta, Z)` such that for every `C > 0` there is `n_0(C)` with
  `#{d | n : X <= d <= X(1 + C n^{-1/4})} <= K` for all `n >= n_0(C)` and **all** `X > 0`.
  * This gives Erdős #887 `parts.ii` (absolute `K`), applied with `C + 1` to match the Lean integer window.
  * It also gives Ruzsa #886 for every `eps > (15/64)(1+2 delta)`, a range below `1/4`.
  * Canonical boundedness of every split `LM_k` would give #886 for every `eps > 0`. For `k`-subsets the averaging holds
    up to `+O(1/M)`, absorbed by the slack.
  * Mechanism: a uniform random 5-subset, profile averaging by thresholds (hypergeometric bound
    `>= 2^{-4}(1 - 15/(M-4))`), Northcott applied to **distinct-coordinate** tuples (referee's fix), and
    Schwartz-Zippel against one equation of `Z`.
  * Both problems are open: formal-conjectures `886.lean`/`887.lean` were fetched raw on 2026-10-28 by the referee and by
    me (J3). Hence the split statement is **not known**.
* **Theorem 2 (proved; the hypersurface claim upgraded to exact by the referee).**
  * *The family.* For `a, b, c ∈ Q` put `d = a+b+c` and `m = abcd`. Let `S, e_2, e_3` be the power sum and the
    elementary symmetric functions of `a², b², c², d²`, and `B = (e_3 - S m)/4`. On the conic `z² = m x² + B` put
    `F_0 = (x-a)(x-b)(x-c)(x+d)`, `F_1 = x⁴ - (S/2)x² + m + 2xz`, `F_2 = (x²-b²)(x²-c²)`, `F_3 = (x²-a²)(x²-c²)`,
    `F_4 = (x²-a²)(x²-b²)`. Then `C_{abc} = [F_0:...:F_4]` is a rational curve.
  * *Intersection profile.* `B.C = 15`, one constant place, and each of the 15 cuts is hit exactly once.
  * *Contact with `e`.* There are two conjugate real-quadratic branches of order 2 at `e`, so `mult_e = 4`.
  * *Liouville constant.* For non-square `m`, `alpha_C = 15/4` is attained: `dist^{15/4} H` stays bounded along the
    convergents (exact, up to `log H ≈ 1900`).
  * *The hypersurface `W`.* The family sweeps an irreducible hypersurface `W ⊂ T`. This is now exact: the map is
    invariant under `(λa, λb, λc, λ²t)` and has exact Jacobian rank 3 at 4 rational points. So **every admissible `Z` in
    split `A_5(delta < 1/30)` or `(L_5^split)` contains `W`**.
  * *The m-czech #887 example* (formal-conjectures PR #7218) is the member `(1,2,4)`: `z² = 56x² + 65`. Its offsets are
    `0, 2√56, 15, 18, 30` times `x²`. split.md printed `2√56 - 35` (typo).
* **Parity lemma (proved).** For a curve over `Q` on the twisted form, `E_A.C = E_{A^c}.C` for every `A` (Galois
  orbits; `v = -v` forces `v = 0`). So a cut-balanced twisted curve has `B.C >= 30`, and at `alpha = 15/4` it has
  `mult_e >= 8`, whatever the branch type. The weight-1 curves of Theorem 2 have no twisted counterpart.
* **Truth value: undecided.**
  * One-parameter subgroups give `alpha >= 8`.
  * Covering families give `alpha >= 4`. This uses the record's inequality; its two-point form was re-derived by neither
    the author nor the referee.
  * Polynomial blocks are overdetermined (**heuristic**, and split.md itself shows that the count can collapse under
    symmetry).
  * A refutation needs Zariski-dense A-points, e.g. infinitely many Pell-type hypersurfaces `W_i`. PR #7218's second
    five-divisor claim (`C = 7/2`) is unexamined and may be a `W'`.
  * A disproof of #887 (absolute `K`) would refute the split statement.
  * Vojta, via McKinnon-Satriano's erratum for toric varieties, predicts it true.
  * Exact data, `lcm <= 2·10^8`: min `gamma = 4.277` (cross-checked independently for `lcm < 2^20`). This is
    uninformative: W-points first reach `gamma < 4` at `log H ≈ 106`.
* **Killed.** §5.3 item 2, "group/coordinate products never improve `alpha`" (**major**): the inequalities given are
  upper bounds and give no lower bound on `alpha`. The claim is false: every point is a product of coordinate-1-PS
  points (`alpha >= 8`), including the `W`-points with `alpha = 15/4`. "Coordinate products do not refute it" is
  withdrawn.
* **Minor fixes.** The `F_1` offset; Northcott wording; `C+1` for the Lean window; the `+O(1/M)` in Remark (2);
  "`log H ≳ 60`" should read `≈ 106`; missing citations: Chan 1406.2230 (IJNT 2015) and 1303.2069, ErRo97 volume 79,
  and the `0.04742` constant (from a KSU abstract).

### 1.2 width (`width.md`; `referee_width.md`): survives in full; one false numerical remark (minor)

* **Theorem 2.5 (proved; referee found no gap).** Let `T_n(a,b) = {y ∈ Z^n : Σy = a-b, Σ y_+ <= a}`, the lattice points
  of `aΔ - bΔ`. For all `a, b >= 0` and `n = 2..5`, `reg T_n(a,b) <= F_n(a,b) = a + b + H_n(a,b)`, where `H_2 = 0`,
  `H_3 = min`, `H_4 = min(floor((a+b)/2), 2 min)` and `H_5 = min(min + floor(max/3), 3 min)`. Consequences:
  * `t_max(m) = floor(10m/3)` for **every** `m >= 1` on `LM_5` and on `Lbar_5^tw`. It is attained by `V^s V_123^r`, and
    for the twisted form by a `Q`-rational section via `sigma`-descent (toric §3; referee's clarity point). So
    **`tau_e(-K) = 10/3`**.
  * `t_max(m) = 3m` on `LM_4` and dP6.
  * (W) holds on every `T_5(a,b)`, with equality exactly on `S_{3s}` (Cor 3.1).
* **Mechanism.**
  * The slices of `T_n(a,b)` by a coordinate are exactly the sets `T_{n-1}(a',b)` and `T_{n-1}(a,b')`.
  * The peeling lemma (`construct.md` Thm B) turns this into a closed recursion.
  * The recursion closes by `min(X,Y) <= X/3 + 2Y/3 <= a + b/3`. That step uses the non-width bound `H_4 <= 2 min`.
* **Consequence (proved, with toric Thm 3.3).** `inf_F delta_F = 2/tau_e - 1/2 = 1/10` over finite families `F` of
  auxiliary Laurent polynomials over `Q(i)` with toric forced divisibility, with equality for `F = {V}`.
  * The cap is uniform in `F`. The witnessing profile may depend on `F`.
  * Toric Liouville never improves Cilleruelo-Córdoba for arbitrary profiles. On balanced profiles the cap stays `15/4`
    (`delta = 1/30`).
  * McKinnon-Roth base-locus Liouville at `e` gives exactly `alpha >= 10/3`.
* **Open.**
  * General (W), i.e. Conjecture C_bal at `M = 5`, reduced to pair-box polytopes. Evidence: 150 certified instances and
    60 pair-box peelings with no violation.
  * `tau_e(-K_{LM_6}) = 10/3`. Evidence: `rho_6(m,m) = floor(10m/3)` for `m <= 30`.
  * `reg T_5(a,b)` off the diagonal (`F_5` is not sharp at `(2,3)`).
* **Killed / corrected.**
  * §3.2: with the conjectural C_bal(4) slice bound, slicing and peeling gives **`7m/2`**, not `15m/4`. `15m/4` comes
    from the Theorem B(4) bound `(7/4)(2m-|c|)`. The conclusion that the mean-width template cannot reach `10m/3` stands,
    since `7m/2 > 10m/3`. I confirmed this independently (J1).
  * Minor: the twisted attainment needs descent; `(0,0)` was silently excluded from "14 cases".

## 2. The split case: known, false, easier, or dead end? Transfer to T_5?

| question | answer | basis |
| --- | --- | --- |
| **known?** | **No.** Split `A_5(delta < 1/30)` would settle Erdős #887 (absolute `K`, at every position on the hyperbola) and a new range of #886. Both are listed open as of 2026-10-28. | split Thm 1 (proved); formal-conjectures fetched raw by the referee and by me (J3) |
| **false?** | **Undecided.** No refutation is known. The split exceptional set must contain an explicit hypersurface `W`. Refutation would need infinitely many such `W_i` (Zariski-dense A-points), or a disproof of absolute #887. | split Thm 2; referee §5 |
| **easier?** | **No.** It has the same geometry and caps (`eps_e = 2`, `tau_e = 10/3`, balanced `15/4`, covering `4`, Vojta `4`), it is at least as hard as an open Erdős problem, and it carries **more** Liouville-extremal structure than the twisted form (the `W`-curves, excluded in `T_5` by parity). | split §7; parity lemma |
| **dead end?** | **Yes, as a source of mechanism; useful as a calibration.** Every unconditional split tool stops at `tau_e`. The only genuinely split tools are phase-0 integrality (Erdős-Rosenfeld `1 + C²`, Chan), which is quotiented out in `T` and whose use for a drifting phase would set the phase to zero, and fixed-template Pell/Schmidt arguments, whose thresholds depend on the template. Both are excluded by the owner's constraints. | split §6.1 |

**What transfers to `T_5` (map value only).**
1. *The parity lemma.* A cut-balanced twisted curve at `15/4` has `B.C >= 30` and `mult_e >= 8`. This explains, rather
   than merely records, why the round-22 twisted searches (`mu_0 = 8`, 30 roots) are so constrained.
2. *A cleaner reduction.* Random 5-subset plus Schwartz-Zippel, with no unit class and no balanced selection. *Judge
   (J4, unrefereed):* in the Gaussian setting this gives the **absolute** form `sup_C M_inf(C) < inf` from
   `A_5(delta < 1/30)`; see §5.
3. *A new twisted search ansatz.* Gaussian blocks on a Pell conic with two conjugate real-quadratic branches,
   `mu_0 = 4`, `B.C = 30`. Such curves cannot cover `T_5` (two-branch covering bound), so success would only enlarge
   `Z`. Map value only.

**Calibration consequence (logical remark, split §6.3).** A proof of `alpha^ess_e(-K) > 15/4` valid for every
`Q`-form of `LM_5` at once (as Vojta would be) solves #887 in its absolute form. So a proof of Question B is either
Erdős-#887-strong or uses the anisotropy of `T_5`. The parity lemma shows that the anisotropy is visible to the
Liouville-extremal curves. It does not show that any method can exploit it.

## 3. Status of (W) and tau_e

* **`tau_e(-K) = 10/3`: proved** (width Thm 2.5, refereed with no gap; J1 re-checks the recursion with my own code on
  `[0,70]^2`). It holds for split `LM_5` and twisted `Lbar_5^tw`. Round 22's interval `[10/3, 15/4]` collapses to its
  lower end. `t_max(m) = 3m` on `LM_4` and dP6 for all `m`.
* **(W) on `aΔ - bΔ`: proved** (Cor 3.1), with equality exactly on `S_{3s}`.
* **(W) in general (C_bal at `M = 5`): open.** It is reduced to pair-box polytopes, with 150 certified instances and 60
  pair-box peelings without violation (evidence). It is no longer needed for `tau_e`. Its only remaining value would
  be a single profile (uniform weight on the ten 2|3 cuts) witnessing the `1/10` cap for all `F` at once.
* **The barrier map changes as follows.** Row 20: best exponent for arbitrary profiles `delta = 1/10` (exactly
  Cilleruelo-Córdoba); for balanced profiles `delta = 1/30`. Both are proved. The "intermediate-exponent lane"
  (`delta ∈ (1/30, 1/10]`, arbitrary profiles) is **closed** for toric Liouville.
* **Effect on the target: none.** The window `alpha^ess ∈ (15/4, 4)` is untouched. The theorem moves every section
  argument further from it, for arbitrary profiles.

## 4. Growth and the uniform bound

**Nothing this round moves either.** Row by row:
* *Growth `M = o(log R/loglog R)`.* Neither write-up contains a growth statement. `tau_e = 10/3` is a cap on a method,
  and the split theorems are statements about the split torus. `(2/3 + eps)` stands.
* *Uniform bound.* The only routes remain those of round 22 §2: `(CB_e^norm1)`, or `alpha^ess > 15/4` on balanced
  one-unit-class sequences, then gauss Thm 3.4/3.4′. Both are open and Vojta-strength. This round adds two
  **difficulty calibrations**, not progress:
  * the split twin of the target implies an open Erdős problem (refereed);
  * the target itself implies the absolute form `sup_C M_inf(C) < inf` (judge J4, unrefereed). This is strictly
    stronger than the owner's goal `M(C) < inf` for each `C`.
* *Fakes.* No input in this round is fake-violated. Theorem 2.5 is about sections, which fakes satisfy. The split
  theorems are about the split torus, where there are no circles.

## 5. Judge checks (unrefereed)

Scripts and outputs are in `judge_checks/`, run under `prlimit --as=3500000000`. Notes are in `judge_checks/notes.md`.
**None of this has been refereed.**

**J1. width recursion, independent code** (`j1_recursion.py` → `j1_recursion.out`, exact integers).
* The exact recursive peeling numbers `rho_n` (with `rho_2 = a+b`) equal `F_n` on `[0,70]^2` for `n = 3,4,5`, with 0
  mismatches.
* `rho_5(m,m) = floor(10m/3)` and `rho_4(m,m) = 3m` for `m <= 70`.
* I re-derived the one non-obvious step by hand: `X <= 2(a-p) + q` and `Y <= (a+b-q)/2 + p` give
  `X/3 + 2Y/3 <= a + b/3`.
* Referee's §3.2 correction, confirmed by slicing `S_m` by `x_1` and peeling:

  | slice bound | result for `m = 2..48` |
  | --- | --- |
  | C_bal(4) bound `(3/2)(2m-\|c\|)` | `7m/2` exactly (`7, 10, 14, 21, 28, 42, 84, 168`) |
  | Theorem B(4) bound `(7/4)(2m-\|c\|)` | `15m/4` |
  | true slice values `F_4` | `floor(10m/3)` |

**J2. split Theorem 2 at new parameters** (`j2_pell.py` → `j2_pell.out`; sympy exact).
* The norm identity `F_1 F̄_1 = Π(x² - r²)` holds on the conic.
* The profile is: one constant place, 15 distinct cuts, `B.C = 15`, at `(1,2,4)` and at four parameter points not used
  by the author or the referee: `(2,5,9)`, `(1,4,6)`, `(3,4,10)` (where `B < 0`) and `(2,7,8)`. In all of them `m` is not
  a square.
* m-czech Pell points `x = 37, 1109, 33233, 995881, 29843197` (exact integers):
  * all five `F_j` divide `n` and lie above `√n`;
  * `C → 30.0000`;
  * the offsets divided by `x²` tend to `0, 14.9666 = 2√56, 15, 18, 30`, which confirms the referee's typo fix;
  * the tuple gcd alternates between 1 and 9.

**J3. Literature status.**
* I fetched `formal-conjectures/.../ErdosProblems/887.lean` raw myself (`judge_checks/dl/`, treated as data):
  * `parts.i` and `parts.ii` are tagged `research open`;
  * `rosenfeld_4` is `research solved` via `z² = 56x² + 65`;
  * the `C = 7/2` claim is not in the file. The PR body could not be fetched here (GitHub REST is not enabled for this
    repository in this session).
* One WebSearch (2026-10-28):
  * Cilleruelo-Granville's Conjecture 14 concerns arcs of length `R^{1-eps}`;
  * at the `√R` scale the known individual bound is `≪ log R`;
  * Oganesyan's claimed `R^alpha` disproof (arXiv 2107.09991) is **withdrawn**, as the record already notes
    (lit18 `primes.md` #26).
  * No source settles whether `sup_C M_inf(C) < inf`. This is search-negative.

**J4. Absolute form of the reduction (judge's proposition; proof below; exact identity check).**

> **Proposition J4 (unrefereed).** Suppose `A_5(delta)` (expert note §3) holds for one `delta ∈ (0, 1/30)`, with
> `Z ⊂ T_5` and `R_0`. Then there is `K = K(delta, Z)` such that for every `C > 0` there is `R_1(C)` with: every arc
> of length `C√R` of every circle `x² + y² = R²` with `R >= R_1(C)` contains at most `K` lattice points. Hence
> `M_inf(C) <= K` for all `C`, and `M(C) < inf` for all `C`.

*Proof.*
1. *Setup.* Let `ζ_1, ..., ζ_M` (`M >= 5`) be the distinct lattice points on the arc. Their mutual distances are at most
   `C√R`. Write `N = R²` and factor
   `ζ_j = ε_j (1+i)^{e_0} Π_q q^{f_q} Π_p π_p^{a_j(p)} π̄_p^{e_p - a_j(p)}`. Here the ramified and inert parts are common
   to all points, because `|ζ_j|² = N`.
2. *Exact reduced radius.* For a 5-subset `S`, the Gaussian gcd is
   `g_S ~ (1+i)^{e_0} Π q^{f_q} Π_p π_p^{min_S a} π̄_p^{e_p - max_S a}`. The reduced tuple `ζ_S/g_S` is primitive, has
   distinct points, and lies on radius `R''_S` with
   `R''_S² = Π_p p^{max_S a(p) - min_S a(p)}`. This is an exact identity, checked on all 252 five-subsets of the
   record's 10-point cluster in `j4_gauss_avg.out`.
3. *Averaging.* `max_S a - min_S a = Σ_{s=1}^{e_p} 1[S meets {a >= s} and {a < s}]`. By the refereed hypergeometric
   bound of split Thm 1(ii), `E log R''_S² <= (15/16 + 15/(16(M-4))) log N` for `M >= 9`. Fix `η > 0` with
   `(15/16 + η)(1 + 2 delta) < 1`; this is possible since `delta < 1/30`. If `M >= 9 + 2/η`, at least
   `(η/2) C(M,5)` subsets have `R''_S <= R^{15/16 + η}`; this is step (iii) of split Thm 1, verbatim.
4. *The good subsets are A-points.* For such `S`, the reduced distances are at most
   `C√R/|g_S| = C R''_S R^{-1/2} <= R''_S^{1/2 - delta}`. That inequality is equivalent to
   `C² R''_S^{1+2 delta} <= R`, and it holds once `R >= R_1(C)`.
5. *The reduced radius is large.* There are finitely many primitive 5-tuples of distinct lattice points with
   `R'' < R_0`. Their classes are `≠ e`, so `dist(·, e) >= c_0 > 0` on them. But
   `dist([ζ_S], e) = max|ζ_j - ζ_1|/R <= C R^{-1/2} < c_0` for `R >= R_1(C)`. So `R''_S >= R_0`, and `A_5(delta)` gives
   `[ζ_S] ∈ Z`.
6. *Schwartz-Zippel.* Over `C`, `T_5 ≅ (C^*)^5/C^*` via `ζ ↦ (ζ_j)`, since the norm-one torus becomes `G_m` through
   `u = x + iy`. So some nonzero homogeneous `F` of degree `D` vanishes on the cone over `Z_C`. Each good `S`, listed
   in angular order, is a zero of `F` in `{ζ_1, ..., ζ_M}^5`; `F` is homogeneous, so dividing by `g_S` does not matter.
   Hence `(η/2) C(M,5) <= D M^4`, which forces `M <= M_1(D, η)`.
7. *Conclusion.* `K = max(9 + 2/η, M_1)` depends only on `delta` and `Z`. `C` enters only through `R_1(C)`, and the
   position on the circle does not enter at all. `□`

*Status and meaning.*
* This is the Gaussian twin of split Thm 1. split.md §6.2 item 1 suggested it as an unchecked assessment. It needs no
  unit class and no balanced selection.
* The record already has `A_5(delta < 1/30) ⇒` the uniform bound (gauss Thm 3.4′, via `(L_5^circ)`). What is new is
  **only the `C`-independence of `K`**.
* So Question B is at least as strong as an *absolute* Cilleruelo-Granville statement, `sup_C M_inf(C) < inf`. This is
  the circle analogue of Erdős #887 `parts.ii`.
* *Contrapositive:* a construction with `M_inf(C) -> inf` as `C -> inf` would refute `A_5(delta)` for every
  `delta < 1/30`, and with it `V_E(5)`.
* The record's families give `M_inf(C) >= 8` for `C > 5.89771`, and nothing unbounded. I make no claim that the
  absolute form is true or false.
* *Constraint check:* the only exceptional set is the hypothesis `Z`; `g_S` is divided out exactly and `R''` is exact;
  no phase occurs; `R_1(C)` depends on `C` and is not claimed uniform.

## 6. Scores

Rubric of `round10/README.md` §1. Each sub-score is out of 10:
* (a) correctness after refereeing;
* (b) fake-violated input;
* (c) movement of growth or uniformity.

The overall score is a judgement that weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **split** | 7 | 0 | 2 | **4** | **Strengths.** Theorem 1 is a correct, clean reduction tying split Question B to the open Erdős #887 (absolute `K`, phase-uniform) and #886. That answers "is it known?" with a proof, not only a search. Theorem 2 recognises the 2026 #887 construction as a balanced `15/4` curve, generalises it to a family, and proves (with the referee's exact upgrade) that every split `Z` contains a hypersurface. The parity lemma is a real structural difference between the split and twisted forms. Literature was handled honestly at snippet level. **Deductions.** One major false claim ("products never improve `alpha`"), which contradicted the write-up's own Theorem 2. The Jacobian rank was labelled as feeding a proved claim while only numerical. Further minor items: an offset typo, the Northcott wording, an unreproduced "`log H ≳ 60`", missing citations, and the unexamined `C = 7/2` family. **(b) = 0:** the split torus has no fakes in play. **(c) = 2:** a proved calibration and a barrier-map item (split ≥ #887; parity), with no movement of the target. |
| **width** | 9 | 0 | 2 | **4** | **Strengths.** A complete, elementary proof of `t_max(m) = floor(10m/3)` for all `m` (also `3m` on `LM_4` and dP6). It settles round 22's open point 1 and turns the round-22 overclaim (kill_toric F1) into a theorem: toric Liouville is capped at exactly Cilleruelo-Córdoba for arbitrary profiles. The quantifiers are careful (the witnessing profile is not claimed uniform in `F`). Labels are clean: computations are certified only in the valid direction. **Deductions.** One false numerical remark (`15m/4` should be `7m/2`) and a descent step left implicit. **(b) = 0:** sections are fake-satisfied. **(c) = 2:** a proved correction to row 20 that closes the intermediate lane. By its own statement it moves the target away from every section argument, not toward it. |

*Referee reports.* Both were thorough and independent.
* referee_split re-derived Theorem 1 and Theorem 2 with fresh code. It upgraded the Jacobian claim to exact, caught the
  one major logical error, and found the omitted `C = 7/2` claim.
* referee_width re-derived every lemma, cross-checked with a different prime and exact sympy, and caught the `7m/2`
  slip.

For comparison: round 22's toric lens scored 4 and lit 3, and round 20's gauss lens 5. This round does not move
(b) or (c) beyond the barrier map.

## 7. Proposed additions to expert_note.md §3

These are proposed text only; the repository was not edited. Items 1-2 are refereed. Item 3 is the judge's and must
be marked unrefereed if adopted. It should be adopted only after a referee pass, or with its label.

**(1) Replace the `tau_e` bullet in "In the approximation-constant literature's own terms".** Replace
"$10/3\le\tau_e(-K)\le15/4$ (**proved**; exact data … suggest $\tau_e=10/3$)" by the following paragraph.

> $\tau_e(-K)=10/3$ (**proved**, round 23 `width.md` Thm 2.5, refereed). For every $m$, the maximal order at $e$ of a
> section of $m(-K)$ is $\lfloor10m/3\rfloor$, for the split and the twisted form, and it is attained by Vandermonde
> products. The proof slices the lattice points of $a\Delta-b\Delta$ coordinatewise, which gives the same family again,
> and peels the slices. Consequently every finite family of auxiliary Laurent polynomials over $\mathbb Q(i)$, used
> through Liouville with toric forced divisibility, gives for arbitrary profiles exactly the Cilleruelo–Córdoba exponent
> $\delta=1/10$ and never more; for balanced profiles it gives $\delta=1/30$. McKinnon–Roth's base-locus Liouville at $e$
> is therefore exactly Cilleruelo–Córdoba. The open window $\alpha^{\rm ess}\in(15/4,4)$ is untouched.

**(2) New paragraph after "Question B, literature form": the split calibration.**

> **Split calibration** (round 23 `split.md`, refereed). The same window occurs on the *split* Losev–Manin 4-fold at
> $e$. There a primitive positive 5-tuple $u$ has $H_{-K}=(\max u/\min u)\,\mathrm{lcm}(u)$ and
> $\mathrm{dist}\asymp(\max u-\min u)/\min u$. Split Question B asks whether five distinct integers with gcd 1 in
> $[U,U(1+\mathrm{lcm}^{-(1+2\delta)/4})]$, $\delta<1/30$, have their ratio class in a fixed proper closed set once the
> lcm is large. Equivalently, it concerns five factorisations of $N$ within hyperbolic angle $N^{-(1+2\delta)/4}$ on
> $xy=N$. It is **not known**: it implies Erdős Problem #887 in its absolute form ($\exists K\,\forall C$: eventually at
> most $K$ divisors of $n$ in $[X,X(1+Cn^{-1/4})]$, for every $X$), and Ruzsa's #886 for
> $\varepsilon>(15/64)(1+2\delta)$. Both are open (formal-conjectures, October 2026). The proof uses a random 5-subset,
> averaging over threshold layers, and Schwartz–Zippel. The 2026 five-divisor construction for #887
> ($n=(x^2-1)(x^2-4)(x^2-16)(x^2-49)$, $z^2=56x^2+65$) is a rational curve hitting each of the 15 cuts once, with
> $\alpha=15/4$ exactly. It lies in a 2-parameter family sweeping a hypersurface $W$, which every admissible split
> exceptional set must contain. Such curves do not exist on the twisted form: for a curve over $\mathbb Q$ there,
> $E_A\cdot C=E_{A^c}\cdot C$, so a balanced twisted curve at $15/4$ has $B\cdot C\ge30$ and $\mathrm{mult}_e\ge8$. A
> proof of Question B valid for every $\mathbb Q$-form of the 4-fold would therefore solve #887; a proof for $T_5$ alone
> may need the anisotropy. Split Question B is undecided. It would fail if infinitely many such hypersurfaces existed,
> or if absolute #887 failed.

**(3) Optional, judge J4 (unrefereed). Add to "Theorems" under "Reduction".**

> *Absolute form (round 23 judge, **unrefereed**).* The same random-5-subset argument gives: if $A_5(\delta)$ holds for
> one $\delta<1/30$, there is $K=K(\delta,Z)$ such that for every $C$, every arc of length $C\sqrt R$ on every circle of
> radius $R\ge R_1(C)$ contains at most $K$ lattice points, i.e. $\sup_CM_\infty(C)<\infty$. The proof uses exact
> reduced radii $R''^2=\prod_pp^{\max_Sa(p)-\min_Sa(p)}$, needs no unit class and no balanced selection, and
> Schwartz–Zippel replaces the counting lemma. So any construction with $M_\infty(C)\to\infty$ as $C\to\infty$ would
> refute Question B. The record's families give $M_\infty(C)\ge8$ and nothing unbounded; Oganesyan's $R^\alpha$
> disproof is withdrawn.

**Also (barrier map row 20 and expert note §5 item 3).** Replace "`tau_e ∈ [10/3, 15/4]`" by
"`tau_e = 10/3` (round 23)". Read the best exponent as "`delta = 1/10` for arbitrary profiles, `1/30` balanced". In
round-22 README §2, the row "toric sections" becomes `delta > 1/10` for arbitrary profiles. In §3's "Question B"
paragraph, add the twisted search ansatz as a *suggestion*: two conjugate real-quadratic branches with `mu_0 = 4`,
`B.C = 30`.

## 8. Next step, probabilities, stop or continue

**Decision: stop internal rounds on growth and uniformity.** This confirms rounds 19, 20 and 22. This round tested the
last two internal items that round 22 listed: the split calibration and (W)/`tau_e`.
* The first shows that the split twin of Question B is Erdős-#887-hard and gives no mechanism.
* The second settles a cap whose value cannot affect the target.

Every unconditional tool still stops at `tau_e`, which is now exactly `10/3` (`15/4` balanced). The record has no
non-section mechanism for the window `(15/4, 4)`.

**Single most promising next step (external).** Fold §7 (1)-(2) into `expert_note.md` §3, and (3) after a referee
pass. Then put Question B to two audiences, each in its own language:
1. *Approximation of rational points* (McKinnon-Roth, McKinnon-Satriano, Huang). Is the twisted toric Fano `X'` of
   the `A_4` root polytope canonically bounded at `e` on `T_5(Q)`, or at least does it have `alpha^ess_e(-K) > 15/4`?
   Add the new facts: `tau_e = 10/3` exactly, and the split form is Erdős-hard and carries the hypersurface `W`. So any
   proof for `T_5` alone must use the anisotropy, or must solve #887.
2. *Erdős problems / close factorisations* (the #887 thread; Chan, Cilleruelo-Jiménez-Urroz, Granville-Jiménez-Urroz
   lineage). Ask the **W-question**: are the Liouville-extremal five-divisor patterns (`alpha = 15/4`) contained in
   finitely many Pell-type families? Is PR #7218's `C = 7/2` family another such hypersurface? A "finitely many"
   answer is the split analogue of what Question B needs; an "infinitely many" answer refutes split `A_5`.

*Cheap internal item, if any:* a referee pass on J4 (one short proof, reusing refereed lemmas). Value: it would let
the expert note state the absolute form as proved.

**Honest probabilities.**

| outcome | probability |
| --- | --- |
| unconditional growth `M = o(log R/loglog R)` from a further internal round | **< 0.1%** |
| an unconditional uniform bound from a further internal round | **essentially 0** |
| a proof of `(L_5^circ)` or `A_5(delta < 1/30)` internally | **< 0.5%** (unchanged; the split calibration shows that an all-forms proof would also solve #887) |
| the literature already contains a theorem implying `(L_5^circ)` | **about 1%** (lowered from 1-2%: such a theorem, if form-independent, would also settle #887, which is listed open) |
| a refutation of Question B, e.g. Zariski-dense twisted A-points, or `M_inf(C) -> inf` via J4 | **about 1-2%** (it would not count as progress) |
| J4 survives an independent referee | **about 85-90%** (every ingredient is a refereed lemma; the new parts are the exact reduced-radius identity and the Schwartz-Zippel coordinates) |
| a search on the two-branch twisted ansatz (`mu_0 = 4`, `B.C = 30`) finds a balanced twisted `15/4` curve | **about 10-20%** (map value only; it cannot cover `T_5`) |
| split Question B is false through infinitely many `W_i` | **about 5-10%** (author's judgement; Vojta says no) |
| an expert recognises Question B or the W-question as known or tractable | small (a few percent); still the only path with non-negligible value |

**Stop.** The record's lasting value:
* the `(2/3 + eps)` theorem;
* the conditional forms of the uniform theorem, the smallest being `alpha^ess_e(-K) > 15/4` on a 4-fold;
* the refereed barrier map, with row 20 now exact (`1/10` arbitrary, `1/30` balanced);
* two calibrations of the target's strength: split ⇒ Erdős #887 (refereed), and twisted ⇒ absolute Cilleruelo-Granville
  (judge, unrefereed).

## 9. Constraint audit

The owner's four constraints, checked for every claim this README relies on:
* **No uniform exceptional set from a fixed-`S` Subspace theorem.** The only exceptional sets are:
  * the hypothetical `Z` of `A_5`/`A_5^split`;
  * the explicit `W`;
  * the zero sets of auxiliary sections (toric Thm 3.3).

  Chan/Temur-type Pell-Schmidt arguments are named and not used.
* **No Gaussian gcd or denominator cost discarded.**
  * split Thm 1 and J4 divide the tuple gcd `g_S` exactly and use the exact toric height (in J4, the exact
    `R''² = Π p^{max-min}`, checked exactly).
  * Theorem 2's rational points have their denominators cleared exactly (author and referee).
  * width Thm 2.5 is pure interpolation, and its consequence keeps `N'^m` (toric Thm 3.2).
* **No drifting common phase replaced by zero.**
  * The phase is quotiented out in `T` and `T_5`.
  * split Thm 1's `n_0(C)` is independent of the position `X` (referee-verified).
  * J4's `R_1(C)` does not depend on the position on the circle.
  * Phase-0 integrality (Erdős-Rosenfeld `1 + C²`) is identified as such and **not** used.
* **No template-dependent threshold treated as uniform in the radius.**
  * `K(delta, Z)` depends only on the hypothesis.
  * `n_0(C)` and `R_1(C)` are explicitly `C`-dependent and are not claimed uniform.
  * The `1/10` cap is an infimum over all finite families, and its witnessing profile is **not** claimed uniform in `F`.
  * The Pell family's "`gamma < 4` from `log H ≈ 106`" is a property of one family, used only to explain why the scan is
    uninformative.

## 10. Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round23/`.
* Write-ups:
  * `split.md` (`split_checks/`, `split_notes.md`);
  * `width.md` (`width_checks/`, `width_notes.md`).
* Referee reports:
  * `referee_split.md` (`referee_split_checks/`, `referee_split_notes.md`);
  * `referee_width.md` (`referee_width_checks/`, `referee_width_notes.md`).
* Judge: this `README.md`, plus `judge_checks/`:
  * `notes.md`;
  * `j1_recursion.py` → `j1_recursion.out` (J1);
  * `j2_pell.py` → `j2_pell.out` (J2);
  * `dl/887.lean` (J3, downloaded data);
  * `j4_gauss_avg.py` → `j4_gauss_avg.out` (J4 identity check).
* The repository, the overview and the expert note were not edited.
