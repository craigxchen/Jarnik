# Round 22 judge: `(L_5^circ)` through toric approximation constants (two lenses, two kill reports)

> **Correction notice (2026-10-29, full-text check).** `../papers_check/README.md` §2.2 lists 12 corrections to this
> README after reading seven of the cited papers in full. In particular, the table row grouping McKinnon–Satriano,
> Monahan–Satriano and Lehmann–McKinnon–Satriano as "conditional on Vojta; split" is wrong: the erratum's Prop 0.1 and
> LMS Prop 1.4 apply to every smooth projective $X$ with $-K$ big. The bottom line is unchanged.

| key | write-up | kill report | lens |
| --- | --- | --- | --- |
| lit | `lit.md` (`lit_checks/`, `lit_notes.md`) | `kill_lit.md` (`kill_lit_checks/`, `kill_lit_notes.md`) | literature on approximating rational points (McKinnon, McKinnon-Roth, Huang, McKinnon-Satriano, Batyrev-Tschinkel), re-examining the lit18 "drifting phase" dismissal |
| toric | `toric.md` (`toric_checks/`, `notes.md`) | `kill_toric.md` (`kill_toric_checks/`, `kill_toric_notes.md`) | Seshadri constant, Liouville threshold and approximation constants of `(Lbar_5^tw, -K, e)` |

All paths are relative to this directory (`scratchpad/round22/`). The repository was read only. Labels are those of the
overview: proved / conditional / exact / evidence / heuristic / unrefereed. **Nothing written this round has been refereed by
anyone other than the killers.** §1 gives the killer-corrected statements. Remarks marked *judge* are the judge's own and
are **unrefereed**; the judge's computations are in `judge_checks/` (§3).

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is still the best bound for every `R`.
**No known theorem and no toric argument reaches `A_5(delta)` for any `delta < 1/30`, or `(L_5^circ)`.** Both lenses
survive their killers with corrections. The round leaves the target where round 20 left it, and changes how it is stated.

What changed:
1. **Literature form of Question B.** In McKinnon-Roth normalisation, with `H = H_{-K}([zeta]) = R''^2` and
   `dist ≍ diam/R''`, an A-point at `delta` has exponent `gamma = 4/(1+2 delta)`. `A_5(delta < 1/30)` is then exactly
   **`alpha^ess_e(-K) > 15/4` off a proper closed set**, and `(L_5^circ)` is the same over balanced sequences. The values
   are:
   * `2` is the best constant `alpha_e = eps_e`; it is attained only on the diagonal conics, which carry repeated points;
   * `10/3` is Cilleruelo-Córdoba;
   * `15/4` is the balanced Liouville constant and an upper bound for every section constant `tau_e`;
   * `4 = dim X` is the counting and Vojta value.

   The open window is `alpha^ess ∈ (15/4, 4)`. In the words of McKinnon-Satriano (TAMS 2021), the strongest form is
   *"`X` is canonically bounded at `e`"*, restricted to norm-one points. That form implies `A_5(delta)` for every
   `delta > 0`, hence the uniform bound. It is a **hypothesis** in their theorems, not a result (killer F4; judge search §3, J3).
2. **The lit18 dismissal is replaced.** In `T_5` the common phase is quotiented out, not set to zero, and the target `e`
   is fixed and rational. Huang, McKinnon-Roth and McKinnon-Satriano are therefore *applicable*. They are simply *too weak*:
   * the best constant they bound is `2`, computed inside `Z_diag`;
   * section and base-locus (Liouville) bounds are capped at `tau_e <= 15/4`;
   * Roth-type bounds satisfy `beta <= tau`;
   * counting works at the natural scale;
   * the strong statements are split, or conditional on Vojta.
3. **The section cap, now for all profiles.** Every finite family of auxiliary Laurent polynomials used through Liouville
   with toric forced divisibility gives exactly `A_5(delta)` for `delta > 2/tau_e - 1/2` (proved, by minimax). Here
   `10/3 <= tau_e(-K) <= 15/4` (proved). Exact data: `t_max(m) = floor(10m/3)` for `m <= 8`, so `tau_e = 10/3` is
   likely (evidence).
4. **Geometry fixed.** `Lbar_5^tw` is only weak Fano. `-K` is nef and big but not ample, and equals `phi^*(-K_{X'})` for
   the small crepant contraction `phi` onto the twisted Gorenstein toric Fano 4-fold `X'` of the `A_4` root polytope.
   `eps_e(-K) = 2` and `(-K)^4 = 70`.

**Scores** (§5): toric 4, lit 3. **Recommendation** (§6): stop internal rounds. Fold the literature form into the expert
note (§4).

## 1. Killer-corrected outcomes

### 1.1 lit (`lit.md`; `kill_lit.md`): the verdict `no_new_result` survives; one major error, repaired

**Survives (killer-corrected).**
1. **Search-negative verdict (labelled, not a proof).** Neither lit.md nor the killer's independent searches found a
   theorem in the literature that implies `A_5(delta)` for some `delta < 1/30`, or `(L_5^circ)`. Neither found an
   unconditional theorem, for any variety, giving `alpha^ess > tau_x` at a rational point. No full text could be read
   (arXiv, ar5iv, uwaterloo, queensu, springer and smf do not resolve in this session). Every statement about the
   *mechanisms* of Huang's and McKinnon-Roth's proofs is therefore **assessment**.
2. **Dictionary (derived; killer-checked).**
   * Work on `X'`, the Gorenstein toric Fano 4-fold of the root polytope, with `L = -K_{X'}`; equivalently use the nef
     and big `-K_X = phi^*(-K_{X'})` on `X = Lbar_5^tw`.
   * Then `X(Q) = X'(Q) = T_5(Q)`, `h_{-K}([zeta]) = log R''^2` exactly, and `dist ≍ diam/R''`.
   * An A-point at `delta` has `gamma = 4/(1+2 delta)`: `15/4` at `delta = 1/30`, `10/3` at `1/10`, `4` at `0`.
   * The target is `alpha^ess_e(-K) > 15/4` on `eta_0`-balanced points of one unit class, off `Z ⊂ T_5`.
   * Vojta `V_E(5)` gives `alpha^ess >= 4 - eps`. This needs only `-K` big, not ample.
3. **Prop 1.3 (correct; a restatement of row 20, `construct.md` Thm B / gauss §4.4).**
   * Every fixed-section, forced-divisibility (Liouville, Seshadri, effective-threshold) bound at `e` gives at most
     exactly `15/4` on balanced points. This holds for **every Laurent polynomial over `Q(i)`**, not only Galois-invariant
     toric `L`, because widths are translation- and sign-invariant.
   * New here: `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}`.
   * Per-point auxiliary sections and multi-point Thue-Roth arguments are covered only through `beta <= tau`. That is
     correct for Ru-Vojta / Evertse-Ferretti with the single irreducible `D = E` and `S = {inf}`.
4. **Exact ranks (reproduced twice, independently).** Maximal `ord_e` on `H^0(m(-K))` is `3, 6, 10, 13, 16` for
   `m = 1..5`. The Vandermonde attains `10` at `m = 3`. So McKinnon-Roth's base-locus Liouville at `e`, at these levels,
   gives exactly Cilleruelo-Córdoba (`10/3`). Also `(-K)^4 = 70`.
5. **Literature verdicts.**
   * McKinnon's conjecture (Coba and its weak form) and all its proved cases concern `alpha_x`, not `alpha^ess`. They are
     split or homogeneous. The conditional cases rest on Vojta for *those* varieties.
   * Roth-type bounds (McKinnon-Roth 2015, Ru-Vojta, Grieve, Heier-Levin) are `<= tau` at a rational target.
   * Counting and equidistribution work at the natural scale with relative errors: Huang (split), Huang-de Saxcé, and
     Batyrev-Tschinkel (which *does* apply to anisotropic `T_5`, with `rank Pic^Gal = 15` and `~ c B (log B)^{14}`).
     They cannot see a sparse Zariski-dense set.
6. **lit18 re-examination: correct.** The drifting-phase objection is void in `T_5`. The theorems apply and are too weak.
7. **Split calibration (correct; proof supplied by the killer for dP6).**
   * On `P^1 x P^1` at `(1,1)`, `alpha^ess = tau = 4`.
   * On split dP6 `= Lbar_3` at `e`, `alpha^ess = tau = 3`: the family `x = (P_1P_2, P_2P_3, P_1P_3)`, `P_i = n + c_i`,
     gives `H ≍ n^3` and `dist ≍ 1/n`.

   At split rational points Liouville is the truth, and the points are *repelled*: `alpha^ess > dim`.

**Killed / corrected.**
* **F1 (major, false statement).** The claim that `X` is a toric Fano 4-fold with `-K_X` ample is false for every
  `k >= 4`. On the chamber of `s`, `m_s = e_{s_1} - e_{s_k}` is not strictly convex: there are 480 violations at
  `k = 5` and 24 at `k = 4`, and the literature confirms "not Fano unless `n <= 3`, always weak Fano". The
  Voskresenskii-Klyachko attribution is wrong. Repair: pass to `X'`. No exponent changes.
* Minor items:
  * "new" on Prop 1.3;
  * "any admissible `Z` contains the big diagonal" (diagonal points are never A-points);
  * "rest on `V_E(5)` itself";
  * `alpha^ess = tau` on dP6 was asserted without proof (now proved);
  * twisted `k = 3` rests on an unchecked corollary;
  * "`(L_5^circ)` for that `delta`";
  * two citation details are unconfirmed: the page range of Huang ANT 15, and the "`alpha >= eps/2` (Thm 6.3)" form of
    McKinnon-Roth 2015.

**Outcome.** The survey is correct and useful as a translation. As mathematics it adds one formula (`tau <= 4 - 2^{3-k}`)
and the rank table. It contains one false geometric statement, now repaired. No route.

### 1.2 toric (`toric.md`; `kill_toric.md`): survives as `restricted_theorem`; one major overclaim outside the proved core

Setting: `X = Lbar_5^tw`, `B = -K_X`, `e` the identity, and `Omega_dist = T_5(Q) \ Z_diag` (the distinct tuples). The
killer re-verified every claim below with independent code.

**Survives (killer-corrected).**
1. **Geometry (proved).**
   * `P_B = R_5 = conv(e_a - e_b) = {sum x = 0, ||x||_1 <= 2}`.
   * `B = H + H'`, the pullbacks of `O(1)` under the two toric blow-downs to `P^4`.
   * `B` is nef and big but not ample, and `B^{k-1} = binom(2k-2, k-1)` (`= 70` at `k = 5`).
   * `width_v(R_5) = 2(max v - min v)`.
   * Every 1-PS orbit through `e` descends to a `Q`-curve, since `sigma = -1` on `N`.
2. **Seshadri constant (proved).** `eps_e(B) = 2`, for the split and twisted forms. The lower bound is
   `H.C, H'.C >= mult_e C`. The constant is attained only by the conics `C_{1_A} ⊂ Z_diag`.

   `alpha_e(B) = 2` on all of `T_5(Q)`. The lower bound is the killer's: `f_ab = 2 - x_a/x_b - x_b/x_a ∈ H^0(B)` has
   order 2, and on norm-one points `|f_ab(y)| = |u_a - u_b|^2 >= H_B^{-1}`.

   So every Seshadri-proportional theorem (`alpha >= c eps`, `c <= 1`) gives `<= 2`, which is vacuous. *Qualification
   (killer):* the restricted constant `eps^diag_e` lies in `[10/3, 8]`, so an "essential Seshadri" statement is not capped
   at 2 by this argument.
3. **Dictionary (proved, Lemma 1.3).** `A_5(delta)` for some `delta < 1/30` holds iff `alpha^ess_e(B; Omega_dist) > 15/4`.
   `(L_5^circ)` holds iff the same is true over `eta_0`-balanced, one-unit-class sequences, with
   `(4 - mu)w - O(eta w) <= H_Sigma <= 10(4 - mu)w + O(eta w)`.
4. **Curves (proved).**
   * Along a rational curve with contact `mu_0` at a rational parameter, `alpha = B.C/mu_0` (Prop 4.1).
   * Balanced curves off `Z_diag` have `alpha >= 15/4`, with equality iff the curve is `V`-null (Prop 2.3,
     Cor 4.2). Here `D_V.C = 80 kappa`.
   * Covering rational families through `e` have `alpha >= k - 1 = 4` (Prop 4.3). This is gauss Thm 4.3(ii) restated.
   * **Thm 4.4 (exact, reproduced at 5 parameter points).** The `k = 4` family of gauss Thm 4.1 consists of balanced
     `V_4`-null rational curves: `E_A.C = 1` on all 14 rays, `B.C = 14`, contact 4, `alpha = 7/2`.
   * The `k = 4 -> 5` switch is the comparison of the covering threshold `k - 1` with
     `alpha_bal(k) = 4(1 - 2^{1-k})`: `3 < 7/2`, but `4 > 15/4`.
5. **Toric Liouville is governed by `tau_e` (proved, Thms 3.2-3.4(a)).** Take finitely many auxiliary Laurent
   polynomials, used point by point through Liouville with toric forced divisibility, off their zero sets.
   * For arbitrary profiles they give exactly `A_5(delta)` with `delta > 2/tau_e(B) - 1/2`. The proof is by minimax
     over profiles, together with `f f^*` (whose order doubles and whose Newton polytope is symmetric).
   * `10/3 <= tau_e(B) <= 15/4`, so toric Liouville **never passes `delta = 1/30`**, for any profile. For balanced
     profiles this is the record's `construct.md` Cor C3 / row 20. The new part is the all-profile, all-finite-family
     form.
6. **Exact data.** `t_max(m) = floor(10m/3)` for `m <= 7` (author, mod `2^61-1`; reproduced mod `2^31-1`) and for `m = 8`
   (killer, new: rank 15421 = `|S_8|`). At `k = 4`, `t_max(m) = 3m` for `m <= 9`, while the balanced constant is `7/2`,
   so the all-profile and balanced thresholds really differ.
7. **Conditional (Prop 5.1).** If `-K` is big, Vojta for `(Bl_e X, E, {inf})` gives `alpha^ess >= dim X`. This is
   gauss Thm 4.2 and, at snippet level, McKinnon-Satriano's corrected Prop 1.3. **Not new.**
8. **Evidence.** No nondegenerate real solution of the 0-dimensional minimal-degree `k = 5` system has been found. That
   system is four coupled complex ideal PTE problems of size 8 sharing 15 roots. The author ran 5 methods with about 1700
   starts; the killer found 0 of 44 converged solutions nondegenerate. The same codes find `k = 4` solutions readily
   (11-45%, 24/100). Nothing depends on this: finitely many curves can be put into `Z`.

**Killed / corrected.**
* **F1 (major overclaim).** The claim "toric Liouville with polynomials of degree `<= 7` cannot improve Cilleruelo-Córdoba
  (`delta > 1/10`)" is not supported. The data cap only a **single** `f` with `max_A width_f(1_A) <= 8`, via
  `f f^* ∈ H^0(8B)`. A section of `7B` can have width up to 14. Finite families are governed by products, which have
  large widths. For families, or for larger widths, the claim **is** the open conjecture `tau_e = 10/3`.
* **F4 (literature).** McKinnon-Satriano is not cited.
  * Their *canonical boundedness* is exactly the essential `alpha >= dim X` statement.
  * "(McK^ess) is not in the literature" is inaccurate.
  * "Not McKinnon's conjecture" holds only for the 2007 infimum form.
* **F3.** "No toric argument can prove `A_5`" is a meta-claim. What is proved is a barrier for three named classes:
  * sections, which is a proved cap;
  * `T(R)`-invariant archimedean data, which is a tautology (an invariant of a transitive action is constant);
  * fake-blind tropical data, which is the record's barrier, killer-checked but not refereed.
* **F2, F5, F6 (minor).**
  * The lower bound `alpha_e >= 2` was missing; it is supplied above.
  * 1-PS orbit closures can be cuspidal at the boundary (smooth at `e`, which is all that is used).
  * The width form of the toric height is an upper bound, with equality on squarefree profiles; that suffices for the
    minimax.
  * Coincident roots in §7 are allowed only for nested blocks.

**Outcome.** A correct and careful placement of `A_5` among the approximation constants of `(X', -K, e)`. It proves the
all-profile section cap and supplies exact data for `tau_e`. It changes no route. The window `alpha^ess ∈ (15/4, 4)` is
the whole content of `(L_5^circ)`, and nothing toric enters it.

## 2. Does any known theorem or toric argument reach A_5(delta), delta < 1/30, or (L_5^circ)?

**No.** In `alpha`-language the target is `alpha^ess_e(-K) > 15/4`, and every candidate stops at or below `15/4`. Each
row below has been checked against the owner's constraints. The exceptional sets are fixed by `(X', -K, e)`, with
`S = {inf}`. There is no phase, and no gcd is discarded.

| candidate | what it gives at `e` (w.r.t. `-K`) | reaches `> 15/4`? | why not | status |
|---|---|---|---|---|
| McKinnon-Roth Roth-type (Invent. 2015); any `alpha >= c eps_e`, `c <= 1` | `<= eps_e = 2` (infimum over all sequences) | no | `eps_e = alpha_e = 2` is computed on `C_{1_A} ⊂ Z_diag` | proved (toric Thm 2.1, killer F2) |
| McKinnon-Roth Liouville with stable/asymptotic base locus (EJM 2016); Seshadri / effective-threshold bounds; any fixed section with forced divisibility, Laurent polynomials over `Q(i)` | `alpha^ess >= tau_e(-K)`, with `tau_e ∈ [10/3, 15/4]`; balanced: exactly `15/4` | no | mean-width Thm B (refereed); minimax over profiles | proved (lit Prop 1.3, toric Thm 3.3-3.4(a); row 20) |
| Ru-Vojta, Evertse-Ferretti, Grieve, Heier-Levin expectation constants with `D = E`, `S = {inf}` | `beta(L,E) <= tau` | no | rational target, `E` geometrically irreducible | proved (gauss §4.4; lit §3) |
| Huang torsor lower bounds (ANT 2021, Acta Arith 2019) | split only; integer non-vanishing in torsor coordinates | no (assessment) | section-type; a twisted extension would reprove `<= tau` unless it had a new mechanism | **assessment**, proofs unread |
| McKinnon conjecture, Coba, weak form; all proved cases (McKinnon 2007, MR16, Huang, McKinnon 2310.01640, MMS 2504.05213) | statements about `alpha_e` (`= 2` here) | no | infimum, not essential; split or homogeneous | snippet-level |
| McKinnon-Satriano (TAMS 2021), Monahan-Satriano, Lehmann-McKinnon-Satriano | canonical boundedness (`alpha^ess >= dim`) | **would** reach it, but it is **assumed**, or derived from Vojta | conditional on Vojta; split | snippet-level (judge J3) |
| Batyrev-Tschinkel (anisotropic tori), Huang equidistribution, Huang-de Saxcé | counts at the natural scale, relative error | no | a sparse Zariski-dense set changes no asymptotic; the need is absolute error `< 1` in balls of radius `B^{-1/4-delta/2}` | derived (lit §5 item 3) |
| toric sections | `delta > 2/tau_e - 1/2 >= 1/30` | no | as row 2 | proved |
| `T(R)`-invariant archimedean data (moment map, `Log_inf`, toric metrics) | nothing | no | constant on `T_5(R)`, a tautology | proved, trivial |
| fake-blind tropical / cocharacter data | nothing | no | balanced 5-row fakes at `R''^{7/15}` | record barrier (killer-checked) |
| covering rational families (refutation direction) | `alpha >= 4` | n/a | gauss Thm 4.3(ii) | proved |

**The only route to the uniform bound in this round's language, written out with every gap flagged.** No step is new
mathematics except the elementary step 1. The chain is conditional.

* **Hypothesis (CB_e^norm1).** For every sequence `y_i` of pairwise distinct points of `T_5(Q)`, Zariski dense in `T_5`,
  with `dist(y_i, e) -> 0`: `dist(y_i,e)^gamma H_{-K}(y_i)` is unbounded for every `gamma < 4`. This is McKinnon-Satriano
  canonical boundedness of `X'` at `e`, restricted to `T_5(Q)`.
  * **GAP 1 (definitions).** Their exact definition of "Zariski-dense sequence" and of `alpha` is known only from snippets.
    So is whether they require `X` smooth, `Q`-factorial or terminal; `X'` is Gorenstein and singular, with `e` a smooth
    point. The hypothesis above is stated self-containedly, so the chain does not depend on their wording.
  * **GAP 2 (the hypothesis itself).** It is **open**. It is implied by Vojta for `(Bl_e X, E, {inf})` (gauss Thm 4.2,
    toric Prop 5.1). No unconditional instance of `alpha^ess > tau` at a rational point is known to the record, for any
    variety.
* **Step 1: (CB_e^norm1) ⇒ `A_5(delta)` for every `delta > 0`.** The argument is the killer's; the judge re-derived it.
  * Fix `delta > 0` and let `S_delta ⊂ T_5(Q)` be the set of classes of primitive tuples of 5 distinct lattice points
    with `diam <= R''^{1/2-delta}`. Then `S_delta ∩ Z_diag = ∅`.
  * If `S_delta` is not Zariski dense, its closure `Z` is a proper closed subset, and `A_5(delta)` holds with this `Z`
    and `R_0 = 1`.
  * Otherwise, enumerate `S_delta` by increasing `R''`. There are finitely many classes of bounded `R''`, since a circle
    of radius `R''` has finitely many lattice points. So `R'' -> inf` along the enumeration.
  * `dist(y,e) <= diam/R'' <= R''^{-1/2-delta} -> 0`, and `dist^gamma H_{-K} = dist^gamma R''^2 <= 1` for
    `gamma = 4/(1+2 delta) < 4`.
  * A Zariski-dense set enumerated as a sequence is a Zariski-dense sequence. So this contradicts (CB). `□`
* **Step 2: `A_5(1/64)` ⇒ `M(C) < inf` for every `C`.** This is gauss Thm 3.4: elementary, killer-checked, and it needs no
  unit class. With `delta <= 2^{-6}` only the convexity averaging of Lemma 3.2 and the counting Lemma 3.3 are used, and
  `Z ⊂ T_5` is similarity-invariant.

So `(CB_e^norm1) ⇒` the uniform bound, with no further gap. The weaker hypothesis "`alpha^ess > 15/4` on
`eta_0`-balanced one-unit-class sequences" gives `(L_5^circ)` by toric Lemma 1.3(c), and then the uniform bound by gauss
Thm 3.4′. That route uses the balanced selection of item 532 within one unit class (`kill_gauss` §2). **Nothing in the
chain is proved except steps 1 and 2.** The hypothesis is Vojta-strength and open.

## 3. Judge's own checks (unrefereed)

Scripts and outputs are in `judge_checks/`, run under `prlimit --as=3500000000`; notes are in `judge_checks/notes.md`.
**None of this has been refereed.**

**J1. Reconciling the two killers on McKinnon-Satriano.** Both killers confirm the bibliographic data (TAMS 374 (2021)
3557-3577, DOI 10.1090/tran/8318, arXiv 2004.05212) and a 2023 erratum. They paraphrase the erratum differently:
* kill_lit: "Vojta does not imply canonical boundedness in general, but does for toric";
* kill_toric: "if `-K` is big, Vojta implies canonical boundedness".

The two readings agree if the toric case is the big-`-K` case, which it is: `-K` of a complete toric variety is big. Both
are snippet-level. Nothing in this README depends on the erratum's exact wording. The record's own derivation (gauss
Thm 4.2, toric Prop 5.1, killer-checked) gives the implication Vojta ⇒ `alpha^ess >= dim` for big `-K` independently.

**J2. Split calibration: the regime `tau_e(-K) < dim X` is geometric** (exact: `split.py` -> `split.out`).
* `eps_e`, `tau_e`, `t_max(m)`, the balanced constant `15/4` and the covering threshold `4` depend only on the fan and
  the polytope. So they are the same for the **split** Losev-Manin 4-fold `LM_5` over `Q` at `e`, with points
  `[u] ∈ G_m^5/G_m(Q)` for primitive integer 5-tuples `u`.
* Exact check (300 random primitive tuples, `k = 3..6`): the toric `-K` height is `H_{-K}([u]) = (max u/min u)·lcm(u)`.
* Balanced split blocks `u_j = prod_{A ∋ j} n_A`, over the 30 proper subsets `A`, give the following, verified on actual
  primes:
  * `gcd(u) = 1`;
  * `lcm = prod n_A`;
  * `gcd(u_a,u_b) = prod_{A ⊇ {a,b}} n_A`;
  * hence the balanced Liouville constant `alpha = 30/(15-7) = 15/4`, and `4(1 - 2^{1-k})` in general.
* So the split integer question has the same shape. Take five distinct integers near `U` with `gcd 1`. Do they satisfy
  `diam <= U·lcm(u)^{-(1+2 delta)/4}` only on a proper Zariski-closed set of ratio classes? That is the split
  `alpha^ess_e > 15/4`.
* It has the same window `(15/4, 4)`. It is predicted by the same Vojta instance, and it is open as far as the record
  knows.
* It is **only a calibration**. It neither implies nor is implied by `(L_5^circ)`. Its value is for the expert: if
  `alpha^ess > tau` were known at `e` for split `LM_5`, the mechanism would be the first candidate to transfer.
* *Contrast (from lit/kill_lit, derived):* at `e` on split `P^1 x P^1` and split dP6, `alpha^ess = tau >= dim`. Points are
  repelled, and Liouville is the truth. `LM_5` (split or twisted) is the first permutohedral case where `tau_e < dim X` is proved (`tau_e <= 15/4 < 4`; at `k = 4`, `tau_e >= 3 = dim`).
  There Liouville is predicted **not** to be the truth.

**J3. Literature search (two WebSearch queries, 2026-10-28; snippets only).**
* Query `McKinnon Satriano "Approximating rational points on toric varieties" canonically bounded abstract`: the snippets
  say their split-toric-surface theorem holds *in the presence of the canonical boundedness condition*, conditional on
  Vojta. Their higher-dimensional statement (`Q`-factorial terminal **split**) compares rational curves with Zariski-dense
  sequences. So canonical boundedness is a **hypothesis** there, not a theorem.
* Query `"canonically bounded" rational point approximation constant Vojta toric unconditional dimension`: no
  unconditional canonical-boundedness theorem surfaced. Pasten, arXiv 2205.07841, describes unconditional progress on
  arithmetic Vojta in higher dimension as very limited.
* This agrees with both killers. It is search-negative, not a proof of absence.

Sources: [AMS page](https://www.ams.org/tran/2021-374-05/S0002-9947-2021-08318-X/),
[arXiv 2004.05212](https://arxiv.org/pdf/2004.05212), [arXiv 2205.07841](https://arxiv.org/pdf/2205.07841) (snippets only;
the hosts do not resolve here).

**J4. Evidence for the weighted mean-width inequality (toric open point 1)** (`wmw.py` -> `wmw.out`, `wmw2.py` -> `wmw2.out`).
* The proposed inequality is
  `(W): ord_e f <= (1/6) sum_{|A| = 2} width_{Q_f}(1_A)`.
  It is equivalent to `reg(S) <= (1/6) sum_{|A|=2} width_S(1_A)` for every finite `S ⊂ M`. For `f ∈ H^0(mB)` all widths
  are `<= 2m`, so (W) gives `tau_e(-K) <= 10/3`, hence `tau_e = 10/3`.
* (W) is an equality for every product of root differences (Vandermondes, 3-point Vandermondes, root binomials). This is
  confirmed on all 70 random root zonotopes.
* Tests: 200 random small sets, plus 270 lattice-point sets of random polytopes, random zonotopes and half-space pieces of
  `mR_5` (`|S|` up to 394). `reg` was computed by rank mod `2147483629`, which gives upper bounds for the rational `reg`.
  * **0 violations.** The largest ratio off root zonotopes is `0.93`.
  * The refereed Thm B bound (sum over all 32 subsets `/16`) was checked as a control: never violated.
* **Evidence only.** Mod-`p` `reg` is an upper bound, so a non-violation is not a certificate.
* Even if (W) were proved, it would settle only the intermediate-exponent question, that toric Liouville cannot improve
  Cilleruelo-Córdoba. It would show that this lane is closed, not open, and it has no bearing on `delta < 1/30`.

No other new mathematical claim is made in this README.

## 4. Additions to the barrier map and the expert note

These are proposed edits; the repository was not edited. All of them are **unrefereed**. Where an item says "proved",
it means proved in a write-up and re-derived by its killer.

**Barrier map (overview §5.4): row 20 extended, and one new row.**

| # | input class | why capped | kind; status | source |
|---|---|---|---|---|
| 20 (extended) | gcd-reduced `k`-point circle statements, attacked by **any finite family of auxiliary Laurent polynomials over `Q(i)`** with toric forced divisibility, used point by point, **for arbitrary profiles**; equivalently McKinnon-Roth Liouville with stable base locus, Seshadri and effective-threshold bounds at `e` | best exponent `delta = 2/tau_e(-K) - 1/2`, with `tau_e ∈ [10/3, 15/4]` (mean width, minimax). Balanced profiles: exactly `15/4` (`delta = 1/30`). Seshadri-proportional bounds: `<= eps_e = 2`, computed in `Z_diag` | V; proved (toric Thms 2.1, 3.3, 3.4(a); lit Prop 1.3) | `round22/toric.md`, `lit.md`, kill reports |
| 21 | **approximation-constant theorems at a fixed rational point** in their unconditional form (McKinnon-Roth 2015/2016, Roth-type and expectation constants with `S = {inf}`, Huang torsor bounds, natural-scale counting and equidistribution including Batyrev-Tschinkel for anisotropic tori), applied at `e ∈ T_5` | they bound `alpha_e = 2`; or they are `<= tau_e <= 15/4`; or they work at relative error at the natural scale; or they are split or conditional on Vojta | Red to row 20, or inapplicable; Huang **assessment** (proofs unread); literature search-negative | `round22/lit.md`, `kill_lit.md`, README J3 |

**lit18 (`lit18/tools.md` §2.8): replace the dismissal reason.** The new reason: in the `T_5` formulation the phase is
quotiented out exactly, so Huang, McKinnon-Roth and McKinnon-Satriano are applicable to the fixed rational point `e`.
They are too weak, for the reasons in row 21. The old reason, "a fixed-point theorem would replace a drifting phase by a
constant", remains correct for the `X_M` formulation only.

**Expert note `expert_note/expert_note.md` §3 (Question B): proposed edits.**
1. *Geometry correction.* `Lbar_k^tw` is smooth projective and **weak Fano** for `k >= 4`; it is not Fano.
   `-K = sum_rho E_rho` is nef and big, with polytope the root polytope `R_k`. `-K = phi^*(-K_{X'})` for the small crepant
   contraction onto the twisted Gorenstein toric Fano `X'`, and `X(Q) = X'(Q) = T_k(Q)`. At `k = 5`: `(-K)^4 = 70`.
2. *Approximation reading, upgraded.* Replace the last two bullets of "Approximation reading" with the following.
   * `eps_e(-K) = alpha_e(-K) = 2` (proved). This is attained only on the diagonal conics, which carry no distinct tuples.
   * `10/3 <= tau_e(-K) <= 15/4` (proved). Exact data `t_max(m) = floor(10m/3)`, `m <= 8`, suggest `tau_e = 10/3`.
   * Every section, base-locus or Seshadri argument gives at most `A_5(delta)` for `delta > 2/tau_e - 1/2 >= 1/30`, for
     every profile (proved).
   * *Question B in the literature's terms:* is `alpha^ess_e(-K_{X'}) > 15/4` on norm-one points (for `A_5`), or on
     balanced norm-one points (for `(L_5^circ)`)? The strongest form is McKinnon-Satriano **canonical boundedness of `X'`
     at `e`** on `T_5(Q)`. It implies the uniform bound by §2 above (elementary). It is predicted by Vojta, and the record
     has found it nowhere proved.
   * Replace "Both statements are the killer's own and **not literature-checked**" with: "Literature-checked at snippet
     level in round 22 (two independent searches, plus the offline arXiv mirror): no unconditional theorem with
     `alpha^ess > tau_x` at a rational point was found for any variety. In the split examples checked (`P^1 x P^1`, dP6)
     `alpha^ess = tau >= dim`."
3. *Optional calibration item (J2).* "The same window `(15/4, 4)` occurs for the **split** `LM_5` at `e`, which is a
   statement about five integers with small lcm in a short interval. An unconditional result there would be the first
   candidate mechanism to transfer."
4. *Section 5, item 3 (auxiliary polynomials) of the expert note.* Add: "on the Gaussian model, for arbitrary profiles
   and finite families: `delta >= 2/tau_e - 1/2 >= 1/30` (round 22)".

**Constraint audit.** Each killer checked the four owner constraints item by item, and the judge re-checked the chain of
§2:
* no uniform exceptional set from a fixed-`S` Subspace theorem: the only exceptional sets are `Z(f)`, the hypothetical
  `Z`, and Vojta's `Z_eps` with `S = {inf}`;
* no Gaussian gcd or denominator cost is discarded: points are primitive tuples, `h_{-K} = log R''^2` exactly, and
  `N'^m` is kept in toric Thm 3.2;
* no drifting phase is set to zero: the phase is quotiented by the similarity action, exactly;
* no template-dependent threshold is treated as uniform: `t_max(m)` is used only as a cap on a method; the PTE searches
  are refutation attempts; the constants of §2 Steps 1-2 are absolute.

Killer F1 of toric was an overclaim about the *reach of a method*, not a uniformity violation.

## 5. Scores

Rubric of `round10/README.md` §1:
* (a) correctness after the kill report;
* (b) fake-violated input;
* (c) movement of growth or uniformity.

Each sub-score is out of 10. The overall score is a judgement that weighs (b) and (c) above (a).

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **toric** | 7 | 0 | 2 | **4** | No error in any statement labelled proved. It is a clean placement of `A_5` among `eps_e = 2`, `tau_e ∈ [10/3, 15/4]`, the balanced `15/4` and the Vojta value `4`. It has two genuinely new proved items: the all-profile, all-finite-family section cap via `tau_e` and minimax, and the identification of gauss Thm 4.1's family as `V`-null balanced curves. Exact `t_max` data to `m = 7` (8 with the killer). The `k = 5` PTE search was calibrated. Deductions: F1, a major overclaim about the reach of toric Liouville ("degree `<= 7`"); F4, McKinnon-Satriano's canonical boundedness was missed and "not in the literature" was asserted; the meta-claim "no toric argument"; the missing `alpha_e >= 2`. (b) = 0: sections are fake-satisfied and capped. (c) = 2: a proved strengthening of row 20 and an exact frame for the expert question; no movement of the target. |
| **lit** | 5 | 0 | 2 | **3** | A correct search-negative verdict with careful hypothesis and exceptional-set checks. It gives the exponent dictionary, replaces the lit18 dismissal correctly, and adds the formula `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}` and the rank table. Deductions: F1, the false headline "toric Fano, `-K` ample", misattributed to Voskresenskii-Klyachko; Prop 1.3 called new although it is row 20; the scope of Prop 1.3 stated too narrowly; several unproved asides (dP6, the twisted `k = 3`); two unconfirmed citation details; no full text read (environmental). It missed, as a *named* target, the canonical-boundedness formulation that the killers found, although it cites the paper. (b) = 0. (c) = 2: a barrier-map row (21) and a corrected lit18 entry. |

*Killer reports.* Both were decisive and independent:
* kill_lit found the non-ampleness by direct convexity check and repaired the dictionary through `X'`; it supplied the
  dP6 proof and confirmed 11 of 13 citations;
* kill_toric caught the F1 overclaim, which is a real logical gap (width versus degree), extended the certificate to
  `m = 8`, supplied `alpha_e >= 2`, and found the canonical-boundedness link, which is this round's most useful
  literature item.

For comparison, round 20's gauss lens scored 5, with the same (b) = 2 standing for the target statement itself. Neither
lens here moved (b) or (c) beyond the barrier map.

## 6. Next step, probabilities, stop decision

**Stop internal rounds on growth and uniformity**, as rounds 19 and 20 recommended. This round tested the last named
literature lane: approximation constants at a fixed rational point, which lit18 had dismissed for a reason that does not
apply in `T_5`. The lane is legitimate, and its unconditional content stops at `tau_e <= 15/4`. That is exactly the
record's row 20, now in the literature's language.

**Single most promising next step (external).** Fold §4's edits into `expert_note.md` §3. Then put Question B to experts
in the McKinnon-Roth / McKinnon-Satriano / Huang circle, in their own terms:

> *Is the twisted toric Fano 4-fold `X'` (the root polytope `A_4`, twisted by `-1`, so that `X'(Q) = T_5(Q)` for the
> anisotropic norm-one torus `T_5`) canonically bounded at the identity on its rational points? Less ambitiously: does
> every Zariski-dense sequence of rational points tending to `e` have approximation constant `> 15/4` with respect to
> `-K`, at least along sequences whose 15 toric heights are asymptotically equal? Here `15/4` is both the balanced
> Liouville constant and an upper bound for `tau_e(-K)`, and Vojta predicts `>= 4`.*

A yes gives the uniform Cilleruelo-Granville bound: §2, with gauss Thm 3.4.

*Owner-side enabler (cheap).* Allow `arxiv.org`, `ar5iv.labs.arxiv.org` and `math.uwaterloo.ca` in the session network
settings. That would let a short follow-up read McKinnon-Roth (EJM 2016, Invent. 2015), McKinnon-Satriano (TAMS 2021
plus the erratum) and Huang (arXiv 1809.02001) in full. It would also confirm the canonical-boundedness definition and the
unconfirmed citation details (§1.1). This affects citation accuracy, not the verdict.

**Map-only internal items (optional; no effect on the target):**
1. prove the weighted mean-width inequality (W), hence `tau_e(-K_{LM_5}) = 10/3` (J4 evidence: 0 violations in 470
   sets). This would close the intermediate-exponent lane for toric Liouville;
2. a certified homotopy or monodromy solve of the 0-dimensional `k = 5` PTE system (calibration only);
3. the split integer calibration of J2: is anything unconditional known at `e` on split `LM_5`?

**Honest probabilities.**
* Unconditional growth `M = o(log R/loglog R)` from a further internal round: **< 0.1%**.
* An unconditional uniform bound: **essentially 0**.
* A proof of `(L_5^circ)` or `A_5(delta < 1/30)` in a further internal round: **< 0.5%**. Every known unconditional
  tool stops at `tau_e <= 15/4`, at a rational point, and the record has no non-section mechanism.
* That the literature already contains a theorem implying `(L_5^circ)`: **about 1-2%**. Two independent snippet-level
  searches and the offline mirror found nothing, and canonical boundedness appears only as a hypothesis. Full texts
  were unread, so this is not lower.
* A refutation of `(L_5^circ)` (which would refute `V_E(5)`): **about 1%**. It would not count as progress.
* `tau_e = 10/3` provable internally via (W): **about 25-35%**, of map value only.
* An expert recognising the canonical-boundedness question for `X'` as known or tractable: unknown, and small. It is still
  the only path with non-negligible value.

**Decision: stop.** The record's lasting value is unchanged:
* the `(2/3 + eps)` theorem;
* four conditional forms of the uniform theorem, the smallest on a 4-fold;
* the refereed barrier map, now with row 20 extended to all profiles and finite families, and with a literature row
  (21).

Question B now has a literature-native statement: canonical boundedness, or `alpha^ess > 15/4`, at `e ∈ X'(Q)`.

## 7. Files

* Write-ups: `lit.md` (`lit_checks/`, `lit_notes.md`) and `toric.md` (`toric_checks/`, `notes.md`).
* Kill reports: `kill_lit.md` (`kill_lit_checks/`, `kill_lit_notes.md`) and `kill_toric.md` (`kill_toric_checks/`,
  `kill_toric_notes.md`).
* Judge: this `README.md`, plus `judge_checks/`:
  * `notes.md`;
  * `split.py` -> `split.out` (J2);
  * `wmw.py` -> `wmw.out`, `wmw2.py` -> `wmw2.out` (J4).
* The repository and the overview were not edited.
