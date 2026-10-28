# Round 19 judge: the abc-gain lane for pencils on Mbar_{0,8}

There were two write-ups this round. An adversarial referee attacked each one, and each report sits next to its
write-up:
* `search.md` (the stop-rule test of `round18/README.md` §7.3), with `referee_search.md`;
* `theory.md` (can any pencil on `Mbar_{0,8}` have positive abc-gain?), with `referee_theory.md`.

§1 gives the referee-corrected statements. Remarks that this README adds are marked *README remark*. Only the judge
has checked them, so they are **unrefereed**. Labels used: proved / conditional / exact / evidence (data) /
unrefereed. Paths are relative to this directory (`round19/`, scratchpad). `round18/...` means
`docs/claude_review/round18/...` in the repository.

## 0. Verdict

**No growth improvement and no uniform bound were proved this round.** That holds both unconditionally and
conditionally on any standard conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the
best bound for actual circles. `(L_8)` is still open, and the uniform theorem is still known only conditionally on
Vojta (`ptolemy.md` Thm D).

**The abc lane is closed, by proof.** "The abc lane" here means the unit-pencil lane of round-18 Prop 2.6, which
is the only part of round 18's "abc + auxiliary pencil" gap that carried a proved conditional implication.
* **Theorem 3.3 of `theory.md`** (refereed, proved): every unit pencil `lambda A/C` on `Mbar_{0,8}` has
  `gain <= -2^{8-k(f)} <= -1`. This holds for every degree `D` and every `lambda in C^*`. The same bound holds for
  every composite `h(A/C)`. So Prop 2.6 never applies.
* **The literal stop rule also fires.** `search.md` shows `gain <= -32` exactly for `D = 4, 5, 6` and every rational
  `lambda`. This is exact, complete by proved reductions, and independently reproduced by the referee. The maximum
  is attained, and only at Belyi pull-backs.
* **The two results agree.** All 140 extremal pencils of the search satisfy theory's Corollary 2.5 with equality
  (judge check J1, §4, unrefereed).

**What is still open is not a lane.** Round 18's wider gap also contained non-unit pencils whose base locus lies in
codimension-2 boundary strata. This part is **undecided**. It also has **no proved conditional implication**: a gcd
extension of Prop 2.6 would be needed, and none is proved. No candidate pencil is known. And even a positive answer
would give only `(L_8^pc)` until the nested-chain extension (J2) is done. **No pencil with positive gain was found
anywhere.** The nested-chain extension is therefore moot for this round.

**Scores** (rubric of `round10/README.md` §1, each out of 10; details in §5):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| theory | 8 | 1 | 2 | **5** |
| search | 8 | 1 | 1 | **4** |

**Recommendation (§6).** The round-18 stop rule has fired twice over, once by proof and once by exact search, so
internal rounds on growth and uniformity stop here. The single next step is external: send the consolidated overview
(`overview.md`, now updated through round 19) to experts. It states `(L_8)` in its lattice form `(L_8^pc)`, and
the abc-gain question is now narrowed to non-unit pencils.

## 1. Referee-corrected statements

The full versions are in `referee_theory.md` §6 and `referee_search.md` §7. The condensed versions below keep
every correction.

### 1.1 theory (`theory.md`; report `referee_theory.md`; survives, 7 minor/cosmetic items)

Setting: `X = Mbar_{0,8}`, `beta` the balanced class (`D_T.beta = 1` for all 119 `T`), and
`k(f) = min{|S| : f = g ∘ pi_S}`. In the anchor-8 chart, `k(f) = |I| + [some vertex sum sum_j e_ij != 0]`
(correction m1). `Lambda_k.beta = 32 - 2^{8-k}` for `k = 5..8`, and `Lambda_4.beta = 16`, giving 16, 24, 28, 30, 31.

1. **Theorem 2.4 (proved): beta-weighted log Riemann-Hurwitz.** For every non-constant rational function `f` on `X`,
   `sum_{p in P^1} (L.beta - R_p^int) <= 2 L.beta + Lambda_{k(f)}.beta`. Equality holds for the cross-ratio and for
   `4z(1-z)`. The proof has four steps:
   * `df` spans the line subsheaf `O(div df)` of `Omega^1_X(log ∂)`, and `div df` is computed from the fibres
     (Lemma 2.3);
   * the forgetful filtration of `Omega^1_X(log ∂)` has quotients `pi_k^* psi_k` (Lemma 2.1; the referee checked
     the class identity on all 1701 F-curves);
   * every line subsheaf of `Omega^1_X(log ∂)` has beta-degree `<= Lambda_k.beta` (Lemma 2.2);
   * the only deep input is positivity of `beta` on effective divisors, from Keel-McKernan via the refereed
     `threshold/L8.md` Thm A(a).

   **Corollary 2.5:** `gain(f) <= Lambda_{k(f)}.beta - n'(f) <= 31 - n'(f)`, where `n'(f)` counts the boundary
   divisors on which `f` is constant.
2. **Lemma 3.2 (proved).** Every non-constant unit `u` of `M_{0,8}` satisfies
   `N(u) = #{T : ord_{D_T} u != 0} >= 32`. This is Schwartz-Zippel (Reed-Muller distance) for the quadric
   `sum e_ij x_i x_j` on `{0,1}^7`. It is sharp for cross-ratios.
3. **Theorem 3.3 and Corollary 3.4 (proved).**
   * For every non-constant unit `u` and every `lambda in C^*`:
     `gain(lambda u) <= Lambda_{k(f)}.beta - N(u) - n_1 <= -2^{8-k(f)} <= -1`.
   * When `k(f) = 4`, `gain <= -32`.
   * For every non-constant `h in C(z)`: `gain(h ∘ u) <= Lambda_{k(u)}.beta - N(u) <= -1`.
4. **Consequence (proved), with the referee's scope corrections m2 and m3.**
   * Round-18 Prop 2.6 never applies. **The unit-pencil abc lane is closed by proof.**
   * The *purpose* of the stop rule, "no positive gain", is met for all `D` and all `lambda`.
   * Theory proves the stop rule's literal threshold `-32` only when `k(f) = 4`. For `k(f) >= 5`, theory alone gives
     `-2^{8-k}`, and `-32` is evidence there. `search.md` makes `-32` exact for `D <= 6` (§1.2).
   * The headline "the abc lane as round 18 defined it is closed" overstates the result (m2). Round 18 also listed
     non-unit pencils with boundary base locus, and those remain open.
5. **Not closed (§4 of theory, as corrected).** An abc route through non-unit pencils in
   `P° = {f regular off the codimension-2 strata}` needs both of the following.
   * (i) A gcd extension of Prop 2.6. It is **unproved**: only an unrefereed sketch of the obstruction exists.
   * (ii) A pencil with `n'(f) <= 31 - 2^{8-k(f)}`, `k(f) >= 5`, not composite with a unit, whose three non-reduced
     fibres satisfy `sum 1/mu_p < 1`. These necessary conditions are proved. **No example is known**, and
     Theorem 2.4 alone allows `gain <= 31`.
6. Minor and cosmetic corrections:
   * m4: the Campana remark must be stated on the resolved pencil `X~`. It is not used.
   * m5: the smallest `N` found at `k(f) = 8` is 44, not 50 (evidence).
   * m6: the 14-`lambda` run is a consistency test, not a special-`lambda` sweep.
   * m7: the "any movable `gamma`" remark needs the maximum over orderings. It is not used.

### 1.2 search (`search.md`; report `referee_search.md`; survives, 4 minor items + 1 remark)

A unit pencil here is `lambda A/C`, where `A` and `C` are coprime products of `D` differences in Kapranov's chart
(point 8 at infinity) and `lambda in Q^*`.

1. **Proved for every `D`:**
   * the excess identity `gain = -N0(A,C) + sum_{Y ⊂ F_1} (m_Y - 1) Y.beta`, with `N0 = N(u)`;
   * `N0 >= 32`, proved independently of theory's Lemma 3.2 by inclusion-exclusion on 32 disjoint quadruples;
   * for `lambda != ±1`, the fibre over 1 has no multiple boundary component. So the gain exceeds `-N0` only through
     a repeated interior factor of `B_lambda`. Every such `lambda` is a rational root of
     `Disc_D(C|_l - lambda A|_l)` on any one line `l` where this discriminant is not identically 0. This makes the
     candidate set **provably** complete (round 18's was complete only with probability 1);
   * a fast exact path for `lambda = ±1`.
2. **Exact (complete by these reductions; independently reproduced).**
   * For `D = 4, 5, 6` and every `lambda in Q^*`, the gain is `<= -32`, and `-32` is attained.
   * The genuinely special `lambda != ±1` number 19, 22 and 31 at `D = 4, 5, 6`. All of them have gain `-32`.
   * At `D = 5` the second value is `-44`.
   * The referee re-ran the complete `D = 4, 5, 6` searches with independent code: different lines and an
     interpolated discriminant. Histograms and special sets are identical entry by entry. This closes the author's
     caveat that `D = 6` was not re-verified pencil by pencil (m3).
   * Together with round 18 (`D <= 3`), this covers all `D <= 6`.
3. **Data.** Every recorded pencil with gain `-32` is a forgetful pull-back `h ∘ chi` of a genus-0 Belyi map.
   * The remark that the maximum "shows no drift upward with `D`" is six data points and has no bearing on `D >= 7`
     (m2).
   * m5: the result covers the `S_8`-saturation of the searched class.
4. Minor corrections: m1, the flint bivariate discriminant needs the universal-discriminant normalisation (harmless,
   independently confirmed); m4, a python-flint `factor()` overflow raises an error and never returns a wrong value.
5. **Decision, not theorem.** The stop rule fires. The search itself proves nothing for `D >= 7` or for non-unit
   pencils. For unit pencils, theory's Theorem 3.3 now covers every `D`.

## 2. Status of the abc lane

| part of the lane | status after round 19 | source |
| --- | --- | --- |
| unit pencils `lambda A/C`, any `D`, any `lambda in Q^*` (indeed `C^*`): the hypothesis of Prop 2.6 | **closed by proof**: `gain <= -2^{8-k(f)} <= -1` | theory Thm 3.3 (refereed) |
| composites `h(A/C)` | **closed by proof**: `gain <= -1` | theory Cor 3.4 (refereed) |
| stop rule, literal threshold `-32` | **fires.** Proved for `k(f) = 4`; exact for all unit pencils with `D <= 6` | theory Thm 3.3; search (refereed, reproduced) |
| non-unit pencils with base locus in codimension-2 boundary strata | **undecided**. No conditional implication is proved (the gcd extension of Prop 2.6 is missing), and no candidate exists. Theorem 2.4 caps the gain at `<= 31 - n'(f)` and gives proved necessary conditions | theory §4 (refereed as corrected, m2) |
| nested-chain extension `(L_8^pc) -> (L_8)` (J2) | **moot** for unit pencils. Still needed after any future positive non-unit pencil | `round18/README.md` §3.3 |

**Classification asked for by the task:** **closed**, by proof (theory) and also by the stop rule (search). It is
**not open**, because no positive-gain pencil was found. The residual non-unit question is **undecided**, but it is
not a lane: it has neither a proved implication nor a candidate.

*README remark (unrefereed; a reading of refereed results).* The margin in Theorem 3.3 is exactly
`32 - 31 = 1`. The two independent write-ups prove `N(u) >= 32` by two different methods: Schwartz-Zippel, and
inclusion-exclusion on quadruples. The data's `-32` comes from `k(f) = 4`, where theory is sharp. The remaining
non-unit question asks for a pencil with three heavily non-reduced fibres and almost no boundary in them. The record
contains no such pencil on `Mbar_{0,8}` other than the composites `h ∘ chi`.

## 3. Consequences for the uniform bound and growth

**Unchanged.** `M <= (2/3 + eps) log R/loglog R` is still the best bound for actual circles. The uniform theorem is
still conditional on Vojta, through `(L_8)` and `ptolemy.md` Thm D. `(L_8)` is still open.

**Effect on conditional theorems.** Round 18 listed one conditional row whose hypothesis was open: "abc + a unit
pencil with gain `> 0` => `(L_8^pc)`". That hypothesis is now **refuted**: no such pencil exists (theory
Thm 3.3). The row therefore yields nothing. Even if it had held, it would have given only an abc-conditional
`(L_8^pc)`, which is not growth progress. So the round removed a (low-probability) path to an abc-conditional uniform
theorem. It removed no path to unconditional growth, because none ran through this lane.

**What does count.** Only unconditional growth `M = o(log R/loglog R)` or a uniform bound counts as progress.
Neither moved. Rubric (b) is also unchanged: the lane concerned `(L_8)`, which is fake-blind by design.

## 4. Judge's own checks (unrefereed)

Scripts are in `judge_checks/`. All arithmetic is exact (python-flint `fmpz_mpoly`, integers). Jobs ran under
`prlimit --as=3500000000`. **None of this has been refereed.**

* **R19-J1** (`j1_bound_vs_records.py` -> `j1_out.txt`, 0.4 s): tests theory Thm 3.3 on search data. The input is
  every one of the 140 pencils with gain `> -40` recorded in `search_checks/out/D{4,5,6}.jsonl`: 43, 42 and 55
  pencils, all with gain `-32`. For each pencil the script computes the following quantities from scratch, without
  using the author's or the referees' code:
  * `k(f)` by the anchor rule of referee m1;
  * `N(u)` by enumerating the 119 `T`;
  * `n_1 = #{T : nu_T(B) > mu_T}`, with `nu_T(B)` computed exactly by the substitution `x_i -> x_{i0} + t x_i`.

  Result:
  * Every one of the 140 pencils has `k(f) = 4`.
  * Thm 3.3's bound holds on all 140.
  * It is an equality on 106 of them: 34 with `(N, n_1) = (32, 16)` and 72 with `(48, 0)`. The 72 are all the
    special `lambda`.
  * On the other 34, `(N, n_1) = (32, 0)` and the bound is `-16`. Of these, 28 have `lambda = -1` (Ptolemy at
    `lambda = -1` is of this type) and 6 have `lambda = 1`. For these 34 the slack comes from 16 further boundary
    divisors on which `f` is constant at a value other than `0`, `1` or `infinity` (J1b).
* **R19-J1b** (`j1b_nprime.py` -> `j1b_out.txt`): tests Corollary 2.5 in its `n'` form on the same 140 pencils.
  * `n'(f)` counts the `T` with `e(T) != 0`, plus the `T` with `e(T) = 0` whose initial forms `in_T(A)` and
    `in_T(C)` are proportional. These are exactly the `D_T` on which `f` is constant.
  * Result: every pencil has `n'(f) = 48`, and `gain = Lambda_4.beta - n'(f) = 16 - 48 = -32`. **Corollary 2.5
    holds with equality on all 140 extremal pencils of the search.**
  * Reading (unrefereed): the two independent write-ups agree exactly at the extremum. The search's maximisers are
    precisely the equality cases of theory's inequality.
* **No other new mathematical claim.** Everything else in this README is a referee-corrected statement or a reading
  of the record.

## 5. Scores

The rubric is that of `round10/README.md` §1: (a) correctness after refereeing, (b) fake-violated input,
(c) movement of growth or uniformity. Each sub-score is out of 10. The overall score is a judgement that weighs (b)
and (c) above (a). Nothing this round moves growth or uniformity. (c) credit is given only for a proved change to
the barrier map. (b) is 1 throughout, because the lane is fake-blind.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **theory** | 8 | 1 | 2 | **5** | The cleanest result in several rounds: a new, sharp, proved inequality (Thm 2.4), and a proof that closes the only uncovered conditional mechanism for every `D` and every `lambda`. The referee re-derived every step and found no fatal or major defect. Deductions: the definition of `k` was ambiguous (m1, which would have doubled the claimed bound); the headline overstated the scope (m2: non-unit pencils were in round 18's lane); "stop rule in a stronger form" holds only for the sign (m3); four cosmetic items. (c) = 2 for a proved correction to the barrier map: the gap is closed. It is a negative result for every conditional route through it, and it does not move growth. |
| **search** | 8 | 1 | 1 | **4** | It did exactly the assigned task and went beyond it. It found proved reductions (excess identity, `N0 >= 32`, a provably complete special-`lambda` set) that made the search exact rather than probability-1. It covered `D = 4, 5, 6`, and the referee reproduced it entry by entry. Deductions: minor items only (m1-m4), plus a trend remark over six data points. Its sign conclusion is superseded by theory Thm 3.3. Its lasting content is that the threshold `-32` is exact for `D <= 6`, which theory proves only for `k(f) = 4`. |

## 6. Recommendation for the program

**Stop internal rounds on growth and uniformity.** Round 18's recommendation stands, and its only exception, the
abc-gain search, has now been decided. Its stop rule fired, and theory went further and closed the lane by proof.

1. **Done this round (§7).** The round-19 lines in the consolidated overview (`overview.md` §7.1, §6.6, Summary
   item 5, the abstract, the route table and §5.5) now give the referee-corrected outcome.
2. **External step (unchanged from `round18/README.md` §7.3, now final).** Put `(L_8)` before experts in Diophantine
   approximation, with three items:
   * `(L_8)` in its lattice form `(L_8^pc)`;
   * the narrowed abc-gain question: does any non-unit pencil in `P°` on `Mbar_{0,8}` have positive gain, and does
     Prop 2.6 extend to such pencils?
   * theory Thm 2.4 as a possibly reusable tool.

   The relevant expertise is in gcd inequalities and Vojta's conjecture for blow-ups of projective space, and in
   abc.
3. **Not recommended:**
   * any further unit-pencil search, at `D >= 7` or otherwise, which Thm 3.3 makes pointless;
   * the nested-chain extension, which is moot without a positive pencil;
   * a search for positive-gain non-unit pencils without first proving the gcd extension. Even with both pieces, the
     result would be an abc-conditional `(L_8^pc)`, not growth;
   * everything on the round-18 "not recommended" list.

**Honest assessment.**
* Unconditional growth `M = o(log R/loglog R)` from any further internal round: essentially 0.
* An unconditional uniform bound: essentially 0.
* An abc-conditional uniform theorem through the residual non-unit lane in one round: well below 1%, because it
  needs three unproved steps (a pencil, the gcd extension and J2). Even then it would be conditional, so it would not
  count.
* The record's lasting value is unchanged: the `(2/3 + eps)` theorem, the conditional theorems and the refereed
  barrier map. The barrier map now has one gap fewer.

## 7. Overview edits and files

**Edits to `overview.md` (scratchpad copy only; the repository is read-only):**
* the main §7.1 "Round 19" placeholder is replaced by a 5-line referee-corrected summary;
* the short placeholders in §6.6 and in Summary item 5 are replaced by one-line pointers;
* the "round-19 search in progress" phrases in the abstract, the route table (§6, abc row) and §5.5 are updated;
* the compile line now says "through round 19".

**Files (this directory):**
* Write-ups: `theory.md` + `theory_checks/` (notes `theory_notes.md`); `search.md` + `search_checks/` (notes
  `search_notes.md`).
* Referee reports: `referee_theory.md` + `referee_theory_checks/`; `referee_search.md` + `referee_search_checks/`.
* Judge: this `README.md`; `judge_checks/j1_bound_vs_records.py` -> `j1_out.txt`; `judge_checks/j1b_nprime.py` ->
  `j1b_out.txt`; `judge_checks/notes.md`.
