# Kill report on round22/lit.md

Status: complete (round 22, adversarial killer of `lit.md`). Scripts and outputs: `kill_lit_checks/`; working
notes: `kill_lit_notes.md`. Labels: [D] derived here, [E] exact computation here, [S] web-search snippet (query
recorded in §5), [A] abstract in the offline arXiv RSS mirror (`scratchpad/lit18/tools/db.json`), [R] recalled.
No full text of any cited paper was reachable (springer, uwaterloo, arXiv hosts fail DNS in this session).

## 0. Verdict

**The verdict `no_new_result` survives. The write-up contains one major mathematical error, which can be repaired, and
about ten minor defects.** No literature theorem implying `A_5(delta)` (`delta < 1/30`) or `(L_5^circ)` turned up in my
own searches or in the offline mirror. All exact numbers were reproduced by independent code.

* **Major (F1, false statement in the headline and in Lemma 1.1).** "`X = Lbar_5^tw` is a toric Fano 4-fold; `-K_X` is
  ample (the permutohedral variety is Fano, Voskresenskii-Klyachko)." This is **false for every `k >= 4`**. `-K_X` is
  nef and big but not ample ([E] `fan.py`: 480 strict-convexity violations at `k = 5`, for example the chamber `id` and
  the ray `1_{{0,2}}` with `<m_sigma, 1_A> = 1`; the literature agrees, [S] Q2). The verdict survives the repair (§2.2):
  `X -> X'` is a small crepant toric contraction onto the Gorenstein toric Fano `X'` of the reflexive root polytope `A_4`
  ([E] `facets.py`). It is an isomorphism on `T_5`, and `X(Q) = X'(Q) = T_5(Q)`. Every literature statement that needs an
  ample `L` must be applied to `(X', -K_{X'})`, not to `(X, -K_X)`. Heights and distances on `T_5(Q)` are unchanged, so
  every exponent in lit.md stays the same.
* **Prop 1.3 survives** (correct; [D] re-derived, §3). It is a transcription of the record's own `gauss.md` §4.4 bullet 2
  (`construct.md` Thm B; row 20), not new mathematics. The only genuinely new items are its specialisation
  `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}` and the exact ranks.
* **Exact ranks survive** ([E] `reg2.py`, `reg2_k5.py`, `vander.py`). The maximal `ord_e` of a section of `m(-K)` is
  `3, 6, 10, 13, 16` for `m = 1..5` at `k = 5`. These are mod-`p` upper bounds at two primes, and the value `10` at
  `m = 3` is certified over `Q` by the Vandermonde. So `max ord/m = 10/3` for `m <= 5`. At `k = 3` the value is `3m`
  (`m <= 5`), and at `k = 4` it is `3m` (`m <= 4`).
* **Citations.** 11 of the 13 bibliographic items I checked are confirmed (§5). The other two are not confirmed: the
  page range of Huang ANT 2021, and the "`alpha_x(L) >= epsilon_x(L)/2` (Thm 6.3)" form of McKinnon-Roth 2015. The
  precise base-locus hypothesis of McKinnon-Roth 2016 is unverified; snippets say "asymptotic" in some places and
  "stable" in others. Everything that lit.md says about the *mechanisms* of Huang's and McKinnon-Roth's proofs is
  assessment, as lit.md itself admits, and it stays assessment.
* **Owner's constraints**: no violation (§6). The re-examination of the lit18 "drifting phase" dismissal is correct.

## 1. What lit.md claims (restated)

1. A dictionary: `L = -K_X`, `H = R''^2`, `dist ≍ diam/R''`. An A-point at `delta` has McKinnon-Roth exponent
   `gamma = 4/(1+2 delta)`. The target is `alpha_ess,e(-K) > 15/4` on balanced points. Liouville gives `15/4`
   (balanced) and `10/3` (all points), Vojta `V_E(5)` gives `>= 4`, and the diagonal conics give `alpha_e <= 2`.
2. Prop 1.3: every section, Seshadri or effective-threshold argument at `e` is capped at `15/4` on balanced points, and
   `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}`.
3. Exact ranks: `ord/m <= 10/3` for `m <= 5`, attained by the Vandermonde.
4. The remaining literature does not apply: McKinnon's conjecture and its proved cases (they concern `alpha`, not
   `alpha_ess`), Roth-type theorems, counting results, and Batyrev-Tschinkel.
5. The lit18 drifting-phase dismissal is void in the `T_5` formulation.

## 2. Normalisation and exponent arithmetic

### 2.1 Dictionary: checked, correct ([D], against `gauss.md` Lemma 1.1, Prop 3.1, Thm 4.2)

* `h_{-K}([zeta]) = W' = log R''^2`. Each Galois pair `E_A + E_{A^c}` has height exactly `log n_A`, because the local
  heights at `pi` and `pibar` carry weight `1/2` each. On `T_5(Q)` every Galois-invariant toric height is
  `sum_c a_c log n_c >= 0`, so all toric heights with every `a_c > 0` are mutually comparable there.
* A-point: `diam <= R''^{1/2-delta}` gives `dist <= H^{-(1+2delta)/4}`, so `gamma = 4/(1+2delta)`. At `delta = 1/30`,
  `gamma = 15/4`; at `delta = 1/10`, `gamma = 10/3`; at `delta = 0`, `gamma = 4`. Correct.
* A counterexample to `(L_5^circ)` with constants `(c, eta_0)` has diameter exponent `<= 7/15 + O(eta_0 + c)`
  (gauss Prop 3.1), hence `alpha <= 15/4 + O(eta_0 + c)`. So the §8 form in lit.md (some `gamma > 15/4`, `eta_0`, `Z`,
  `B_0`) implies `(L_5^circ)`: take `eta_0' <= eta_0` and `c` small enough that `15/4 + O(eta_0' + c) < gamma`.
  Correct. Its wording "gives `(L_5^circ)` for that `delta`" is loose, since `(L_5^circ)` has no `delta` (minor).
* Vojta: `lambda_E + h_{K_Y} >= n lambda_E - h_{-K} - O(1)`, because `h_E >= lambda_{E,inf} - O(1)`. So `V_E` gives
  `alpha_ess >= n - eps`. Correct. It needs only `h_A <= c h_{-K} + O(1)` on `T_5(Q)`, which holds for big `-K` (§2.2).
  It does **not** need `-K` ample, contrary to lit.md's phrase "any variety with `-K_X` ample".
* 1-PS conics: `-K.C_v = width_v(P) = 2(max v - min v)` ([E] `fan.py`, 2000 random `v`). Diagonal conics give
  `gamma = 2`. Correct. But the phrase "the big diagonal, which any admissible `Z` contains" is wrong as stated: their
  points are not tuples of *distinct* points, so they are never A-points, and no `Z` needs to contain them (minor).
* `P^1 x P^1` example (§5 item 4): `tau = 4`, the dense set `(1+1/q, 1+1/s)` has `dist ≍ 1/q` and `H ≍ q^4`, and every
  point off the two lines has `dist >= H^{-1/4}`. Correct ([D]).
* Split dP6 (`= Lbar_3`): `tau_e(-K) = 3` ([E] `reg2.py`; Vandermonde in `H^0(-K)` with `ord 3`). lit.md §8 also asserts
  `alpha_ess = tau` for it without proof. Here is a proof ([D]). `H_{-K}` of a primitive integer triple `x` near `e` is
  `≍ lcm(x_1, x_2, x_3)`. Take `x = (P_1P_2, P_2P_3, P_1P_3)` with `P_i = n + c_i`. Then `dist ≍ 1/n` and `H ≍ n^3`, so
  `gamma = 3`. The union over `c ∈ Z^3` is Zariski dense. So `alpha_ess = 3 = tau`, and the claim is true.
* Twisted `k = 3`: "Liouville-sharp" rests on the round-20 killer's projection corollary, which is labelled "not
  independently checked" (minor labelling defect).

### 2.2 F1: `-K_X` is not ample. Repair ([E] `fan.py`, `facets.py`; [S] Q2)

* On the chamber of the permutation `s`, the linear function of `-K` is `m_s = e_{s_1} - e_{s_k}`. It depends only on
  the first and last entries of `s`. Strict convexity fails whenever `k >= 4` (480 violations at `k = 5`, 24 at `k = 4`,
  none at `k = 3`). There are 20 distinct `m_s` at `k = 5`, namely the vertices `e_i - e_j` of the root polytope.
* All 30 inequalities `<u, 1_A> <= 1` are facets, and every vertex lies on 8 of them, so the polytope is reflexive and
  not simple.
* Hence `-K_X = phi^*(-K_{X'})`, where `phi: X -> X'` is the small crepant toric contraction onto the toric variety `X'`
  of the root polytope. `X'` has the same 30 rays, is Gorenstein Fano and is singular. The twisted form is defined over
  `Q` because the polytope is symmetric under `u -> -u`.
* No nonzero cone of either fan is Galois-stable. Galois maps the chain of a cone to the complementary chain, and for
  `X'` it maps a cone `sigma` to `-sigma`. So `X(Q) = X'(Q) = T_5(Q)`, matching the abstract of Batyrev-Tschinkel
  alg-geom/9411009 [S Q8].
* Literature: "the permutohedral variety ... is not Fano unless `n <= 3`, but always weak Fano" ([S] Q2, snippet
  attributed to *Fano and weak Fano Hessenberg varieties*, Michigan Math. J. 73(3), arXiv 2003.12286). The citation of
  Voskresenskii-Klyachko for "the permutohedral variety and its twisted forms are Fano" is therefore a misattribution of
  a false statement. Q1 shows only that VK classify centrally symmetric toric Fano varieties.
* **Consequence.** `alpha_x(L)` (McKinnon-Roth), McKinnon's conjecture and Huang's constants are all stated for ample
  `L`. The dictionary must therefore read "`L = -K_{X'}` on `X'`" (with `e` a smooth point of `X'`), or use an ample
  toric `L` on `X` with a rescaled exponent. Every number in lit.md stays the same, because `phi` is an isomorphism near
  `e` and on all rational points. So the error does not change the verdict, but the statement as written is false.

## 3. Prop 1.3 (section/Seshadri cap 15/4): re-check

**(a) Correct** ([D]).
* `construct.md` Thm B (read: `reg(S) <= 2^{1-n} sum_{P ⊆ [n]} w_P(S)` for any finite `S ⊂ Z^n`) applies with `n = k`
  to supports in `M = {sum u = 0}`. A degree-0 homogeneous `f` is constant along the diagonal `G_m`, so its order at
  `(1,..,1)` in `A^k` equals its order on the quotient torus.
* `Q_L` has width `<= 2 a_c` in the directions `1_A` and `1_{A^c}`, because `<u,1> = 0`. So
  `ord <= 2^{1-k} · 4m sum_c a_c = 2^{3-k} m sum_c a_c`.
* For `-K` the widths are all exactly `2` ([E] `fan.py`), which gives `4 - 2^{3-k}`: `3, 7/2, 15/4, 31/8` for
  `k = 3..6`.

**(b) Correct, with a scope gap in the statement (minor).**
* The forced-divisibility (tropical) local height of a support `S` at a balanced point is
  `sum_layers (1/2) w_A(S) log p`.
* Combining this with Thm B gives exactly `-log dist <= (w/4) sum_A w_A / ord + O(1)`, and the right-hand side is
  `>= 2^{k-3} w = 4w` at `k = 5`. So the cap is exactly `15/4`, never more.
* Statement (b) quantifies only over Galois-invariant toric `L`. The decisive sections, however, are non-symmetric
  Laurent polynomials over `Q(i)`; the pair difference `w_a - w_b` is the one that realises the balanced Liouville
  bound. The proof covers them verbatim, because widths are invariant under translation and sign. The statement should
  say "every Laurent polynomial over `Q(i)`".
* Extra finite divisibility from contact with `e` or with the diagonals is bounded by the residues `|t_ab|`, which are
  automatically small on the relevant points.

**Novelty: overstated (minor).** (a) and (b) are `gauss.md` §4.4 bullet 2 ("forced-divisibility Liouville excludes
balanced points iff `2^{k-2} ord_e f > sum_cuts width`"; `gamma_hat` movable), i.e. row 20. lit.md's own §1.3 last line
says so, but its §0 calls Prop 1.3 "new here".

**Scope (not an error, but should be stated).**
* "Every section argument" means a *fixed* section with *forced* (profile) divisibility.
* Point-dependent auxiliary sections (Siegel's lemma per point) and multi-point Thue-Roth arguments are covered only
  through the separate `beta <= tau` remark. That remark is correct for Ru-Vojta / Evertse-Ferretti-type bounds with
  one irreducible `D = E` and `S = {inf}`.
* "Covers Huang-type torsor bounds" is assessment only (unread proofs). lit.md labels it so, and it must stay labelled.

## 4. Rank computation (ord/m <= 10/3 for m <= 5): independent recomputation

The method differs from `lit_checks/reg.py` in five ways:
* lattice points are enumerated from the full facet list, not from `||u||_1`;
* the coordinates are the *last* `k-1` entries;
* the basis is the integer-valued binomial basis `prod binom(u_i + m, a_i)`;
* ranks are computed by python-flint `nmod_mat` at the primes `1000003` and `998244353`;
* "full rank of `Poly_{<=d}` on `S` mod `p`" certifies `ord <= d` over `Q`.

| `k` | `|mP ∩ M|`, `m = 1..` | max `ord_e` mod `p` (both primes agree), `m = 1..` | `ord/m` |
|---|---|---|---|
| 3 | 7, 19, 37, 61, 91 | 3, 6, 9, 12, 15 | 3 (= Thm B bound, so `tau_e = 3` exactly with the Vandermonde) |
| 4 | 13, 55, 147, 309 | 3, 6, 9, 12 | 3 (< 7/2) |
| 5 | 21, 131, 471, 1251, 2751 | 3, 6, 10, 13, 16 | 3, 3, 10/3, 3.25, 3.2 |

* **Vandermonde** ([E] `vander.py`, exact over `Z`). Its shifted support (120 points) lies in `3P` and not in `2P`. All
  moments of degree `< 10` vanish and some moment of degree 10 does not, so `ord_e = 10` exactly. This confirms
  `tau_e(-K_{X_5}) >= 10/3`.
* **Ehrhart.** The fourth difference of `1, 21, 131, 471, 1251` is `70`, so `(-K)^4 = 70`. This is the volume of a nef
  and big divisor, and it agrees with the degree of `X'`.
* **Caveat (lit.md already states it).** Levels `m <= 5` do not determine `tau_e(-K) ∈ [10/3, 15/4]`. Level 6 needs
  `ord >= 21` to beat `10/3`, and the square of the Vandermonde reaches only `20`.

## 5. Citations: independent searches (queries and what the snippets support)

| # | query (abridged) | result | lit.md item |
|---|---|---|---|
| Q1 | Voskresenskii Klyachko toroidal Fano varieties root systems permutohedral | VK, Math. USSR-Izv. 24 (1985) 221-244: finiteness of toric Fano varieties; classification of centrally symmetric ones; no "permutohedral is Fano" | **misattributed** (F1) |
| Q2 | "permutohedral variety" "is Fano" if and only if dimension | snippet (Michigan Math. J. 73(3), arXiv 2003.12286): "not Fano unless `n <= 3`, always weak Fano" | contradicts lit.md; agrees with `fan.py` |
| Q3 | McKinnon Roth "An analogue of Liouville's theorem..." | EJM 2016, arXiv 1306.2977: Liouville analogue on arbitrary varieties plus an extension with the asymptotic/stable base locus; cubic surfaces, all nef `L`, `x` off lines; on `P^1`, `alpha_x(O(1)) >= 1/d` and Roth `>= 1/2` (same normalisation as lit.md) | confirmed; exact hypotheses unseen |
| Q4 | McKinnon Roth Invent "epsilon_x(L)/2" | Invent. Math. 200 (2015) 513-583 confirmed; the `eps/2` (Thm 6.3) form **not** found | partly unconfirmed |
| Q5 | Huang "Rational approximations on toric varieties" ... Pagelot | arXiv 1809.02001: split toric, Salberger torsor, best approximations on minimal-degree rational curves, McKinnon's conjecture; "essential constant (Definition 4.7) first introduced by Pagelot [Pag08]" | confirmed; ANT 15 pages 461-512 not confirmed |
| Q6 | McKinnon JAG 16 (2007) Conj 2.7 | JAG 16 (2007) 257-303, DOI 10.1090/S1056-3911-06-00458-9; conjecture as paraphrased; snippet: Castañeda's thesis extended it to split smooth rational surfaces | confirmed ("San19" item still unidentified) |
| Q7 | McKinnon Satriano TAMS 374 3557 erratum | TAMS 374 no. 5 (2021) 3557-3577, DOI 10.1090/tran/8318; erratum 2023: Vojta does not imply canonical boundedness in general, but does for toric | confirmed |
| Q8 | Batyrev Tschinkel anisotropic tori IMRN 1995 | IMRN 1995 no. 12, 591-635 (arXiv alg-geom/9411009): Poisson summation, Batyrev-Manin conjecture, "all `K`-rational points of `X` are contained in `T(K)`"; JAG 1998 paper alg-geom/9510014 covers "smooth projective toric varieties", Galois-invariant fans, no Fano hypothesis | confirmed; applies to `X` |
| Q9 | offline mirror lookup | 2403.02480, 2504.05213, 2509.03998, 2308.07911, 2111.01509 (v3), 2510.09954, 2509.05300 (math.GM): abstracts match lit.md; a keyword scan (approximation/essential constant, anisotropic torus, best approximation, McKinnon) finds nothing on anisotropic tori | confirmed |
| Q10 | arXiv 2308.11847 Monahan Satriano | Vojta implies a curve through `P` approximating better than any Zariski-dense sequence, split horospherical; Vojta implies canonical boundedness | confirmed |
| Q11 | arXiv 2310.01640 McKinnon | cubic hypersurfaces with a rational line; approximation constant computed; Coba | confirmed |

Not re-checked: Kleinbock-Merrill, Fishman-Kleinbock-Merrill-Simmons, Grieve 1705.01023 (peripheral; lit.md does not
rely on them), and Huang IJNT 2017 / Acta Arith. 2019 (bibliographic only). Nothing in §5 changes the verdict. No
search turned up an unconditional theorem with `alpha_ess > tau_x` at a rational point. I cannot certify that none
exists, so the gap statement of lit.md §8 remains a search-negative, labelled claim.

## 6. Owner's constraints

* *Drifting phase.* `[zeta]` is a ratio vector, and the common phase is quotiented by similarity. gauss Thm 3.4' (killer-
  checked) passes from arcs in any position to `A_5` using only `Z ⊂ T_5`. The target `e` is genuinely fixed and
  rational. lit.md's §7 is correct: the lit18 reason is void here, and the theorems are applicable but too weak.
* *Gaussian gcd / denominators.* The height is that of the primitive reduced tuple, and the cost is paid in Lemma 3.2 /
  the balanced selection. Nothing is discarded.
* *Fixed-`S` exceptional sets.* The literature results invoked have exceptional sets fixed by `(X, L, e, gamma)` with
  `S = {inf}`: base loci, accumulating subvarieties, Vojta's `Z_eps`. None is a Subspace set over the primes of `N`.
* *Template-dependent thresholds.* None. The constants depend on the fixed section or variety only.

One thing for the record: if any literature result were used through `(X, -K_X)` as written, its ampleness hypothesis
would fail (F1). It must go through `X'`. That is a hypothesis check, not a constraint violation.

## 7. What survives

**Survives (corrected):**
1. The verdict: no theorem found in the literature implies `A_5(delta)` (`delta < 1/30`) or `(L_5^circ)`.
2. The dictionary (`gamma = 4/(1+2delta)`; target `alpha_ess > 15/4` on balanced points; Liouville `15/4` / `10/3`;
   Vojta `>= 4`; diagonal conics `2`), with `-K_X` read as `phi^*(-K_{X'})` on the Gorenstein Fano `X'`.
3. Prop 1.3: correct, a restatement of row 20 / gauss §4.4. New only in the formula `tau_e(-K_{Lbar_k}) <= 4 - 2^{3-k}`.
4. The exact ranks (`3, 6, 10, 13, 16`; Vandermonde `10/3`).
5. The re-examination of the lit18 dismissal.

**Killed / must be corrected:**
* F1 (major): "toric Fano", "`-K_X` ample", and the VK attribution.
* Minor:
  * the "new" label on Prop 1.3;
  * the scope of Prop 1.3(b), which should cover all Laurent polynomials over `Q(i)`;
  * "any admissible `Z` contains the big diagonal";
  * "rest on Vojta, i.e. on `V_E(5)` itself", which is Vojta for *other* varieties;
  * `alpha_ess = tau` for split dP6 asserted without proof (true, §2.1);
  * twisted `k = 3` resting on an unchecked corollary;
  * "`(L_5^circ)` for that `delta`";
  * the unconfirmed MR15 `eps/2` form and the ANT page range.

**Corrected statement of the main result.**
* No theorem in the literature that I or lit.md could locate implies `A_5(delta)` for any `delta < 1/30`, or
  `(L_5^circ)`.
* Translation: `X' =` the (twisted, Gorenstein, singular) toric Fano 4-fold of the `A_4` root polytope. It is the small
  crepant contraction of `Lbar_5^tw`, which is only weak Fano. Then `X'(Q) = T_5(Q)`, `H_{-K_{X'}} = R''^2` and
  `dist ≍ diam/R''`. In these terms the target is: the essential approximation constant of the rational point `e`
  w.r.t. `-K_{X'}`, restricted to `eta_0`-balanced points in one unit class, exceeds `15/4`.
* Every fixed-section forced-divisibility argument is capped at exactly `15/4` on balanced points: this is row 20,
  restated through `construct.md` Thm B. Also `tau_e(-K) <= 4 - 2^{3-k}`, and exactly, at levels `m <= 5`,
  `max ord_e/m = 10/3` (Cilleruelo-Córdoba).
* The proved cases of McKinnon's conjecture concern `alpha`, not `alpha_ess`, and are split or homogeneous. Roth-type
  bounds satisfy `beta <= tau`. Counting results (Batyrev-Tschinkel apply: `rank Pic(X')^Gal = 15`,
  `~ c B (log B)^{14}`) work at the natural scale.
* The lit18 drifting-phase dismissal is void in the `T_5` formulation. The theorems apply and are too weak.

**Files.** `kill_lit_checks/`:
* `fan.py`, `fan.out`: ampleness, Picard rank, widths, 1-PS degrees;
* `facets.py`, `facets.out`: reflexive, non-simple root polytope;
* `reg2.py`, `reg2_k34.out`, `reg2_k5.py`, `reg2_k5.out`: ranks;
* `vander.py`, `vander.out`: exact Vandermonde order.
