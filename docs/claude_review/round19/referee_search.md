# Referee report on round19/search.md (adversarial)

Status: complete. All scripts and outputs are in `referee_search_checks/` (absolute path
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/referee_search_checks/`).

## 0. Verdict

**The result survives. The stop rule fires.**
* I re-derived every proof step: the excess identity, `N0 >= 32`, the special-`lambda` reduction, the provably
  complete discriminant candidate set and the `±1` fast path. All are correct for every `D`.
* I recomputed the complete `D = 4`, `D = 5` and `D = 6` searches with independent exact code: my own orbit
  representatives (checked by Burnside), my own lines, a different discriminant routine, and my own literal
  implementation of the Prop 2.6 formula.
  * The `lambda = ±1` histograms and the generic histograms agree with the author's **entry by entry** at every
    `D`.
  * The sets of genuinely special `(A, C, lambda)` are **identical**: 19, 22 and 31.
  * Every special `lambda` has gain exactly `-32`.
* **Maximum gain: `-32` at every `D <= 6`, for all `lambda in Q^*`.** At `D = 5` the second value is `-44`.
* No fatal or major issue. There are four minor issues (m1-m4) and one remark (m5), §6.
* Evidence is labelled as evidence. No growth or uniformity claim is made.

## 1. Object being computed (Prop 2.6 as corrected)

I checked the quantity in search.md §1 against `round18/theory.md` §2.4 and `round18/referee_theory.md` §2.6 / §3.4
(item F), line by line.

* **Pencils.** Units of `M_{0,8}` in Kapranov's chart (point 8 at infinity) are, up to constants, the products
  `prod (x_i - x_j)^{e_ij}` with `sum e_ij = 0`. Writing the positive and negative parts as `A` and `C` gives
  `f = lambda A/C` with `A, C` coprime products of `D` differences each. So "coprime pairs of `D`-multisets of edges of
  `K_7` with disjoint supports, times `lambda in Q^*`" is exactly the set of unit pencils of degree `D`. ✓
* **All rational `lambda`** (the item-F correction). Covered: `lambda = ±1` explicitly, and every other `lambda`
  through Prop 2.1 / Thm 2.3 (see §2). ✓
* **Multiplicities of every component of the three fibres.** The gain formula of theory.md §2.4 needs, for each
  visible `T`, the orders `nu_T(A), nu_T(B), nu_T(C)`, the order `mu_T`, and for every `Q`-irreducible interior
  factor `Y` of `B` its degree and all 98 orders `nu_T(Y)`. Only `B` can have interior factors, because `A` and `C`
  are units. search.md uses exactly this formula (via the author's `exact.py`, the round-18 referee's formula) and
  never a reduced version of it, except where Lemma 1.1 / Prop 2.1 / Prop 3.1 prove that the two agree.
  I re-derived those (§2). ✓
* **`beta`-degrees.** `H.beta = 31` and `E_T.beta = 1` for all 98 `T` with `|T| >= 3` (including the seven
  points, `|T| = 6`). Every boundary divisor (strict transform of `x_i = x_j`) has `beta`-degree `31 - 30 = 1`. An
  interior `Y` has `31 deg Y - sum_T nu_T(Y)`. `K.beta = -186 + 196 = 10`. All ✓, recomputed.
* **Exceptional divisors of the base-locus blow-up.** After the fixed part `sum mu_T E_T` is removed, the
  base locus of a unit pencil has codimension 2. The exceptional divisors of its blow-up get `beta~ = pi^* beta`
  degree 0 by the definition used in theory.md §2.4 ("exceptional degrees 0 on `Bal^pc`"). More to the point, the
  arithmetic proof of Prop 2.6 never uses them: on `Bal^pc`, `gcd(A(x), C(x)) = prod d_T^{mu_T}` exactly (referee
  round 18, point 3), so no accidental gcd primes appear. The quantity that Prop 2.6 actually needs is the
  displayed formula, and that is what the search computes. ✓ (This is a reading of the refereed Prop 2.6, not a new
  claim.)
* **The `lambda = p/q` normalisation.** `B = qC - pA` has the same divisor as `C - lambda A`, and constants change
  `rad` by `O(1)`. ✓

**Conclusion of §1: the search computes the right object.**

## 2. Re-derivation of the proof steps

### 2.1 Excess identity (1): correct

Here is a cleaner form of the author's proof that I checked independently. The fibre over 1 is
`F_1 = sum_e b_e D_e + sum_{T in Ts} (nu_T(B) - mu_T) E_T + sum_Y m_Y Y~`. Its class is `L`, so `F_1.beta = L.beta`.
The fibres over 0 and infinity are pure boundary. Their reduced supports are the visible `T` with
`nu_T(A) > mu_T` or `nu_T(C) > mu_T`, that is, with `a_T ≠ c_T`, each of `beta`-degree 1. So
`F_{0,red}.beta + F_{inf,red}.beta = N0`. A component is counted in at most one fibre: two of `A, B, C` sharing a
component would force the third to share it too. Moreover `nu_T(B) = mu_T` when `a_T ≠ c_T`. Together these give
`gain = (F_1.beta - F_{1,red}.beta) - N0 = -N0 + sum_{Y ⊂ F_1} (m_Y - 1) Y.beta`, which is (1.1). The author's
bookkeeping proof, including the count "30 sets of `Ts` contain a pair" (`2^5 - 2`), is also correct.

The `D <= 3` pencils (all 21020, both signs) were recomputed with my own implementation of the **literal** formula
(no use of (1.1)), and the D = 4, 5, 6 samples were recomputed with it as well (§4). All agree with (1.1).

### 2.2 N0 >= 32 (2): correct

`w(T' + i + j) - w(T' + i) - w(T' + j) + w(T') = w_ij`, because the cross terms `sum_{k in T'} w_ik` and
`sum_{k in T'} w_jk` cancel. Since `w = 0` on `|T| <= 1` and on `T = [7]` (`deg A = deg C`), every `T` with
`w(T) ≠ 0` is one of the 119 visible sets. The 32 quadruples `{T', T'+i, T'+j, T'+i+j}` partition `2^[7]`. ✓ The
bound is sharp (cross-ratio). Exact data agree: `min N0 = 32` at every `D <= 6` in my run (§4).

### 2.3 Special lambda reduction (3): correct

* **Prop 2.1(1).** If `x_i - x_j | B_lambda`, then `{i,j}` is not an edge of `A` or `C`. Restricting to
  `x_i = x_j` gives `C| = lambda A|`, two nonzero products of primitive linear forms whose (lex) leading
  coefficients are `±1`. Comparing leading coefficients gives `lambda = ±1`. ✓
* **Prop 2.1(2).** Same argument for the initial forms along `L_T`, which are products of forms `y_i - y_j`, `y_i`,
  `x_{i0} - x_k` and `x_k - x_l`, all with leading coefficient `±1`. ✓ Also `nu_T(B) >= mu_T` always, with
  `nu_T(B) > mu_T` iff `in_T(C) = lambda in_T(A)`, because both are homogeneous of the same `eps`-degree
  `mu_T = a_T = c_T`. ✓
* **Thm 2.3.** If `B_lambda = Y^2 G` with `deg Y >= 1`, then on any line `Y(sp+tq)` is either `≡ 0` or a nonzero
  binary form of degree `>= 1`. Either way `Disc_D(c - lambda a) = 0`, where `Disc_D` is the discriminant of formal
  degree `D` (so a root at infinity counts). ✓
  * *Implementation caveat.* The author computes the affine `t`-discriminant in `Z[t, lambda]`. This equals
    `Disc_D(c - lambda a)` as a polynomial in `lambda` only if the leading `t`-coefficient is a nonzero polynomial in
    `lambda` and the discriminant is the universal one (exact division by the leading coefficient). The author
    checks the first condition (`deg_drop_line1 = 0`), and the second holds for flint. It is in any case confirmed
    by my independent computation, which interpolates `Disc_D` from `2D - 1` exact univariate discriminants taken
    at `lambda` with nonvanishing leading coefficient, on lines of my own. My code finds **the same genuinely
    special sets** (§4).
* **The discard step.** "`q c - p a` squarefree as a binary form of formal degree `D` on some line implies `B`
  squarefree" follows by the same argument. The test "degree `>= D - 1` and `gcd(b, b') = 1`" is the correct
  criterion for a binary form of formal degree `D`. ✓
* **Consequence.** For `lambda ∉ {±1} ∪ Λ_ℓ`, `gain = -N0 <= -32`. Since `Λ_ℓ` is finite and explicit, the
  candidate list is **provably** complete, which improves on the probability-1 list of the round-18 referee. ✓

### 2.4 Fast path for lambda = ±1 (4): correct

Prop 3.1(a): the inside-`T` edges map injectively to the forms `y_i - y_j` (with `y_{i0} = 0`), and the outside
forms involve some `x_k`. So `in_T(C) = s in_T(A)` forces the inside edges of `A` and `C` to coincide, and since the
supports are disjoint they are empty. Hence `T in S` implies `a_T = c_T = 0`. ✓
Prop 3.1(b) is (1.1) with `b_e = m_Y = 1`. ✓

My implementation of this path is independent of the author's. It does not use the author's combinatorial initial
forms or hash. For each `T` with `a_T = c_T = 0`, it tests `C(P_T) = s A(P_T)` at a random integer point
`P_T in L_T`; this is necessary for `B|_{L_T} = 0`, so the filter loses nothing. It then computes `nu_T(B)` exactly
by symbolic substitution at every surviving `T`. Every filtered candidate in every run was a true member of `S`.

## 3. Enumeration completeness

* **Reduction to `A` up to `S_7`.** `S_7` permutes the 119 visible sets and fixes the class `beta`, so the formula is
  `S_7`-invariant. A relabelling sends `x_i - x_j` to `±(x_{σi} - x_{σj})`. That changes `lambda` by a sign, which
  is harmless because all `lambda in Q^*` are covered (both signs at `lambda = ±1`). ✓
* **Orbit representatives, recomputed independently** (`referee_search_checks/rorbits.py`). The canonical form is
  the minimum over all 5040 permutations of a base-7 code of the edge-count vector, and the orbit count is checked
  against Burnside's lemma (cycle index on edges):

  | `D` | multisets | orbits (mine = Burnside = author) | coprime pairs (mine = author) |
  | --- | --- | --- | --- |
  | 1 | 21 | 1 | 20 |
  | 2 | 231 | 3 | 590 |
  | 3 | 1771 | 8 | 9900 |
  | 4 | 10626 | 22 | 128460 |
  | 5 | 53130 | 60 | 1278340 |
  | 6 | 230230 | 173 | 11313931 |

* **The author's representatives are a complete system** (`chk_author_reps.py`). For `D = 4, 5, 6` the author's
  `reps_D*.json` hit 22 / 60 / 173 distinct orbits, the same set as mine. In particular the augmentation-built
  `D = 6` list is complete.
* **The author's chunks cover every `C`.** For every orbit, the chunk intervals in `out/D{4,5,6}.jsonl` are
  contiguous from 0, and the last one ends at `binom(21 - #supp A + D - 1, D)`. The totals are
  128460 / 1278340 / 11313931 pairs.
* **The special-`lambda` list** is complete by Thm 2.3 (proved; §2.3), given that `Delta_ℓ1 ≢ 0` for every pair,
  which the author logs and my run checks independently. My candidate sets come from different lines and a
  different discriminant routine (interpolation). They give the same genuinely special `(A, C, lambda)` (§4).

## 4. Independent recomputation

All of this code is my own, in `referee_search_checks/`, and none of it imports the author's code:
* `rcore.py`: the literal Prop 2.6 formula, flint factorisation over `Z` with a sympy fallback, exact `nu_T` by
  symbolic substitution, and the line tools;
* `rorbits.py`;
* `rsearch.py`: one chunk;
* `rtodo.py` and `rrun.sh`: the resumable driver, at most 2 processes, each under `prlimit --as=3500000000`.

The search uses its own random lines (seed 4242, coordinates in `[-1000, 1000]`, all `q_i` distinct, so every
difference keeps `t`-degree 1). It computes the candidate polynomial `Delta(lambda) = Disc_D(c - lambda a)` by
**exact interpolation** from `2D - 1` univariate discriminants, taken at integer `lambda` where the leading
coefficient does not vanish. That is a different route from the author's bivariate `fmpz_mpoly.discriminant`.

**`lambda = ±1`.** If `C - sA` restricted to some line is a squarefree form of degree `D`, the gain is
`-N0 + sum_{T in S} (nu_T(B) - 1)`. `S` is found by testing `C(P_T) = s A(P_T)` at random points `P_T in L_T` for
the `T` with `a_T = c_T = 0`, followed by exact `nu_T(B)`. Otherwise the gain comes from the literal formula.
A random fraction of the pencils (`VERIFY`) is recomputed with the literal formula.

**`lambda ≠ ±1`.** Every rational root of `Delta`. A root is discarded if `qC - pA` is squarefree on one of two
further lines; then the gain is `-N0`, and with `VERIFY > 0` it is also recomputed with the literal formula.
Otherwise the gain comes from the literal formula.

### 4.1 D <= 3 (every pencil through the literal formula)

| `D` | pairs | `±1` pencils, all re-done literally | `±1` histogram top | special `lambda ≠ ±1` | mismatches |
| --- | --- | --- | --- | --- | --- |
| 1 | 20 | 40 | `-32` (20), `-48` (20) | none | 0 |
| 2 | 590 | 1180 | `-32` (24), `-40` (12), `-48` (152), `-52` (144), `-56` (272) | 6: `-1/4` (5), `-4` (1); all `-32` | 0 |
| 3 | 9900 | 19800 | `-32` (20), `-40` (4), `-44` (2), `-48` (208), `-52` (188), `-56` (876) | 11: `±4/27` (5 each), `27/4` (1); all `-32` | 0 |

These are identical to round 18 (referee §3.3 and §3.4: 17 special values, all `-32`) and to search.md §4. All
222 + 407 squarefree candidates were confirmed literally at `-N0`.

### 4.2 D = 5, complete (all 60 orbits, 1,278,340 pairs; `out/R5.jsonl`, 60 chunks, 0 failures, 527 s CPU)

* `Delta ≢ 0` for every pair, and no `lambda = ±1` fallback was needed: every `C ∓ A` is squarefree on a line.
* **Generic value:** histogram top `-32` (10), `-44` (1), `-48` (167), `-52` (44), `-54` (14), `-56` (998). The
  full histogram is **identical** to the author's.
* **`lambda = ±1`:** 2,556,680 pencils, top `-32` (20), `-44` (2), `-48` (334), `-52` (88), `-54` (28),
  `-56` (1996). The full 50-value histogram is **identical** to the author's `hist_lam_pm1`. 24,223 randomly chosen
  pencils (about 1%) were re-done literally, with 0 mismatches. All 13,951 filtered `S`-candidates were true.
* **Special `lambda ≠ ±1`:** 842 candidates from my lines (the author had 848 and 840 from theirs; these counts
  depend on the line), of which 820 are squarefree. All 820 were confirmed literally at `-N0`.
  **22 genuinely special**, `±256/3125` (5 each), `±108/3125` (5 each), `3125/256` (1), `-3125/108` (1), **all
  with gain exactly `-32`.** As a set of `(A, C, lambda)` this is identical to the author's 22.
* **Maximum gain at `D = 5`: `-32`. The second value is `-44`.** Confirmed.

### 4.3 D = 4, complete (22 orbits, 128,460 pairs; `out/R4.jsonl`; 5% of the `±1` pencils re-done literally)

* `lambda = ±1`: 256,920 pencils, a 43-value histogram **identical** to the author's and to round 18:
  `-32` (24), `-40` (14), `-48` (318), `-50` (10), `-52` (238), `-54` (6), `-56` (1636). 12,288 literal
  re-checks, 0 mismatches.
* Generic histogram: identical to the author's; top `-32` (12 pairs).
* Special `lambda ≠ ±1`: 648 candidates on my lines (the author had 654), 629 squarefree, all confirmed
  literally. **19 genuinely special**, `-27/256` (10), `1/16` (5), `-256/27`, `16`, `-1/4`, `-4` (1 each),
  **all `-32`**. The same set of `(A, C, lambda)` as the author's.
* Three chunks first failed on the python-flint `factor()` overflow (m4). After the sympy fallback they reran
  cleanly. The failure raised an exception and could not have produced a wrong value.

### 4.4 D = 6, complete (173 orbits, 11,313,931 pairs; `out/R6.jsonl`, 334 chunks, 0 failures, 5151 s CPU)

* `Delta ≢ 0` for every pair, and no `±1` fallback was needed.
* `lambda = ±1`: 22,627,862 pencils, a 58-value histogram **identical** to the author's: `-32` (24), `-40` (16),
  `-44` (2), `-48` (508), `-50` (10), `-52` (440). 12,692 random literal re-checks, 0 mismatches. All 37,614
  filtered `S`-candidates were true.
* Generic histogram: identical to the author's; top `-32` (12), `-40` (8), `-44` (1).
* Special `lambda ≠ ±1`: 1067 candidates on my lines (the author had 1076), 1036 squarefree, all confirmed
  literally. **31 genuinely special**:
  * `-3125/46656` (10), `16/729` (10), `-1/64` (5);
  * their reciprocals `-46656/3125`, `729/16`, `-64` (1 each);
  * the 4-point values `-4/27`, `4/27`, `27/4` (1 each).

  **All 31 have gain `-32`.** The same set as the author's.
* This is an independent exact `D = 6` result. It closes the author's caveat that `D = 6` was not re-verified
  pencil by pencil.

### 4.5 Random literal samples and extreme cases

* **Random `lambda`.** I drew 300 random pairs each at `D = 5` and `D = 6`, each with a random `lambda = p/q`
  (`|p|, q <= 50`, `lambda ≠ ±1`), and evaluated them with the **literal** formula (`rsample.py`). The gain
  equals `-N0` in 600 / 600 cases. This tests Prop 2.1, Thm 2.3 and (1.1) directly.
* **Every claimed extreme case.** All 140 of the author's records (gain `> -40` and all special) were recomputed
  literally: 0 mismatches, all `-32` (§5).
* **`min N0 = 32` at every `D <= 6`** (generic histograms). This agrees with Prop 2.2.

## 5. Evidence vs proof

* **Labels are used correctly.**
  * "Exact" is reserved for finite exact computations made complete by proved reductions.
  * The forgetful-pull-back structure of the maximisers is labelled Data.
  * "Closing the abc lane" is called a decision rule, not a theorem.
  * The non-unit pencils of Prop 2.6 (base locus in codimension-2 boundary strata) and `D >= 7` are listed as
    not searched.
* **Proved for every `D`** (correct as stated): a unit pencil can beat `-32` only at a special `lambda` (`±1` or
  a rational root of `Delta_ℓ`), and only through a non-reduced fibre over 1. When `lambda ≠ ±1`, the non-reduced
  component must moreover be a repeated interior factor.
* **The "Data" claim about the maximisers is reproduced.** I recomputed all 140 of the author's recorded pencils
  (gain `> -40` plus every genuinely special `lambda`; 43 / 42 / 55 at `D = 4 / 5 / 6`) with my literal formula and
  found 0 mismatches. Every one has gain `-32`, and `A ∪ C` is supported either on 3 labels or on 4 labels with
  `A/C` invariant under `x -> 1/x` (sympy check, `chk_author_records.py`). Such a function factors through a
  forgetful map to `Mbar_{0,4}`, so the pencils are forgetful pull-backs. The gain of a pull-back `h ∘ chi_Q` is
  `16 deg h - 16 #h^{-1}{0,1,∞}` by the formula directly. This does **not** need Bruno–Mella (which is cited from
  the abstract only), so search.md's remark on morphism extremals is sound as data.
* **Overreach to correct (minor).** In §8, "the maximum is isolated and shows no drift upward with `D`" is a reading
  of six data points. It must not be taken as evidence about `D >= 7`. The gain can exceed `-32` only at special
  `lambda`, and the counts of special `lambda` grow with `D` (17, 19, 22, 31 at `D = 2+3, 4, 5, 6`). Nothing here
  bounds the gain at `D >= 7`.

## 6. Issues

No fatal or major issue was found. Every proof step re-derives, every number I recomputed agrees exactly, and the
enumeration is complete.

**Minor.**
* **m1 (Thm 2.3 implementation).** search.md §2 identifies the bivariate `fmpz_mpoly.discriminant` with `Disc_D`
  without arguing it. The identification needs the universal discriminant (exact division by the leading
  coefficient) and a leading `t`-coefficient that is not identically 0 in `lambda`. The second condition is logged,
  but the first is not argued. It is harmless: my interpolation of `Disc_D`, on other lines, gives the same genuinely
  special sets at `D = 4, 5, 6`. A sentence should be added.
* **m2 (overreach in §8).** "The maximum is isolated and shows no drift upward with `D`" reads six data points as a
  trend. Label it data with no bearing on `D >= 7`. The number of special `lambda` grows with `D`
  (19, 22, 31), and the possible excess at a special `lambda` is not bounded by anything proved here.
* **m3 (D = 6 provenance).** search.md is candid that its `D = 6` gains rest on the fast path, without
  pencil-by-pencil literal re-verification. My independent `D = 6` run (§4.4) closes this.
* **m4 (python-flint 0.9 bug, for the record).** `fmpz_mpoly.factor()` raises `OverflowError` in its sort key
  for some polynomials with large coefficients. I hit it on 3 `D = 4` squarefree-candidate verifications. It
  raises an exception and never returns a wrong factorisation, and the author's runs logged no failure. I added a
  sympy fallback, which checks that the product of the factors reproduces `B`.
* **m5 (remark, unrefereed).** `D` is the degree in the chart with point 8 at infinity, so it is not
  `S_8`-invariant. The gain `10 - K^log_rel.beta~` is intrinsic, because `beta` is `S_8`-symmetric. The result
  therefore holds for the `S_8`-saturation of the searched class, which is slightly more than stated. Nothing
  depends on this.

## 7. Corrected statement

**(Exact; complete by proved reductions; independently reproduced.)** Let `f = lambda A/C` be a unit pencil on
`Mbar_{0,8}`: `A, C` coprime products of `D` differences `x_i - x_j` (Kapranov chart, point 8 at infinity), and
`lambda in Q^*`. Let `gain(f)` be the theory.md Prop 2.6 quantity. Then:

1. **(Proved, every `D`.)**
   `gain = -N0(A,C) + sum_{Y ⊂ F_1} (m_Y - 1) Y.beta`, and `N0(A,C) >= 32`.
   * For `lambda ≠ ±1`, the fibre over 1 has no multiple boundary component. So `gain > -N0` only if `B_lambda` has
     a repeated interior factor, and every such `lambda` is a rational root of `Disc_D(C|_ℓ - lambda A|_ℓ)` for any
     line `ℓ` on which this discriminant is not identically 0.
   * Hence `gain <= -32` for every `lambda` outside an explicit finite set `{±1} ∪ Λ_ℓ(A,C)`.
2. **(Exact.)** For `D = 4, 5, 6` and every `lambda in Q^*`, `gain(f) <= -32`, with equality attained.
   * The genuinely special `lambda ≠ ±1` number 19, 22 and 31, all of gain `-32`.
   * At `D = 5` the next value is `-44`.
   * Together with round 18 (`D <= 3`, re-checked here with the literal formula on every pencil), this covers
     every `D <= 6`.
3. **(Data.)** Every pencil of gain `-32` with `D <= 6` that was recorded (all with gain `> -40`) is a forgetful
   pull-back `h ∘ chi_Q` of a genus-0 Belyi map. Nothing is claimed about `D >= 7`.
4. **The round-18 stop rule fires:** every unit pencil with `D = 5` (all `lambda`) and with `D = 4`
   (`lambda ≠ ±1`) has gain `<= -32`. This is a decision rule, not a theorem.
   * The abc-pencil question stays **open** for unit pencils with `D >= 7` and for the non-unit pencils of
     Prop 2.6 (base locus in codimension-2 boundary strata).
   * A positive answer would still give only `(L_8^pc)`, which is not enough for Theorem D without the unproved
     nested-chain extension.
5. **No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` is unchanged.
