# Referee report on round18/census.md

Referee scripts (independent; none of the author's code is imported): 
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/referee_census_checks/`
(`rlib.py`, `test_rlib.py`, `rcensus.c` with `rc5.out`, `rc678.out`, `rsym.py`, `remptiness.py`,
`rlattice.py`, `rk6fam.py`, `notes.md`). The only foreign code that was run is the repository's own
`threshold/L8/k6fam` generator, read in place, to produce the `k = 6` points. It was used only to
generate points; every contact was recomputed by the referee.

## 0. Verdict

**The result survives, in the weak form in which it is stated: outcome (c), weak evidence, no
growth or uniformity claimed.** Every formal proposition I re-derived is correct: Lemma 1.1, Lemma 1.2,
Prop. 1.3, Theorem E, Prop. 7.1, Cor. 7.2, Cor. 7.3, Prop. 7.4, Prop. 7.5(1)-(4) and Lemma 7.7. Every
computational claim I re-ran with independent code is reproduced exactly, or within sampling
(§4, §5).

**One major problem, which is interpretive.** In §6.3, at the end of §7.2 and in §8(a), the text
draws consequences from Prop. 7.5 that are false as stated:

* "no family of curves of genus `<= 1` and bounded degree [...] produces balanced points densely";
* "a counterexample [...] lives on no bounded-degree family of curves of genus `<= 1`";
* "a dense supply would need infinitely many rigid curves of class `∝ beta`, or sporadic points on
  no curve of genus `<= 1`".

Every `Q`-point of `M_{0,8}` lies on a `Q`-line of Kapranov's `P^5`, and these lines form a
bounded-degree covering family of rational curves. So the alternative "sporadic points on no
curve of genus `<= 1`" is the empty set. Prop. 7.5(4) gives `(L_8)` on each member only with
member-dependent constants `O_C(.)`. Reading it as a statement about the family treats a
curve-dependent threshold as uniform, which the owner's constraints forbid. What Prop. 7.5 does
prove is narrower: no bounded-degree covering family of genus-`<= 1` curves can have members that
each carry infinitely many balanced points with small `H_Sigma`. That rules out the `(L_5)` and
`(L_6)` mechanisms and Thm 2(a)-(c) of `round4/L7.md` with bounded degree. It does not rule out a
dense supply made of finitely many points per member.

There are also several minor issues (§7). The census studies the right object: §1 below.

## 1. Is the object exactly (L_8) of ptolemy.md §7 as used in Steps 0-5?

Yes. I checked it item by item against `ptolemy.md` §1, §6.2 and §7.

* **Anchoring.** `P_0 = (1,0)` (`P_0 = 1` in §1), `Y_j = det(P_0,P_j)`, and residues
  `t_ij = det/G_ij`. The census uses the same. Step 2 needs the configuration of the 8 circle
  points to lie outside `Z`. This is `Phi = [P_0:...:P_7]` (Thm B(i)), and it is what the census
  quantifies over.
* **Balance.** Step 1 delivers `|log n_T - w| <= eta w` for **every nonempty** `T ⊆ [7]`. The
  census's §1.1 statement imposes it on all blocks, and then argues that the 8 singletons and the
  top block can always be filled freely. This argument is correct:
  - singletons enter no `G_ij`;
  - the top block can be supplied by `P_j -> (n X_j, 1)`, which changes neither `Phi` nor any
    `t_ij` or `Y_j`, and adds a one-block chain at a new prime.

  So only the 119 visible blocks matter. `eta_0` is the balance constant of `(L_8)`, and Step 4
  uses it through `eta = 2 sqrt2 · 255/sqrt(M-7) -> 0`. That is consistent.
* **Nested chains.** "Two nested chains on disjoint row sets at every prime" (item 532) is used
  verbatim. The consequence "at most `m - 2 = k - 3 = 5` visible blocks per prime" is correct:
  - a chain of visible sets inside `L` has at most `|L| - 1` members;
  - if `U = ∅` the chain stops at size `m - 1`;
  - a one-element `U` contributes no visible block.
* **`H_Sigma`** is `sum_j log|Y_j| + sum_{i<j} log|t_ij|`, with 7 + 21 terms. Same in the census.
* **`Z`** is one fixed proper closed subset of `Mbar_{0,8}`, independent of `w` and of the data.
  Same in the census. Prop. 1.3(2) has a small labelling issue: see §7, item m1.
* **Scale.** `w = W/2^m = W/128` (`ptolemy.md` §1, Step 0). The census's remark is correct.
* **Typo finding.** The finding is correct. The visible `T` are `T ⊆ [m]` with `2 <= |T| <= m-1`, so there are
  `2^m - 1 - m - 1 = 2^m - m - 2` of them (119 at `m = 7`, which is the number of boundary divisors
  `2^(k-1) - k - 1`). `ptolemy.md` §1 prints `2^m - m - 1`. The error is harmless.

## 2. Prop. 1.3 and the reduction to contact profiles

* **Lemma 1.1** is correct. It is `ptolemy.md` Lemma 5.2 with anchor at infinity, plus a
  normalisation.
  - One cosmetic fix is needed. Dividing the `X_j` by `g = gcd(X_i - X_j)` does not in general
    keep them integral. A translation `X_j -> (X_j - X_1)/g` does, and it is still a Möbius map
    fixing `infinity`.
  - `gcd_{i<j}(X_i - X_j) = prod p^{delta_p}` holds because, at every `p`, some pair `i, j` lies
    in no positive visible `T`. The maximal positive `T` are pairwise disjoint proper subsets of
    `[7]`.
  - **Residue-class rule (independently checked).** Three independent implementations agree on
    150 random configurations with `k = 5..8`, in 14,618 comparisons with 0 mismatches
    (`test_rlib.py`):
    - the `ptolemy.md` Thm B(iv) min-formula, by factoring;
    - the residue-class rule, in the unanchored form where a split counts `#{t : S is a class}`
      plus `#{t : S^c is a class}`;
    - the anchored gcd inversion `g_T = delta · prod_{T' ⊇ T} d_{T'}`.

    The gcd identity holds for the same reason as the gcd claim above.
* **Lemma 1.2** is correct for general `(L_8)` data, not only for circle data. The proof of Thm B(iv)(2) uses only
  - `g_ij = v_p(t_ij) + sum_{T ∋ i,j} v_p(n_T)`, with `g_0j = v_p(Y_j)`;
  - the compatibility of the blocks at `p`, which two chains on disjoint row sets guarantee;
  - the fact that the min-formula is 1-Lipschitz.
* **Prop. 1.3(1)** is correct with `eta_0 = eta_1/2` and `c = eta_1/4`. If `H <= cw`, then
  `|ell_T - w| <= 2cw + eta_0 w = eta_1 w`, so `w < w_1`. Hence `H >= cw - cw_1` always.
* **Prop. 1.3(2)** is correct for admissible `Phi`.

  Note that the summary's phrase "`(L_8)` is, up to constants, a statement about the 119 finite
  boundary heights" is a two-sided **sandwich**, not an equivalence:
  - `(L_8^c)` over **all** `Phi` is sufficient for `(L_8)`;
  - `(L_8^c)` over **admissible** `Phi` is necessary.

  The gap is real. A datum with small `H_Sigma` can have non-two-chain trees at the primes that
  divide residues. The body states (1) and (2) correctly.

## 3. Theorem E (emptiness)

The proof is correct.
* Each balanced visible block is `> 1`, and each prime lies in at most `k - 3` of them.
* So if `n_(j) < p_r` with `r = ceil(j/(k-3))`, the `j` smallest blocks are covered by `r - 1`
  primes, which forces `j <= (k-3)(r-1) < j`.
* Part (3) uses the `2^(k-3) - 1` visible `T ⊇ {i,j}`, `i, j >= 1`. It is a statement about the
  21 minors `det(P_i,P_j)` with `i, j >= 1`, not about the `Y_j`, which may be 1. The summary says
  "all 21 minors", which is right.

`remptiness.py` reproduces every entry of the table, including:

| `k` | min # primes `r_k` | `p_{r_k}` | class empty for `e^w <` (`eta = 0.1`) | minor `>=` |
|---|---|---|---|---|
| 5 | 5 | 11 | 8.8 | `10^2.6` |
| 6 | 9 | 23 | 17.3 | `10^7.8` |
| 7 | 14 | 43 | 30.5 | `10^20.0` |
| 8 | 24 | 89 | 59.2 | `10^49.4` |

At `eta = 1/2`, `k = 8`: `e^w >= 19.9` and minors `>= 10^20.1`. The sum bound is always weaker than
the max bound, as the census implicitly uses.

## 4. Calibration data, re-run independently

* **Exhaustive census** (`rcensus.c`, C, residue-class rule, own enumeration). Configuration
  counts and coverage histograms are reproduced **exactly** for:
  - `k = 5`, `B = 100`: 3,646,069 configurations, full coverage 775,434;
  - `k = 6`, `B = 60`: 5,302,911 configurations, maximum coverage 18/25, attained 198 times;
  - `k = 7`, `B = 40`: 3,797,877 configurations, maximum 23/56, attained 42 times;
  - `k = 8`, `B = 34`: 5,359,838 configurations, maximum 29/119, attained once, by
    `(0,4,11,13,21,23,30,34)`.

  The best `k = 5` points are also reproduced: `eps = 0.1333` at `w* = 2.767`, for
  `(0,22,35,46,85)`, with contacts `11,11,13,13,17,20,20,21,21,23`; and `eps = 0.1345` with
  `w* >= 3`, for `(0,17,29,48,99)`. At `k = 5` every tree is admissible, because a prime lies in
  at most 2 visible blocks. So `H_Sigma = 0` there is legitimate.
* **`k = 6` family** (`rk6fam.py`). I used the repository's generator unchanged; the census's copy
  is byte-identical to `threshold/L8/k6fam`. The contact integers were recomputed with my own
  gcd inversion. `w*` and `eps` match the census table **exactly** for `N = 2..14`:
  `eps = 0.3880, 0.2107, 0.1490, 0.1189, 0.0951, 0.0784` at `N = 4, 6, ..., 14`.
  - `eps` is not monotone (`N = 7` gives `0.2355` and `N = 13` gives `0.1031`); only the trend
    decreases.
  - The `max log|X|` values differ by the trivial normalisation `delta = 2`.
* **Lattice upper bounds** (`rlattice.py`). Lattice basis by my own successive sublattice
  intersection, index checked exactly against `prod n_T^{|T|-1}`, `flint` LLL, and the same
  `[-2,2]^{m-1}` search. Results:
  - `k = 5`, `A = 1000`: `H/w = 8.5, 9.1, 9.5`;
  - `k = 8`, `A = 20000`: `H/w = 465.2` and `467.9`, with exact divisibility verified.

  This is consistent with the census (`464 w`). The values are genuine **upper** bounds for the
  census minimum, and they are correctly not presented as evidence.

## 5. Section 7, re-derived

* **Prop. 7.1: correct.** The closure in Knudsen's model meets `𝒟_T` only over finitely many
  primes, so its structure sheaf is killed by some `N_Y`. Pulling back along the section gives
  `N_Y ∈ (s^* f_T)`, and `l_T(p) = (Phi.D_T)_p` (Thm B(iv)(1)).
* **Cor. 7.2: correct.**
  - Incompatible boundary divisors are disjoint over `Z`.
  - A `sigma`-fixed point has a `sigma`-fixed integral section.
  - `rsym.py` (my own code, all of `S_k`) gives minima **3, 6, 25, 48** at `k = 5, 6, 7, 8`. At
    `k = 8` the minimum is attained exactly by cycle type `2^4`; the next type, a transposition,
    gives 62.
* **Cor. 7.3: correct.** The counts are `(2^(q-1) - q - 1) · 2^(8-q)`, which gives 48, 80, 100, 112
  for `q = 4..7`. This checks the census's numbers by formula.
* **Prop. 7.4: correct as stated**, under its hypotheses `eta_n -> 0` and `H = o(w)`. On a curve of
  genus `<= 1`, a degree-0 class has height `O(sqrt h_A + 1)`. Lemma 6.2 holds for general `(L_8)`
  data.
* **Prop. 7.5: correct.**
  - (1) `7K = -2D_2 + D_3 + 2D_4`, with 28, 56 and 35 divisors, so `K.beta = 10`. The bound
    `-56(1+eta) + 126(1-eta) = 70 - 182 eta` holds, and `eta >= 5/13` follows.
  - (2) This is the standard normal-sheaf argument: the normal sheaf is generically generated
    modulo torsion, so `deg N >= 0`.
  - (3) This is a finiteness argument for the Chow variety.
  - (4) This is Step 3 plus the curve's own height theory. **The constants depend on `C`.**
* **Cor. 7.6: correct for the criteria of `round4/L7.md` Thm 2.**
  - Criterion (b) needs `beta` in the cone of ruling classes. Each ruling class has
    `-K.gamma >= 2` (free curves), and `K.beta = 10 > 0`, so criterion (b) never applies at
    `k = 8`.
  - Criteria (a) and (c) need classes near the `beta`-ray. With bounded degree those classes lie
    in `Z_d` by 7.5(3).
  - One wording fix: "lie in no covering family, i.e. infinitely many rigid curves" should read
    "lie in no covering family". Not covering does not mean rigid; such curves could move in
    divisor-sweeping families.
* **Lemma 7.7: correct.**
  - `G'_ab = G_ab`, because `T -> T \ {j}` is 2-to-1 onto the blocks containing `a, b`.
  - The residues are unchanged, so `H' <= H`.
  - The merged blocks are balanced at scale `2w` with the same `eta`.
  - The chains stay nested on disjoint supports.

  Hence `(L_7)` with `(c, eta_0, Z_7)` implies `(L_8)` with `(2c, eta_0, pi_j^{-1}Z_7)`.

## 6. Constraints of the project owner

* **Fixed-`S` Subspace exceptional sets.** Not used anywhere.
* **Gaussian gcd or denominator costs.** Not discarded. The census works with the form-blind
  `(L_8)`, whose sufficiency (Theorem D) already includes the item-532 extraction.
* **Drifting phase.** Not applicable.
* **Template-dependent threshold treated as uniform. Violated in the interpretive text.**
  - Prop. 7.5(4) gives `(L_8)` on each fixed member `C`, with error terms `O_C(sqrt(w + H) + 1)`
    that depend on the member.
  - From this, §6.3 concludes that no bounded-degree family "produces balanced points densely",
    and §8(a) concludes that a dense supply needs rigid curves or "sporadic points on no curve of
    genus `<= 1`".
  - These conclusions would need the constants to be uniform over the family. They are not.

  This is issue M1 below. The formal propositions themselves comply with the constraints.
* **Evidence vs proof.** The Data and Heuristic labels are used consistently. The §7.3 count over
  `E_t` is labelled heuristic. Its "`O(sqrt W)` points" assumes rank 1, and the census only says
  "positive rank". The sum still converges for any rank (a theta series), so the heuristic
  conclusion is unaffected.

## 7. Issues

**M1 (major, interpretive, wrong as stated).** The following sentences are wrong:
* §6.3, bullet 3: "No family of curves of genus `<= 1` and bounded degree, and no product chart,
  produces balanced points densely";
* §7.2, final paragraph: "a counterexample [...] lives [...] on no bounded-degree family of
  curves of genus `<= 1`";
* §8(a): "A dense supply would need infinitely many rigid curves of class `∝ beta`, or sporadic
  points on no curve of genus `<= 1`".

*Counterexample to the dichotomy.* `Mbar_{0,8}` is the iterated blow-up of `P^5`, and `M_{0,8}` is
the complement of the 21 hyperplanes `<p_A>`, `|A| = 5`. Take a `Q`-point `x`. A general `Q`-line
through `x` avoids the centres, which have codimension `>= 2`. Its strict transform is a smooth
rational curve of class `h`, with `h.H_A = 1` and `h.E_B = 0` (class imbalance `eta = 1`). These
lines form a covering family of rational curves of bounded degree, defined over `Q`, and **every**
rational point of `M_{0,8}` lies on members of it.

Likewise, `Mbar_{0,8}` is birational to `(P^1)^5`, so every set of rational points off a closed set
is the image of rational points of a product chart. "No product chart produces balanced points
densely" is therefore, up to a closed set, the non-existence of a dense balanced supply, which is
essentially `(L_8)` itself.

What is proved is the following. For each fixed member `C` with `K.C <= 0`, only finitely many
points have `eta <= eta_0 < 5/141` and `H_Sigma <= c w` with `c < (10 - 282 eta_0)/60`. That is
Prop. 7.5(4) plus Northcott. The bound on `w` depends on `C`. So a refutation made of finitely many
balanced points on each of infinitely many members of a covering bounded-degree family, for
instance the lines, is **not excluded**. The correct reading of Prop. 7.5 is: "mechanisms in which
each member of a bounded-degree covering family of genus-`<= 1` curves carries infinitely many
balanced small-`H` points, and criteria (a)-(c) of `round4/L7.md` Thm 2 with bounded degree, cannot
refute `(L_8)`".

The §0 summary sentence about "mechanisms" is acceptable under this reading. The §6.3, §7.2 and
§8(a) sentences are not.

**m1 (minor).** Prop. 1.3(2) says "admissible for some anchor and labelling". Changing which label
is the anchor replaces `Phi` by `sigma Phi`. So `(L_8)` then gives the conclusion only off
`sigma^{-1} Z`. Fix: either use the anchor 0 only, or replace `Z` by `∪_sigma sigma^{-1} Z`, which is
still proper and closed.

**m2 (minor).** "`(L_8)` is equivalent, up to constants, to a statement about finite boundary
heights" overstates the result. It is a sandwich: all `Phi` is sufficient, admissible `Phi` is
necessary (§2).

**m3 (minor).** §7.2 / §8(b) say "`Z` ⊇ every curve of genus `<= 1` and class `d beta` with
infinitely many rational points (Prop. 7.4)". Prop. 7.4 proves the converse: a curve carrying
balanced points has class `∝ beta`. What `Z` must contain is every curve carrying infinitely many
balanced small-`H` points. For genus 0 curves with a `Q`-parametrisation, of class `d beta`, that
satisfy (i)-(iv), this is `round4/L7.md` Thm 2(a), which is proved. For genus 1 it is Thm 2(c),
which is only sketched. The citation should be to those results, with their hypotheses.

**m4 (minor).** §4.2: "bounded cost `<= 2·79`". Moving the `Sigma`-parts `s_T` out of the blocks
multiplies each `t_ij` by `prod_{T ∋ i,j} s_T`. That costs `sum_T C(|T|,2) log s_T <= 6 · 79` at
`k = 6`, not `2 · 79`. The cost is still bounded, so the conclusion is unchanged.

**m5 (minor).** §4.1 says "at `k = 5` balanced bounded-`H` points have positive density among small
integer configurations". The data show 21% **full coverage**, not 21% balanced. The best `eps` is
0.133, and only at `w* <= 3.3`. The CRT comparison rows (`w = 3.9-8.5`) lie at scales where the
exhaustive data show no balanced point. There, abundance rests on Prop. C and the count, not on
the data shown.

**m6 (minor).** Lemma 1.1: the normalisation needs a translation (§2).

## 8. Corrected statement

The census's study of `(L_8)` (ptolemy.md §7, normalisations as in Steps 0-5) proves:
* **(i) Sandwich reduction (Prop. 1.3).**
  - `(L_8^c)(eta_1, Z)` for all `Phi` implies `(L_8)` with `eta_0 = eta_1/2` and `c = eta_1/4`.
  - `(L_8)` implies `(L_8^c)(eta_0, Z)` for `Phi` admissible from anchor 0. With a symmetrised
    `Z`, it holds for any anchor.
* **(ii) Theorem E.** Balanced data need at least 24 primes. The hypothesis class is empty for
  `e^w < 59.2` (`eta = 0.1`), and all 21 minors `det(P_i,P_j)`, `i, j >= 1`, are `>= 10^49.4`.
* **(iii) `(L_8)` holds on special loci, with locus-dependent constants:**
  - every subvariety missing some `D_T`, including all symmetric loci (at least 48 trivial
    contacts), all forgetful fibres (48/80/100/112 bounded contacts) and the boundary;
  - each fixed curve of genus `<= 1` with `K.C <= 0`, for `eta_0 < 5/141`.
* **(iv) Families.**
  - Covering families of genus-`<= 1` curves have class imbalance `>= 5/13`. Hence no
    bounded-degree covering family can have members carrying infinitely many balanced small-`H`
    points.
  - Criterion (b) of `round4/L7` Thm 2 never applies at `k = 8`. Criteria (a) and (c) require
    unbounded degree.
  - This does **not** exclude a Zariski-dense balanced supply with finitely many points per member
    of a covering family. For example, every rational point lies on a `Q`-line of Kapranov's
    `P^5`.
* **(v) Lemma 7.7.** `(L_k)` implies `(L_{k+1})`.
* **(vi) Exact data, reproduced independently.**
  - Box censuses at `k = 5..8`.
  - The `k = 6` family: `eps -> 0` with `H_Sigma = O(1)`.
  - Lattice upper bounds `H <= ~465 w` at `k = 8`.
  - Block-sampling searches return `H ≍ w` also at `k = 5, 6`, so they cannot test `(L_8)`.

No counterexample to `(L_8)`, no candidate `Z` and no case of `(L_8)` on a Zariski-dense set was
produced. The evidence for `(L_8)` is weak (outcome (c)). No growth or uniformity is claimed.
