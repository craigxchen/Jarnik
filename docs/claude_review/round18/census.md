# Round 18: exact census of small-height balanced anchored configurations (k = 6, 7, 8)

Round-18 task 1 of the round-17 judge (§7.3): census of small-height balanced anchored
8-configurations, with `k = 6, 7` as calibration. Scripts: `census_checks/` (§9).

## 0. Summary

**Verdict: outcome (c), nothing small exists, and the evidence for `(L_8)` is weak.** In
addition, the census shows exactly why it cannot be strong, and it rules out every explicit candidate
locus for `Z` that was proposed.

1. **What is being counted (§1, proved).** By the anchoring lemma (`ptolemy.md` Lemma 5.2;
   `round4/L7.md` Thm 1 at `k = 7`), every point of `M_{0,k}(Q)` has an exact Boolean
   representation with `Y_j = ±1`, `t_ij = ±1` and blocks equal to its contact integers
   `d_T = prod p^{l_T(p)}`. With the Lipschitz bound `|log n_T - log d_T| <= 2 H_Sigma`, this makes
   `(L_8)` equivalent, up to its unspecified constants, to a statement about the 119 finite boundary
   heights `ell_T = log d_T` of points `Phi` (Prop. 1.3). So **minimum `H_Sigma` bounded `<=>`
   intrinsic imbalance `eps(Phi) -> 0`**, and the census is a census of contact profiles. This is
   the formulation used throughout; it was checked against `ptolemy.md` §7 and Steps 0-5.
2. **Theorem E (emptiness, proved).** A prime divides at most `k - 3` visible blocks (two nested
   chains), so a balanced configuration needs at least `ceil(V_k/(k-3))` primes: 9, 14, 24 at
   `k = 6, 7, 8`. At `k = 8` some block is `>= 89`; the hypothesis class of `(L_8)` is **empty for
   `e^w < 59`** at `eta = 0.1` (`< 20` at `eta = 1/2`), and every realisation has all 21 minors
   `>= 10^49` (`>= 10^20` at `eta = 1/2`). A literal small-height census at `k = 8` is empty by
   theorem.
3. **Calibration (Data).** `k = 5`: an exhaustive census of all 3.6 million configurations of
   height `<= 100` finds 21% touching every boundary divisor, with `eps` down to `0.133`, so
   `min H_Sigma = 0`. `k = 6`: the `(L_6)` family `E_{t_0}` of `L8.md` Thm B is reproduced with
   independent code: `H_Sigma = 0` (bounded after the fixed bad primes), and `eps` decreases from
   `0.39` to `0.078` as `w` grows from 32 to 494. Its points have `|X_i| ≈ 10^99` to `10^1490`.
4. **Blindness of search (Data plus count).** The natural CRT/lattice census over prescribed
   blocks gives `H_Sigma ≈ 8w, 42w, 150w, 464w` at `k = 5, 6, 7, 8`, linear in `w` for **every**
   `k`, including `k = 5, 6` where `H_Sigma = 0` points are abundant at the same scales. The reason
   is that a sampled block vector supports a bounded-`H` point with probability `e^{-5w}`,
   `e^{-17w}`, `e^{-49w}`, `e^{-129w}`. The `k = 8` numbers are exact **upper bounds** for the
   census minimum (e.g. a balanced 8-configuration at `w = 9.9` with `H_Sigma <= 4591`), not
   evidence. Exhaustive box searches at `k = 6, 7, 8` (heights `<= 60, 40, 34`) reach at most
   18/25, 23/56, 29/119 of the boundary divisors.
5. **Candidate loci for `Z` (§7, proved).**
   * Every locus missing a boundary divisor carries no balanced points (Prop. 7.1). This covers
     **all symmetric loci** (at least 48 of 119 contacts are identically 1; Cor. 7.2), **all
     forgetful fibres** (any configuration containing a fixed sub-configuration; Cor. 7.3) and
     boundary strata.
   * On a curve of genus `<= 1`, balanced bounded-`H` points force the class to be proportional to
     `beta`, so `K.C > 0` (Prop. 7.4).
   * **Families threshold (Prop. 7.5):** `7K.gamma >= a(70 - 182 eta(gamma))`. Covering families
     of genus-`<= 1` curves have `K.gamma <= 0`, so their classes have imbalance
     `eta >= 5/13`. Consequently `beta` is not in the cone of such classes, and **none of the
     mechanisms that refuted `(L_5)` and `(L_6)`, nor the product-chart criterion (b) of
     `round4/L7.md` Thm 2, can work at `k = 8` with families of bounded degree.** On each member
     of a covering family, `(L_8)` holds for `eta_0 < 5/141`.
   * `(L_k)` implies `(L_{k+1})` (Lemma 7.7). Balanced 8-points lie over balanced 7- and
     6-points. Over the `(L_6)` curves `E_t` the expected number of balanced 8-points is finite
     (heuristic, §7.3).
6. **Decision among (a)/(b)/(c) (§8).** (a) is not supported: no dense supply exists among the
   searched or excluded classes. For (b) there is no positive candidate: the only shape `Z` can
   have that the analysis produces is "curves of genus `<= 1` and class `d·beta`" (none known; their
   expected dimension is negative). (c) holds in the sense of the data, but because of item 4 it
   is weak evidence. No growth or uniformity claim is made.

## 1. The object being counted (formulation check against ptolemy.md §7 and §6.2)

Conventions as in `L8.md`: *Theorem/Proposition/Lemma* = proved here with all quantifiers (or cited
with the place of proof); *Data* = exact finite computation (Python integers, `flint`); *Heuristic*
= parameter or volume count, never used as proof. Scripts are in `census_checks/` (§9).

### 1.1 The statement tested

`(L_8)` is used verbatim from `ptolemy.md` §7 (= `L8.md` §1, = round-17 §7.3):

> There are `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z ⊂ Mbar_{0,8}` such that for every
> datum `D = (P_0 = (1,0), P_1, ..., P_7; (n_T)_{∅≠T⊆[7]}; (t_ij))` of primitive integer vectors
> with `det(P_i,P_j) = t_ij prod_{T∋i,j} n_T` (`1 <= i < j <= 7`, `t_ij ∈ Z\{0}`),
> `Y_j = det(P_0,P_j) ≠ 0`, blocks lying in two nested chains on disjoint row sets at every prime
> (item 532), `|log n_T - w| <= eta_0 w`, and `Phi(D) = [P_0 : ... : P_7] ∈ M_{0,8}(Q) \ Z`:
> `H_Sigma(D) := sum_j log|Y_j| + sum_{i<j} log|t_ij| >= c w - O(1)`.

Normalisations, checked against `ptolemy.md` §1 and Steps 0-5 of §6.2:

* `k = 8`, `m = 7`, anchor label `0`, rows `1..7`; `w` is the scale (`w = W/128` in Step 0;
  in `(L_8)` it is any real with the balance property).
* **Visible blocks** are `T ⊆ [7]` with `2 <= |T| <= 6`: there are `2^7 - 7 - 2 = 119`, one for each
  boundary divisor `D_T` of `Mbar_{0,8}` (the split `T | T^c ∪ {0}`). (`ptolemy.md` §1 prints the
  count as `2^m - m - 1`; its own table in §4 has the correct `2^(k-1) - k - 1 = 2^m - m - 2`. Typo
  only; nothing depends on it.)
* **Invisible blocks**: singletons `{i}` never occur in a minor `det(P_i,P_j)`; the top block `[7]`
  occurs in every minor with `i, j >= 1` and in no Plücker relation. Both can be prescribed freely
  for a given point `Phi` (top: `P_j -> (n X_j, Y_j)` with a prime `n ≈ e^w` dividing no `Y_j`,
  a Möbius rescaling that does not move `Phi`; singletons: any new primes `≈ e^w`), exactly as in
  the proof of `L8.md` Thm B. So only the 119 visible blocks matter, and "balanced" below refers to
  them.
* **Nested chains.** At a prime `p`, the visible blocks divisible by `p` are `T_1 ⊂ T_2 ⊂ ...` on
  a row set `L` and `S_1 ⊂ S_2 ⊂ ...` on a disjoint row set `U` (item 532: `L = {a_i(p) < a_0(p)}`,
  `U = {a_i(p) > a_0(p)}`). A chain of visible sets inside `L` has at most `|L| - 1` members, so
  **at most `|L| + |U| - 2 <= m - 1 - 1 = k - 3 = 5` visible blocks contain a given prime.** Call a
  laminar family of visible sets *admissible* if it is such a union of two chains.
* `H_Sigma` counts the 7 anchor residues `Y_j` and the 21 residues `t_ij`, each with `|.| >= 1`.

### 1.2 Reduction: residues vanish, the census object is the contact profile

For `Phi ∈ M_{0,k}(Q)` and a visible `T` let `l_T(p)` be the `p`-adic tree edge length
(`ptolemy.md` Thm B(iv)) and `d_T(Phi) = prod_p p^{l_T(p)}` the **contact integer**;
`ell_T(Phi) = log d_T(Phi)` is the finite part of the boundary height `h_{D_T}(Phi)`.

**Lemma 1.1 (exact Boolean form; `ptolemy.md` Lemma 5.2, `round4/L7.md` Thm 1 for `k = 7`).**
Every `Phi ∈ M_{0,k}(Q)` is represented by integers `X_1, ..., X_m` (anchor at `infinity`,
`P_i = (X_i, 1)`, `gcd_{i<j}(X_i - X_j) = 1`), and then for all `1 <= i < j <= m`

```text
|X_i - X_j| = prod_{T visible, T ∋ i,j} d_T(Phi),     Y_j = 1,   t_ij = ±1,   top block 1.
```

Moreover `l_T(p) = #{t >= 1 : T is exactly one residue class of {X_i mod p^t}}`.

*Proof.* Lemma 5.2 with anchor `P_0` gives `P_j = (D x_j/|y_j|, ±1)` with
`v_p(det(P_i,P_j)) = delta_p + sum_{T∋i,j} l_T(p)`. Replacing `P_j` by `(±X_j, 1)` changes nothing,
and dividing all `X_j` by `prod p^{delta_p} = gcd(X_i - X_j)` is the Möbius map `x -> x/delta`,
which fixes `Phi`. The class description: the `p`-adic Bruhat-Tits tree spanned by the ends
`X_1..X_m, infinity` has, at depth `t >= 1` below the standard vertex, one vertex for each residue
class mod `p^t` that contains some `X_i`; the edge entering it separates that class from the rest,
so it lies on the edge of `T` iff the class is `T`. `[]`

`census_lib.py` implements this; `test_lib.py` compares it with the Gromov-product min-formula of
Thm B(iv) on 300 random configurations (`k = 5..8`): 284,752 comparisons, 0 mismatches (Data).

**Lemma 1.2 (two-sided transfer).** Let `D` be any datum (as in §1.1) representing `Phi`. Then for
every visible `T`, `|log n_T - ell_T(Phi)| <= 2 H_Sigma(D)`.

*Proof.* `ptolemy.md` Thm B(iv)(2) (refereed, `referee_ptolemy.md` §1): at each `p`,
`|l_T(p) - v_p(n_T)| <= 2 max_{i<j} v_p(t_ij)` (with `t_0j = Y_j`), the blocks at `p` being
compatible since they form two chains. Multiply by `log p` and sum, using
`max_{i<j} v_p(t_ij) <= sum_{i<j} v_p(t_ij)`. `[]`

Put `ell_min = min_T ell_T`, `ell_max = max_T ell_T`, and the **intrinsic imbalance**
`eps(Phi) = (ell_max - ell_min)/(ell_max + ell_min)` (`= 1` if some `d_T = 1`), attained at
`w*(Phi) = (ell_max + ell_min)/2`. So `Phi` has an `eta`-balanced contact profile iff
`eps(Phi) <= eta`.

**Proposition 1.3 (what the census must look for).** Consider

> `(L_8^c)(eta_1, Z)`: there is `w_1` such that no `Phi ∈ M_{0,8}(Q) \ Z` has
> `|ell_T(Phi) - w| <= eta_1 w` for all 119 visible `T` and some `w >= w_1`.

1. If `(L_8^c)(eta_1, Z)` holds, then `(L_8)` holds with the same `Z`, `eta_0 = eta_1/2`,
   `c = eta_1/4`.
2. If `(L_8)` holds with `(c, eta_0, Z)`, then `(L_8^c)(eta_0, Z)` holds for all **admissible**
   `Phi` (two-chain trees at every prime for some anchor and labelling); in particular for every
   `Phi` whose 119 contact integers are pairwise coprime.

*Proof.* (1) Let `D` satisfy the hypotheses of `(L_8)` with `eta_0 = eta_1/2`. If
`H_Sigma(D) <= c w` then by Lemma 1.2 `|ell_T - w| <= eta_0 w + 2cw = eta_1 w` for all `T`, so
`w < w_1`. Hence always `H_Sigma(D) >= c w - c w_1`. (2) Lemma 1.1 represents an admissible `Phi`
with `H_Sigma = 0` and `n_T = d_T`; the invisible blocks are filled as in §1.1. `(L_8)` then gives
`0 >= c w - O(1)`. Pairwise coprime contact integers give one edge per prime, which is a chain. `[]`

So, up to constants that `(L_8)` does not fix, **a counterexample to `(L_8)` is the same thing as a
Zariski-dense supply of (admissible) points of `M_{0,8}(Q)` whose 119 finite boundary heights are
`(1 ± o(1)) w` with `w -> infinity`; residues never have to be controlled.** The census object is
therefore the contact profile `(ell_T(Phi))_T`, a function of the point `Phi` alone. The quantity
asked for, "minimum of `H_Sigma` at scale `w`", satisfies by Lemma 1.2

```text
max(0, max_T |ell_T(Phi) - w| - eta w)/2  <=  min_{D represents Phi, eta-balanced at w} H_Sigma(D),
```

with equality-type upper bounds (shrinking an oversized block `n_T` by `f` costs
`binom(|T|,2) log f` in the `t_ij`) when the trees are admissible. In particular **minimum
`H_Sigma` bounded along a sequence `<=>` `eps(Phi) -> 0` along it** (for admissible points).

## 2. Reduction of the census to a finite search; emptiness below an explicit scale

**Theorem E (emptiness).** Let `k >= 5`, `V_k = 2^(k-1) - k - 1` and `r_k = ceil(V_k/(k-3))`. Let
`D` be a datum (§1.1, any `k`) whose visible blocks satisfy `|log n_T - w| <= eta w` with
`0 <= eta < 1`, `w > 0`, and the two-chain condition. Then, with `p_j` the `j`-th prime:

1. at least `r_k` distinct primes divide `prod_T n_T`, and the `j`-th smallest visible block is
   `>= p_{ceil(j/(k-3))}`; in particular `max_T n_T >= p_{r_k}`;
2. hence `(1+eta) w >= log p_{r_k}` and `119(1+eta)w >= sum_j log p_{ceil(j/5)}` (for `k = 8`);
3. every minor satisfies `log|det(P_i,P_j)| >= (2^(k-3) - 1)(1-eta) w`.

The same holds for the contact integers `d_T(Phi)` of any `Phi` with `eps(Phi) <= eta`
(at a prime the positive edges form a tree on `k` leaves: at most `k-3` internal edges).

*Proof.* `n_T >= e^{(1-eta)w} > 1`, so each visible block has a prime factor; a prime divides at
most `k-3` visible blocks (§1.1). If the `j` smallest blocks all had all prime factors `< p_r`,
they would be covered by `r-1` primes, so `j <= (k-3)(r-1)`; this gives (1), and (2) follows from
`n_T <= e^{(1+eta)w}`. (3): `|t_ij| >= 1` and `2^(k-3) - 1` visible blocks contain `{i,j}`. `[]`

*Data* (`emptiness.py`, exact apart from printing logs):

| `k` | `V_k` | blocks/prime | min #primes `r_k` | `p_{r_k}` | class empty for `e^w <` (`eta = 0.1`; `0.25`; `0.5`) | then `|det| >=` for every pair (`eta = 0.1`) |
|---|---|---|---|---|---|---|
| 5 | 10 | 2 | 5 | 11 | 8.8; 6.8; 4.9 | `10^2.6` |
| 6 | 25 | 3 | 9 | 23 | 17.3; 12.3; 8.1 | `10^7.8` |
| 7 | 56 | 4 | 14 | 43 | 30.5; 20.3; 12.3 | `10^20.0` |
| 8 | 119 | 5 | 24 | 89 | 59.2; 36.3; 19.9 | `10^49.4` |

So **the `k = 8` census is provably empty below `e^w ≈ 59` (`eta = 0.1`), and every balanced
8-configuration has all 21 minors above `10^49`** (above `10^20` even at `eta = 1/2`). An exhaustive
search over integer configurations of bounded height can therefore never reach a balanced
8-point configuration: the smallest candidates have `|X_i| > 10^49`, in a 6-dimensional space.
Any census must be structured (lattices, curves, families). §§3-6 do this.

## 3. Algorithms and scripts

Three exact searches, all reporting the intrinsic statistics `eps(Phi)`, `w*(Phi)` of §1.2 (and
coverage `cov = #{T : d_T > 1}`), plus a heuristic count used only for interpretation.

1. **Exhaustive small-height census** (`census_small.py k B`). All `k`-point configurations
   `0 = x_0 < x_1 < ... < x_{k-1} <= B` (`gcd = 1`, infinity unmarked; every configuration with
   a representative of this shape is visited at least once). Each is anchored at `x_0`
   (`X_i = L/(x_i - x_0)`, divided by the gcd of differences) and its `V_k` contact integers are
   computed exactly by the residue-class rule of Lemma 1.1, with the two-chain test at every prime.
2. **Families** (`k6_family_census.py`). The rational points `Q_1 + N(Q_2 - Q_1)` of the `(L_6)`
   family `E_{t_0}` (`L8.md` Thm B, `t_0 = (2,3,-5,1)`), produced with the repository's
   `threshold/L8/k6fam/` scripts (copied to `census_checks/k6fam/`), then anchored and profiled
   with this census's own code. Contact integers of thousands of digits are obtained without
   factoring: in the integral chart `gcd_{i,j∈T}(X_i - X_j) = delta prod_{T' ⊇ T} d_{T'}`, inverted
   from the largest `T` down (`census_lib.contact_integers_gcd`).
3. **Lattice (CRT) census** (`lattice_census.py k A s seed`). For prescribed pairwise coprime
   prime blocks `n_T ∈ [A, 2A]` (so `eta <= log 2/log A`, two chains automatic), the configurations
   realising them are the lattice `Lambda(n) = {X ∈ Z^(m-1) : X_i ≡ X_j mod n_T (i, j ∈ T)}` of index
   `prod n_T^{|T|-1}`; for `X ∈ Lambda(n)`, `H_Sigma(X) = sum_{i<j} log(|X_i - X_j|/G_ij)` exactly.
   The script LLL-reduces (`python-flint`) and evaluates `H_Sigma` on all `5^(m-1) - 1`
   combinations with coefficients in `[-2, 2]`. The output is an **upper bound** for
   `min_{X ∈ Lambda(n)} H_Sigma`, hence for the census minimum at scale `w = log A`.

**Heuristic count (not used in proofs).** Per prescribed block vector, the lattice has
volume-to-target ratio `exp(((m-1) b_k - I_k) w)`, `b_k = 2^(k-3) - 1` blocks per pair,
`I_k = sum_T (|T| - 1)`; there are `≈ exp(V_k w)` block vectors. Expected number of bounded-`H`
balanced configurations at scale `w`, in units of `w` in the exponent:

| `k` | `V_k` | `(m-1) b_k` | `I_k` | per block vector `(m-1)b_k - I_k` | total `V_k + (m-1)b_k - I_k` | `ptolemy.md` `1 - K.beta` |
|---|---|---|---|---|---|---|
| 5 | 10 | 9 | 14 | -5 | +5 | 6 |
| 6 | 25 | 28 | 45 | -17 | +8 | 9 |
| 7 | 56 | 75 | 124 | -49 | +7 | 8 |
| 8 | 119 | 186 | 315 | -129 | **-10** | **-9** |

(The two last columns differ by the one invisible top block.) Two consequences matter for reading
the data:

* A census that samples block vectors (CRT, method 3) finds a bounded-`H` configuration only with
  probability `exp(-5w), exp(-17w), exp(-49w)` per sample at `k = 5, 6, 7`, although such
  configurations are abundant (`k = 5, 6`) or expected (`k = 7`). **Method 3 cannot see the
  difference between `k <= 7` and `k = 8`**; this is confirmed at `k = 5` in §4.1.
* Adding one point to a balanced `k`-configuration at scale `2w` costs `exp((b_{k+1} + 1 - V_k) w)`
  (`b_{k+1}` sizes of the new point's minors, one scale choice, `V_k` contacts that must split
  in two): `-2, -9, -24` for `5 -> 6, 6 -> 7, 7 -> 8`. These reproduce the totals
  (`5·2 - 2 = 8`, `8·2 - 9 = 7`, `7·2 - 24 = -10`), a consistency check of the count.

## 4. Calibration: k = 5 and k = 6

### 4.1 k = 5 (where `(L_5)` is false by `L8.md` Prop. C): the exhaustive census works, CRT does not

*Data* (`census_small.py 5 100`, 151 s): all 3,646,069 five-point configurations
`0 < x_1 < ... < x_4 <= 100` with `gcd = 1`. Coverage histogram (`#{T : d_T > 1}` of 10):

| cov | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| configurations | 23 | 1,614 | 23,999 | 164,555 | 537,322 | 919,346 | 1,223,776 | **775,434** |

21% of all configurations touch every boundary divisor, all admissible ones have `H_Sigma = 0`
(Lemma 1.1), and the best are already well balanced: `eps = 0.133` at `w* = 2.77`
(`(0, 22, 35, 46, 85)` and its reflection), `eps = 0.135` at `w* = 3.27` (`(0, 17, 29, 48, 99)`).
Frontier: `min eps = 0.133` over full-coverage points with `w* >= 2`, `0.135` with `w* >= 3`,
`0.32` with `w* >= 4` (the box `B = 100` then runs out). This reproduces Prop. C: at `k = 5`
balanced bounded-`H` points have positive density among small integer configurations, and the
census minimum of `H_Sigma` at these scales is `0`.

*Data* (`lattice_census.py`, 40 random prime-block vectors per row): the CRT method at the same
scales finds only `H_Sigma ≈ 8 w`:

| `k` | `A = e^w` | `w` | min `H_Sigma` found | `min H/w` | mean `H/w` | Gaussian-heuristic `H_typ/w` |
|---|---|---|---|---|---|---|
| 5 | 50 | 3.91 | 18.8 | 4.80 | 8.59 | 12.6 |
| 5 | 200 | 5.30 | 40.0 | 7.56 | 9.14 | 12.1 |
| 5 | 1000 | 6.91 | 59.5 | 8.61 | 9.37 | 11.5 |
| 5 | 5000 | 8.52 | 66.7 | 7.83 | 9.30 | 11.3 |

So **a census that samples block vectors returns `H_Sigma ≍ w` even where `H_Sigma = 0` points are
abundant** (the probability per sample is `≈ e^{-5w}`, §3). Its `k = 8` output (§6) is therefore
not evidence for `(L_8)`; it is reported only as an upper bound for the census minimum.

### 4.2 k = 6 (where `(L_6)` is false by `L8.md` Thm B): the known family is reproduced

*Data* (`k6_family_census.py`): points `Q_1 + N(Q_2 - Q_1)` on `E_{t_0}`, anchored at `x_1`, exact
contact integers by the gcd method (no factoring), `H_Sigma = 0` representation (Lemma 1.1).
`eps_Sigma` is the imbalance after the structural bad primes `Sigma` of `L8.md` §3.3 are moved into
the residues (bounded cost `<= 2 · 79`).

| `N` | `w*` | `eps` (intrinsic) | `eps_Sigma` | log of `Sigma`-part | `log max|X_i|` |
|---|---|---|---|---|---|
| 4 | 32.4 | 0.388 | 0.460 | 78.3 | 229 |
| 6 | 82.2 | 0.211 | 0.224 | 46.3 | 590 |
| 8 | 152.4 | 0.149 | 0.157 | 62.3 | 1092 |
| 10 | 243.7 | 0.119 | 0.132 | 73.3 | 1743 |
| 12 | 358.1 | 0.095 | 0.097 | 61.7 | 2485 |
| 14 | 493.9 | 0.078 | 0.080 | 52.6 | 3431 |

(all `N = 2..14` in `k6_family_census.out`; `N = 1` is degenerate.) The `Sigma`-free scale agrees
with `L8.md` §3.3 (e.g. `30.06` vs `30.2` at `N = 4`; `L8` uses the mean, here the midpoint).
`eps -> 0` like `O(1/N)` while `H_Sigma` stays bounded: **the `k = 6` census minimum of `H_Sigma`
is `O(1)` along `w ≈ 0.25 N^2 -> infinity`, as it must be.** Note the heights: `log max|X_i| ≈ 7 w*`,
i.e. the first balanced point (`eps < 0.4`) has `|X_i| ≈ 10^99`. No box search reaches it; the
family is found only because it is a curve with positive-rank arithmetic.

## 5. k = 7

* **Emptiness (Theorem E).** Balanced 7-configurations need at least 14 primes, a block `>= 43`,
  `e^w >= 30.5` at `eta = 0.1` (`>= 12.3` at `eta = 1/2`), and every minor `>= 10^20` at
  `eta = 0.1`.
* **Exhaustive census** (`census_small.py 7 40`, §6.1 table): no configuration of height `<= 40`
  touches more than 23 of the 56 boundary divisors.
* **Lattice census** (Data, 40 samples per row): `H_Sigma/w ≈ 150` at `A = 500, 2000, 10^4`
  (min 151.1, 150.0, 148.6; Gaussian heuristic 170, 165, 162). By §3 and §4.1 this says nothing
  about `(L_7)`.
* **What would test the heuristic "`(L_7)` is false".** The count predicts `≈ e^{7w}` balanced
  bounded-`H` configurations per shift class at scale `w`, but they sit at heights `>= 10^20` in a
  5-dimensional space, a fraction `≈ e^{-49w}` of the CRT block vectors, and no algebraic family
  carrying them is known (the dense good curves `V_A` of `round7/L7.md` have class
  pair 6 / triple 12, so most pair divisors are missed: `eps = 1`). Two-step lifting from the
  abundant `k = 5` supply costs `e^{-4w}·e^{-9w}` per 5-point (§3), i.e. `≈ 10^{-17}` at the smallest
  admissible scale. **No feasible census tests the `k = 7` heuristic**; it remains a heuristic.
  The structural facts are unchanged: `6K = -2 sum_{|S|=2} D_S`, so Prop. 7.5 is empty at `k = 7`
  and covering families of balanced class are not excluded (as `round4/L7.md` Thm 2 requires).

## 6. k = 8

### 6.1 Exhaustive small-height census (`k = 6, 7, 8`)

*Data* (`census_small.py`; every configuration `0 = x_0 < ... < x_{k-1} <= B`, `gcd = 1`, anchored at
`x_0`; contact integers exact):

| `k` | `B` | configurations | `V_k` | max coverage `#{T: d_T > 1}` | # attaining it | full coverage | seconds |
|---|---|---|---|---|---|---|---|
| 5 | 100 | 3,646,069 | 10 | 10 | 775,434 | 775,434 (21%) | 151 |
| 6 | 60 | 5,302,911 | 25 | 18 | 198 | 0 | 292 |
| 7 | 40 | 3,797,877 | 56 | 23 | 42 | 0 | 281 |
| 8 | 34 | 5,359,838 | 119 | 29 | 1 | 0 | 554 |

Full coverage histograms are in `census_small.out`; at `k = 8` the median coverage is 20 of 119.
At `k = 6, 7, 8` no configuration in these boxes touches every boundary divisor (so `eps = 1`
throughout), in line with Theorem E: a balanced configuration needs at least 9, 14, 24 distinct
primes in its contacts. The census shows how far small points are from the balanced regime. It is
data about small heights only and says nothing about `(L_k)`.

### 6.2 Lattice upper bounds for the census minimum

*Data* (`lattice_census.py`, 40 prime-block vectors per row, `n_T ∈ [A, 2A]` pairwise coprime,
`eta <= log 2/w`):

| `k` | `A` | `w` | min `H_Sigma` found | `min H/w` | mean `H/w` | heuristic `H_typ/w` |
|---|---|---|---|---|---|---|
| 6 | 200 | 5.30 | 220.4 | 41.6 | 42.8 | 51.0 |
| 6 | 1000 | 6.91 | 289.0 | 41.8 | 42.6 | 48.5 |
| 6 | 5000 | 8.52 | 348.1 | 40.9 | 42.6 | 47.5 |
| 7 | 500 | 6.22 | 939.0 | 151.1 | 152.2 | 170.3 |
| 7 | 2000 | 7.60 | 1139.8 | 150.0 | 151.1 | 165.5 |
| 7 | 10000 | 9.21 | 1369.1 | 148.6 | 150.5 | 162.3 |
| 8 | 1500 | 7.31 | 3426.1 | 468.5 | 470.6 | 509.8 |
| 8 | 5000 | 8.52 | 3966.7 | 465.7 | 467.8 | 501.3 |
| 8 | 20000 | 9.90 | 4591.2 | 463.6 | 465.2 | 493.7 |

These are exact **upper** bounds: there are balanced (pairwise coprime prime blocks, `eta < 0.1`)
8-configurations at `w = 7.3, 8.5, 9.9` with `H_Sigma <= 3426, 3967, 4591`; the census minimum is
at most `≈ 464 w` there, and the minimum over the sampled block vectors is concentrated (spread
`< 1%`), as the Gaussian heuristic predicts for random lattices of dimension 6. The same method
gives `≈ 42 w` at `k = 6` and `≈ 8 w` at `k = 5`, where the true minimum is `O(1)` (§4). The
slopes are the dimension count `binom(m,2)(I_k/(m-1) - b_k)` (`= 11.5, 47.5, 162, 493` for large
`A`, minus the gain from searching small combinations), not information about `(L_k)`.

### 6.3 What the k = 8 census establishes

* Nothing below `e^w ≈ 59` (`eta = 0.1`) exists at all, and nothing below `|det| ≈ 10^49`
  (Theorem E); a literal small-height census is empty by theorem.
* No symmetric configuration, no configuration with a fixed sub-configuration (forgetful fibre)
  and no boundary-stratum point is balanced (Cor. 7.2, 7.3).
* No family of curves of genus `<= 1` and bounded degree, and no product chart, produces balanced
  points densely (Prop. 7.5, Cor. 7.6); on each member of a covering family `(L_8)` holds
  (Prop. 7.5(4)).
* Every balanced 8-point lies over balanced 7- and 6-points (Lemma 7.7); over the `(L_6)` curves the
  heuristic expectation is finite (§7.3).
* Constructed balanced 8-configurations have `H_Sigma ≈ 464 w` (upper bound; §6.2).

No small `H_Sigma` was found at `k = 8`, but the searches that could have found one are
demonstrably blind at `k = 5, 6` (§4.1), so this absence is weak evidence (as the round-17 judge
anticipated).

## 7. Where small-height points could lie: candidate loci, tested by theorems

Theorem E makes a literal census of `k = 8` points of bounded height useless, so the question "do
the small points lie on explicit subvarieties?" is answered here locus by locus, rigorously. Notation:
`X = Mbar_{0,8}`, `Phi ∈ M_{0,8}(Q)`, `ell_T`, `eps`, `w*` as in §1.2.

### 7.1 Loci missing a boundary divisor carry no balanced points

**Proposition 7.1.** Let `Y ⊂ X` be an irreducible closed subvariety over `Q` meeting `M_{0,8}`,
and let `T` be visible with `Y ∩ D_T = ∅`. Then there is an integer `N_Y >= 1` with
`d_T(Phi) | N_Y` for every `Phi ∈ Y(Q) ∩ M_{0,8}`. Consequently every datum `D` representing a
point of `Y` with `eta`-balanced blocks at scale `w` has `H_Sigma(D) >= ((1-eta) w - log N_Y)/2`;
so `(L_8)` holds on `Y` (with `c = (1-eta_0)/2`).

*Proof.* Let `𝒴` be the closure of `Y` in Knudsen's smooth projective model `Mbar_{0,8}/Z` and
`𝒟_T` the boundary divisor over `Z`. `𝒴 ∩ 𝒟_T` has empty generic fibre, so it lies over finitely
many primes, and its structure sheaf is killed by some integer `N_Y` (a coherent sheaf on a
noetherian scheme is killed by a power of the ideal of its support). The section
`s_Phi: Spec Z -> 𝒴` extending `Phi` (properness) pulls the ideal of `𝒴 ∩ 𝒟_T` back to the ideal
generated by `s_Phi^* f_T` (`f_T` a local equation of `𝒟_T`), which therefore contains `N_Y`. Hence
the local intersection number `(Phi.D_T)_p = v_p(s_Phi^* f_T) <= v_p(N_Y)`, and
`(Phi.D_T)_p = l_T(p)` (`ptolemy.md` Thm B(iv)(1)). The last claim is Lemma 1.2:
`log n_T <= ell_T + 2H_Sigma <= log N_Y + 2 H_Sigma` and `log n_T >= (1-eta) w`. `[]`

**Corollary 7.2 (symmetric loci).** For a label permutation `sigma ≠ id` of `{0, ..., 7}` let
`Y_sigma` be the closed set of configurations fixed by `sigma` (some projectivity induces `sigma`).
Then `d_T(Phi) = 1` for **at least 48** of the 119 visible `T`, for every `Phi ∈ Y_sigma(Q) ∩ M_{0,8}`.
So `eps(Phi) = 1` on every symmetric locus.

*Proof.* If `T` and `sigma T` are incompatible splits, `D_T ∩ D_{sigma T} = ∅` already over `Z`
(boundary divisors of incompatible splits are disjoint in Knudsen's model). The `Z_p`-point
`s_Phi` is `sigma`-fixed (uniqueness of the extension), so its reduction lies on `D_T` iff on
`D_{sigma T}`; hence `l_T(p) = 0` for all `p`. `check_symmetric_loci.py` (Data, all of `S_8`)
finds, for every `sigma ≠ id`, at least 48 visible `T` with `T`, `sigma T` incompatible (minimum
attained by the fixed-point-free involutions, cycle type `2^4`; e.g. for `(12)(34)(56)(78)` the
24 splits of type `{a, sigma a, b}` and the 24 of type `{a, sigma a, b, c}`). At `k = 5, 6, 7` the
minima are 3, 6, 25. `[]`

**Corollary 7.3 (forgetful fibres).** Let `Q ⊆ {0..7}`, `4 <= |Q| <= 7`, `pi_Q` the forgetful map and
`c ∈ M_{0,Q}(Q)`. For `Phi ∈ pi_Q^{-1}(c)` and each boundary split `S` of `Q`,
`sum_{T restricts to S} ell_T(Phi) = ell_S(c)`, a constant. Exactly 48, 80, 100, 112 visible `T`
restrict to a boundary split of `Q` for `|Q| = 4, 5, 6, 7` (`check_forgetful_fibres.py`; `3·16 = 48`
by hand for `|Q| = 4`), and all of them have `ell_T(Phi) <= max_S ell_S(c)`. So balanced points on a
fibre have `(1-eta) w <= max_S ell_S(c)`.

*Proof.* `pi_Q^* D_S = sum_{T ↦ S} D_T` (reduced, over `Z`), and local intersection numbers are
compatible with pull-back. `[]`

These cover the natural "special" loci one would first suspect: configurations with extra
projective symmetry, configurations containing a fixed sub-configuration (harmonic quadruples,
equianharmonic ones, any fixed cross-ratio), and boundary strata (which contain no points of
`M_{0,8}`).

### 7.2 Curves: balanced points force the balanced class, and no covering family has it

**Proposition 7.4 (curves).** Let `C ⊂ X` be an irreducible curve meeting `M_{0,8}`, of geometric
genus `<= 1`, with class `gamma`. Suppose `Phi_n ∈ C(Q)` are distinct, represented by data `D_n`
with balance `eta_n -> 0`, `w_n -> infinity` and `H_Sigma(D_n) = o(w_n)`. Then
`gamma.D_T = gamma.D_{T'}` for all visible `T, T'`, i.e. `gamma ∈ Q_{>0} beta`; then
`K.gamma = 10 gamma.D_T > 0`.

*Proof.* `ptolemy.md` Lemma 6.2 (refereed) gives, for every visible `T`,
`(1-eta)w - 2H - c_0 <= h_{D_T}(Phi) <= (1+33 eta) w + 3H + c_0`. On the normalisation `C~`
(genus `<= 1`), `h_{D_T} = (gamma.D_T/gamma.A) h_A + O_C(sqrt(h_A) + 1)` (heights of a degree-0
class on a curve of genus `<= 1`), and `h_A ≍ w` by the same lemma. Divide by `w_n` and let
`n -> infinity`: `(gamma.D_T/gamma.A) lim h_A/w_n = 1` for every `T`. `[]`

**Proposition 7.5 (families threshold at `k = 8`).**

1. *(Class inequality.)* `7K = -2 sum_{|S|=2} D_S + sum_{|S|=3} D_S + 2 sum_{|S|=4} D_S`
   (28, 56, 35 divisors; `ptolemy.md` §6.2 Step 3). If `gamma.D_S ∈ [(1-eta)a, (1+eta)a]` for all
   `S` then `7 K.gamma >= a(70 - 182 eta)`. So `K.gamma <= 0` forces
   `eta(gamma) := (max_S gamma.D_S - min_S gamma.D_S)/(max + min) >= 5/13`.
2. *(Covering families.)* If irreducible curves of geometric genus `<= 1` and class `gamma` cover
   a dense subset of `X`, then `K.gamma <= 0` (and `K.gamma <= -2` in genus 0). Hence
   `eta(gamma) >= 5/13`, and `beta` is not in the closed cone spanned by such classes.
3. *(Bounded degree.)* For every `d`, the union of all curves of geometric genus `<= 1`, degree
   `<= d` and `eta(gamma) < 5/13` lies in a proper closed subset `Z_d ⊂ X`.
4. *(Points.)* If `C` is a fixed curve of genus `<= 1` with `K.C <= 0` (e.g. any member of a covering
   family), every datum representing a point of `C` satisfies
   `60 H_Sigma >= (10 - 282 eta) w - O_C(sqrt(w + H_Sigma) + 1)`. So `(L_8)` holds on `C` for every
   `eta_0 < 5/141` and `c < (10 - 282 eta_0)/60`.

*Proof.* (1) `7K.gamma >= -2·28(1+eta)a + (56 + 70)(1-eta)a`. (2) After shrinking and base change, the family is
`F: 𝒞 -> X` with `𝒞 -> B` smooth projective of relative genus `g <= 1` over an irreducible `B`,
members `f_b: C_b -> X` birational onto their images, and `F` dominant. In characteristic 0, `dF`
is surjective at a general point `(b, c)` (generic smoothness). Its image is
`df_b(T_{C_b,c})` plus the values at `c` of the first-order deformations along `B`, which are
global sections of the normal sheaf `N_b = coker(T_{C_b} -> f_b^* T_X)`. So `N_b` modulo torsion is
generically generated by global sections, and its determinant has a nonzero section. Hence
`deg N_b >= 0`, i.e. `-K.gamma - (2 - 2g) >= 0`, that is `-K.gamma >= 2 - 2g >= 0`. For `g = 0` this
is the bound `-K.gamma >= 2` for free curves (Kollár, *Rational curves on algebraic varieties*,
II.3.11; `round4/L7.md` Prop. 3(i) uses it at `k = 7`). Then apply (1). If `beta = sum alpha_k gamma_k`,
`alpha_k >= 0`, then `10 = K.beta = sum alpha_k K.gamma_k <= 0`, impossible. (3) Curves of bounded
degree lie in finitely many components of the Chow variety; the sublocus of genus `<= 1` is
constructible, with finitely many components; each covering component has `eta >= 5/13` by (2)
(the class is constant on a component); the others sweep out proper closed subsets. (4) `ptolemy.md`
Step 3 gives `h_K(Phi) >= 10w - 282 eta w - 60 H_Sigma - c`; on `C`,
`h_K = (K.gamma/gamma.A) h_A + O_C(sqrt(h_A) + 1) <= O_C(sqrt(h_A) + 1)` and `h_A = O(w + H_Sigma)`. `[]`

**Corollary 7.6.** At `k = 8` none of the three refutation criteria of `round4/L7.md` Thm 2 can
apply with families of bounded degree: criterion (b) (a dominant product chart with `beta` in the
cone of its ruling classes) never applies; criteria (a) and (c) (rational, resp. positive-rank
genus-one curves with classes near the `beta`-ray and dense union) need curves of unbounded degree
that lie in no covering family, i.e. infinitely many rigid curves. Every mechanism that refuted
`(L_5)` and `(L_6)` (Prop. C and the linear family C1: covering rational families of class
`beta_5`, `K.beta_5 = -5`; `L8.md` Thm B: a covering genus-one family of class `beta_6`,
`K.beta_6 = -8`) is of the excluded kind. For comparison, at `k = 7` the same inequality is empty:
`6K = -2 sum_{|S|=2} D_S`, so every class has `K.gamma <= 0`.

So a counterexample to `(L_8)`, if one exists, lives on no special locus of §7.1 and on no
bounded-degree family of curves of genus `<= 1`; by Prop. 7.4 the only curves of genus `<= 1` that
can carry one have class proportional to `beta` with `K.C = 10 d > 0`. **The set `Z` that any proof
must discard therefore contains every curve of geometric genus `<= 1` and class `d·beta` with
infinitely many rational points (if such curves exist); no such curve is known, and none was found
(§6).** That is the only concrete shape of `Z` the analysis produces.

### 7.3 Forgetful monotonicity and the pull-backs of the `(L_6)` families

**Lemma 7.7.** Let `pi_j: Mbar_{0,8} -> Mbar_{0,7}` forget a label `j ≠ 0`. For every visible `S` of
the 7-point configuration, `ell_S(pi_j Phi) = ell_S(Phi) + ell_{S ∪ {j}}(Phi)`, both terms visible.
Hence if all `ell_T(Phi) ∈ [(1-eta)w, (1+eta)w]` then all `ell_S(pi_j Phi) ∈ [(1-eta)2w, (1+eta)2w]`.
At the level of data: deleting row `j` from an 8-datum `D` and merging `n'_S = n_S n_{S ∪ j}`
(top: `n'_{[7]\j} = n_{[7]\j} n_{[7]}`) gives a 7-datum `D'` with the **same** `G_ab`, `t_ab`, `Y_b`
(so `H_Sigma(D') <= H_Sigma(D)`), blocks balanced at scale `2w` with the same `eta`, and still two
nested chains at every prime. Consequently `(L_7)` with `(c, eta_0, Z_7)` implies `(L_8)` with
`(2c, eta_0, pi_j^{-1} Z_7)`; the same holds for every `k -> k+1`.

*Proof.* `pi_j^* D_S = D_S + D_{S ∪ j}` over `Z`, `|S ∪ j| <= 6`. For data: `T ↦ T \ {j}` is a
2-to-1 map from the blocks containing `a, b` (`a, b ≠ j`) onto the 7-point blocks containing
`a, b`, so `G'_ab = G_ab`; deleting `j` keeps chains nested and supports disjoint. `[]`

So a balanced 8-point configuration has 8 balanced 7-point and 28 balanced 6-point
sub-configurations (scales `2w`, `4w`). Every `k = 8` census point lies over the `k = 6` census,
which is Zariski dense (`L8.md` Thm B); the question is only whether `k = 8` points concentrate over
**special** 6-points, e.g. over the curves `E_t`. The 3-folds `pi^{-1}(E_t)` meet every `D_T`, so
Prop. 7.1 does not apply. *Heuristic* (§3): over a balanced 6-point at scale `4w`, adding two
points costs `exp(-9·2w)·exp(-24w) = exp(-42w)`, while `E_t` has only `O(sqrt W)` points of scale
`<= W` (they have `w_6 = alpha N^2`). The expected number of balanced bounded-`H` 8-points on
`pi^{-1}(E_t)` is therefore finite (`sum_N exp(-42 alpha N^2/4) < infinity`). The pull-backs of the
`(L_6)` families are not, heuristically, where small 8-points concentrate.

## 8. Verdict among (a), (b), (c), and limits of the evidence

**Formulation re-checked first, as instructed.** `(L_8)` was used verbatim (§1.1). All
normalisations were checked against `ptolemy.md` §1, §7 and Steps 0-5: the anchor, visible versus
invisible blocks, the two-chain condition with at most 5 visible blocks per prime, `w`, `eta` on
the visible blocks, and `H_Sigma` with 7 + 21 residues. Prop. 1.3 shows that this object is the
contact profile of `Phi`, and Lemma 7.7 shows the statements are monotone in `k`. One typo was found:
`ptolemy.md` §1 gives `2^m - m - 1` for the number of visible blocks, where `2^m - m - 2` is
correct. It is harmless. Nothing in the census suggests that `(L_8)` is mis-stated.

* **(a) Bounded-`H` points Zariski dense: not supported.** No bounded-`H` balanced 8-point was found.
  Every mechanism known to produce dense bounded-`H` points at `k <= 6` is provably unavailable at
  `k = 8` with bounded degree: covering curve families of genus `<= 1` and product charts
  (Prop. 7.5, Cor. 7.6). A dense supply would need infinitely many rigid curves of class `∝ beta`,
  or sporadic points on no curve of genus `<= 1`. Either would contradict Vojta (`D = 0`) on
  `Mbar_{0,8}`.
* **(b) Concentration on explicit subvarieties: no positive candidate. The natural candidates are
  excluded by theorem, not just untested.**
  * Symmetric loci: `eps = 1` identically, since at least 48 contacts are trivial (Cor. 7.2).
  * Forgetful fibres and fixed sub-configurations (Cor. 7.3).
  * Boundary and Kapranov strata (no interior points). Points near a stratum over-contact only
    the stratum's divisors, so they are unbalanced unless all 119 are hit equally.
  * The pull-back 3-folds `pi^{-1}(E_t)`: these are not excluded by theorem, but heuristically
    they carry finitely many balanced points (§7.3).

  The only shape of `Z` the analysis produces is **`Z ⊇` the curves of geometric genus `<= 1`
  and class `d·beta`, `d >= 1`, with infinitely many rational points (Prop. 7.4)**. In Kapranov's
  model (`Mbar_{0,8}` = iterated blow-up of `P^5` in 7 general points `p_i` and the spans of at
  most 4 of them; `D_{{8} ∪ A} = E_A` for `|A| <= 4`, the strict transform of the hyperplane `<p_A>`
  for `|A| = 5`), class `beta` means: a curve of degree 31 through `p_1, ..., p_7`, meeting each
  of the 21 lines, 35 planes and 35 solids `<p_A>` once (strict transforms), and meeting each of the
  21 hyperplanes `<p_A>` (`|A| = 5`) in exactly one point outside those flats (`31 - 30 = 1`). In
  the integral chart of §1.2, with the anchor fixed at infinity, a genus-0 member is
  `X_i(s) - X_j(s) = c_ij prod_{T ∋ i,j} (s - r_T)` with 119 distinct `r_T` and degree
  `31 = 2^5 - 1`. Their expected dimension is `-10d + 2 < 0`, and none is known. This is
  the first concrete candidate in equations, but it may well be empty.
* **(c) Nothing small exists: this is what the census shows.** The evidence for `(L_8)` is weak, for
  two reasons.
  1. The only searches that can reach balanced 8-configurations (Theorem E puts them above
     `10^49`) are structured. The structured search that samples block vectors returns
     `H_Sigma ≍ w` also at `k = 5, 6`, where `H_Sigma = 0` is abundant (§4.1, §6.2).
  2. The known bounded-`H` supplies at `k = 5, 6` were found by algebraic geometry (integral
     charts, positive-rank genus-one curves), not by search. At `k = 8` the corresponding
     algebraic geometry is exactly what Prop. 7.5 excludes.

  So the census agrees with Vojta's prediction and with the heuristic count (`-10` per shift class
  at `k = 8`), but it adds no quantitative support beyond them.

**What the census does add.**
1. The reduction of `(L_8)` to finite boundary heights (Prop. 1.3), with explicit constants.
2. The explicit scale below which the problem is vacuous (Theorem E).
3. A rigorous list of loci that any refutation must avoid (§7).
4. The families threshold (Prop. 7.5). It is the geometric counterpart of Theorem D's `K.beta > 0`:
   at `k = 8` bounded-degree families cannot refute `(L_8)`, just as at `k = 8` (and not before)
   Vojta's inequality says something. This bears on the judge's estimate in round 17 (§7.2): a
   refutation of `(L_8)` cannot be found the way `(L_5)` and `(L_6)` were refuted. A proof of
   `(L_8)` still needs a Vojta-strength input, and none is available (round 17 §7.1, items 3 and 5).

**Not claimed.**
* No growth improvement.
* No uniform bound.
* No unconditional case of `(L_8)` beyond the loci of §7.1 and the individual curves of
  Prop. 7.5(4). On those loci `(L_8)` holds trivially or by the curve's own height theory, and
  together they are far from Zariski dense.
* No claim about `(L_7)`. Its heuristic falsity is untestable by any feasible census (§5).

## 9. Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round18/census_checks/`; run with `python3` from that
directory (memory-limited with `prlimit --as=3500000000`).

| file | content | type |
|---|---|---|
| `census_lib.py` | contact integers by residue classes (Lemma 1.1); factoring-free gcd version; admissibility (two chains); Boolean identity check; integral chart | exact |
| `test_lib.py` | 284,752 comparisons with the Thm B(iv) min-formula, 0 mismatches | exact |
| `emptiness.py`, `emptiness.out` | Theorem E numbers | exact (logs printed) |
| `check_forgetful_fibres.py` | Cor. 7.3 counts 48, 80, 100, 112 | exact |
| `verify_lattice_example.py` | independent exact divisibility check of one `k = 8` lattice output (`H = 3442.5`, `w = 7.31`, `eta = 0.095`, `|X_i| <= 10^175`) | exact |
| `check_symmetric_loci.py` | Cor. 7.2: for all `sigma ≠ id` in `S_k`, `k = 5..8`, number of `T` with `T`, `sigma T` incompatible (min 3, 6, 25, 48) | exact |
| `census_small.py`, `census_small.out` | exhaustive small-height census, `k = 5..8` | exact |
| `k6_family_census.py`, `k6_family_census.out`; `k6fam/` (copied from `threshold/L8/k6fam/`) | `k = 6` calibration on `E_{t_0}` | exact |
| `lattice_census.py`, `lattice_census.out` | CRT/lattice census, upper bounds for min `H_Sigma`, `k = 5..8` | exact values; logs printed |
| `../census_notes.md` | working notes | - |
