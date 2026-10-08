# Constructions and data: how large is M(C)?

Route: constructions and data. Could the conjecture fail, and what is the true `M(C)`?
All files referenced below are under `scratchpad/growth/families/`; the independent exact
checker is `families/check_families.py` (output in `families/check_families_output.txt`).

Normalisation throughout: a cluster of lattice points of one modulus `R` is first divided by its
Gaussian gcd (this keeps the points distinct and can only decrease the constant). Its constant is
`C = (arc length of the shortest arc containing it)/sqrt(R)`. `M(C)` is the supremum, over **all**
radii and arc positions, of the number of points in such an arc. A single cluster therefore gives
a lower bound for `M(C)`. An infinite family is needed only for statements of the form "for every
`C > C_0`" or "as `R -> infinity`".

## 0. Summary

**Rigorous new results.**

1. *Golden eight-point family* (Section 2). With `J_r = f_{r+1} + i f_r` (Fibonacci) and `F = 1+2i`,
   the eight odd-parity words in the four blocks `J_n, ..., J_{n+3}`, each with prefactor `F` or `conj F`,
   give 8 distinct points with
   `C -> C_8 = sqrt2 (4+sqrt5)/5^(1/4) = 5.897708962894...`
   (along `n = 1 mod 3`). Hence **`M(C) >= 8` for every `C > 5.8978`**, and `M(C) >= 7` for every
   `C > sqrt2 phi^2 5^(1/4) = 5.5364678...`. The eight-point family recorded in the notes (translated
   Pell, unit `9+4 sqrt5`) has `C = 44,652,598.0016` (computed exactly here). The new family improves
   that constant by a factor of about `7.6*10^6`.
2. *Axis six-point Pell family* (Section 3.2, elementary). For every solution of `y_2^2 - 2y_1^2 = -146`
   with `y_1` odd, put `a = (y_1^2-143)/2`. Then `(a,±12)`, `(a-1,±y_1)`, `(a-2,±y_2)` lie on
   `x^2+y^2 = a^2+144`, with `C = 2 sqrt(4+140/a)(1+o(1)) -> 4` from above. Hence **`M(C) >= 6` for every
   `C > 4`**. This beats the golden six-point template (`sqrt2 phi^3/5^(1/4) = 4.0062`) and the
   Cilleruelo–Granville six points. The latter have `C -> 4 sqrt2` once the common factor `(1+i)^2` is
   removed (the notes record 8 without the gcd).
3. *Single record clusters* (exact, finite verification). There are 8 points with `C = 5.673925`
   (`R^2 = 32988780417125 = 5^3 13^2 17^2 89 109 557`) and 7 points with `C = 5.269927`
   (`R^2 = 235370525`). So `M(5.6740) >= 8` and `M(5.2700) >= 7`.
4. *Theorem A* (Section 5). In every fixed-shift quadratic-unit template (the class of
   `quadratic_unit_template_bound.md`), a sub-cluster whose normalised arc constant tends to 0 has
   **at most 4 points**. The bound is sharp: the notes' canonical four-point family attains it. Five
   points with `C -> 0` therefore need growing rates, as in the cyclotomic five-row family.
5. *Symmetric clusters* (Section 3, finite congruence computations). No circle has lattice points
   at four consecutive abscissae (mod 8). For a cluster closed under a lattice reflection, with
   `R -> infinity`:
   - 6 points need `C >= 4 - o(1)`;
   - 8 points need `C >= 4 sqrt2 - o(1)`, and the only extremal depth pattern is `(0,1,3,4)`;
   - 10 points, or 9 with a point on the axis, need `C >= 8 - o(1)`;
   - 12 points need `C >= 6 sqrt2 - o(1)`.

   The golden family is such a symmetric cluster, with depth pattern `(0,1,3,4)`.

**Data (exhaustive, all circles `x^2+y^2 = n`, `n = R^2 <= 10^12`, about `2*10^10` primitive circles).**
- The largest cluster on an arc of length `(1/2)sqrt R` has **3** points (first at `n = 276108509`).
- On an arc of length `sqrt R` it has **4** points (first at `n = 412333805`).
- No 5-point cluster has `C < 1.996`. No 6-point cluster has `C < 4.0001`; all 6-point records are
  axis clusters converging to 4. No cluster of 9 or more points has `C <= 6`.
- Beyond the scan, targeted searches of symmetric clusters up to `R^2 <= 4*10^20` find exactly six
  8-point clusters of depth type `(0,1,3,4)` (four golden, two with `C = 5.674` and `5.754`), and
  no 10-point cluster of the minimal depth types (`C` up to about 10).

**Answer to "what is the true M(C)?"** Rigorously, `M(C) >= 5` for every `C > 0` (the notes'
cyclotomic family, re-run here, at radii around `10^3428`); `M(C) >= 6` for `C > 4`; `M(C) >= 7` for
`C >= 5.27`; `M(C) >= 8` for `C >= 5.674`. No family or configuration with 9 or more points at
`C <= 6` is known; all known mechanisms stop at 8 (bounded `C`) and at 5 (`C -> 0`). The heuristics
(Section 7) make a uniform bound plausible and predict **`M(1/2) = M(1) = 5`**. The fifth point is
supplied only by the cyclotomic family at astronomically large radii; for all `R <= 10^6` the true
values are 3 and 4.

**Status.** These are constructions, a restricted theorem (Theorem A) and data. They give no new
general upper bound: the general growth rate `(1+eps) log R/loglog R` is unchanged. The missing lemma
is the one stated in the handoff (Section 8).

## 1. Normalisation and the three mechanisms that beat randomness

For `k` points of one circle in an arc of length `C sqrt R`, randomness would require `k-1`
independent phase coincidences of size `R^(-1/2)`. Three algebraic mechanisms avoid this.

- **Phase-locked blocks (quadratic-unit templates).** Products of blocks `H_r = a lambda^r + tau(a) lambda'^r`
  all lie within angle about `|H_r|^(-2)` of one ray. The `k-1` phase conditions then hold identically
  along a one-parameter family. The cost is rigidity: at most 8 rows (notes), and at most 4 rows with
  `C -> 0` (Theorem A).
- **Resonant cyclotomic layers** (notes): growing rates `G_{m d}` with shared cyclotomic factors.
  These give 5 rows with `C -> 0`, and at most 2251 in that model.
- **Reflection symmetry (new here).** A cluster closed under `z -> conj z` consists of pairs
  `(x, ±y)` on the vertical lattice lines near the tangent point `(sqrt n, 0)`. Each line costs one
  arithmetic condition but yields two points, so `2m` points cost only `m-1` conditions.

## 2. Golden-unit fixed-shift templates (the improved eight-point family)

**Construction.** Let `f_r` be the Fibonacci numbers, `J_r = f_{r+1} + i f_r` (so `|J_r|^2 = f_{2r+1}`),
`phi = (1+sqrt5)/2`, `q = -phi^(-2)`, `F = 1+2i`. For a word `w in {0,1}^4` with `p = |w|` odd, put

    Z_w(n) = F^[p=1] conj(F)^[p=3] * prod_{s=0}^{3} J_{n+s}^{w_s} conj(J_{n+s})^{1-w_s}.

**Proposition G.** For `n = 1 mod 3` the eight points `Z_w(n)` are distinct and have equal norm.
They are all divisible by `(1+i)^2`, and their normalised constant tends to
`C_8 = sqrt2(4+sqrt5)/5^(1/4) = 5.8977089628...`. Their best 7-point sub-cluster tends to
`C_7 = sqrt2 phi^2 5^(1/4) = 5.5364678128...`.

*Proof.* Binet's formula gives `J_r = a phi^r (1 - i q^r/phi)` with `a = (phi+i)/sqrt5`. Hence
`arg J_r = theta - delta_r`, where `theta = arg a`, `2 theta = arg F`, and
`delta_r = atan(q^r/phi) = q^r/phi + O(phi^(-6r))`. Note `tau(a)/a = -i/phi` is purely imaginary.

The leading phase of `Z_w` is `(2p-4) theta ± 2 theta = 0`: the prefactor `F` cancels `-2 theta` when
`p = 1`, and `conj F` cancels `+2 theta` when `p = 3`. So

    arg Z_w = -sum_s (2w_s - 1) delta_{n+s} = -(q^n/phi) sum_s (2w_s - 1) q^s + O(phi^(-6n)).

Over the odd words, `sum_s w_s q^s` has maximum `1 + q^2 + q^3` (at `w = 1011`) and minimum `q`
(at `w = 0100`). Their difference is `1 - q + q^2 + q^3 = 2 sqrt5 - 3`. The angular spread is therefore
`2|q|^n (2 sqrt5 - 3)/phi (1 + o(1))`.

Since `3 | 2n+1` and `3 | 2n+7`, the norms `f_{2n+1}` and `f_{2n+7}` are even with 2-adic valuation 1.
So `(1+i)` divides `J_n`, `J_{n+3}` and their conjugates exactly once, and `(1+i)^2` divides every row.
Dividing by it,

    R^2 = 5 prod_s f_{2n+2s+1} / 4 = phi^(8n+16)/20 (1 + o(1)).

Hence

    C = R * spread / sqrt R = 2 phi^3 (2 sqrt5 - 3)/20^(1/4) = sqrt2(4+sqrt5)/5^(1/4),

using `phi^3 (2 sqrt5 - 3) = 4 + sqrt5`. Any further common factor would only decrease `C`;
numerically the gcd norm is exactly 4.

Dropping an extreme word changes the spread to `1 - q = sqrt5/phi`, which gives `C_7`. Distinctness
holds because different words differ in the orientation of a block of norm greater than 5. ∎

**Verification.** `check_families.py` computes the actual clusters at `n = 100, 160, 220` and
matches all four closed forms to `10^(-59)`. The member `n = 4` has `R^2 = 537606725 = 5^2*17*61*89*233`;
it is exactly the exhaustive scan's 8-point record (`C = 5.897763`).

**Systematic search.**
- `golden_search_E4.py` covers every golden template of total width 4: four binary shifts up to 10,
  and widths `(2,1,1)`, `(2,2)`, `(3,1)`, `(4)`, in both parity classes, over a full period of `n`.
  Best constants per `k`:
  - `k = 5`: `sqrt2 5^(1/4) phi = 3.4217`;
  - `k = 6`: `sqrt2 phi^3/5^(1/4) = 4.0062` (shifts 0,1,2, widths 1,2,1);
  - `k = 7`: `5.5365`;
  - `k = 8`: `5.8977`.
- Silver unit (`P_{r+1} + i P_r`, Pell numbers), `silver.py`: `k = 5`: `2 + 2 sqrt2`;
  `k = 6`: `4 + sqrt2`; `k = 7`: `11.657`; `k = 8`: `11.899`.
- `units_search.py` covers 2-term recurrences `H_r = (A u_{r+1} + B u_r) + i(C u_{r+1} + D u_r)` for
  the units of `X^2 - tX ∓ 1` (`t <= 4`) and `|A|, ..., |D| <= 2`, with prefactor alignment searched
  automatically. Nothing beats the golden values.
- For a Lucas-type sequence `u_{r+1} + i u_r` the general formula is
  `C = 2|c|^(1/2) lambda^(S/2) spread / sqrt(t^2+4)`, divided by `sqrt|g|`. The golden unit
  minimises it.
- Widths `E > 4` would need contact order 3 for every pair (`eta` is imaginary, so even orders vanish).
  The searches `golden_nu3*.py` find at most 4 such rows in width at most 12.

By the template theorem of the notes, 8 is the maximum for fixed-shift quadratic templates.
**No template with 9 or more points was found**, consistent with that theorem.

## 3. Reflection-symmetric ("axis") clusters

**3.1 Reduction.** Let a cluster be closed under `z -> conj z` (after a unit rotation; the diagonal
reflections `z -> ±i conj z` are treated in 3.4). Its points are pairs `(x_j, ±y_j)` with distinct
`x_j`, near `(sqrt n, 0)`. Write `a = max x_j` and `d_j = a - x_j` (the *depths*, `D = {d_j}`). Then

    y_d^2 = y_0^2 + 2da - d^2   (d in D),   C = 2 y_{dmax}/n^(1/4) * (1 + o(1)),
    C^2 -> 8 dmax + 4 y_0^2/a.

Eliminating `a` leaves `m - 2` quadrics in the `m = |D|` unknowns `y_d`. The projective closure is a
complete intersection of `m-2` quadrics in `P^m`:

| `m` | surface | consequence |
|---|---|---|
| 3 | quadric | infinitely many points: Pell orbits |
| 4 | quartic del Pezzo; affine part = complement of an elliptic curve (log Calabi–Yau) | |
| 5 | K3 | |
| 6 | general type | finitely many under Bombieri–Lang, apart from special curves |

This is a circle analogue of Büchi's problem (second difference `-2`).

**3.2 Six points: the Pell family (elementary proof).** Take `D = (0,1,2)` and `y_0 = 12`. The
conditions reduce to `y_2^2 - 2y_1^2 = -146` with `a = (y_1^2 - 143)/2`, and directly

    (a-1)^2 + y_1^2 = a^2 + 144,   (a-2)^2 + y_2^2 = a^2 + 144.

`(y_2, y_1) = (1656, 1171)` is a solution. Multiplying `y_2 + y_1 sqrt2` by `(3+2 sqrt2)^2` preserves
the norm and the parity of `y_1`, so there are infinitely many solutions. The six points are distinct.
Since `y_2^2 = 4a + 140` and `R = sqrt(a^2+144)`,

    arc = 2R atan(y_2/(a-2)),   C = arc/sqrt R -> 4   (from above).

Exact values along the orbit (`check_families.py`):

| `a` | `C` |
|---:|---:|
| 685549 | 4.000106 |
| 23290241 | 4.0000031 |
| 791184349 | 4.00000009 |
| 26876979329 | 4.0000000027 |

So `M(C) >= 6` for every `C > 4`. The search `axis/axis.c` finds 12,500,186 six-point axis clusters
with `y_1 <= 10^8`. These are the exhaustive scan's 6-point records (`C = 4.0001` at
`n = 469977431545`).

**3.3 Four consecutive abscissae are impossible.** Among four consecutive integers the squares
represent all of `{0, 1, 4}` mod 8. But `n`, `n-1`, `n-4` cannot all be squares mod 8, since
`{0,1,4} ∩ {1,2,5} ∩ {4,5,0} = ∅` (checked in `check_families.py`). This explains the absence of `m = 4`
solutions with `D = (0,1,2,3)`: none among the 12.5 million six-point clusters.

*Local solubility of depth sets* (`axis/local.py`, moduli `16, 32, 9, 27, 5, 25, 7, 49, 11, 13, 17`;
every insolubility is a proof of non-existence):

| `m` | locally soluble depth sets with smallest `dmax` | limiting `C` | points |
|---:|---|---:|---:|
| 3 | `(0,1,2)` | 4 | 6 |
| 4 | `(0,1,3,4)` only | `4 sqrt2 = 5.657` | 8 |
| 5 | `(0,1,2,5,8)`, `(0,1,4,7,8)`, `(0,3,6,7,8)` | 8 | 10 |
| 6 | `(0,1,2,5,8,9)`, `(0,1,4,7,8,9)` | `6 sqrt2 = 8.485` | 12 |

For the four-depth sets with `dmax = 4`, the others die mod 8 or 16, e.g. `(0,1,2,4)` dies mod 16.
Every 5-set with `dmax <= 7` is insoluble; for example `(0,1,3,4,x)` dies mod 8 or 16 for every
`x <= 7`.

Since `y_dmax^2 >= 2 dmax x_max - dmax^2` and `x_max = sqrt n - O(C^2)`, this gives the following.

**Corollary S.** A conjugation-symmetric cluster with `R -> infinity` has `C >= 4 - o(1)` for 6 points,
`C >= 4 sqrt2 - o(1)` for 8 points, `C >= 8 - o(1)` for 10 points (or 9 points including the axis
point), and `C >= 6 sqrt2 - o(1)` for 12 points.

**3.4 Eight points.** The golden family of Section 2 is conjugation-symmetric: complementing a word
maps the odd class to itself and `F` to `conj F`. Its depth set is `(0,1,3,4)`, so its points lie on
the del Pezzo surface `X_{0134}`:

    3y_1^2 - y_3^2 - 2y_0^2 = 6,   4y_1^2 - y_4^2 - 3y_0^2 = 12.

`axis/axis3.c` searches all integral points with `a <= A` and `y_0 <= 3 sqrt a` (results in
`axis/s_134.txt`). For `a <= 2*10^10` it finds exactly six (Section 9):

| `a` | `y_0` | `C` | source |
|---:|---:|---:|---|
| 74 | 7 | 5.92 (exact arc) | golden |
| 23186 | 127 | 5.8977 | golden |
| 7465178 | 2279 | 5.8977 | golden |
| 2403763490 | 40895 | 5.8977 | golden |
| 5743586 | 527 | **5.6739** | non-golden |
| 590749322 | 12809 | **5.7542** | non-golden |

The golden points recur with ratio `phi^12`. The non-golden points have smaller constants, approaching
the bound `4 sqrt2` of Corollary S.

Whether `X_{0134}` has infinitely many integral points with `y_0^2/a -> 0` is open. It would give
`M(C) >= 8` for every `C > 4 sqrt2`. A fixed `y_0` gives only a genus-1 curve, so finitely many points.

For the diagonal reflection `z -> i conj z` on a primitive (odd) circle, the lines `x+y = u` must have
`u` odd. The same argument then gives `C >= sqrt24 * 2^(1/4) - o(1) = 5.83 - o(1)` for 8 points.

**3.5 Ten points.** Integral-point searches on the three K3 surfaces with `dmax = 8` are in
`axis/s_1258.txt`, `s_1478.txt`, `s_3678.txt`; results are in Section 9. The heuristic count of
Section 7 predicts at most finitely many. Any such point would be a 10-point cluster with `C` near 8,
beating 8 points.

## 4. Polynomial templates in Z[i][t]

**Lemma P.** Let `z_j in Z[i][u]` be monic of degree `D` with `z_j conj(z_j) = N(u)`, and let `S_j` be
the root multisets. Real parts of the power sums `p_r(S_j)` do not depend on `j`. The cluster
`{z_j(u)}` (`u -> infinity`) has bounded `C` iff `Im p_r(S_j)` agree for `r < D/2`. Then

    C = (2/D) * spread_j Im p_{D/2}(S_j) / sqrt|g|,

with `g` the eventual cluster gcd. `C -> 0` iff the agreement extends to `r <= floor(D/2)`; for odd
`D` that is automatic once the bounded-`C` condition holds.

*Proof.* Newton's identities: `z_j - z_l` has degree `D - r_0`, where `r_0` is the first index with
`p_r` differing. The normalised chord is the leading coefficient of `z_j - z_l` at degree `D/2`, which
is purely imaginary by the norm identity. ∎

**Exhaustive search** (`poly/pte.c`): roots `x+iy`, `0 <= x <= X`, `1 <= y <= Y`, multisets of size `D`,
translation-normalised. Largest fibres:

| `D` | box | bounded `C` | `C -> 0` |
|---:|---|---|---:|
| 2 | `X=6, Y=3` | 4 | — |
| 3 | `X=6, Y=3` | — | 3 |
| 4 | `X=6, Y=3` | **6** (base `C` = 8; CG type, all `y = 1`) | — |
| 5 | `X=14, Y=5` (6.4M multisets) | — | 3 |
| 6 | `X=14, Y=4` (27M multisets) | 5 (base 48) | 3 |
| 7 | `X=6, Y=3` | — | 3 |
| 8 | `X=7, Y=2` | 4 | 3 |

After the gcd `(1+i)^2`, the CG family has `C -> 4 sqrt2` (verified at `u = 10^6, 10^12+1`). So
polynomial templates are the weakest mechanism found: 6 points at bounded `C` and 3 with `C -> 0`.
This is consistent with the notes' fixed-shift argument (at most 6 even, at most 3 odd).

## 5. Theorem A: fixed-shift quadratic templates have at most 4 points with C -> 0

**Theorem A.** Take a template as in `quadratic_unit_template_bound.md`, with an unbounded index
sequence and `R_n -> infinity`. If its distinct eventual points lie in arcs of length `o(sqrt R_n)`,
there are at most 4 of them. Equality occurs (canonical four-point Pell family of the notes).

*Proof.* By (9) and (14) of that note, `R_n ≍ lambda^(K_0 n)`, and the first nonzero order `nu` of a
pair gives angle `≍ lambda^(-2 nu n)`. So the pair's normalised chord is `≍ lambda^((K_0/2 - 2nu) n)`,
and `C -> 0` forces `K_0 < 4 nu_jl` for every pair.

*Case 1: some pair has `nu = 1`.* Then `K_0 <= 3`, and (24) of the note gives at most
`2^(K_0-1) <= 4` rows.

*Case 2: every pair has `nu >= 2`.* By (21), `nu <= 3`, so `K_0 <= 11`, and every pair difference
`P` satisfies `P(q) = 0`. Section 4 of the note rules out `lambda > sqrt11`, which covers
`Q(sqrt3)`. Otherwise `||P||_1` is even (26) and at least `lambda^2 + 1` (27). Suppose there are 5 rows.
Each integer threshold of a coordinate separates at most 6 pairs, so

    sum over pairs of ||P||_1 <= 6 K_0 <= 66.

- If `lambda != phi`: every unit has `lambda^2 + 1 > 6.8`. That is `1+sqrt2` (norm `-1`), or
  `lambda >= 3.30` (norm `-1`, trace `>= 3`), or `lambda >= phi^2` (norm 1). So every `||P||_1 >= 8`
  and the sum is at least 80.
- If `lambda = phi`: distance-6 pairs are `±X^s T` by (33). They form an integer-grid graph with at
  most `f(5) = 5` edges (34), and all other distances are at least 8. The sum is at least
  `30 + 40 = 70`.

Both contradict 66. ∎

So the cyclotomic five-row phenomenon (`C -> 0`) genuinely needs growing (resonant) rates. Fixed
shifts reach at most 4 points with `C -> 0`, 6 points for `C < 2 sqrt2` (notes), and 8 points
overall (notes, attained at `C = 5.8977` by Section 2).

## 6. Data: exhaustive scan of all circles with R^2 <= 10^12

**Method** (`families/scan2`, Rust; adapted from `families/scan`).
- *Reduction.* Dividing a cluster by its Gaussian gcd decreases `C` and removes 2 and all inert
  primes. So it suffices to scan `n` whose prime factors are all `1 mod 4`.
- *Representation.* For such `n` the `d(n)` point classes are stored as angles mod `pi/2`, in `u64`
  with `2^64 = pi/2`, built multiplicatively over the factorisation.
- *Part A:* a DFS over `n` with all prime factors at most `10^6`.
- *Part B:* `n = m p` with `p > 10^6` prime. Here `m = 1` and `m` prime are skipped: for `n = qp`, any
  3 of the 4 classes contain a coprime pair, and the pair inequality gives
  `C^4 >= |u_1 - u_2|^2 |u_2 - u_4|^2 >= 4`, so `C > sqrt2`.
- *Measurement.* For each `n` and each `k <= 24`, the scan records `C_k(n) = n^(1/4) * (least span of
  k consecutive angles)` whenever it is at most 6.
- *Validation and accuracy.* Output agrees with a brute-force enumeration of every lattice point for
  `n <= 10^6`. The worst angle error (about `10^(-14)`) is negligible against spans of at least
  `3*10^(-4)`.
- *Scope caveat.* Entries with `k <= 4` and `C > sqrt2` may omit `n = qp` with `p > 10^6`; this
  affects no row of the tables.
- *Cost.* `1.96*10^10` circles in 851 s on 4 threads.

**Largest cluster versus N** (all `n <= N`):

| `N` | C <= 0.5 | C <= 1 | C <= sqrt2 | C <= 2 | C <= 3 | C <= 4 |
|---:|---:|---:|---:|---:|---:|---:|
| 2.6e2 | 2 | 2 | 2 | 3 | 4 | 4 |
| 6.6e4 | 2 | 2 | 3 | 4 | 4 | 5 |
| 1.0e6 | 2 | 3 | 3 | 4 | 4 | 5 |
| 1.7e7 | 2 | 3 | 4 | 4 | 5 | 5 |
| 2.7e8 | 2 | 3 | 4 | 4 | 5 | 5 |
| 4.3e9 | 3 | 4 | 4 | 4 | 5 | 5 |
| 6.9e10 | 3 | 4 | 4 | 4 | 5 | 5 |
| 1.0e12 | **3** | **4** | 4 | 5 | 5 | 5 |

The scan records only `k >= 3`. Pairs (`k = 2`) are attained at every threshold once `N >= 85`: the pair `(x, x+1)`, `(x+1, x)` on `n = 2x^2+2x+1` has chord `sqrt2`. An entry 2 means only pairs are present.

First occurrences:
- 3 points: `C <= 1/2` at `n = 276108509`, `C <= 1` at `128605`, `C <= sqrt2` at `4205`.
- 4 points: `C <= 1` at `412333805`, `C <= sqrt2` at `1762645`.
- 5 points: `C <= 2` at `683617125425`, `C <= 3` at `10414625`.

**Best constant per k over n <= N:**

| `N` | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 | k>=9 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e6 | 0.8115 | 1.6273 | 3.0061 | 4.1323 | 5.5479 | 5.9148 | none <= 6 |
| 1e8 | 0.5669 | 1.0100 | 2.4170 | 4.0067 | 5.5479 | 5.9148 | none |
| 1e10 | 0.3759 | 0.7001 | 2.2584 | 4.0002 | 5.2699 | 5.8978 | none |
| 1e12 | 0.2524 | 0.5060 | 1.9962 | 4.0001 | 5.2699 | 5.8978 | none |

The record minima of `C_3` decay like `n^(-1/12)` (Jarník's exponent). The `C_4` records also
decrease, and `C_5` decreases slowly (3.94 at `n = 325`, 1.996 at `6.8*10^11`). `C_6` is pinned just
above 4 by the axis family. The `k = 7, 8` records are a sporadic cluster (5.2699, `n = 235370525`,
not symmetric) and the golden member `n = 537606725`. Below `C = 4` the scan finds 640 circles with
5-point clusters and none with 6.

## 7. Heuristics: when is the expected number of configurations finite?

**Random model.** On a circle `n` with `d = d(n)` classes, the expected number of `k`-sets in an arc of
constant `C` is about `d (2Cd/pi)^(k-1) n^(-(k-1)/4)/(k-1)!`. Summing over `n <= N` (the means of
`d(n)^k` are polylogarithmic) gives `N^(1-(k-1)/4)` up to logarithms. So:
- `k <= 4`: infinitely many, polynomially;
- `k = 5`: log-divergent;
- `k >= 6`: finitely many for every fixed `C`.

**Balanced simplex profile** (handoff normalisation, `m = k-1` non-anchor rows, `2^m - 1` blocks of
log-norm `w`, `W = (2^m - 1)w`). There are `e^((2^m-1)w)` block choices and the `m` independent phase
windows have probability `(C e^(-W/4))^m`. The expected number per scale is
`C^m exp((2^m-1)(1-m/4) w)`. This gives the same threshold: finite for `k >= 6`, borderline for
`k = 5`.

**Symmetric model.** For `2m` points, a pair `(a, y_0)` with `y_0 <= K a^(1/2)` gives `N^(3/4)`
candidates, and each further depth is a square with probability about `N^(-1/4)`. The count is
`N^((4-m)/4)`:
- 6 points: about `N^(1/4)`, matching the linear growth in `y_1` of the 12.5M solutions;
- 8 points: log-divergent, matching the handful of points on `X_{0134}` up to `n ≈ 10^20`;
- 10 points: finite;
- 12 or more: finite, and under Bombieri–Lang plus the absence of special curves, finitely many
  altogether.

**How algebraic families beat the count.**
- Symmetry halves the number of conditions (6 and 8 points).
- Phase locking makes all conditions hold identically along a family, at the price of row-count
  rigidity (at most 8; at most 4 with `C -> 0`).
- Resonance gives 5 with `C -> 0`.

Every known mechanism has a uniform cap: 8 at bounded `C` and 5 with `C -> 0`. The only infinite
sources of 8-point clusters are log-divergent, and their constants are at least `4 sqrt2 - o(1)`
(Corollary S) or at least 5.27 (sporadic records).

**Plausibility and likely values.** A uniform bound is plausible. For `C <= 4` the predicted picture
is:
- infinitely many 5-point clusters at every `C` (random log-divergence plus the cyclotomic family);
- finitely many 6-point clusters, none found;
- `M(C) = 5` for `C <= 4`, in particular `M(1/2) = M(1) = 5`.

For `4 < C < 5.27`: `M(C) >= 6`, with 7 only from (not yet found) sporadic clusters. For
`5.67 <= C < 8`: `M(C) >= 8`; 9 or more would need sporadic configurations or asymmetric families,
and none with `C <= 6` occur up to `10^12`. These values are conjectural beyond the stated lower
bounds.

## 8. What is rigorous, what is missing, and why the method stalls

**Rigorous here:**
- Proposition G (golden 8 and 7);
- the axis six-point Pell family (`M(C) >= 6` for `C > 4`);
- the exact record clusters (`M(5.2700) >= 7`, `M(5.6740) >= 8`);
- Theorem A;
- the mod-8 lemma and Corollary S (finite congruence computations plus an elementary inequality);
- Lemma P;
- the finite searches and the exhaustive scan, as stated computations.

**Not proved:**
- an upper bound for `M(C)`;
- infinitude of the integral points of `X_{0134}` with `y_0^2/a -> 0`;
- non-existence of 9- or 10-point families at bounded `C`. For 10 points, the K3 surfaces of
  Section 3.5 and the one-class cyclotomic model with nullity at least 5 remain open; the notes'
  bound there is 2251 rows.

**Missing lemma for the uniform theorem.** Unchanged from the handoff: a lower bound `H >= c w` (or
merely `H_sum >= f(w)` with `f(w)/log w -> infinity`) for the residual height of extracted
balanced-profile tuples. The constructions here sharpen what such a lemma must tolerate:

- a Pell family with six points at `C -> 4`, and phase-locked families with eight points at `C = 5.9`
  and residual height `O(1)`;
- infinitely many symmetric 8-point clusters are expected (log-divergent), all at `C >= 4 sqrt2`.

The method of this route (families plus data) cannot prove an upper bound. The surfaces controlling
symmetric clusters (del Pezzo for 8 points, K3 for 10) are exactly where integral points are neither
provably finite nor provably infinite. Asymmetric clusters have no comparable algebraic model, and
their finiteness for `k >= 6` is only heuristic.

## 9. Ten-point and eight-point surface searches (final status)

`axis/axis3.c` enumerates every integral point with `a <= 2*10^10` (so `R^2 = n <= 4*10^20`) and
`y_0 <= 3 sqrt a`, i.e. every symmetric cluster of these depth types with limiting `C <= 2 sqrt(9 + 2 dmax)`.
It uses about `1.5*10^10` candidate pairs per depth set, all four runs completing. Logs are in
`axis/s_*.err`.

**`D = (0,1,3,4)`, 8 points:** exactly 6 points.

| `a` | `y_0` | `C` | source |
|---:|---:|---:|---|
| 74 | 7 | 5.9148 (exact arc) | golden |
| 23186 | 127 | 5.8978 | golden |
| 7465178 | 2279 | 5.8977 | golden |
| 2403763490 | 40895 | 5.89770896 | golden |
| 5743586 | 527 | 5.673925 | non-golden |
| 590749322 | 12809 | 5.754210 | non-golden |

All are verified exactly (primitive, 8 distinct points) in `check_families.py`.

**`D = (0,1,2,5,8)`, `(0,1,4,7,8)`, `(0,3,6,7,8)`, 10 points:** **no** integral points in the same
range. So there is no symmetric 10-point cluster with `C <= 2 sqrt(9+16) = 10` (asymptotically) and
`R^2 <= 4*10^20` of these minimal depth types. This is consistent with the convergent heuristic count
of Section 7.

## 10. Files

- `families/check_families.py`: independent exact verification of Sections 2–3 and the records.
- `families/scan2/`: exhaustive scanner (`src/main.rs`), output `out_1e12.txt`, log `err_1e12.txt`.
  The earlier scanner and outputs to `10^11` are in `families/scan/`, with a brute-force cross-check
  in `families/py/`.
- `families/poly/`:
  - `gpoly.py`: exact `Z[i]` helpers;
  - `pte.c`: polynomial-template fibres;
  - `golden.py`, `golden_search_E4.py`, `golden_nu3*.py`: golden templates;
  - `silver.py`, `units_search.py`: other units;
  - `pell8.py`, `fib8.py`: family evaluations;
  - `cluster_at.py`: extracts the best `k`-cluster of a given circle and its prime-orientation pattern.
- `families/axis/`:
  - `axis.c`: three consecutive abscissae;
  - `axis2.c`, `axis3.c`: general depth sets;
  - `local.py`: local solubility;
  - result files.
- `families/notes_checks/`: re-run of the notes' five-row checker
  (`k=19`, `d=1`: `7 chord^4/R^2 < 1`; `d=3`: the same).
