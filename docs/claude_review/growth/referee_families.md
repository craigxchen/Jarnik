# Referee report: "Constructions and data: how large is M(C)?" (growth/families.md)

Status of the submission: data_or_construction (no general upper bound is claimed).
Verdict: **the claimed results survive.** Several statements need small corrections:
- a proof gap in Corollary S, which I closed by computation;
- mislabelled rows in a data table;
- a search restriction missing from the summary;
- an incomplete symmetric search. It missed a primitive 10-point cluster with `C = 8.2393`, which I found; see the
  addendum to §6. Novelty
relative to the project notes is genuine for items (1), (2), (4) and (5). I could not complete the
literature check because arXiv and the authors' home pages are blocked from this sandbox.

All referee scripts are under `growth/referee_families_checks/`. They were written from scratch and do not
import `families/check_families.py`:

- `ref_check.py`: exact Gaussian arithmetic. It checks the golden family, using a Decimal arcsine with argument reduction
  at 110 to 1000 digits, and verifies distinctness *after* the unit rotation into the arc.
- `ref_check2.py`: the notes' translated Pell eight-point family.
- `ref_check3.py`: the axis six-point family and all record clusters.
- `ref_local.py`: local solubility of depth sets, for both reflection types.
- `brute.c`: brute-force enumeration of every lattice point of every circle `n <= N`, with no factorization and no gcd
  reduction. Outputs are `brute_1e9.txt` and `brute_1e10.txt`.
- `cyc5_copy.py`: a verbatim copy of the notes' cyclotomic five-row checker.

---

## 1. Golden eight-point family: CONFIRMED

**Construction.** Blocks `J_r = f_{r+1} + i f_r`. Take the 8 odd words on `J_n..J_{n+3}`, with prefactor `F = 1+2i` when `p = 1`
and `conj F` when `p = 3`.

**Exact independent recomputation.** For every `n = 1..13`, plus `n = 100` and `n = 301`:
- the eight points are distinct lattice points after unit rotation into one arc;
- they have equal norm;
- the cluster gcd has norm exactly 4 when `n = 1 (mod 3)`, and 2 otherwise.

Values of `C_8`:

| `n` | `C_8` |
|---|---|
| 1 | 5.914765 (`n = 5525`) |
| 4 | 5.897763 (`N = 537606725`, the scan record) |
| 7 | 5.8977091305 |
| 10 | 5.89770896341 |
| 100 | `C8 + 3e-85` |
| 301 | `C8 + 2e-109` (limited by working precision) |

The sequence decreases to `C8 = sqrt2 (4+sqrt5)/5^{1/4} = 5.8977089628938...` **from above**. So "M(C) >= 8 for every
C > C8" is correct. It already follows from the explicit member `n = 7`, which has `C < 5.8978`. For `n` not `1 (mod 3)`
the limit is `2^{1/4} C8 = 7.0136`, as predicted (only one block has even norm). The 7-point sub-cluster converges to
`C7 = sqrt2 phi^2 5^{1/4} = 5.536467812830`.

**Proposition G, checked line by line.**
- `tau(a)/a = -(psi+i)/(phi+i) = -i/phi`. This uses `phi^2+1 = sqrt5 phi`.
- `a^2` is proportional to `phi(1+2i)`, so `2 arg a = arg F`.
- The leading phases cancel in both parity classes.
- The odd-word extremes are `1+q^2+q^3` (word 1011) and `q` (word 0100), with spread `2 sqrt5 - 3`.
- `R^2 ~ phi^{8n+16}/20`.
- `C = 2 phi^3 (2 sqrt5 - 3)/20^{1/4}`.

All correct.

**Comparison constant.** The notes' family (`primitive_eight_point_translated_pell_family.md`, n = 61) has gcd 1 and
`C = 44652598.0015640`. This reproduces the author's value, so the ratio of about `7.6e6` is right.

**Remark on novelty and mechanism.** The Pell blocks satisfy `Norm(H_r) = f_{12r+1} = Norm(J_{6r})` (checked at r = 5).
So the notes' family is the same golden template with block spacing 12, and the new family uses spacing 1. That is a
natural refinement, but no explicit golden-unit eight-point family, and no bounded-C eight-point cluster with a small
constant, appears in the notes. I grepped for the record radii, for `5.897`, golden/Fibonacci eight-point families,
items 511–512, `pell_product_endpoint_obstruction.md` and `quadratic_unit_template_bound.md`. The notes never optimised
C, so the factor `7.6e6` is a comparison with an unoptimised constant, not a measure of difficulty.

**Correction of emphasis.** For the *uniform* `M(C)` of the handoff (supremum over all R), the single exact clusters
of item (3) supersede the family statements for 7 and 8 points:
- `M(C) >= 8` already holds for `C >= 5.6740`;
- `M(C) >= 7` already holds for `C >= 5.2700`.

The family statements (`C > 5.8978` and `C > 5.5365`) matter only for the asymptotic quantity
`M_inf(C) = limsup_{R -> inf}`. There they are the best known bounds, because the non-golden (0,1,3,4) points are isolated.
The write-up should state which of the two quantities each bound refers to.

## 2. Axis six-point Pell family: CONFIRMED (elementary, exact)

The identities `(a-1)^2 + y1^2 = a^2 + 144 <=> a = (y1^2-143)/2` and `(a-2)^2 + y2^2 = a^2+144 <=> y2^2 - 2y1^2 = -146`
are correct.

The parity of `y1` is preserved already by one application of `3 + 2sqrt2`. The new value is `y1' = 2y2 + 3y1`, so the
write-up's use of the square is harmless overkill. I checked 8 consecutive orbit members from `(1171,1656)` up to
`a = 3.6e16`. In each, all six points lie on the circle, the cluster gcd is 1 (primitive), and `C - 4` takes the values
`1.06e-4, 3.12e-6, ..., 2.03e-15`, always positive.

Asymptotically, `C = 4(1 + 18.17/a + ...)`. This matches the first member: `4 * 18.17/685549 = 1.06e-4`.

The fundamental solutions with odd `y1` include `(9,4)`, `(19,24)`, `(35,48)`, `(105,148)`, `(201,284)` and `(611,864)`.
The member at `n = 405176785` (`C = 4.003609`) is the best 6-point cluster below `1e9` in my brute-force scan.

**Verdict.** `M(C) >= 6` for every `C > 4` is rigorous, with clusters in `M_inf` as well.

**Novelty.**
- Relative to the notes: the notes' best six-point constant is `4 sqrt2`. Two sources give it:
  - Cilleruelo–Granville's six points, recorded as 8 in `codex_uniform_bound_research.md` §7, which is `4 sqrt2` after
    the gcd `(1+i)^2`;
  - `full_clipping_six_point_pell_sharpness.md`.

  The write-up does not cite the second family. It is itself a conjugation-symmetric cluster, with depth set (0,1,3) and
  `y0 ~ sqrt(2a)`, so `C^2 -> 24 + 8 = 32`. The new family improves `4 sqrt2` to 4, and that is new for the project.
- Relative to the literature: I could not access the Cilleruelo–Granville papers. The construction (three consecutive
  abscissae via a Pell conic) is elementary and may be folklore. No searchable source states it.

## 3. Single exact clusters: CONFIRMED

Each cluster was recomputed from scratch; all are primitive (gcd norm 1).

| Cluster | Circle `n = R^2` | `C` | Notes |
|---|---|---|---|
| 8 points, a = 5743586, depths (0,1,3,4) | 32988780417125 | 5.673924770 | |
| 8 points, a = 590749322 | 348984761607530165 | 5.754209881 | |
| 7 points | 235370525 | 5.269927439 | circle has 72 first-quadrant points |
| 5 points | 683617125425 | 1.996248844 | |

Hence `M(5.6740) >= 8` and `M(5.2700) >= 7` hold rigorously.

## 4. Theorem A (fixed-shift quadratic templates: at most 4 points with C -> 0): CONFIRMED

I re-derived the proof against `quadratic_unit_template_bound.md`.

**Reduction.** By (9) and (14), the normalised chord of a pair is `≍ lambda^{(K0/2 - 2nu) n}`. So `C -> 0` iff
`K0 < 4 nu` for every pair. Measuring C with a sub-cluster's own reduced radius only lowers it, so "sub-cluster" is
harmless.

**Case nu = 1 for some pair.** Then `K0 <= 3`, and (24) gives at most 4 rows.

**Case nu >= 2 for all pairs.** First-order vanishing means `(eta - conj eta) P(q) = 0`, so `P(q) = 0`. By (21),
`nu <= 3`, hence `K0 <= 11`. Section 4 gives `lambda^2 + 1 <= ||P||_1 <= K0 <= 11`, which excludes `Q(sqrt3)`
(`(2+sqrt3)^2 = 13.9`). So parity (26) applies.

**Five rows give a contradiction.**
- A threshold separates at most `2*3 = 6` pairs, so the sum of pair distances is at most `6 K0 <= 66`.
- If `lambda != phi`: the smallest units are `1+sqrt2`, with `lambda^2 + 1 = 6.83`, and `phi^2`, with 7.85. Parity then
  forces `||P||_1 >= 8`, so the sum is at least `80 > 66`.
- If `lambda = phi`: (33) classifies distance-6 differences as `±X^s T`. These form an integer-grid graph with at most
  `f(5) = 5` edges (34), and every other distance is at least 8. So the sum is at least `30 + 40 = 70 > 66`.

**Sharpness.** The canonical four-point Pell family (`canonical_four_point_counterexample.md`, `K0 = 3`, C -> 0)
attains 4.

**Novelty.** The notes have at most 6 rows for `C < 2 sqrt2` (`multiplicative_rectangle_separation.md`) and at most 5
rows for `lambda > sqrt11` and `10 C^2 <= 8` (`small_width_affine_circuits.md`). Theorem A (at most 4 as C -> 0, for
every unit) is new.

**Scope.** It is a statement about one template class only. The informal sentence "five points with C -> 0 need growing
(resonant) rates" should read "... cannot come from a fixed-shift quadratic-unit template".

## 5. Reflection-symmetric clusters (Corollary S): TRUE, but the write-up proves only part of it

**Mod-8 lemma.** Four consecutive x give squares `{0,1,4}` mod 8. Then `{0,1,4} ∩ {1,2,5} ∩ {4,5,0} = ∅`. Correct.

**Geometric inequality.** It is exact up to `o(1)`. The arc is at least the chord `2 y_dmax`, and
`y_dmax^2 = y0^2 + 2a dmax - dmax^2`. Since `a = sqrt n - O(C^2)`, this gives `C^2 >= 8 dmax (1 - O(C^2/sqrt n)) - O(1/sqrt n)`.

**Axis type.** My own local-solubility test used moduli `16, 32, 64, 128, 9, 27, 81, 5, 25, 125, 7, 49, 11, 13, ..., 41`.
It reproduces `local.py` exactly:

| `m` | smallest locally soluble `dmax` | depth sets | bound |
|---|---|---|---|
| 3 | 2 | (0,1,2) | 4 |
| 4 | 4 | (0,1,3,4) only | `4 sqrt2` |
| 5 | 8 | (0,1,2,5,8), (0,1,4,7,8), (0,3,6,7,8) | 8 |
| 6 | 9 | (0,1,2,5,8,9), (0,1,4,7,8,9) | `6 sqrt2` |

**Gap.** Corollary S is stated for clusters closed under *a lattice reflection*. Removing the cluster gcd can turn an
axis-symmetric cluster into a diagonal-symmetric one, since `conj(1+i) = -i(1+i)`. But §3.4 treats the diagonal
reflection `z -> ±i conj z` only for 8 points. For 10 and 12 points the axis argument does **not** transfer
automatically. In coordinates `u = x+y`, `v = x-y` on `u^2 + v^2 = 2n` (n odd, u and v odd),
`C_xy = C_uv/2^{1/4}`. The trivial spans give only:
- 10 points: `dmax_uv >= 8`, so `C >= 6.73`, which is less than 8;
- 12 points: `dmax_uv >= 10`, so `C >= 7.52`, which is less than `6 sqrt2`.

**Closing the gap** (`ref_local.py`, diagonal case):

| points | smallest locally soluble `dmax_uv` | depth sets | resulting bound |
|---|---|---|---|
| 6 | 4 | (0,2,4) | `C >= 4.757` |
| 8 | 6 | (0,2,4,6) | `C >= 5.826`, the author's value |
| 10 | 12 | six sets, e.g. (0,2,4,6,12) | `C >= 8.239 > 8` |
| 12 | 14 | (0,2,6,8,12,14) only | `C >= 8.899 > 6 sqrt2` |

Every insolubility is a finite congruence certificate. A diagonal axis point would need `x = y`, hence `n = 2x^2` even,
which is impossible after reduction. So **Corollary S holds as stated for both reflection types**, but the proof needs
this extra case. Add it to the write-up.

**Further remarks on §3.**
- The surface table (quadric, quartic del Pezzo, K3, general type) is the canonical-class count `K = (m-5)H` for a
  complete intersection of `m-2` quadrics in `P^m`. It presumes smoothness or mild singularities, which was not checked.
  For integral points the relevant pair is `(X, D = hyperplane at infinity)`, with `K+D = (m-4)H`. This is consistent
  with the heuristic exponent `(4-m)/4`.
- The search `axis3.c` is complete for `a <= 2*10^10` **only under its restriction `y0 <= 3 sqrt a`**, i.e. limiting
  `C <= 2 sqrt(9+2 dmax)`: about 8.2 for (0,1,3,4) and 10 for `dmax = 8`. `families.md` §9 says so. The submitted
  summary ("Up to R^2 = 4e20 there are exactly six (0,1,3,4) points") omits the restriction and should include it.
  I checked the search logic: the factorisation `(y_{d1} - y0)(y_{d1} + y0) = 2 d1 a - d1^2`, the `u` and `v` loops, and
  the break condition.

## 6. Data: CONFIRMED, with one mislabelled table

**Independent brute force.** `brute.c` enumerates every lattice point directly, uses no reduction, and wraps the angle
list as often as needed, so small circles are handled.

Up to `n <= 10^9` it reproduces every first occurrence:

| Event | First `n` |
|---|---|
| 3 points, `C <= 1/2` | 276108509 |
| 3 points, `C <= 1` | 128605 |
| 3 points, `C <= sqrt2` | 4205 |
| 4 points, `C <= 1` | 412333805 |
| 4 points, `C <= sqrt2` | 1762645 |

It also reproduces the records:
- `k = 7`: 5.269927 at 235370525;
- `k = 8`: 5.897763 at 537606725;
- `k = 9`: none below 7.99, the best being at `n = 5525`.

So there is no cluster of 9 or more points with `C <= 6` up to `10^9`.

**Mislabelled table.** The "best constant per k over n <= N" table in §6 has wrong labels.
- The row "1e8" is in fact the minimum over `n < 2^26 ≈ 6.7e7`. This can be read off the scanner's dyadic `blockmin`
  lines. For `n <= 1e8` the true `k = 3` value is **0.551921**, at `n = 88027829 = 229*269*1429` with points
  (9305,1202), (9302,1225), (9298,1255). The table gives 0.5669. The other `k` in that row agree with my scan.
- The row "1e10" is likewise `n < 2^33 ≈ 8.6e9`. See the addendum below for the `1e10` scan.
- The rows "1e6" and "1e12" are correct.

This error is cosmetic: no conclusion depends on the intermediate rows.

**Scope of the 1e12 scan.** The Part-B skip of `n = p` and `n = qp` is justified. For `n = qp`, any three classes
contain two pairs sharing one prime, and the unit-adjusted chords are at least `sqrt(2q)` and `sqrt(2p)`. So
`C > sqrt2`, which affects no table entry. I did not re-run the full `1e12` scan, which would take about a day with my
brute force. My independent check stops at `1e10`.

**Cyclotomic family (notes).** The verbatim checker passes. It shows `k = 19`, `d = 1, 3` with
`7 chord^4/R^2 < 1`, i.e. a 5-point instance with `C < 0.62` at `R ≈ 10^3428`. "M(C) >= 5 for all C > 0" rests on the
notes' asymptotic proof that this template has `C -> 0`. I did not re-audit that proof here.

### Addendum: brute force to 1e10, and a 10-point cluster the submission missed

**Scan results.** `brute_1e10.txt` (9 minutes, 4 threads) gives the following best constants for `n <= 10^10`:

| `k` | best `C` | at `n` |
|---|---|---|
| 3 | **0.371418** | 9520005625 |
| 4 | 0.700105 | |
| 5 | 2.258390 | |
| 6 | 4.000197 | |
| 7 | 5.269927 | |
| 8 | 5.897763 | |
| 9 | **7.880380** | 1176852625 |
| 10 | **8.239327** | 1176852625 |

- The `k = 3` value at `n = 9520005625 = 5^4*13*673*1741` confirms the second mislabelled row: the table's "1e10"
  entry 0.3759 is the minimum over `n < 2^33`.
- The scanner keeps only one minimum per dyadic block. When that minimum lies above `N`, building the "n <= N" rows from
  block minima overstates the true value.
- The `k = 4..8` entries agree with the table.

**New cluster: 10 points with C = 8.2393, missed by the submission.** It lies on `n = 1176852625 = 5^3*13^2*17*29*113`.
The ten points are:

| `x` | 24791 | 24745 | 24636 | 24567 | 24260 |
|---|---|---|---|---|---|
| `y` | 23712 | 23760 | 23873 | 23944 | 24255 |

together with their reflections `(y, x)`.

- The cluster is primitive (gcd norm 1).
- It contains no further lattice points in its arc.
- Its exact constant is `C = 8.239326866`.
- Dropping one end point gives 9 points with `C = 7.880379573`.

Hence, rigorously, **`M(7.8804) >= 9` and `M(8.2394) >= 10`**.

**Why it was missed.** The cluster is symmetric under the *diagonal* reflection `z -> i conj z`. Its `u = x+y` depth set
is `(0,4,6,10,12)`, which is one of the minimal locally soluble diagonal sets from §5 above, and it has `v0 = 5`. So it
essentially attains the diagonal 10-point bound `sqrt(96)/2^{1/4} = 8.239069`.

The submission searched only the axis-type 10-point K3 surfaces, and its scan records only `C <= 6`. So its summary
sentences need correcting:
- "There are no 10-point symmetric clusters of the minimal depth types" is true for the axis type only.
- "All known mechanisms stop at 8 (bounded C)" is false as a statement about configurations: the symmetric mechanism
  produced 10 points at `C = 8.24`.

**Is it sporadic?** I ran the author's own `axis3.c` (copied) on all six minimal diagonal depth sets with `dmax_uv = 12`,
for `a <= 10^10` and `v0 <= 3 sqrt a`. That covers `n` up to about `5*10^19`. Exactly one integral point appears, this
one (`diag_*.txt`). So it looks sporadic, consistent with the convergent heuristic count for 10-point symmetric
clusters. No claimed theorem is affected: Corollary S predicts at least 8.239 for diagonal 10-point clusters, and this
cluster meets that bound.

## 7. Heuristics and the "prediction" M(1/2) = M(1) = 5

These are clearly labelled as heuristics in the write-up and are not claimed as proofs. One minor point: the random
model at `k = 5` gives `sum d(n)^5/n`, which diverges like a *power* of `log N`, not like `log N`. The prediction
`M(1/2) = M(1) = 5` is conjectural. It is consistent with all data, but the data cover only `R <= 10^6`, while the fifth
point appears at `R ~ 10^3428`.

## 8. Corrected statement

- **(1)** For `n = 1 (mod 3)`, the eight odd golden words with prefactor `F` or `conj F` are primitive after dividing by
  `(1+i)^2`. Their constant decreases strictly to `C8 = sqrt2 (4+sqrt5)/5^{1/4} = 5.89770896...`, and the 7-point
  sub-cluster decreases to `5.5365...`.
  - Hence `limsup_{R -> inf}`-versions `M_inf(C) >= 8` for `C > C8` and `M_inf(C) >= 7` for `C > 5.5365`.
  - For the uniform `M(C)`, single clusters give more: `M(C) >= 8` for `C >= 5.6740` and `M(C) >= 7` for `C >= 5.2700`.
- **(2)** `M(C) >= 6`, and `M_inf(C) >= 6`, for every `C > 4`, by the elementary Pell family.
- **(3)** Unchanged.
- **(4)** Theorem A is unchanged. It is a restricted statement about fixed-shift quadratic-unit templates and depends on
  `quadratic_unit_template_bound.md`.
- **(5)** For a cluster invariant under any lattice reflection, as `R -> inf`, the bounds are:
  - 6 points: `C >= 4 - o(1)`;
  - 8 points: `C >= 4 sqrt2 - o(1)`, with extremal pattern axis type (0,1,3,4) only; diagonal type needs `5.826`;
  - 10 points: `C >= 8 - o(1)`; diagonal type needs `8.239`;
  - 9 points with an axis point: `C >= 8 - o(1)`;
  - 12 points: `C >= 6 sqrt2 - o(1)`; diagonal type needs `8.899`.

  The diagonal cases for 10 and 12 points need the extra congruence certificates computed here.
- **Data.** As stated, except for three points:
  - the intermediate rows of the best-C table are dyadic;
  - the `X_{0134}` and K3 searches are complete only for `y0 <= 3 sqrt a`;
  - the axis-only search misses a diagonal-symmetric 10-point cluster.
- **Added (referee).** `M(C) >= 9` for `C >= 7.8804` and `M(C) >= 10` for `C >= 8.2394`, from the primitive diagonal
  cluster on `n = 1176852625`. No 9- or 10-point configuration with `C <= 6` is known.

No general upper bound is claimed. None follows, and the write-up says so correctly.
