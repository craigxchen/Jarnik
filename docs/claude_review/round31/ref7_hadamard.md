# Referee report, round 31, Part VII, block "hadamard" (VII.4)

Source: Astra, `astra_free_VII/answer.md`, section VII.4 (lines 1615-1650), summary items 9-10; scripts
`vii_double_core.cpp/.out` (read only after my own runs). Binary model: marks z_i in F_4[t] (theta^2+theta+1=0),
one norm N in F_2[t], equal degree D, critical: deg(z_i - z_j) <= D/2 (c = 0), K = ceil(D/2), U_K = 1 + t^{-K}F_4[[1/t]].

My scripts (own F_4[t]/F_2[t] arithmetic, nothing imported) are in `ref7_hadamard/`:
* `brute_hadamard.cpp` (`bh`): all assignments of monic nonreal blocks f_v of degree <= dmax to a set of Walsh
  columns; full gcd removed; critical test; distinctness; disjoint vs shared base-prime support. Outputs `small.out`,
  `m8_*.out`.
* `twocore.cpp` (`tc`): the eight-row first-core test, in two independent modes: `direct` (every 4-subset of primes,
  derivative of the actual product) and `key` (my own derivation of the A^2+tB^2 key); line and plane counts.
* `twocore_big.cpp` (`tcb`): same, memory-lean, faster line/plane search; used for e = 18 and positive controls.
  Outputs `tc*_*.out`.

## 0. Verdicts in brief

| claim | verdict |
|---|---|
| (a) VII.4 theorem: m >= 8, all rows, proper column subset, disjoint supports => no critical cluster | **CORRECT** (line by line; brute force m = 8 agrees). Threshold m >= 8 is exact: **false at m = 4** (explicit cluster below). Holds in a stronger form without disjointness (§1.9). |
| (b) scope "with VI.2, only the full m-1 columns remain for a single complete set of m >= 8 rows" | **CORRECT**, narrow class (§2) |
| (c) two-core computation e = 10, 12, 14, 16 | **CORRECT**, all counts reproduced exactly (67/10/0, 319/34/0, 4382/226/0, 40546/1453/0); **extended by me to e = 18: also impossible** (840521/44233/0). Not a structural obstruction: a random model fitting the line counts within 3-10% predicts the first admissible planes near e = 22-24 (§3.4, heuristic). |
| (c') criterion (deg H' <= e/2-1; key floor(t^{2e-b-1}B/A)) | **CORRECT** (derived independently; key mode = direct mode at e = 10, 12) |
| (d) growth scale over F_2[t] | **unchanged**: 3M log_2 M - O(M) <= W(M) <= 2M^2 log_2 M |
| (d') Gram/Laplacian step over Z | the step **does not transfer** (no exact sign of the off-diagonal over Z) and its content (rank >= M/2) is **weaker than what the record already has** for integer profiles. Nothing new over Z. |

Novelty: VII.4 is new; it closes the record's open class "incomplete Walsh cores" (README_VI §3(d), ref6_bounds §1.4.1,
§2 open item 1) in the all-rows, block-form case. Significance for the integer problem: **none** (0-1/10). It is
an exclusion of one construction class over F_2[t], and its decisive step (depth multiplication by 2^a) has no
archimedean counterpart.

## 1. (a) The theorem, line by line

Setting after reduction. Marks z_i = prod_{v in R} g_{v,i}, g_{v,i} = f_v if chi_v(i) = +1, conj(f_v) otherwise.
Write g_{v,i} = conj(f_v) * r_v^{(1+chi_v(i))/2}, r_v = f_v / conj(f_v).

1.1 **Reduction.** With disjoint base-prime supports, gcd_i(z_i) = G * prod_v gcd(f_v, conj f_v) (common factor G
of the marks times the self-conjugate gcd in each block), since primes of different blocks cannot interact.
Reducing f_v by g_v = gcd(f_v, conj f_v) leaves f_v/g_v coprime to its conjugate; r_v is unchanged, d_v drops.
Reduced degree D and reduced differences: deg(z_i - z_j) <= D_0/2 - g <= (D_0 - g)/2, so the reduced cluster is
critical (the standing caution "do not discard gcd costs" is respected: division only tightens). A reduced block
that becomes constant is dropped, which only enlarges the missing set; distinctness of marks survives division.
A reduced nonconstant block is nonreal (coprime to its conjugate). Correct.

1.2 **The identity** prod_i z_i^{(chi_v(i) - chi_w(i))/2} = r_v^{m/4} (w missing, v retained). Exponents
e_i in {-1,0,1}, sum zero by balance, so the factor prod conj(f_u) cancels and the product is
prod_u r_u^{sum_i e_i chi_u(i)/2} = prod_u r_u^{(<chi_v,chi_u> - <chi_w,chi_u>)/4}. For u = v: (m - 0)/4. For any
other retained u: u != v, and u != w because w is missing, so both inner products vanish. Hence exactly r_v^{m/4};
4 | m by orthogonality of two balanced columns. Correct. (The missing column is essential: with w retained one gets
(r_v/r_w)^{m/4}, VI.2's identity.)

1.3 **Depth doubling.** U^1 of F_4((1/t)) is torsion-free and (1+x)^{2^a b} - 1 = (x^{2^a})(b + ...) for b odd:
depth multiplied by exactly 2^a. With m = 2^A u: depth(r_v) * 2^{A-2} >= K, i.e. depth(r_v) >= 4uK/m >= 2uD/m.
Correct.

1.4 **d_v >= depth(r_v)**: r_v - 1 = (f_v + conj f_v)/conj f_v, numerator nonzero (nonreal), so
depth = d_v - deg(f_v + conj f_v) <= d_v. Summing, D >= n * 2uD/m, so n <= m/(2u). Correct.

1.5 **Gram matrix.** S_{iv} = chi_v(i) sqrt(d_v). (SS^T)_{ii} = D; row sums sum_v chi_v(i) d_v sum_j chi_v(j) = 0
(balance). Off-diagonal = 2 Agr_ij - D, where Agr_ij = sum of d_v over retained v with chi_v(i) = chi_v(j). The
product P_ij of the common choices divides z_i and z_j literally (it is a sub-product of both), hence divides
z_i - z_j != 0, so Agr_ij = deg P_ij <= deg(z_i - z_j) <= D/2. **No coprimality or disjointness is used here**, so
the question "when reduced blocks share factors" does not arise: divisibility of a common sub-product is formal.
Off-diagonals are <= 0 and rows sum to 0: SS^T = sum_{i<j} w_ij (e_i - e_j)(e_i - e_j)^T, w_ij = -L_ij >= 0, a
weighted Laplacian, nullity = number of components of the graph {w_ij > 0}. Each vertex has diagonal D > 0 and row
sum 0, so it has a positive edge; components have size >= 2, nullity <= m/2, rank >= m/2. Rank = n because the
columns of S are nonzero and orthogonal. Correct. (Remark: rank >= m/2 is the classical "at most 2n vectors with
pairwise nonpositive inner products in R^n"; the Laplacian form additionally yields the pairing.)

1.6 **Consequences.** n >= m/2 and n <= m/(2u): u = 1, n = m/2, nullity = m/2, every component is a pair {i,j}
with L_ij = -D, i.e. Agr_ij = 0: opposite signs in every retained column. Equality throughout 1.4: D even, K = D/2,
d_v = 2D/m = D/n for all v, depth(r_v) = d_v, so deg(f_v + conj f_v) = 0: writing f_v = A_v + theta B_v,
f_v + conj f_v = B_v = 1. Correct.

1.7 **Restriction.** Choose from each pair the row with chi_{w0} = +1. For retained v != v', chi_v chi_v' takes the
same value on both rows of a pair, so its sum over the n chosen rows is half of 0. The n retained columns
restricted to the chosen rows are orthogonal: an n x n Hadamard matrix with w0 as all-ones column. n = m/2 is a power
of two. Correct.

1.8 **Final step.** On the chosen rows, exponents (chi_v - chi_v')/2 (v, v' != w0) are integral with zero sum
(orthogonal to w0), and the same computation as 1.2 gives (r_v / r_v')^{n/4} (w0 contributes <chi_v, 1> = 0). The
chosen marks are still pairwise in U_K with the original K. n/4 is a power of two **only if n >= 4, i.e. m >= 8**;
then depth(r_v/r_v') >= 4K/n = 2D/n = 2d. Numerator f_v conj f_v' + conj f_v f_v' = A_v + A_v' (I recomputed:
(A+theta)(A'+theta^2) + (A+theta^2)(A'+theta) = A + A'), denominator degree 2d; depth = 2d - deg(A_v + A_v')
(valid with or without cancellation). So deg(A_v + A_v') <= 0, A_v in {A, A+1}. Since A + 1 + theta = A + theta^2 =
conj(A + theta), at most two distinct blocks, and the second is the conjugate of the first. n - 1 >= 3 nonconstant
blocks with disjoint supports are pairwise distinct (indeed even two would already be impossible with disjoint
supports, but n - 1 >= 3 is automatic). Contradiction. Correct.

1.9 **Threshold and hypotheses.**
* **m = 4 is a genuine exception** (the theorem's m >= 8 is exactly right). At m = 4, n = 2 and step 1.8 has no
  integral exponent. Explicit critical cluster with a missing column (found by `bh 2 3 1 2`, `small.out`):
  rows of H_4, retained columns chi_1 = (+,+,-,-), chi_2 = (+,-,+,-), missing chi_3;
  f_1 = t^2 + t + theta (norm t^4+t+1, irreducible), f_2 = t^2 + theta = (t + theta^2)^2 (norm (t^2+t+1)^2);
  marks f_1 f_2, f_1 conj f_2, conj f_1 f_2, conj f_1 conj f_2, D = 4, gcd 1, pair differences of degree 2, 2, 1, 1, 2, 2
  (z_1 - z_4 = A_1 + A_2 + 1 = t + 1, z_2 - z_3 = A_1 + A_2 = t). By the same algebra, for m = 4 **every** pair
  f_v = A_v + theta of equal degree with disjoint supports gives a critical 4-cluster. 96 disjoint-support
  assignments with deg <= 3 (plus 24 shared-support ones).
* **Common-factor remark.** Not needed for the logic of 1.2-1.8, which only use (i) the identity, (ii) the cluster
  condition at the degree D = sum d_v of the displayed products, (iii) pairwise distinct nonconstant nonreal blocks.
  It is what makes the theorem cover non-reduced inputs, and with disjoint supports it is harmless because reduction
  preserves the block form (1.1). **Strengthened statement (referee; same proof):** VII.4 holds whenever the sign
  matrix plus all-ones column is Hadamard of order m >= 8, a proper subset of the nonconstant columns is used, the
  f_v are monic, nonconstant, f_v != conj f_v, pairwise distinct with f_v' != conj f_v, and the *undivided* block
  products are a critical cluster at D = sum d_v. Disjointness is used only for distinctness and to make the
  reduction preserve block form — the exact analogue of ref6_bounds' (H1)-(H4) for VI.2. It does **not** cover
  shared supports where dividing by the gcd cuts into blocks.
* Astra's one-line summary (final7 item 9, Part VII summary item 9) drops the disjoint-support/block hypotheses. The
  statement in VII.4 itself is correct; the summary line should carry "in block form with disjoint (or distinct,
  undivided) blocks".
* Standing cautions: points equal-degree and monic (critical forces equal leading coefficients, so unit twists of
  marks are excluded automatically; unit multiples of blocks are absorbed); one norm in F_2[t] (prod of N(f_v));
  conjugates: a row negation is not allowed (it would conjugate a mark), and the proof never uses one; infinity is
  inert (degree-1 place, F_4/F_2); thresholds: m >= 8 is uniform, no m-dependent threshold. All respected.

1.10 **Brute force (m = 8, Sylvester H_8; all normalized H_8 are equivalent to it under row/column permutations and
column negations, the latter covered by enumerating f and conj f).** Retained-column orbits under GL(3,2):

| retained columns | rows distinct | block degrees | assignments | critical clusters (disjoint / shared supports) |
|---|---|---|---|---|
| {1,2,4} (3, independent) | yes | 1..3 | 343,000 | 0 / 0 |
| {1,2,3} (a line) | no | 1..3 | 343,000 | 0 / 0 |
| {1,2,4,7} (complement of a line) | yes | 1..3 | 24,010,000 | 0 / 0 |
| {1,2,4,3} (complement of an independent triple) | yes | 1..3 | 24,010,000 | 0 / 0 |
| {1,...,5} | yes | 1..3 | 1,680,700,000 | 0 / 0 |
| {1,...,6} | yes | 1..2 | 7,529,536 | 0 / 0 |
| control: complete {1,...,7} | yes | 1..2 | 105,413,504 | 0 / 0 |
| control m = 4: {1,2} | yes | 1..3 | 4,900 | 96 / 24 |
| control m = 4: {1,2,3} | yes | 1..3 | 343,000 | 144 / 144 |

Full gcd removed in every case, so these tests include non-common-factor situations with shared primes that the
theorem does not cover. No violation. (Degrees are small, D <= 15; this is a consistency check, not a proof.)

## 2. (b) What VII.4 + VI.2 exclude, and what not

Excluded over F_2[t] (block form, each column a sign-pattern of one block, distinct/disjoint blocks):
* any **single** Hadamard matrix of order m >= 8 with **all m rows** used and a proper subset of its columns
  (VII.4), at any degrees, any per-prime orientations inside a block, any common factor;
* complete cores of non-power-of-2 order, and complete power-of-2 cores below D >= m(m-1)(log_2(m-1)-3)/6
  (VI.2; vacuous for m <= 8, so "Omega(m^2 log m)" is asymptotic).

So for m >= 8 rows forming one Hadamard matrix the only surviving block profile is the complete one, which costs
Theta(m^2 log m) (VI.2 with VI.5). This closes, in this form, the record's open class "incomplete Walsh cores" —
note that ref6_bounds §1.4.1 only had D >~ k log k there; VII.4 shows they do not exist at all for m >= 8.

Not excluded:
1. **Row subsets**: m marks on a subset of rows of a larger Walsh matrix (balance and orthogonality fail on the
   subset; identity 1.2 picks up all blocks).
2. **Several Hadamard matrices on the same rows** (two Walsh cores on x and x^3, or on permuted rows): the columns
   of the union are not mutually orthogonal, so neither 1.2 nor rank = n holds.
3. **Extended/flipped profiles**: a Hadamard core plus any extra column (even one prime outside the core breaks
   1.2 and the Gram balance).
4. **Shared supports with non-common division** (Part V-type designs once the gcd cuts into blocks); Part V itself
   satisfies the undivided hypothesis and is complete, so it is consistent with both theorems.
5. **Multi-level exponents** (columns that are not ±1 patterns: marks taking pi^a conj(pi)^b with varying a),
   independent per-prime choices not organised as columns of one matrix, unions of designs, sparse hyperplane
   subfamilies, non-hyperplane transversals (README_VI §3(d) list).
6. m = 4 (clusters exist, §1.9).

## 3. (c) The two-core computation

3.1 **Criterion.** For the first Walsh core of order 8 (d = 3) inside the two-core profile, D = 14 * e/2 = 7e,
K = 7e/2 (e even). For each nonzero word y the product H_y of the four primes with lambda(y) = 1 satisfies
ord(H_y'/H_y) >= K + 1 (VI.1; d/dt of a product of ratios in U_K, char 2 makes (pi/conj pi)'/(pi/conj pi) = p'/p),
i.e. deg H_y' <= 4e - 7e/2 - 1 = e/2 - 1. Correct. The seven supports are the complements of the Fano lines; with
the empty set they form the [7,3,4] simplex code, and conversely every 3-dimensional subspace of 4-sets is such a
system (a constant-weight-4 linear code of dimension 3 without zero coordinates is the simplex code). So "no
3-dimensional subspace of admissible 4-sets" is exactly the impossibility of the first core. Correct.

3.2 **Key.** For a pair product P = A^2 + t B^2 (deg A = e, monic; deg B <= e - 1) and Q = C^2 + t E^2:
P' = B^2, (PQ)' = B^2 C^2 + A^2 E^2 = (AE + BC)^2. So deg (PQ)' <= e/2 - 1 iff deg(AE + BC) <= b =
floor((e/2-1)/2). Dividing by AC (degree 2e): iff ord_inf(B/A + E/C) >= 2e - b, iff
t^{2e-b-1}(B/A + E/C) has negative degree, iff the polynomial parts of t^{2e-b-1}B/A and t^{2e-b-1}E/C agree.
Correct (I derived it before reading her code). Each admissible 4-set is found from its three pairings.

3.3 **Reproduction.** Primes: own irreducibility test (t^{2^e} = t mod f and gcd conditions).

| e | primes | admissible 4-sets | lines (2-dim) | planes (3-dim) | method |
|---|---|---|---|---|---|
| 10 | 99 | 67 | 10 | 0 | direct over all 3,764,376 four-sets **and** key |
| 12 | 335 | 319 | 34 | 0 | direct over all 515,421,285 four-sets **and** key |
| 14 | 1161 | 4382 | 226 | 0 | key (each admissible set re-verified directly) |
| 16 | 4080 | 40546 | 1453 | 0 | key (re-verified directly, 128-bit products) |
| **18** | **14532** | **840521** | **44233** | **0** | key (re-verified directly) — **new, referee** |

All of Astra's numbers are reproduced exactly. Positive control for the plane search (e = 10, derivative bound
raised): bound 16 -> 1964 sets, 4848 lines, 116 planes; bound 20 -> 3859 planes; bound 24 -> 230604 planes. So the
zero at the true bound is not an artefact of the search.

3.4 **Does it settle anything beyond?** Only e <= 18 (my run) for d = 3, i.e. m = 8 with one prime of equal degree
per column and exponent one. It says nothing for e >= 20, for d >= 5, for unequal degrees, or for several primes per
column. Moreover, two structural remarks:
* The derivative conditions are **linear** in the word (char 2), and the row span of the two cores is a direct sum
  (VI.1). Hence the full VI.1 test of the x, x^3 profile is exactly two independent first-core tests on disjoint
  prime sets: the x, x^3 coupling plays **no role** in this test. Two disjoint admissible planes at some e would
  make the whole profile pass VI.1.
* **Heuristic (referee; random model).** One admissible condition costs c = 2e - 1 - b bits; a line needs 2c, a plane
  3c (linearity). Predicted admissible sets C(P,4)*2^{-c}: 29, 246, 4489, 42949, 864930 for e = 10..18 against observed 67,
  319, 4382, 40546, 840521. Predicted lines C(P,6)*15*2^{-2c}: 1.0, 6.4, 179, 1329, 42501 against observed
  10, 34, 226, 1453, 44233 — within 3-10% from e = 16 on. Predicted planes C(P,7)*30*2^{-3c}:
  e = 16: 0.006, e = 18: 0.08, e = 20: 0.16, e = 22: 2.6, e = 24: 5.8. The observed zeros up to e = 18 are what
  the random model expects; it predicts the first admissible first cores at **e about 22-24**, and then the
  two-core profile would pass the derivative test (expected count about the square). So the computation is not
  evidence of a structural obstruction; the case "e = 4d - 2" that VI.1 permits is just far below where solutions
  are expected. Passing VI.1 would still be only a necessary condition for an actual cluster.

## 4. (d) Growth over F_2[t]; what transfers to Z

* **Status unchanged:** 3M log_2 M - O(M) <= W(M) <= 2M^2 log_2 M (Part V; README_VI §3(d)). VII.4 removes one
  construction class whose surviving member (complete cores) already sits at M^2 log M; it neither improves the
  lower bound nor offers a construction.
* **Which steps exist over Z.** Identity 1.2 (pure algebra of Gaussian-integer block products): yes. Depth
  doubling (1.3) and the equality analysis (1.6, 1.8): **no** — |r_v| = 1 on the archimedean unit circle, which has
  torsion; r^{m/4} close to 1 only says r is close to an (m/4)-th root of unity. Gram step: the divisibility
  "common sub-product divides z_i - z_j" exists over Z[i], but the size bound becomes
  log N(P_ij) <= W/2 + 2 log C, so off-diagonals are <= 4 log C (log-weighted), **not <= 0** for C >= 1. A Laplacian
  argument with positive off-diagonal slack gives no rank bound without a spectral gap (vectors with pairwise
  inner products <= delta > 0 can be exponentially many). For C < 1 one gets strictly negative inner products and a
  trivial bound, already superseded by Cilleruelo-Cordoba type results for short arcs.
* **Record comparison (grep).** "Laplacian" occurs only in round20/kill_arakelov (an unrelated chain Laplacian).
  The Gram identities <f_x,f_y> = W' - 2d_xy and Plotkin balance are in the record (round8/referee_stability.md §3,
  Prop 4.1(b); round5/exact.md Lemma 2.2(e), with rank >= M - 1 under (AI)); overview §3: uniform theorem <=>
  bounded rank of the exponent lattice, and integer fakes have full rank M - 1. The Gram step's conclusion
  (rank >= M/2, pairing structure) is therefore weaker than what is already known/assumed for integer profiles, and
  its exact-sign hypothesis is unavailable over Z. **The Gram/Laplacian step gives nothing new over Z.**

## 5. Novelty and significance

* Novel in the record: the VII.4 theorem (and the Laplacian device in this context); the two-core computation; my
  e = 18 extension and the observation that the derivative test of the x, x^3 profile decouples into two single-core
  tests.
* Significance for the integer goal: **essentially none.** It closes one F_2 construction class (incomplete
  single-matrix Hadamard cores with all rows) using characteristic-2 depth doubling. Astra says so herself ("It
  supplies no corresponding integer bound"); I concur. Binary growth scale remains open; the most informative next
  step for (b) would be the e = 20-24 search for an admissible order-8 first core (predicted to exist), which would
  decide whether the x, x^3 profile can pass VI.1 at W = O(m log m) scale for m = 8 — still only a necessary
  condition.
