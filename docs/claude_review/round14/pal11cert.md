# Round 14 — pal11 product-certificate search (route 5, first overdetermined data point)

## 0. Bottom line

**No growth improvement and no uniform bound.** `M <= (2/3 + eps) log R/loglog R` stands. This is the route-5 item of
README-13 §7.1: the first overdetermined data point for the barrier test, on the fake pattern box pal11 (`T = 1`, `L = 2`,
`m = 8`, `n = 12`, `|A| = 12093`, 436 tangent octics in 11 variables).

| # | result | status |
|---|---|---|
| 1 | **The "`Z_11` symmetry" does not exist for pal11.** `Aut(A)` among signed coordinate permutations is `Z_5 x {±1}` (multipliers by squares). `Z_11` is the symmetry of Paley-12. (Lemma 1.1) | exact, exhaustive |
| 2 | **Certification set-up.** `rank_Q K_{<=8} = rank_p K_{<=8} = |A| = 12093` (two primes, `Z_5`-isotypic). So every mod-`p` non-membership is a rigorous rational non-membership (referee criterion). Validated against symmetry-free direct computations. (Prop 3.1) | exact |
| 3 | **The support-4 search is empty for all entries, with nothing to test.** Every balanced generator of support `<= 5` (any real entries) touches `>= 6` of the 11 Paley core columns. Each core column has width 4, so at most 7 such generators fit. Exact table of `min t` by support. (Theorem 4.1) | exact |
| 4 | **No product certificate exists at all.** No 8 real linear forms satisfy the width counts of pal11, so `J(A)_8^perp` (and even the fat-point system `L_8`) contains no product of linear forms, zonotopal or not. Proof: an exact profile reduction (25 profiles) plus 25 exhaustive exact searches. (Theorem 4.3) | exact, computer-assisted |
| 5 | **Non-zonotopal certificates:** none, by 4. `zeta(pal11) <= 7` already excluded zonotopes. | exact |
| 6 | **Round-13 correction.** fatpoint's `prodsearch.py` capped base-row incidence at 2, but `tau_A(1)` is 4-8. Corrected tables: `M8'` at `m = 10` **has** 1092 zonotope certificates (round 13: none). `M8` at `m = 8` has 42 count-feasible products, none a certificate, so that conclusion survives. Every certificate found in every exact model is a zonotope. (§4.5) | exact |
| 7 | **(2.2) at fixed `n` is nondegeneracy** (Prop 6.1). Whether `TC(pal11, 8)` has a nondegenerate real point is **not decided**. New exact structure: `TC` contains the 12 five-dimensional coordinate spaces `span(e_P)`, `P` a twin pattern (Prop 6.2). Evidence (200 minimisations): all real zeros found lie on rational subspaces or within `3.4e-4` of collision strata. (§6) | open; exact structure + evidence |

**Answer to the task.**
* *Support-4 and larger product-certificate search:* complete, for **all** generators. No product passes even the
  width counts, so no polarisation test is needed and none is reached.
* *Non-zonotopal certificate:* none exists.
* *Real tangent cone / (2.2):* reduced to the existence of a nondegenerate real point (Prop 6.1). That is not
  decided. It is not decidable by certificates (Theorem 4.3), and a saturation computation of the R11 J4 type (436
  octics in 11 variables) is out of reach. The numerical evidence (§6.3) is on the "degenerate" side.

**Fake-circle test.** Every statement is about the box and holds verbatim for fakes and for actual clusters with
this pattern set. Nothing here distinguishes them, and rigidity (Step 0) remains the only input fakes violate.

## 1. The model pal11 and the criteria used

**Labels.** *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted with an exact source.
*Exact*: integer or rational computation, or a mod-`p` computation used only in a direction in which it is rigorous
(stated each time). *Evidence*: anything else. Scripts are in `pal11cert_checks/` (`x.py` writes `x_out.txt` unless
said otherwise); working notes in `notes.md`. Sources: `round12/fake.md` (§1, §3, Thm 2.1, (2.2)),
`round13/fatpoint.md` (Props 2.5, 2.6, Thm 3.3) and `round13/referee_fatpoint.md` (§2.5, the rank criterion).

**1.1 The box.** `pal11` is the round-5 Paley-I fake profile at `q = 11` (`construct.md` §2), as built by
`fatpoint_checks/fplib.py::pal_patterns(11)` (copied, not imported): `n = 12` points, base row `0` (= `inf`), rows
`1 + i` for `i in F_11`, `d = 11`, `T = 1`, half-width `L = 2`:

```text
A = { k' in {-1,0,1}^11 : |k(P)| <= 2 for each of the 23 relative patterns P },   |A| = 12093,   m = 8.
```

The 23 patterns are the 11 Paley core columns (relative size 6) and 12 twins (relative size 5). All 23 slabs are
active (`k(P)` takes every value in `[-2,2]` on `A`), and the base coordinate `k_0 = -sum k'` ranges over `[-2,2]`
(`explore.py`). Its Hilbert function (fake §3.1, recomputed here, §3) is
`H_A = (1, 11, 66, 275, 880, 2253, 4664, 3507, 436)`, `rho = 8`, so `m = rho = 8` is the top order and the system is
**overdetermined**: `436` tangent octics in `11` variables. `zeta <= 7 < 8` (fake §3.1, LP certificate): `A` contains no
zonotope of length `8`.

**1.2 The symmetry group is Z_5, not Z_11.** The task (and README-13 §7.1) speaks of "the `Z_11` symmetry". That is
the symmetry of **Paley-12** (`fplib.paley12_patterns`: columns `i + D`, shift-invariant), not of `pal11`.

> **Lemma 1.1 (exact, `autgroup.py`).** The group of signed permutations of the 11 coordinates of `k'` that map
> `A` onto itself has order 10. It is `{±1} x Z_5`, where `Z_5` is the multiplier group `x -> a x` (`a` a nonzero square
> mod 11) acting on the rows `1 + x`. In particular no shift `x -> x + b` (`b != 0`) preserves `A`.

*Proof.* `autgroup.py` runs an exhaustive backtracking search over all `2^11 * 11!` signed permutations. It prunes
only with the invariant `#{k in A : k_x = s, k_y = t}`, which an automorphism must preserve, and it tests every
surviving complete map by exact set equality `g(A) = A`. Exactly 10 maps survive, listed in the output. `explore.py`
checks independently that among all 110 affine maps `x -> a x + b` only the five multipliers by squares preserve
the pattern set and `A`. QED

The shifts fail because the twin of the base row (`a* = 1`, `construct.md` §2) is not shift-invariant. `-1` acts on
a product of linear forms by a sign only. So the effective symmetry of the certificate search is `Z_5`, and it
is used below only to organise computations, never to restrict them.

**1.3 What a certificate is.** Notation as in fatpoint §1.2. `J(A)` is the ideal of leading forms of the vanishing
ideal of `A`. For `z_1..z_8 in Z^11 \ 0`, let `P(u) = prod_i <z_i,u>`.

* *Product certificate* (fatpoint Prop 2.5, refereed): `P in J(A)_8^perp = span{Phi_c : c in V(A,8)}`. Equivalently,
  the polarisations of `J(A)_8` vanish at `(z_1..z_8)`. If one exists, `TC(A,8) ⊂ ∪_i {<z_i,u> = 0}`, so `TC(A,8)` has
  no nondegenerate point.
* *Coordinates.* Write `P = sum_a P_a u^a`. Then `P in J(A)_8^perp` iff `P_a = 0` for every non-reduced `a` (some
  exponent `>= 3`) and the vector `(P_a / multinom(8;a))_a` over the reduced monomials of degree 8 lies in
  `T_8 = { (sum_k c_k k^a)_a : c in V(A,8) }`. Membership is solvability of the linear system
  `K_{<=8} c = (0; P_a/multinom(8;a))`, with `K_{<=8}` the `31824 x 12093` evaluation matrix of all reduced monomials of
  degree `<= 8` on `A`.
* *Width counts* (fatpoint Prop 2.6(c), refereed). A necessary condition is `N_v := #{i : <z_i,v> != 0} <= tau_A(v)` for
  every `v in Z^11 \ 0`, where `tau_A(v) = #{<v,k'> : k' in A} - 1`. It binds only when `tau_A(v) < 8`. For a family
  `F = (v_1, v_2)`, Prop 2.6(b) adds `#{i : z_i not orthogonal to span F} <= rho(B_F)`.
* *Rank criterion for negatives* (referee_fatpoint §2.5). If `rank_p K_{<=8} = rank_Q K_{<=8}`, then a system that is
  unsolvable mod `p` is unsolvable over `Q`. The right-hand sides are `p`-integral, since `multinom(8;a) | 8!`. Mod-`p`
  kernel containment alone is **not** a valid certificate. The referee's counterexample is `A = {0,p}`.
* *Zonotope.* `P` is zonotopal if `k_0 + sum_i [0,1] z_i ⊂ A` for some `k_0` (fake Lemma 5.1(1)). Since `zeta(pal11) <= 7`,
  **any** product certificate for `pal11` at `m = 8` would be non-zonotopal.

## 2. Search design: why enumeration fails, and the core-column count

**2.1 Naive enumeration is hopeless.** fatpoint's pool is the balanced `z in Z^12` (row 0 = base) with support `<= 4`
and entries `±1`: 66 pairs plus 1485 four-term vectors, 1551 generators. A product of 8 of them must pass the width
counts at the 1519 test directions with `tau_A < 8`. Those directions come from all of `{-1,0,1}^11` and from all of
`{-2..2}^11` with support `<= 4`, 129,713 directions in all; their widths are exact (`dirs.py`). A depth-first
enumeration in C (`csearch2.c`, bitsets, `Z_5` orbit representatives as first element) finds:
* 289,269 count-feasible multisets of size 2;
* 101,176,023 of size 3.

Size 8 is out of reach this way. The round-13 run, restricted to the 66 pairs, already took 23 minutes in Python.
The tree has to be cut by an argument.

**2.2 The width counts that matter.** The 1519 directions with `tau_A < 8` have `tau in {2, 4, 6}`:
* 11 coordinate directions `e_x` (`tau = 2`);
* 134 directions with `tau = 4`: the 11 Paley **core** columns `1_P` (`|P| = 6`), the 12 **twin** columns (`|P| = 5`),
  the all-ones direction (the base coordinate, since `<1,k'> = -k_0`), and the 110 directions `e_x ± e_y`;
* 1374 directions with `tau = 6`.

For a balanced `z in Z^12` let `sigma(z) = |supp z|` (base row included) and
`t(z) = #{core columns P : z(P) != 0} = #{c : <z, beta_c> != 0}`, where `beta_c = 2*1_P - 1`. The 11 vectors `beta_c`,
together with the all-ones vector, are the columns of the `12 x 12` Paley Hadamard matrix (asserted in
`hadflats2.py`). The width counts give:

```text
(core)       sum_i t(z_i) <= 11 * 4 = 44 ,        (incidence)   sum_i sigma(z_i) <= 11 * 2 + 4 = 26 .
```

The incidence bound uses `tau(e_x) = 2` for the 11 non-base rows and `tau(1) = 4` for the base row.

**2.3 The LP shadow.** The LP relaxation of the width counts is `max sum x_i` subject to `sum_i x_i [z_i touches v] <= tau(v)`
and `x >= 0` (`lpcount.py`, `lpdual.py`, scipy/HiGHS; floating point, used only for orientation). Its value is:
* exactly `22/3 = 7.333` for the pools {support `<= 4`, entries `±1`} (1551), {support `<= 5`, entries `±1`} (1551; no
  balanced `±1` vector has odd support), {support `<= 3`, entries `<= 2`} (726) and {support `<= 4`, entries `<= 2`} (8151);
* `8.259` for {support `<= 6`, entries `±1`} (10,791).

The optimal dual of the small pools is `1/6` on each of the 11 core columns and `0` elsewhere. It certifies the
value `44/6 = 22/3` as soon as every generator touches at least 6 core columns. That is the content of Theorem 4.1
below, and it is proved there **for all real generators**, not only for the pools. Adding all 28,444 informative
2-direction family counts of fatpoint Prop 2.6(b), `#{i : z_i touches v or v'} <= rho(B_{v,v'})` with `rho(B)` exact
(`fam2.py`), leaves the LP values unchanged, `7.333` and `8.259` (`lpfam.py`, `lpfam2.py`).

## 3. Exact certification method (rank_p = rank_Q)

**Proposition 3.1 (exact).** `rank_Q K_{<=8} = 12093 = |A|`. Hence, for every product `P` with `p`-integral reduced
coefficients, mod-`p` unsolvability of `K_{<=8} c = (0; P_a/multinom(8;a))` for `p in {p_1, p_2}` implies
`P notin J(A)_8^perp`. Here `p_1 = 1152921504606847081` and `p_2 = 576460752303423761`, both `= 1 mod 5`.

*Proof.* `isotan.py` computes `rank_p K_{<=8}` for `p = p_1, p_2` in the `Z_5`-isotypic basis. The change of basis is
invertible over `F_p` because `p = 1 mod 5`, so `F_p` contains the 5th roots of unity and `5` is a unit. The five
blocks have ranks `2421, 2418, 2418, 2418, 2418`, summing to `12093` (output `isotan_p1_out.txt`, `isotan_p2_out.txt`).
Since `rank_p M <= rank_Q M` for integer matrices and `rank_Q K_{<=8} <= #columns = |A|`, we get
`rank_Q K_{<=8} = rank_p K_{<=8} = |A|`. The second statement is the referee_fatpoint §2.5 criterion. If the system
is solvable over `Q`, then `rank_Q [K|b] = rank_Q K = |A|`, so `rank_p [K|b] <= |A| = rank_p K`, and the system is
solvable mod `p`. The right-hand side is `p`-integral because `multinom(8;a) | 8!` and `p > 8!`. QED

The same runs reproduce, as mod-`p` lower bounds for the ranks, the cumulative values `1, 12, 78, 353, 1233, 3486,
8150, 11657, 12093` of fake §3.1 (`H_A = 1, 11, 66, 275, 880, 2253, 4664, 3507, 436`). They also give
`dim_{F_p} T_8 = 436` (blocks `88, 87, 87, 87, 87`). Over `Q` only `dim T_8 = dim V(A,8) = 12093 - rank_Q K_{<8} <= 436` is
proved this way. No negative below uses the value 436: every negative uses only Proposition 3.1.

**Bug found and fixed before any use.** A first version of `isotan.py` scaled the rows of the `Z_5`-fixed monomials
inconsistently: free columns got `E(a,k)` where `sum_t E(a, g^t k) = 5 E(a,k)` is correct. Only the trivial block was
affected. The fix is marked in the source, and every number in this write-up comes from the fixed code.

**Validation** (`validate.py`, prime `p_2`, output `validate_out.txt`), all passed:
1. Each of the 436 isotypic kernel vectors is lifted back to a function `c` on `A`. For 15 of them (3 per block),
   `E_{<8} c = 0` is checked **directly**, on all 17,469 reduced monomials of degree `< 8`, without using the symmetry.
2. For the same `c`, the moment vector `mu_8(c)`, computed directly, is declared a member by the isotypic membership
   test (`memlib.py`).
3. Every perturbed vector (`mu_8(c)` plus one unit in a random coordinate) and a random vector are declared
   non-members.
4. The 436 lifted vectors are linearly independent mod `p_2`.

So the isotypic coordinates, the projection formulae and the membership test agree with the untransformed
definitions.

## 4. Results of the certificate search

### 4.1 The support-4 search is complete, and empty, for every entry size

**Theorem 4.1 (exact).** Let `z_1..z_8` be nonzero balanced real vectors on the 12 rows (equivalently, nonzero
`z'_i in R^11`, `z_{i,0} = -sum z'_i`), let `sigma_i = |supp z_i|` with the base row counted, and let `t_i` be the
number of core columns touched by `z_i`. Then:

(a) `t(z) >= f(sigma(z))` for every nonzero balanced `z`, with

| `sigma` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `f(sigma) = min t` | 6 | 9 | 6 | 7 | 2 | 5 | 4 | 3 | 4 | 4 | 1 |

(b) If `prod_i <z'_i, u> in J(A)_8^perp` (a product certificate for pal11 at `m = 8`), then `sum_i t_i <= 44` and
`sum_i sigma_i <= 26`. In particular **some `z_i` has `sigma_i >= 6` and `t_i <= 5`**.

(c) Consequently no product certificate has all generators of support `<= 5`, whatever the entries. This covers
fatpoint's support-4 pool, all of its entry-bounded extensions, and every pool of support `<= 5`. There is nothing to
test for membership.

*Proof.* (a) Fix a support `R`. A balanced `z` with `supp z = R` is orthogonal to the all-ones column. Let `Z(z)` be
the set of Hadamard columns `h` (all-ones included) with `<z,h> = 0`. Restricted to `R`, `Z(z)` is a flat of the
vector configuration `{h|_R}`, and `z` is a generic point of the orthogonal complement of `span{h|_R : h in Z(z)}`.
Conversely, for every flat `Z` whose complement contains a vector of full support `R`, a generic such vector `z` has
`Z(z) = Z` and `t(z) = 12 - |Z|`. So `min t` over support `R` is `12 - max|Z|` over the flats with a full-support
complement. `hadflats2.py` enumerates, for each of the `2^12 - 13` supports `R` with `|R| >= 2`, all flats containing the
all-ones column (closure by exact integer rank, `fmpz_mat`). It tests full support of the complement exactly, using
an exact nullspace basis (a generic element has full support iff no coordinate vanishes on the whole basis). The
table is the output (`hadflats2_out.txt`). The Hadamard property `H H^T = 12 I` is asserted in the script.
(b) Each core column `1_P` is a test direction with `tau_A = 4` (`k(P)` takes the 5 values `-2..2` on `A`). By
fatpoint Prop 2.6(c) (refereed; necessity of the width counts for every product certificate), at most 4 of the
`z_i` touch it, so `sum_i t_i <= 44`. In the same way `tau_A(e_x) = 2` and `tau_A(1) = 4` (the base row) give
`sum_i sigma_i <= 26`. If every `t_i >= 6`, then `sum t_i >= 48 > 44`. By (a), `t_i <= 5` forces `sigma_i >= 6`.
(c) is immediate from (a) and (b). QED

**Remark.** Theorem 4.1(b) sharpens fatpoint Thm 3.3 for this box. The uncertainty principle alone gives `t >= 12/sigma`
and allowed `m <= 9`. The exact table shows that small supports cost **at least 6 core touches**: the Paley structure
gives more than the uncertainty bound at `sigma = 4, 5`. The LP value `22/3` of §2.3 is exactly `44/6`.

### 4.2 Profiles: what a certificate would have to look like

Combining (a), (b) and the exact list of **attainable** pairs `(sigma, t)` with `t <= 5` (`excflats.py`; see below),
`profiles.py` enumerates every multiset of eight `(sigma_i, t_i)` with `sum sigma <= 26`, `sum t <= 44`, `t >= f(sigma)`.
The exceptional generators (`t <= 5`) that occur are:

| `(sigma, t)` | how they arise (exact; `excflats_out.txt`) |
|---|---|
| `(6,2)` | `(beta_a ± beta_b)/2`: 110 vectors |
| `(6,4)` | 1320 explicit vectors (complement dimension 1) |
| `(7,5)` | 9240 explicit vectors |
| `(8,4)` | 1980 planes (complement dimension 2) |
| `(8,5)` | 2310 explicit vectors and 13,860 planes |
| `(12,1)` | `beta_c`: 11 vectors |
| `(12,2)` | `a beta_c + b beta_d`, `a != ±b`: 55 planes |

No support-6 generator has `t = 3` or `t = 5`. The possible exceptional parts of a certificate are therefore:
`{(6,2)}`, `{(6,2),(6,2)}`, `{(6,2),(6,4)}`, `{(6,2),(7,5)}`, `{(6,2),(8,4)}`, `{(6,2),(8,5)}`, `{(6,4),(6,4)}`, `{(6,4),(8,4)}`,
`{(12,1)}` and `{(12,2)}`. In the last two cases the other seven generators are pairs. In most others the core
budget is exactly tight, so every non-exceptional generator touches exactly 6 core columns.

**Lemma 4.2 (exact attainable types, `typeflats.py`, `excflats.py`).** For `sigma <= 8`, the attainable values of `t` are:

| `sigma` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| attainable `t` | 6 | 9 | 6, 8, 9, 10, 11 | 7..11 | 2, 4, 6..11 | 5..11 | 4..11 |

For `sigma = 9..12` only the minimum `f(sigma)` is used. For each type the flats are listed together with the dimension of
their complement. The types `(2,6)`, `(4,6)`, `(4,9)`, `(6,2)`, `(6,4)`, `(7,5)`, `(12,1)` have **only** flats with
1-dimensional complement, so the explicit lists (66, 990, 495, 110, 1320, 9240, 11 vectors, entries `|z_x| <= 2`) contain
**every** generator of that type, up to scale. *Proof:* as for Theorem 4.1(a), with all flats recorded instead of the
maximal one. A random test is in `sanity_types.py`: 295,447 random balanced integer vectors, no violation. QED

With these exact type sets, `profiles2.py` finds exactly **25** feasible profiles (listed in `univrun_out.txt`). Every
one of them has at least 4 pairs and at most one generator of a *parametric* type, i.e. one of `(3,9)`, `(4,8)`, `(4,10)`, `(6,6)`,
`(7,6)`, `(8,4)`, `(8,5)`, `(8,6)`, `(12,2)`.

### 4.3 Exhaustive count search on explicit pools

`csearch4.c` enumerates every multiset of 8 pool generators satisfying all 1519 width counts. Exactness:
* integer arithmetic only;
* pruning only by consequences of the counts: saturation bitsets, and knapsack bounds `sum_i w_q(z_i) <= sum_v Y_q(v)(tau(v) - N_v)` for 10 nonnegative integer weight vectors `Y_q` (`mkweights.py`), using that each
  generator occurs at most twice (it touches a coordinate direction with `tau = 2`);
* the first element is restricted to a `Z_5` orbit representative, which is complete by Lemma 1.1;
* the pool is ordered with exceptional generators first, so the bound becomes tight once they are placed.

| pool | generators | count-feasible 8-multisets | certificates |
|---|---|---|---|
| support `<= 4`, entries `±1` (fatpoint's pool) | 1551 | **0** (also by Theorem 4.1) | 0 |
| support `<= 4`, entries `<= 2` in absolute value | 8151 | **0** (Theorem 4.1) | 0 |
| support `<= 6`, entries `±1` | 10,791 | **0** (`csearch4`, 3 parts; nodes per level `1, 22, 2836, 91815, 1.0e6, 1.8e6, 1.8e5, 0, 0`) | 0 |
| all balanced `±1` vectors (any support) | 36,894 | run stopped once Theorem 4.3 (all real generators) was complete; it contains this pool | 0 |

### 4.4 The universal result: no product of 8 linear forms satisfies the width counts

**Theorem 4.3 (exact; computer-assisted).** For pal11 (`T = 1`, `L = 2`, `|A| = 12093`), there are no nonzero real
vectors `z'_1..z'_8 in R^11` such that `#{i : <z'_i, v> != 0} <= tau_A(v)` for all 1519 test directions `v` with
`tau_A(v) < 8`. Consequently:
* (a) `J(A)_8^perp` contains **no** product of 8 real linear forms. There is no product certificate of any kind, with
  any supports, any entries, rational or irrational.
* (b) The same holds for the fat-point (width) system `L_8 = V_{F_1}(8)` of refJ §2.2, and for every local system
  `V_F(8)` whose family set contains these 1519 directions (fatpoint Prop 2.6(c)).

*Proof.* Suppose `z_1..z_8` satisfy the counts. Their profile `{(sigma_i, t_i)}` is then one of the 25 profiles of
§4.2, by Theorem 4.1(b) and Lemma 4.2. Each profile `p` is excluded by one exhaustive run of `csearch7` (profiles 0-17) or `csearch9` (profiles 18-24)
(driver `univrun.py`; outputs `univrun_out_0_17.txt`, `univrun_out2.txt`, `univrun_summary.txt`):

1. *All types explicit* (13 profiles). By Lemma 4.2 the generators lie in the explicit lists. The run enumerates every
   multiset with exactly `p`'s composition, using `Z_5` orbit representatives for the first element (Lemma 1.1;
   compositions are `Z_5`-invariant). It finds 0 count-feasible multisets.
2. *`p` = 7 pairs + one `(12,2)` generator `z = a beta_c + b beta_d`*. Write `Pi = span(beta_c, beta_d)`. A direction
   `v` restricts to a linear functional on `Pi`. Every nonzero restricted functional vanishes on exactly one line of
   `Pi`. Take one integer point on each such line: these are the 1001 representatives (`reps_12_2.pkl`).
   * If `z` lies on one of these lines, `z` is a representative up to scale.
   * Otherwise `z` touches every direction whose functional on `Pi` is nonzero. Any representative touches a subset
     of these, so replacing `z` by it preserves all counts.

   The run enumerates 7 pairs together with each representative and finds 0 count-feasible multisets.
3. *Exactly one generator of a parametric type `H`* (11 profiles; `H` is one of `(3,9)`, `(4,8)`, `(4,10)`, `(6,6)`,
   `(7,6)`, `(8,4)`, `(8,5)`, `(8,6)`). The other 7 generators have explicit types. The run enumerates them exactly
   (FREE mode). For each count-feasible 7-multiset it computes the saturated set `S = {v : N_v = tau_A(v)}`. The 8th
   generator must be orthogonal to every `v in S`, since otherwise that count is exceeded. So `rank S <= 10` is
   necessary. The rank is computed mod `P = 2147483629`, and `rank_P(S) = 11` implies `rank_Q(S) = 11`. Every
   count-feasible 7-multiset found has `rank_P(S) = 11`, so none can be completed.

Every run is exact: integer counts and pruning only by consequences of the counts (§4.3). Its bounds are valid for the
free slot too: the free generator's weights are either not subtracted, or (`csearch8/9`) subtracted only by their exact
values `t_f` and `sigma_f` for the profile's free type, and `0` for the twin and `tau = 6` vectors. QED

**Run data** (`univrun_summary.txt`). The 25 profiles, numbered 0-24 in the order printed by `profiles2.py`:
* *Full mode* (`p` = 0, 1, 3, 4, 5, 7, 12, 15, 16, 17, 19, 22, 24; and `p = 2` with the `(12,2)` representatives): **0
  count-feasible 8-multisets** in every run. The longest run was `p = 19`, at 1917 s on 4 cores.
* *Free-slot mode*: all count-feasible 7-multisets have saturated rank 11, so 0 are completable:

  | `p` | free type | rest | count-feasible 7-multisets |
  |---|---|---|---|
  | 6, 8, 9, 10, 11 | `(6,6)`, `(7,6)`, `(8,4)`, `(8,5)`, `(8,6)` | 6 pairs + `(6,2)` | 14,442 each |
  | 13 | `(8,4)` | 6 pairs + `(6,4)` | 11,799 |
  | 14 | `(3,9)` | 5 pairs + 2 `(6,2)` | 7230 |
  | 18 | `(6,6)` | 5 pairs + `(4,6)` + `(6,2)` | 5042 |
  | 20 | `(4,8)` | 5 pairs + 2 `(6,2)` | 2779 |
  | 21 | `(4,8)` | 5 pairs + `(6,2)` + `(6,4)` | 246 |
  | 23 | `(4,10)` | 5 pairs + 2 `(6,2)` | 0 |

The counts differ between rows with the same rest because of the reserved weights of `csearch8/9` (profiles 18-24).
Profiles 0-17 ran with `csearch7`, profiles 18-24 with `csearch9`. `csearch9` adds an exact column-requirement prune:
a 0/1 weight class with slack `s` forces at least `R_j - s - FREE` further touches of each of its directions `j`.
Cross-check: on profile 3, `csearch9` reproduces the `csearch7` result (0 leaves) with about 8 times less work
(level-6 nodes: 311 against 176,000).

### 4.5 Correction to the round-13 certificate tables (fatpoint §5.2)

**Finding (exact).** `fatpoint_checks/prodsearch.py` prunes with `any(inc[x] >= 2 for x in S)` over **all** rows,
the base row included, and with `capacity = 2n`. In effect it imposes **base incidence `<= 2`**. Its docstring
justifies this as "`v = -1` (base incidence `<= 2`)". But the width of the base direction is
`tau_A(1) = #{k_0 = -sum k'} - 1`, and that is `4` (sdc7, `L = 2`), `6` (pal7, `L = 2`), `6` (`M8`), `8` (`M8'`) and `4`
(pal11), not 2. So the round-13 searches were **not** exhaustive even within their pool.

The positive control `control_pal7.py` reproduces fatpoint's 57 exactly when base incidence is capped at 2. Without
the cap there are 2134. Rerunning the same pool (support `<= 4`, entries `±1`, test directions `{-2..2}^7`) with the
correct counts (`control_model.py`, `csearch4`), followed by membership tests (`memsmall.py`), gives:

| model | `m` | round 13: pass counts / certificates | corrected: pass counts | certificates (all zonotopes) | `rank_Q K_{<=m} = rank_p` |
|---|---|---|---|---|---|
| sdc7 | 8 | 0 / 0 | 0 | 0 | — |
| sdc7 | 7 | 24 / 24 | 45 | 24 | 756 = 756 |
| pal7 | 7 | 57 / 57 | 2134 | **300** | 861 = 861 |
| `M8` | 8 | 0 / 0 | **42** | **0** (rigorous) | 1037 = 1037 |
| `M8` | 7 | 696 / 696 | 23,355 | **1746** | 995 = 995 |
| `M8'` | 10 | 0 / 0 | **1092** | **1092** | 1870 = 1870 |
| `M8'` | 11 | 0 / 0 | 0 | 0 | — |

How the columns are certified:
* `rank_Q` is an exact `fmpz_mat` rank. Since `rank_Q = rank_p`, every non-membership is rigorous.
* Every member found mod `p` passes the exact vertex check: a translate of the zonotope lies in `A`. So by fake
  Lemma 5.1(1) every member is an **exact** certificate. The first five of each model were also re-checked by exact
  `fmpq` rank.

**Consequences for the record.**
* README-13 §1.4 item 7 says "non-existence for `M8` (`m = 8`) and `M8'` (`m = 10`) holds only within that pool".
  * For `M8` at `m = 8` the conclusion survives (now rigorous within the pool): 42 products pass the counts and none
    is a certificate. "Certificates are not necessary" (`M8` at `m = 8`) therefore still stands. Its justification
    changes from "0 pass the counts" to "42 pass, 0 are members".
  * For `M8'` at `m = 10` it is **false**. The pool contains 1092 zonotope certificates of length 10. One is
    `-e_1, -e_1, -e_2, -e_2, -e_3, -e_4, -e_5, -e_6, e_3 - e_4, e_5 - e_6` (in `k'` coordinates). It uses the base row 8 times
    (`= tau_A(1)`). So `zeta(M8') = 10`, consistent with the earlier bound `zeta <= 10` (fake §3.1), and
    `TC(M8', 10)` has no nondegenerate point. That is consistent with J4 at `m = 11`. The fatpoint reading "the search
    adds `M8'` at `m = 10`" to the zonotope-free cases must be withdrawn.
* The qualitative pattern of round 13 survives and is strengthened. In the exact `n = 8` models, **every**
  certificate found (24 + 300 + 1746 + 1092) is a zonotope.
* Theorem 4.3 does not depend on this bug: its counts use the true `tau_A(1) = 4` of pal11.

## 5. Non-zonotopal certificates?

The task asked for any certificate that is **not** a zonotope. For pal11 at `m = 8` there is nothing to find, for
two independent reasons:
* `zeta(pal11) <= 7 < 8` (fake §3.1). `A` contains no zonotope of length 8, so **every** product certificate
  would be non-zonotopal. This is why pal11 is the right first overdetermined test.
* Theorem 4.3 (§4.4) shows that **no** product of 8 real linear forms satisfies even the width counts. So there is no
  product certificate of any kind: not zonotopal, not non-zonotopal, not with irrational coefficients.

The membership machinery of §3 (validated, rigorous for negatives) was therefore never needed on a search leaf: no
leaf reached the membership test. In the corrected round-13 data (§4.5), the width counts alone kill every pool
product only for `sdc7` at `m = 8` and `M8'` at `m = 11`. For `M8` at `m = 8`, 42 pool products pass the counts and
fail membership. pal11 is the first **overdetermined, faithful, fake** box in which the counts kill every product,
and the statement covers all real generators, not a pool.

What this does **not** say: it does not say that `TC(A,8)` has nondegenerate points. `M8` (round 11 J4) has no
nondegenerate real tangent-cone point and also no product certificate. Product certificates are sufficient, not
necessary (fatpoint §5.2). §6 deals with the tangent cone directly.

## 6. Real tangent cone and the Diophantine condition (2.2)

### 6.1 At fixed n, (2.2) is nondegeneracy

**Proposition 6.1.** Use the pal11 pattern data of fake §1: `n = 12`, `r = 9n - 8 = 100`, `rho_1 = 1`, `k = 0`. Then
every nondegenerate `u* in R^11` satisfies the Diophantine condition (2.2) of fake.md with some finite `a = a(u*)`.
The hypotheses of fake Thm 2.1 together with (2.2), applied to this one box, therefore reduce to "`u*` is a
nondegenerate point of `TC(A,8)`".

*Proof.* The condition quantifies over the relevant `0 != lambda in L_T` with `||lambda||_1 <= r/rho_1 + 1 = 101`. The
lattice `L_T` is discrete, so this set is finite. A nonzero zero-sum `lambda` has a nonzero rational restriction
`lambda'` to the 11 non-base rows. With `u*_0 = 0`, `lambda.u* = <lambda', u*>`, which is nonzero because `u*` is
nondegenerate. Take `a >= max_lambda ln(1/|lambda.u*|)/(||lambda||_1 ln 12)` over the finite set. QED

This was already observed for sdc7 in fatpoint §5.2. So the question for pal11 is exactly: **does `TC(A,8)` contain a
nondegenerate real point?** The quantitative content of (2.2), a bounded `a` as `n -> infinity`, cannot be tested
on a single box.

### 6.2 What is proved about TC(pal11, 8)

1. *Lower bound (refJ Prop 2.2, fatpoint Prop 2.3; quoted).* `TC(A,8)` contains every cheap span `span(v_1..v_s)` with
   `sum tau_A(v_i) < 8`. For pal11, from §2.2 (`tau in {2,4,6}` on the test set), these include:
   * the 165 coordinate 3-spaces `span(e_x, e_y, e_z)`;
   * the planes `span(e_x, v)` with `tau(v) = 4` (core, twin, all-ones, `e_y ± e_w`);
   * the lines `R v` with `tau(v) = 6`;
   * and joins of family pieces with total `nu < 8`.

   All of these are rational subspaces, hence degenerate.
2. *No product certificate* (Theorem 4.3). This also excludes every certificate of the form
   `prod_{i<=8} <z_i,u>` in the local systems `V_F(8)` that contain the width conditions.
3. *Span principle (fatpoint Prop 2.4; quoted).* A nondegenerate real point of `TC(A,8)` needs a `Q`-irreducible
   component of `TC_C` whose linear span is `C^11`. If the point is algebraic, it has degree `>= 11` over `Q`.


**Proposition 6.2 (exact; `tcpost.py`, output `tcpost_out.txt`).** For a set `S` of coordinates of `k'`, let
`B_S = pi_S(A)`. If `rank_p K_{<=7}(B_S) = |B_S|`, then `span(e_x : x in S) ⊂ TC(A,8)`. This holds for:
* all `S` with `|S| <= 3` (165 three-sets);
* 115 four-sets;
* exactly 12 five-sets, namely the **12 twin patterns** `P` (relative size 5).

No set with `|S| = 6, 7, 8` passes the test. The maximal certified coordinate subspaces are the 12 five-dimensional
`span(e_P)` (`P` a twin) and 55 four-dimensional ones.

*Proof.* `rank_p <= rank_Q <= |B_S|`, so `rank_Q K_{<=7}(B_S) = |B_S|`. Hence `H_{B_S}(e) = 0` for `e >= 8`,
`J(B_S)_8` contains all forms of degree 8, and `TC(B_S, 8)` is everything. By fatpoint Prop 2.2(b), the base locus of the
coordinate family `F_S` is all of `span(e_S)`, and fatpoint Lemma 2.1 puts it inside `TC(A,8)`. QED

The tangent cone therefore has rational linear components of dimension at least 5. These are not cheap spans in the
width sense, since `sum tau(e_x) = 10 > 8`. They are projection-cheap: the twin slab `|k(P)| <= 2` cuts the projection.
All of them are degenerate.

Beyond items 1-3 and Proposition 6.2, nothing is proved. In particular it is **not known** whether `TC_R(A,8)` is a
finite union of rational subspaces, which would mean no nondegenerate points. Its negation, a nondegenerate point,
is not known either.

### 6.3 Numerical evidence (EVIDENCE ONLY)

`tcnum.py` computes a floating orthonormal basis of `V(A,8)` from the isotypic blocks. Check:
`max |E_{<8} Q| = 3.1e-15` on 400 random monomials of degree `< 8`, and the rank is 436. It minimises
`F(u) = ||Q^T y(u)||^2/||y(u)||^2`, where `y(u)_k = <k',u>^8`, over the unit sphere from 200 random starts (L-BFGS,
seed 1); `F = 0` iff `u in TC(A,8)`. `tcpost.py` and `tcpost2.py` classify the minimisers (outputs `tcnum_out.txt`,
`tcpost_out.txt`, `tcpost2_out.txt`):
* **All 200 runs reach `F < 2e-18`.** The cone is large: it contains 5-dimensional linear pieces.
* **93 of the 200** are zeros (`F <= 1e-23`) whose support lies inside one of the certified coordinate subspaces of
  Proposition 6.2: 79 with support 5 (twin patterns) and 14 with support 4. At 87 of these 93 points the numerical
  Jacobian of the 436 forms has rank exactly the codimension of the containing coordinate space: 6 at 73 of the
  support-5 points, 7 at all 14 support-4 points (singular values `~1e-2` against `~1e-14`). So near them the cone is
  that coordinate space, with no extra branch. The other 6 support-5 points report rank 7 at the `1e-6` threshold.
  That is a numerical artefact: containment in a 5-dimensional linear piece of the cone forces rank `<= 6`.
* **The other 107** all lie within `3.4e-4` of a collision stratum. At tolerance `1e-3`, their 12 positions
  `(0, u_1..u_11)` form only 3-8 clusters (counts 9, 68, 21, 6, 2, 1 for 3, 4, 5, 6, 7, 8 clusters), and `F` ranges
  over `1e-29..2e-18`. No minimiser has its 12 positions pairwise separated by `1e-3`. The best minimal gap found is
  `3.4e-4`, against `2/11 = 0.18` for equally spaced points.

Reading. This agrees with, and sharpens, fake §4.2 (pal11 at `s >= 0.02`: `kappa2* > 0`). Every numerically found
point of the real cone is on, or extremely close to, a rational subspace: a twin coordinate space, a coordinate
4-space, or a collision stratum. The cone shows no sign of a nondegenerate component. **This is evidence, not
proof.** Local minimisation can miss small components, and near-strata points with `F ~ 1e-20` are not
distinguished from zeros by this method.

## 7. What this does and does not show

**Growth / uniformity.** Nothing here improves `M <= (2/3 + eps) log R/loglog R`. pal11 is a toy box (`n = 12`,
LP-subcritical, `m/(2L) = 2`), and every statement is about the box alone. It holds verbatim for fakes and for actual
clusters with this pattern set.

**For route 5's barrier test (README-13 §3, §7.1).**
1. *The first overdetermined data point is a complete negative for product certificates* (Theorem 4.3). It covers
   all generators, not a pool. In particular the "non-zonotopal certificate" that would have been the first sign that
   route 5 is cappable does **not** exist in pal11. README-13 §7.1 anticipated that a complete negative would be
   weak evidence the other way. It is weaker still, because the negative is forced already at the level of the width
   (`b = 1`) counts. The obstruction is the **core budget** `sum t_i <= 44` together with the Paley uncertainty table
   (Theorem 4.1). Pairs cost 6 core touches each, so only "Hadamard-aligned" generators (`(beta_a ± beta_b)/2`,
   `beta_c`) are cheap, and those exhaust the incidence budget.
2. *The mechanism is the regime one.* Theorem 4.1(b) is the exact small-`n` analogue of fatpoint Thm 3.1 / Prop 3.2.
   Pairs pay `n/2` core touches, so certificates need heavy generators, and the heavy ones are too wide for the
   incidence budget. In pal11 this closes completely. In the regime, fatpoint Thm 3.3 shows that the corresponding
   budget is **not** binding (factor `1.15` at `n = 2^12`). So the pal11 negative does not extrapolate. It confirms
   only that in toy boxes the certificate route fails for the same counting reason in fake and genuine profiles.
3. *The tangent-cone question for pal11 stays open.* No product certificate exists, so the barrier question cannot
   be settled for pal11 by the only exact upper-bound device the record has. §6 gives numerical evidence that the
   real tangent cone consists of degenerate points. That is the behaviour ZE-TC predicts, and the opposite of what a
   cap of route 5 would need. But it is not a proof. M8 shows that degenerate tangent cones can exist without
   certificates; there the proof went through a saturation argument (R11 J4) that is out of reach here: 436 octics in
   11 variables, against 42 in 7.
4. *Methodological correction.* The "`Z_11` symmetry" of the task belongs to Paley-12, not to pal11, whose symmetry
   group is `Z_5 x {±1}` (Lemma 1.1). Mod-`p` negatives are certified by `rank_p K_{<=8} = |A|` (Prop 3.1), as the
   referee required. In the end no membership test was needed, because no product passed the counts.

**Recommendation.** README-13 §7.1 was right to call this the one cheap route-5 item, and it is now done. It gives no
reason to reopen route 5. The certificate device is dead in pal11 at the count level, and the tangent-cone question
is not computable by the tools at hand. Next work should go to README-13 §7.2 (Rig(K), P3-all\* in regime (1.4)).

## 8. File index (`pal11cert_checks/`)

| script | output | what | status |
|---|---|---|---|
| `fplib.py` | — | copy of `round13/fatpoint_checks/fplib.py` (models, boxes, monomials) | library |
| `explore.py` | (stdout) | pal11 box, patterns, affine symmetries | exact |
| `autgroup.py` | (stdout) | full signed-permutation automorphism group of `A` (order 10) | exact (Lemma 1.1) |
| `dirs.py` | `dirs_pal11.npz`, `dirs_pal11.txt` | widths `tau_A(v)` of 129,713 test directions; the 1519 with `tau < 8` | exact |
| `isotan.py` | `isotan_p1_out.txt`, `isotan_p2_out.txt`, `T8_p*.pkl` | `Z_5`-isotypic ranks of `K_{<=e}`; `rank_p K_{<=8} = 12093`; `T_8` mod `p` | exact (Prop 3.1) |
| `memlib.py` | — | isotypic membership test in `T_8` | library |
| `validate.py` | `validate_out.txt` | direct (symmetry-free) validation of kernel, membership, non-membership | exact checks |
| `csearch2.c` | (stdout) | naive bitset enumeration (tree sizes, §2.1) | exact counts |
| `lpcount.py`, `lpdual.py`, `lp6.py`, `lpfam.py`, `lpfam2.py` | (stdout) | LP relaxations, duals, 2-family constraints | floating point, orientation only |
| `fam2.py` | `fam2_out.txt`, `fam2_pal11.npy` | `rho(B_F)` for 2-direction families | exact |
| `hadflats2.py` | `hadflats2_out.txt` | `min t` by support (Theorem 4.1(a)) | exact |
| `typeflats.py`, `excflats.py` | `typeflats_out.txt`, `excflats_out.txt`, `*.pkl` | all flats/types for `sigma <= 8` and `sigma = 12` (Lemma 4.2) | exact |
| `sanity_types.py` | (stdout) | random cross-check of the type tables | evidence (consistency) |
| `profiles.py`, `profiles2.py`, `freeslot.py` | (stdout) | feasible `(sigma,t)` profiles (25) | exact enumeration |
| `reps_by_type.py`, `univpool.py` | `reps_12_2.pkl`, `univpool_out.txt` | representatives of parametric types | exact |
| `mkpool.py`, `mkpool2.py`, `mkpool7.py`, `mkweights.py` | `pool_*.txt`, `weights.txt` | pools, knapsack weight vectors | — |
| `csearch3.c` ... `csearch8.c` (all exact; 3 = knapsack bounds, 4 = stronger candidate filter, 5/6 = free slot + composition, 7/8 = per-profile + reserves) | `out4_b61_*.txt`, `univrun_*` | exact count searches (§4.3, §4.4) | exact |
| `univrun.py` | `univrun_out_0_17.txt`, `univrun_out2.txt`, `univrun_summary.txt` | driver for Theorem 4.3 | exact |
| `tcnum.py` | `tcnum_out.txt` | real zeros of the 436 octics, Jacobian ranks (its LLL 'rational-span dim' column is unreliable; use `tcpost2.py`) | **evidence** |
| `tcpost.py` | `tcpost_out.txt`, `tc_coord_maximal.pkl` | coordinate subspaces inside `TC` (Prop 6.2) | exact (positives) |
| `tcpost2.py` | `tcpost2_out.txt` | classification of the 200 minimisers | **evidence** |
| `control_pal7.py`, `control_model.py` | `c_*_pool.txt`, `c_*_leaves.txt` | round-13 models with correct base widths (§4.5) | exact |
| `csearch4w.c` | (stdout) | `csearch4` with 4608-bit sets (`M8'`, `m = 11`) | exact |
| `memsmall.py` | `mem_*_out.txt` | membership plus exact rank and zonotope checks in the `n = 8` models | exact (negatives via `rank_Q = rank_p`; positives are zonotopes) |
| `csearch9.c` | — | `csearch8` plus exact column-requirement pruning | exact |
