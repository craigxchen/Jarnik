# Round 28, lane compute: reach of the multi-shift criterion vs zonotope order on growing slab boxes

Status: complete (unrefereed). Labels: **exact** (over Q / Z), **mod-p** (ranks mod two primes > 2^30; evidence),
**evidence**, **proved**. Everything here is unrefereed.

## 0. Verdict and decision-rule reading

**No growth statement is proved or claimed.** `M <= (2/3+eps) log R/loglog R` stands. This is a data report on the
round-27 pre-registered test; the decision rule (README §11) is applied as fixed in advance.

1. **Reproduction (own pipeline, validated against the killer's exact-`Q` code): all round-27 values reproduce.**
   sdc7 `L = 2`: reach_pair 6, `zeta` 7; syl8 `L = 2`: reach_pair 6, reach_cube 7, `zeta` 7; M8' = syl8 `L = 3`:
   reach_pair 11 = `rho`, `zeta` 10 (§3). The killer's 27 random boxes give the same 20/7/0 split, now with exact
   `zeta` (all 27 are zonotopal, `zeta = rho`).
2. **Extension: 105 boxes in total** (`d = 5..9`; 28 structured, 77 random), reach mod two primes (identical on
   every box), `zeta` exact on every box:
   * `d = 8`: 18 structured boxes (Sylvester-16 / sdc11 / Paley-Hadamard-12 subcluster profiles, `L = 1, 2, 3`) and
     30 random boxes (killer law, seeds 1–30);
   * `d = 9`: 4 structured (`L = 1, 2`) and 20 random (`L = 1`); `d = 9` at `L >= 3` and random `L = 2` are beyond
     the memory cap (`|A| > 4900`, `H_A` needs a `|A| x N_{<=rho}` elimination).
3. **Outcome.** `reach > zeta` on **2 of 105** boxes: M8' and its one-point extension syl16|J9a (`d = 8`, `L = 3`:
   reach 11 = `rho` = truth, `zeta = 10` exact although the LP bound is 11). Every other box has `reach <= zeta`:
   ties on the other four Sylvester `L = 3` restrictions and on the Hadamard `L = 2` boxes (`d = 8, 9`), one short on
   the Hadamard and sdc `L = 3` boxes, `reach < zeta` on all 30 random `d = 8` boxes and on 19 of 20 random `d = 9`
   boxes. The gap `zeta - reach` **grows** with `rho/L` (mean 0.6 at
   `rho/L in (3,4]` to 2.2 at `rho/L > 6`) and, for random boxes, with `d` (0.73, 1.17, 1.87 at `d = 5, 6, 8`).
4. **Decision rule: UNFAVOURABLE.** The favourable branch needs `reach > zeta` on a growing share as `d` goes
   7 → 8 → 9 and `rho/L` increases; the share is 1/6, 1/48, 0/24, and is 0 at every `rho/L > 11/3`. The unfavourable
   branch ("`reach <= zeta` on essentially every box beyond M8'") holds: the single exception beyond M8' is M8''s own
   point configuration plus one point, and the four other `d = 8` Sylvester restrictions at the same `L` tie.
   Per the pre-registered rule, Cor 4.1 is to be treated as presumptively vacuous and the internal route-5 rounds stop;
   Questions A and B, ZE^all_H, the dP2 instance and the `Mbar_{0,7}` reformulation go to experts.
5. **Caveats** (§7): mod-`p` evidence for every reach; cube shift sets truncated at `d = 8` (random) and `d = 9`;
   all models are LP-subcritical (`rho/L <= 8`, with large `rho/L` only at `L = 1`), far from the regime
   (`m'/j = 8` at `j -> infinity`). The data cannot exclude a change of behaviour at regime size; they show no sign of
   one.


## 1. Definitions and conventions (checked against the record)

All definitions are those of `round27/interp.md` §§1, 3 and `round27/kill_interp.md` §5, in the box convention of
round 11 (`round11/compute_checks/bigger.py`, `core.py`) and of the round-27 killer (`kill_interp_checks/klib.py`):

* **Slab box.** `n` points, `d = n - 1`, frequencies `k' in {-1,0,1}^d`, `kfull = (k', -sum k')` (base point `n-1`).
  A *profile* is a set of canonical subsets `S_p` of `{0..n-1}` (a subset and its complement are identified; the
  canonical one is the shorter, then lexicographically smaller). The box is
  `A = { k' in {-1,0,1}^d : |kfull(S_p)| <= L for every p }`.
  For a sign column `beta_p` (entries `+-1`) one has `<k, beta_p> = +-2 kfull(S_p)`, so this is the task's
  `|<k, beta_p>| <= 2L`; the record's "`L`" (and ours) is the bound on `|k(S_p)|`, the analogue of `j = floor(w/2)`
  in the regime. All boxes here have one common `L` for all cuts (`T = 1`).
* `V(S,e)`, `M_e(S)` (degree-`e` least space), `H_S(e) = dim M_e(S)` (graded Hilbert function), as in interp §1.
  `rho(A) = max{e : H_A(e) > 0}` (the top degree; record's `rho`, "regularity"). The record's LP bound is a bound
  for `zeta`, not for `rho`; it is reported separately as `zeta_LP`.
* `(SHIFT_e)` w.r.t. a shift set `L`: `dim sum_{l in L} M_e(A cap (A - l)) = H_A(e)` (interp (3.2)).
  **reach(A, L)** = the first `e` at which `(SHIFT_e)` fails (= the largest `m` with `(SHIFT_e)` for all `e < m`);
  Thm 3.1 then certifies ZE-braid at order `m = reach`.
  *Fact (used to stop at `rho - 1`):* `reach <= rho`, since `(SHIFT_e)` for all `e <= rho` would, by Thm 3.1(1)
  with `m = rho + 1`, put every `Y` with `Y^l != 1` outside `Bs(V(A, rho+1)) = Bs({0}) = G_m^d`.
* Shift sets: `L_pair = {+-e_x, +-(e_x - e_y)}` (`d(d+1)` shifts); `L_cube = {-1,0,1}^d \ 0`. Since
  `A cap (A + l) = (A cap (A - l)) + l` and leading forms are translation invariant, `M_e(A_{-l}) = M_e(A_l)`; only
  one shift per class `{l, -l}` is computed (`d(d+1)/2` resp. `(3^d - 1)/2` classes).
* **zeta.** Order of a lattice zonotope form `c = prod_i (T_{v_i} - 1)^{a_i} delta_{k0}` (order `t = sum a_i`,
  `P_c = Y^{k0} prod (Y^{v_i} - 1)^{a_i}`), all vertices in `A`. Two generator classes:
  `zeta_all` (all primitive `v in {-2..2}^d`, the class of round 11 Prop 2.1, the killer's `k6_zono.py` and of the
  value `zeta(M8') = 10` quoted in the task), and `zeta_pair` (generators in `L_pair`, the task's wording).
  Both are reported. They coincide on every structured box; on 8 of the 77 random boxes `zeta_pair = zeta_all - 1`.
  The comparison below uses `zeta_all` (the larger, i.e. the harder bar); against `zeta_pair` the counts change only
  where noted (against `zeta_pair`: two more ties, `n = 9` seed 27 and `n = 10` seed 25; no new `reach > zeta`).
* **truth** = the largest `m` with no nondegenerate point in `Bs(V(A,m))`. Always
  `max(reach, zeta) <= truth <= rho` (a zonotope form of order `zeta` is nonzero at every nondegenerate point;
  `V(A, rho + 1) = 0`). It is not computed beyond these bounds except where the record decides it.
  *Refinement (proved here, used in `boxes.csv`):* if `H_A(rho) = 1` and `zeta_all < rho`, then `truth <= rho - 1`.
  *Proof.* `V(A, rho) = <c>`. If every zero of `P_c` were degenerate, Thm 2.2's classification (interp (2.3)) gives
  `P_c = c' Y^a prod_j (Y^{v_j} - zeta_j)^{m_j}` with `v_j` primitive in `{-2..2}^d` and
  `ord_e P_c = sum_{zeta_j = 1} m_j = rho`. By Ostrowski `Newt(P_c) = a + sum_j m_j [0, v_j]`, which lies in `conv(A)`;
  so the lattice zonotope `a + sum_{zeta_j = 1} m_j [0, v_j]` (a Minkowski summand containing the vertex `a`) has its
  vertices among the lattice points of `conv(A)`, i.e. in `A` (`A = conv(A) cap Z^d` for a slab box), and is a zonotope of
  order `rho` in `A`, contradicting `zeta_all < rho`. So `Bs(V(A, rho))` contains nondegenerate points. QED
  (Here `H_A(rho) = 1` is a mod-`p` value, both primes.)
* `rho/L` is the model analogue of the regime ratio `m'/j ~ 8` (regime: `m' = ceil(4w) = 8j`, zonotopes capped at
  `2w(1+o(1)) = 4j`, so the decisive window is `e/L in (4, 8)`).

## 2. Pipeline and its validation

**Least spaces (own code, `compute_checks/slib.py`).** On the ternary cube `k_x^3 = k_x`, so every monomial equals,
on `S`, a standard monomial (exponents `<= 2`) of no larger degree. Hence `V(S,e)` is cut out by the standard
monomials of degree `< e`, and for `c in V(S,e)` the degree-`e` moment of a non-standard monomial is a moment of
order `<= e - 2`, i.e. zero; so `M_e(S)` lives in the span of the `N_e` standard monomials of degree `e`
(`N_e = [x^e](1+x+x^2)^d`, at most 1107 at `d = 8`, 3139 at `d = 9`).
Let `E` be the `|S| x N_{<=E}` evaluation matrix with columns grouped by degree and `R = rref(E)` (mod `p`).
*Lemma (row-echelon extraction).* `M_e(S)` is spanned by the degree-`e` blocks of the rows of `R` whose pivot lies
in degree block `e`, and these are independent. *Proof.* Row space of `R` = `{c^T E}` (moment vectors). A combination
of rows of `R` vanishes on the blocks of degree `< e` iff it gives coefficient 0 to every row with pivot in such a
block (look at the pivot columns, which are unit columns); the remaining rows with pivot in a block `> e` vanish on
block `e`. Distinct pivots give independence. QED
So one `rref` per set gives `H_S(e)` and `M_e(S)` for all `e`. Sums `sum_l M_e(A_l)` are accumulated per degree
(lazy re-reduction), shift classes in decreasing order of `|A_l|`, and a degree is closed as soon as the sum reaches
`H_A(e)` (only the failing degree needs all classes). Shifts with `|A_l| <= e` are skipped (`H_{A_l}(e) = 0`).
For the cube shift set the computation continues from the `L_pair` accumulators (`L_pair` is a subset of `L_cube`)
on the degrees from the `L_pair` failure on.

**Primes.** `P1 = 2147483629`, `P2 = 2147483587` (both prime, `> 2^30`). Every reach value in this report was
computed independently mod both primes; all mod-`p` statements are **evidence** (a mod-`p` rank can differ from the
rational one for finitely many `p`). They agreed on every box (§§3-5).

**Zonotopes (`compute_checks/zono.py`).** Candidates: all primitive `v in {-2..2}^d` (sign-normalised) with
`|lam(v)| <= 2 b_lam` for every box functional (`lam` = a coordinate, `b = 1`; or a slab, `b = L`).
*Upper bound (exact):* the record's LP (`round11/compute_checks/lpzeta.py`): any zonotope in `A` has
`sum_i a_i |lam(v_i)| <= 2 b_lam`; a rational dual `y >= 0` with `sum_lam y_lam |lam(v)| >= 1` for **every**
candidate (checked in exact integer arithmetic over the full candidate list) gives `zeta <= sum 2 b_lam y_lam`.
*Exact optimum (`zono.dfs_feasible`).* Every zonotope of order `T` satisfies
`sum_v x_v c_v <= B` with `c_v = sum_lam y_lam |lam(v)| >= 1` (exact rationals), hence `sum_v x_v (c_v - 1) <= B - T`:
only generators with `c_v - 1 <= B - t` can occur in a zonotope of order `>= t`. For `t = floor(B), floor(B) - 1, ...`
an exhaustive depth-first search over these generators (multiplicities 1, 2; the set of admissible base points `k0`
is carried along and shrinks monotonically; capacity and excess pruning) decides whether a zonotope of order `>= t`
exists. The first feasible `t` is `zeta`; the witness is then re-verified **exactly**: support in `A`, all moments of
order `< t` vanish (all standard monomials of degree `< t`, integer arithmetic), one order-`t` moment is nonzero (the
coefficient of the lex-leading monomial of `prod <v_i,u>^{a_i}`). Every `zeta` in this report is **exact** in this sense
(upper bound = exact LP bound, or exact LP bound plus a completed exhaustive search; no timeouts occurred).
*Cross-checks:* the first version used a HiGHS MILP over all candidates (it ran out of memory on some `d = 8` boxes);
on the 39 boxes where it finished, MILP optimum = DFS value for both classes. Where `floor(LP)` is not attained (the
`d = 8` Sylvester restrictions at `L = 3`, §4) the "no zonotope of order 11" verdict was re-checked for J9a, J9b and
one random-`J` box by an independent HiGHS feasibility MILP on the dual-reduced candidate set (`xcheck_zeta.py` →
`xcheck_out.txt`): infeasible in all six (box, class) cases, as the DFS says.

**Validation.**
* `v0_validate.py` → `v0_validate_out.txt`: on the cube `d = 3`, the hexagon and 10 random slab boxes (`d = 3, 4`),
  every degree, `slib` (mod `P1` and `P2`, one sign per class) agrees with the round-27 killer's independent library
  `klib.shift_dims_Q` (**exact over `Q`**, both signs of every shift), including the non-spanning dimensions
  (e.g. 10 vs 11 and 3 vs 4). 0 disagreements. This also checks `M_e(A_{-l}) = M_e(A_l)` in practice.
* Boxes: `|A|` = 757 (sdc7 `L = 2`), 1037 (syl8 `L = 2`), 1877 (M8'), 1765 (sdc7 `L = 3`), as in the record; the
  Hilbert functions reproduce round 27/round 11 exactly (§3).
* `t1_zeta_small.py` → `t1_zeta_small_out.txt`: `zeta` reproduces the record/killer values (§3).

## 3. Reproduction of the round-27 values

Own pipeline (`slib.py`, `zono.py`), run by `run_box.py` / `run_zeta.py` (outputs `results_<p>.jsonl`, `zeta.jsonl`,
logs in `compute_checks/logs/`). Profiles: `models.py` (Sylvester `H_8`; skew double core `q = 7`, the record's
construction; both rebuilt here and checked against the record's `|A|`). Reach columns: mod `P1` / mod `P2`.

| box | `|A|` | `H_A` (both primes) | `rho` | reach `L_pair` | first failure (sum/`H`) | reach `L_cube` | `zeta_all` = `zeta_pair` (LP bound) | round 27 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sdc7 `L = 2` | 757 | 1 7 28 77 161 251 188 43 1 | 8 | **6 / 6** | `e = 6`: 166/188 | **6 / 6** | **7** exact (LP 7) | 6, 6, 7 ✓ |
| syl8 `L = 2` (M8) | 1037 | 1 7 28 77 161 259 287 175 42 | 8 | **6 / 6** | `e = 6`: 286/287 | **7 / 7** (`e = 7`: 141/175) | **7** exact (LP 7) | 6, 7, 7 ✓ |
| syl8 `L = 3` (M8') | 1877 | 1 7 28 77 161 266 357 386 329 195 63 7 | 11 | **11 / 11** | none below `rho` | 11 / 11 | **10** exact (LP 21/2) | 11, 10 ✓ |
| sdc7 `L = 3` | 1765 | 1 7 28 77 161 266 357 379 301 153 35 | 10 | 10 / 10 | none below `rho` | 10 / 10 | 10 exact (LP 21/2) | (record: `zeta = 10`) |
| syl8 `L = 1` | 57 | 1 7 28 21 | 3 | 3 / 3 | — | 3 / 3 | 3 exact (LP 7/2) | 3 ✓ |
| sdc7 `L = 1` | 57 | 1 7 28 21 | 3 | 3 / 3 | — | 3 / 3 | 3 exact (LP 7/2) | — |

**All three requested round-27 values reproduce**: sdc7 `L = 2` (reach_pair 6, `zeta` 7), syl8 `L = 2` (reach_pair 6,
reach_cube 7, `zeta` 7), M8' (reach_pair 11, `zeta` 10); the failure dimensions agree with the killer's
(`k3_*`, `k10_syl8_e7`). Additions:
* M8' contains an order-10 zonotope **with all generators in `L_pair`**: `k0 = (-1,-1,-1,0,-1,0,1)`, generators
  `e_0, e_1, e_2 (twice), e_3, e_4, e_0 - e_6, e_1 - e_3, e_4 - e_5, e_5 - e_6` (exact: 540 support points in `A`,
  moments of order `< 10` vanish, an order-10 moment is nonzero, `|c|_inf = 2`). So `zeta_pair(M8') = 10` too (the
  killer's witness used the non-pair generator `e_0 - e_1 - e_2`). LP bound `21/2` makes `zeta = 10` exact.
* sdc7 `L = 3` (`rho = 10`): reach 10 = `rho` = `zeta`; the criterion ties with zonotopes there.
* **The 27 random boxes of round 27** (killer's law and seeds: `n = 6`, seeds 100–114; `n = 7`, seeds 1–12),
  recomputed with exact `zeta` (the killer had only `zeta >= reach + 1`): reach w.r.t. `L_pair` and `L_cube`
  coincide on all 27; `zeta_all > reach` on **20**, `=` on **7**, `<` on **0** — exactly the killer's 20/7/0.
  The gaps are `zeta - reach = 1` (16 boxes) and `2` (4 boxes: `n=6` seed 101; `n=7` seeds 7, 8, 10). Full rows in
  `boxes.csv`.
* `zeta_pair` can be smaller than `zeta_all`: on `n = 7` seeds 7 and 8, `zeta_pair = 5 < zeta_all = 6` (the order-6
  witnesses need a non-pair generator). Even against `zeta_pair` the criterion loses there (reach 4).

## 4. Larger structured boxes (Sylvester, sdc) at d = 8, 9

**What "Sylvester / sdc at `d = 8, 9`" means here.** Neither family exists at `n = 9, 10` (Hadamard orders are
`1, 2, 4k`; the skew double core needs `q = n - 1 = 3 mod 4` prime). The computable members are `n = 8` (`d = 7`) and
`n = 12` (`d = 11`: sdc11 / Paley-Hadamard 12, `|A| = 9373` at `L = 2`, `N_{<=7} ~ 1.7*10^4` standard monomials; out
of reach of this pipeline, and `L = 1` gives the 133-point box with `rho = 3`). The faithful substitute is the
**subcluster profile**: if an actual cluster has a Sylvester-16 / sdc11 / Paley-Hadamard-12 profile, each 9- or
10-point subcluster `J` has the restricted profile `S_p -> S_p cap J` (trivial columns dropped, canonicalised,
deduplicated; `models.restrict`). Boxes used (point indices refer to the record's indexing: `F_2^4` in binary order
for Sylvester; `(inf, 0, 1, ..., q-1)` for Paley/sdc; the last element of `J` is the base point):

| name | profile | `J` | cuts |
| --- | --- | --- | --- |
| syl16\|J9a | Sylvester `H_16` | `0..8` (an affine hyperplane `F_2^3 x {0}` plus one point) | 15 |
| syl16\|J9b | Sylvester `H_16` | `0,1,2,4,8,15,3,5,6` (not containing a hyperplane) | 15 |
| sdc11\|J9 | skew double core `q = 11` | `inf,0..7` | 22 |
| had12\|J9 | Paley-I Hadamard `q = 11` | `inf,0..7` | 11 |
| syl16\|J10, sdc11\|J10, had12\|J10 | as above | `0..9` resp. `inf,0..8` | 15, 22, 11 |

All restricted profiles are full rank. At `L = 1` all four `d = 8` profiles give the same 73-point box
(`H = 1 8 36 28`, `rho = 3`): uninformative (reach = `zeta` = `rho` = 3).

**Results at `d = 8`** (reach mod `P1` / mod `P2`; `zeta` exact; `zeta_all = zeta_pair` on every row):

| box | `L` | `|A|` | `H_A` | `rho` | `rho/L` | reach `L_pair` (fail dims) | reach `L_cube` | `zeta` (LP bound) | reach − `zeta` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| syl16\|J9a | 2 | 1837 | 1 8 36 112 266 489 580 309 36 | 8 | 4 | 6 / 6 (572/580) | 7 / 7 | 7 (15/2) | −1 (cube: 0) |
| syl16\|J9b | 2 | 1841 | 1 8 36 112 266 490 577 310 41 | 8 | 4 | 6 / 6 (576/577) | 7 / 7 | 7 (15/2) | −1 (cube: 0) |
| sdc11\|J9 | 2 | 1367 | 1 8 36 112 266 479 381 84 | 7 | 3.5 | 6 / 6 (336/381) | 6 / 6 | 7 = `rho` (22/3) | −1 |
| had12\|J9 | 2 | 2049 | 1 8 36 112 266 492 680 412 42 | 8 | 4 | 7 / 7 (12/412) | 7 / 7 | 7 (22/3) | 0 |
| **syl16\|J9a** | **3** | **4861** | 1 8 36 112 266 504 784 1001 1029 756 315 49 | **11** | 3.67 | **11 / 11** (none below `rho`) | — | **10** (LP 11; no order-11 zonotope) | **+1** |
| sdc11\|J9 | 3 | 4457 | 1 8 36 112 266 504 784 998 984 606 158 | 10 | 3.33 | 9 / 9 (604/606) | not run | 10 = `rho` (1316/121) | −1 |
| had12\|J9 | 3 | 4997 | 1 8 36 112 266 504 784 1007 1044 817 368 50 | 11 | 3.67 | 10 / 10 (342/368) | not run | 11 = `rho` (LP 11) | −1 |

**The Sylvester `L = 3` box at `d = 8` repeats M8' exactly.** syl16|J9a at `L = 3`: `(SHIFT_e)` w.r.t. `L_pair`
holds for every `e <= 10` (both primes) and the reach is `rho = 11`, so Thm 3.1 certifies ZE-braid at the top order and
truth = 11. The exact LP bound is exactly `11`, but **no zonotope of order 11 lies in the box**: the dual reduction
leaves only the 36 `L_pair` classes (`c_v = 1` exactly; every other primitive generator has `c_v > 1`), and both the
exhaustive DFS (3 506 942 nodes, complete) and an independent HiGHS feasibility MILP on that reduced set
(`xcheck_zeta.py`; infeasible) find none. An order-10 `L_pair` zonotope exists (exact witness in `zeta.jsonl`), so
`zeta = 10` exactly and **reach − zeta = +1**, as in M8'.
**Every other faithful structured box is at or below `zeta`**: the `L = 2` Sylvester restrictions reproduce the M8
pattern (pair reach `zeta - 1`, cube reach `= zeta`), had12|J9 `L = 2` ties, and sdc11|J9 (`L = 2, 3`) and
had12|J9 `L = 3` are zonotopal (`zeta = rho`), so the criterion cannot beat `zeta` there and falls one short.

**Further `d = 8`, `L = 3` restrictions** (added to test whether the J9a success is specific to M8''s configuration;
J9b fixed in advance, the other six drawn with `random.Random(28).sample`; `jobs5.txt`):

| box | `|A|` | `H_A` | `rho` | reach `L_pair` (fail dims) | reach cube | `zeta` (LP) | reach − `zeta` | truth |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| syl16\|J9b `0,1,2,4,8,15,3,5,6` | 4787 | 1 8 36 112 266 504 784 1004 1026 735 281 30 | 11 | 10 / 10 (215/281) | 10 (300 largest classes, `P1`; dims unchanged) | 10 (LP 11, no order 11) | 0 | 10 or 11 |
| syl16\|`3,11,2,8,9,13,15,10,7` | 4809 | 1 8 36 112 266 504 784 1005 1034 746 284 29 | 11 | 10 / 10 (218/284) | — | 10 (LP 89/8) | 0 | 10 or 11 |
| syl16\|`13,3,14,2,6,12,10,15,8` | 4809 | same | 11 | 10 / 10 | — | 10 (89/8) | 0 | 10 or 11 |
| syl16\|`0,6,3,1,12,14,5,2,8` | 4809 | same | 11 | 10 / 10 | — | 10 (89/8) | 0 | 10 or 11 |
| had12\|`8,11,4,3,7,1,6,0,2` | 4981 | 1 8 36 112 266 504 784 1006 1045 815 359 45 | 11 | 10 / 10 (334/359) | — | 11 = `rho` (LP 11) | −1 | 11 |
| had12\|`10,11,4,1,9,5,0,7,3` | 4997 | 1 8 36 112 266 504 784 1007 1044 817 368 50 | 11 | 10 / 10 | — | 11 = `rho` | −1 | 11 |
| sdc11\|`2,3,5,1,0,10,6,11,8` | 4465 | 1 8 36 112 266 504 784 997 986 611 160 | 10 | 10 / 10 (= `rho`) | — | 10 = `rho` | 0 | 10 |

(The three random Sylvester 9-subsets give the same Hilbert function and LP bound; they are presumably equivalent under
`AGL(4,2)`, which was not checked. Their `zeta = 10` is exact: exhaustive DFS for order 11, 423 297 nodes each.)
So among the five `d = 8` Sylvester `L = 3` restrictions only the hyperplane-containing one beats `zeta`; the
other four tie at `zeta = 10 = rho - 1`. (Checked directly: J9a contains exactly one affine hyperplane of `F_2^4`, `{x_0 = 0}`; J9b and the three random
sets contain none.) The `reach > zeta` phenomenon at `d = 8` is therefore attached to M8''s point configuration
(Sylvester `H_8` on 8 of the points), not to the Sylvester-16 profile as such.

**Results at `d = 9`** (`n = 10`; `J = 0..9` for Sylvester, `inf,0..8` for Paley/sdc; reach mod `P1` / `P2`):

| box | `L` | `|A|` | `H_A` | `rho` | `rho/L` | reach `L_pair` (fail dims) | reach cube (300 largest classes) | `zeta` (LP) | reach − `zeta` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| syl16\|J10 | 1 | 91 | 1 9 45 36 | 3 | 3 | 3 / 3 | 3 | 3 (15/4) | 0 |
| syl16\|J10 | 2 | 3847 | 1 9 45 156 414 867 1251 891 213 | 8 | 4 | 6 / 6 (1250/1251) | 7 / 7 | 7 (15/2) | −1 (cube 0) |
| sdc11\|J10 | 2 | 2709 | 1 9 45 156 414 858 895 331 | 7 | 3.5 | 6 / 6 (893/895) | 6 / 6 | 7 = `rho` (22/3) | −1 |
| had12\|J10 | 2 | 4525 | 1 9 45 156 414 870 1446 1227 357 | 8 | 4 | 7 / 7 (1112/1227) | 7 / 7 | 7 (22/3) | 0 |

The `d = 9`, `L = 2` structured boxes repeat the `d = 7, 8` pattern exactly (Sylvester: pair `zeta - 1`, cube `zeta`;
sdc: zonotopal, one short; Hadamard: tie). `L = 3` at `d = 9` (`|A|` = 11561–13419) is out of reach (§0 item 2).


## 5. Random multi-cut slab boxes at d = 8 (and d = 9)

**Sampling law** (the round-27 killer's `k8_random.py`, reproduced exactly by `models.random_profile`; seed `s`,
Python `random.Random(s)`): `n` points; the number of cuts `r` is uniform on `{n-1, ..., 2n}`; cuts are drawn one at a
time, each point of `{0..n-1}` put in `S` independently with probability `1/2`, rejected if `S` is empty or everything,
canonicalised (`S` or its complement), duplicates rejected, until `r` distinct cuts; `L` uniform on `{1, 2}`; `T = 1`.
`d = 8`: `n = 9`, seeds `1..30` (all 30 boxes, no selection). `d = 9`: `n = 10`, the seeds in `1..40` that draw `L = 1`
(the law conditioned on `L = 1`; the `L = 2` boxes at `n = 10` have 4905–7909 points and are out of reach of the
`H_A` computation within the memory cap), first 20 such seeds.
Cube shifts: at `d = 8` the 300 largest non-pair classes (by `|A_l|`) after `L_pair`, degree by degree from the
`L_pair` failure (a lower bound for `reach(L_cube)`; see §7).

**Results, `d = 8` (30 boxes; reach identical mod `P1` and mod `P2` on all 30, `H_A` identical; `zeta` exact on all
30).** Columns: seed, cuts `r`, `L`, `|A|`, `rho`, `rho/L`, reach `L_pair`, reach with cube shifts (largest 300
classes for seeds 1–5, 100 for 6–30; lower bound), `zeta_all` (exact LP bound), `zeta_pair`, reach_pair − `zeta_all`.

| seed | `r` | `L` | `|A|` | `rho` | `rho/L` | reach pair | reach cube (≥) | `zeta_all` (LP) | `zeta_pair` | reach − `zeta` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 10 | 1 | 351 | 6 | 6.0 | 4 | 4 | 5 (40/7) | 5 | -1 |
| 2 | 8 | 2 | 3803 | 11 | 5.5 | 8 | 8 | 11 (11) | 11 | -3 |
| 3 | 11 | 2 | 2907 | 10 | 5.0 | 7 | 7 | 10 (10) | 10 | -3 |
| 4 | 11 | 2 | 3187 | 10 | 5.0 | 8 | 8 | 10 (21/2) | 10 | -2 |
| 5 | 17 | 2 | 2281 | 9 | 4.5 | 7 | 7 | 9 (19/2) | 9 | -2 |
| 6 | 17 | 2 | 2293 | 9 | 4.5 | 7 | 7 | 9 (9) | 9 | -2 |
| 7 | 13 | 2 | 2851 | 10 | 5.0 | 7 | 7 | 10 (10) | 10 | -3 |
| 8 | 11 | 1 | 363 | 6 | 6.0 | 4 | 4 | 6 (6) | 5 | -2 |
| 9 | 15 | 2 | 2377 | 9 | 4.5 | 7 | 7 | 9 (28/3) | 9 | -2 |
| 10 | 17 | 2 | 2187 | 8 | 4.0 | 7 | 7 | 8 (53/6) | 8 | -1 |
| 11 | 15 | 1 | 249 | 5 | 5.0 | 4 | 4 | 5 (11/2) | 5 | -1 |
| 12 | 15 | 2 | 2521 | 9 | 4.5 | 7 | 7 | 9 (19/2) | 9 | -2 |
| 13 | 12 | 1 | 303 | 5 | 5.0 | 4 | 4 | 5 (6) | 5 | -1 |
| 14 | 9 | 1 | 633 | 8 | 8.0 | 6 | 6 | 8 (8) | 8 | -2 |
| 15 | 11 | 1 | 365 | 6 | 6.0 | 4 | 4 | 6 (44/7) | 5 | -2 |
| 16 | 13 | 1 | 391 | 6 | 6.0 | 4 | 4 | 6 (19/3) | 6 | -2 |
| 17 | 16 | 1 | 121 | 4 | 4.0 | 3 | 3 | 4 (19/4) | 4 | -1 |
| 18 | 10 | 1 | 501 | 6 | 6.0 | 5 | 5 | 6 (20/3) | 6 | -1 |
| 19 | 18 | 1 | 261 | 5 | 5.0 | 4 | 4 | 5 (17/3) | 5 | -1 |
| 20 | 18 | 2 | 2233 | 9 | 4.5 | 6 | 6 | 8 (9) | 8 | -2 |
| 21 | 10 | 2 | 3307 | 10 | 5.0 | 7 | 7 | 10 (21/2) | 10 | -3 |
| 22 | 10 | 1 | 355 | 6 | 6.0 | 4 | 4 | 6 (6) | 5 | -2 |
| 23 | 12 | 1 | 453 | 7 | 7.0 | 4 | 4 | 7 (15/2) | 7 | -3 |
| 24 | 14 | 1 | 285 | 5 | 5.0 | 4 | 4 | 5 (17/3) | 5 | -1 |
| 25 | 14 | 1 | 257 | 5 | 5.0 | 4 | 4 | 5 (39/7) | 5 | -1 |
| 26 | 11 | 2 | 2859 | 10 | 5.0 | 8 | 8 | 10 (10) | 10 | -2 |
| 27 | 18 | 1 | 129 | 5 | 5.0 | 4 | 4 | 5 (26/5) | 4 | -1 |
| 28 | 9 | 2 | 3175 | 10 | 5.0 | 8 | 8 | 10 (41/4) | 10 | -2 |
| 29 | 16 | 2 | 2651 | 10 | 5.0 | 7 | 7 | 10 (10) | 10 | -3 |
| 30 | 16 | 1 | 173 | 5 | 5.0 | 3 | 4 | 5 (21/4) | 5 | -2 |

* **reach < `zeta` on 30 of 30** (pair and cube); `=` on 0; `>` on 0. Gap `zeta - reach`: 1 on 10 boxes, 2 on 14,
  3 on 6 (mean 1.87).
* 28 of 30 boxes are zonotopal (`zeta_all = rho`); the two faithful ones (seed 1: `zeta = 5 < rho = 6`; seed 20:
  `zeta = 8 < rho = 9`) have reach 4 and 6, i.e. still below `zeta`.
* Cube shifts changed the reach on 1 box of 30 (seed 30: 3 → 4), never up to `zeta`.
* The largest `rho/L` reached is 8 (seed 14, `L = 1`, `|A| = 633`): reach 6 = `rho - 2`, `zeta = rho = 8`.

**Results, `d = 9` (`n = 10`, `L = 1`, the 20 seeds in `1..38` drawing `L = 1`; reach and `H_A` identical mod both
primes on all 20; cube = 100 largest non-pair classes; `zeta` exact).** Columns: seed, `r`, `|A|`, `rho` (= `rho/L`),
reach pair, reach cube (≥), truth, `zeta_all` (LP), `zeta_pair`, reach_pair − `zeta`.

| seed | `r` | `|A|` | `rho` | reach pair | reach cube (≥) | truth | `zeta_all` (LP) | `zeta_pair` | reach − `zeta` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 9 | 1255 | 7 | 5 | 6 | 7 | 7 (15/2) | 7 | -2 |
| 5 | 18 | 177 | 4 | 3 | 3 | 4 | 4 (34/7) | 4 | -1 |
| 6 | 18 | 239 | 5 | 4 | 4 | 5 | 5 (11/2) | 5 | -1 |
| 8 | 12 | 419 | 6 | 4 | 4 | 6 | 6 (6) | 6 | -2 |
| 9 | 16 | 271 | 5 | 4 | 4 | 5 | 5 (163/30) | 5 | -1 |
| 11 | 16 | 265 | 5 | 4 | 4 | 5 | 5 (104/19) | 5 | -1 |
| 14 | 10 | 815 | 7 | 4 | 4 | 7 | 7 (7) | 7 | -3 |
| 15 | 12 | 329 | 5 | 4 | 4 | 5 | 5 (11/2) | 5 | -1 |
| 16 | 14 | 325 | 5 | 4 | 4 | 5 | 5 (11/2) | 5 | -1 |
| 17 | 17 | 377 | 5 | 4 | 4 | 5 | 5 (17/3) | 5 | -1 |
| 19 | 19 | 291 | 5 | 3 | 4 | 5 | 5 (28/5) | 5 | -2 |
| 20 | 20 | 255 | 5 | 4 | 4 | 5 | 5 (28/5) | 5 | -1 |
| 21 | 11 | 707 | 6 | 5 | 5 | 6 | 6 (20/3) | 6 | -1 |
| 22 | 11 | 711 | 6 | 5 | 5 | 6 | 6 (84/13) | 6 | -1 |
| 23 | 13 | 307 | 5 | 4 | 4 | 5 | 5 (16/3) | 5 | -1 |
| 25 | 15 | 269 | 5 | 4 | 4 | 5 | 5 (16/3) | 4 | -1 |
| 27 | 19 | 205 | 4 | 4 | 4 | 4 | 4 (98/19) | 4 | 0 |
| 30 | 17 | 273 | 5 | 4 | 4 | 5 | 5 (37/7) | 5 | -1 |
| 31 | 9 | 833 | 7 | 5 | 5 | 7 | 7 (29/4) | 7 | -2 |
| 38 | 19 | 309 | 6 | 4 | 4 | 6 | 6 (55/9) | 6 | -2 |

* **reach < `zeta` on 19 of 20, `=` on 1 (seed 27, `rho = 4`), `>` on 0.** All 20 are zonotopal (`zeta = rho` =
  truth). Cube shifts raise the reach on 2 boxes (seeds 2, 19), never to `zeta`.


## 6. Tables: reach - zeta, reach/rho vs rho/L and box size; trend

All boxes of §§3–5 (`boxes.csv`, one row per box; `trend.csv` holds the same groupings). "best" = the larger of the
`L_pair` reach and the (possibly truncated) cube reach. `zeta` = `zeta_all`. Generated by `md_tables.py`.

**Table 6.1 (by `d`, all boxes).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d = 5 | 15 | 0 | 0 / 5 / 10 | 0 / 5 / 10 | 0.00 | -0.73 | -0.73 | 0.869 | 1.000 | 4.40 |
| d = 6 | 12 | 0 | 0 / 2 / 10 | 0 / 2 / 10 | 0.00 | -1.17 | -1.17 | 0.819 | 1.000 | 5.00 |
| d = 7 | 6 | 3 | 1 / 3 / 2 | 1 / 4 / 1 | 0.17 | -0.17 | +0.00 | 0.917 | 0.943 | 3.50 |
| d = 8 | 48 | 10 | 1 / 10 / 37 | 1 / 12 / 35 | 0.02 | -1.29 | -1.23 | 0.807 | 0.977 | 4.59 |
| d = 9 | 24 | 2 | 0 / 3 / 21 | 0 / 4 / 20 | 0.00 | -1.17 | -1.04 | 0.785 | 0.990 | 5.10 |

**Table 6.2 (by `d` and kind).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d = 5, random | 15 | 0 | 0 / 5 / 10 | 0 / 5 / 10 | 0.00 | -0.73 | -0.73 | 0.869 | 1.000 | 4.40 |
| d = 6, random | 12 | 0 | 0 / 2 / 10 | 0 / 2 / 10 | 0.00 | -1.17 | -1.17 | 0.819 | 1.000 | 5.00 |
| d = 7, structured | 6 | 3 | 1 / 3 / 2 | 1 / 4 / 1 | 0.17 | -0.17 | +0.00 | 0.917 | 0.943 | 3.50 |
| d = 8, random | 30 | 2 | 0 / 0 / 30 | 0 / 0 / 30 | 0.00 | -1.87 | -1.83 | 0.741 | 0.991 | 5.23 |
| d = 8, structured | 18 | 8 | 1 / 10 / 7 | 1 / 12 / 5 | 0.06 | -0.33 | -0.22 | 0.916 | 0.954 | 3.53 |
| d = 9, random | 20 | 0 | 0 / 1 / 19 | 0 / 1 / 19 | 0.00 | -1.30 | -1.20 | 0.767 | 1.000 | 5.40 |
| d = 9, structured | 4 | 2 | 0 / 2 / 2 | 0 / 3 / 1 | 0.00 | -0.50 | -0.25 | 0.871 | 0.938 | 3.62 |

**Table 6.3 (faithful boxes only, `zeta < rho`, by `d`).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d = 7 | 3 | 3 | 1 / 0 / 2 | 1 / 1 / 1 | 0.33 | -0.33 | +0.00 | 0.833 | 0.886 | 3.89 |
| d = 8 | 10 | 10 | 1 / 5 / 4 | 1 / 7 / 2 | 0.10 | -0.40 | -0.20 | 0.834 | 0.889 | 4.08 |
| d = 9 | 2 | 2 | 0 / 1 / 1 | 0 / 2 / 0 | 0.00 | -0.50 | +0.00 | 0.812 | 0.875 | 4.00 |

**Table 6.4 (by `rho/L`, all boxes).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rho/L <= 3 | 7 | 0 | 0 / 7 / 0 | 0 / 7 / 0 | 0.00 | +0.00 | +0.00 | 1.000 | 1.000 | 3.00 |
| 3 < rho/L <= 4 | 39 | 13 | 2 / 15 / 22 | 2 / 19 / 18 | 0.05 | -0.56 | -0.46 | 0.880 | 0.964 | 3.82 |
| 4 < rho/L <= 5 | 38 | 1 | 0 / 1 / 37 | 0 / 1 / 37 | 0.00 | -1.45 | -1.39 | 0.781 | 0.997 | 4.93 |
| 5 < rho/L <= 6 | 15 | 1 | 0 / 0 / 15 | 0 / 0 / 15 | 0.00 | -1.73 | -1.73 | 0.715 | 0.989 | 5.97 |
| rho/L > 6 | 6 | 0 | 0 / 0 / 6 | 0 / 0 / 6 | 0.00 | -2.17 | -2.00 | 0.696 | 1.000 | 7.17 |

**Table 6.5 (by box size `|A|`, all boxes).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| |A| < 100 | 18 | 0 | 0 / 11 / 7 | 0 / 11 / 7 | 0.00 | -0.39 | -0.39 | 0.917 | 1.000 | 3.89 |
| 100 <= |A| < 500 | 44 | 1 | 0 / 4 / 40 | 0 / 4 / 40 | 0.00 | -1.25 | -1.20 | 0.778 | 0.996 | 5.07 |
| 500 <= |A| < 2000 | 15 | 5 | 1 / 1 / 13 | 1 / 4 / 10 | 0.07 | -1.13 | -0.87 | 0.800 | 0.961 | 5.23 |
| |A| >= 2000 | 28 | 9 | 1 / 7 / 20 | 1 / 8 / 19 | 0.04 | -1.32 | -1.29 | 0.829 | 0.966 | 4.23 |

**Table 6.6 (random boxes by `d` and `L`).**

| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| d = 5, L = 1 | 10 | 0 | 0 / 3 / 7 | 0 / 3 / 7 | 0.00 | -0.70 | -0.70 | 0.856 | 1.000 | 4.70 |
| d = 5, L = 2 | 5 | 0 | 0 / 2 / 3 | 0 / 2 / 3 | 0.00 | -0.80 | -0.80 | 0.896 | 1.000 | 3.80 |
| d = 6, L = 1 | 8 | 0 | 0 / 2 / 6 | 0 / 2 / 6 | 0.00 | -1.12 | -1.12 | 0.804 | 1.000 | 5.38 |
| d = 6, L = 2 | 4 | 0 | 0 / 0 / 4 | 0 / 0 / 4 | 0.00 | -1.25 | -1.25 | 0.850 | 1.000 | 4.25 |
| d = 8, L = 1 | 16 | 1 | 0 / 0 / 16 | 0 / 0 / 16 | 0.00 | -1.50 | -1.44 | 0.727 | 0.990 | 5.62 |
| d = 8, L = 2 | 14 | 1 | 0 / 0 / 14 | 0 / 0 / 14 | 0.00 | -2.29 | -2.29 | 0.756 | 0.992 | 4.79 |
| d = 9, L = 1 | 20 | 0 | 0 / 1 / 19 | 0 / 1 / 19 | 0.00 | -1.30 | -1.20 | 0.767 | 1.000 | 5.40 |

**What the trend is (plain statement).**
1. **`reach > zeta` does not occur on a growing share of boxes.** It occurs on exactly two boxes: M8' (`d = 7`) and
   syl16|J9a at `L = 3` (`d = 8`), whose 9-point set contains an affine hyperplane of `F_2^4`, i.e. M8''s own point
   configuration plus one point. The other four `d = 8` Sylvester-16 restrictions at `L = 3` (J9b and three random
   9-subsets) all **tie** (`reach = zeta = 10 < rho = 11`; cube shifts on J9b, 300 largest classes, change nothing:
   215/281 at `e = 10` before and after). Share with `reach > zeta`: `d <= 6`: 0/27; `d = 7`: 1/6; `d = 8`: 1/48;
   `d = 9`: 0/24. Among faithful boxes (`zeta < rho`, the only ones where `reach > zeta` is possible): `d = 7`: 1/3;
   `d = 8`: 1/10; `d = 9`: 0/2.
2. **`reach - zeta` decreases with `rho/L`** (Table 6.4): mean `-0.6` for `3 < rho/L <= 4`, about `-1.5` for
   `4 < rho/L <= 5`, `-1.7` for `5 < rho/L <= 6`, `-2.2` for `rho/L > 6`; every box with `rho/L > 4` has
   `reach < zeta` except one tie. `reach/rho` falls from `~0.9` to `~0.7` over the same range, while `zeta/rho` stays
   `~1`. The only `reach > zeta` boxes have `rho/L = 11/3`.
3. **With box size and `d`, the random boxes get worse, not better** (Tables 6.1, 6.2, 6.6): mean `reach - zeta` is
   `-0.73` (`d = 5`), `-1.17` (`d = 6`), `-1.87` (`d = 8`; `-2.29` at `L = 2`), `-1.30` (`d = 9`, `L = 1` only).
   0 of 77 random boxes have `reach > zeta`; 69 have `reach < zeta`.
4. **Random boxes are almost all zonotopal** (`zeta = rho` on 75 of 77): for them `reach > zeta` is impossible a
   priori, and the criterion stops 1–3 degrees below the top. This is a property of the sampling law (few cuts relative
   to `2^n`), not of the criterion; but the faithful random boxes (`n = 9` seeds 1 and 20) also have `reach < zeta`,
   and there `truth = zeta` exactly (refinement lemma, `H_A(rho) = 1`): the criterion misses true statements that a
   single product form proves.
5. **Structured boxes** (all 28): the criterion is within one degree of `zeta` (pair reach `zeta - 1` or `zeta`; cube
   reach `= zeta` on every Sylvester `L = 2` box, `d = 7, 8, 9`), above it only on the two Sylvester `L = 3` boxes above.


## 7. Caveats

* **Mod-`p` evidence.** Every reach value is a mod-`p` rank statement (two primes `> 2^30`, independent runs, always
  in agreement). It is not a proof over `Q`. All `zeta` values are proved (exact witness; exact rational LP dual,
  plus an exhaustive search where the LP bound is not attained). The mod-`p` caveat also applies to `H_A`, hence to
  `rho` and to the truth bracket that uses `H_A(rho) = 1`.
* **Truncated cube shift sets.** At `d = 8` (random) and `d = 9` the cube reach was computed with the largest
  `K` shift classes only (`K` stated per row; ordering by `|A_l|`). A truncated value is a **lower bound** for
  `reach(L_cube)` when the sum did not span; when it spans, the value is exact for that degree.
* **LP-subcritical models.** No computable box reaches the regime ratio: here `rho/L <= 8` only at `L = 1` (where the
  boxes are tiny or zonotopal) and `rho/L ~ 4` at `L = 2, 3`, against `m'/j = 8` with `zeta/j <= 4 + o(1)` in the
  regime. Nothing here is a statement about regime boxes; the decision rule was fixed in advance to read model data.
* **Box convention.** `L` bounds `|k(S_p)|` (record convention); in the `beta` normalisation the slab half-width is
  `2L`.
* **Sampling.** The random law is the round-27 killer's (stated in §5); it is not the law of actual cluster profiles.
  At `d <= 6` it produces only zonotopal boxes (`zeta = rho`), where `reach > zeta` is impossible a priori (§3, §6).
* `truth` is only bracketed (`max(reach, zeta) <= truth <= rho`) except where reach or `zeta` equals `rho`.

## Appendix: scripts, outputs, CSV

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/`.
Heavy jobs ran under `prlimit --as=3500000000`, two at a time (one per prime); `zeta` jobs under `--as=2000000000`/
`2500000000`.

| file | content |
| --- | --- |
| `boxes.csv` | one row per box: spec, family, `d`, `n`, cuts, `L`, `|A|`, `rho`, `rho/L`, `H_A`, reach (pair, cube) mod `P1` and `P2`, first-failure dims, cube completeness, `zeta_all`/`zeta_pair` with exact LP bounds, reach − `zeta`, reach/`rho`, truth bracket, `H_A(rho)` |
| `trend.csv` | grouped statistics (by `d`, family, kind, `rho/L`, `|A|`, `L`, faithfulness) |
| `compute_checks/slib.py` | boxes, standard monomials, least spaces via one `rref`, shift accumulation (mod `p`) |
| `compute_checks/zono.py` | zonotope candidates, exact LP-dual bound (constraint generation), exact DFS, exact witness verification |
| `compute_checks/models.py` | Sylvester, sdc (`q = 7, 11`), Paley-Hadamard, restrictions to subclusters, the random law |
| `compute_checks/run_box.py` | driver: `H_A`, reach w.r.t. `L_pair`, then cube shifts degree by degree → `results_<p>.jsonl` |
| `compute_checks/run_zeta.py` | driver: `zeta_pair`, `zeta_all` → `zeta.jsonl` |
| `compute_checks/v0_validate.py` → `v0_validate_out.txt` | validation against the killer's exact-`Q` library |
| `compute_checks/t1_zeta_small.py` → `t1_zeta_small_out.txt` | first `zeta` reproduction (MILP version) |
| `compute_checks/xcheck_zeta.py` | independent MILP check of "no order-`t` zonotope" on the dual-reduced set |
| `compute_checks/table.py`, `trend.py`, `md_tables.py` | CSV and tables |
| `compute_checks/worker.sh`, `chain2.sh`, `jobs*.txt` | job queues (box specs exactly as in `boxes.csv`) |
| `compute_checks/results_2147483629.jsonl`, `results_2147483587.jsonl`, `zeta.jsonl` | raw results (last line per spec is authoritative) |
| `compute_checks/logs/` | run logs (per-shift progress, timings) |
| `compute_notes.md` | working notes |

**Box specs** (reproducible from `run_box.parse`): `syl:k:L` (Sylvester `H_{2^k}`), `sdc:q:L`, `had:q:L`;
`rsyl:4:J:L`, `rsdc:11:J:L`, `rhad:11:J:L` (restriction to the ordered point list `J`, last = base);
`rand:n:seed` (killer law). Every cut list is also stored verbatim in `results_*.jsonl` (`pats`).

