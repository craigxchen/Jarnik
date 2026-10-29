# Round 28 verification (independent, adversarial)

## 0. Verdict

**No growth statement is proved or claimed by either lane, and none follows from anything here.**
`M <= (2/3+eps) log R/loglog R` stands. All code is mine (`verify_checks/`). Reach is computed mod two primes that
differ from the primary's, and zeta exactly.

1. **Data confirmed (no numerical discrepancy).** I recomputed the following, and every `|A|`, `H_A`, reach,
   failure dimension, `zeta_pair`, `zeta_all` and LP bound agrees exactly with `compute.md` / `boxes.csv`:
   * the three round-27 reference values;
   * both boxes with reach > zeta;
   * 9 further boxes of the primary's set, including the largest computed ones at `d = 8` (4997, 4787 points) and
     `d = 9` (4525, 3847), plus one new box (J9a');
   * reach on the 27 round-27 random boxes.
   
   Details in §§2–5.
2. **Upgrades to exact.** Two values are now **certified over `Q`** (§6), not only mod `p`:
   * reach(M8', `L_pair`) = 11. So the killer's M8' value and the theory lane's Cor 2.2 are theorems.
   * reach(M8) = 6, with a defect that is exactly one-dimensional and spanned by `sigma^6`.
   
   J9a's reach 11 stays mod-`p` evidence, now at four independent primes.
3. **zeta certificates verified** (§5).
   * `zeta(M8') = zeta(J9a) = 10` for both generator classes. Exact witnesses; exact LP duals over the full
     primitive candidate set.
   * For J9a (LP bound exactly 11), my own exhaustive search over the dual-reduced set (exactly the 36 `L_pair`
     classes) finds no order-11 zonotope.
4. **New adversarial box: J9a'.** It is the `d = 8` box that contains the M8' box as a coordinate face: the same
   9 points as J9a, with the base moved into the hyperplane. Its reach is **9 < zeta = 10** (both primes; 740/741 at
   `e = 9`, unchanged with the 150 largest cube-shift classes). The primary's J9a does **not** contain the M8' box as a face; its face is `M8' cap {|k_base| <= 1}`.
   So the one `d = 8` success depends on where the base point is placed, and it does not propagate from M8' to its
   natural one-point extension. This makes the success more fragile than `compute.md` says. It supports the
   unfavourable reading.
5. **Decision rule.** I agree with the reading **unfavourable**, with stated caveats (§8 items 3–5).
   * The favourable branch (a growing share of reach > zeta as `d` = 7 → 8 → 9 and `rho/L` grow) is clearly not met.
   * Read literally, the unfavourable branch ("reach <= zeta on essentially every box beyond M8'") has exactly one
     exception, J9a, a `d = 8` restriction containing M8''s point configuration. Its base-moved twin J9a' fails.
   * Caveat: the `d = 9` data contain no `L = 3` box, the only `L` at which the phenomenon was ever seen.
6. **Theory lane.** Every proof I checked is correct; there are minor gaps, none affecting a conclusion (§7). Every
   numerical claim I re-ran reproduces (rigidity data at `n = 8, 16`, the 3-of-27 random BCO cases, `sigma^6`, the
   subset reach 7, the linkage identity, the product theorem on samples). It correctly states that it does not
   settle the rule.

## 1. Definitions as implemented (own code, `verify_checks/vlib.py`, `zono_v.py`)

Nothing is imported from round 27 or from the primary/theory lanes. The definitions below were taken from
`round27/interp.md` §§1, 3 and the round-27 killer's random law (`k8_random.py`), and checked against the two lanes.

* **Box.** `n` points `0..n-1`, base `n-1`, `d = n-1`, `k' in {-1,0,1}^d`, `kfull = (k', -sum k')`. A profile is a set
  of subsets `S` of the points (a set and its complement identified). `A = {k' : |kfull(S)| <= L for every S}`. As
  both lanes say, the bound `L` is on `|k(S)|`, i.e. half-width `2L` in the `beta` normalisation. The task's wording
  "`|<k, beta_p>| <= L`" is the same box up to this factor 2. No discrepancy, but the factor matters when comparing
  `rho/L` with the regime ratio `m'/j = 8`.
* **Profiles.** Sylvester `H_{2^k}`: `S_u = {x : <u,x> = 1}`, `u != 0`. Paley-I matrix `K` on `(inf, 0..q-1)`
  with `K + K^T = 2I`. Hadamard profile: rows `2..` of `K` with its columns normalised by the first row. sdc: the
  record's skew double core. Restriction to an ordered point list `J`: `S -> S cap J`, trivial sets dropped, last
  element of `J` = base. Random law: `r ~ U{n-1..2n}` distinct canonical cuts, each point with probability 1/2,
  `L ~ U{1,2}`, `random.Random(seed)`. All 16 box sizes `|A|` that I rebuilt agree with the record and with
  `boxes.csv`. These are 1877 (M8'), 4861 (J9a), 4787 (J9b), 4997 (had12|J9 `L=3`), 4457 (sdc11|J9 `L=3`),
  4525 / 3847 (had12 / syl16 |J10 `L=2`), 3803 / 2233 / 351 / 363 (random `n=9` seeds 2/20/1/8), and 1255 (random
  `n=10` seed 2).
* **Least spaces.** On the cube every monomial agrees with a standard monomial (exponents `<= 2`) of no larger
  degree, so `M_e(S)` lives in the degree-`e` standard monomials. Mod `p`: one rref of the degree-ordered evaluation
  matrix. The rows whose pivot lies in block `e`, restricted to block `e`, span `M_e(S)`. `H_S(e)` = their number.
* **(SHIFT_e), reach.** `dim sum_{l} M_e(A cap (A-l)) = H_A(e)`, with one representative per `+-l`
  (`M_e(A_{-l}) = M_e(A_l + l) = M_e(A_l)`, checked exactly in §2). reach = first failing `e`, and it is `<= rho`.
  The containment `sum_l M_e(A_l) subset M_e(A)` is re-checked numerically at every degree of every run.
* **Primes.** `Q1 = 2^31 - 1 = 2147483647` and `Q2 = 2^30 + 3 = 1073741827`, both prime (sympy) and both
  **different** from the primary's `P1, P2`. So every mod-`p` value below comes from two primes independent of the
  primary's.
* **zeta.** `zeta_pair`: generators in `L_pair` (one per sign class), with multiplicity. `zeta_all`: all primitive
  `v in {-2..2}^d`. Both use exhaustive depth-first search. The state is the set `K` of admissible base points
  (`K <- K cap (K - v)` when `v` is added). Pruning uses an exact LP-dual bound on the widths of `K`. The global
  bound is an LP dual over coordinates, the base coordinate `sigma` and the slabs. It is solved by HiGHS, rationalised,
  and rescaled so that `min_v c_v = 1` exactly over **all** candidates, all in integer arithmetic. To show that no
  zonotope of order `t > zeta_pair` exists, only generators with `c_v - 1 <= B - t` can occur, with total excess
  `<= B - t`. Every witness is re-verified exactly: all `2^t` subset sums lie in `A`, every moment of order `< t`
  vanishes (integer arithmetic), and the order-`t` moment equals `t! prod <v_i,u>`, which is nonzero.
* **Exact over Q** (`exact_shift.py`, `exact_count.py`). Exact integer bases of `V_Q(S,e)` come from flint fmpz
  nullspaces, re-verified by exact multiplication. A mod-`p` rank of an *integer* moment matrix built from them is a
  lower bound for its `Q`-rank. Together with containment and `sum_e H_Q(e) = |A|`, this certifies `(SHIFT_e)` over
  `Q` (details in §6).

## 2. Engine validation and the round-27 reference values

**Validation against exact arithmetic** (`v0_exact_small.py` → `v0_exact_small_out.txt`). The test set is 16 small
slab boxes: 14 random boxes with `d = 3, 4`, the full cube `d = 3`, and syl8 `L = 1`. On each:
* `H_A`, and `dim sum_l M_e(A_l)` with **both signs** of every `L_pair` shift, were computed exactly over `Q`
  (fmpq rref);
* these were compared with my mod-`Q1` and mod-`Q2` engine at every degree;
* the check `M_e(A_{-l}) = M_e(A_l)` was made exactly.

There were **0 disagreements**, including non-spanning dimensions (6 vs 9, 5 vs 7, 4 vs 6) and the cube's
reach `2d = 6`.

**Round-27 reference values** (reach mod `Q1` / mod `Q2`; `zeta` exact):

| box | `|A|` | `H_A` | reach `L_pair` (fail dims) | reach cube (fail dims) | `zeta_pair` / `zeta_all` (LP) | primary / round 27 |
| --- | --- | --- | --- | --- | --- | --- |
| sdc7 `L = 2` | 757 | 1 7 28 77 161 251 188 43 1 | 6 / 6 (166/188) | not rerun | 7 / 7 (7) | 6, 7 ✓ |
| syl8 `L = 2` (M8) | 1037 | 1 7 28 77 161 259 287 175 42 | 6 / 6 (286/287); **exact over `Q`** (§6) | 7 / 7 (`e = 6`: 287 = 287 after 155 classes; `e = 7`: 141/175) | 7 / 7 (7) | 6, 7, 7 ✓ |
| syl8 `L = 3` (M8') | 1877 | 1 7 28 77 161 266 357 386 329 195 63 7 | 11 / 11 (none below `rho`); **exact over `Q`** (§6) | — | 10 / 10 (21/2) | 11, 10 ✓ |
| sdc7 `L = 3` | 1765 | 1 7 28 77 161 266 357 379 301 153 35 | 10 / 10 (= `rho`) | — | 10 / 10 (21/2) | 10, 10 ✓ (tie) |

**All three requested round-27 values reproduce**, with primes different from those of every earlier run.

**Witnesses.** The primary's M8' `L_pair` witness verifies exactly: `k0 = (-1,-1,-1,0,-1,0,1)`, generators
`e_0, e_1, e_2, e_2, e_3, e_4, e_0-e_6, e_1-e_3, e_4-e_5, e_5-e_6`, support 540, `|c|_inf = 2`. I also found an
independent `L_pair` witness, `k0 = (-1,-1,-1,0,-1,0,0)`, generators `e_0, e_0, e_1, e_1, e_2, e_3, e_4, e_5,
e_2-e_3, e_4-e_5`, verified exactly (support 324, `|c|_inf = 4`). So **`zeta_pair(M8') = 10` is confirmed**, and
the killer's non-pair generator was not needed. On M8' the pair LP bound is already 21/2, so `zeta_pair <= 10`
needs no search.

## 3. Every box with reach > zeta (both of the primary's)

| box | definition | `|A|` | `H_A` | reach `L_pair` | `zeta_pair` / `zeta_all` | label |
| --- | --- | --- | --- | --- | --- | --- |
| M8' | Sylvester `H_8`, `n = 8`, base 7, 7 cuts, `L = 3` | 1877 | 1 7 28 77 161 266 357 386 329 195 63 7 | **11** (= `rho`; `Q1`, `Q2`; **exact over `Q`**) | **10 / 10** | reach exact; zeta exact |
| syl16\|J9a `L = 3` | Sylvester `H_16` (`x in F_2^4` by binary index), restricted to `J = (0,1,...,8)`, base 8, 15 cuts | 4861 | 1 8 36 112 266 504 784 1001 1029 756 315 49 | **11** (= `rho`; `Q1`, `Q2`; per-degree sums 8, 36, ..., 315 = `H_A(e)` for every `e <= 10`; at `e = 10` the sum spans after 33 of 36 classes mod `Q2`) | **10 / 10** | reach mod `p` (4 primes in all); zeta exact |

Both reach > zeta claims are **confirmed**: reach − zeta = +1 on both. Since reach = `rho` on both, Thm 3.1 gives
ZE-braid at the top order, and truth = 11 (exactly for M8'; mod `p` for J9a).

**Is J9a "M8''s own configuration plus one point"?** As a point configuration, yes. As a box, only partly. The
check is `t5_slices.py` → `t5_slices_out.txt` (exact set comparison):
* `J9a = {0..7} u {8}`, and `{0..7}` is the affine hyperplane `{x_3 = 0}`. It is the only affine hyperplane in
  `J9a`. J9b and the three random 9-subsets contain none, which confirms the primary's statement.
* But the extra point 8 is the **base**. So the face `k_8 = 0` of the J9a box is **not** the M8' box. It is
  `M8' cap {|k_base| <= 1}`, 1093 of M8''s 1877 points, because the M8' base coordinate is unconstrained in M8'
  (it ranges over `[-4, 4]`) but is a cube coordinate in J9a.
* The box that **does** contain M8' as a coordinate face is J9a': the same 9 points with the base moved into the
  hyperplane, `J = (1,2,...,8,0)`, base 0, `|A| = 4791`. Its face `k_8 = 0` equals the M8' box exactly (1877
  points). This is the natural "one-point extension of M8'". The primary did not compute it; I did (§4).
* *Discrepancy (descriptive, minor):* `compute.md` calls J9a "M8''s own point configuration plus one point". That is
  true of the 9-point set, but the box is not an extension of the M8' box. Both base choices are needed before the
  `d = 8` success can be attributed to M8''s box.

## 4. Further boxes (my choice: the largest computed ones, faithful random boxes, `zeta_pair < zeta_all` cases, and one new adversarial box)

Box specs as in `boxes.csv`: `rsyl:4:J:L` = Sylvester-16 restricted to the ordered list `J` (last = base);
`rhad:11:J:L` = Paley-I Hadamard 12 on `(inf,0,...,10)` (index 0 = inf); `rand:n:seed` = killer's law. Reach mod
`Q1` / `Q2`; zeta exact (§5). The table is generated by `summarize.py` and saved as `verify_checks/summary_table.md`.
`agree` compares `|A|`, `H_A`, reach, failure dimensions and both zetas with `boxes.csv`.

| box | d | L | `|A|` | `rho` | reach (fail dims) | `zeta_pair` / `zeta_all` | reach − `zeta_all` | primary | agree |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| had12\|J9 (`rhad:11:0..8:3`), largest `d = 8` | 8 | 3 | 4997 | 11 | 10 / 10 (342/368) | 11 / 11 | −1 | 10; 11/11; 342/368 | yes |
| syl16\|J9b (`rsyl:4:0,1,2,4,8,15,3,5,6:3`) | 8 | 3 | 4787 | 11 | 10 / 10 (215/281) | 10 / 10 | 0 | 10; 10/10; 215/281 | yes |
| had12\|J10 (`rhad:11:0..9:2`), largest `d = 9` | 9 | 2 | 4525 | 8 | 7 / 7 (1112/1227) | 7 / 7 | 0 | 7; 7/7; 1112/1227 | yes |
| syl16\|J10 (`rsyl:4:0..9:2`) | 9 | 2 | 3847 | 8 | 6 / 6 (1250/1251) | 7 / 7 | −1 | 6; 7/7; 1250/1251 | yes |
| rand `n=9` seed 2 (largest random) | 8 | 2 | 3803 | 11 | 8 / 8 (722/723) | 11 / 11 | −3 | 8; 11/11; 722/723 | yes |
| rand `n=9` seed 20 (faithful) | 8 | 2 | 2233 | 9 | 6 / 6 (655/657) | 8 / 8 | −2 | 6; 8/8; 655/657 | yes |
| rand `n=9` seed 8 (`zeta_pair < zeta_all`) | 8 | 1 | 363 | 6 | 4 / 4 (140/143) | 5 / 6 | −2 | 4; 6/5; 140/143 | yes |
| rand `n=10` seed 25 (`zeta_pair < zeta_all`) | 9 | 1 | 269 | 5 | 4 / 4 (25/89) | 4 / 5 | −1 (tie vs `zeta_pair`) | 4; 5/4; 25/89 | yes |
| rand `n=10` seed 2 (largest random `d = 9`) | 9 | 1 | 1255 | 7 | 5 / 5 (413/415) | 7 / 7 | −2 | 5; 7/7; 413/415 | yes |
| **J9a' (`rsyl:4:1,2,3,4,5,6,7,8,0:3`), new** | 8 | 3 | 4791 | 11 | **9 / 9 (740/741)** | 10 / 10 | **−1** | not in primary set | — |
| the 27 round-27 random boxes (`n = 6` seeds 100–114, `n = 7` seeds 1–12) | 5, 6 | 1, 2 | 39–603 | — | equal to the CSV on all 27 | (CSV) | 20 `<`, 7 `=` | 20/7/0 | yes |

J9a' details:
* `H_A = 1 8 36 112 266 504 784 1002 1023 741 286 28`.
* `zeta = 10` exactly (LP 43/4).
* Its face `k_8 = 0` is exactly the M8' box (§3), so by the subset criterion (theory Thm 3.1 / §2.4)
  `reach_chain(J9a') >= 11`. Only the constant-chain reach falls below zeta.
* The one-dimensional defect at `e = 9` is **not** of the implied-slab type of theory Thm 2.3. No ternary `phi` has
  range 10 on `A` and at most 9 on every `A_l` (`t6_J9ap_out.txt`; `sigma` ranges over `[-4, 4]`).
* So this failure is not a design/BCO effect in the theory lane's sense. It is one of the "non-univariate residues"
  that theory §3.4 leaves unexplained.
* Cube shifts do not rescue it. With the 150 largest classes of `{-1,0,1}^8` (all 36 pair classes included), the
  sum at `e = 9` is still 740/741 mod `Q1` and mod `Q2` (`logs/J9ap_cube_Q*.log`). That is a lower bound for the
  full cube set, the same truncation the primary uses.

**Not rerun by me** (accepted from the primary, with the mod-`p` label it gives them): the other 64 boxes of
`boxes.csv`, including sdc11|J9/J10, the remaining `d = 9` random boxes and the cube-shift reaches.

## 5. zeta certificates and the "no larger zonotope" claims

All values are from my own LP and DFS (`zono_v.py` → `zeta.jsonl`, logs in `verify_checks/logs/zeta_*.log`).
* "LP" is my exact rational LP-dual bound, over all primitive candidates (`zeta_all`) or over `L_pair`
  (`zeta_pair`).
* "search" is the exhaustive DFS that rules out order `floor(LP)` when it is not attained.
* Every witness passed the exact check of §1.

| box | `zeta_pair` (LP, DFS nodes) | `zeta_all` (LP; candidates; search for `zeta+1`) | primary | agree |
| --- | --- | --- | --- | --- |
| M8' | 10 (21/2) | 10 (21/2; 34 399) | 10 / 10, LP 21/2 | ✓ |
| syl8 `L=2` | 7 (7) | 7 (7) | 7 / 7 | ✓ |
| sdc7 `L=2` | 7 (7) | 7 (7) | 7 / 7 | ✓ |
| sdc7 `L=3` | 10 (21/2) | 10 (21/2) | 10 | ✓ |
| **syl16\|J9a `L=3`** | **10** (11; exhaustive, 873 100 nodes) | **10** (LP **11** exactly; 155 263 candidates; the dual leaves exactly the 36 `L_pair` classes with `c_v = 1`; exhaustive search for order 11, 873 085 nodes: none) | 10 / 10, LP 11, 36 classes, 3 506 942 nodes, MILP infeasible | ✓ |
| syl16\|J9b `L=3` | 10 (11; 33 699 nodes) | 10 (LP 11; my dual leaves 59 generators; search for 11: none, 35 206 nodes) | 10, LP 11 | ✓ |
| syl16\|`3,11,2,8,9,13,15,10,7` `L=3` | 10 (89/8) | 10 (89/8; 24 generators; none of order 11) | 10, LP 89/8 | ✓ |
| had12\|J9 `L=3` | 11 (11) | 11 (11) | 11 = `rho` | ✓ |
| syl16\|J10 `L=2` | 7 (15/2) | 7 (15/2) | 7 (15/2) | ✓ |
| had12\|J10 `L=2` | 7 (22/3) | 7 (22/3) | 7 (22/3) | ✓ |
| rand `n=9` s2 | 11 (11) | 11 (11) | 11 (11) | ✓ |
| rand `n=9` s20 | 8 (9) | 8 (LP 9; 37 generators; none of order 9) | 8 (9) | ✓ |
| rand `n=9` s8 | 5 (5) | 6 (6; order-6 witness uses a non-pair generator) | 6 / 5 | ✓ |
| rand `n=10` s25 | 4 (41/8) | 5 (16/3; witness found) | 5 / 4 | ✓ |
| rand `n=10` s2 | 7 (7) | 7 (15/2) | 7 (15/2) | ✓ |
| *new:* syl16\|`1,...,8,0` `L=3` (J9a', base inside the hyperplane) | 10 (43/4) | 10 (43/4) | not in the primary's set | — |

**Certificate checks requested by the task.**
* *Vertices in `A`:* every one of the `2^t` subset sums of every witness lies in `A`, checked as a condition of the
  search and again in the exact verifier.
* *Generators in `L_pair`:* every `zeta_pair` witness is in `L_pair` by construction. The primary's M8' witness
  is all-`L_pair` (re-checked, §2).
* *Order:* checked exactly (vanishing moments below `t`, nonzero order-`t` moment).
* *Non-existence of a larger zonotope:* where `floor(LP) = zeta` the LP dual is the proof. I checked it in exact
  integer arithmetic over every candidate: all primitive `v in {-2..2}^d` that fit the widths, which is the complete
  generator class. Non-primitive generators never help, because `m v'` can be replaced by `m` copies of `v'` in a
  slab box. Where `floor(LP) = zeta + 1` (J9a, J9b, the random Sylvester subset, random `n=9` s20) the exhaustive
  search over the dual-reduced generator set is the proof. The reduction is valid because `sum_i (c_{v_i} - 1) <= B - t`
  for any order-`t` zonotope, with `c_v >= 1` exact.

So **`zeta(J9a) = zeta(M8') = 10` is confirmed exactly, for both generator classes**, by independent code. One
point on the definition: the zeta of `boxes.csv` is `zeta_all`, and the task's definition is `zeta_pair`. On every
box where reach > zeta occurs they coincide, so the reading does not depend on the choice.

## 6. Exact-over-Q checks

**Method.**
* `dim V_Q(S,e)` is the exact nullity of `E_{<e}(S)^T`, from flint `fmpz_mat.nullspace`. The returned integer basis
  is re-verified by exact multiplication.
* Lower bound for `dim_Q sum_l M_e(A_l)`: take the exact integer kernel vectors of `E_{<e}(U)` for subsets
  `U subset A_l` and form their integer degree-`e` moment vectors. Their mod-`p` rank is at most their `Q`-rank.
* Upper bound: `sum_l M_e(A_l) subset M_e(A)`, which is a theorem.
* `exact_shift.py` uses exact `H_Q(e) = dim V_Q(A,e) - dim V_Q(A,e+1)` degree by degree.
* `exact_count.py` uses a counting argument. If the exact lower bounds `s_e` equal `H_p(e)` for `1 <= e < rho`,
  and `dim V_Q(A,rho) = H_p(rho)`, `dim V_Q(A,rho+1) = 0` and `1 + sum s_e + H_p(rho) = |A|`, then
  `sum_e H_Q(e) = |A|` together with `H_Q(e) >= s_e` forces `H_Q = H_p` and `(SHIFT_e)` over `Q` for every `e < rho`.

| statement | result | script / log |
| --- | --- | --- |
| small boxes (16), all degrees, both signs | exact = mod `Q1` = mod `Q2`, 0 disagreements | `v0_exact_small_out.txt` |
| **M8': reach(`L_pair`) = 11 over `Q`** | **certified** (`exact_count`: `s_e = H_p(e)` for `e = 1..10`; `dim V_Q(A,11) = 7`, `dim V_Q(A,12) = 0`; sum = 1877). Independently, `exact_shift` (e = 6..11): exact `H_Q = 357 386 329 195 63 7`, `(SHIFT_e)` certified for `e = 6..10`, and at `e = 11` every `V_Q(A_l, 11)` is exactly 0, so `(SHIFT_11)` fails exactly | `logs/exactc_M8p.log`, `logs/exact_M8p.log` |
| M8: `(SHIFT_e)` for `e <= 5`, exact `H_Q(1..6)` = 7 28 77 161 259 287, `dim_Q sum_l M_6(A_l) >= 286` | certified | `logs/exact_M8.log` |
| M8: `sigma^6 notin J(A)_6`, `sigma^6 in J(A_l)_6` (56 signed shifts), `sigma^7 in J(A)_7` | exact | `t1_sigma_out.txt` |
| hence **reach(M8, `L_pair`) = 6 over `Q`, with an exactly one-dimensional defect spanned by `sigma^6`** | exact | (combination of the two rows above) |

Consequences:
* the round-27 killer's M8' value is a theorem, not just mod-`p` evidence;
* so is the theory lane's Cor 2.2 (`(M8')^t`: reach `11t`, `zeta 10t`);
* the "mod-`p` caveat" on M8' in `compute.md` §7 and `theory.md` §0 can be dropped.

**Not done exactly: J9a (`d = 8`, 4861 points).** The same certificate needs exact nullspaces of matrices up to
about 4850 x 4300. Extrapolating from M8' (218–247 s for all degrees), that would take many hours of fraction-free
elimination, so I stopped it. J9a's reach 11 is therefore **mod-`p` evidence at four independent primes**: the
primary's `P1, P2` and my `Q1, Q2`, with identical `H_A` and identical per-degree dimensions.

## 7. Theory lane (`theory.md`), line by line

Verdict: **every proof I checked is correct up to minor gaps (listed). Every numerical claim I re-ran reproduces.**
One correction upgrades a claim: Cor 2.2 no longer depends on mod-`p` data, because reach(M8') = 11 is now exact
over `Q` (§6).

**§1 (definitions, elementary inequalities).** Correct. The convention (`L` bounds `|k(S)|`) and `zeta_L` with
generators in `L` both match. `reach <= truth_L <= rho` and `zeta_L <= truth_L` are correct. The reference table
reproduces: `H_A` for syl8 and sdc7 at `L = 2, 3`, defects 1 (M8, `e = 6`), 22 (sdc7 `L = 2`, `e = 6`) and 7 (M8',
`e = 11`), and `zeta_pair` = 7, 10, 7, 10 for syl8 `L = 2, 3` and sdc7 `L = 2, 3` (my exhaustive DFS). I did not
rerun `L = 1`. The linkage identity `H_A(e) = T_d(e) - H_B(2d-e)` holds for syl8 `L = 2, 3` and sdc7 `L = 2` at
every `e` (`t4`, mod `Q1`).

**Thm 2.1 (products).** Proof checked:
* (i) `gr` of a tensor product of filtered spaces.
* (ii) `A_l = (A_1)_{l_1} x (A_2)_{l_2}`. The in-block and cross-root contributions are as stated, the quotient
  identity `(M (x) M')/(N (x) M' + M (x) N') = (M/N) (x) (M'/N')` holds, and the first failure is at `r_1 + r_2`.
* (iii) Projection. *Minor gap:* for `zeta_all` the block-1 projection of a primitive generator can be
  non-primitive (e.g. `(2, 1) -> 2`). This is harmless: in a slab box (`A = conv(A) cap Z^d`) a generator `m v'` can
  be replaced by `m` copies of `v'`, which does not lower the order. The proof should say so.
* Numerically, on 4 products of small boxes the Hilbert functions convolve, reach adds (both primes) and
  `zeta_pair` adds (`t3_product_out.txt`).

**Cor 2.2.** Correct. Its one mod-`p` input, reach(M8') = 11, is now **certified exactly over `Q`**
(`exact_count.py`, §6). `zeta_pair(M8') = zeta_all(M8') = 10` exactly (mine: DFS plus LP bound 21/2). So
`(M8')^t` has `reach = 11t`, `zeta = 10t` **unconditionally**, and `reach <= zeta + c` (and
`reach <= (1+o(1)) zeta`) fails on slab boxes. I agree with the theory lane that this is vacuous for the route:
block profiles violate the character gap.

**Thm 2.3 (implied-slab obstruction).** Correct: `omega_l(phi) in I(S_l)` has leading form `phi^{R_l}`, `J(S_l)` is
an ideal, and the shift-by-shift discussion and the zonotope certificate for `(N_phi)` (`t! prod phi(v_i)`) are
correct.

**Cor 2.4 (BCO).** Correct for slab boxes. *Minor:* "`R - 1 <= 2h`" is stated for any `S` in the cube but needs
`min sigma >= -h`, e.g. `S = -S` (all slab boxes are symmetric). For non-symmetric `S` only `R - 1 = h - min sigma`
holds.

**Prop 2.5.** Correct as stated (balanced, `f = 1/2`). For Sylvester, `f = n/(2(n-1)) > 1/2`, which only
strengthens it. Whether regime profiles are balanced at the base in this exact sense is a hypothesis, as stated.

**§2.3 (rigidity data).** Reproduced by exact enumeration:
* syl8 `L = 1, 2, 3`: `h = 1, 3, 4`; (R) no / yes / no. sdc7 `L = 1, 2, 3`: (R) fails, with `h = 4` at `L = 3`.
* Sylvester `n = 16` (`3^15` points, `t2_syl16_out.txt`): `L = 2, 3, 4` give `h = 3, 4, 7`, (R) yes / no / yes,
  max `sigma` on `A_{e_0-e_1}` = 2, 4, 6, and the top levels have 35, 420, 15 points. This is exactly the theory
  table's `n = 16` row.
* The killer's 27 random boxes (`t4`): (R) holds on exactly 3 (`n=6` seeds 105, 106; `n=7` seed 9), and on each
  reach `= 2h = 4` (both primes), as claimed.
* I did not rerun the IP rows `n >= 32` or the random-profile IPs (`t16`).

**Thm 3.1 (chains).** Proof correct; it is round-27 Lemma 3.0 applied to `S_{e+1}`. The zonotope chain satisfies
`(CH_e)` because `S_e subset (S_{e+1})_{v_{e+1}}`. So `reach_chain >= max(zeta_L, reach(A'))` holds. On M8,
`A cap {|sigma| <= 2}` (1023 points) has reach 7 with defect 141/169 at `e = 7` (both primes, `t4`), as claimed.

**Lemma 3.2.** Correct: `I(S cap S') = I(S) + I(S')` for finite point sets (`K^S (x)_R K^{S'} = K^{S cap S'}`), and
the leading form of `tau_l f - f` is `d_l F`.

**Prop 3.3 (linkage form).** Labelled "proved given the linkage identity". I checked the identity numerically on
three boxes but **did not re-derive** the duality step (the inclusion `(R_e/J(C_l)_e)^* -> (R_e/J(C)_e)^*` as
multiplication by `mu_l`), so it is not independently verified. It is used nowhere else.

**§3.4, §3.5.**
* M8': `h = 4`, 28 top points, 252 legal pair moves (exact enumeration, `t1_sigma_out.txt`).
* sdc7 `L = 3`: `h = 4`, 126 legal moves, not rigid.
* Not rerun: the power-span counts (6 of 8 defects) and the staircase result (`e <= 7` on M8').

**Thm 4.1 (the M8 obstruction is `sigma^6`).** Re-verified **exactly over `Q`** (fmpz ranks, `t1_sigma.py`):
* `sigma(A) = [-3, 3]`.
* The top level consists of 7 0/1 points whose supports are exactly the Fano lines (complements of the 7 patterns).
* `sigma^7 in J(A)_7` and `sigma^6 notin J(A)_6`.
* `sigma^6 in J(A_l)_6` for all 56 **signed** `L_pair` shifts.
* The defect at `e = 6` is 1 mod `Q1` and `Q2` (286 vs 287), and `dim_Q sum_l M_6(A_l) >= 286` with
  `H_Q(6) = 287` is certified exactly (§6). Together with the exact non-membership, the defect is **exactly**
  one-dimensional and spanned by `sigma^6` (item 4 is upgraded from mod `p` to exact).
* The proof of item 3 (range 7 to 6 on every `A_l`) is correct.

**§5.2.** Item 4 ("the plain count `H_A(e) <= sum_l H_{A_l}(e)` always holds in the regime") is asserted, not
proved. It is plausible (`~d^2` shifts), but it should be labelled heuristic. Items 1–3 follow from the results above.

**§5.3, §6.** These are judgements, and I agree with them with one caveat. The advice to read failures at
Sylvester `L = 2^m` as design effects cuts both ways: by the same logic the two successes (M8', J9a; vertex cuts
`L = |P| - 1 = 3`, no design) are also artefacts of small symmetric profiles (theory §3.5 says so for M8'). Removing
both kinds of artefact leaves the random and non-design boxes, where reach is never above `zeta` (0 of 77 random
boxes).

## 8. Discrepancies and caveats (complete list)

**Numerical discrepancies with the primary lane: none.** Every `|A|`, `H_A`, reach, failure dimension, `zeta_pair`,
`zeta_all` and LP bound that I recomputed agrees exactly with `compute.md` / `boxes.csv`. On the 27 round-27 random
boxes my reach (two new primes) equals the CSV on all 27. The CSV's own counts are consistent with the text: 2 of 105
boxes with reach > zeta, P1 and P2 agreeing everywhere, every zeta flagged exact, 69 of 77 random boxes below zeta.

**Numerical discrepancies with the theory lane: none** among the claims I re-ran (§7).

**The primary lane's own proofs** (`compute.md` §§1–2), checked:
* *Row-echelon extraction:* one rref per set gives every `M_e(S)`. Correct; my engine uses the same fact,
  independently implemented, and it is validated exactly in §2.
* *`reach <= rho`:* correct.
* *Refinement lemma* (`H_A(rho) = 1` and `zeta_all < rho` imply `truth <= rho - 1`): correct, given the record's
  classification (interp Thm 2.2) and Ostrowski. The zonotope `a + sum_{zeta_j = 1} m_j [0, v_j]` is a Minkowski
  summand of `Newt(P_c)` containing its vertex `a`, so it lies in `conv(A)`, and its lattice vertices lie in `A`.
  Non-primitive `v_j` only increase the order after splitting. It applies to random `n = 9` seed 20, where I confirm
  `H_A(9) = 1` and `zeta_all = 8`.
* *Translation invariance* `M_e(A_{-l}) = M_e(A_l)`: correct, and checked exactly in §2.

**Corrections, upgrades and caveats:**
1. *(upgrade)* reach(M8') = 11 and reach(M8) = 6 are **exact over `Q`** (§6), not only mod `p`. Theory Cor 2.2 is
   therefore unconditional. The phrase "mod p" attached to these values in `compute.md` §§0, 7 and `theory.md` §§0,
   2.1 can be removed.
2. *(definition)* The task defines `zeta` with generators in `L_pair`. `compute.md` compares against `zeta_all`.
   On the two reach > zeta boxes the two coincide (10), so no conclusion changes. Against `zeta_pair` the random
   boxes give two more ties and still no reach > zeta (CSV).
3. *(coverage gap, material for the rule)* reach > zeta was seen only on Sylvester restrictions at `L = 3`. The
   `d = 9` data contain **no `L = 3` box at all** (out of memory reach). So the `d = 9` entry "0/24" says nothing
   about the class in which the phenomenon lives. The data for that class are: `d = 7`, 1 of 1 Sylvester `L = 3`
   (M8'); `d = 8`, 1 of 5 Sylvester `L = 3` restrictions (J9a). Here J9a' is not counted: it is not in the primary's
   set (§4).
4. *(denominator)* The headline shares (1/6, 1/48, 0/24) are dominated by random boxes, 75 of 77 of which are
   zonotopal (`zeta = rho`), where reach > zeta is impossible a priori. The faithful-box shares (1/3, 1/10, 0/2) rest
   on 3, 10 and 2 boxes. The direction (no growth) is supported; its strength is modest.
5. *(confounding)* "reach − zeta falls as `rho/L` grows" (`compute.md` Table 6.4) compares different families.
   * Every one of the 53 boxes with `rho/L > 4` is a **random** box, and 51 of these 53 are zonotopal.
   * Every structured box has `rho/L <= 4`.
   * So the decline is mainly random versus structured (and mostly `L = 1`), not a `rho/L` effect within a family.
   * The trend is in the data, but it does not show that the gap widens as boxes grow towards the regime.
6. *(theory, minor)* Thm 2.1(iii): non-primitive projected generators are not mentioned (harmless). Cor 2.4: needs
   `min sigma >= -h` (true for slab boxes). §5.2 item 4 is a heuristic, not a proof. Prop 3.3's duality step was not
   independently re-derived.
7. *(process, mine)* While reorganising my job queue I briefly had more than two heavy processes for about a minute,
   because a `pkill -f` pattern matched my own shell. The extra jobs were killed and rerun. No result is affected.

## Appendix: scripts and outputs (`verify_checks/`, all own code)

| file | content | label |
| --- | --- | --- |
| `vlib.py` | profiles (Sylvester, Paley-Hadamard, sdc, restriction, random law), boxes, standard monomials, mod-`p` least spaces via one rref | — |
| `run_reach.py` → `reach_2147483647.jsonl`, `reach_1073741827.jsonl`, `logs/*_pair_.log` | reach w.r.t. `L_pair` or cube shifts, lazy per-degree sums, containment re-check | mod `p` |
| `zono_v.py` → `zeta.jsonl`, `logs/zeta_*.log` | `zeta_pair`, `zeta_all`: exact LP dual, exhaustive DFS, exact witness check | exact |
| `v0_exact_small.py` → `v0_exact_small_out.txt` | engine vs exact `Q` on 16 small boxes, both signs | exact |
| `exact_shift.py` → `exact_shift.jsonl`, `logs/exact_M8p.log`, `logs/exact_M8.log` | exact `H_Q`, exact lower bounds for shift sums | exact |
| `exact_count.py` → `exact_count.jsonl`, `logs/exactc_M8p.log` | counting certificate reach_Q(M8') = 11 | exact |
| `exact_count2.py` | faster variant (moment product mod `p`, Freivalds kernel check); re-certified M8' in 218 s | exact (kernel check probabilistic) |
| `t1_sigma.py` → `t1_sigma_out.txt` | `sigma` ranges, rigidity (R), moves; exact `sigma^6`/`sigma^7` memberships on M8 | exact |
| `t2_syl16.py` → `t2_syl16_out.txt` | Sylvester `n = 16`, `L = 2,3,4`: `h`, (R), max `sigma` on `A_{e_0-e_1}` | exact enumeration |
| `t3_product.py` → `t3_product_out.txt` | Thm 2.1 on 4 products | mod `p` + exact DFS |
| `t4_theory_misc.py` → `t4_theory_misc_out.txt` | 27 random boxes: (R), reach; linkage identity; M8 minus `|sigma| = 3` | mod `p` / exact enumeration |
| `t5_slices.py` → `t5_slices_out.txt` | J9a and J9a' faces vs the M8' box | exact |
| `t6_J9ap_obstruction.py` → `t6_J9ap_out.txt` | implied-slab (Thm 2.3 range) search on J9a' | exact enumeration |
| `summarize.py` → `summary_table.md` | my results vs `boxes.csv` | — |
| `notes.md` | working notes | — |
| `runq.sh`, `chain*.sh`, `jobs*.txt`, `zjobs.sh` | job queues (at most two heavy jobs, `prlimit --as=3500000000`) | — |
