# Round 13 — Locality Loc(b0) of the leading-form ideal J(A): exact computations

Status: complete. Not completed: Paley `L_col(7)` and Syl16 degree `>= 8`. Every claim carries its status (proved /
exact / mod `p` evidence / heuristic). Scripts are in `loccompute_checks/`; working notes are in `loccompute_notes.md`.

## 0. Summary of results

**Question (README round 12, §7.1).** `Loc(b_0)` asks whether, for some `b_0` independent of `n`, `J(A)_i` is the sum
of the pullbacks `pi_F^* J(pi_F A)` over families `F` of at most `b_0` integer directions, in all degrees `i <= m`.

**Answer on exact models.** `Loc(b_0)` holds with small `b_0` in most models, but it is **false at small `b_0` in an
exact model** that is zonotope-free, overdetermined and indecomposable (Paley-12, which is LP-subcritical). It is
therefore not a structural property of slab boxes.

| model (`T = 1`) | `m` | `m/(2L)` | `b_0` (columns) | `b_0` (columns + coordinates) | status |
| --- | --- | --- | --- | --- | --- |
| `Syl(8,1)` | 3-4 | 1.5-2 | 1 | 1 | exact |
| `M8 = Syl(8,2)` | 8 | 2 | **3** (exactly the 7 lines) | **3** | exact, both directions (`b = 2` fails over `Q`) |
| `M8' = Syl(8,3)` | 11 | 1.83 | **1**: `J(A) = J_w` in all degrees `<= 11` | 1 | exact |
| two-cut blocks (400 of 402 with `d <= 7`) | `<= rho + 1` | up to 4.5 | 1 | 1 | exact |
| three-cut blocks, `d <= 6` (953 types) | `<= rho` | up to 5 | 1 (827), 3 (27), none (99) | `<= 4` always | equalities exact; failures mod `p` |
| `Syl(16,1)` | 3 | 1.5 | — | 3 | exact (`b = 2` failure mod `p`) |
| **Paley-12, `L = 2`** | 8 | 2 | **`> 5` over `Q`; `> 6` mod `p`** | `> 4` mod `p` | exact (Proposition 4.4) |
| **`Syl(16,2)`** | `<= 7` | `<= 1.75` | **`>= 3` over `Q`** (lines give exact equality through degree 6); **`>= 7` mod `p`** at `m = 7`: whole 7-column planes are needed (`L_col(6)_7 = 6226`, `+ planes = 9691`, `J(A)_7 = 9769`) | `>= 3` over `Q` | see §3.6 |

**Main results.**
1. **Proposition 4.4 (exact; the key negative outcome).** In Paley-12 (`T = 1`, `L = 2`, `m = 8`), over `Q`,
   `J(A)_7 ⊄ L_col(5)_7`. In the `Z_11`-invariant part, `dim_Q P_0(J(A)_7) = 137`, against `82` for `L_col(5)` and `72`
   for `J_w`.
   * The certificate uses 55 exact invariant polynomials of degree 7, vanishing at all 20153 points.
   * Mod `p`, the deficit is 725 (`J_w`), 615 (`<= 5` or `<= 6` columns) and 725 (`<= 4` columns or coordinates).
   * The non-local generators are **H-basis completions** of the slab system (Lemma 4.3): leading forms from non-Koszul
     syzygies of the width powers `y_p^5` in `R = S/(x^3)`. There are 440 such syzygies in degree 9 (mod `p`), where the
     Fröberg-type count predicts none.
   * They couple at least 3 slabs on 9-11 of the 11 coordinates. Projections onto `<= 4` columns are grids
     (Lemma 4.1), and pair cut boxes are exactly width-generated.
2. **README §7.1's family restriction `rho(pi_F A) < ceil(4w)` is wrong as a criterion.** In Paley every family that
   contributes below `m` has `rho(pi_F A) >= m`, which is exact by counting. So the route to
   `b_0 <= m'/(2j) + O(1)` via projection costs (README §7.3(3)) fails at its first step. The Sylvester data (item 5)
   point the same way: `b_0 >= 7` at `m/(2L) = 1.75` (mod `p`).
3. **Two mechanisms.**
   * (a) Arithmetic congruences among columns (Lemma 4.2) make projections non-grid. These are projection-local, with
     `b` = support of the congruence: Sylvester lines `b = 3` (mod 2), Paley `b = 5` (mod 3).
   * (b) Cut-box completions need coordinate directions, or are global as in Paley.
4. **No dependence on `m/(2L)`.**
   * Two-cut blocks are width-local up to ratio 4.5.
   * Three-cut blocks have `b_0(mixed) <= 4` at every ratio up to 5, with only a mild drift from 3 to 4.
   * Among the Hadamard models at ratio about 2, `b_0` ranges from 1 to `> 6`.
5. **Sylvester profiles: `b_0` grows with the subspace dimension.**
   * In `Syl(16,2)`, the first Sylvester model whose line images are full as in the regime, the 35 line families are
     **needed and exactly sufficient** through degree `6 = 3L`. This is exact: an exhaustive scan shows that no width
     power or 2-dimensional family supplies them. Unlike `M8'`, nothing is absorbed.
   * In degree 7 the 15 plane families (7 columns each) add `3465` generators beyond the lines (mod `p`).
   * No proper subfamily of a plane, and no family of `<= 6` columns outside a plane, produces any of them
     (`L_col(6)_7 = 6226`). So column `b_0 >= 7` already at `m = 7 = 1.75 * (2L)`.
   * The abstract large-coset images agree (mod `p`): lines have non-width generators exactly at `3j`, and planes have
     generators that their lines do not produce.
   * A dimension-count heuristic puts subspace families up to 15 columns below `m' ~ 8j`. Evidence therefore points to
     a `b_0` for Sylvester regime boxes that is **bounded in `n` but at least 7, plausibly 15**: not the hoped-for
     `m'/(2j) + O(1)`.

**What this does not do.** It proves nothing about growth or uniformity. `M <= (2/3+eps) log R/loglog R` stands.
`Loc` is neither proved nor refuted for any regime box. All exact models are LP-subcritical, and the Hadamard ones are
ball-bound (§6.1, §6.3).

## 1. Definitions, conventions, family sets

**Boxes.** All models are `T = 1` slab boxes `A = {k in {-1,0,1}^d : |<c_p, k>| <= L_p for every column p}`. Here `d = n - 1`, and the
coordinates are the non-base points; the base coordinate `k_base = -sum k_x` is constrained only through the columns.
* `Syl(2^k, L)`: points `F_2^k`, base `0`, columns `k(S_a)` with `S_a = {x : a.x = 1}` (`a != 0`), `|k(S_a)| <= L`.
  This is zetc §4's model (`hat k(a) = -2 k(S_a)`). `M8 = Syl(8,2)` and `M8' = Syl(8,3)`.
* Paley-12: coordinates `Z_11`, base `infinity`, columns `i + D` with `D = {0,1,3,4,5,9}`, `|k(i+D)| <= L` (zetc §5.5).
* Cut blocks: columns are indicator vectors of the sets `S_j`.

In this notation the width of a column is `2L` and its width power has degree `2L + 1`. The regime analogue is `2j` and
`2j + 1`, with `m' = ceil(4w) ~ 8j`, so the regime ratio is `m'/(2j) ~ 4`. The ratio quoted for each model below is `m/(2L)`.

**The reduced ring.** `k_x^3 - k_x` vanishes on the cube, so `x^3 in J(A)` for every `A`. Every local ideal below contains
the coordinate width powers `x^3` (`b = 1` families). All ideals are therefore compared in
`R = Q[x_1..x_d]/(x_1^3, ..., x_d^3)`, which has monomial basis `{x^e : e in {0,1,2}^d}` and `N_i = [t^i](1+t+t^2)^d`.
In `R`, `dim J(A)_i = N_i - H_A(i)`. Here `H_A` is the Hilbert function of `A`, computed from reduced monomials evaluated
on `A`.

**Local ideals.** For a family `F = (v_1, ..., v_b)` of integer directions, write `pi_F(k) = (<v_1,k>, ..., <v_b,k>)`.
The family contributes `pi_F^* J(pi_F(A))`. Pulled-back vanishing polynomials vanish on `A`, so this lies in `J(A)` over
`Q`. Two remarks:
* (span invariance) `J(T P) = (T^{-1})^* J(P)` for `T in GL_b(Q)`, so the pullback depends only on the `Q`-span of `F`.
  Hence "`b`" below always means the **dimension of the span**.
* (`b = 1`) For a single direction `v`, `J(pi_v A) = (y^s)` with `s = |<v, A>|`. The `b = 1` part is therefore the width
  ideal `J_w` generated by `<v,x>^{|<v,A>|}`.

The family sets used are:

| name | definition |
| --- | --- |
| `J_w(R)` | `<v,x>^{|<v,A>|}` for all primitive `v in [-R,R]^d` with `|<v,A>| <= m` (`R = 2` for `d <= 7`; `R = 1` for `d >= 8`, stated per model) |
| `L_col(b)` | `J_w(R)` plus `pi_F^* J(pi_F A)` for **all** families `F` of `<= b` profile columns |
| `L_mixed(b)` | `J_w(R)` plus pullbacks for all families of `<= b` directions drawn from columns and coordinate vectors (span dimension `b`) |
| `L_col(b) + extra` | as stated per run (e.g. column pairs together with line half-sums, `§3.1`) |

**`Loc` on a model.** `b_0(model, m)` is the least `b` with `L_col(b)_i = J(A)_i` for all `i <= m`. Whenever a column
family's image is a full box `[-L,L]^b`, its pullback is generated by the `y_q^{2L+1}`, which are already in `J_w`. So only
families with **non-box images** can matter, and the tables record image sizes.

## 2. Certification method (mod-p sandwich)

All mod-`p` work uses `p = 2147475331` (prime, `= 1 mod 2310`). The kernels of the projection ideals use
`P62 = prevprime(2^62)`.

1. **Generators are exact.** Width powers are integer forms. The projection ideals `J(pi_F A)` are computed in one of two
   ways:
   * by exact `fmpz` nullspaces (`loclib.jproj_gens`), used for all `b <= 2`, and for `b <= 4` in `colb.py` and
     `cyc_run.py`;
   * by a nullspace mod `P62` followed by rational reconstruction **and exact verification over Z on every point of
     `pi_F(A)`** (`jfast.py`, `jfast2.py`). Unverified vectors are discarded; no failures occurred.

   Every generator is therefore an exact element of `J(A)` over `Q`, and only its reduction mod `p` is used.
2. **Equality (sandwich).** Cumulatively in the degree,
   `dim_p L_{<=i} <= dim_Q L_{<=i} <= dim_Q J(A)_{<=i} = N_{<=i} - rank_Q E_{<=i} <= N_{<=i} - rank_p E_{<=i}`.
   * Here `E` is the evaluation matrix of reduced monomials on `A`.
   * Any set of rows of `L` gives a lower bound for `dim_p L`. That covers the isotypic computations below, which use
     only some of the rows.
   * So **cumulative equality of the outer two numbers proves `L = J(A)` over `Q` in all degrees `<= i`**.
3. **Non-equality** needs an exact **upper** bound on `dim_Q L` and an exact element of `J(A)` outside `L`. It is
   certified in §4 by:
   * exact `fmpz` ranks of the invariant projection `P_0(L_i)`, which is rational;
   * exact invariant polynomials that vanish on `A`, verified over `Z` at **every** point of `A`;
   * a mod-`p` rank increase, which is a lower bound and so is valid in the right direction.

   Non-equalities that are **only** mod `p` are labelled *(mod p, evidence)*.
4. **Isotypic decomposition** (`cyclib.py`). For Paley-12 (`Z_11`) and `Syl(16)` (the Singer cycle `Z_15`, with
   coordinates `alpha^l` in `F_16 = F_2[t]/(t^4+t+1)`), the coordinate shift preserves `A` and every family set used.
   * Over `F_p` (`p = 1 mod 11, 15`) every ideal splits into character components.
   * Equivariance gives `P_s(R_1 J_{i-1}) = sum_{s'} P_s(x_0 J^{s'}_{i-1})`, so one variable suffices.
   * Generators are needed only up to the shift.
   * Short orbits (for `q = 15`) are handled by periodic orbit rows.
   * The decomposition is validated by reproducing zetc's Paley Hilbert function `[1,11,66,275,880,2265,4785,7525,4345]`
     exactly.

## 3. Models and results

Dimensions are those of `J(A)_i` **in `R`**, i.e. modulo `(x^3)`. The judge's full-ring numbers are larger by
`dim (x^3)_i`.

### 3.1 Sylvester 8 (`d = 7`): `M8`, `M8'`, `L = 1`

Scripts: `colb.py syl8` (`colb_syl8_out.txt`), `m8_lines.py`, `mixed_syl8.py` (`mixed_syl8_out.txt`), `m8_cert.py`
(`m8_cert_out.txt`).

| model | `m` | `m/(2L)` | `H_A` | `dim J(A)_i` (in `R`), `i = 0..m` | `J_w(R=2)` | `b_0` (columns) | `b_0` (columns + coordinates) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Syl(8,1)` | 3, 4 | 1.5, 2 | `[1,7,28,21]` | `[0,0,0,56,161]` | equal | **1** | 1 |
| `M8 = Syl(8,2)` | 8 | 2 | `[1,7,28,77,161,259,287,175,42]` | `[0,0,0,0,0,7,70,218,315]` | `[..,7,49,162,294]` (fails at `i = 6`) | **3** | **3** |
| `M8' = Syl(8,3)` | 11 | 1.83 | `[1,7,28,77,161,266,357,386,329,195,63,7]` | `[0,..,0,7,28,71,98,70]` | **equal in all degrees `<= 11`** | **1** | 1 |

**Findings.**
1. **`M8'` is local with `b_0 = 1`.** `J(M8') = J_w(R=2)` in every degree `<= 11 = rho`. This is exact, by the sandwich.
   So the round-12 locality hypothesis holds in `M8'` in its strongest form: width powers alone suffice.
2. **`M8`: `b_0 = 3`, and the needed families are exactly the 7 lines.**
   * All 21 two-column images are full `5 x 5` boxes. So `L_col(2) = J_w` exactly (§4.1, Lemma 4.1).
   * The 28 independent triples have full `5^3` images.
   * Each of the 7 **line triples** `{a, b, a+b}` has a 57-point image with `rho = 6`. Its new generators are 6 in
     degree 5 and **3 in degree 6**, and `7 x 3 = 21` is exactly the degree-6 gap (`70 - 49`).
   * `J_w + (7 line pullbacks) = J(A)` in all degrees `<= 8`. This is exact.
3. **Mechanism: a parity congruence.** `1_{S_a} + 1_{S_b} + 1_{S_{a+b}} = 2 * 1_{complement of a.x = b.x = 0}`, so on `A`
   * `y_a + y_b + y_{a+b} = 2z`, with `z` an integer direction (`|<z,A>| = 7`);
   * the line image is the box `[-2,2]^3` cut to the even sublattice.

   Replacing the line by 2-dimensional families does **not** work. The 21 pairs `(column, z)` give only 2 of the 3
   degree-6 generators per line (`63 < 70`). Pairs `(z, z')` give nothing.
4. **Exact over `Q`** (`m8_cert.py`).
   * Let `L_mixed(2)` be `J_w(R=2)` plus the pullbacks of **all** 2-dimensional families from columns and coordinate
     vectors (84 families). Then `dim_Q L_mixed(2)_6 = 49` exactly (fmpz rank).
   * On the other side, `dim_Q J(A)_6 >= 70`.
   * So `b_0(M8) = 3` holds over `Q` for both family sets. The judge's "pairs add nothing" is now exact, and extended to
     mixed pairs.

### 3.2 Paley-12, `L = 2`, `m = 8` (`d = 11`, `|A| = 20153`, `rho = 8`, `m/(2L) = 2`): **not local for `b <= 5` (exact) or `b <= 6` (mod `p`)**

Scripts: `cyc_run.py paley 2 8 4` (`cyc_paley_L2_out.txt`), `paley_images.py`, `paley_cong.py`, `paley_b5.py`,
`paley_bgen.py` (`paley_bgen_out.txt`, `paley_bgen67_out.txt`), `paley_mixed.py` (`paley_mixed_out.txt`),
`paley_pair_block.py`, `paley_syz.py`, `paley_cert.py` (`paley_cert_out.txt`).

**Hilbert function and `J(A)`.**
* `H_A = [1,11,66,275,880,2265,4785,7525,4345]`. This is mod `p` and agrees with zetc.
* `dim J(A)_i = [0,0,0,0,0,12,132,1517,10010]` (`i <= 8`). These are mod-`p` upper bounds; the lower bounds come from §4.
* Per `Z_11`-component in degree 7, `dim J(A)^s_7 = 137` (`s = 0`) and `138` (`s != 0`), against `72` for `J_w`.

**Width ideal.**
* `J_w(R=1)`, with all `88573` primitive `{-1,0,1}` directions, has `12` quintics: the 11 columns and the all-ones
  direction `sigma = -k_base` (`|k_base| <= 2`).
* It has `17` orbits of septic width powers. These lie in `R_2 * (quintics)` mod `p`.
* `dim J_w = [..,12,132,792,3300]`. **It fails at `i = 7`**, with a deficit of 725 (65 or 66 per component).

**Column families, by image (`paley_images.py`).**

| `b` | image sizes (number of shift orbits) | full box `5^b`? |
| --- | --- | --- |
| 2, 3, 4 | `25 (5)`, `125 (15)`, `625 (30)` | **all full** |
| 5 | `3113 (30)`, `3123 (6)`, `1031 (5)`, `1041 (1)` | none full |
| 6 | `5073 (30)`, `5105 (6)`, `4971 (5)`, `5023 (1)` | none |
| 7..10 | `7917..7981`, `11781..11833`, `16325`, `20153` | none (`b = 10` is injective on `A`) |

* Hence `L_col(4) = J_w` **exactly** (§4.1, Lemma 4.1).
* The column matrix has Smith form `[1^5, 3^5, 6]`, i.e. a 6-dimensional space of **mod-3 congruences** among the columns
  (minimal support 5), plus `sum_p y_p = 6 sigma`.
* The six orbits with `|P| = 1031, 1041` are exactly the 5-column mod-3 congruences: `index-3` sublattices, SNF
  `[1,1,1,3]`.
* The near-full images (`3113`, `3123`) miss only corners cut by the ball `sum_p y_p^2 = 3 |k|^2` (zetc Lemma 4.4).

**Locality results (mod `p` lower bounds on `dim L`; deficits are evidence unless stated).**

| family set | `dim L_7` | `dim L_8` | deficit at `i = 7` |
| --- | --- | --- | --- |
| `J_w(R=1) = L_col(4)` | 792 | 3300 | 725 (**exact**, §4.2) |
| `L_col(5)` (all 462 five-column families) | 902 | 7546 | 615 (**exact** non-containment, §4.2) |
| `L_col(6)` (contains `L_col(5)`: pullbacks are monotone under refinement of the family) | 902 | 7612 | 615 |
| `L_col(7)` | not completed (estimated `> 10 h`) | | |
| columns + coordinates, `b <= 3` (all images full products) | 792 | 3300 | 725 |
| columns + coordinates, `b <= 4` | 792 | 5335 | 725 |
| `J(A)` | 1517 | 10010 | — |

* Among the 5-column families, only the all-ones mod-3 congruence `(0,1,2,6,9)` (1 orbit, image 1041) has degree-7
  generators: 11 per family. Its 11 shifts together raise `dim L_7` by 110.
* The other congruence families have generators only in degrees 5 and 8 (51 each in degree 8).
* Near-full images contribute only `y_q^5`, up to degree 8.
* Families with `rho(pi_F A) >= m = 8` do contribute generators below `m`: every 5-column family has
  `rho(pi_F A) >= 8`. So README §7.1's restriction to `rho(pi_F(A)) < ceil(4w)` would lose them.

**Cut boxes do not explain it either.**
* The box cut by two Paley columns is the two-cut block `(3,3,3)`, `L = (2,2)`, `d = 9`, `|A_T| = 12657`.
* It satisfies `J(A_T) = J_w(A_T)` in all degrees `<= 8` (exact, by the sandwich).
* Widths on `A_T` are at least those on `A`, so `J(A_T)_{<=8} ⊂ J_w(A)`. **Pair cut boxes contribute nothing.**

### 3.3 Two- and three-cut blocks (`d <= 7` and `d <= 6`; all `L`, resp. `L in {1,2}`)

Scripts: `cuts_run.py` (`cuts_two_d7_out.txt`, `cuts_three_d6_out.txt`), `three_mixed.py` (`three_mixed_out.txt`,
`three_mixed_summary.txt`), `mixed.py`.

The test runs in all degrees `<= rho + 1`, with `J_w(R=2)` for `d <= 6` and `R = 1` for `d = 7`. Free coordinates
factor off, so blocks are enumerated by Venn-region sizes.

**Two cuts.**
* 400 of the 402 blocks with `a + b + c <= 7` (all binding `L_1, L_2`) were completed. The two missing ones are the
  nested `(0,7,0)`, `L = (5,6), (6,6)`, cut off when the run was stopped. **Every completed block has `J(A) = J_w` in
  every degree `<= rho + 1`** (exact, by the sandwich).
* `rho` reaches 13. The ratio `rho/(2 L_max)` reaches 4.5, for example `(1,2,4)`, `L = (1,1)`, `rho = 9`.
* So two cuts are local with `b_0 = 1` at every ratio tested, including ratios above the regime's 4.
* The same holds for the Paley pair block `(3,3,3)`, `L = (2,2)`, `d = 9` (§3.2).

**Three cuts.** All 953 Venn types up to `S_3` (the same count as zetc §5.2), block size `<= 6`, `L_j in {1,2}` binding:

| outcome | number of blocks |
| --- | --- |
| `J(A) = J_w` (`b_0 = 1`) | 827 |
| needs the 3-column projection (`b_0(col) = 3`; never `2`) | 27 |
| **not generated by `J_w` + all column projections** (`b_0(col) = None`) | 99 |

* For the 126 non-width blocks, the first failing degree satisfies `rho - i_fail in {0,1,2,3,4}`.
* With families from columns and coordinates (`three_mixed.py`, all 126 non-width blocks), the minimal span dimension
  is `b_0(mixed) = 2` (6 blocks), `3` (81) or `4` (39) (`three_mixed_summary.txt`). It is **at most 4 in every
  three-cut block with `d <= 6`**.
  * The equalities at `b_0(mixed)` are exact (sandwich).
  * That `b_0 - 1` fails is mod `p`.
* In these blocks the projection-invisible generators therefore need coordinate directions: they are cut-box
  generators (Lemma 4.3) supported on a few coordinates.
* Example (`mixed.py`): `S_1 = {0,1,4}`, `S_2 = {0,2,3,4}`, `S_3 = {1,2,3,4}`, `L = (1,1,1)`, `d = 5`, `|A| = 87`.
  * `J(A)_4 = 18 > 15`, the value for `J_w`, for `L_col(2)` and for `L_col(3)` alike.
  * Mixed families of span 3 restore equality.
  * **Exact** (`three_cert.py`, `three_cert_out.txt`): `dim_Q L_col(3)_4 = 15` (fmpz), against
    `dim_Q J(A)_4 >= 18`, the mod-`p` dimension of the exact mixed local ideal. So no column family, not even the full
    3-column projection, generates `J(A)` here.

### 3.4 `Syl(16, 1)` (`d = 15`, `|A| = 241`, `rho = 3`; ball-bound) (`syl16_L1.py`, `syl16_L1_out.txt`)

* `H_A = [1,15,120,105,0]` and `dim J(A) = [0,0,0,560,2835]`.
* Width powers give 15 in degree 3; with the line triples, 155.
* Mixed families (columns and coordinates) of span `<= 2` give `330`. Span `<= 3` gives **equality in all degrees
  `<= 4`**, which is exact.
* So `b_0(mixed) = 3`; that `b = 2` fails is mod `p` evidence.
* Here the ball `|k|^2 <= 15/4` (zetc Lemma 4.4) allows at most 3 nonzero coordinates, so coordinate-triple images are
  non-grid. This is a feature of the small Hadamard models; in the regime the ball is slack.

### 3.5 Abstract Sylvester congruence images (the regime limit of line and plane families)

Scripts: `line_lattice.py` (`line_lattice_out.txt`), `plane_lattice.py` (`plane_lattice_L2_R1_out.txt`,
`plane_lattice_L2_R2_out.txt`).

When every coset is large (`Syl(2^k, j)` with `k >> s`), the image of the `2^s - 1` columns of an `s`-dimensional
subspace is exactly `Lambda_s ∩ [-j,j]^{2^s-1}`. Here `Lambda_s` is the `Z`-span of the simplex-code vectors
`(a.x)_a`. Its index is `2, 32, 131072` for `s = 2, 3, 4`, matching `2^{(s-2)2^{s-1}+1}`.

The test is `J(P)` against the width powers of `P`'s own dual-lattice directions, in the polynomial ring of the family.
All results are mod `p` (evidence).

* **Lines (`s = 2`, `b = 3`).**
  * `rho(P_L) = 3L`.
  * `J(P_L) = J_w(P_L)` in all degrees `< 3L`, and **fails exactly at `3L`**. This holds for `L = 2, ..., 6`
    (e.g. `L = 6`: `126` against `63` in degree 18).
  * So every line family has non-width generators in degree `3j`. That is below the regime order `m' ~ 8j`.
* **Planes (`s = 3`, `b = 7`, `L = 2`).**
  * `|P| = 2605`, `H_P = [1,7,28,84,210,455,812,756,252]`, `rho = 8`.
  * `dim J(P) = [..,7,112,960,2751,5005]`.
  * `J_w` (`R = 1` or `2`) fails at degree 6 (`49`).
  * **`J_w + the 7 line pullbacks` restores degree 6 (`112`) but fails at degree 7 (`582` against `960`).** So the plane
    family has generators that its lines do not produce.
* **Count heuristic** (labelled heuristic). By the dimension count `C(D+b, b) >= |P|`,
  `rho(P_s) ~ (b!/index_s)^{1/b} (2j+1)`. That gives about `2.9j, 4.1j, 5.9j, 8.3j, 11.8j` for `s = 2..6`.
  * If the non-width generators of `P_s` sit near `rho(P_s)`, as they do for lines (and, at `L = 2`, for planes), then
    Sylvester subspace families up to `s = 4` (**15 columns**) can matter below `m' ~ 8j`, and `s = 5` (31 columns) is
    borderline.
  * Whether the actual box's width powers absorb them, as happens in `M8'` but not in `M8`, is not decided.

### 3.6 `Syl(16, 2)` (`d = 15`, `|A| = 213937`, `T = 1`, degrees `<= 7`; Singer cycle `Z_15`) (`syl16.py`, `syl16_L2_m7_out.txt`)

**Data.**
* `H_A = [1,15,120,665,2835,9813,28125,62186]`. This is mod `p`, from `hilbert_lean`. The per-component column subsets
  only lower ranks, so `dim J(A)` below is an upper bound, which is the sound direction.
* `dim J(A)_i = [0,0,0,0,0,15,540,9769]` in `R`.
* Family set: the width powers of the columns (`|image| = 5`), of `sigma = -k_base` (`|image| = 7`), of the 35 line
  half-sums and of the 15 plane quarter-sums (all `|image| = 7`). Column families are added on top.

| family set | `i = 5` | `i = 6` | `i = 7` |
| --- | --- | --- | --- |
| `J(A)` | 15 | 540 | 9769 |
| width powers as above | 15 | **225** | 1851 |
| + 35 line triples | 15 | **540 (equal)** | 6226 |
| + 15 plane 7-column families (1 Singer orbit) | 15 | 540 | **9691** (deficit 78) |

1. **Lines are needed and are exactly enough in degree 6, i.e. `b_0 >= 3`.**
   * Every line image is the **full** abstract `P_2 = [-2,2]^3 ∩ {sum even}` (63 points; regions of 4 points). It
     contributes 9 non-width generators in degree 6 (§3.5).
   * `35 x 9 = 315 = 540 - 225`: they are independent, and **none is absorbed by the box's width powers**.
   * Exact scan (`syl16_widthscan.py`, `syl16_widthscan_out.txt`, integer arithmetic): among all `7174228` primitive
     `v in {-1,0,1}^15` of support `>= 3`, the only ones with `|<v,A>| <= 6` are the 15 columns.
   * So `J_w(R=1)` agrees in degrees `<= 6` with the width ideal used here. Smaller supports give powers that vanish
     in `R`.
   * The `R = 1` width ideal therefore misses all 315 line generators in degree 6.
   * All 435 two-dimensional families from columns and coordinates have **grid** images (`syl16_mixed2.py`,
     `syl16_mixed2_out.txt`). By Lemma 4.1, `L_mixed(2)_6` is spanned by the products `R_1 *` (15 column quintics).
   * **Exact:** `dim_Q L_mixed(2)_6 <= 15 * 15 = 225 < 540 <= dim_Q J(A)_6`. The lower bound comes from the mod-`p`
     dimension of the exact line ideal.
   * So `b_0(Syl(16,2)) >= 3` over `Q`, for columns and for columns + coordinates, and `b_0 = 3` suffices through
     degree 6.
   * Equality through degree 6 is exact, by the sandwich.

   This is the first Sylvester model with regime-type (full) line images. In it the M8'-type absorption does **not**
   happen.
2. **Degree 7.** Lines leave a deficit of `3543` (mod `p`). The plane images (`|P| = 2507`, not full; regions of 2
   points) have non-width generators `[(5,14),(6,63),(7,232)]`. **Adding the plane families raises degree 7 from 6226 to
   9691**, so 7-column families are needed in degree 7 (`b_0 >= 7` for column families, mod `p` evidence).
   * The plane contribution is not generated by lines, matching the abstract plane image (§3.5).
   * **The full plane is needed.** All 4-, 5- and 6-column subfamilies of a plane, added to the lines, give nothing in
     degree 7: it stays at `6226` (`syl16_planesub.py`, `syl16_planesub_out.txt`).
   * The mechanism is the mod-4 congruence `sum_{a in V} y_a = 4 z_V` (Lemma 4.2), which involves all 7 columns.
   * Families outside planes add nothing on top of the planes. Of the `b = 3, 4, 5` families outside planes, `0`, `56`
     (all grids) and `180` (non-grid) orbits respectively change nothing (`syl16_rest_out.txt`,
     `syl16_rest56_out.txt`). The run that added the 328 non-grid 6-column orbits on top of the planes was OOM-killed;
     see the incremental run below.
   * Without the planes, the 180 non-grid 5-column families outside planes, added to the lines, also give **nothing**:
     degree 7 stays at `6226` (`syl16_noplane5_out.txt`).
   * **Incremental run** (`syl16_inc.py 6`, `syl16_inc6_out.txt`). Since `L_6 = J(A)_6` already, only degree-7
     generators matter. Adding to width + lines all 328 non-grid 6-column orbits outside planes leaves degree 7 at
     **`6226`**. Adding the planes then gives `9691`.
   * Hence, mod `p`, **`L_col(6)_7 = 6226 < 9769`**. The pieces are: lines, in-plane subfamilies of size `<= 6`
     (nothing new), non-plane families of size 4 (grids), 5 and 6 (nothing new).
   * The plane families (`b = 7`) add `3465`.
   * A residual of 78 remains in degree 7 after planes (`syl16_planes_out.txt`). Its source (families of `>= 8`
     columns, or mixed families) is not identified.
   * **Mod `p`, then, column-`b_0(Syl(16,2), m = 7) >= 7`.** The 7-column plane families are genuinely needed, and even
     they leave 78.

## 4. Non-local generators

### 4.1 Three elementary lemmas (proved)

**Lemma 4.1 (grids are invisible).** Let `F` be a family with `pi_F(A) = I_1 x ... x I_b`, a product of finite sets.
Then `pi_F^* J(pi_F A) = (<v_q, x>^{|I_q|} : q <= b)`. Every generator is a width power, so the family contributes
nothing beyond `J_w`.

*Proof.*
1. `I(I_1 x ... x I_b) = (p_q(y_q))_q`, with `p_q = prod_{t in I_q}(y_q - t)`. This is the standard grid fact: reduce
   modulo the `p_q` to degree `< |I_q|` in each `y_q`, then use induction on `b`.
2. Fix a degree-compatible monomial order `<`. The `p_q` have pairwise coprime leading monomials `y_q^{|I_q|}`, so by
   Buchberger's first criterion they form a Gröbner basis, and `in_<(I) = (y_q^{|I_q|})`.
3. For a degree-compatible order, `in_<(J(I)) = in_<(I)`. So `J(I) ⊇ (y_q^{|I_q|})` are two homogeneous ideals with
   the same initial ideal. Hence they are equal: `J = (y_q^{|I_q|})`.
4. Pull back along `y_q = <v_q, x>`. A pulled-back form is the leading form of the pulled-back polynomial, or it is
   `0`.

QED

**Lemma 4.2 (congruences give integer directions).** Suppose `A` contains every unit vector `e_x` (true whenever all
`L_p >= 1`), and `sum_q c_q <v_q, k> = 0 mod M` for all `k in A`. Then `z = (sum_q c_q v_q)/M` is in `Z^d`.

*Proof.* Evaluate at `k = e_x`. QED

So a non-grid image produced by a congruence is the image of a family whose span contains the integer direction `z`.
By span invariance (§1), the family may be rewritten using `z`.

**Lemma 4.3 (the cube-and-slab ideal).** For integer columns,
`I(A) = (x_l^3 - x_l (all l), g_p(y_p) (all p))`, where `g_p(y) = prod_{|t| <= L_p}(y - t)`.

*Proof.*
1. Let `B = Q[x]/(x^3 - x)`. It is the algebra of functions on the cube `C`, a product of copies of `Q`, one per point.
2. In such an algebra the ideal generated by the images of the `g_p` is the ideal of functions vanishing on their
   common zero set `C ∩ ∩_p S_p = A`.

QED

Hence `J(A) ⊇ (x^3, y_p^{2L_p+1})`, and **equality holds iff `{x^3 - x, g_p}` is an H-basis.** Every generator of `J(A)`
outside the width ideal is therefore an H-basis completion: the leading form of an S-polynomial `sum_p c_p g_p`, reduced
by `x^3 -> x`, built from a syzygy `sum_p c_p y_p^{2L+1} = 0` in `R` that is not Koszul. Since `x^3 -> x` lowers the
degree by 2, a degree-`D` syzygy produces candidates in degree `D - 2` (then `D - 4`, and so on).

### 4.2 The Paley non-local generator: exact certificate (`paley_cert.py`, `paley_cert_out.txt`)

**Proposition 4.4 (exact).** Paley-12, `T = 1`, `L = 2`, `m = 8`. Over `Q`:
* (i) `J(A)_7 ⊄ L_col(4)_7 = J_w(R=1)_7`;
* (ii) `J(A)_7 ⊄ L_col(5)_7` (all five-column families);
* (iii) in the `Z_11`-invariant part, `dim_Q P_0(J(A)_7) = 137`, `dim_Q P_0(L_col(4)_7) = 72` and
  `dim_Q P_0(L_col(5)_7) = 82`.

So column-`Loc(b_0)` fails in degree `7 < m` for every `b_0 <= 5`. The width ideal uses every `{-1,0,1}` direction.

*Certificate.*
1. Let `P_0` denote orbit-sum coordinates (the `Z_11`-invariant part, which is rational).
2. The 149 integer rows `P_0(mu * g)` span `P_0(J_w(R=1)_7)`. Here `g` runs over the exact width generators (2 quintic
   and 17 septic orbit representatives), and `mu` over the monomials of `R_{7 - deg g}`.
3. Their **exact** rank is `72` (fmpz).
4. A nullspace mod `P62` of the invariant evaluation matrix (1589 orbit-sum monomials of degree `<= 7` against 1834
   point orbits), followed by rational reconstruction, yields an invariant polynomial `f` of degree 7 with integer
   coefficients (max 48).
5. `f` vanishes at **all 20153 points of `A`** (checked over `Z`). Its top form `T_f` (384 nonzero coefficients) lies in
   `J(A)_7`.
6. `rank_p([rows; P_0(T_f)]) = 73 > 72`, which forces the rational rank above 72 as well. So `T_f ∉ J_w(R=1)_7` over
   `Q`.
7. By Lemma 4.1, `L_col(4) = J_w(R=1)`, because all images of `<= 4` columns are full boxes (`paley_images.py`).

The coefficients of `f` are in `paley_cert_poly.txt`.

For (ii) and (iii) (`paley_cert3.py cong`, `paley_cert3_cong_out.txt`):
1. The rows now also contain the exact (verified) pullbacks of the six congruence families. Their exact rank is `82`.
   * **Completeness, exact.** For each of the 6 congruence orbits, the ideal generated mod `p` by the exact generators
     has the same cumulative dimension in degrees `<= 7` as `J(P)` mod `P62`, which is an upper bound for `dim_Q J(P)`
     (`paley_cong_check.py`). The same holds for the 36 near-full orbits, where the comparison is with the grid ideal
     `(y_q^5)` (`paley_nearfull_check.py`). By the sandwich, the rows span `P_0(L_col(5)_7)` over `Q`.
2. Pick a nonsingular pivot block `E[pr, pc]` of the invariant evaluation matrix: rank 1438 mod `p`, hence nonsingular
   over `Q`.
3. Exact Dixon solves give 55 invariant polynomials of degree `<= 7`, each verified to vanish at all 20153 points over
   `Z`. Their top forms raise the mod-`p` rank of the rows from 82 to 137.
4. This yields (ii), and `dim_Q P_0(J(A)_7) >= 137`.
5. Equality in (iii): the mod-`p` kernel of the invariant evaluation in degrees `<= 7` has dimension 151, an upper bound
   for the rational one. The rational one is at least `2 + 12 + 137 = 151` (degrees 5, 6, 7).

The first Dixon polynomial is in `paley_cert_poly_cong_dixon.txt`.

**Caveat on (ii).** The saved polynomial of the `b <= 4` certificate has its top form **inside** `L_col(5)`: it comes from
the all-ones mod-3 congruence family. The certificate for (ii) uses different polynomials. QED

**Size of the phenomenon (mod `p`, evidence for the exact sizes).**
* In each `Z_11`-component `J(A)^s_7` has dimension `137`/`138` against `72` for `J_w`: a deficit of `65`/`66`
  (`= dim R_2`) per component, 725 in total.
* The 12 quintics have **440 non-Koszul syzygies in degree 9 in `R`** (`paley_syz.py`): the ideal they generate has
  `dim = [..,12,132,792,3300,10120,23353]`, against the free count `10560` in degree 9.
* By Lemma 4.3 these syzygies produce degree-7 S-polynomial candidates. Since `440 < 725`, second-generation
  completions (degree-11 syzygies, a drop of 4) are needed as well.
* The syzygies mix the cube relations with **several slabs at once**. Projections cannot see them, because every image of
  `<= 4` columns, and of `<= 3` columns-or-coordinates, is a grid.

## 5. Trend in m/(2L)

**Cut blocks** (`ratio_table.py`). Columns: `rho/(2 L_max)` bins. Entries: number of blocks by outcome.

| `rho/(2L_max)` | two cuts: `b_0 = 1` | three cuts: `b_0 = 1` | needs 3-col projection (`b_0(mixed)`) | needs coordinates: `b_0(mixed) = 3` | `b_0(mixed) = 4` |
| --- | --- | --- | --- | --- | --- |
| 1.0-1.5 | 100 | 12 | 2 (2) | 0 | 0 |
| 1.5-2.0 | 42 | 236 | 11 (2 or 3) | 20 | 6 |
| 2.0-2.5 | 23 | 411 | 14 (3) | 31 | 18 |
| 2.5-3.0 | 3 | 31 | 0 | 0 | 1 |
| 3.0-3.5 | 4 | 30 | 0 | 3 | 7 |
| 3.5-4.0 | 5 | 69 | 0 | 5 | 7 |
| 4.0-4.5 | 5 | 34 | 0 | 1 | 0 |
| 4.5-5.0 | — | 4 | 0 | 0 | 0 |

(`L_max` is the largest slab half-width of the block; there are 19 further two-cut blocks in 0.5-1.0.)

**Hadamard models.**

| model | `m/(2L)` | first non-width degree `i_fail` (`rho`) | `b_0` (columns) | `b_0` (columns + coordinates) |
| --- | --- | --- | --- | --- |
| `Syl(8,1)` | 2 | none (3) | 1 | 1 |
| `M8 = Syl(8,2)` | 2 | 6 (8) | 3 (lines) | 3 |
| `M8' = Syl(8,3)` | 1.83 | none (11) | 1 | 1 |
| `Syl(16,1)` | 1.5 | 3 (3) | — | 3 |
| Paley-12, `L = 2` | 2 | 7 (8) | `> 5` (mod `p`); `> 4` exact | `> 4` (mod `p`) |

**Reading.**
1. **`b_0` is not a function of `m/(2L)`.** At ratio 2, the four Hadamard models give `b_0 = 1, 3, 1` and `> 5`. In the
   cut blocks the share of non-width blocks does not grow with the ratio: it is about 10% around ratio 2, 18% in
   3.0-4.0, and 1 of 39 in 4.0-5.0. Two-cut blocks are width-local at every ratio up to 4.5, which includes the regime
   value 4.
   * Among the non-width three-cut blocks there is a mild drift towards `b_0(mixed) = 4` as the ratio grows: 14 of 22
     at ratio 3-4, against 6 of 37 at ratio 1.5-2.
   * `b_0(mixed)` stays `<= 4 < d` throughout.
2. What does determine `b_0`:
   * (a) **arithmetic congruences among columns** (Lemma 4.2). Sylvester lines are mod 2 and need `b = 3`. Paley's
     5-column relations are mod 3.
   * (b) **H-basis completions across several slabs** (Lemma 4.3), which projections onto columns do not see. These
     appear in every non-width three-cut block that needs coordinates, and they make up the bulk of Paley's deficit.
3. **Position of the failures.** Non-width generators start at `i_fail >= rho - 4` in every model. In the abstract
   Sylvester line image the first non-width generators sit at degree `3L = rho(P_L)` (`line_lattice.py`, `L = 2..5`).
   * The regime has `m'/rho -> 0`, so a failure tied to the top degrees of the *full* box would be harmless there.
   * Failures tied to `rho` of a *projection* (like the line's `3L ~ 3j < m' ~ 8j`) are not harmless.

## 6. What this says about Loc(b0) in the regime

**6.1 `Loc(b_0)` as stated in README §7.1 is refuted as a model-independent principle.** It fails on an exact model,
Paley-12 (`L = 2`, `m = 8`):
* in degree `7 < m`, for every `b_0 <= 5`, over `Q` (Proposition 4.4);
* for `b_0 <= 6` mod `p` (`L_col(6)_7 = 902 < 1517`). `b = 7` was not completed.

The model is zonotope-free, overdetermined and indecomposable, but **LP-subcritical, and ball-bound** (zetc Lemma 4.4).
Two further points:
* README §7.1's restriction to families with `rho(pi_F A) < ceil(4w)` makes `Loc` weaker still: in Paley the only
  families with non-width generators below `m` have `rho(pi_F A) >= 8 = m`. This is exact, by counting:
  `|pi_F A| >= 1031 > 792`, the number of monomials of degree `<= 7` in 5 variables.
* So the projection-cost bound proposed in README §7.3(3), "a family contributes generators in degree `<= m'` only if
  `rho(pi_F(A)) < m'`", is **false** as a general statement.

What survives:
* `Loc` with `b_0 = 1` holds in `M8'`, `Syl(8,1)`, every two-cut block and 827/953 three-cut blocks.
* `b_0 = 3` holds in `M8` (lines) and in `Syl(16,1)` (mixed).
* `b_0 <= 4` holds in every three-cut block once coordinates are allowed.

**6.2 The non-local generators have one algebraic source** (Lemma 4.3). They are H-basis completions of the slab system
`{x^3 - x, g_p(y_p)}`, i.e. leading forms produced by **non-Koszul syzygies of the width powers `y_p^{2L+1}` in
`R = S/(x^3)`**. A completion is "local" in README's sense only when it already appears in the leading-form ideal of a
projection onto a few directions. That happens when:
* (a) the syzygy comes from an **arithmetic congruence** among few columns (Lemma 4.2), whose projection image is a
  sublattice-in-box, not a grid. Examples: Sylvester lines (mod 2, `b = 3`), Paley 5-column relations (mod 3, `b = 5`);
* (b) the completion lives on a cut box `C ∩ S_T` with **small coordinate support** (the three-cut blocks with
  `b_0(mixed) <= 4`).

Paley's degree-7 deficit is neither:
* its column projections of size `<= 4` are grids;
* its pair cut boxes are width-generated (exact);
* its 440 degree-9 syzygies, and the further completions, couple `>= 3` slabs whose supports cover 9-11 of the 11
  coordinates.

In the language of the round-12 README, these are **global low-degree relations coupling many columns**.
* They can only occur where the width powers are *not Fröberg-generic*. In Paley the expected count
  `N_9 - 12 N_4 = 9295 > 0` predicts no degree-9 syzygy, yet there are 440.
* Non-genericity is necessary but **not sufficient**, and it is common in small models (`froberg_check.py`, mod `p`):
  * in `M8` the 7 quintics span `154 < 196` in degree 7 and `287 < 357` in degree 8;
  * in `M8'` the 7 septics span `28 < 49` in degree 8, yet `J(M8') = J_w`, because the completions are absorbed by
    other width powers.

**6.3 Regime implications (heuristic, labelled as such).**
1. *Congruence families.* Admissible integer profiles (Sylvester, the J5 regime box) carry congruences.
   * The abstract Sylvester line image `[-j,j]^3 ∩ {sum even}`, which is the exact image once cosets are large, has
     non-width generators exactly at its `rho = 3j` (computed for `j = 2..6`, `line_lattice_out.txt`). Since
     `3j < m' ~ 8j`, line families are needed unless the width powers of `A` absorb them.
   * They are absorbed in `M8'` (small images) and not in `M8`.
   * **They are not absorbed in `Syl(16,2)`**, the first Sylvester model whose line images are full (§3.6, exact). All
     35 x 9 line generators are needed in degree `6 = 3L`, and no width power or 2-dimensional family supplies them.
   * Planes (7 columns) add generators beyond their lines, both in the abstract image (`582 < 960`, §3.5) and in
     `Syl(16,2)` itself (degree 7: `6226 -> 9691`, §3.6; mod `p`).
   * So along Sylvester profiles **`b_0` grows with the subspace dimension `s` (`b = 2^s - 1`) as the degree rises**:
     `b = 3` at degree `3L`, `b = 7` at degree `3L + 1` for `L = 2`.
   * The count heuristic of §3.5 puts `s = 4` (15 columns) below `8j`. This is evidence for a **bounded but large**
     `b_0`, at least 7 and plausibly 15, for Sylvester regime boxes. It is not a proof in either direction.
2. *Generic real profiles.* For profiles with **real** columns, `<beta_p, k>` takes exponentially many values on `A`.
   * No column width power exists below degree `m'`, so every element of `J(A)_{<= m'}` beyond `x^3` is an H-basis
     completion or an integer-direction width power.
   * `Loc(b_0)` then asks whether these completions are visible to projections onto few **integer** directions.
   * Nothing computed here bears on that case: every exact model has integer columns.
3. *Where the failures sit.*
   * In every exact model except `Syl(16,2)` the first non-width degree satisfies `rho(A) - i_fail <= 4`.
   * `Syl(16,2)` has `i_fail = 6` and `8 <= rho <= 11`. The upper bound is the referee's; the lower bound holds because
     `sum_{i<=7} H_A <= N_{<=7} = 114084 < 213937 = |A|`.
   * In the regime `m'/rho(A) -> 0`, so these failure *positions* do not transfer.
   * But the congruence mechanism attaches failures to the `rho` of a **projection** (`3j` for a line), and that does
     transfer.

**6.4 Bottom line for route 5.** `Loc(b_0)` cannot be proved as a general structural fact about slab boxes: Paley-12 is
an exact counterexample at small `b_0`. Any locality statement for the regime must use regime-specific input.
* One option is a Fröberg-type genericity of the width powers `{x^3, y_p^{2j+1}}` in degrees `<= m' + 2`.
* Another is an explicit list of the arithmetic congruence families, with their projection-local completions.

The round-12 hope that `Loc` reduces route 5 to few-cut inverse systems survives only in that conditional form. The key
negative structure (§4.2) is the type of global relation README §7.2 singled out as the possible source of
nondegenerate tangent-cone points on fake pattern boxes.

## 7. Scripts and reproduction

All scripts are in `loccompute_checks/`. Each writes `<name>_out.txt` (or the file named below) next to itself. Working
notes are in `loccompute_notes.md`.

| script | what | status of its claims |
| --- | --- | --- |
| `loclib.py` | reduced ring `R`, Hilbert functions, exact `J(P)` (fmpz), pullbacks, mod-`p` / exact ideal dimensions | library |
| `bigring.py`, `fullring.py` | `R` for `d = 15`; full polynomial ring | library |
| `cyclib.py` | `Z_q` isotypic ideal dimensions (short orbits handled), isotypic Hilbert function | library; validated on Paley |
| `hilbert_lean.py` | memory-lean isotypic Hilbert function (column subsets, which keep the sound direction) | library; validated on Paley |
| `jfast.py`, `jfast2.py`, `echelon.py` | `J(P)` via nullspace mod `P62`, rational reconstruction, exact verification | library (sound) |
| `crtker.py` | exact verification helper (int64 fast path only when overflow is impossible) | library |
| `models.py`, `locrun.py`, `mixed.py` | models, drivers, mixed families | library |
| `colb.py syl8`, `m8_test.py`, `m8_lines.py`, `mixed_syl8.py` | §3.1 | exact (sandwich) for equalities; mod `p` for `M8` pairs |
| `m8_cert.py` | §3.1: `dim_Q L_mixed(2)_6 = 49 < 70` | **exact** |
| `cyc_run.py paley 2 8 4`, `paley_images.py`, `paley_cong.py`, `paley_b5.py`, `paley_bgen.py`, `paley_mixed.py` | §3.2 | mod `p` (deficits = evidence) |
| `paley_pair_block.py` | §3.2: pair cut block is width-generated | exact (sandwich) |
| `paley_syz.py` | §4.2: 440 degree-9 syzygies | mod `p` (evidence) |
| `paley_cert.py`, `paley_cert3.py cong`, `paley_cong_check.py`, `paley_nearfull_check.py` | Proposition 4.4 | **exact** |
| `cuts_run.py two/two_a1/three`, `three_mixed.py`, `ratio_table.py` | §3.3, §5 | exact for `b_0 = 1` equalities; failures mod `p` |
| `syl16.py` (Hilbert function + first run; OOM-killed after the plane `J(P)` step), `syl16_planes.py`, `syl16_rest.py 3,4` / `5,6`, `syl16_planesub.py`, `syl16_noplane.py 5`, `syl16_inc.py 6` | §3.6 | equalities exact; degree-7 deficits mod `p` |
| `syl16_widthscan.py`, `syl16_mixed2.py` | §3.6: only the columns have width `<= 5`; all mixed pairs are grids | **exact** (integer scan) / grid check exact |
| `syl16_L1.py` | §3.4 | as stated |
| `three_cert.py` | §3.3 example: `dim_Q L_col(3)_4 = 15 < 18` | **exact** |
| `froberg_check.py` | §6.2: column width powers versus the Fröberg count | mod `p` |
| `line_lattice.py`, `plane_lattice.py` | §3.5 | mod `p` (evidence) |
