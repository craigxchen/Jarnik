# Referee report on loccompute.md (round 13, Route 5, locality Loc(b0))

Status: complete. Referee: adversarial, independent code (`referee_loccompute_checks/`), notes in
`referee_loccompute_notes.md`.

## 0. Verdict

**The exact computations are sound, but the headline conclusions overreach, because every `b_0` in the write-up is
measured on a restricted family set (columns, or columns plus non-base coordinates) while README §7.1's `Loc(b_0)`
counts arbitrary integer directions.**

* **Survives (re-verified independently):** the sandwich method; Lemmas 4.1-4.3 as lemmas; `M8` (`b_0 = 3`, exact),
  `M8'` (`J = J_w` through degree 11, exact), two-cut and three-cut equalities spot-checked; **Proposition 4.4 for the
  family set `L_col(5)`**, including an independent exact `dim_Q P_0 L_col(5)_7 = 82 < 137`.
* **Does not survive as stated:**
  1. "`Loc(b_0)` ... fails in two [models], so it is not a structural property of `T = 1` slab boxes" / "refuted as a
     model-independent principle" (§0, §6.1, §6.4). Only **column** locality is refuted (Paley, `b <= 5`, exact).
     Wider families recover generators the write-up calls non-local (Issue 1): pairs give README-`b_0 = 2` in every
     three-cut block tested (the §3.3 example with `b_0(mixed) = 3`, and 10 of the 39 blocks with `b_0(mixed) = 4`); pairs give at least 280 of the 315 `Syl(16,2)`
     line generators; base-coordinate families `{k_base, x_a, x_b, x_c}` give invariant elements of the Paley `J(A)_7`
     outside `L_col(5)` (exact).
  2. `Syl(16,2)` degree 7 (Issue 3): "whole 7-column planes are needed" holds only for column subfamilies: in the
     plane's ring, families of `<= 3` integer directions give 147 of the 231 plane-only generators (mod `p`). (The
     degree-7 width set also omits 105 width-7 directions, Issue 2, but their powers lie in the line pullbacks.)
  3. The regime heuristic "`b_0` bounded in `n` but at least 7, plausibly 15" counts columns of subspace families, not
     directions; it is unsupported.
* Nothing here proves or refutes README `Loc(b_0)` in any model. No growth or uniformity claim is made by the write-up
  or by this report; `M <= (2/3 + eps) log R / loglog R` stands.

## 1. What was checked and how

All checks below use my own code in `referee_loccompute_checks/` (no import of the author's scripts). Library:
`rlib.py` (exact leading-form ideals `J(P)` of finite point sets via `fmpz` nullspaces), `reng.py` (the reduced ring
`R = Q[x]/(x^3)` as dense vectors mod `p = 2147483629`, multiplication by linear forms, pullbacks, and the `Z_q` orbit
projection `P_0`). Rigour conventions are the author's (§2 of the write-up), which I re-derived and accept:

* a mod-`p` rank of integer rows is a **lower** bound for the rational rank, so a mod-`p` dimension of an ideal generated
  by exact elements of `J(A)` is a lower bound for its rational dimension;
* a mod-`p` rank of the evaluation matrix is a lower bound for the rational rank, so `N_i - H_A^{(p)}(i)` is an
  **upper** bound for `dim_Q J(A)_i`;
* equality of the two (cumulatively) proves equality over `Q`; a strict inequality between them proves nothing.

Models rebuilt from scratch: Paley-12 `L = 2` (`paley_base.py`), `M8 = Syl(8,2)` (`m8.py`), `Syl(16,2)`
(`syl16_base.py`), and the abstract Sylvester line image (`t1_line.py`).

## 2. Re-derivation of definitions and lemmas

**Family sets versus README §7.1.** README `Loc(b_0)` quantifies over families of `<= b_0` **arbitrary integer
directions** (with `rho(pi_F A) < ceil(4w)`). The write-up's `b_0` is always computed for the restricted sets
`L_col(b)` (columns) and `L_mixed(b)` (columns plus the `d` non-base coordinate vectors), plus `J_w(R)` with
`R in {1,2}`. The base coordinate `k_base = -sigma`, the line/plane half- and quarter-sums, the coset sums, the
mod-3 quotient directions and every other integer direction are **not** in either family set (except as width powers).
So every lower bound "`b_0 >= ...`" in the write-up is a lower bound for **column** or **column+coordinate**
locality only. The write-up mostly says so in its tables, but the summary, §0 item 5, §6.1 and §6.3 draw README-level
conclusions from it (Issues 1-3 below).

**Span invariance, `b = 1`, the sandwich.** Correct as stated: `J(TP) = (T^{-1})^* J(P)` for `T in GL_b(Q)`; the
pullback depends only on the rational span; `J(pi_v A) = (y^{|<v,A>|})`; the sandwich inequalities run in the
directions claimed.

**Lemma 4.1 (grids).** Correct. `I(I_1 x ... x I_b) = (p_q(y_q))`; the `p_q` have coprime leading monomials, so they
are a Gröbner basis for every degree-compatible order; `in(J(I)) = in(I)` gives `J(I) = (y_q^{|I_q|})`. "Contributes
nothing beyond `J_w`" requires each `v_q` to be in the direction set of `J_w`; for columns and coordinates it is.

**Lemma 4.2 (congruences).** Correct (evaluate at `e_x`, which lies in `A` when all `L_p >= 1` and the columns are
0/1). Consequence worth recording: the quotient directions can be very short. In Paley,
`y_5 + y_6 - y_7 + y_8 - y_10 = 3 (x_6 + x_9)` on `Z^11` (`paley_dirs.py`): the mod-3 congruence families contain a
**two-coordinate sum** in their span. The mod-3 space has `3^6 - 1 = 728` nonzero elements with support histogram
`{5:132, 6:132, 8:330, 9:110, 11:24}` (confirmed).

**Lemma 4.3 (cube-and-slab ideal).** The identity `I(A) = (x^3 - x, g_p(y_p))` is correct for integer columns (a finite
product of fields; `g_p(y_p(k)) = 0` iff `|y_p(k)| <= L_p` because `y_p(k)` is an integer). The text after "Hence" is an
interpretation, not part of the lemma: "every generator outside the width ideal is the leading form of an
S-polynomial built from a non-Koszul syzygy" is the standard H-basis (Möller-Sauer) criterion and is fine as a
heuristic description, but the quantitative statements attached to it (440 syzygies, "couple at least 3 slabs on 9-11
coordinates", "can only occur where the width powers are not Fröberg-generic") are mod-`p` evidence or unproved and
should be labelled so. In particular a mod-`p` ideal dimension is a lower bound, so "440 syzygies" is an **upper**
bound for the number of rational syzygies.

## 3. Independent computations

### 3.1 Reproduced (agree with the write-up)

| claim | my number | script | status |
| --- | --- | --- | --- |
| Paley-12 `L = 2`: `|A| = 20153`, `sigma = -k_base` has width 5 | same | `paley_base.py` | exact |
| Paley widths: coordinates 3; 11 columns + `sigma` 5; 17 orbits of width-7 `{-1,0,1}` directions (supports 5, 7, 10) | same | `paley_base.py` | exact (full scan of 88573 directions) |
| Paley invariant part: `dim P_0 J(A)_i <= 2, 12, 137` (`i = 5, 6, 7`) | same | `paley_inv.py` | mod-`p` upper bound |
| `dim P_0 (J_w(R=1))_7 = 72`, `dim P_0 L_col(5)_7 = 82`; only the orbit of `(0,1,2,6,9)` (image 1041) adds (10) | same | `paley_run1.py` | mod-`p` lower bound (author's 82 is an exact upper bound) |
| `M8`: `|A| = 1037`, `dim J(A)_{5..8} = 7, 70, 218, 315`, `J_w(R=2) = 7, 49, 162, 294`, `J_w + 7 lines = J(A)` | same | `m8.py`, `m8_run.py` | equality exact (sandwich) |
| `Syl(16,2)`: `|A| = 213937`, every line image is the full 63-point `P_2`, plane image 2507 points | same | `syl16_base.py` | exact |
| three-cut example (`S_1={0,1,4}`, `S_2={0,2,3,4}`, `S_3={1,2,3,4}`, `L = 1`): `|A| = 87`, `dim J(A)_4 = 18`, `J_w = L_col(3) = 15` in degree 4 | same | `three_ex.py` | exact |

Minor discrepancy: Paley has **1833** `Z_11`-orbits of points (`20153 = 11 * 1832 + 1`, the only fixed point is `0`), not
1834 as stated in §4.2 step 4. It does not affect any rank.

The logic of Proposition 4.4 is sound: invariant forms `g` satisfy `P_0 g = g`, so an exact invariant element of
`J(A)_7` outside `P_0(L_7)` lies outside `L_7`; `82` is an exact upper bound only if the rows span `P_0(L_col(5)_7)`,
and the completeness check (mod-`p` ideal of exact generators = mod-`P62` upper bound for `J(pi_F A)` in degrees
`<= 7`) gives that. My independent mod-`p` value 82 is the matching lower bound. **Proposition 4.4 is correct as a
statement about the family set `L_col(5)` (columns only, plus `J_w(R=1)`).**

### 3.2 New: families of non-column integer directions (README §7.1 allows them)

**(a) Abstract Sylvester line image `P_2 = [-2,2]^3 ∩ {sum even}` (63 points; this is the exact line image in
`Syl(16,2)`)** (`t1_line.py`, `t1_line_detail.py`). In `Q[a,b,c]`, degree 6:

| ideal | `dim_6` | status |
| --- | --- | --- |
| `J(P_2)` | 18 | exact (`fmpz`) |
| width powers of all primitive directions with coefficients `<= 4` | 9 | mod `p`, = exact `R_1 * (a^5,b^5,c^5)` |
| + pullbacks of **all** 2-dimensional integer projections (kernel vectors with entries `<= 4`) | **17** | mod `p` lower bound, generators exact |

The contributing pairs are, e.g., `(y_a, y_b + y_c) ~ (y_a, z)` (2 each), `(y_c, y_b - y_a)` (2 each) and many pairs
with images of 39 or 51 points (1 each). So **8 of the 9 non-width line generators are produced by families of two
integer directions**; only one per line needs the full 3-dimensional family (within line spans, coefficients `<= 4`;
mod-`p` evidence that 17 is the maximum).

Transfer to `Syl(16,2)`: the line image is the full `P_2` (`syl16_base.py`), so these pairs are families of `A` with the
same images. The write-up's exact data (`dim_Q (J_w)_6 <= 225`, `dim_Q (J_w + lines)_6 >= 540`, each line contributing
`<= 9` modulo `J_w`) force the 35 line contributions to be independent and 9-dimensional modulo `J_w`. Hence
**`dim_Q (J_w + all 2-dimensional families inside line spans)_6 >= 225 + 35 * 8 = 505`**, rigorously. The write-up's
"no width power or 2-dimensional family supplies them" (§3.6 item 1, §6.3) is false; it is true only for 2-dimensional
families of **columns and coordinates**.

**(b) `M8`** (`m8_run.py`). All 2-dimensional sublattice families inside the 7 line spans (kernel entries `<= 3`)
give `63, 211, 315` in degrees `6, 7, 8` against `J(A) = 70, 218, 315`: again one generator per line in degree 6 (and 7)
needs the 3-dimensional family. `b_0(M8) = 3` survives within line spans (mod-`p` evidence for the failure part).

**(c) The three-cut example of §3.3** (`three_ex.py`). With `J_w(R=2)` plus pullbacks of **all pairs of primitive
`{-1,0,1}` directions**, the local ideal equals `J(A)` in **every** degree `<= 7` (dims `3, 18, 40, 44, 30` in degrees
`3..7`; `J(A)` exact, local ideal mod `p`, hence equality exact by the sandwich). The degree-4 generators come from the
pairs `(x_4, x_2 + x_3 - x_4)`, `(x_4, x_0 + x_1 + x_2 + x_3 - x_4)`, `(x_2 + x_3 - x_4, x_1 + x_4)`. So for this block
README-`b_0 = 2`, against the write-up's `b_0(col) = None`, `b_0(mixed) = 3`.

**(d) Paley-12, families containing `sigma = -k_base`** (`paley_run2.py`, `paley_run3.py`, `paley_run4.py`).
Invariant part of degree 7, starting from `P_0 L_col(5)_7 = 82`:

| added families (all orbit representatives, mod-`p` lower bounds) | `dim P_0 L_7` |
| --- | --- |
| all pairs from 320 directions (columns, `sigma`, coordinates, `x_a +- x_b`, all width-7 directions) | 82 (from `J_w`: 72 -> 77, via `(y_0, sigma - x_j)`, `j notin D`, image 29 points) |
| all pairs from 210 directions incl. `sigma - x_a - x_b` | 82 |
| all triples from 45 directions (columns, `sigma`, `sigma - x_j`, coordinates, `sigma - y_i`) | 82 |
| all quadruples from 23 directions (columns, `sigma`, `sigma - x_j`) | **87** |
| `sigma` + 3 coordinates (all orbits) | **87** (the source of the 5) |
| `sigma` + 4 coordinates | 87 |

The 5 new invariant dimensions come from the families `{sigma, x_a, x_b, x_c}` (equivalently the base coordinate
`k_base` and three non-base coordinates; image 125 of the 135 points of `[-2,2] x {-1,0,1}^3`, missing exactly 10
sign patterns at `sigma = +-2`; `rho = 8`). These are **4-direction coordinate families**, i.e. exactly the
cut-box mechanism (b) of §6.2, with the base coordinate as one of the coordinates. The write-up's
"columns + coordinates, `b <= 4`: 792 / deficit 725" omits `k_base` and therefore misses them.

## 4. Issues

**Issue 1 (major; scope of every negative claim).** README §7.1's `Loc(b_0)` counts **arbitrary integer directions**.
Every lower bound on `b_0` in the write-up holds only for columns or for columns plus non-base coordinates, yet the
summary, §0 ("false at small `b_0` in an exact model ... not a structural property of slab boxes"), §6.1 ("`Loc(b_0)`
as stated in README §7.1 is refuted as a model-independent principle"), §6.3 and §6.4 read them as statements about
README `Loc`. My computations (§3.2) show that the restricted family sets do overstate `b_0`:
* three-cut example: README-`b_0 = 2` exactly in all degrees `<= 7`, against `b_0(col) = None` and `b_0(mixed) = 3`;
* `Syl(16,2)` lines: two-direction families give at least 280 of the 315 degree-6 line generators (exact lower bound);
* Paley-12: the families `{k_base, x_a, x_b, x_c}` add 5 invariant dimensions beyond `L_col(5)` (82 -> 87) and are not
  in `L_mixed(4)` either, because the write-up's "coordinates" exclude the base coordinate;
* `Syl(16,2)` planes: see Issue 3.

None of this proves README `Loc` in any model: Paley still has `137 - 87 = 50` invariant degree-7 dimensions that no
family I tried produces. But Paley-12 is a counterexample only to **column** locality (`b <= 5`, exact) and, as
evidence, to locality with families drawn from the finite direction sets I and the author tried. It does not refute
README `Loc(b_0)` "as stated", even setting aside that Paley is not a regime box.

**Issue 2 (minor; the degree-7 width set of `Syl(16,2)` is incomplete, but harmlessly).** The write-up's degree-7
width set for `Syl(16,2)` is: 15 columns, `sigma`, 35 line half-sums, 15 plane quarter-sums. My full scan of the
7174453 primitive `{-1,0,1}^15` directions (`syl16_w7scan2.py`; filter on 4751 points of `A`, survivors checked on all
of `A`) finds 2546 directions of width `<= 7`. Besides the listed ones and those whose 7th powers vanish in `R` or lie in
`(y_S^5)` (`1_S +- e_x`: `(y +- x)^7 = y^7 +- 7 y^6 x + 21 y^5 x^2`), there are **105 directions of support 4**: the
indicators `u` of the affine 2-flats `{a.x = 1, b.x = 1}` missing 0, of width 7, with `u^7 != 0` in `R`. They are not in
the write-up's width set, so the row "width powers: 1851" is not `dim J_w(R=1)_7`. However
`u = (y_a + y_b - y_{a+b})/2` lies in a line span (`syl16_affspan.py`), and `((a+b-c)/2)^7` is a width power of the
line image `P_2` (width 7), so `u^7` lies in the line pullback. Every row from "+ lines" on is therefore unaffected.
Confirmed in the Singer-invariant component (`syl16_inv_fix.py`, mod `p`): `J_w + all 35 lines` has dimension 416 there, and the 7 orbits of `u^7` add 0.

**Issue 3 (major; "whole 7-column planes are needed", "`b_0 >= 7`, plausibly 15").** Only column subfamilies of a plane
were tested. The plane span equals the span of the 7 **two-coordinate sums** `kappa_c = x_u + x_{u+v}` over the
nonzero cosets of `V^perp = {0, v}` (the incidence matrix is invertible over `Q`), and it also contains `z_V` and the 7
line half-sums `z_l`. In the plane's own ring `Q[y_1..y_7]` (image of `Syl(16,2)`, 2507 points; `syl16_plane.py`;
generators exact, dimensions mod `p`):

| ideal in `Q[y]_7` | `dim` |
| --- | --- |
| `J(P)_7` | `<= 1002` (mod-`p` upper bound) |
| width powers of all `{-1,0,1}` directions of `y`-space and of the 22 candidates | 421 |
| + the 7 line families | 771 (plane-only part: `1002 - 771 = 231`, matching the write-up's `3465 / 15`) |
| + all non-grid pairs from 22 candidates (`y_a`, `2 z_l`, `4 z_V`, `4 kappa_c`) | 784 (`(y_a, z_V)` and `(z_l, z_V)`, images 29 / 27 points) |
| + all non-grid triples from the same 22 | **918** |

So, mod `p`, **147 of the 231 plane-only generators** per plane come from families of at most 3 integer directions;
4-direction families were not run. (Rigorous part: `dim_Q` of the `<= 3`-direction local ideal is `>= 918`; the
comparison with `771` and `1002` is mod-`p` evidence, like the write-up's own degree-7 numbers.)

Anything generated in the plane ring pulls back to `R`, so these are genuine local generators of `J(A)` from families
of `<= 3` integer directions. The degree-7 evidence therefore does not show column `b_0 >= 7` to be the right measure,
and the regime heuristic "`b_0` at least 7, plausibly 15" (summary; §0 item 5; §6.3(1)) counts **columns of a subspace
family**, not directions, and is unsupported.

**Issue 4 (moderate; three-cut statistics measure the family set, not `Loc`).** The table "`b_0(mixed) <= 4` in every
three-cut block", the 827/27/99 split and the §5 trend table describe `L_col` and `L_mixed`. With pairs of arbitrary
`{-1,0,1}` directions (`three_ex.py`, `three_b4.py`; `J(A)` exact, local ideal mod `p`, so equalities are exact):
* the §3.3 example (`b_0(col) = None`, `b_0(mixed) = 3`): pairs suffice in all degrees `<= 7`;
* the blocks with `b_0(mixed) = 4`, sorted by `(d, rho)`, first 10 of 39 (`three_b4_out*.txt`; e.g. regions
  `[0,0,2,1,1,1,0]`, `[0,0,0,1,1,1,3]`, `[0,1,1,1,1,1,1]`, all `L = (1,1,1)`): pairs suffice in **every** block tested,
  in all degrees `<= rho + 1`.
So the statement "the projection-invisible generators need coordinate directions: they are cut-box generators
supported on a few coordinates" is fine, but README-`b_0` in these blocks is 2, not 3 or 4.

**Issue 5 (moderate; Paley "columns + coordinates" rows of the §3.2 table).** "Columns + coordinates, `b <= 4`: 792"
and "`> 4` mod `p`" exclude `k_base = -sigma`. Including it, `{k_base, x_a, x_b, x_c}` adds generators already at
`b = 4` (Issue 1; exact). The table entries are correct for the stated set but the reading "projections cannot see
them, because every image of `<= 4` columns and `<= 3` columns-or-coordinates is a grid" (§4.2) does not extend to the
base coordinate: `{sigma, x_a, x_b, x_c}` has a non-grid image.

**Issue 6 (minor; labels).**
* "440 non-Koszul syzygies" (§4.2, §6.2) is a mod-`p` count: since mod-`p` ideal dimensions are lower bounds, 440 is an
  **upper** bound for the rational number of syzygies. The write-up labels it "mod `p`", but "Paley has 440 such
  syzygies" in the summary reads as exact.
* "They couple at least 3 slabs on 9-11 of the 11 coordinates" (summary, §4.2, §6.2) has no certificate in the
  write-up; label as evidence.
* "can only occur where the width powers are not Fröberg-generic" (§6.2) is a heuristic reading of the H-basis
  criterion, not a lemma.
* Paley has 1833 point orbits, not 1834 (§4.2).

**What survives unchanged.**
* The sandwich method, Lemmas 4.1-4.3 (the lemma parts), and every **equality** claim I re-ran: `M8` (`J_w + lines`),
  `M8'` (`J(A) = J_w(R=2)` in all degrees `<= 11`, `m8p_run.py`), the two-cut block `(1,2,4)`, `L = (1,1)` (all degrees
  `<= 10`, `twocut.py`), the three-cut example.
* Proposition 4.4 for the stated family sets, with an independent exact `dim_Q P_0 L_col(5)_7 = 82`
  (`paley_exact.py`: exact `fmpz` rank of 4920 integer rows; completeness of the 36 near-full orbits checked by the
  mod-`p` upper bound, the 6 congruence orbits by exact `J(P)`).
* `M8`: `b_0 = 3` for columns and coordinates (exact), and also within line spans for arbitrary 2-dimensional
  sublattices (mod-`p` evidence, `63 < 70`).
* The write-up's critique of README's `rho(pi_F A) < ceil(4w)` restriction: every non-column family I found that
  contributes in degree 7 also has `rho(pi_F A) = 8 = m` (Paley pair image 29 points, `{sigma, x_a, x_b, x_c}` image
  125 points; `paley_pairimg.py`, `paley_quadimg.py`). So this point is supported, not contradicted, by the wider
  family sets.

## 5. Corrected statement

All statements below are at `T = 1`, on the exact toy models named (none is a regime box; all are LP-subcritical; the
Hadamard ones are ball-bound), with `m` the model's order (`m = rho(A)` or as stated) and degrees `i <= m`.

1. **Column / column-plus-non-base-coordinate locality** (family sets `L_col(b)`, `L_mixed(b)` with `J_w(R)`):
   * holds with `b_0 = 1` in `Syl(8,1)`, `M8' = Syl(8,3)` (all `i <= 11`), the 400 completed two-cut blocks with
     `d <= 7` and 827 of the 953 three-cut blocks with `d <= 6` (exact);
   * holds with `b_0 = 3` and fails at `b = 2` in `M8` (exact);
   * **fails in Paley-12** (`L = 2`, `m = 8`) in degree 7 for every `b <= 5` (exact, Proposition 4.4; independently
     re-verified `dim_Q P_0 L_col(5)_7 = 82 < 137 = dim_Q P_0 J(A)_7`); for `b = 6` and for `L_mixed(4)` the failure
     is mod-`p` evidence;
   * in `Syl(16,2)` needs `b >= 3` in degree 6, with the 35 line families sufficient through degree 6 (exact); degree 7
     needs 7-column plane families **among column families** and leaves a residual of 78 (mod-`p` evidence).
2. **README `Loc(b_0)` (arbitrary integer directions) is neither proved nor refuted in any model.** Allowing
   non-column directions lowers the observed `b_0`:
   * three-cut blocks: pairs of `{-1,0,1}` directions give `J(A)` exactly in all degrees `<= rho + 1` in every block
     tested (the §3.3 example, `b_0(mixed) = 3`, and 10 of the 39 blocks with `b_0(mixed) = 4`);
   * `Syl(16,2)`, degree 6: two-direction families inside line spans give at least `225 + 35 * 8 = 505` of the 540
     dimensions (exact lower bound); only one generator per line was not produced by any pair tried;
   * `M8`: the same pairs give `63` of `70` in degree 6 (mod `p`);
   * Paley-12, degree 7: the base-coordinate families `{k_base, x_a, x_b, x_c}` give invariant elements outside
     `L_col(5)_7` (exact: rank 82 -> 84 with two of them; 87 mod `p` with all orbits); 50 invariant dimensions remain
     unexplained by every family I tried (pairs from 320 directions, triples from 89, quadruples from 23, `k_base` with
     up to 5 coordinates, `k_base` with up to 3 coordinates and a column);
   * `Syl(16,2)` planes: in the plane's ring, families of `<= 3` directions (columns, `z_l`, `z_V`, coset sums
     `kappa_c`) give 147 of the 231 plane-only degree-7 generators (mod `p`; `>= 918` is rigorous).
3. Consequently the regime heuristic "`b_0` at least 7, plausibly 15" and the reading "`Loc` is not a structural
   property of slab boxes" are not supported; what is supported is "projection locality onto **columns** (and
   non-base coordinates) is not a structural property". The README §7.3(3) route (projection costs of column families
   bounding `b_0` by `m'/(2j) + O(1)`) is refuted at its first step for column families, as the write-up says. In
   Paley the non-column families found here that contribute below `m` also have `rho(pi_F A) = m = 8`, so the
   `rho < ceil(4w)` restriction would exclude them too; in `Syl(16,2)` (where `m >= rho(A) >= 8`) the contributing
   line pairs, e.g. `(y_a, z)` with a 23-point image, have `rho = 6 < m` and would be admitted.
4. No growth or uniformity statement follows; `M <= (2/3 + eps) log R / loglog R` stands.

## 6. Referee scripts (`referee_loccompute_checks/`)

| script | what | status |
| --- | --- | --- |
| `rlib.py`, `reng.py`, `rengc.py` | exact `J(P)` (`fmpz`); `R` mod `p`, pullbacks, `Z_q` projection | library |
| `paley_base.py`, `paley_inv.py` | Paley `A`, width scan, invariant Hilbert function | exact / mod-`p` upper bound |
| `paley_fam.py`, `paley_run1.py`..`paley_run4.py` | invariant `P_0 L_7` for column, pair, triple, quadruple, `k_base`+coordinate families | mod-`p` lower bounds |
| `paley_exact.py` | exact `dim_Q P_0 L_col(5)_7 = 82`; `+ {sigma, x_a, x_b, x_c}` -> 84 | **exact** |
| `paley_pairimg.py`, `paley_quadimg.py`, `paley_dirs.py` | images, `rho`, congruences | exact |
| `m8.py`, `m8_run.py`, `m8p_run.py` | `M8`, `M8'` | equalities exact; pair deficits mod `p` |
| `t1_line.py`, `t1_line_detail.py` | abstract line image vs 2-dimensional projections | `J(P)` exact; local ideal mod-`p` lower bound |
| `three_ex.py`, `three_b4.py` (`three_b4_out*.txt`) | three-cut blocks with arbitrary pairs | equalities exact |
| `twocut.py` | two-cut block `(1,2,4)`, `L = (1,1)` | exact |
| `syl16_base.py`, `syl16_w7scan2.py`, `syl16_affspan.py` | `Syl(16,2)`: `A`, width-`<= 7` scan, affine flats in line spans | exact |
| `syl16_plane.py` (`syl16_plane7_b3_out.txt`) | plane ring: `<= 3`-direction families | mod `p` |
| `syl16_inv_fix.py` | Singer-invariant degree 7: `J_w + 35 lines = 416`, affine 4-sums add 0 | mod `p` |
| `syl16_inv.py` | Singer-invariant degree 7 | **superseded**: its line orbit representatives used the coordinate permutation on column indices (wrong action), which produced a spurious "+1" for the affine 4-sums; `syl16_inv_fix.py` uses all 35 lines |
