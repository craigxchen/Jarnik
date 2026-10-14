# Referee report on refJ.md

Referee of `round13/refJ.md` (round-13 referee of judge findings J5-J7). All checks use my own code in
`referee_refJ_checks/` (library `rlib.py`). I read the author's `m8lib.py` only to audit how it uses mod-`p` versus
exact arithmetic. The repository was not modified.

## 0 Verdict

**The result survives, with one major error removed.**

| item | refJ claim | referee verdict |
| --- | --- | --- |
| J5 uniform `Q_0`, `2^12..2^20` | proved | **confirmed** (all numbers reproduced, exact R1 for `k <= 20`) |
| J5 per-subspace, `2^12..2^14`, near-extremal `W` | proved | **confirmed**; my thresholds are lower (3.89 / 5.42 / 5.42) |
| Lemma 1.3 `rho(P_j) <= 3j`; Lemma 1.4 Fano blocks | proved | **confirmed** (proof re-derived) |
| J6(a) `TC ⊇ Bs(L_m)`; Prop 2.3 Fano spans in the base locus | proved | **confirmed** (independent exact rank 45 of 45, plus a direct `TC` membership test) |
| J6(b) SLP proof, Thm 3.3 (disjoint cuts) | proved | **confirmed**, with 4 new boxes; **scope: decomposable toy boxes** |
| J6(c) `M8` dims; Fano triples suffice; pairs add nothing; `b_0 = 2` fails on the tested planes | exact | **confirmed** (`(J_w)_{6,7,8} = 616, 1485, 2940` exact); **scope: `M8` is a toy model, not LP-supercritical, `m/j = 4` instead of `~8`** |
| J7(a) hierarchy | proved | **confirmed** (typo: `C'^{-d}`, not `C'^{-2d}`) |
| J7(b) **Prop 5.2(2)**: "MIN^pair fails whenever `max eps_xy >= 2 ln(n/2)`", and the "window of width `< 2/(n-1)`" | proved | **FALSE as stated** (inequality reversed; explicit counterexample with empty `G`). The correct version needs `min eps_xy >= 2 ln(n/2)`, so `abar >= 2 ln(n/2)`, i.e. `W >~ 2 W_min`. The "no escape on non-empty gap regions" conclusion is **withdrawn**; README §4.3(b)'s "open" stands. |

**Issues.**
* **Major:** Prop 5.2(2) and its consequences (§5.2 and the summary row J7(b)).
* **Minor:**
  - the `M8` and disjoint-cut scope is not flagged;
  - the `C'` exponent typo;
  - the step-7 label "None = rho > 6" should read "rho >= 6";
  - "exactly the pencil" is checked only for normals `<= 5`.

**Bearing on Loc(b_0).** Nothing in refJ proves or refutes `Loc(b_0)` for any regime profile. The proved instances
are decomposable toy boxes, and the `M8` evidence comes from a non-regime model. I add one point (§6.1):
XOR-closed blocks of growing dimension are a concrete candidate for non-local generators that would make `b_0`
grow with `n`.

## 1 J5 (refJ §1)

**Verdict: survives.** I re-derived Lemmas 1.1-1.4 and recomputed every table entry with my own script
(`referee_refJ_checks/j5_mine.py`, output `j5_mine_out.txt`: 60-digit `mpmath`, exact `Fraction` collision counts, exact
big-integer R1 for `k <= 20`).

* **Lemma 1.1 (case A).** Correct. A partition with more within-block pairs than `P(L)` cannot refine `pi(L)`.
* **Lemma 1.2 (case B).** Correct, including the closed form of `g*`. The derivative sign is
  `o(1 - ln(c_n Q)) - d/Q - a`, which I re-derived. The side conditions `eta <= 1/4` and `sigma <= 3/10` follow as
  stated. For `sigma <= 3/10` I used `Q_0 abar >= 2 ln(5/3)`, which holds with a wide margin at `Q_0 >= Q_4`.
* **Lemma 1.3 (`rho(P_j) <= 3j`).** The peeling proof is correct. I checked the index bound `index + 2j - t <= 3j - 2m`
  in both cases (`t <= j - 2m` and `t > j - 2m`), including slices with `2m > j`, where only the second case occurs.
* **Lemma 1.4.** Correct. The parity identity holds because `a.x + b.x + (a+b).x = 0 mod 2`, so each `x` lies in 3 or
  1 of the three hyperplanes. Using `sum k = 0` and `|k(S_a)| <= floor(w/2) = j`, the image lies in `P_j`.
  - `rho` is monotone under inclusion and subadditive on products, so `rho(L_{s,c}) <= (3s+2c)j`.
  - `8j < m'` holds because `4w` is not an integer (distance at least `0.07` for `k = 12..21`, recomputed).
* **Uniform-`Q_0` table.** The values `o_1 Q_5 - d` (50526, ..., 37249, -527080) and `g*(Q_5)` (0.357, ..., 9.681) are
  reproduced to all printed digits. R1 is exact `True` for `k <= 20`.
* **Per-subspace table.** Reproduced: (1,1) at `k = 12, 13` has `o Q - d = 2724.7` and `1753.5`, (2,0) at `k = 14` has
  `3010.8`, and the coordinate-span rows 14.9028 and 17.9347 match. With the sharper mixing condition
  `(delta - eta)/(1+eta) >= e^{-abar/2}`, which is the correct one, the thresholds are lower than refJ's "+1" bound:

| `k` | span | refJ threshold | referee threshold |
| --- | --- | --- | --- |
| 12 | (1,1) | 4.89 | **3.89** |
| 12, 13, 14 | (2,0) | 6.42 | **5.4161** (= 2 ln 15) |

  So refJ's numbers are valid upper bounds on the failure level, and the conclusion holds with more margin.
* **Monotonicity in `abar`.** The test points refute the hypothesis at every level at or above the threshold. Raising
  `abar` lowers `eta`, and `o Q > d`, so both conditions only get easier.

**Remarks (minor).**
- (a) The `C`-dependence correction (refJ §1.4.1) is right. Since `C' = max(C,1)`, at `C <= 1` we have
  `abar = W/(2d)`, and case B's threshold moves by `2 ln(1/C)/(o_1 Q - d)`.
- (b) The scope is correct: these are finite-`n` refutations. As refJ itself says, every construction has
  `o N_2/P(L) = O(j) = o(d)`, so nothing asymptotic follows.
- (c) The Sylvester boxes used here are genuine regime boxes: `d >= 316 ln 4d` for `k >= 12`, and R1 is exact. The
  refutation is therefore not a toy-model artefact.
- (d) Twin and duplicate columns do not affect any of this. The `a`'s of `L_{s,c}` must be linearly independent, which
  the lemma requires; otherwise the class count `2^{2s+c}` would be wrong.

## 2 J6(a) (refJ §2)

**Verdict: survives.**

**Proofs re-derived.**
* Lemma 2.0 (duality). The dimension count is consistent.
* Prop 2.1:
  - (i) Leading form of `prod_t (l_v - t)`.
  - (ii) Taylor expansion along `v`: `t`-degree `<= tau` is equivalent to multiplicity `>= m - tau` at `[v]`.
  - (iii) Apolarity in degree `m - tau - 1`.
* Prop 2.2. Every degree-`m` monomial in `s` variables has some exponent `> tau_i` once `sum tau_i < m`.

**Prop 2.3, checked independently** (`referee_refJ_checks/prop23.py`, `prop23_out.txt`).
* Exact widths: `h(p_i) = h(p_i + p_j) = 4` and `h(z_l) = 6`, with `tau = h`.
* The ternary octics with multiplicity `>= 4` at the six points `P_i`, `Q_ij` and `>= 2` at `Z` form a system of
  45 monomials and 70 conditions of rank 45. The system is therefore zero (sympy, exact).
* I re-ran the six peeling sums `12, 11, 10, 7, 6, 5` and the final conic step by hand. They are correct.
* Independent of the base-locus argument, I tested membership directly: for three random integer `u` in `U_l`,
  `l_u^8` lies in `J(A)_8`. The test is that `ev(l_u^8)` lies in `ev(R_{<=7})` on the 1037 points of `A`, by exact
  `fmpz` ranks. A control point outside `U_l` fails.

**Remarks.**
1. The re-reading of "cheap spans = linear base cycles" is a correction of wording, not of mathematics. It is
   correct, and the README sentence should be changed as refJ proposes.
2. Prop 2.1(iii) gives only an inclusion `TC ⊇ Bs(L_m)`. In `M8`, `(J_w)_8 = 2940 < 2961`, so `Bs(L_8)` may be
   strictly smaller than `TC`. refJ says so.
3. The extra base-locus subspaces are exactly the projection-cheap Fano spans. That is a statement about `M8` only.

## 3 J6(b) (refJ §3)

**Verdict: survives.**

**Proofs re-derived.**
* The `sl_2` proof of Theorem 3.1. I checked the commutators on `1, x, x^2` and the two weight cases. In the
  injective case, `mu - 2 < N`, which is what is needed.
* Lemma 3.2. The telescoping gives `sum_{i=s-L}^{s+L} h_i`, and the top degree is strict.
* Theorem 3.3. The two-ideal dimension argument holds; the case `L >= s` is trivial since `l^{2L+1} = 0` in `B`.

**Independent exact checks** (`referee_refJ_checks/thm33.py`, `thm33_out.txt`). These are 4 boxes not in refJ's list.
The Hilbert function of `R/J(A)` comes from exact `fmpz` evaluation ranks, and it equals the predicted product of
Lefschetz truncations in every degree.

| `d` | cuts and signs | `|A|` |
| --- | --- | --- |
| 6 | 3-cut `(+,-,+)`, `L=1` and 2-cut `L=0` | 171 |
| 7 | 4-cut `(+,-,-,+)`, `L=1` and 3-cut `(-,+,+)`, `L=0` | 357 |
| 7 | three disjoint cuts | 171 |
| 6 | full 6-cut, mixed signs, `L=2` | 573 |

**Scope remark (needed).** These are exact statements about **decomposable** boxes. A disjoint-cut box is a product of
one-cut boxes, and such boxes are far from the regime: they are not overdetermined, not indecomposable and not
LP-supercritical. So "`Loc(1)` holds for disjoint cuts" is a proved toy statement, with no bearing on whether `Loc`
holds for regime profiles. refJ §6.1 calls these "proved instances of locality", which is accurate only with this
qualifier.

## 4 J6(c) (refJ §4)

**Verdict: all numbers survive; the scope of the model needs to be stated.**

### 4.1 Independent exact recomputation

Scripts: `referee_refJ_checks/rlib.py`, my own library (no import of `m8lib`), and `m8_a.py` to `m8_d.py`.

| quantity | value | method |
| --- | --- | --- |
| `|A|` | 1037 | enumeration |
| `H(R/J(A))` | `[1,7,28,77,161,259,287,175,42]` | exact `fmpz` cumulative evaluation ranks |
| `dim J(A)_{0..8}` | `[0,0,0,7,49,203,637,1541,2961]` | from `H` |
| width set | 358 primitive directions (7 / 49 / 302 with `h` = 2 / 4 / 6), `tau = h`, `max |v_i| = 2` | enumeration |
| `(J_w)_6` | **616, exact** | `fmpz` rank |
| `(J_w)_7` | **1485, exact** | `fmpz` rank |
| `(J_w)_8` | **2940, exact** | `fmpz` rank, 9464 x 3003 (§4.3) |

* **Width set.** The completeness argument (`M >= 4` forces 9 distinct values) is correct.
* **Images.** Re-derived exactly:
  - column pairs: full grid, 25 points, `rho = 8`;
  - non-line triples: full grid, 125 points, `rho = 12`;
  - Fano lines: 57 points, the parity box minus exactly the six mixed-sign points `(±2,±2,±2)`, `rho = 6`;
  - `dim J(Y)_{5..8} = 6, 21, 36, 45` for the Fano lines.
* **Locality.**
  - `J_w` plus the 7 Fano pullbacks reaches `[203,637,1541,2961]` in degrees 5..8.
  - Adding the pairs or the non-line triples leaves `[203,616,1485,2940]`.
* **Rigour of the sandwich.** The pullback generators are exact rational forms in `J(A)`, computed by exact `fmpq`
  kernels. The mod-`p` evaluation rank is at most their `Q`-span dimension, which is at most `dim J(A)` (exact). So
  equality is rigorous.
* **Pairs and non-line triples add nothing over `Q`.** This needs the grid argument, which is correct: products of
  1-D sets, so `J(Y) = (y_i^5)`.
* **No misuse of mod-`p` arithmetic.** I read the author's `m8lib.py`. The upper bounds (616, 630 and the like) use
  `fmpz` ranks, and the lower bounds use mod-`p` ranks.

### 4.2 Minor defects

* `m8_step7_out.txt` prints "rho histogram ... None = rho > 6". `proj_ideal_gens(..., 6)` detects `rho` only when
  `rho <= 5`, so "None" means `rho >= 6`. The planes `span(z_l, c_a)` have `rho = 6` (step 9) and are among them.
  Harmless, but the label is wrong.
* "The contributing planes inside `U_l` are exactly the pencil through `z_l`" is verified only for normals with entries
  `<= 5`.

### 4.3 Adversarial search for non-local planes, and `(J_w)_8`

**Exact `(J_w)_8`.** `(J_w)_8 = 2940` exactly over `Q` (`m8_jw8_out.txt`: 9464 x 3003 `fmpz` rank, 29 s). So refJ's
dims of `J_w` are confirmed exactly in every degree `<= 8`.

**Plane search** (`m8_planes.py`, `m8_planes_out.txt`). I looked for planes outside refJ's family (refJ used planes
spanned by two of the 358 low-width directions).
* Rows `P_x in Z^2` were drawn from small random palettes, so the planes have coincident and lattice-related rows.
  There were 75,545 distinct planes.
* Of these, 4,038 pass the necessary admissibility filter `|pi_P(A)| <= 36`.
* Each pullback ideal (all degrees `<= 6`, exact `fmpq` generators) was added to `J_w` plus the Fano-pencil planes
  `span(z_l, c)`, whose mod-`p` dimension is 630, matching refJ's exact 630.
* **No plane raised the mod-`p` dimension above 630.**

This is **evidence only**. An increase mod `p` would have been rigorous; the absence of an increase is not. It
supports refJ's "`Loc(2)` fails in `M8`, evidence in general". I found no non-local plane that refJ missed.

### 4.4 Scope of `M8` as a model (not stated in refJ)

* `M8` is **not** LP-supercritical: the smallest LP-supercritical Sylvester model has `n >= 32` (zetc). It is not a
  regime box.
* Its ratio of order to slab is `m/j = 8/2 = 4`. Regime boxes have `m'/j ≈ 8`, since `m' = ceil(4w)` and
  `j = floor(w/2)`.
* So in `M8` only single Fano blocks (`3j = 6 < 8`) are admissible, while in the regime products of two Fano blocks and
  a column (`8j < m'`) and possibly whole XOR-closed blocks are admissible.
* Hence "Loc(3) holds in `M8`" is a toy statement. It gives a *lower-bound* heuristic (`b_0 >= 3` is plausible for
  Sylvester-type profiles). It says nothing about whether `b_0` is bounded in the regime. refJ §4.5 and §6.3 are
  worded as suggestions, which is acceptable; the summary table should add "toy model, `m/j = 4`, not
  LP-supercritical".

## 5 J7 (refJ §5)

### 5.1 Prop 5.1 (hierarchy): survives, with one typo

**Re-derived.** Using `Gamma = (2/N_2) sum ln^+(1/|u_x-u_y|)` (nearstrata l.105):
* `G ⊂ Gam`: on `G`, `Gamma <= (1/N_2) sum eps_xy <= abar`.
* `R_per ⊂ R_unif(Q_0)`: on the subfamily `P(pi) >= max(1, N_2/Q_0)`, `sigma_pi >= (1/2)e^{-Q_0 abar/2}`.

**Typo in `MIN* => (MIN)`.** `e^{-d abar/2} = e^{-W/4} C'^{-d}`, not `C'^{-2d}`. With `C' = max(C,1)` the condition is
`C C'^d <= 1`, which holds iff `C <= 1`, so the conclusion is unchanged.

**Convention `P(pi) >= 1`.** refJ is right that it is needed to make the uniform family a subfamily.

### 5.2 Prop 5.2 (dichotomy): part (2) is FALSE AS PROVED. This is a major error.

The proof of (2) says: "every pair has `|u_x - u_y| >= 2/n >= e^{-max eps/2} >= e^{-eps_xy/2}`, so `u in G`".

**The last inequality is reversed.** Since `eps_xy <= max eps`, we have `e^{-eps_xy/2} >= e^{-max eps/2}`. The gap
region `G` (README §2.6, nearstrata l.437) puts a **per-pair** constraint `|u_x - u_y| >= e^{-eps_xy/2}`. The binding
pairs are those with **small** `eps_xy`, so `max eps` is the wrong quantity.

**Counterexample to the statement (2) itself.**
* Fix four points `x_1..x_4` and set `eps_{x_i x_j} = 0` for `i < j <= 4`. Set every other `eps_xy = E` with
  `E >= 2 ln(n/2)`.
* Then `max eps = E >= 2 ln(n/2)`.
* Membership in `G` needs four reals in `[-1,1]` with pairwise distances `>= e^0 = 1`. That is impossible: the
  extreme two would be `>= 3` apart. So `G` is empty and MIN^pair holds **vacuously** (empty minimum `= +infinity`).
* This contradicts "(2) MIN^pair fails".
* The data are also compatible with the cluster budget `sum eps_xy <= N_2 abar` as soon as `abar >= E`.

So Prop 5.2(2) is false as stated, for every `k`.

**Counter-data to the proof alone** (`referee_refJ_checks/prop52.py`, `prop52_out.txt`), where `G` is non-empty:
* Take `eps_xy = 2 ln(n/2) + 1` for all pairs except one pair of adjacent values of `u`, given `eps = 0`.
* Then `max eps >= 2 ln(n/2)`, but that pair needs gap `>= 1` and has gap `2/n`. So `u notin G`, and the proof gives
  nothing.
* `G` is non-empty for these data. For this single pair the failure can be rescued by an affine relabelling of
  `F_2^k`, since all columns are present. With several small `eps_xy`, no argument is given.

**Corrected statement (proved by the same argument).**
> (2') If `min_{x<y} eps_xy >= 2 ln(n/2)`, then at every order in `(2jk, m_max]` MIN^pair fails.

Part (1) is correct as stated: if `max eps < 2 ln((n-1)/2)` then `G` is empty, and that direction uses `max`
correctly.

**Consequence: the "window" conclusion collapses.**
* (2') says nothing about data in the region `min eps < 2 ln(n/2) <= max eps`.
* Every actual cluster has `sum eps_xy <= N_2 abar`, hence `min eps_xy <= abar`. So (2') needs `abar >= 2 ln(n/2)`.
* At `k = 16..20`, `2 ln(n/2) = 20.79..26.34`, while the near-extremal level is `abar = 9.64..13.79`. (2') therefore
  needs `W >= 1.91-2.16 W_min` (`prop52_out.txt`).

So the claims below are **unproved**:
* "MIN^pair fails at every order in `(2jk, m_max]` for every actual cluster unless `max eps` falls in a window of width
  `< 2/(n-1)`";
* "MIN^pair gives no escape on non-empty gap regions";
* the proposed sharpening of README §4.3(b).

The README's original wording, "whether it is weaker on non-empty gap regions is open", stands. It is now settled
only for data with `min eps >= 2 ln(n/2)`, which actual near-extremal clusters cannot have.

### 5.3 Empty-gap witness: survives

`(n-1) e^{-max eps/2} = 104.94 > 2`, and the argument correctly uses `max`. Confirmed.

## 6 Other checks of refJ §6, and the corrected statement

### 6.1 Checks of refJ §6

* **§6.4 (strata half).** The deduction is correct.
  - A `t`-dimensional XOR-closed block has `2^t - 1` columns. Its span is all functions of `t` bits, with `2^t` equal
    classes.
  - Its image lies in `Y_t(j)`, so `o >= m' - rho_t`.
  - Lemma 1.2 with per-subspace radii then fails the hypothesis once `(m' - rho_t) Q_t > d`, with `Q_t ~ 2^t`. The
    mixing threshold is only `2 ln(2^t - 1) + o(1)`.
  - The `rho_3` bounds come from mod-`p` evaluation ranks. Full mod-`p` rank implies full `Q` rank, so these are
    rigorous upper bounds on `rho`, which is the direction refJ uses.
* **A possible non-local mechanism refJ does not draw out.** If `rho_t <= c j` with `c < 8` for `t` growing with `n`,
  then whole `t`-dimensional XOR blocks (`2^t - 1` columns) are admissible families. Their pullbacks are candidates
  for generators of `J(A)` in degree `<= m'` that need `b_0 >= 2^t - 1`, i.e. **`Loc(b_0)` with `b_0` independent of
  `n` could fail through XOR blocks of growing dimension.**
  - At `j = 2` refJ's own data give `rho_3 <= 8 = 4j < 8j ~ m'`, so 7-column blocks are admissible in the regime.
  - Whether they contribute non-local generators in degree `<= m'` is not computed. This is the same open question as
    refJ §6.4, and it bears on `Loc` directly, not only on the strata half.

### 6.2 Corrected statement

All items are at `T = 1`, Sylvester profile `n = 2^k`, `r = d = n-1`, `w = sqrt(8d ln 4d)`, `j = floor(w/2)`,
`m' = ceil(4w)`, `C = 1`, and `abar = W/(2d)`.

1. **(J5, proved.)**
   - ZE-TC^coll_eps`(A, m'; Q_0)` is false for every `Q_0 > 0`, for `12 <= k <= 20`, at every
     `abar >= max(2 ln 15, g*(Q_5))`. This includes all `W >= 3n ln n - 14n`.
   - ZE-TC^coll_sigma is false for `k = 12, 13, 14` at every `abar >= 3.89 / 5.42 / 5.42` respectively, by the
     sharper mixing bound; refJ's 4.89 / 6.20 / 6.42 are also valid. This includes all near-extremal `W`.
   - For every `k >= 5`, `L_{2,1}` (32 classes) lies in `TC(A, m')`.
   - All of these are finite-`n` statements; nothing asymptotic follows.
2. **(J6(a), proved.)** `TC_C(A,m) ⊇ Bs_C(L_m)` and `TC(A,m) ⊇ Bs_R(L_m)`, where `L_m` has multiplicities
   `m - tau_A(v)` at the integer `[v]` with `tau_A(v) < m`. Cheap spans are the positive-multiplicity linear cycles
   `L_I`. In `M8` (a toy model) the base locus also contains the 3-dimensional Fano spans.
3. **(J6(b), proved.)** For pairwise disjoint signed cuts on `{-1,0,1}^d`, `J(A) = (x_i^3, l_j^{2L_j+1})` in all
   degrees. Hence `Loc(1)` holds there. These boxes are decomposable toy boxes, not regime boxes.
4. **(J6(c), exact, toy model `M8`: `m/j = 4`, not LP-supercritical.)**
   - `J(A) = J_w + sum_{7 Fano lines} pi_F^* J(pi_F(A))` in all degrees `<= 8`.
   - `dim (J_w)_6 = 616 < 637`.
   - Column pairs and non-line triples add nothing over `Q`.
   - `Loc(2)` fails for column pairs and for the tested planes (630 < 637 in degree 6). For all planes it is open
     (evidence: fails; see §4.3).
5. **(J7(a), proved.)** The hierarchy holds with the convention `P(pi) >= 1`. `MIN* => (MIN)` needs `C <= 1`, and the
   exponent is `C'^{-d}`.
6. **(J7(b), corrected.)** The empty-gap witness is confirmed. Prop 5.2 must read:
   - (1) if `max eps_xy < 2 ln((n-1)/2)`, then `G` is empty;
   - (2') if `min eps_xy >= 2 ln(n/2)`, then MIN^pair fails at every order in `(2jk, m_max]`, `k >= 16`.

   (2) with `max` is **false**: four points with mutual `eps = 0` give an empty `G`. (2') requires
   `abar >= 2 ln(n/2)`, i.e. `W >~ 2 W_min`. So at near-extremal `W` whether MIN^pair is weaker on non-empty gap
   regions **remains open**, as README §4.3(b) said.

**No growth improvement and no uniform bound follow.**
