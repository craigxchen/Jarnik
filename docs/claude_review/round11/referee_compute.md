# Referee report on round11/compute.md (adversarial)

Referee: adversarial review of `round11/compute.md` (2026-10-13). I re-derived every proof and recomputed every
exact input with my own code. Scripts and outputs are in `round11/referee_compute_checks/`.

## 0. Verdict (summary)

**The exact results survive. Two interpretive claims do not.**

The status "evidence, no growth bound, no proof of ZE" is correct.

**What I confirmed independently** (own kernel, own factorisation, flint Groebner bases instead of sympy):
* Prop 2.1 (zonotope certificate), Prop 2.2 (no such certificate in the regime), Prop 3.1 (Hadamard bound) and
  Prop 2.3 are correct.
* Prop 2.4 is correct after a repair (§2.2). I also prove its conclusion for *every* slab-box at `n - 1 = 2`.
* The `n - 1 = 3` scan reproduces exactly: 355 distinct admissible sets, all with `zeta = rho`.
* **Theorem 5.2 (both models) is correct.** The saturations reach the unit ideal on my own basis, with a different
  engine.
* **Prop 6.2 is correct.** I rebuilt the skew double core from the Paley construction and recomputed every exact
  input: `dim V = 1`, `c` even, `Phi = 8! g` with `g` irreducible, a sign change, and `zeta <= 7` by a different
  rigorous argument.

**What fails.**
1. **(Major) The recommended next target ZE-braid is false in the regime** (§2.8).
   * For Sylvester-Hadamard profiles of order `n = 2^k >= 2^16`, Theorem 2.1's slab-box contains, inside `Bs(A,m)`,
     a `k`-dimensional rational subtorus whose generic points have no collisions. This follows from refereed
     `determinant.md` Prop 5.2.
   * The braid-arrangement structure seen at `n = 8` is a small-`m/w` artefact.
2. **(Major) The `n = 8` models are not "faithful to the regime" in the sense that matters** (§2.9).
   * `zeta < m` is only a necessary condition.
   * `M8` and `M8'` violate the README's own model condition (ii), `h_A < m` (`8 = 8`, `12 > 11`), and miss the
     regime's LP margin `h_A <= m/2` by a factor `>= 2`.
   * So Theorem 5.2 is a correct theorem about two LP-subcritical toy models. It is weak evidence for ZE in the
     regime.
3. **(Minor)**
   * Quantifier and scope slips: "ZE holds" means the model statement for ad hoc `(T, L, m)`, not ZE(8), which is
     vacuous.
   * Prop 6.1 needs `A' ⊂ B_T`.
   * The "`m <= 2T` unfaithful" remark needs `L >= T`.
   * Labels `{1..8}` versus `0..7`.
   * The basis-dependent "six" relations.
   * The unproved "nondegenerate real zeros of each single `P_c`".
   * Overgeneralised scope statements.
   See §4.

The corrected statement is in §5.

## 1. What was checked and how

I re-derived every proposition of `compute.md` (§§1-7) and re-ran every exact claim that a theorem rests on with
my own code. I did not import the author's scripts or data. I read two of the author's scripts only to see what
they compute (`cone_sat.py`, `cone_reduce.py`; see §2.6).

My scripts are in `referee_compute_checks/`. Each writes `<name>_out.txt`. They use python-flint 0.9.0
(`fmpz_mat` exact rank and nullspace, LLL, `fmpz_mpoly.factor`, `fmpz_mpoly_vec.buchberger_naive` over `Q`) and
exact Python integers. The saturations use a different Groebner engine from the author's (flint, not sympy), on
a basis I computed myself.

| script | what it checks | status |
| --- | --- | --- |
| `rmodels.py`, `r_had8.py` | Sylvester `H_8` slab-boxes, `T = 1`, `L = 1, 2, 3`: `|A|`, ranges, `rho`, `dim V(A,m)` | exact rank mod `2^61-1`; exact over `Z` in `r_kernel.py` |
| `r_kernel.py` | exact integer kernel of the moment matrix for `M8`, `M8'`; LLL; `M c = 0` checked exactly | exact |
| `r_factor.py` | factorisation over `Q` of every basis polynomial; two-term factor shape; "form in `x = w-1`" test | exact |
| `r_saturate.py` | Lemma 5.1 / Theorem 5.2: Bayer saturations of `(G_j)` by the 28 collision forms (`G_j` = product of the non-binomial factors) | exact (flint Groebner over `Q`) |
| `r_tc_m8.py` | §4 tangent-cone claims for `M8` | exact |
| `r_sdc7.py`, `r_sdc7_c.py` | skew double core `q = 7` from scratch (Paley); `|A|`, `rho`, `dim V`; Prop 6.2 inputs; rigorous `zeta <= 7` | exact |
| `rzono.py`, `r_scan.py`, `r_scan_triples.py` | all sign profiles at `n = 3` (`T <= 8`, `L <= 7`) and `n = 4` (`T <= 4`, `L <= 3`): admissible models and zonotope witnesses for `zeta = rho` | exact (witness search; a witness certifies `zeta >= rho`) |
| `r_scan5.py`, `r_n5_hard.py` | my own random sample at `n = 5` (seed 7, 12 profiles, `T, L <= 3`); the two instances the search did not settle are decided by factorising an exact basis | exact |
| `r_had8_T2.py`, `r_had8_rho.py` | `H_8` with `T = 2`: is condition (ii) reachable at `n = 8`? | rank mod `p` (upper bound for `rho`) |
| `r_single_real.py` | single-form real zeros (§2.10) | floating point, evidence only |
| `r_braid_threshold.py` | exact `m` (big-integer binomials) and slab widths for Sylvester profiles `n = 2^k`; failure of ZE-braid in the regime (§2.8) | exact `m`; `w` in floating point with a large margin |

Notes are in `referee_compute_checks/notes.md`.

## 2. Section-by-section re-derivation

### 2.1 §1: reformulations R1-R3

All three are correct.
* **R1.** `V(A,m)` is the annihilator of `Poly_{<m}|_A`, so `w` lies in `Bs(A,m)` iff `k -> w^k` agrees with a
  polynomial of degree `< m` on `A`.
* **R2.**
  * `dim V(A,m) = |A| - sum_{i<m} dim (S/J)_i = sum_{i>=m} dim (S/J)_i`.
  * `S/J` is standard graded, so its nonzero degrees form an initial segment, and `V(A,m) != 0` iff `m <= rho(A)`.
  * The description of `TC(A,m)` is the annihilator argument given in the text.
* **R3.** Every slab-box is symmetric (`A = -A`).

### 2.2 §2: Props 2.1-2.4

* **Prop 2.1 (zonotope certificate): correct.**
  * `sum_k c_k P(k) = (prod Delta_i^{alpha_i} Q)(0)` with `deg Q < m = |alpha|`, so the sum is 0.
  * `|c|_1 <= 2^m`, and `H = (U T^{m-1})^2 >= 2^m` for `T >= 2`, `m >= 2`. For `T = 1`, `H = U^2 >= 2^m` holds in
    every model used (`2^8 <= 1037^2`, `2^11 <= 1877^2`).
  * The conclusion "`w^{v_i} = 1` for some `i` with `alpha_i > 0`" is degenerate in the ZE sense near `e`.
* **Prop 2.2 (no zonotope with `sum L_i >= m` in the regime): correct.** I re-derived the chain.
  * Widths add on zonotopes. With `sum omega_p <= 1`, `2w <= m/2` and `omega_p |<v,beta_p>| = (2 log p/W)|<v,a_p>|`,
    this gives `sum_i L_i X_{v_i} <= mW/4`.
  * So some `X_{v_i} <= W/4`.
  * The corrected character inequality is `X_v >= W/2 - 2 log(C|v|_1/(2 sqrt 2))`. It follows from
    `2 arcsin(e^{-X/2}/sqrt 2) <= |<v,delta>| <= |v|_1 C e^{-W/4}/2`.
  * The lift of `v_i` to `Z^n_0` has `|v_i|_1 <= 2 sum_{x<n} |v_{i,x}| <= 4n(n-1) = 2D`.
  * Together these give `W <= 8 log(CD/sqrt 2) <= 8 log^+ C + 16 ln n + 4 ln 2`, which contradicts (4.1).
* **The sentence after Prop 2.2 is not correct as stated (minor).** It reads: "every model with `m <= 2T` (the
  coordinate segments have `2T+1` points) is unfaithful".
  * A coordinate segment `{-T..T} e_x` changes `k(S)` by `2T` for every cut `S` that separates `x` from the base point.
  * So the segment lies in `A` only if `L >= T`, or if no cut separates `x` from the base point.
  * In the regime `L = w/2 >> T`, so the remark is harmless there, but it is false for general models.
* **Prop 2.3 (section principle): correct.**
* **Prop 2.4 (`n - 1 = 2`): the conclusion is correct, but the proof covers only condition (ii).**
  * The case "two or more cutting patterns" uses (ii), `h_A(beta_p) < m`, to get `m >= 4L + 1`.
  * The faithful `n = 8` models of §3.1, and the scans at `n - 1 = 3, 4`, use the weaker (ii'), `range < m`.
  * Under (ii') the written proof breaks: `m > 2L` only, against `rho <= 4L`.
  * *Repair (referee).* No condition of type (ii) is needed. Every slab-box at `n = 3` is the set of lattice points
    of `P = {|k_1| <= a, |k_2| <= b, |k_1 + k_2| <= s}` with integers `a, b, s`.
    * If `s >= a + b`, `P` is a rectangle. If `s <= |a - b|`, `P` is a unimodular parallelogram. In both cases two
      independent functionals give a complete intersection, so `rho <= zeta`.
    * Otherwise `P` is a lattice hexagon. Its edges have lattice lengths `s - a + b` (direction `(0,1)`),
      `a + b - s` (direction `(1,-1)`) and `s + a - b` (direction `(1,0)`). So `P` is a lattice zonotope with
      `zeta >= a + b + s`.
    * The author's line-by-line Newton interpolation, with line `sigma` carrying at most `a + b - |sigma| + 1`
      points, gives `rho <= a + b + s`.
    * Hence `zeta = rho` for **every** slab-box at `n - 1 = 2`, whatever condition (ii)/(ii') is imposed.
  * Exact confirmation (`r_scan.py 3 8 7`): 55 distinct admissible sets under (ii') with `T <= 8`, `L <= 7`, all with
    a zonotope witness of total `rho`. None is admissible under (ii), which agrees with the written proof.

### 2.3 §3: scans at `n - 1 = 3, 4`

* **`n - 1 = 3`: reproduced.** `r_scan.py 4 4 3` enumerates all 127 sets of the 7 non-constant patterns, keeps the
  full-rank ones, and takes `T <= 4`, `L <= 3` with conditions (i), (ii'), (iii), (v).
  * It finds **355 distinct admissible `A`**, the author's number.
  * For **all 355** my own search finds a zonotope with `sum L_i = rho`, so `zeta = rho` (exact).
  * 216 of them are admissible also under (ii).
  * Counting admissible triples `(profile, T, L)` (`r_scan_triples.py`) gives **980** with (ii') and **339** with
    (ii): the author's numbers exactly. Counting every admissible `m` as well gives 3664.
* **`n - 1 = 4`.** The author's 12 profiles cannot be reproduced without the seeds. I drew my own 12 full-rank
  random profiles (`r_scan5.py`, seed 7). All 75 distinct admissible sets have `zeta = rho` (§3).
* **The general claim "the README's model sizes are vacuous" is right in substance but overstated in form.**
  * `n - 1 = 2`: proved (with the repair above).
  * `n - 1 = 3`: exhaustive only for `T <= 4`, `L <= 3`, and only for one common `L` for all columns.
  * `n - 1 = 4`: a random sample.
  * The README asked for these sizes under (ii), and under (ii) fewer models qualify. The conclusion should say
    "in every case computed".
* **§3.1 Prop 3.1 (Hadamard): correct.**
  * `sum_p <v,h_p>^2 = n|v|_2^2`, `<v,h_1> = 0` and `|<v,h_p>| <= |v|_1 <= |v|_2^2`.
  * Hence `sum_{p>=2} |<v,h_p>| >= n`, which gives `n sum L_i <= (n-1) max h_A`.
* **The `M8`/`M8'` data are reproduced exactly.**
  * `|A| = 1037`: `rank Poly_{<8} = 995`, `rank Poly_{<9} = 1037`, `dim V(A,8) = 42`.
  * `|A| = 1877`: `rank Poly_{<11} = 1870`, `rank Poly_{<12} = 1877`, `dim V(A,11) = 7`.
  * My LLL bases have sup-norms `6` and `4`, so `V_H = V`.
  * At `L = 1` the model is `|A| = 57`, `rho = 3`, which explains why the table starts at `L = 2`.
* **§3.2 (the hard `n = 5` instance).** The argument "a product of a monomial and two-term factors of total order 13
  supported in `A` is a zonotope certificate" is correct.
  * The vertices of a Newton polytope carry nonzero coefficients, so the Minkowski-sum zonotope lies in
    `conv(supp) ⊂ conv(A)`.
  * `A` is the set of all lattice points of its polytope.
  * I did not re-run this single instance.

### 2.4 §4: tangent cones

* **Unfaithful models.** `Phi_c = m! prod <v_i,u>^{alpha_i}` follows from comparing lowest-order terms. Correct.
* **`M8`: every exact claim is reproduced with my own basis** (`r_tc_m8.py`).
  * All exponents are `<= 2`, and every monomial involves at least 4 variables.
  * All 63 coordinate subspaces of dimension `<= 3` lie in `TC`.
  * Exactly 7 of the 35 coordinate 4-planes lie in `TC`, and they are exactly the 7 cuts avoiding the base point.
  * No coordinate subspace of dimension `>= 5` lies in `TC`.
  * (My maximal coefficient is 161280 against the author's 120960. This depends on the basis and is harmless.)

### 2.5 §5.1: Lemma 5.1 (cone reduction)

Correct, with one wording repair (minor).
* The two-term factors in the dichotomy must be those of the `P_j` with `j in J`.
* The lemma also needs the irreducible two-term factors to be `w^a - w^{a'}` with `a - a'` primitive.
  * Irreducible factors over `Q` automatically have `a - a'` primitive.
  * A factor `w^a + w^{a'}` would give a torsion translate, not a subtorus through `e`.
  * In both models every two-term factor is of the form `w^a - w^{a'}` (checked, §3).
* Bayer's criterion (grevlex, saturating variable last, homogeneous ideal) is applied correctly. It holds for any
  Groebner basis, because for homogeneous `f` in grevlex `x_r | f` iff `x_r | in(f)`.

### 2.6 §5.2: Theorem 5.2, independently re-proved

The theorem survives an independent recomputation. I used my own kernel, my own factorisation and a different
Groebner engine.

* **`M8'` (`L = 3`, `m = 11`)** (`r_kernel.py 3 11`, `r_factor.py`, `r_saturate.py`).
  * My LLL basis of `V(A_3,11)` (`dim 7`, verified `Mc = 0` exactly) has the structure the author reports. Each
    basis polynomial is a monomial, times two-term factors that are all linear in `x = w - 1` (`w_i - 1`, `w_0 - w_1`,
    `w_3 - w_4`), times **one** irreducible form in `x`. There are three cubics with 6 `x`-terms and four quartics
    with 12, and the orders at `e` add up to 11.
  * Example: `x2x3x4 + x2x3x5 - x2x4x5 - x3x4x5 - x2x3x6 + x4x5x6`, identical to the author's.
  * The ideal of the 7 forms has a 30-element Groebner basis (as reported).
  * Saturating successively by `x_0, ..., x_6, x_0-x_1, ..., x_2-x_3` gives the unit ideal, in 2 s with flint.
* **`M8` (`L = 2`, `m = 8`)** (`r_kernel.py 2 8`).
  * My basis also splits 24 / 18: 24 polynomials whose non-binomial factors are all forms in `x`, and 18 with a
    non-homogeneous factor.
  * There are 27 distinct two-term factors, all of the form `w^a - w^{a'}` with `a - a'` primitive.
  * **Seven** of them are not linear in `x` in my basis: `w0w3 - w4`, `w0w5 - w2`, `w0w6 - w1`, `w1w3 - w5`,
    `w1w6 - w0`, `w2w3 - w6`, `w2w5 - w0`. The author has six. The list of non-collision subtori produced by Lemma 5.1
    depends on the basis (the true `Bs` does not), so "six four-term relations" is a property of the author's
    basis, not of the model.
  * With `G_j` the **product** of the non-binomial factors of `P_j` (as Lemma 5.1 requires), the ideal of the 24
    `G_j` saturates to the unit ideal after the same 28 forms.
  * I checked that the author's `cone_reduce.py` also takes products. The danger was that the 11 polynomials with
    two cubic factors would be entered as separate generators, which would compute a smaller variety; this did not
    happen.
* **Verdict on Theorem 5.2: correct as an exact computation.**
  * Every point of `Bs(A_3, 11)` in `(C^*)^7` satisfies `w_x = 1` or `w_x = w_y`.
  * Every point of `Bs(A_2, 8)` lies on `{w_x = 1}`, on `{w_x = w_y}`, or on one of finitely many explicit
    subtori `{w^a = w^{a'}}` through `e`.
* **Quantifiers of "ZE holds in M8 and M8'" (minor, must be stated).** The theorem is not ZE(8), and not ZE for any
  `(T, w, m, H)` of Theorem 2.1.
  * ZE(8) is vacuous: `7 >= 300 ln(4r)` fails for every `r >= 1`.
  * With Theorem 2.1's own parameters at `n = 8` (`T = 8`, `r = 7`), `w = 109.3`. But `max_{B_T} |<k,beta_p>| = 64`, so
    **the slabs do not cut**, `A = B_8`, `m = 39 < 4w`, and ZE holds by referee Prop R1 of round 10.
  * What Theorem 5.2 proves is the model statement "`Bs(A,m)` contains no nondegenerate point" for one profile
    (Sylvester `H_8`, base point 8) and ad hoc parameters:
    * `T = 1`, much smaller than `n`;
    * slab `|k(S)| <= L`, i.e. `|<k,beta>| <= 2L` (`= 4` or `6`);
    * `m = rho(A)`;
    * `V_H = V(A,m)`;
    * all `delta` on the whole torus.
  * The `delta` quantifier is stronger than ZE needs. The profile and parameter quantifiers are much weaker.

### 2.7 §6: Prop 6.1 and the skew double core (Prop 6.2)

* **Prop 6.1: correct up to a missing hypothesis (minor).** The proof is "word for word Prop 2.2", which bounds
  `|v_i|_1 <= 2D` using the box `B_T`. For an arbitrary finite `A' ⊂ Z^d` one needs either `A' ⊂ B_T`, or `log diam(A')`
  in place of `log D` in the condition on `W`. Every support the author has in mind lies in a box, so this is a
  statement fix only.
* **Prop 6.2: all exact inputs reproduced from scratch** (`r_sdc7.py`, `r_sdc7_c.py`).
  * Construction. I built the skew-Hadamard matrix of order 8 from the Paley construction at `q = 7` and normalised
    it as in `round5/exact.md` Thm 4.1. The 14 cuts are distinct and balanced, and every pair of points is separated
    by exactly 8 of them (Plotkin-tight).
  * With `T = 1`, `L = 2`: `|A| = 757` for **every** choice of base point.
  * Exact ranks: `rank Poly_{<8} = 756` and `rank Poly_{<9} = 757`, so `dim V(A,8) = 1` and `rho = 8`.
  * The primitive `c` is even, with sup-norm 42 and support 687, so `V_H = V`.
  * `Phi = 8! g`, where `g` has 196 terms with coefficients in `{-2,-1,1,3}` and is **irreducible over `Q`** (flint).
  * In my coordinates (base point 8), `g` is negative at all 65604 nonzero-valued points of `[-2,2]^7` (12520
    zeros). A positive value is `g(25, 6, -12, 1, -11, -14, -20) = 76047300`.
  * Faithfulness, rigorously: every zonotope generator `v` (with `|v_x| <= 2` off the base) has
    `sum_p |v(S_p)| >= 8`. The minimum over all `5^7 - 1` candidates is exactly 8. The 14 slab functionals have range
    4, so `8 zeta <= 14 * 4` and `zeta <= 7 < 8 = m`.
* **The proof of Prop 6.2 is correct.** The remainder bound for `cos` is
  `|cos x - sum_{j<=4} (-1)^j x^{2j}/(2j)!| <= x^10/10!`. The segments parallel to `e = u_- - u_+` project to single
  points of `e^perp`. A linear map between 6-dimensional spaces sends null sets to null sets. And
  `F|_{H_v} == 0 => <v,u> | Phi`, which irreducibility excludes.
  * The description "a real **hypersurface** through `e`" is informal, since only nondegenerate zeros
    `O(t)`-close to `e` are proved. This is harmless.
* **Interpretation (fair).** The model satisfies neither condition (ii) nor the regime's `dim V_H >= |A|/2`, and
  the author says so. It shows that `zeta < m` alone does not imply ZE-real.

### 2.8 §5.3: ZE-braid is false in the regime (major)

The write-up's recommended next target is **ZE-braid(A, m)**: *the degree-`m` forms `H` whose difference operators
`H(Delta)` fit, after translation, inside `A` have no common zero off the braid arrangement
`{x_x = 0} ∪ {x_x = x_y}`*.

The implication ZE-braid => ZE is correct, and needs no Lemma 5.1.
* If `w` is in `Bs(A,m)` and `H` is realisable via `tau_a H(Delta)`, then `c = tau_a H(Delta)` is in `V(A,m)`, because
  `w^a H(w-1)` vanishes to order `m` at `e`.
* Hence `0 = w^a H(w-1)`, so `w - 1` is in `V_C(Sigma_A)`.

**But ZE-braid is false for the regime slab-box of the Sylvester-Hadamard profiles, so it cannot be the target.**

*Proof (referee).* Let `n = 2^k` and let `P` be the Sylvester-Hadamard sign profile: points `x in F_2^k`, cuts
`S_a = {x : a.x = 0}`, `a != 0`, so `r = n - 1`, `e_p = 1` and squarefree. Let `(A, m)` be the slab-box of
`determinant.md` Theorem 2.1 (`T = n`, `w = n sqrt(8(n-1) ln 4r)`).
1. Take the `k` columns `a = e_1, ..., e_k`. For each of them `<k, a_p> = <k, beta_p>/2` is an integer in `[-w/2, w/2]`
   on `A`. Hence `h_A(a_p) <= 2 floor(w/2)`.
2. By refereed `determinant.md` Prop 5.2 (with `kappa = 0`), the complex subtorus
   `Y = {delta = sum_j t_j a_{e_j} : t in C^k}` lies in `Bs(A,m)` as soon as `k * 2 floor(w/2) < m`.
3. By Theorem 2.1(b), `m >= 0.7 n(n-1)`, so it suffices that `log_2(n) sqrt(8(n-1) ln(4(n-1))) < 0.7(n-1)`.
   * This holds at `n = 2^16`: `40921 < 45874`.
   * It holds for all larger `n`, since the ratio of the two sides grows when `n` doubles (a factor `~1.39` against
     at most `17/16`).
4. With the exact `m` (`r_braid_threshold.py`) the first Sylvester order is `n = 2^16 = 65536`:
   `k * 2 floor(w/2) = 2.682 * 10^9 < m = 3160214892`.
   * This `n` is in the regime: `n - 1 >= 300 ln 4r` and `m >= 4w`.
   * For `k = 12, ..., 15` this argument does not apply.
5. On `Y`, `delta_x - delta_y = sum_j t_j (a_{e_j,x} - a_{e_j,y})`. The map `x -> (a_{e_j,x})_j` is injective on
   `F_2^k`, so for `t` off countably many hyperplanes:
   * all `delta_x` are distinct;
   * `w_x = e^{i(delta_x - delta_n)} != 1`;
   * `w_x != w_y`.
6. These points of `Bs(A,m)` are off the braid arrangement. They exist arbitrarily close to `e` and on the real
   torus. By the implication above, `w - 1` is a common zero of all realisable symbols off the braid arrangement.
   **So ZE-braid(A,m) fails** for every Sylvester order `n = 2^k >= 2^16`, at Theorem 2.1's own parameters in the
   regime. QED

The same argument works for any profile in which some `s` columns separate all points and `s * max h_A(a_p) < m`. For
random sign profiles `s ~ 2 log_2 n` suffices, and `m/w ~ 0.25 sqrt(n/ln 4r)` eventually exceeds it.

**Consequences.**
* These points are *degenerate*: they lie on a `k`-dimensional rational subtorus. So they do not contradict ZE. They
  show that in the regime **`Bs` contains degenerate subtori with no collisions**.
* ZE-braid is strictly stronger than ZE and false. Any algebraic reformulation must allow rational subtori, and
  rational subtori are not cones in `x = w - 1`.
  * Only the braid ones `{w_x = w_y}` are linear in `x`.
  * So the "no exponentials" feature of ZE-braid cannot be kept as stated.
* **Why `n = 8` shows only the braid arrangement.** In `M8` one column direction costs range `4` against `m = 8`, and
  in `M8'` range `6` against `m = 11`. So at most **one** column subgroup is affordable (Prop 5.2 needs
  `sum ranges < m`), and a one-column subgroup moves the 4 points of a cut together: a collision.
  * In the regime `m / w -> infinity` (`>= 0.25 sqrt(n/ln 4r)`), so subgroups spanned by `~ m/w` columns, with
    collision-free generic points, are in `Bs`.
  * The braid structure of Theorem 5.2 is therefore a **small-`m/w` artefact**, not a feature to extrapolate.

### 2.9 Is the `n = 8` evidence "faithful to the regime"? (major, interpretive)

"Faithful" in `compute.md` means `zeta(A) < m <= rho(A)`. Prop 2.2 shows this is a *necessary* condition for a
support of an actual nondegenerate cluster in the regime. It is not sufficient, and the `n = 8` models miss the
regime's two defining features.
1. **LP supercriticality.** In the regime `sum_p omega_p h_A(beta_p) <= 2w <= m/2` (Theorem 2.1(d)). That factor 2
   is what makes Theorem 4.1 force (4.2).
   * The README's model condition (ii) is the weaker `h_A(beta_p) < m`.
   * `M8` has `h_A = 8 = m`, and `M8'` has `h_A = 12 > 11 = m`, with uniform weights `1/7`. **Neither satisfies the
     README's (ii)**, let alone the regime's `h_A <= m/2`.
   * In the language of Theorem 4.1, the forced-divisibility argument would not force `F_c(z) = 0` on these
     supports.
   * The author flags "(ii'; not (ii))" in the table, but §0 item 2, §8 item 3 and the Result line still present
     them as faithful models and as "evidence for ZE in the sense of README §8.5, first bullet". README §8.5
     presupposes models built under §8.4 (i)-(iii).
   * I checked whether (ii) is reachable at `n = 8` for `H_8`:
     * `T = 1`: `L = 1` gives `rho = 3 < h_A = 4`; `L = 2` gives `rho = 8 = h_A`; `L = 3` gives `rho = 11 < 12`.
     * `T = 2`: `L = 1` gives `|A| = 57`, `rho = 3`; `L = 2` gives `|A| = 2507`, `rho <= 8 = h_A` (rank mod `p`, an upper
       bound).
     * So (ii) is unreachable for this profile at `T = 1` (every cutting `L`) and at `T = 2`, `L <= 2`. I did not
       compute `T = 2`, `L = 3` (`|A| = 16857`, `h_A = 12`) or larger boxes. The README's finite test may have no
       admissible instance at `n = 8` at all. That is worth recording, but it is not "evidence for ZE in the regime".
2. **Scale of `m / w`.** In the regime `m/w -> infinity`. At `n = 8` it is about 2. §2.8 shows that this changes the
   qualitative structure of `Bs` (multi-column subtori), so the models do not test the configurations ZE must
   handle.

So Theorem 5.2 is a correct exact theorem about two toy models. As *evidence for ZE* it is weak. The honest reading
is: "in the two computable models that are free of zonotope certificates, the base locus is toric; these models are
LP-subcritical and far from the regime's `m/w`."

### 2.10 §7 (fakes) and the numerical checks

* The statement "no nondegenerate (hence no fake) position lies on `Z_P` for `M8`, `M8'`" follows from Theorem 5.2.
  Correct.
* The 100/100 and 80/100 certifications (`fakes_m8.py`) are labelled as a sanity check, and the 20 uncertified
  cases are disclosed. Acceptable. They are superseded by the theorem, as the author says.
* `explain.py`, `tc_search.py` and `bs_search*.py` are floating point and labelled evidence. Nothing rests on them.
* **Unproved remark (minor).** "Each single `P_c` has non-toric zeros through `e` **with nondegenerate real
  points**" (§5.2 Remarks, §8 item 3) is asserted, not shown.
  * On the real torus, `P_c(e^{iu}) = e^{i<a,u>} H_c(e^{iu} - 1)` is complex-valued for the non-symmetric basis
    vectors, so `P_c = 0` is two real equations.
  * The complex statement (`V(H_c)` is a cone in `x`, not a subtorus) is correct. The real nondegenerate part is
    plausible but unproved. It should be labelled as such, or proved as in Prop 6.2, which needs `c` even.
  * Floating-point evidence (`r_single_real.py`) is consistent with it. For the `M8'` cubic
    `x2x3x4 + ... + x4x5x6`, 26 of 40 least-squares starts give real `u` with `|u| ~ 0.06-0.2` and
    `|H(e^{iu}-1)| < 1e-14`, and no relation `<v,u> = 0` with `|v|_inf <= 3` to within `1e-4 |u|`. This is evidence
    only.
* **§6 "Over sub-clusters" (minor).** "Its own slab-box is zonotopal whenever `|J| <= 5`" is stated as a fact. It is
  proved for `|J| = 3` (Prop 2.4 as repaired), exhaustive in a range for `|J| = 4`, and sampled for `|J| = 5`.

## 3. Independent exact computations (results)

| claim in `compute.md` | referee result | file |
| --- | --- | --- |
| `M8`: `|A| = 1037`, `rho = 8`, `dim V(A,8) = 42`, `V_H = V` | reproduced (exact rank 995; LLL sup-norm 6) | `r_had8_out.txt`, `r_kernel_L2_out.txt` |
| `M8'`: `|A| = 1877`, `rho = 11`, `dim V(A,11) = 7`, `V_H = V` | reproduced (exact rank 1870; LLL sup-norm 4) | `r_had8_out.txt`, `r_kernel_L3_out.txt` |
| `M8'` basis = monomial x collision binomials x one form in `x` (3 cubics, 4 quartics) | reproduced with my own basis (same forms) | `r_factor_L3_out.txt` |
| `M8`: 24 form-type / 18 non-form; 27 two-term factors, all `w^a - w^{a'}` | reproduced; **7** non-linear two-term factors in my basis (author: 6), which depends on the basis | `r_factor_L2_out.txt` |
| Theorem 5.2(a), (b): saturations end in the unit ideal | **reproduced** (flint Groebner, own code, `G_j` = products) after `x_0..x_6, x_0-x_1, ..., x_2-x_3` | `r_saturate_L3_out.txt`, `r_saturate_L2_out.txt` |
| §4 tangent-cone facts for `M8` | reproduced | `r_tc_m8_out.txt` |
| skew double core `q = 7`: `|A| = 757`, `dim V = 1`, `rho = 8`, `c` even, sup 42, support 687 | reproduced from scratch (any base point) | `r_sdc7_out.txt`, `r_sdc7_c_out.txt` |
| `Phi = 8! g`, 196 terms, coefficients `{-2,-1,1,3}`, `g` irreducible, sign change | reproduced (own positive witness) | `r_sdc7_c_out.txt` |
| sdc7: `zeta <= 7` | reproduced by a different rigorous argument (min cut cost 8) | `r_sdc7_c_out.txt` |
| `n - 1 = 3`: 355 distinct admissible sets, all `zeta = rho` | reproduced (355; all witnessed) | `r_scan4_out.txt` |
| `n - 1 = 2`: no faithful model | reproduced under (ii') for `T <= 8`, `L <= 7` (55 sets); proof repaired for all cases (§2.2) | `r_scan3_out.txt` |
| `n - 1 = 4`: `zeta = rho` in the sample | my own sample: 75 distinct sets, all `zeta = rho` (73 by search, 2 by factorisation) | `r_scan5_out.txt`, `r_n5_hard_out.txt` |
| (new) `H_8`, `T = 2`: condition (ii) reachable? | no for `L = 1, 2` (`rho = 3`; `rho <= 8 = h_A`) | `r_had8_T2_L1_out.txt`, `r_had8_T2_L2_out.txt` |
| (new) ZE-braid in the regime | **false** for Sylvester `n = 2^k >= 2^16` (§2.8) | `r_braid_threshold_out.txt` |

All exact claims that a theorem or proposition of `compute.md` uses were confirmed. No numerical claim is
presented as proof, except the informal remark flagged in §2.10.

**My own `n - 1 = 4` sample** (`r_scan5.py 7 12 3 200000`, `r_n5_hard.py`; seed 7, 12 random full-rank sign
profiles with `r` between `n` and `3n` columns, `T, L <= 3`, conditions (i), (ii'), (iii), (v)):
* There are 75 distinct admissible slab-boxes.
* For 73 of them my zonotope search finds a witness with `sum L_i = rho`.
* The other two hit the search limit (both at `T = L = 3`, `rho = 12`). For both, an exact LLL basis of `V(A,12)`
  (`dim 4` and `dim 1`) contains a polynomial that factors over `Q` as a monomial times two-term factors
  `w^a - w^b` of total order 12. An example is `(w2-w3)^3 (w1-w3)^3 (w0-w2)^3 (w0-w1)^3`. This is a zonotope
  certificate (the argument of §3.2).
* So **all 75 have `zeta = rho`**. This agrees with the author's sample, and it is still only a sample.

## 4. Issues (fatal / major / minor)

**Fatal: none.** Every exact statement that carries a proposition or theorem was confirmed independently.

**Major.**
1. **The recommended next target, ZE-braid, is false in the regime** (§2.8).
   * Proof: refereed `determinant.md` Prop 5.2 plus Theorem 2.1(b).
   * For the Sylvester-Hadamard profile of every order `n = 2^k >= 2^16`, the regime slab-box has a
     `log_2 n`-dimensional rational subtorus in `Bs(A,m)` whose generic points have no collisions.
   * ZE-braid implies ZE but is strictly stronger, and false. The braid-arrangement pattern of Theorem 5.2 is an
     artefact of `m/w ~ 2`, where at most one column subgroup is affordable.
   * The "purely algebraic, no exponentials" reformulation must be replaced by one that allows rational subtori,
     which are not cones in `x = w - 1`.
2. **"Faithful" is only a necessary condition, and the `n = 8` evidence is overstated** (§2.9).
   * `M8` and `M8'` fail the README's (ii) (`h_A = 8 = m` and `12 > 11`), and miss the regime's LP margin
     `h_A <= m/2` by a factor `>= 2`.
   * At `n = 8` (ii) is unreachable for `H_8` at `T = 1` and at `T = 2`, `L <= 2` (§2.9).
   * Theorem 5.2 is a correct theorem about two LP-subcritical toy models. It should not be called "evidence for ZE
     in the sense of README §8.5", and the summary's "first faithful models" needs the qualifier "free of zonotope
     certificates; LP-subcritical".

**Minor.**
1. Prop 2.4's proof uses (ii), while the comparison models use (ii'). The conclusion holds for every slab-box at
   `n - 1 = 2`; the general zonotope proof is in §2.2.
2. "Every model with `m <= 2T` is unfaithful" needs `L >= T`, or no cut separating the coordinate from the base
   point (§2.2).
3. Prop 6.1 needs `A' ⊂ B_T`, or `log diam A'` in the condition on `W` (§2.7).
4. Lemma 5.1's two-term factors should be "of the `P_j`, `j in J`"; a primitive `a - a'` is needed for "subtorus
   through `e`" (§2.5).
5. "ZE holds in M8 and M8'" is a model statement with ad hoc `(T, L, m)`. It is not ZE(8), which is vacuous, and not
   ZE at Theorem 2.1's parameters, where at `n = 8` the slabs do not cut (§2.6).
6. Theorem 5.2(a) writes `x, y in {1, ..., 8}` while the computation indexes points `0..7` with base point 7.
   "Six four-term relations" in (b) depends on the basis (my basis gives seven different ones).
7. "Each single `P_c` has non-toric zeros through `e` with nondegenerate real points" is unproved (§2.10).
8. The scope statements are overgeneralised:
   * "the model sizes `n - 1 = 2, 3, 4` are all unfaithful": proved at 2; exhaustive only for `T <= 4`, `L <= 3`
     (one common `L`) at 3; a sample at 4;
   * "zonotopal whenever `|J| <= 5`" (§6).


## 5. Corrected statement

**Status: evidence (exact model theorems). No growth bound and no proof of ZE.**

1. **Zonotope certificates (Prop 2.1; proved).**
   * If the slab-box `A` contains a lattice zonotope `k_0 + sum [0,L_i] v_i` with `sum L_i >= m`, an explicit
     `c in V_H` gives `sum c_k w^k = w^{k_0} prod (w^{v_i}-1)^{alpha_i}`, and `Bs(A,m)` lies in finitely many
     subtori `{w^{v_i} = 1}`.
   * For an actual cluster with nondegenerate positions satisfying (4.1), its own slab-box (Theorem 2.1, `T = n`)
     has `zeta(A) < m` (Prop 2.2; proved).
   * The same holds for any support `A' ⊂ B_T` with `sum_p omega_p h_{A'}(beta_p) <= m'/2` (Prop 6.1, with
     `A' ⊂ B_T` added).
   * Hence `zeta(A) < m <= rho(A)` is a **necessary** condition for a model to resemble a regime support. It is not
     sufficient.
2. **Small sizes.**
   * `n - 1 = 2`: every slab-box from sign profiles is a lattice zonotope with `zeta = rho` (proved; the proof is in
     §2.2 of this report, and it needs no condition (ii)).
   * `n - 1 = 3`: all 355 distinct admissible slab-boxes from full-rank sign profiles with one common `L`, `T <= 4`,
     `L <= 3`, under (ii'), have `zeta = rho` (exact; reproduced).
   * `n - 1 = 4`: every sampled slab-box has `zeta = rho` (exact; sample only).
   * So the README §8.4 test at these sizes, in the ranges computed, is settled by zonotope (character-type)
     certificates.
3. **`n = 8`, Sylvester `H_8`, base point 8, `T = 1`.**
   * `M8 = (L = 2, m = 8)` and `M8' = (L = 3, m = 11)` satisfy `zeta < m = rho`, `V_H = V(A,m)`, and condition (ii').
   * They **fail** the README's (ii) (`h_A(beta_p) = 8, 12 >= m`) and the regime's `h_A <= m/2`.
4. **Theorem 5.2 (exact; independently confirmed).**
   * In `(C^*)^7`, every point of `Bs(A_3, 11)` satisfies `w_x = 1` or `w_x = w_y`.
   * Every point of `Bs(A_2, 8)` lies on `{w_x = 1}`, on `{w_x = w_y}`, or on one of finitely many explicit subtori
     `{w^a = w^{a'}}` through `e`. (The list of the latter depends on the basis used; there are six in the author's
     basis and seven in mine.)
   * Hence no nondegenerate `delta` lies on `Z_P` for these two `(A, m)`.
   * This is a model statement with ad hoc `(T, L, m)`. It is not ZE(8), which is vacuous, and not ZE at Theorem
     2.1's parameters, where at `n = 8` the slabs do not cut.
5. **Prop 6.2 (proved from exact data; independently confirmed).** For the skew double core at `q = 7`, `T = 1`,
   `L = 2`, `m = 8` (`zeta <= 7 < 8 = rho`, `dim V = 1`), `Z_P` contains nondegenerate real points arbitrarily close
   to `e`. So `zeta < m` alone does not give ZE-real.
6. **ZE-braid is false in the regime** (§2.8 of this report, proved).
   * For the Sylvester-Hadamard profile of order `n = 2^k >= 2^16`, Theorem 2.1's slab-box has a `k`-dimensional
     rational subtorus, with collision-free generic points, inside `Bs(A,m)`.
   * ZE-braid is therefore not a valid target. Any algebraic replacement must allow rational subtori, not just the
     braid arrangement.
   * The braid pattern of Theorem 5.2 comes from `m/w ~ 2` at `n = 8`.
7. **Reading.**
   * The two computable models that are free of zonotope certificates have a toric base locus, while a single form
     has non-toric zeros.
   * This is consistent with ZE, but it is weak evidence: the models are LP-subcritical, and they miss the
     multi-column subtori that are present in the regime.
   * No statement about ZE in the regime follows.
