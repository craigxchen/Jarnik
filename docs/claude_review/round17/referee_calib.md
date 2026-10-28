# Referee report on round17/calib.md (CIRCLE(q) calibration)

Referee: adversarial pass, independent code (`referee_calib_checks/`). Date: 2026-10-28.

## 0. Verdict

**calib.md survives as what it says it is: calibration evidence plus one small proved reduction. CIRCLE(q) is
undecided. No growth improvement, no uniform bound. `M <= (2/3 + eps) log R/loglog R` stands.**

* **Prop 2.1 (isotypic reduction): proved, all quantifiers check (§1).**
  * Part (1) needs only translation invariance. I give an independent two-line proof (translation invariance of
    `Psi_A` kills every frequency `2l` with `0 < |l| <= m < q`).
  * The hypothesis `m < q` is used, and it fails exactly at `m = q = 7` (sdc7, `L = 2`).
  * Its consequence (2.1) bounds `kappa_V`, never `kappa_H`.
* **Lemma 7.1: proved (§2).**
* **Numbers: everything I recomputed agrees (§3).**
  * *Exact, symmetry-free, rank certified over `Q`, no multimodular step:* 24 sdc cases `(q, L, m)` plus 5 uncut-cube cases, identical rationals,
    including the 66-digit sdc11 `m = 5` value.
  * *Own float engine:* sdc19, `m = 5, 6`, and sdc23, `m = 5`, to 10 digits, with the same eigenspace ranks.
  * *Also reproduced:* the LP values, the regime bookkeeping, and the box counts.
* **Labelling: mostly honest, four corrections needed (§4).**
  * "rigorous `k_lo`" for multimodular rows, where the unexcluded failure mode goes the wrong way for a lower bound;
  * a `0.28` that is `0.27`;
  * an unsupported "`kappa_H` within factor 4 of `kappa_V`";
  * one datum with no saved output.
* **Extrapolation (§5).**
  * Labelled evidence, not used in any proof. Toys are not passed off as regime boxes.
  * But "`m/rho(A)` governs the decay" is not shown: `m/rho` and `m/(2L+1)` are confounded in the toys.
  * The margin `e^{4.8e5}` is cube-model-specific.
  * "Robust no" must say **`kappa_V`**: CIRCLE is a `kappa_H` question, and the evidence is silent at regime height.
  * A model-free comparison still supports the qualitative reading. A `kappa_V`-cap needs an average decay of
    `>= 466` per degree; the steepest single step in any toy is `7.0`.

`result_survives = true` in the sense of the write-up's own (evidence-level) claims, after the corrections of §6 and
the corrected statement of §7.

## 1. Prop 2.1 (isotypic reduction): re-derived, survives with all quantifiers

Setting: `q` an odd prime, `A ⊂ {-1,0,1}^{F_q}` finite, `mu_A` uniform, `1 <= m < q`, `Pi` the orthogonal projection
of `L^2(mu_A; C)` onto `Poly_{<m}|_A ⊗ C`.

* **(1) Constancy on `C_1`. Correct, and it needs only translation invariance.** `l_{u_theta} = (e^{i theta} khat_1 +
  e^{-i theta} khat_{-1})/2` with `khat_{-1} = conj(khat_1)` (since `k` is real). The translation `g` preserves `A`, `mu_A`
  and the degree filtration, so `Pi` commutes with `f -> f o g`. `F_i o g = zeta^{2i-m} F_i`. The exponents `2i - m`,
  `0 <= i <= m`, are pairwise incongruent mod `q`, because `0 < 2(i'-i) <= 2m < 2q` and `q` is odd. So the `R_i` are
  eigenvectors of a unitary operator for distinct eigenvalues, hence pairwise orthogonal. The complex projection of
  the real function `l^m` is real and equals the real projection, so `Psi_A(u_theta) = 4^{-m} sum_i binom(m,i)^2
  ||R_i||^2`, which does not depend on `theta`. The same holds for `E l^{2m}`.
  *Independent two-line proof (referee).* `Psi_A(u o g) = Psi_A(u)` and `u_theta o g^{-1} = u_{theta + 2 pi/q}`. Also
  `theta -> Psi_A(u_theta)` is a trigonometric polynomial with frequencies `2l`, `|l| <= m` (it is a sum of squares
  of forms `Phi_c(u_theta)` of degree `m`, reduction Thm 2.1(5)). Invariance under the shift `2 pi/q` kills every
  frequency with `q ∤ l`, so every frequency with `0 < |l| <= m < q` vanishes. This is the R16-J3 mechanism applied to
  the `L^2` form.
  *Exact check (referee, symmetry-free code, §3).* `<R_i, R_{i'}> = 0` holds exactly in `Q(zeta)` for every `i != i'`
  in every case with `m < q` that I computed. At `m = q = 7` (sdc7, `L = 2`) it **fails**, so constancy does not hold
  there. Hence the hypothesis `m < q` is used, and it is needed.
* **(2) Rationality. Correct.** `Pi` has a rational matrix (the integer evaluation matrix `M` gives
  `Pi = M (M^T M)^+ M^T`) and `F_i(k) in Z[zeta]`. `sigma_t(F_i(k)) = F_i(pi_t k)` holds with `(pi_t k)_y = k_{t^{-1}y}`;
  I checked this from `khat_1(pi_t k) = sum_y zeta^y k_{t^{-1}y} = sum_x zeta^{tx} k_x`. Since `pi_t` preserves `A`
  (the AGL(1,q) invariance of `A_L(q)` is proved in referee_barrier §4, and I re-derived it: `x -> tx` maps `P_a` to
  `P_{ta}` if `chi(t) = 1` and to `P'_{ta}` if `chi(t) = -1`), the norm `||R_i||^2` is Galois-fixed, hence rational.
  `||R_i|| = ||R_{m-i}||` follows from `F_{m-i} = conj(F_i)`. Note that (2) uses multiplication invariance and (1) does
  not.
* **(3) Correct** (it is the eigenvector statement).
* **(4) Correct.** The `2q` points `theta + pi x/q` are equally spaced with gap `pi/q` (here `q` odd is used), so
  `|u_theta|_inf >= cos(pi/(2q))`. `Psi_A` is homogeneous of degree `2m` in `u`, and `|u_0|_inf = 1`. Then
  (2.1) follows from `Psi_A^{1/2} <= Cheb_A = d^m kappa_V` (reduction Thm 2.1(1)-(2), refereed).
* **Scope (not an error in the proposition, but in how it is used).**
  * (2.1) bounds `kappa_V`, not `kappa_H`. CIRCLE is a `kappa_H` statement and `kappa_H <= kappa_V`, so (2.1) can
    never certify "no". calib says this in Remark (d) and §7.3.
  * (2.1) is rigorous for a given `(q, L, m)` only if `Psi_A(C_1)` is known rigorously **from below**. calib's "E"
    values are multimodular. The residual is computed relative to the span of pivot columns `S` chosen modulo one
    prime. If `rank_Q > rank_p` at that prime, the computed number is the residual relative to a **proper** subspace,
    which **over**-estimates the true `Psi_A`. So the unexcluded failure mode of the multimodular caveat goes in the
    direction that would invalidate a lower bound.
    * calib §3.2 states the caveat. §0 and the §4.1 header nonetheless call `k_lo` "rigorous whenever `Psi_A(C_1)` is
      exact" without carrying the caveat along.
    * The referee's fully exact recomputation (§3) removes the caveat for the rows listed there. For the other E rows
      (`L = 1`, `q = 31..103`; `sdc19`, `sdc23` at `L = 2`, `m <= 4`; `sdc11`, `L = 2`, `m = 6, 7, 8`), `k_lo` is
      exact **modulo the multimodular caveat**, not rigorous. (The values are very probably right; the label is what
      is wrong.)
* **"Psi_A(hat u_theta) = Psi_A(C_1)/|u_theta|_inf^{2m}"** is correct. One wording slip: Remark (a) says Prop 2.1
  "sharpens barrier Thm 3.2". It does so for the `L^2` form only. It says nothing new about the sup form
  `kappa_V`, as calib itself notes. **No error.**

**Verdict on Prop 2.1: proved, with all quantifiers.** It is a short and correct observation. It is not new
mathematics beyond R16-J3 plus Galois invariance, but it is a correct and useful reduction.

## 2. Lemma 7.1 (`rho(A) >= d + 1` when `|A| >= 3^d/2`): correct

`Poly_{<rho}|_A = R^A` needs at least `|A|` reduced monomials `alpha in {0,1,2}^d` with `|alpha| <= rho - 1`, because
these span `Poly_{<rho}` on the ternary cube (`k^3 = k`). The involution `alpha -> 2 - alpha` maps `{|alpha| <= d-1}`
bijectively onto `{|alpha| >= d+1}`, and `{|alpha| = d}` is nonempty, so `#{|alpha| <= d-1} < 3^d/2 <= |A|`. Hence
`rho - 1 >= d`. **Proved.** Its premise `|A| >= 3^d/2` for regime boxes is quoted from t1rigour §2 (refereed record).
`m'/(d+1) = 2057/3252 = 0.6325` at `q = 3251` (checked).

**Caveat on its use.**
* The lemma is used in §7.2 to place the regime at "comparable `m/rho(A)`" to the toys.
* calib §5.4 claims the decay "is governed by `m/rho(A)` ... not by `m/(2L)`". **The data cannot separate these two
  parameters.** Only `L in {1, 2}` occur, with `(2L+1, rho) = (3, 4)` and `(5, 9)`, so the two ratios move together.
* In the regime they differ enormously: `m'/rho <= 0.633`, but `m'/(2j+1) = 2057/515 ~ 4.0`.
* Measured against the slab degree `2L+1`, every toy has `m/(2L+1) <= 1.6` (sdc11, `m = 8`). Once `m > 2L+1`, every
  toy shows `rel_A/rel_cube` falling:
  * sdc11: `0.67, 0.55, 0.25, 0.00066` at `m = 5, 6, 7, 8`;
  * sdc19: `0.88, 0.55` at `m = 6, 7`.
* So the choice of `m/rho` as the extrapolation parameter is an **assumption**, not something the data show (see §5).

## 3. Independent recomputation (referee code only; nothing imported from `calib_checks/`)

All scripts are in `referee_calib_checks/`, with outputs in `referee_calib_checks/out/`. Each run took under 6 minutes.

**3.1 Boxes.**
* `enum_check.py` (numpy brute force over `{-1,0,1}^q`) and `enum2.c` (my own DFS) give `|A_L(q)| = 57, 133`
  (`L = 1`; equal to `{0, +-e_x, e_x - e_y}`, i.e. `1 + q(q+1)`) and `757, 9373, 264291, 818755` (`L = 2`,
  `q = 7, 11, 19, 23`). This agrees with calib §1.
* I re-derived the support bound `|supp k| <= 4L^2 + (2L-1)^2/q` (Parseval with `|chihat(b)|^2 = q`, and
  `sum_a k(P_a) = S(q+1)/2`). It is correct.
* My count-only DFS for `q = 31`, `L = 2` did not finish within its 25-minute limit (fewer than `2^24` points found).
  This is consistent with calib's report that `L = 2` is out of reach for `q = 31, 43`.
* The `L = 1` boxes for `q >= 19` were taken from the closed form. Brute force confirms that form only for
  `q = 7, 11`; for larger `q` I rely on calib's enumeration.

**3.2 Fully exact, symmetry-free `Psi_A(C_1)` (`exact_full.py`).** The method:
* exact rank and pivots from the fraction-free rref of the integer Gram `M^T M` (python-flint `fmpz_mat`);
* exact rational solve, verified by multiplying back;
* the coefficients of `F_i = khat_1^i khat_{-1}^{m-i}` in the group ring `Z[Z_q]`, reduced in `Q(zeta_q)` at the end.

This involves no eigenspace decomposition and no multimodular step. It also checks Prop 2.1 itself: it tests
`<R_i, R_{i'}> = 0` in `Q(zeta)` for `i != i'`, and that each `||R_i||^2` is rational.

| box | `m` | result |
| --- | --- | --- |
| sdc7, `L = 1` | 1-4 | identical rationals to `calib_checks/out/exact_q7_L1.txt` (e.g. `6895/6498`, `735/608`; `dim V(A,4) = 0`) |
| sdc7, `L = 2` | 1-6 | identical rationals (e.g. `m = 6`: `1317448821...587501875/379887844...018688`), `dim V = 756, 749, 721, 644, 483, 232` |
| sdc7, `L = 2` | 7 (`= q`) | **off-diagonal `<R_i,R_{i'}> != 0`**: `Psi_A` is not constant on `C_1` at `m = q`, so Prop 2.1's `m < q` is sharp here |
| sdc7, `L = 2` | 8 | `dim V = 1`; only `R_4 != 0`; `Psi_A(C_1) = 420175/18168 = 23.1272016733` (= calib §6). Also `|c|_2^2 = Phi_c^2/(U Psi) = 3402` exactly for the record's octic (`Phi_c = 7717.5`), consistent |
| sdc11, `L = 1` | 2-4 | identical (`39765/35378`, `363/304`, `0`) |
| sdc11, `L = 2` | 2-5 | identical, including the 66-digit `m = 5` rational (`Psi = 17047.761788`, `dim V(A,5) = 8140`) |
| sdc19, sdc23, `L = 1` | 2, 3 | identical (`1.17254`, `1.18716`, `1.16586`) |
| full cube `q = 7` | 2-6 | identical to calib's `cube_q7.txt` (so the formula of reduction Thm 5.1(a) is confirmed exactly at these points) |

For these rows `k_lo = Psi_A(C_1)^{1/2}/d^m` is now a **rigorous** lower bound for `inf_theta kappa_V(hat u_theta)`.
The multimodular caveat (§1) is removed for them.

**3.3 Floating-point (evidence), own isotypic engine (`iso_float.py`, complex128, SVD rank at three tolerances).**
* sdc7, `L = 2`, `m = 4, 5, 6`, and sdc11, `L = 2`, `m = 5`: agree with the exact values to 10 digits.
* sdc11, `L = 2`, `m = 6`: `1.6740358758e5` (calib E: `167404`).
* sdc19, `L = 2`, `m = 3, 4`: `183.08893733`, `2106.1183932` (calib E).
* sdc19, `L = 2`, `m = 5`: `27393.324041`, rank `446` in every eigenspace (deficiency `13464`); calib F: `27393.3`.
  `||R_i||^2/U = 1.042551e6, 1.822760e5, 8.425931e4`, as in calib §6.
* sdc19, `L = 2`, `m = 6`: `347191.14645`, rank `2024` in every eigenspace (so deficiencies `11887` for `c = 0` and
  `11886` for `c != 0`); calib F: `347191`, ranks `2024`. `||R_i||^2/U = 3.127219e7, 4.283093e6, 1.427698e6,
  1.021760e6`, as in calib §6. The value is stable across SVD tolerances `1e-8, 1e-10, 1e-12`.
* sdc23, `L = 2`, `m = 5`: `31253.905496`, rank `739` in the six needed eigenspaces; calib F: `31253.9`.
* (My first version had a spurious rank-`+1` in eigenspace `c = q - 1`. It came from a numerically-zero constant
  column that was normalised up to unit length. It was caught by the `||R_i|| = ||R_{m-i}||` symmetry and fixed. I
  mention it because the same failure mode would affect any floating "F" row that lacks such an internal check.
  calib's F rows do satisfy the symmetry.)

**3.4 LP (`lp_check.py`, HiGHS, my own code, evidence).** For sdc7, `L = 2`, `m = 3`, at 5 `theta` in
`[0, pi/(2q)]`: `inf kappa_V = 6.5951e-2`, `inf/sup = 0.9970`; `m = 6`: `1.0726e-3`, `0.9946`; sdc7, `L = 1`, `m = 3`
(9 `theta`): `4.9933e-3`, `0.8564`. All equal the table. The fundamental domain `[0, pi/(2q)]` is correct: it is
generated by `theta -> -theta` (`x -> -x`) and `theta -> theta + pi/q`, the latter because
`u_{theta+pi/q} = -u_theta o` (translation by `(q+1)/2`).

**3.5 `V_H` heights (`hred_check.py`, exact integer kernel plus flint LLL).**
* `V(A,m)` has an integer basis with `|c|_inf <= 3` for sdc7, `L = 2`, `m = 5`, and `<= 1` for sdc11, `L = 1`,
  `m = 3`. These match calib.
* For sdc7, `m = 6`, my LLL gives max entry `7`, against calib's `4`. LLL output depends on the implementation, and
  "`V_H = V` for `H >= H_red`" is an existence statement, so this does not refute `H_red = 4`. It does mean that
  `H_red` is an upper bound produced by one reduction, not a minimum.

**3.6 Regime bookkeeping (`regime_check.py`).** I reproduced every number in the calib §7.1 table: `j`, `m'`, `tau_D`,
`W`, `ln T`, the needed `ln Psi_A(C_1)`, and the needed decay `466.1, 613.1, 841.1, 2805.9, 9320.6`. Also the
first regime prime `q = 3251` and the "floor" variant (`tau = 2 ln 4r = 20.33`, needed decay `30.3`).

**3.7 Uncut cube at the regime (`cube_regime.py` at 1500 and 2500 digits, cross-checked by the log-domain product
`cube_logprod.py`).**
* At `(q, m) = (3251, 2057)` on `C_1`: `ln Psi_cube = 27188.36`, `ln E_cube l^{2m'} = 28949.37`. calib's Gaussian
  model gives `29436.5`.
* Per-degree decay: `-ln(rel_cube)/m' = 0.856` (and `0.814` at `(103, 65)`, the same `m/q`).
* calib's statement that the cube decay "tends to `ln 2 = 0.69` per degree" holds only for `m << q`. At the regime
  ratio `m'/q = 0.63` it is `0.86`. This is harmless (it changes the margin by about `350` in `ln Psi`), but it
  should be corrected.
* (A 40-digit run of the same series gave a negative `Psi`, from cancellation. The 1500- and 2500-digit runs and
  the cancellation-free log product agree to 7 digits.)

## 4. Exact vs floating labelling

On the whole the labelling is honest: every F row is marked F, LP values are marked evidence, and §7 is headed
"LABELLED EVIDENCE". These are the places where a label overstates what was shown.

1. **"E" means multimodular.** calib §3.2 says so ("exact modulo the standard multimodular caveat"). But §0 ("Hence
   ... rigorously"), the §4.1 header ("`k_lo` ... a **rigorous** lower bound ... whenever `Psi_A(C_1)` is exact") and
   §5.7 ("Using the rigorous lower bound `k_lo`") drop the caveat.
   * As shown in §1, the failure mode the caveat leaves open (`rank_p < rank_Q` at the pivot prime) would make
     `Psi_A` too **large**, which is exactly the wrong direction for a lower bound.
   * Corrected: `k_lo` is rigorous for the rows recomputed in §3.2. For the remaining E rows it is exact modulo the
     multimodular caveat. For F rows it is evidence.
   * Rational reconstruction accepted after "stable over 3 further primes" is also heuristic unless a height bound is
     used. Here this is moot for the rows I recomputed, since all of them agree.
2. **"The circle is never in `TC_V` below saturation"** (§0, first bullet).
   * This is exact (now rigorous, by §3.2) for the recomputed rows, and exact modulo the caveat for the other E rows.
   * For the F rows (sdc19, `m = 5, 6, 7`; sdc23, `m = 5, 6`) it is floating-point evidence, although robust
     (`rel >= 0.012`).
   * It is a statement about the computed toy `(q, L, m)` only.
3. **§5.6: "`kappa_H(hat u, H = 1) >= 0.28 sqrt(Psi_A)/d^m` on sdc7, `m = 5, 6`".** calib's own output
   (`out/kappaH_7_2_5.txt`, `kappaH_7_2_6.txt`, at `theta = 0`) gives ratios `6.254e-4/2.248e-3 = 0.278` and
   `1.370e-4/5.006e-4 = 0.274`. So the constant is **`0.27`**, and the result is floating point.
   * It is computed at two `theta` only.
   * It does extend to the whole circle, by an argument calib does not give. `V_H` is `AGL(1,q)`-invariant, because
     permuting coordinates preserves integrality and `|c|_inf`. So `Psi_{A,H}` is constant on `C_1` by the proof of
     Prop 2.1(1), and `Psi_{A,1} >= max_theta Psi_span(theta)`.
   * The "span" used (LLL vectors of height 1) is not translation-invariant. That is why the two `theta` values differ.
4. **§7.2: "`kappa_H` on the circle is `kappa_V` up to a factor `<= 4` even at `H = 1`". Not supported.**
   * What was computed are **lower bounds**: `kappa_1 >= 0.27 k_lo`.
   * `kappa_V` itself (LP) is `2.4 k_lo` (`m = 5`) and `2.1 k_lo` (`m = 6`).
   * So the data give `kappa_V/kappa_1 <= ~9` (sdc7, `m = 5, 6`, float), not `<= 4`.
   * "`V_H = V(A,m)` already for `H <= 4`" was computed only for sdc7 (`L = 2`, `m = 5, 6`) and the pair boxes
     (`L = 1`, `q = 11, 19, 23`, `m = 2, 3`). It was not computed for any `L = 2` box with `q >= 11`, and the sdc11
     `m = 8` run was killed (`kappaH_11_2_8.txt` ends in "Killed"). The §0 sentence ("in every computed case") is
     accurate. The task summary's "whenever `m < q`" is not.
5. **"inf over `theta`" in the LP column** is a minimum over a finite grid (5-9 points). That is evidence, and it is
   labelled as such. My LP reproduces it (§3.4).
6. **§6, sdc7 `m = 8`:** "`kappa_h.py`: primitive integer `c` with `|c|_inf = 42`". There is no
   `out/kappaH_7_2_8.txt` in `calib_checks/out/`, so the datum is unsupported by a saved output. It is consistent with
   the exact `Psi` (`|c|_2^2 = 3402` follows from `Phi_c = 7717.5` and my exact `Psi`), but `|c|_inf = 42` was not
   checked.
7. **Summary-level slip.** "`rel_A/rel_cube` ... 1.00 to 1.04 for `m <= 5`" (task summary) is false at sdc23,
   `m = 5`, where it is `0.98` (calib §5.5 itself prints `0.98`).
8. **Cube decay.** "Gaussian-model decay tends to `ln 2 = 0.69` per degree" (§7.2, and `regime_need.txt`'s "cube
   (Gaussian model) decay = 0.691") is the `m << q` asymptotic. The exact cube at the regime point is `0.856` per
   degree (§3.7).

None of these changes a conclusion of calib. Items 1, 3, 4 and 6 are mislabelled or unsupported sentences that need
correcting.

## 5. Extrapolation and toy-vs-regime labelling

**Toys are not presented as regime boxes.** calib §1 says so explicitly ("no toy is a regime box"), §4-§7 are
headed as toy evidence, and §7.3 caveat 2 lists the gaps:
* no toy reaches `m/(2L) = 4` with `dim V >= 2`;
* every toy has density `<= 0.35 < 1/2`.

The extrapolation is labelled evidence and is **not used in any proof**. The only proved statements are Prop 2.1,
Lemma 7.1 and the bookkeeping identity of §7.1. **This part of the brief is satisfied.**

**What is robust and what is model-dependent.**
1. *Proved bookkeeping (modulo the record's `C = 1` and `c_0 = 0.184`).* A "yes" **through `kappa_V`** at
   `q = 3251` forces `ln Psi_A(C_1) <= -9.293e5`. This uses (2.1) and the Lemma 3.4 value of `tau_D`, and I
   re-derived it.
   * The conversion to "decay `>= 466` per degree" uses a heuristic `ln E_A l^{2m'}`. The exact cube value
     (`28949`, against the model's `29437`) changes this by `< 0.3`.
   * The choice of `tau_D` is the one most favourable to "yes". `C <= 1` and `c_0 <= 0.184` (the MILP minimum over
     `q <= 43` caps any uniform `c_0`) give the smallest admissible `tau_D`, unless H_light turns out to hold at regime
     `q` with a larger `c_0`.
   * Even at `c_0 = 1/2` the other two terms of `tau_D` give `tau_D >= 62.3`. The needed decay is then still about
     `97` per degree (`regime_check.py` arithmetic).
2. *Model-dependent.* The margin "`e^{4.8e5}` in `kappa`" assumes that the regime circle decays at the cube's rate
   (`c ~ 1` per degree). That is one model among several the toys allow.
   * **Parameter choice.** As noted in §2, the toys cannot separate `m/rho(A)` from `m/(2L+1)`. Under the second
     parametrization the regime (`m'/(2j+1) = 4.0`) lies far outside the toy range (`<= 1.6`), and the toys show
     `rel_A/rel_cube` decreasing once `m > 2L+1`.
   * **A model-free comparison** therefore reads better. A "yes" through `kappa_V` needs an **average** decay of
     `>= 466` per degree over all `2057` degrees.
     * The steepest **single-degree** drop of `rel` in any toy is `ln(0.006483/6.0e-6) = 7.0` (sdc11, `m = 7 -> 8`,
       the last degree before saturation).
     * Saturation itself cannot occur in the regime, by Lemma 7.1.
     * Even if every regime degree dropped as steeply as the steepest toy step, `ln Psi` would fall by only about
       `1.4e4`, against the `9.6e5` needed.
   * So the **qualitative** reading ("far from a `kappa_V`-cap") survives every parametrization the toys support. The
     **quantitative** margin `e^{4.8e5}` belongs to the cube model and should be labelled that way.
3. *Not addressed by the extrapolation at all.*
   * **`kappa_H`.** A "no" for CIRCLE is a `kappa_H` statement. `kappa_H <= kappa_V`, so every piece of evidence here
     points to "no" **only for the `kappa_V` version**. The record says fakes may live in `TC_H \ TC_V` (README-16
     §2.1 item 5). calib §0 and the task summary say "CIRCLE is a robust 'no' by a factor `e^{4.8e5}` in kappa"; this
     must read "in `kappa_V`". §7.3 has it right ("Through `kappa_V`").
   * **Exact regime identities.** As calib caveat 3 says, toys cannot exhibit a degree-`m'` identity along `psi_1`
     that exists only at regime size.

**L = 1 rows.** `A_1(q) = {0, +-e_x, e_x - e_y}` does not depend on the Paley patterns: any pattern family that
separates every pair gives the same box. So the 39 `L = 1` rows (`q <= 103`) calibrate a generic pair box at
`m <= 3 = rho - 1`. They carry no `sdc`-specific information. calib calls them "pair boxes", but §0 and the task
summary list them as `sdc` calibration data ("exact for all primes `7 <= q <= 103`, which covers `q = 19, 23, 31,
43`"). The judge's request for `q = 31, 43` is met only in this degenerate sense.

## 6. Issues list

| # | severity | issue | where |
| --- | --- | --- | --- |
| 1 | major | `k_lo` is called "rigorous" for every E row. E is multimodular, and the unexcluded failure (`rank_p < rank_Q`) would make `Psi_A` too large, i.e. it would break a **lower** bound. It is rigorous only for the rows recomputed exactly here (§3.2); elsewhere it is exact modulo the caveat | §0, §4.1 header, §5.7 |
| 2 | major | "CIRCLE is a robust 'no' by `e^{4.8e5}` in kappa" must read "in `kappa_V`". CIRCLE is a `kappa_H` statement, and the evidence is silent on `kappa_H` at regime height | §0 (task summary) |
| 3 | major | The extrapolation parameter `m/rho(A)` is not supported over `m/(2L+1)`: they are confounded in the toys (`L in {1,2}`). Under `m/(2L+1)` the regime (`4.0`) is far outside the toy range (`<= 1.6`), where the toys show `rel_A/rel_cube` falling. The qualitative "far from a `kappa_V` cap" survives a model-free comparison (needed average decay `466`/degree against a steepest toy step of `7.0`); the margin `e^{4.8e5}` is cube-model-specific | §5.4, §7.2-7.3 |
| 4 | minor | "`kappa_H` is `kappa_V` up to a factor `<= 4` even at `H = 1`" is unsupported. The data give lower bounds `kappa_1 >= 0.27 k_lo`, hence `kappa_V/kappa_1 <~ 9` (sdc7, float) | §7.2 |
| 5 | minor | "`>= 0.28`" should be "`>= 0.27`" (calib's own outputs give `0.278` and `0.274`). It holds at all `theta` only by the `V_H`-invariance argument, which calib does not give | §5.6 |
| 6 | minor | sdc7 `m = 8` "`|c|_inf = 42`" has no saved output | §6 |
| 7 | minor | The cube decay at the regime is `0.856`/degree (exact series), not `ln 2` | §7.2, `regime_need.txt` |
| 8 | minor | The `L = 1` rows are the pattern-independent pair box, so `q = 31, 43` are "covered" only by a box that does not see the Paley structure | §0, §3 |
| 9 | minor | "1.00 to 1.04 for `m <= 5`" is false at sdc23, `m = 5` (`0.98`) | task summary |
| 10 | minor | "`V_H = V` ... whenever `m < q`" was computed for 7 small cases only (no `L = 2` box with `q >= 11`) | task summary |
| 11 | minor | The random comparison directions are uniform in the cube, not "random points of `G_1`" as the judge asked (`G_1` is undefined for toys without cluster data). This is acceptable but should be stated | §3.4, §5.2 |

No fatal issue. Prop 2.1 and Lemma 7.1 are correct. Every E value I recomputed is identical, and every F value I
recomputed agrees to 10 digits.

## 7. Corrected statement

> **Status: evidence, plus one proved reduction. CIRCLE(q) is undecided. There is no growth improvement and no
> uniform bound. `M <= (2/3 + eps) log R/loglog R` stands.**
>
> **Proved.**
> * **Prop 2.1.** For every finite `A ⊂ {-1,0,1}^{F_q}` invariant under translations, and every `1 <= m < q`:
>   `Psi_A(u_theta) = 4^{-m} sum_i binom(m,i)^2 ||R_i||^2` is independent of `theta`, with the `R_i` pairwise
>   orthogonal translation eigenvectors. If `A` is also invariant under `x -> tx`, every `||R_i||^2` is rational.
>   Hence `inf_theta kappa_V(hat u_theta) >= Psi_A(C_1)^{1/2}/d^m`.
>   * This bounds `kappa_V` only. It cannot certify a CIRCLE "no", which is a `kappa_H` statement.
>   * The hypothesis `m < q` is needed: at sdc7, `L = 2`, `m = q = 7` the cross terms are exactly nonzero.
>   * The same argument shows that `Psi_{A,H}` is constant on `C_1` for every `H`.
> * **Lemma 7.1.** `|A| >= 3^d/2` implies `rho(A) >= d + 1`. So in the regime `m'/rho(A) <= m'/(d+1) <= 0.633`, and
>   saturation (`dim V = 0`) cannot occur at `m'`.
> * **Bookkeeping.** At `q = 3251`, with `C = 1`, `c_0 = 0.184` and `tau = tau_D = 296.1`, a "yes" through `kappa_V`
>   forces `ln Psi_A(C_1) <= -9.293e5`.
>
> **Exact data (toy boxes only; none is a regime box).**
> * *Rigorous: exact rational, symmetry-free, rank certified over `Q`, recomputed by the referee.* `Psi_A(C_1)` and
>   `dim V(A,m)` for:
>   * sdc7: `L = 1`, `m <= 3`; `L = 2`, `m <= 6`;
>   * sdc11: `L = 1`, `m <= 3`; `L = 2`, `m <= 5`;
>   * sdc19 and sdc23: `L = 1`, `m <= 3`.
>
>   For these, `k_lo` is a rigorous lower bound for `inf_theta kappa_V`.
> * *Exact modulo the multimodular caveat.* `L = 1` for the other primes `q <= 103`; sdc11, `L = 2`, `m = 6, 7, 8`;
>   sdc19 and sdc23, `L = 2`, `m <= 4`.
> * *Floating-point evidence.* sdc19, `L = 2`, `m = 5, 6, 7`; sdc23, `L = 2`, `m = 5, 6`. The values sdc19
>   `m = 5, 6`, sdc23 `m = 5` and sdc11 `m = 6` were independently reproduced to 10 digits.
> * *Scope.* The `L = 1` boxes are the pattern-independent pair box. `L = 2` was not reached for `q = 31, 43`.
>
> **Evidence (toys).**
> * Below saturation the circle stays far from `TC_V`: `rel >= 0.0058` for `m <= rho - 2`.
> * The relative residual tracks the uncut cube's: `rel_A/rel_cube in [0.29, 1.31]` for `m <= rho - 2`. It falls
>   once `m` exceeds the slab degree `2L + 1`.
> * The LP `kappa_V` varies by at most 15% along the circle.
> * `V_H = V(A,m)` at `H <= 4` in the 7 small cases where it was computed.
> * `kappa_V` reaches the CIRCLE threshold shape only at `tau* <= 2.22`, against `tau_D = 164 - 221`.
>
> **Extrapolation (evidence, model-dependent, not used in any proof).**
> * A `kappa_V`-cap at `q = 3251` needs an average relative-residual decay of `>= 466` per degree (`>= 30` even at
>   `tau = 2 ln 4r`; `>= 97` at `c_0 = 1/2`).
> * The toys show at most `1.51` per degree on average, and `7.0` in the single steepest step. The exact uncut cube at
>   `(3251, 2057)` shows `0.856`.
> * So every toy-supported model puts the regime circle far above the threshold **for `kappa_V`**. The margin
>   `~e^{4.8e5}` holds in the cube-rate model only.
> * Nothing here bears on `kappa_H` at the regime height `H' = e^{m'W/40}` (the passage (P1)/(P2)), or on an exact
>   degree-`m'` identity along `psi_1` that would exist only at regime size.
>
> **Reading.** Under the judge's rule, route 5 is undecided after this round. The evidence points to the circle test
> **not** capping route 5 through `kappa_V`, but it supplies no certificate in either direction.
