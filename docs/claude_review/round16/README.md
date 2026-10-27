# Round 16 (judge): route 5, discrete anti-concentration reframing (steelman S3)

## 0. Verdict

**No growth improvement and no uniform bound** were proved this round, either unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for actual circles.

**AC / MIN^pair: open.** Nothing was proved, refuted or capped. What changed is the map around them.
* **The anti-concentration reframing (steelman S3) is the wrong frame.**
  * The form MIN^pair needs is an `L^2` statement about one vector: the top orthogonal component of `<k',u>^(m')`
    in `L^2(mu_A)`. It is not a small-ball statement.
  * The steelman's cube form `AC(u) = AC_T(u;1/2)` is strictly stronger as a statement form, by an exact toy and by
    the uncut cube.
  * `AC_T(u;1/2)` is false at explicit points of the averaged region for Sylvester `n = 2^12..2^22` and
    `2^26, 2^30, 2^31` (finite `n` only).
  * By judge check R16-J4 (unrefereed), `AC_T(u;1/2)` is false for **every** `u` at the smallest admissible `W`
    when `2966 <= n <= 3750` (`r = n`).
* **The literature cannot reach the regime, provably.** Every anti-concentration bound depending only on degree,
  variance, influence or rank is vacuous at degree `>= 2j + 1` (reduction Prop 7.1, antitheory Thm 3.1, both
  exact). Inverse Littlewood-Offord theorems are bounded-degree, and Kane's theorem is Gaussian. This does **not**
  refute anti-concentration for the relevant class `{l_u^(m') - P : u in G_1}`; it shows that tools blind to
  `u in G` are useless.

**Route 5: still open, stuck, not capped, not advanced.** The precise missing lemma is unchanged in content and
now stated sharply (§5.1):

> **(R_H)** `min_{u in G_1} Psi_{A,H'}(u) > d^(2m') eps_0^2`, with `H' = e^(m'W/40)`, for the cluster's own regime
> box `A` and gap sphere `G_1`.

* (R_H) implies MIN^pair losslessly, hence `n = O(loglog R)` on all norms (round 12).
* MIN^pair implies (R_H) at level `eps_0 3^(-d/2)`.
* It is a statement about `(J(A), V_H', G)` alone. By fake Thm 2.1 it is the barrier question.

Round 16's one new, sharp handle is on the **barrier side**: the **Fourier-circle test** for the AGL(1,q)-invariant
fake boxes `sdc_q` (barrier Thm 3.5, as corrected). Route 5 is capped if `inf_{C_1} kappa_H <= eps_0/(2 sqrt(m'+1)
e^(pi^2/(4q)))` at one regime `q`, modulo the hypothesis H_light. An exact "circle not in TC" does not settle it,
because the test is quantitative.

**Scores** (rubric of `round10/README.md` §1; details in §6):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| reduction | 6 | 1 | 1 | **4** |
| barrier | 5 | 2 | 1 | **4** |
| antitheory | 6 | 1 | 1 | **3** |

**Recommendation (§7).** Run the quantitative Fourier-circle test on regime `sdc_q`, starting with a proof of
H_light by character sums. It is the only route-5 task in the record where one outcome is decisive at regime size
(a cap). The other outcome would need the first non-zonotopal lower-bound certificate at an irrational direction of
an overdetermined regime box, which is exactly the missing tool for (R_H). Discrete anti-concentration, density
forms, locality and further reformulations of MIN^pair are not recommended.

## 1. Setting and conventions

Unless a line says otherwise, every statement is at:
* box size `T = 1`, `d = n - 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m' = ceil(4w)`;
* regime `d >= 316 ln 4r`, any profile, any norm;
* scale `d^(m') eps_0`, with `eps_0 = C d e^(-W/4)/(m'+1)`, and `H' = e^(m'W/40)` (t1rigour Remark 3.3);
* `C = C' = 1` wherever numbers are given.

Notation, following `round12/README.md` §1 and the round-16 write-ups:
* `U = |A|` and `mu_A` is uniform on `A`.
* `Cheb_A(u) = dist_inf(l_u^(m')|_A, Poly_(<m')|_A) = d^(m') kappa_V(u)` (nearstrata Prop 2.1).
* `Psi_A(u)` is the squared `L^2(mu_A)` norm of the top orthogonal component of `l_u^(m')`, and `Psi_(A,H)` is the
  same with `V_H` in place of `V`.
* `rho_A = mu(A) Psi_A/Psi_cube` is the restriction ratio.
* `AC_A` is the sup form, `AC_2` the `L^2` form, `AC_p` the "mass `p`" small-ball form on `mu_A`, and
  `AC_T(.;theta)` the cube form with threshold `theta` (the steelman's `AC` is `AC_T(.;1/2)`).
* `MIN^pair_inf` is MIN^pair with `kappa_V` (`H = infinity`).

Labels:
* *proved*: all quantifiers, survived refereeing;
* *exact*: integer or rational computation;
* *evidence*: anything else.

**No exact model this round is a regime box.** The regime starts at `n = 2966` (`r = n`). The first regime `sdc_q`
has `q = 3251` (R16-J2).

## 2. Referee-corrected statements

### 2.1 reduction (`reduction.md`, `referee_reduction.md`): core survives, two framing claims withdrawn

The referee re-derived every proof and redid the key checks with independent exact code (63 cube cases, 21 product
cases with a non-slab factor, Prop 7.1 tail probabilities), finding 0 mismatches.

1. **Sup and `L^2` forms (proved, any finite `A` in the ternary cube).**
   * `Cheb_A = d^(m) kappa_V`.
   * `Psi_A <= Cheb_A^2 <= U Psi_A`.
   * `Psi_A = (1/U) sum_i Phi_(c_i)^2` is a sum-of-squares form whose zero set is `TC_V`.
   * For the cluster box, `AC_2` (`min_(G_1) Psi_A > d^(2m') eps_0^2`) implies `MIN^pair_inf` losslessly.
     `MIN^pair_inf` implies `AC_2` only at the lower level `eps_0 3^(-d/2)`.
   * The Gaussian "easy half" is the identity `Psi_(N(0,I)) = m'! |u|_2^(2m')`.
2. **Probabilistic forms (proved in one direction only).**
   * `AC_p` (fixed `p > 1/U`) and `AC_cube` (cube measure, threshold `mu(A)`) each **imply** `AC_2`.
   * The converses are **false in general**. Exact toy: `d = 2`, `m = 3`, one-cut `A` with 7 points,
     `u = (1/2, 1)`.
     * `Psi_A(u) = 3/8` in sum form, and `TC_V(A,3)` lies in the braid arrangement, so `MIN^pair_inf` holds.
     * A 7-point `B` in the cube has `Psi_B(u) = 0`, so `AC_cube` fails.
     * A 6-point subset of `A` has `Psi = 0`, so `AC_p` fails for every `p > 1/7`.
     * The judge reproduced all of this (R16-J1).
   * Only the threshold `p = 1/U` ("some `k`") is equivalent up to `e^(O(d))`.
   * The `e^(O(d))` losses are **not** shown to be negligible for growth. `eps_0` is fixed by the cluster, and
     shifting `W` also enlarges `G`.
   * **Withdrawn:** "all probabilistic versions are the same statement up to `e^(O(d))`", "(R_H) *is* MIN^pair for
     the growth question", and "`e^(O(d))` shifts `W` by `O(d)`".
3. **Graded form and restriction ratio (proved).**
   * `Psi_A(u) = q_A(l_u^m mod J(A)_m)`, with `q_A` positive definite on `R[k']_m/J(A)_m`. The slab polynomials enter
     only through `J(A)`.
   * `AC_2` is equivalent to `rho_A(u) > mu(A) eps_0^2 d^(2m')/Psi_cube(u)` on `G_1`.
   * A sufficient condition is `rho_A > e^(-W/2 + O(m'(ln n + abar_n)))`, which is `e^(-(1/2 - o(1))W)`.
4. **Codimension (Thm 4.1, proved; Bombieri-Vaaler quoted, theorem number unconfirmed).**
   * `codim(V_H in V(A,m)) <= r_0 ln U/(2 ln H)`, with `r_0 = dim Poly_(<m)|_A`.
   * This is relatively `e^(-Omega(d))` only when `ln r = o(d)`. That holds in the growth range, since `r <= W/ln 5`
     and `W = O(n ln n)`, but **not** at the regime edge `ln 4r = d/316`.
   * It improves `theory.md` Lemma 1.1 by removing its `U ln 2/ln H` term. **No codimension bound gives a pointwise
     passage.**
5. **The passage `H = infinity -> V_H'`.**
   * (P1) (span plus `l^1` norm) is **proved for the uncut cube only**, with `K = sum_alpha 2^(n_2(alpha)) =
     e^(O(m ln d))`; the proof is the referee's.
   * For one cut and few cuts (P1) is unproved. There MIN^pair with `kappa_H` is proved directly instead.
   * (P1) and (P2) are open in the regime. Fakes may use `TC_H \ TC_V`, so the passage belongs to the barrier.
6. **Model cases (proved).**
   * Cube: exact `Psi`, and MIN^pair for `m <= theta d` with `theta < 1`.
   * Products: exact product formula. Disjoint signed cuts need sides that avoid the base coordinate.
   * One cut: MIN^pair for every `L`.
   * Few cuts: MIN^pair for `q <= (1/2) log_2(d/ln 4r) - 5.6` and `n >= n_0(C')`. This is asymptotic, not
     "explicit room".
   * Thin cuts: `rho_A >= 1/2` when `mu(A^c) <= 3^(-3m)/4` (O'Donnell Thm 10.21 constant quoted from memory). With
     the weaker secondary-source form the threshold is `120^(-4m)/4`, and the conclusions are unchanged.
   * **None applies to actual or gap-admissible regime data.** Prop 5.7 (proved): pair-zonotope certificates need
     `|M| >= m'/2 >= 2w` disjoint pairs, while `W > 8 ln C'` forces `|M| < 2w`. "Regime slabs are thick" is proved
     only for balanced sign columns; for general data it is evidence.
7. **Prop 7.1 (proved, exact numbers).**
   * At `(n, r) = (2901, n)`, `(5000, n)` and `(6000, n^2)`, a slab polynomial of degree `2j + 1` is nonzero on only
     `10^(-13.18)`, `10^(-13.92)` and `10^(-25.70)` of the ternary cube.
   * So bounds depending only on degree, variance and influence are vacuous at degree `>= 2j + 1`. MNV is vacuous
     for a different reason (rank `<= d/m`).
   * Scope: `Q_S` has leading form `l_(1_S)^(2j+1)`, and `1_S` lies outside `G`. So the proposition does **not**
     refute anti-concentration for `{l_u^(m') - P : u in G_1}`.
   * `n = 2901` is the record's exact `m_max >= 4w` threshold, not a regime size. The regime sizes are 5000 and
     6000.
8. **Minimal open statement (corrected).** `(R_H)`: `min_(G_1) Psi_(A,H') > d^(2m') eps_0^2`.
   * It follows from `(R_inf^lambda)`: `rho_A > e^(lambda) eps_0^2 d^(2m')/Psi_cube`, together with
     `(P2^lambda)`: `Psi_(A,H') >= e^(-lambda) Psi_A`, for **one matched** `lambda >= 0`.
   * **Withdrawn:** "`(R_H) <== (R_inf) + (P2)`". The `e^(-o(n ln n))` loss of (P2) is not absorbed.
   * The judge re-derived the matched form:
     `Psi_(A,H') >= e^(-lambda) rho_A Psi_cube/mu(A) > eps_0^2 d^(2m')/mu(A) >= eps_0^2 d^(2m')`.

### 2.2 antitheory (`antitheory.md`, `referee_antitheory.md`): proofs survive, headline lesson narrowed to finite `n`

The referee re-derived every proof and recomputed every table independently (exact central trinomials, 50-digit
arithmetic) without running the author's code. Everything reproduces once `theta_0` is corrected.

1. **Restatement (proved).**
   * `AC_A(u) <=> kappa_V(u) > eps_0`. The dual form: two probability measures on `A` with equal moments of order
     `< m'` and different `m'`-th moments along `u`.
   * `AC_T(u;theta) => AC_A(u)` for `theta <= mu(A)`.
   * Strictness of the density form is witnessed only **outside** the regime (uncut cube, `u = 1`). For regime
     boxes the converse is merely unproved.
   * Referee addition (proved, circuits plus Hadamard): `kappa_V = kappa_(H*)` with `ln H* = (E/2) ln E`. Together
     with Minkowski-Vaaler counting this leaves the `H`-gap (`ln H` between `2 ln U` and `(E/2) ln E`, against the
     rigidity tolerance `~ wW/2`) **open**, not "real".
2. **Degree barrier (Thm 3.1, proved, exact table).**
   * For `b = floor(d/m') >= 2` there are top-chaos degree-`m'` polynomials with equal influences `1/b` and rank
     `b`, and `mu{Q != 0} <= exp(-0.4 m'^(3/2)/sqrt d)`. The exact values run from `10^(-648.5)` (`n = 10^4`) to
     `10^(-2088)` (`n = 10^6`).
   * At `b = 1` (`n = 2900, 4096`) the example is a single monomial, and no degree-`m'` polynomial is regular at all.
   * Conclusion: no hypothesis of degree, influence, chaos or rank type gives anti-concentration at degree `m'`.
     Only the prescribed leading form can be used.
3. **Density forms on the averaged region (proved, finite `n`; Sylvester `(2^k, 2^k - 1)`, `W = 3n ln n`).**
   * `AC_T(.;1/2)` fails at explicit points of `{|u|_inf = 1, Gamma(u) <= abar}` for `k = 12..22` (author) and
     `k = 26, 30, 31` (referee). The construction stops working at `k = 32, 33, 34`.
   * `AC_T(.;theta_0)` with `theta_0 = 1 - (1/2)(4r)^(-1/2)` (the text's `1/4` is wrong; the code is right) fails for
     `k = 12..16`, and for `k = 17` at `W = 4n ln n`.
   * At the Sylvester-specific density `theta_S`, the exceptional bit strata are exactly the `<= 4`-column cheap
     range (`k <= 60`).
   * Nothing is shown at `W_min`, nothing on `G`, and nothing for `k >= 32`. The necessary condition `K o_q > d`
     dies at `k = 36` (`theta = 1/2`).
   * **Withdrawn:** "a true anti-concentration statement must see the slabs" and "a profile-free statement of the
     steelman type is false at regime sizes". Also withdrawn: "never at `theta_S`" (false for `k = 12, 13`, harmless).
   * Prop 4.2 must read `|X_i| <= N`, not `>= N`.
4. **Fibre restriction (Lemma 5.1, Cor 5.2) and grid model (Lemma 5.3): proved.** Fibre restriction is
   height-preserving. It must be applied with `eta = dist_inf(u, L)` and `dist_inf(v, L) = 1`, because the fibre
   quantity depends only on `v mod L`.
5. **ILO_m (conditional).**
   * ILO_m (a `kappa_V` statement) plus separation gives `MIN^pair_inf`, provided every projection-cheap `L` has an
     adequate collision count.
   * At exponent `m'` the room suffices from `n ~ 5.3 * 10^5` for 16-class spans, and from `n ~ 2.2 * 10^6` for the
     known 32-class spans `L_(2,1)` (referee ledger).
   * It gives no growth, because of the `H`-gap.
   * It quantifies over all regime profiles, fakes included. It is refuted by nondegenerate real **`TC_V`** points of
     regime fake boxes (not `TC_H` points, as the text says).
6. **Literature (WebSearch snippets only; primary sources unreachable).**
   * Confirmed to exist: MNV (ToC 2016), Kane (Ann. Probab. 2017), Kwan-Sauermann (Compositio 2025),
     Jin-Kwan-Sauermann-Wang (arXiv 2505.23335), Grebennikov-Kwan (arXiv 2505.24699).
   * The verbatim MNV quote is **unconfirmed** and should be labelled a paraphrase.
   * Grebennikov 2610.08708 is unverified.
   * Filmus-Mossel needs degree `o(n)` in the journal version but `o(sqrt n)` in the CCC version.
   * The §5.2 meta-claim about how anti-concentration theorems are proved is a heuristic reading.
   * No quoted theorem is a premise of a proved statement.

### 2.3 barrier (`barrier.md`, `referee_barrier.md`): proofs survive with slips, three headline claims corrected

Barrier scale: `W = r tau`, `tau = K ln n`. Fake data: `sdc_q` (`n = q + 1`, `r = 2q`, `k = e_1`). The referee's code
is independent.

1. **AC is box-free (Prop 2.1, proved, elementary).** `AC(u)` holds iff no half-cube `S` has
   `kappa_V^((S))(u) <= eps_0`.
   * It is meaningful only for `W >~ 2m' ln d + 5m'`. Below that, a one-variable Chebyshev witness refutes it for
     every `u` (referee; made rigorous and quantified by the judge in R16-J4).
   * It implies `MIN^pair_inf`.
   * **Withdrawn:** "AC is strictly stronger than MIN^pair" and "a proof of AC would decide the barrier question for
     all fake boxes". AC controls `kappa_V >= kappa_H`. It implies MIN^pair, and decides the `kappa_H` barrier, only
     together with `V_H = V(A,m')`, which is open.
2. **Thm 2.3 (proved): only light characters need Diophantine control.**
   * In fake Thm 2.1, only characters with `sigma_lambda < r/2 + O(n/tau)` need (2.2'), in the form
     `|lambda . u*| >= (13 s*/C) e^(-tau sigma/2)`.
   * Heavy relations, even exact ones, cost nothing. (The first heavy group is bounded by exactly `1/4`, which is
     harmless.)
   * This is a proved correction to the hypothesis of fake Thm 2.1: the barrier test is MIN^char restricted to light
     characters.
3. **Collision strata (Prop 2.4), corrected.**
   * Proved: transversal-vanishing certificates along a stratum carrying a light relation `lambda_0` cannot certify a
     barrier point when `2 o sigma_(lambda_0) < r - 4 ln d/tau`.
   * For `sdc_q`, `sigma_xy = 1` holds exactly for every pair (`q <= 83`, referee), so this covers every collision
     stratum.
   * **Withdrawn:** "collision strata never give barrier points". Near-zeros near collision strata coming from other
     mechanisms (other components of `TC_H`, vanishing transversal coefficients) are not excluded. The argument
     should also use the necessary character inequalities, not (2.2').
4. **Collision-free rational strata (Prop 2.5, proved after correction).**
   * A rational `L ⊂ TC_H` with no pair relation, light relations of slack `>= c_0 r`, and transversal order
     `o > 1/(2c_0)` yields barrier-scale fakes. For `sdc_q` that means `o >= 3`.
   * Corrections: `tau_L` gains a term of about `ln(2r + 3)/c_0`, and `rho_1 >= 1/2` must be assumed.
5. **Products and sections (Prop 3.1, proved).**
   * `nu` is additive on products; this is a self-contained proof.
   * **Correction:** the regime bound is `m' <= 8j + 8`, not `8j + 4`. The `8j + 4` form fails for 4446 regime
     primes, the first being `q = 3271`; the judge confirms (R16-J2). The conclusion survives for `j >= 4`.
6. **Cyclic symmetry (Thm 3.2, Cor 3.3, proved; strengthened by the referee).**
   * The averaging identity holds for `m < N` when `N` is odd, not only for `2m < N`: `|f_c|^2` has only even
     frequencies. The judge re-checked this numerically (R16-J3).
   * This matters because the claim "`2m' < q` throughout the regime" is **false**. The first regime `q = 3251` has
     `m' = 2057`, and `2m' < q` first holds at `q = 5479`. But `m' < q` holds at every regime `q` (R16-J2).
   * Referee addition (exact for `q <= 83`, two-line proof): `sdc_q` is invariant under all of `AGL(1,q)`. Hence
     `kappa_H` is **identical** on all Fourier circles `C_a`, and `sigma_xy` is constant.
7. **Circle test (Thm 3.5, proved modulo H_light; interpretation corrected).**
   * Suppose H_light(`c_0`) holds (every relevant non-pair light character has `sigma >= c_0 r`). If
     `inf_(C_1) kappa_H <= eps_0/(2 sqrt(m'+1) e^(pi^2/(4q)))` at some regime `q`, then route 5 is capped.
   * Evidence for H_light: the referee's global MILP (HiGHS-certified optimal) gives min non-pair `sigma/r` of
     `0.184` to `0.227` for `q <= 43`. H_light is **not proved**.
   * **Withdrawn:** "the circle case is one explicit algebraic question". Only "yes" is decisive. An exact "circle
     not in TC" leaves the quantitative test open. Liouville at the algebraic circle points gives only
     `e^(-Theta(q^2))`, against the needed `e^(-Theta(q log q))`.
   * The exact circle data are exact-zero statements at toy `(L, m)`, not at regime `(j, m')`:
     * `sdc7` (`L = 2`, `m = 8`): `dim V = 1`, and `Phi_c = -7717.5` is constant on every circle (referee,
       independent);
     * `sdc11` at `m = 8` and `sdc19` at `m = 5, 6`: author only, mod-`p` or floating point.
8. **Also corrected or withdrawn.**
   * The Sylvester aside "MIN^pair would fail at every `W`" is unsupported (Sylvester pair margins are zero).
   * The §5 "iff" is only a sufficient trichotomy.
   * The fractional character `(2/(q+1))(1_(F_q) - q e_inf)` has `sigma = r/2` and vanishes identically on every
     circle. It is harmless only because it is irrelevant for `k = e_1`.
   * Barber arXiv 2104.15057 concerns **five** roots of unity.
9. **Net.** No barrier point was constructed, and none was excluded. The referee also tried four further candidate
   sources: multiplicative-character circles, subgroup fixed spaces, the integer ramp, and cheap-plus-perturbation.
   None yields a barrier point.

## 3. Status of AC / MIN^pair

| statement | status | source |
| --- | --- | --- |
| **MIN^pair** (`kappa_H`, cluster's own `G`), any regime box: actual, generic or fake | **open**. Not proved, not refuted, not capped. No evidence either way at regime size | all three |
| `(R_H)`: `min_(G_1) Psi_(A,H') > d^(2m') eps_0^2` | **open**. It implies MIN^pair losslessly, and is implied by it at level `eps_0 3^(-d/2)`. It is the barrier question by fake Thm 2.1 | reduction §8, corrected |
| `MIN^pair_inf` = `AC_A` = `AC_2` (sup and `L^2` forms, up to `U <= 3^d` one way) | **open** | reduction Thm 2.1; antitheory Prop 1.2 |
| `(R_inf^lambda) + (P2^lambda) => (R_H)` with matched `lambda` | **proved** implication. Both hypotheses open | referee_reduction §4 |
| `(R_inf) + (P2)` as stated in reduction §8 `=> (R_H)` | **false as stated** (quantifier) | referee_reduction M2 |
| "Most `k`" forms `AC_p` (`p > 1/U` fixed) and cube form `AC_cube` | **sufficient, not equivalent**. The converse is false in an exact toy and open for regime boxes | referee_reduction §2; R16-J1 |
| Steelman `AC = AC_T(.;1/2)` | **strictly stronger form** (toy; uncut cube) | referee_antitheory §1.1 |
| ... on the averaged region `{Gamma <= abar}`, Sylvester, `W = 3n ln n` | **false** at explicit points for `k = 12..22, 26, 30, 31` (finite `n`). The construction fails at `k = 32..34` | antitheory Thm 4.4; referee §2.3 |
| ... for every `u`, at `W <= 4m' ln sqrt(3d) + 4 ln(d/(2(m'+1)))` | **false** (rigorous one-variable Chebyshev witness). With `r = n` this covers the smallest admissible `W = 3n ln n - 14n` for every regime `n <= 3750` | R16-J4 (judge, unrefereed) |
| ... on `G`, asymptotically, or at `W_min` | **open**. Nothing is shown | — |
| `AC_T(.;theta_0)` (profile-free density) | **false** for `k = 12..16` (and `17` at `W = 4n ln n`), finite `n` | antitheory Thm 4.4 |
| Discrete anti-concentration from degree / variance / influence / rank hypotheses | **provably vacuous** at degree `>= 2j + 1` | reduction Prop 7.1; antitheory Thm 3.1 |
| Inverse Littlewood-Offord theorems (Kwan-Sauermann, JKSW, Grebennikov-Kwan) | right genre, but **bounded degree** only | both literature sections |
| MIN^pair on model boxes: cube, one cut, few cuts, thin cuts, products | **proved**. **Inapplicable** to actual and gap-admissible regime data (Prop 5.7) | reduction §5 |
| Passage `kappa_V -> kappa_H'`: codimension | **proved**: relative codim `e^(-Omega(d))` when `ln r = o(d)` | reduction Thm 4.1 |
| Passage, pointwise (P1)/(P2) | **open** in the regime. (P1) is proved for the uncut cube only | referee_reduction §4 |
| `H`-gap | **open**: `kappa_V = kappa_(H*)`, `ln H* = (E/2) ln E`, against a tolerance of `~ wW/2` | referee_antitheory §1.2 |
| Barrier points (nondegenerate near-zeros in `G`) of regime fake boxes | **open**. None constructed, none excluded. Collision strata: the transversal ceiling cannot certify them. Collision-free rational strata with `o >= 3`, or one Fourier circle at the quantitative level, would cap route 5 (H_light assumed) | barrier §2-§3, corrected |

**Reading.** Round 16 settles what kind of statement MIN^pair is: an `L^2` restriction statement about the top
component of one fixed power, on one box, at an absolute scale. It is not an anti-concentration statement, and no
probabilistic weakening of it is known to be equivalent. Every positive tool located, from the literature or from
the record, is either blind to `u in G` (vacuous) or needs twin classes or thin slabs, which regime data lack. The
barrier side has one new, sharp, symmetric candidate (the Fourier circles of `sdc_q`) and one proved simplification
(only light characters matter).

## 4. Judge's own checks (unrefereed)

All scripts are in `judge_checks/`, and each writes `*_out.txt`. **None of this has been refereed.**

* **R16-J1 (exact, sympy rationals): the referee's toy against probabilistic equivalence.**
  `j1_probforms.py` uses independent code: exact projection residual onto the column space of the degree-`< m`
  evaluation matrix. It reproduces the referee's toy:
  * one-cut `A` (`d = 2`, 7 points), `m = 3`: symbolically `Psi_A^Sigma(u) = 6 u_1^2 u_2^2 (u_1 - u_2)^2`. The zero
    set is exactly the braid arrangement (with `u_3 = 0`), and the value at `(1/2, 1)` is `3/8`;
  * `B` = cube minus `{(1,1), (-1,0)}`: `Psi_B(1/2, 1) = 0`;
  * `A` minus `(-1,0)`: `Psi = 0`;
  * cube at `m = 4`: `Psi^Sigma = 16 u_1^4 u_2^4`, which is 1 at `(1/2, 1)`;
  * cube minus a point at `m = 4`: `Psi = 0` identically.

  So `MIN^pair_inf` holds while `AC_cube` and `AC_p` (`p > 1/7`) fail. The probabilistic-equivalence claim is false
  in general.
* **R16-J2 (50-digit arithmetic): sdc regime arithmetic.** `j2_regime.py` covers primes `q = 3 mod 4` below
  `2 * 10^5`, with `d = q`, `r = 2q`:
  * 8777 regime primes; the first is `q = 3251`, with `w = 514.20`, `j = 257`, `m' = 2057`;
  * `2m' < q` first holds at `q = 5479`;
  * `m' < q` holds at **every** regime prime;
  * `m' > 8j + 4` at 4446 of them, and `max(m' - 8j) = 8`.

  This confirms referee_barrier on all three points. With the referee's `m < N` relaxation, Thm 3.2, Cor 3.3 and
  Thm 3.5 apply at every regime `q`. Sylvester's regime starts at `k = 12` (`m' = 2256`).
* **R16-J3 (50-digit numerics; the statement has a two-line proof): the `m < N` relaxation of Thm 3.2.**
  * Proof: `f` real and homogeneous of degree `m` in `(cos t, sin t)` makes `f^2` a trigonometric polynomial with
    frequencies `2l`, `|l| <= m`. The sum over the `N` equispaced points kills frequency `2l` unless `N | 2la`, i.e.
    `N | l` (`N` odd, `gcd(a, N) = 1`), which forces `l = 0` when `m < N`.
  * Check (`j3_orbitavg.py`, random integer coefficients, step `a = 2`): the orbit average divided by the mean is
    exactly 1 for `N = 7`, `m <= 6` and `N = 11`, `m <= 10`. It differs from 1 at `m = N` and `m = N + 1`.
* **R16-J4 (rigorous, 50-digit arithmetic): the steelman's cube form is false for every `u` at the regime start.**
  * Claim. Let `|u|_inf = 1`, `X = <k',u>`, `k'` uniform on `{-1,0,1}^d`.
    1. `Var X = (2/3)|u|_2^2 <= 2d/3`, so `Pr(|X| < s) >= 1/2` for `s = sqrt(4d/3)` (Chebyshev's inequality).
    2. `Q = 2^(1-m') s^(m') T_(m')(X/s)` has leading form `l_u^(m')` and satisfies `|Q| <= 2^(1-m') s^(m')` on
       `{|X| < s}`.
    3. Hence `AC_T(u;1/2)` fails for **every** such `u` whenever
       `W <= W_rig := 4m' ln sqrt(3d) + 4 ln(C d/(2(m'+1)))`.
  * Numbers (`j4b_chebscale_rig.py`, `r = n`, `C = 1`). The record's smallest admissible `W` is
    `3n ln n - 14n` (two_thirds; README-13 J5). It lies **below** `W_rig` for every `n` in `[2966, 3750]`, e.g.
    `29615 < 34317` at `n = 2966`.
  * Asymptotically `W_rig = O(m' ln d) = o(n ln n)`, so this is a finite-`n` effect only.
  * It concerns the cube form only. The sup form `AC_A` (MIN^pair) is untouched: on `A` the same witness gives only
    `|Q| <= 2(d/2)^(m')`.
  * Consequence: at the regime start, near-extremal `W` makes the steelman's `AC` false everywhere, and not only at
    structured points. This is a second, independent reason (after the referee's toy) to drop density forms as
    targets.
  * Using the Gaussian median `0.55|u|_2` in place of Chebyshev's inequality (evidence-level constant;
    `j4_chebscale.py`) extends the range to `n <= 4741`.

## 5. What it means for route 5 and the growth rate

**Growth rate: unchanged.** `M <= (2/3 + eps) log R/loglog R` remains the best bound for actual circles. Nothing this
round improves the growth rate or touches the uniform bound, unconditionally or from a standard conjecture.

### 5.1 Route 5: still open, with a precise missing lemma; not capped; not advanced

* **Missing lemma.** It is `(R_H)` of §0, or equivalently (up to the one-way factor `3^(d/2)`) MIN^pair itself. The
  lossless sufficient split, with one matched `lambda`, is:
  * `(R_inf^lambda)`: the top component of `l_u^(m')` keeps more than an `e^(lambda) eps_0^2 d^(2m')/Psi_cube(u)`
    fraction of its cube `L^2` mass on the slab event, at every `u in G_1`. This is `e^(-(1/2 - o(1))W + lambda)`.
  * `(P2^lambda)`: `Psi_(A,H')(u) >= e^(-lambda) Psi_A(u)` on `G_1`.
* **Not capped.** No nondegenerate near-zero of `kappa_H` in `G` is known for any regime fake box.
* **Not advanced.** No new route-5 input was proved for actual boxes.
  * Every proved positive case needs twin classes or thin slabs, and neither actual nor gap-admissible regime data
    have them.
  * The literature is provably vacuous at our degree for any tool that does not see `u in G`.
* **Arithmetic is spent, and the barrier question is unchanged in substance.** `(R_H)` depends only on
  `(J(A), V_H', G)`. A proof must therefore separate actual from fake pattern boxes through their algebra alone.
  * Judge remark (heuristic, not used anywhere): a sign-column profile with `e_p = 1` is realised by lattice points
    on a circle, by choosing `pi_p` or its conjugate per point. What actual clusters add is the short-arc
    condition, i.e. the angles, and Step 0 has already consumed it.
  * So a box-level proof can only work if **no** admissible regime box, fake ones included, has a near-zero in `G`.
    That is the barrier question itself.

### 5.2 What the round did achieve (barrier-map sharpening only; no growth credit)

1. **The target's type is settled.**
   * MIN^pair is a sup or `L^2` statement about a fixed top component.
   * The density and anti-concentration forms are sufficient but possibly strictly stronger (exact toy).
   * The steelman's cube form fails at finite regime sizes, at structured points (antitheory) and, at near-extremal
     `W` for `n <= 3750`, everywhere (R16-J4).
   * This closes local_smoothing S3 as a route: the "easy continuous half" is the identity `Psi = m'! |u|^(2m')`, and
     the discrete half is the whole problem.
2. **The literature verdict is now a theorem.** Exact degree-barrier examples (Prop 7.1, Thm 3.1) show that no
   bound depending only on degree, variance, influence or rank can work. Only the prescribed leading form, with
   `u in G`, can be used.
3. **The `H = infinity -> V_H` passage is localised.**
   * Codimension is no obstacle (Thm 4.1).
   * The pointwise passage is open.
   * `kappa_V = kappa_(H*)` with `ln H* = (E/2) ln E` bounds the `H`-gap from the other side.
4. **The barrier test is sharpened.**
   * Only light characters need the Diophantine condition (Thm 2.3).
   * For `sdc_q`, collision strata cannot be certified by transversal ceilings (Prop 2.4, corrected).
   * Collision-free rational strata with `o >= 3` would cap route 5 (Prop 2.5).
   * The `AGL(1,q)`-symmetry makes `kappa_H` constant on all Fourier circles, and Thm 3.5 turns a single quantitative
     near-zero on one circle into a cap, modulo H_light.

### 5.3 Conditional theorems on record after round 16 (changes against round 13 only)

| hypothesis | conclusion | implication | hypothesis status |
| --- | --- | --- | --- |
| `(R_H)`, cluster's own regime box | `n = O(loglog R)`, all norms | proved (via MIN^pair, R12) | **open** |
| `(R_inf^lambda) + (P2^lambda)`, matched `lambda` | `(R_H)` | proved (referee_reduction §4) | both **open**; neither has a regime instance |
| `AC_2`, `AC_p`, `AC_cube`, `AC_T(.;1/2)` | `MIN^pair_inf` only (not MIN^pair: `H`-gap) | proved | `AC_T(.;1/2)` **false** at finite `n` (antitheory Thm 4.4; R16-J4); others **open** |
| ILO_m + adequate collision counts for every projection-cheap `L` | `MIN^pair_inf` once the room covers every such `L` (`n >~ 5.3 * 10^5` for 16-class spans, `>~ 2.2 * 10^6` for the known 32-class spans) | proved (antitheory §5.4, corrected) | **open**; refuted by any regime fake box with a nondegenerate real `TC_V` point |
| H_light(`c_0`) + `inf_(C_1) kappa_H <= eps_0/(2 sqrt(m'+1) e^(pi^2/(4q)))` for one regime `sdc_q` | **route 5 capped** (MIN^pair fails on gap-admissible fakes) | proved (barrier Thm 3.5 with the `m < N` relaxation) | H_light: **evidence** only (`c_0 ~ 0.184`, `q <= 43`); circle near-zero: **open** |
| collision-free rational `L ⊂ TC_H` of a regime fake box, transversal order `o > 1/(2c_0)` | route 5 capped | proved (barrier Prop 2.5, corrected) | **open**; no candidate known |

**Rubric (b), fake-violated input: unchanged.** The input that fakes violate is still rigidity at the actual point,
spent in Step 0. The round shows from three sides (reduction §8, antitheory §5.4, barrier Prop 2.1) that the
residual is a box statement, which fakes are entitled to test.

**Verdict on route 5.** Still open, with the precise missing lemma `(R_H)`. Not capped. Not advanced on the growth
side. The anti-concentration handle proposed by the steelman is closed: it does not unblock route 5, and its
probabilistic forms are partly false. Route 5 should stay **paused** for growth purposes. The only route-5 work
with a decisive possible outcome is the barrier test of §7.

## 6. Scores

Rubric of `round10/README.md` §1:
* (a) correctness after refereeing;
* (b) fake-violated input;
* (c) movement of growth or uniformity.

The overall score is a judgement, and it weighs (b) and (c) above (a). No write-up moves the growth rate, so no (c)
score exceeds 1. Partial (c) credit is only for proved sharpenings of the route-5 map.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **reduction** | 6 | 1 | 1 | **4** | Correct core (Thm 2.1, Prop 3.1/Cor 3.2, Thm 4.1, §5 model theorems, Prop 7.1), independently re-derived and exactly re-checked. It turns the steelman's reframing into an exact `L^2` restriction statement, `(R_H)`, and proves the literature vacuity with exact numbers. The codimension theorem improves theory.md L1.1. Deductions: two major framing errors (the false "all probabilistic versions are equivalent", refuted by an exact toy; the quantifier error in the §8 split); "negligible `e^(O(d))`" heuristics presented as results; (P1) claimed beyond the cube; the hypercontractivity constant unconfirmed. |
| **barrier** | 5 | 2 | 1 | **4** | The round's only new handle on the decisive question. Thm 2.3 is a proved simplification of the fake-test hypothesis (light characters only). The cyclic-symmetry dichotomy plus the referee's `AGL(1,q)` invariance yields a single, symmetric, quantitative cap test (Thm 3.5). (b) = 2 because it pins down exactly what a fake must exhibit. Deductions: three logically overstated headline claims (AC vs MIN^pair, "collision strata never", "one algebraic question"); two false regime-arithmetic statements (`2m' < q`, `m' <= 8j + 4`), repaired by the referee; H_light unproved. |
| **antitheory** | 6 | 1 | 1 | **3** | Solid proofs: the exact dual form, the degree-barrier Thm 3.1 (stronger than Prop 7.1 when `b >= 2`), the finite-`n` refutation of the steelman's density form, fibre restriction, and the grid model. The tables reproduce. Deductions: the headline lesson ("must see the slabs"; "false at regime sizes") is a finite-`n` phenomenon of one family and fails from `k = 32`; ILO_m's fake test names the wrong cone; several minor slips (`theta_0`, Prop 4.2, "never at `theta_S`", "maximally regular" at `b = 1`). Its main message overlaps with reduction §7. |

## 7. Single most promising next step (whole record)

**Recommended: the quantitative Fourier-circle test on the regime fake boxes `sdc_q`, starting with a proof of
H_light.**

### 7.1 Statement to decide

> **CIRCLE(q).** Take `q = 3 mod 4` prime with `q >= 3251` (regime), the `sdc_q` box `A_q` at `T = 1`, `j`, `m'`, and
> the barrier scale `W = 2q tau` with `tau = K ln n >= tau_D` (barrier Lemma 3.4, Thm 2.3). Decide whether
>
> ```text
> inf_theta kappa_H(hat u_theta) <= eps_0 / (2 sqrt(m'+1) e^(pi^2/(4q))),   u_theta(x) = cos(theta + 2 pi x/q), u_inf = 0.
> ```
>
> **Yes** at one regime `q`, together with H_light, means route 5 is **capped** (barrier Thm 3.5, with the
> `m < N` relaxation, which covers every regime `q` by R16-J2).
>
> **No** needs a lower bound for `kappa_H` on `C_1`. By `AGL(1,q)`-invariance, this holds on all circles at once.

* For the cap direction `H` is irrelevant, since `kappa_H <= kappa_V`. It suffices to exhibit one polynomial `P` of
  degree `< m'` with `|l_(u_theta)^(m') - P| <= d^(m') eps_0/(2 sqrt(m'+1) e^(pi^2/(4q)))` on all of `A_q`.
* For the other direction the height matters, and one needs a certificate `c in V_H'(A_q, m')`.

### 7.2 Why this step, from the whole record

1. **Every growth route is blocked at a named statement.**
   * Vojta: conditional, unchanged.
   * Route 3 / P3-all\*: by round 15 (`rigSLtheory` §4.2), no outcome of Rig^SL moves the growth rate by itself. The
     growth inputs there are H^gen + L1, statements about genuine Gaussian integers with no new handle since round 9.
   * Route 5: the residual `(R_H)` is the barrier question (§5.1). Every analytic, algebraic and literature handle
     proposed for it has been closed: Lojasiewicz (R12), locality (R13), product and zonotope certificates (R11,
     R13, R14 pal11), and now anti-concentration (R16).
2. **CIRCLE is the only route-5 question in the record that is decisive at regime size in one direction, by theory
   rather than by toy computation.**
   * Round 13 §3 found that "neither direction of the barrier test has a viable method". CIRCLE is the first test
     that is one explicit object, with a symmetry group of order `q(q - 1)`.
   * The Diophantine side is already handled by Lemma 3.4 (given H_light), and "somewhere" equals "everywhere" on the
     circle (Thm 3.2).
3. **Both outcomes produce something the record lacks.**
   * *Yes:* route 5, in its MIN^pair form, is capped at the barrier scale. That completes the barrier map for the
     determinant method at `T = 1` and redirects the program to the genuine-prime inputs (H^gen, L1) and Vojta.
   * *No, with a proof:* this is the first non-zonotopal lower-bound certificate at an irrational, collision-free
     direction of an overdetermined regime box. That is precisely the tool missing for `(R_inf)` (antitheory §5.2
     gap; barrier open point 4; round 13 §3). The circle, with its full affine symmetry, is the most favourable place
     to build it.
4. **It is cheap to start.** H_light is a character-sum statement about Paley patterns, and the first exact
   calibration is feasible with the translation-isotypic reduction.

### 7.3 First concrete tasks

1. **Prove H_light(`c_0`) for `sdc_q`** (theory). The statement: every relevant integral non-pair `lambda` that is
   light has `sigma_lambda >= c_0 r`.
   * Large norms are automatic from `sigma_lambda >= rho_1 ||lambda||_1` (barrier Lemma 3.4's proof).
   * Small norms should follow from Weil-type bounds for quadratic-residue patterns. Barrier §3.3 states the heuristic
     `sigma_lambda = r(E|lambda(S)| - 1/2) + O(||lambda||_1 r/sqrt q)` [M]; this needs checking.
   * The intermediate range `sqrt q <~ ||lambda||_1 <~ 0.4 r` is the real content.
   * Evidence so far: the referee's global MILP gives `c_0 ~ 0.184` for `q <= 43`.
2. **Cap attempt (upper bound).**
   * Use the slab coordinates `s_a = k(a + QR)`, with `|s_a| <= j` on `A`. The Fourier coefficient `hat k(1)` is
     `hat s(1)` divided by a Gauss-sum factor of modulus about `sqrt q/2` [standard identity; verify for the exact
     `sdc` patterns].
   * Look for degree-`< m'` approximants of `Re(e^(i theta) hat k(1))^(m')` on `A` that use `k_x^3 = k_x` and
     `|s_a| <= j` jointly.
   * Any success at one `theta`, at the stated scale, caps route 5 (given H_light).
3. **Survival attempt (lower bound).**
   * Look for an `AGL(1,q)`-equivariant dual certificate `c in V_H'(A_q, m')`. It lies in the translation eigenspace
     of character `zeta^(m' - 2i)` matching the monomials `l_(psi_1)^i l_(psi_(-1))^(m'-i)` (Cor 3.3).
   * It needs `|c|_1 <= e^(o(n ln n))` and `|Phi_c(u_theta)| >= |c|_1 d^(m') eps_0/(2 sqrt(m'+1))`.
4. **Calibration (exact; evidence only; state `(q, L, m, H_A(m))` for every number).**
   * Models: `sdc_q` with `q in {19, 23, 31, 43}`, `L in {1, 2}`, every `m < q` with `dim V(A,m) >= 2`.
   * Compute `Psi_A` on `C_1` exactly in the `Z_q`-isotypic component, and compare it with `Psi_cube(C_1)` and with
     `Psi_A` at random points of `G_1`.
   * Certify non-membership with `rank_p = rank_Q` (README-13 conventions).
   * Report any overdetermined degree with `C_1 ⊂ TC_V`. That would be the first such point on an overdetermined
     fake box.

### 7.4 Reading of the outcomes

* **Cap (yes).** Record route 5 (MIN^pair and `(R_H)` forms) as capped at the barrier scale by `sdc_q` fakes.
  * The only route-5 variant left would use the actual point itself: a lower bound for degree-`m'` forms evaluated
    at the vector of Gaussian-prime angles.
  * Judge remark (heuristic, unverified, from memory): that is a measure of algebraic independence for many
    logarithms, which known transcendence theory does not supply at `~ n` variables.
* **No, with a proof.**
  * Next, test collision-free rational strata with `o >= 3` (barrier Prop 2.5), and `pal_q`, which lacks translation
    symmetry.
  * Then try to transplant the equivariant certificate to actual-type Sylvester boxes through the Singer cycle.
* **Undecided after one round.** Stop route 5 for good. Its residual is then the barrier question with no method in
  either direction.

**Not recommended:**
* discrete anti-concentration and inverse Littlewood-Offord literature at degree `m'` (provably vacuous for
  hypothesis-type bounds: Prop 7.1, Thm 3.1);
* density or cube forms `AC_T` (false at the regime start; R16-J4, antitheory Thm 4.4);
* further reformulations of `(R_H)`;
* Loc;
* LP-subcritical toys as regime evidence, beyond the calibration in task 4;
* Rig^SL as a growth step (round 15: it changes only the P3-all\* barrier map).

## 8. Records to correct, and files

### 8.1 Records to correct

* **`local_smoothing/steelman.md` S3.**
  * `AC(u)` (the cube form `AC_T(.;1/2)`) is a sufficient condition, and it is strictly stronger as a form (exact
    toy, R16-J1). It is false for every `u` at the regime start at near-extremal `W` (R16-J4), and at structured
    points for Sylvester `n = 2^12..2^31` (antitheory Thm 4.4, finite `n`).
  * The statement MIN^pair needs is the sup / `L^2` form `(R_H)`.
  * "sdc7 and sdc11 ... so `AC` fails there exactly": those toys have density `< 1/2`, so `AC`'s premise
    `|A| >= |T|/2` does not hold there.
  * "Any proof of `AC` must decide the barrier question": only its `H = infinity` (`kappa_V`) version.
* **`round12` `fake.md` Thm 2.1.** The Diophantine condition (2.2) is needed only for light characters (barrier
  Thm 2.3, refereed).
* **`theory.md` Lemma 1.1.** The codimension of `V_H` in `V` is at most `r_0 ln U/(2 ln H)`; the `U ln 2/ln H` term
  is not needed (reduction Thm 4.1; Bombieri-Vaaler quoted).
* **Top-level `README.md`, growth-rate status.** Add: round 16 shows that MIN^pair is an `L^2` restriction statement,
  not an anti-concentration statement. Discrete anti-concentration tools are provably vacuous at its degree. The
  barrier test is sharpened to the `sdc_q` Fourier circles (a cap if one quantitative near-zero exists, modulo
  H_light). There is no growth movement, and route 5 is still open and not capped.

### 8.2 Files (all under `round16/`)

| file | content |
| --- | --- |
| `reduction.md`, `referee_reduction.md` | exact reduction of MIN^pair to `Cheb_A` / `Psi_A`; graded form; Thm 4.1; model cases; Prop 7.1; report with the exact toy, the (P1) cube proof, and the matched-margin fix |
| `reduction_checks/`, `referee_reduction_checks/` | `lib.py`, `check_cube_product.py`, `check_onecut_zero.py`, `cond_lambda.py` (evidence), `slab_concentration.py`; `rlib.py`, `check_cube_product_ref.py`, `check_prob_forms.py`, `check_slab_conc_ref.py` |
| `antitheory.md`, `antitheory_notes.md`, `referee_antitheory.md` | literature survey, degree barrier, density-form failures, fibre restriction, ILO_m; report with the `k = 26..34` extension and the circuit bound |
| `antitheory_checks/`, `referee_antitheory_checks/` | author's scripts; referee's `tables.py`, `bits.py`, `necK.py`, `ledgerQ.py` |
| `barrier.md`, `barrier_notes.md`, `referee_barrier.md` | box-free AC, light characters, rational strata, cyclic symmetry and circle test; report with `AGL(1,q)` invariance, the H_light MILP, and the `m < N` relaxation |
| `barrier_checks/`, `referee_barrier_checks/` | author's `sdc_circles.py`, `sdc_fourier.py` (sdc7, sdc11, sdc19), `m8_tc.py`, `m8_refine*.py` (M8), `hlight.py`, `sigma_pairs.py`; referee's `sdc_sym.py`, `hlight_ref.py`, `hlight_milp.py`, `thm32_check.py`, `sdc7_octic.py`, `regime_arith.py` |
| `judge_checks/` | `j1_probforms.py` (exact), `j2_regime.py`, `j3_orbitavg.py`, `j4_chebscale.py` (evidence constant), `j4b_chebscale_rig.py` (rigorous), each with `*_out.txt` |
| `notes.md` | round working notes |

**Repository note.** Only files under 1 MiB are committed; pickles, arrays, caches and compiled binaries are
not. The scripts next to them regenerate them.
