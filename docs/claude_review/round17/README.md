# Round 17 - judge's report (CIRCLE(q) test)

## 0. Verdict

**No growth improvement and no uniform bound** were proved this round, unconditionally or from a standard
conjecture. `M <= (2/3 + eps) log R/loglog R` (`two_thirds/proof.md`) is still the best bound for actual circles.

* **CIRCLE(q): undecided** at every regime `q`. All four attacks agree, and all four referees uphold that verdict.
  * No cap (approximant or `TC_{H'}` membership) was built.
  * No quantitative lower bound for `kappa_H` on `C_1` was proved at any regime `q`.
* **H_light: proved** for `sdc_q`, every prime `q = 3 mod 4`, `q >= 7` (hlight Thm 2.8, two independent referee
  passes):
  * `c_0 = 1/6` for every such `q`, by elementary means;
  * `c_0 >= 0.2349` at every regime `q >= 3251`, using Weil's bound for 4-term characters;
  * the constant is sharp up to `O(q^(-1/2))`: dual pairs have `sigma = (q-1)/2` exactly.

  So barrier Thm 3.5, as corrected, is an **unconditional** cap test.
* **Route 5: STOPS.** This follows the round-16 rule (R16 §7.4): CIRCLE(q) is undecided after one round. Its residual
  `(R_H)` is the barrier question, and the record has no method for it in either direction.
  * It is **not capped**: no barrier point of a regime fake box is known.
  * It is **not unblocked**: no lower-bound certificate exists at regime size.
* **What the round did establish** (the barrier side only; no growth credit):
  1. H_light, which removes the last hypothesis from the cap test.
  2. survive Thm 3.1 / 3.4: on the Fourier circles of `sdc_q`, `kappa_H` has **no small positive values**. It is
     either identically 0 or at least `2^(-m)/(|A| H q^m)`, and any integral witness gives a Liouville-free,
     `theta`-uniform lower bound. The quantitative test therefore becomes two **exact** membership questions, with
     an undecided **height window** of `ln H` in `[3.5*10^5, 7.7*10^7]` at `q = 3251` (J1).
  3. An exact `L^2` splitting. `Psi_A` is constant on circles and is a positive combination of `m'+1` rational
     numbers. This was found independently by cap (Prop 4.1/4.2) and calib (Prop 2.1).
  4. Evidence (`kappa_V` only, toy boxes, model-dependent extrapolation) that the circle is far from `TC_V`.
     This points to "no" for `kappa_V`. It says nothing about `kappa_H` at the regime height.

**Scores** (rubric of `round10/README.md` §1, each out of 10; details in §5):

| write-up | (a) | (b) | (c) | **overall** |
| --- | --- | --- | --- | --- |
| hlight | 9 | 1 | 1 | **5** |
| survive | 7 | 2 | 1 | **5** |
| cap | 5 | 1 | 1 | **3** |
| calib | 6 | 1 | 1 | **3** |

**Recommendation (§7).** Stop route 5, and do not reopen the box question on toy evidence. The program has no
growth route with a credible method. If it continues, the single most promising target is the Vojta-type
arithmetic lemma `(L_8)` (Ptolemy route), treated as a self-contained problem about integer `2 x 8` matrices. It is
the only open statement in the record that has all three of these properties:
* it implies the uniform theorem;
* a standard conjecture predicts it is true;
* fake circles do not test it.

Honest chance that a further round proves it: well below 1%. The first task is a cheap, evidence-only census of its
small-height configurations. In parallel, the proved results should be written up for outside experts.

## 1. What was asked; conventions

**CIRCLE(q)** (R16 §7.1):
* *Setting.* `q` is a prime, `q = 3 mod 4`, with `q >= 3251` (the regime: `n - 1 = q >= 316 ln 4r`, `r = 2q`).
  * `w = sqrt(8 q ln 8q)`, `j = floor(w/2)`, `m' = ceil(4w) < q`.
  * The box is the `sdc_q` box `A` at `T = 1`, at the barrier scale `W = 2q tau`, `tau >= tau_D`, `tau = O(log q)`,
    together with the thresholds of barrier Thm 2.3.
* *Question.* Decide whether `inf_theta kappa_H(hat u_theta) <= T* := eps_0/(2 sqrt(m'+1) e^(pi^2/(4q)))` with
  `H = H' = e^(m'W/40)`.
* *Outcomes.*
  * **Yes**, together with H_light: route 5 is capped.
  * **No**, proved: the first non-zonotopal certificate at an irrational direction of an overdetermined regime box.
  * **Undecided:** route 5 stops.

**Labels.**
* *Proved*: all quantifiers, survived refereeing.
* *Exact*: integer or rational computation.
* *Evidence*: anything else.
* *J-n*: the judge's own checks, **unrefereed**.

`C = 1` wherever numbers are given. The first regime prime is `q = 3251`, with `j = 257` and `m' = 2057`
(`m'/q = 0.633`). **No exact model in this round is a regime box.**

**Round structure.** There were four attacks with adversarial referees:
* `hlight` had two independent referee passes;
* the other three attacks had one pass each.

Two container restarts killed running agents. The results of the first six agents were recovered from
`cached_results.json`, and `calib` and its referee were rerun. I read every write-up and report in full where it
carries a verdict, and the cached summaries everywhere else.

## 2. Referee-corrected statements
### 2.1 hlight (`hlight.md`; `referee_hlight_pass1.md`, `referee_hlight_pass2.md`): survives, minor corrections only

Both passes re-derived every step by hand. Their exact code is independent of the author's and of each other's.
* Pass 2 rebuilt `sdc_q` from the record's Paley-I definition (`exact.md` Thm 4.1). It matches the author's patterns
  up to complement, which leaves `sigma` unchanged.
* Pass 2 also proved, with a Smith-normal-form check, that `L_T/Z^n_0 = Z/((q+1)/2)`. So the relevant characters
  are exactly the integral ones, and the theorem covers all of H_light's domain.

**Theorem (H_light for `sdc_q`; proved).** Let `q = 3 mod 4` be prime, `q >= 7`, `r = 2q`. Let `lambda` be integral
on `F_q ∪ {inf}`, with `sum lambda = 0`, and not an integer multiple of a pair. Then:
1. `sigma_lambda >= q/3 = r/6`. The proof is elementary. At `q = 7`, the case `F = 0`, `n >= 6` uses
   `sigma >= n - 1 >= 5`.
2. For `q >= 11`, `sigma_lambda >= LB(q) = min(max(q/3, q/2 - 11 - 1.5 sqrt q), (q-3)/2)`. The term
   `q/2 - 11 - 1.5 sqrt q` uses Weil's bound for a squarefree monic quartic with the quadratic character
   ([L] Lidl-Niederreiter Thm 5.41, located at citation level by pass 2).
3. Hence `sigma_lambda >= 0.23515 r` at `q = 3251` (`LB = 1528.97`), and `>= 0.2349 r` at every regime `q`.
4. Sharpness: the dual pairs `f = (chi(x) - chi(x-1))/2` (with `f(0) = 0`, `f(1) = 1`) have `sigma = (q-1)/2`
   exactly, so `c_0 <= (q-1)/(4q)`. Some base 4-terms `e_x + e_y - e_z - e_inf` attain `(q-1)/2`, and the others
   have `(q+3)/2`.

The proof has five steps:
* the pattern-sum formula `sigma = sum_a max(|F + f(a)|, |(f*chi)(a)|) - q`;
* a reduction to ternary `lambda`;
* the defect identity from the 2-design property;
* a **new zero bound**: `f*chi` has at most `(q-1)/2 - j_0` zeros. This uses Euler's criterion
  `chi(u) = u^((q-1)/2)` and a Vandermonde argument;
* Weil, used for 4-terms only.

The zero bound closes the "intermediate norm range" that R16 §7.3 expected to be the real content.

**Exact data** (corroboration, not used in the proof):
* the global minimum of `sigma` is 3 at `q = 7` and 5 at `q = 11` (exhaustive);
* for `19 <= q <= 1091` and for `q = 3251, 3259, 3271, 5479`, it equals the exact 4-term minimum. This was computed
  by two independent C enumerations (author, pass 2), which together with the theorem's case split determine it.
  The values are 1571, 1575, 1581, 2667 at the four regime primes, so **`c_0(3251) = 1571/6502 = 0.2416` exactly**;
* at those four primes, the Lemma 2.7 bound is verified for every 4-term without Weil;
* the HiGHS MILP values for `q <= 43` are solver-optimal corroboration, not exact certificates.

**Corrections.**
* M1: the displayed `LB(7) = 2 < 7/3`. State `LB` for `q >= 11`, and `sigma >= 3` at `q = 7`.
* M2: not every base 4-term has `sigma = (q-1)/2`; some have `(q+3)/2`.
* M3: the MILP results are labelled "exact"; they are not certificates.
* M4: the Weil citation is [L], with its hypotheses checked.
* M5: §4 item 2 omits barrier Thm 2.3's thresholds `sigma*(tau) <= r+1`, `||lambda||_1 Delta <= 1/2` and `n >= 11`.
  These are all `O(log q)` or automatic. Also, "inf_theta" should read "some point of `C_1`"; the two agree because
  `kappa_H` is continuous on a compact circle.
* M6: the "proved LB/r" table rounds `sqrt q` up. This is conservative.

**Consequences (proved).**
* Barrier Lemma 3.4 holds unconditionally, with
  `tau_D = max((2/c_0)(2 ln(8q+5) + ln(1000/C)), 4 ln(8q+5) + 8, 2 ln(1000 q(q+1)^2/C))`.
* Barrier Thm 3.5 (with the R16 corrections `m' < q` and `AGL(1,q)`) is an **unconditional** cap test.
* Barrier Prop 2.5 (collision-free rational strata) needs transversal order `o >= 3` at regime `q`, and `o >= 4`
  for all `q`.
* With the proved `c_0`, the regime numbers change (J1, §6). At `q = 3251`:
  * `tau_D = 231.68` (Weil `c_0`), or `225.48` (exact `c_0`);
  * `ln T* = -376,602`;
  * the survival budget is `ln|c|_1 < 110.3 q`.

  The R16 evidence value `c_0 = 0.184` (used by cap, survive and calib) gave `296.09`, `-481,302` and `142.5 q`.
### 2.2 cap (`cap.md`, `referee_cap.md`): survives as "undecided"; one major logical inversion corrected

The referee reproduced every number and re-derived every proof, using independent code:
* the Gauss-sum identity, exact in `Z[zeta]` for every `q = 3 mod 4` up to 139;
* the exact `N_i` on `sdc7`, `L = 2`, `m = 6`;
* the margin MILP for `q <= 43`.

**Corrected statements.**
1. **Scale (proved arithmetic).** At `q = 3251`: `ln T* <= -101,333` for every `c_0` (floor `tau = 62.337`). With
   the proved `c_0`: `ln T* = -376,602` (J1); with `c_0 = 0.184`: `-481,302`.
   * A cap needs `kappa_H(hat u_theta) <= T*` for some `theta`.
   * A cap **built from an approximant on `A`** (a "`kappa_V`-cap") needs `kappa_V(hat u_theta) <= T*`.
   * **Major correction (M1).** cap.md says throughout (summary, §1, Cor 3.3, §8.4) that a cap *needs*
     `kappa_V <= T*`. That is false: since `kappa_H <= kappa_V`, the condition is only sufficient. A cap can come from
     `TC_{H'} \ TC_V`, and cap.md never examines that route. By survive Thm 3.4 (§2.3), such a cap means either
     `C_1 ⊂ TC_{H'}`, or that every integral witness of height `<= H'` has `ln|c|_1 >= -ln T* - m' ln 2q` (§3).
2. **Gauss-sum identity (proved).** `(1 - g) hat k(1) = 2 sum_a zeta^a S_a`, with `S_a = k({a} ∪ (a + QR))` and
   `g^2 = -q`. Hence `l_(u_theta) = sum_a c_a S_a` with `|c_a| <= 2/sqrt(q+1)`.
3. **Single relation families fail for approximants (proved).** Any degree-`< m'` approximant within `T* d^(m')` on
   `A` has error at least:
   * `e^(-4507) d^(m')` somewhere in `{-1,0,1}^q \ A` (ternary relations only);
   * `e^(-11312) d^(m')` somewhere on the `S`-grid outside `A` (one slab family only).

   Rewriting the slab powers touches at most `e^(-2770)` of the `l^1` coefficient mass of `(sum c_a S_a)^(m')`. The
   reading "no traction" from this is an **interpretation** (m2).
4. **Fat-point form (proved).** `l_psi^(m') in (J_gen)_(m')` iff `[psi]` is a common zero of the degree-`m'` forms
   with multiplicity `>= m' - 2` at the `q` coordinate points and `>= m' - 2j` at the `2q` pattern points.
   * `(J_gen)_(m')` has relative dimension `<= e^(-367)` modulo `(k^3)`.
   * `H_A(e) > 0` for all `e <= q`.
   * The dimension count says nothing about the special point `[psi]` (m3).
5. **`L^2` splitting (proved; also calib Prop 2.1).** For `m' < q`:
   * `Psi_A(u_theta) = 4^(-m') sum_i binom(m',i)^2 N_i`, constant on every circle;
   * `N_i = N_(m'-i)`, and each `N_i` is **rational** (Galois plus `AGL(1,q)`).

   So survival at the `kappa_V` level follows from one bound, `ln N_0 > -166,546` (at the floor; about `-7.2*10^5`
   needed in `ln Psi` at the proved `c_0`). A `kappa_V`-cap needs every `N_i` small. **Neither statement decides
   `kappa_H`.**
6. **Zonotope bound (proved from the refereed `sdc` margin).** Every zonotope certificate in `A` has order
   `<= floor(4jq/(q+1)) < 4j <= m'/2`.
7. **Toys (evidence).**
   * On `sdc11`, `L = 2`, the circles have the largest restriction ratio among the tested directions.
   * On `sdc7` they are typical: beaten by controls at `m = 5, 6`, and at `m = 8` by a uniform control at `5.4e-5`.
   * **Withdrawn:** "for `m < q` the circle has the largest ratio of all tested directions" (m5).
   * Also withdrawn: "height is not the obstacle (`ln H' ~ 10^8`)". This was a heuristic; on the circle the true
     statement is survive Thm 3.1.
   * Remark 6.2's "proved negatives" rest partly on an Alon-Furedi bound quoted from memory [M] (m4).
### 2.3 survive (`survive.md`, `referee_survive.md`): survives; one false number in a remark, conclusion restored

The referee:
* re-derived Lemma 2.1, (2.2), Thm 3.1, Thm 3.4, Cor 3.2/3.3 and Props 4.1/4.2;
* checked that the write-up's patterns equal the record's `sdc_q` for every `q <= 383`, and are `AGL`-invariant for
  `q <= 83`;
* checked Thm 3.1, including its item-2 uniformity mechanism on a 61-point `theta` grid, with exact `Z[zeta]`
  witnesses at `(q, L) = (7,1), (7,2), (11,1)`;
* reproduced the `sdc7` datum: `L = 2`, `m = 8`, `|c|_1 = 1008`, `<c, M_4> = 28224`, value `7.65625`.

Setting: `sdc_q`, `q = 3 mod 4` prime, `A = A_q(L)`, `1 <= m < q`, `M_i = hat k(1)^i hat k(-1)^(m-i)`.

1. **Thm 3.1 (equivariant witness bound; proved).** Let `c in V(A,m)` be integral with `alpha = <c, M_i> != 0`.
   Then for every `theta` and every `H >= |c|_inf`:

   ```text
   kappa_H(hat u_theta) >= 2^(-m) binom(m,i) max_s |sigma_s alpha| / (|c|_1 q^m) >= 2^(-m) binom(m,i)/(|c|_1 q^m).
   ```

   * *Mechanism.* The multipliers `mu_s` are box automorphisms that realise `Gal(Q(zeta_q)/Q)`. So one evaluates at
     the largest conjugate, which has modulus `>= 1` by the norm, instead of dividing by the other conjugates.
     Translation averaging, with exponents `2i - m` distinct mod `q`, gives uniformity in `theta`.
   * *What it removes.* This removes the R16 Liouville obstruction `e^(-Theta(q^2))` **for a given witness**. It does
     not bound the height of the minimal witness.
2. **Thm 3.4 (gap principle; proved).** For every `H >= 1`, either `C_1 ⊂ TC_H(A,m)` (and then every `C_b` is too),
   or `inf_theta kappa_H(hat u_theta) >= 2^(-m)/(|A| H q^m)`.
3. **Cor 3.2 / 3.3 (proved).** Put `L = j`, `m = m'`, `W = 2q tau`. Let `lambda_q` be the minimal `|c|_1` of an
   integral witness (`Phi_c` not identically 0 on `C_1`). If
   `ln lambda_q < q tau/2 - m' ln 2q - ln(Cq/(m'+1))`, then `inf_theta kappa_(H')(hat u_theta) > eps_0 > T*` on every
   circle, so the circle test does not cap route 5 at that `tau`, or at any larger `tau`.
   * The budget at `q = 3251` is `110.3 q` at the proved `c_0` (J1), `142.5 q` at `c_0 = 0.184`, and `>= 25.6 q` for
     every `c_0` at the `tau` floor.
   * A cap with `C_1` not in `TC_V` is possible only in the **height window** `ln H in [~(1/2) q tau, (m'/20) q tau]`.
     At the proved `c_0` and `q = 3251` this is `[3.5*10^5, 7.7*10^7]` (J1).
   * Scope: a pass of the circle test does not mean route 5 survives, because other barrier points of `G` are not
     excluded.
4. **Props 4.1, 4.2 (proved).**
   * `V(A_q(j), m')` contains no product-type measure: no zonotope, no sub-cube tensor, and no convolution of these.
     The order of such measures is `<= 8jq/(q+1) < m'`, by an exact Paley-Parseval identity and an extreme-point
     argument.
   * It contains no measure supported on a union of coordinate sub-cubes inside `A`; their order is `<= 4j`.
   * The witness must therefore be a non-Prouhet PTE-type configuration of order exactly `m'`.
5. **Remark 4.3 (proved numbers; scope narrowed).** One Bombieri-Vaaler pigeonhole puts a vector outside
   `V' = {c : <c, M_i> = 0 for all i}` only at height `exp(e^(1073))`, because `ln r_0 >= 1073.3`. "Counting cannot
   produce the witness" refers to this argument only. Fukshansky (arXiv math/0409374) was read at abstract level.
6. **Remark 4.4 (corrected; major).** The author's `min_p G = +1823` is a grid artefact: with the author's atom
   `p/2`, `G -> -infinity` as `p -> 1`. With the correct minimal atom `min(p/2, 1 - p)`, an interval-certified bound
   gives `G >= 2848 > -ln 2` on `(0,1)`. So the `p`-biased `L^2`/hypercontractive chain cannot close for any `p`,
   and the conclusion stands.
7. **Minor corrections.**
   * The equality in (2.2) is unjustified; only `>=` holds (counterexample with the cube roots of unity). Nothing
     downstream uses the equality.
   * The calibration factors "7-33 of LP" measure the bound before `max_s |sigma_s alpha| >= 1` is applied, and they
     depend on the witness: the referee found shorter witnesses at `(7,2)`.
   * "In practice at about `4j`" rests on sampled data.
   * Cor 3.3(1) "every barrier `tau`" should read `tau >= tau_D` plus the thresholds of Thm 2.3.

**Judge's note.** This is the round's most useful structural result for future work. On a box whose symmetry group
realises the Galois group of the test direction's field, the quantitative test reduces to an exact membership
question plus **one height**. The question it leaves is sharp and explicit:

> **(W_q)** Is there an integral `c in V(A_q(j), m')` with `Phi_c` not identically 0 on `C_1` and
> `ln |c|_1 < 110 q` (at `q = 3251`)?

Nobody knows how to construct such a `c`, or how to rule one out.
### 2.4 calib (`calib.md`, `referee_calib.md`): survives as evidence plus a small reduction; four labelling corrections

The referee's checks used only the referee's own code:
* an exact recomputation of 24 `sdc` cases and 5 uncut-cube cases (symmetry-free, rank certified over `Q`, solves
  checked by multiplying back): every rational agrees, including a 66-digit value;
* a floating-point check with an independent isotypic engine: `sdc19` `m = 5, 6` and `sdc23` `m = 5` agree to
  10 digits;
* the LP values, the box counts (757, 9373, 264291, 818755) and the §7.1 regime table, all reproduced.

**Proved.**
1. **Prop 2.1.** Let `A ⊂ {-1,0,1}^(F_q)` be translation-invariant and `1 <= m < q`.
   * `Psi_A(u_theta) = 4^(-m) sum_i binom(m,i)^2 ||R_i||^2` is independent of `theta`, and the `R_i` are pairwise
     orthogonal translation eigenvectors.
   * If `A` is also invariant under `x -> tx`, every `||R_i||^2` is rational.
   * Hence `inf_theta kappa_V(hat u_theta) >= Psi_A(C_1)^(1/2)/d^m`. This bounds **`kappa_V` only**.
   * The hypothesis `m < q` is needed: at `sdc7`, `L = 2`, `m = q = 7`, the cross terms are exactly nonzero.
   * The same argument makes `Psi_(A,H)` constant on `C_1` for every `H`.
   * This is the same result as cap Prop 4.1/4.2, found independently, and is credited once.
2. **Lemma 7.1.** If `|A| >= 3^d/2` then `rho(A) >= d + 1`. In the regime, `m'/rho(A) <= 0.633`.
3. **Bookkeeping.**
   * At `c_0 = 0.184`, a `kappa_V`-cap at `q = 3251` forces `ln Psi_A(C_1) <= -9.293*10^5`.
   * At the proved `c_0` the bound is `<= -7.20*10^5` (J1). The qualitative reading is unchanged.

**Exact data (toy boxes only; none is a regime box).**
* *Rigorous* (exact rational, recomputed by the referee): `sdc7` `L = 1` (`m <= 3`) and `L = 2` (`m <= 6`); `sdc11`
  `L = 1` (`m <= 3`) and `L = 2` (`m <= 5`); `sdc19` and `sdc23` `L = 1` (`m <= 3`). For these rows `k_lo` is a
  rigorous lower bound for `inf kappa_V`.
* *Exact modulo the multimodular caveat* (pivots chosen mod `p`; if `rank_p < rank_Q`, `Psi_A` comes out too large,
  which would break a lower bound): the other `L = 1` primes `q <= 103`; `sdc11` `L = 2` `m = 6..8`; `sdc19` and
  `sdc23` `L = 2` `m <= 4`.
* *Floating point* (evidence): `sdc19` `L = 2` `m = 5..7`; `sdc23` `L = 2` `m = 5, 6`. The `m = 7` and `m = 6` rows
  rest on 2 random directions, not 6.
* The `L = 1` boxes are the pattern-independent pair box. `L = 2` was not reached for `q = 31` or `q = 43`; the
  referee's DFS also failed for `q = 31`.

**Evidence (toys).**
* Below saturation the circle stays far from `TC_V`.
* `rel_A/rel_cube` lies in `[0.29, 1.31]` for `m <= rho - 2`, and it **decreases once `m` exceeds the slab degree
  `2L+1`**.
* The LP `kappa_V` varies by at most 15% along the circle.
* `V_H = V` at `H <= 4` in the 7 small cases computed.
* `tau* <= 2.22`, against `tau_D = 164-221` in the toys.

**Extrapolation (evidence, model-dependent, used in no proof).**
* A `kappa_V`-cap at `q = 3251` needs an average relative-residual decay of `>= 466` per degree at `c_0 = 0.184`.
  At the proved `c_0` the figure is about `364` per degree (J1).
* The toys show at most 1.51 per degree on average and 7.0 in the steepest single step. The exact uncut cube at
  `(3251, 2057)` decays at `0.856` per degree; the referee computed this at 1500 and 2500 digits, and it is not the
  `ln 2` quoted in calib.
* The margin `e^(4.8*10^5)` holds only in the cube-rate model.

**Corrections.**
* (major) "Rigorous `k_lo`" is overstated for the multimodular rows.
* (major) "A robust no in kappa" must read **`kappa_V`**. CIRCLE is a `kappa_H` question, and the evidence says
  nothing about `V_(H')` at `ln H' = m'W/40`.
* (major) The extrapolation parameter `m/rho` cannot be separated from `m/(2L+1)` in the toys. Measured by
  `m/(2L+1)`, the regime (4.0) is far outside the toy range (`<= 1.6`).
* (minor) "kappa_H within a factor 4 of kappa_V" is unsupported; the data give about 9.
* (minor) The constant is 0.27, not 0.28.
* (minor) `|c|_inf = 42` had no saved output when calib was written. The survive referee has since reproduced it
  independently.
* (minor) "1.00-1.04" is false at `sdc23`, `m = 5`, where the value is 0.98.
* (minor) The random controls are uniform in the cube, not on `G_1`.
* **Judge's addition.** calib §7.3 says a proved "no" needs a lower bound on one `||R_i||` plus (P1)/(P2). That is
  true for the `L^2` route only. survive Thm 3.1 shows that a "no" through an **integral witness** needs no (P1)/(P2)
  on the circle, only the witness height `(W_q)`.

## 3. Status of CIRCLE(q), H_light, route 5

| statement | status | source |
| --- | --- | --- |
| **H_light(`c_0`)** for `sdc_q` | **proved**: `c_0 = 1/6` for all `q >= 7`; `c_0 >= 0.2349` for `q >= 3251` (Weil); `c_0(3251) = 0.2416` exactly; `c_0 <= (q-1)/(4q)` | hlight Thm 2.8, Prop 2.6; two referee passes |
| Barrier Thm 3.5 (circle cap test) | **unconditional**, with the R16 corrections (`m' < q`, `AGL(1,q)`) and the Thm 2.3 thresholds | hlight §4 |
| **CIRCLE(q)**, any regime `q` | **undecided**: no cap, and no lower bound at regime size | all four write-ups |
| `kappa_V`-cap (an approximant on `A` within `T* d^(m')`) | not found. Each single relation family provably fails, by `>= e^(90000)` | cap Prop 3.2 |
| Cap through `C_1 ⊂ TC_(H')` (exact) | open; unexamined by cap (referee M1) | survive Thm 3.4 |
| Cap through `C_1 ⊄ TC_(H')` with `kappa_(H') <= T*` | possible only if **every** integral witness with `|c|_inf <= H'` has `ln |c|_1 >= -ln T* - m' ln 2q` (`~ 3.6*10^5` at `q = 3251`) | J2 (from survive Thm 3.1) |
| Survival ("no") | follows from `(W_q)`: one integral witness with `ln |c|_1 <` budget (`110.3 q` at `q = 3251`). No candidate exists. Product-type and sub-cube witnesses are provably excluded | survive Cor 3.2, Props 4.1/4.2; J1 |
| `kappa_V` on the circle | **evidence** for "far from `TC_V`" (toys; model-dependent extrapolation). A `kappa_V`-cap needs a per-degree decay of `~ 364`; toys show `<= 1.51` on average and `7.0` at most | calib, referee_calib; J1 |
| `H`-gap (`kappa_(H') << kappa_V`) on the circle at regime height | **open**; no toy shadow (`V_H = V` at `H <= 4` in 7 small cases) | calib, referee_calib |
| Collision-free rational strata (barrier Prop 2.5) | hypothesis 3 now free; needs `o >= 3` at regime `q`. **Open**, no candidate known | hlight §4 |
| **Route 5** (MIN^pair / `(R_H)`) | **STOPPED** (R16 §7.4 rule): open, **not capped**, **not unblocked** | this README |

**Route 5 reading.** The decisive test was well posed and its last hypothesis is now proved. Both directions
reduce to one explicit question about integral moment-matching measures on one symmetric box:
* is `C_1 ⊂ TC_(H')`? If so, the test caps route 5;
* does a witness with `ln |c|_1 < 110 q` exist? If so, the test passes.

Neither has a method:
* every certificate type the record can analyse stops at order `<= 8jq/(q+1)`, about half of `m'`;
* counting is `e^(1073)` orders of magnitude too weak;
* toy computation cannot reach `m'/(2L+1) = 4`.

So route 5 now has a sharply stated and well-understood residual, and no tool for it. It **stops**. Two points are
recorded for anyone who reopens it:
1. **"Yes"** would cap route 5 in its MIN^pair / `(R_H)` form. It would leave only the actual-point variant, which
   needs a measure of algebraic independence for about `n` logarithms. No such measure exists in the literature
   (lit18 §1.2).
2. **"No"** would not prove route 5. It would show only that this one family of candidate barrier points fails. MIN^pair
   for actual boxes would still be open. Route 5 cannot give growth in either outcome without a new input.

## 4. What it means for the growth rate

**Unchanged.** `M <= (2/3 + eps) log R/loglog R` remains the best bound for actual circles. No result of this round
touches actual clusters. Every statement is about the fake boxes `sdc_q`, which is the barrier side.

**Conditional theorems on record after round 17.** Only the rows that changed since R16 §5.3 are listed.

| hypothesis | conclusion | implication | hypothesis status |
| --- | --- | --- | --- |
| `inf_(C_1) kappa_(H') <= T*` at one regime `q` (no other hypothesis) | route 5 (MIN^pair / `(R_H)`) capped | proved (barrier Thm 3.5 + hlight Thm 2.8) | **open** (CIRCLE undecided) |
| `C_1 ⊂ TC_(H')(A_q(j), m')` at one regime `q` | same | proved (survive Cor 3.3(1) + hlight) | **open** |
| `(W_q)`: integral witness with `ln |c|_1 < q tau_D/2 - m' ln 2q - ln(q/(m'+1))` | the circle test does not cap route 5 (at every `tau >= tau_D`) | proved (survive Cor 3.2) | **open**; product-type and sub-cube witnesses excluded |
| collision-free rational `L ⊂ TC_H`, transversal order `o >= 3`, regime `sdc_q` | route 5 capped | proved (barrier Prop 2.5 + hlight) | **open**; no candidate |

**Rubric (b), fake-violated input: unchanged.** The input that fakes violate is still rigidity at the actual point,
which is spent in Step 0. This round's results are all about fakes.

**What the record now says about growth** (whole record, R4-R17 and lit18):
* *Capped at `log R/loglog R`* by refereed barrier theorems:
  * every local input (pairs, gcds, cut, residues);
  * character, GAP-Fourier and energy arguments;
  * single-place auxiliary polynomials (`tau* < 2`);
  * Baker/Matveev;
  * single-form Lang-Waldschmidt;
  * robust analytic inputs (round 10 Thm B);
  * the round-8 Subspace route `QS` (lit18 Props A, B).
* *Open, with no method:*
  * route 5's box statement `(R_H)`, now **stopped**;
  * route 3's genuine-prime inputs H^gen + L1;
  * `QS^prod`, at least as strong as the target;
  * composite-supported analytic inputs (**unclassified**);
  * P3-all;
  * full Lang-Waldschmidt with `eps -> 0`.
* *Conditional on a standard conjecture:* only Vojta's conjecture for one fixed variety (the blow-up of `X_7`, or
  `Mbar_(0,8)` through `(L_8)`). It gives the **uniform** bound.

## 5. Scores

Rubric of `round10/README.md` §1, each sub-score out of 10:
* (a) correctness after refereeing;
* (b) fake-violated input;
* (c) movement of growth or uniformity.

The overall score is a judgement, and it weighs (b) and (c) above (a). Nothing this round moves growth, so no (c)
exceeds 1. Partial (c) credit is only for proved sharpenings of the route-5 / barrier map.

| write-up | (a) | (b) | (c) | **overall** | reason |
| --- | --- | --- | --- | --- | --- |
| **hlight** | 9 | 1 | 1 | **5** | It completes its assigned task outright: a sharp theorem (up to `O(q^(-1/2))`) with all quantifiers, upheld by two independent passes. The new Euler-criterion/Vandermonde zero bound replaces the character-sum work R16 expected. It makes barrier Thm 3.5 unconditional. Deductions: the `q = 7` slip, an overstated "all base 4-terms", MILP values labelled "exact", and incomplete threshold bookkeeping. All are minor. It moves nothing on growth, and it is about fakes only. |
| **survive** | 7 | 2 | 1 | **5** | The round's best structural result. Thm 3.1 is Liouville-free: Galois is realised by box automorphisms. Thm 3.4 is a gap principle. Together they turn the quantitative test into exact membership plus one height, and they make the `kappa_V -> kappa_H` passage free on the circle. Props 4.1/4.2 exclude every certificate type the record knows. (b) = 2 because it pins down exactly what a survival certificate must be. Deductions: a false number in Remark 4.4 (major, though not load-bearing; the referee repaired the conclusion), an unjustified equality in (2.2), witness-dependent calibration factors, and wording overclaims. |
| **cap** | 5 | 1 | 1 | **3** | Correct and useful pieces: the Gauss-sum identity, single-family lower bounds with exact numbers, the fat-point form, the `L^2` splitting with rationality, and the zonotope order bound. Deductions: the major necessary/sufficient inversion (`kappa_V` vs `kappa_H`) runs through its headline, so the `TC_(H') \ TC_V` cap route is never examined; a false toy headline; interpretations presented as proofs (m2, m4). The `L^2` splitting duplicates calib's. |
| **calib** | 6 | 1 | 1 | **3** | Careful evidence: 68 rows, and every rational the referee recomputed agrees. Prop 2.1 (shared with cap) and Lemma 7.1 are proved. Deductions: three major labelling or scope errors ("rigorous" for multimodular rows, "no in kappa" for `kappa_V`, an extrapolation parameter the toys cannot identify) and several summary-level overstatements. As planned, it is evidence and cannot decide CIRCLE. |

## 6. Judge's own checks (unrefereed)

Scripts are in `judge_checks/`, and each writes `*_out.txt`. **None of this has been refereed.**

* **R17-J1 (60-digit arithmetic, closed forms): regime numbers under the proved H_light.**
  * *Script.* `j1_regime.py` evaluates barrier Lemma 3.4's `tau_D`, `T*`, the survival budget of survive (3.2), the
    `kappa_V`-cap requirement of calib (2.1), and `ln H'`. All at `C = 1`, `d = q`, `n = q + 1`, `r = 2q`.
  * *Reproduction.* At `c_0 = 0.184` it reproduces the write-ups exactly: `tau_D = 296.09`, `ln T* = -481,301.5`,
    budget `142.49 q`, `ln Psi` needed `<= -929,334`, decay needed `466.1` per degree, `ln H' = 9.90*10^7`. At the
    floor `tau = 62.337` it gives `ln T* = -101,332.9` and budget `25.61 q`.
  * *Updated values at `q = 3251` under the proved constants:*

    | `c_0` | `tau_D` | `ln T*` | survival budget | `kappa_V`-cap needs `ln Psi_A(C_1) <=` | decay needed per degree | `ln H'` |
    | --- | --- | --- | --- | --- | --- | --- |
    | `1/6` (elementary) | 326.89 | -531,356 | 157.9 q | -1,029,444 | 514.8 | 1.09e8 |
    | `0.23515` (Weil) | 231.68 | -376,602 | 110.3 q | -719,936 | 364.3 | 7.75e7 |
    | `0.24162` (exact) | 225.48 | -366,528 | 107.2 q | -699,787 | 354.5 | 7.54e7 |

  * *Notes.*
    * The cap test may use **any** admissible `tau`. The cheapest cap scale is therefore at the largest proved
      `c_0`.
    * The exact `c_0(3251) = 0.2416` rests on the two independent exhaustive 4-term enumerations plus Thm 2.8's case
      split.
    * At `q = 5479` the Weil values are `tau_D = 236.82`, `ln T* = -648,781` and budget `113.8 q`.
    * The height window of survive Thm 3.4 at the Weil `c_0` is `ln H in [3.55*10^5, 7.75*10^7]`, using
      `|A| <= 3^q`.
  * *Reading.* Nothing qualitative changes. The proved `c_0` makes the cap threshold about `e^(10^5)` less severe
    than R16's evidence value did, and the toys are still `~ 10^2` per degree short of a `kappa_V`-cap.
* **R17-J2 (two-line consequence of survive Thm 3.1; proof here).**
  * *Claim.* Let `lambda'_q = min { |c|_1 : c in V(A_q(j), m') integral, |c|_inf <= H', Phi_c not identically 0 on C_1 }`.
    1. If `inf_theta kappa_(H')(hat u_theta) <= T*` with `C_1 ⊄ TC_(H')`, then
       `ln lambda'_q >= -ln T* - m' ln(2q)`.
    2. Conversely, `ln lambda_q < -ln T* - m' ln(2q) - ln 2 - (1/2) ln(m'+1) - pi^2/(4q)` implies "no". This is
       survive's (3.2), rewritten.
  * *Proof of 1.* Thm 3.1 with `H = H'` and `binom >= 1` gives
    `T* >= kappa_(H') >= 2^(-m')/(lambda'_q q^(m'))`.
  * *Reading.*
    * Every non-exact cap must avoid **all** integral witnesses of `l^1` norm below `e^(3.6*10^5)` (at
      `q = 3251`, Weil `c_0`).
    * The thresholds for "yes is impossible" and for "no holds" differ by `ln 2 + (1/2) ln(m'+1) ~ 4.5`.
    * This is **not** an equivalence. If all witnesses are long, that does not make `kappa_(H')` small, because
      `V_(H')` may still equal `V` through real combinations. The undecided set is exactly survive's height window.
* **No other new mathematical claims.** Everything else in this README is a referee-corrected statement or a reading
  of the record.

## 7. Single most promising next step (whole record), now that route 5 stops

### 7.1 What any growth or uniform argument must avoid (from refereed barrier theorems)

Read together, the record says that an argument giving `o(log R/loglog R)` or a uniform bound must satisfy all of
the following:
1. **No local input.** Inputs satisfied by fake circles are capped: pairs, gcds, cut, residues, and characters
   (rounds 4-7, 9).
2. **No averaged analytic statement.** Robust ensemble inputs are capped (round 10 Thm B). "Almost all"
   exceptional sets swallow any `O(M)`-point selection (Prop C). Pointwise statements about Gaussian primes
   reconstitute the problem (Thm D).
3. **No single-place auxiliary polynomial on `X_M`** (`tau* < 2`, threshold round), and **no auxiliary polynomial
   on `Mbar_(0,k)`** for any `k` (threshold `L8.md` Thm A: the balanced class lies in the interior of the movable
   cone).
4. **No single linear form in logarithms** (Baker/Matveev, single-form Lang-Waldschmidt).
5. **No `S`-integer argument whose threshold depends on the prime set only through `|S|`** (lit18 Props A, B). It
   must use equality of norms, or a threshold tied to `prod_(p in S) p`.
6. **Route 5's determinant method at `T = 1`** reduces, after Step 0, to a box statement that fakes may test. Its
   decisive test was posed this round and is undecided, with no method in either direction (§3).

### 7.2 The candidates

| candidate | would give | fake-tested? | predicted by a standard conjecture? | method known? | judge's estimate of success in a further round |
| --- | --- | --- | --- | --- | --- |
| Route 5 reopened via `(W_q)` or `C_1 ⊂ TC_(H')` | a cap or a pass of one test only; **no growth either way** (§3) | it *is* the fake test | no | no: certificate order stalls at `~ m'/2`; counting is `e^(1073)` orders too weak | `~0` for growth |
| Route 3, H^gen + L1 | `W >~ M log^2 M/loglog M`, i.e. growth (round 14 transfer theorem with rigidity, after the associate-twist correction; round 15: rigidity alone moves nothing) | no (genuine primes) | heuristically (Poisson models) | no. Worst-case statements over adversarial `O(M)`-prime selections are needed (Prop C), and those reconstitute the problem (Thm D). Ford/Gorodetsky have the right shape but are short by an unbounded power (lit18) | `< 0.1%` |
| `(L_7)` (Ptolemy, `Mbar_(0,7)`) | the uniform bound | no | **no**: `K.beta = -7`, so Vojta is vacuous there; heuristics predict `(L_7)` **false** (expected dimension of balanced rational curves `> 0`) | refutation only (round 7 Q1-Q3, finite computations) | refutation plausible; a proof `~0` |
| **`(L_8)`** (Ptolemy, `Mbar_(0,8)`) | **the uniform bound** (ptolemy Steps 0-5, refereed) | **no**: a statement about integer `2 x 8` matrices and divisibility of minors, with no circle | **yes**: Vojta (`D = 0`) on `Mbar_(0,8)`, `K.beta = 10 > 0` | no: items 3 and 5 exclude every unconditional technique in the record | **`< 1%`** for a proof; useful evidence likely |
| `QS^prod` (lit18 §2.3) | growth or uniform | no | no (non-standard) | no; at least as strong as the target | `~0`; a refutation attempt is cheap |
| composite-supported analytic inputs | unknown | unclassified | — | none found (lit18) | `~0` for growth; a classification would be map-only |
| low-height window `P_+(N) < M^(2-delta)` (no fakes known) | restricted theorems only. Adversarial `N` can put primes above `M^2`, where fakes exist | no fakes known there | — | none | not a growth route |

By items 1-6, a "genuinely new mechanism" must be global, multi-place and tied to the actual arithmetic of `N`.
The only known statements of that type are Vojta-type inequalities, and in the record their sharpest arithmetic
form is `(L_8)`. **No candidate from this round's ideas** qualifies:
* the Galois-automorphism trick of survive Thm 3.1 needs a box symmetry realising the Galois group of the test
  direction, while actual angle vectors are transcendental;
* the Euler-criterion zero bound of hlight is specific to Paley patterns.

### 7.3 Recommendation

**Target: `(L_8)` as a self-contained Diophantine problem. The Vojta-conditional route stays the only route from a
standard conjecture.**

The statement (ptolemy §7): there are `c > 0`, `eta_0 > 0` and a proper Zariski-closed `Z ⊂ Mbar_(0,8)` such that the
following holds for anchored primitive integer vectors `P_0 = (1,0), P_1, ..., P_7`. Suppose
`det(P_i, P_j) = t_ij prod_(T contains i,j) n_T`, where the `n_T` lie in two nested chains at each prime and satisfy
`|log n_T - w| <= eta_0 w`. Then `sum log|Y_i| + sum log|t_ij| >= c w - O(1)` off `Z`.

**First tasks** (cheap; evidence only; neither decides `(L_8)`):
1. **Census of small-height configurations at `k = 8`.**
   * Enumerate balanced anchored 8-configurations with nested-chain minors exactly, at increasing `w`.
   * Record the minimum of `H_Sigma = sum log|Y_i| + sum log|t_ij|`, and the locus where it stays bounded.
   * The aim is to locate the exceptional set `Z` that any proof must remove. Candidates are pull-backs of the
     `(L_6)` genus-one families (`L8.md` Thm B), and forgetful and Kapranov curves.
   * Possible outcomes:
     * bounded-`H_Sigma` points Zariski-dense in `Mbar_(0,8)`: then either `(L_8)` is mis-stated or Vojta fails for
       `Mbar_(0,8)`. This is implausible but decisive for the record;
     * they concentrate on explicit subvarieties: this gives the first concrete `Z`;
     * nothing small is found: weak evidence for `(L_8)`.
2. **Lifting test.** Do the `(L_6)`-refuting families `E_t ⊂ Mbar_(0,6)` lift to balanced 8-point configurations
   with bounded residues? An exact rank count, like the overdetermination count in `L8.md` §5, would settle the
   cheapest version.
3. **In parallel, and independent of the above:** consolidate the proved results into one document for outside
   experts in Diophantine approximation. These are the `(2/3 + eps)` bound, the two Vojta-conditional theorems, the
   restricted Hadamard/Walsh theorems and the barrier map. The record now shows that further internal rounds on
   reformulations have very low yield.

**Honest assessment.**
* A proof of `(L_8)` would be an unconditional case of Vojta's conjecture with `D = 0` on a 5-dimensional rational
  variety, for a structured but unbounded family of points. Every unconditional method in the record is excluded by
  a refereed barrier theorem (items 3 and 5). Outside the record I know of no such method [M]. Chance in a further
  round: **well below 1%**.
* The census and the lifting test are likely to produce usable evidence: an explicit candidate `Z`, or a
  confirmation of the formulation. They will **not** move the growth rate.
* No other candidate in §7.2 does better. Route 3 and `QS^prod` are at least as hard as the problem. `(L_7)` is
  predicted false. Route 5 cannot give growth in either outcome of its test.
* **The judge's overall view:** the program has reached the point where the expected return of further internal
  rounds on growth is low. The most valuable action is item 3.

**Not recommended:**
* reopening route 5 (any form of `(W_q)`, `TC_(H')`, or toy calibration at larger `q`);
* further reformulations of MIN^pair / `(R_H)`;
* H^gen + L1 through analytic number theory;
* proof attempts on `(L_7)`;
* further literature sweeps (lit18 §3 lists the only full texts worth reading).

## 8. Records to correct, and files

### 8.1 Records to correct

The repository is read-only for this round. Apply these corrections on merge.
* **Top-level `README.md`, growth-rate status.** Add the following for round 17:
  * H_light is proved for `sdc_q` (`c_0 = 1/6`; `0.2349` in the regime), so the circle cap test is unconditional;
  * CIRCLE(q) is undecided: on the circles the quantitative test reduces to exact `TC_(H')` membership versus one
    witness height (survive Thm 3.1/3.4);
  * **route 5 is stopped**, neither capped nor unblocked;
  * there is no growth movement.
* **`round16/README.md` §5.3.**
  * In the H_light row, "evidence only (`c_0 ~ 0.184`, `q <= 43`)" becomes **proved** (round 17 hlight).
  * In §7.1, "No needs a lower bound for `kappa_H` on `C_1`" should add: it **suffices** to have one integral witness
    with `ln |c|_1 < q tau/2 - m' ln 2q - ln(q/(m'+1))` (survive Cor 3.2), and no (P1)/(P2) is needed on the circle.
* **`round16/barrier.md`.**
  * Lemma 3.4 and Thm 3.5: drop the H_light hypothesis, and take `tau_D` with `c_0 = 1/6` (all `q`) or `0.2349`
    (regime).
  * §3.3: the heuristic `sigma = r(E|lambda(S)| - 1/2) + O(...)` is superseded by hlight Thm 2.8.
  * Prop 2.5: `o >= 3` at regime `q`.
* **The partial snapshot `round17/README.md` in the repository** is replaced by this file. Its line on cap should
  also mention survive Thm 3.4. A cap through `TC_(H') \ TC_V` means either `C_1 ⊂ TC_(H')`, or that every integral
  witness of height `<= H'` is long (J2).
* **Within round 17:**
  * cap.md "a cap needs `kappa_V <= T*`" should read "a `kappa_V`-cap needs ...";
  * calib.md §0 "robust no in kappa" should read "in `kappa_V`, model-dependent";
  * survive Remark 4.4: `+1823` should be `>= 2848` (corrected atom);
  * hlight `LB(7)`: see §2.1;
  * the regime numbers at `c_0 = 0.184` in cap, survive and calib are superseded by J1.

### 8.2 Files (all under `round17/`)

| file | content |
| --- | --- |
| `hlight.md`, `hlight_notes.md`; `referee_hlight_pass1.md`, `referee_hlight_pass2.md` (= `referee_hlight.md`) | H_light theorem; two independent referee passes |
| `hlight_checks/` | `check_formula.py`, `check_lemmas.py`, `exhaust_F1.py`, `dualpair.py`, `weil4.py`, `milp_tern.py`, `c0table.py`, `min4.c` (exhaustive 4-term minima, `run_min4.sh`, `summarize_min4.py`), `regime_spot.py`, each with `*_out.txt` |
| `referee_hlight_checks/` | pass 1: `check1-5.py`, `r4.c` (4-term enumeration), `rlib.py`; pass 2 (`pass2/`): `c1_definition.py` (pattern rebuild from `exact.md`, SNF of `L_T`), `c2_exhaust.py`, `c3_misc.py`, `c4_regime.py`, `r4b.c` |
| `cap.md`; `referee_cap.md` (= `referee_cap_pass1.md`), `referee_cap_notes.md` | cap attempt; report with the `kappa_V`/`kappa_H` correction |
| `cap_checks/`, `referee_cap_checks/` | `gauss_identity.py`, `regime_numbers.py`, `slab_tail.py`, `toy_calibration.py`, `sdclib.py`; `regime_ref.py`, `gauss_ref.py`, `toy_ref.py`, `exact_Ni.py`, `margin_milp.py`, `r2_*.py` |
| `survive.md`, `survive_notes.md`; `referee_survive.md` (= `referee_survive_pass1.md`) | equivariant witness bound, gap principle, witness obstructions; report |
| `survive_checks/`, `referee_survive_checks/` | `sym_check.py`, `sdc7_check.py`, `zono_cost.py`, `r0_bounds.py`, `biased_l2.py`, `regime_budget.py`, `calib.py` (toy witnesses), `slib.py`; referee's `pats.py`, `parseval.py`, `toy.py`, `regime.py`, `g44.py`, `g44_cert.py` (certified Remark 4.4 bound) |
| `calib.md`, `calib_notes.md`; `referee_calib.md` | calibration (68-row table) and report |
| `calib_checks/`, `referee_calib_checks/` | `exact_iso*.py`, `float_iso*.py`, `cheb_lp.py`, `kappa_h.py`, `enumA.c`, `regime_need.py`, `insert_table.py`; referee's `exact_full.py`, `iso_float.py`, `enum2.c`, `lp_check.py`, `hred_check.py`, `regime_check.py`, `cube_regime.py`, `cube_logprod.py` |
| `judge_checks/` | `j1_regime.py`, `j1_regime_out.txt` (R17-J1). J2 is proved in §6 |
| `judge_notes.md`, `notes.md` | working notes |
| `cached_results.json` | the structured results of the six agents that finished before the container restarts |

**Repository note.** Commit only files under 1 MiB. Generated box enumerations, arrays and compiled binaries
(`calib_checks/data`, `referee_calib_checks/data`, `enumA`, `enum2`, `min4`) are regenerated by the scripts next to
them.

## Repository note

Copied from the session scratch directory. Excluded: `*.npz` arrays and files over 1MB (10 files, calibration data).
Because of two container restarts, the first six agents' return values are kept in `cached_results.json`. The H_light
referee ran twice (`referee_hlight_pass1.md`, `referee_hlight_pass2.md`); `referee_cap.md` and `referee_survive.md` are the
first passes (also saved as `*_pass1.md`).
