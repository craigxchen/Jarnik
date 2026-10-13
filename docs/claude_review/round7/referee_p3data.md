# Referee report: `round7/p3data.md` (genuine residue demands and the discreteness question (Q))

Referee: adversarial, round 7. Scripts and logs: `round7/referee_p3data_checks/` (own code; nothing
of the author's code is imported; the author's files were read only to learn data formats and
conventions such as the cut-mask encoding and the column order of the Paley profile).

## 0. Verdict

**The draft survives as a data paper, with corrections.** It claims no growth result and no barrier
for actual circles, and that is correct. Every exact computation I re-derived or recomputed
independently is right:
* the genuine collision moduli;
* the certified `W_cut` values;
* the equal-weight bound `LB_eq`;
* the actual-prime profiles, both the pair checks and the separation (placement) checks;
* the character minima.

The proofs of Lemma 1.1, Prop. 3.1, Prop. 3.2, Prop. 4.1, Remark 4.2, Example 6.2A (as a prime-floor
statement), Thm 5.1, Cor. 5.2, Lemma 6.1, Prop. 6.2, Lemma 6.3 and Cor. 6.4 are correct as stated.

Three interpretive claims do not survive, and two of them are in the Verdict.

1. **(major) "The fitted `M = 8` profiles are infeasible for their own genuine residues"** (Verdict
   item 11, §5.1 last paragraph, §6.3 item 2). This is an artefact of fixing the windings at `k = 0`.
   In P3-gen the windings are free: they come from the free angles `phi` through (1.1) of
   `outside.md`; see `sharp.md` Def. 1.1(b) and `p3lower.md` Def. 1.1.
   * At `k = 0` I reproduce the author's result on the author's 16 profiles: 0 of 400 feasible
     conjugation choices each.
   * Once a winding vector `k in (Z/4)^8` may be chosen, 23 of the 24 realised `M = 8` profiles become
     pair-feasible for their own genuine residues for some conjugation choices. The profiles with
     separation do so in a median of 23% of 400 choices, and in up to 90% (§2, issue 1).

   The same remark limits how Thm 5.1 may be read: it is a theorem for windings fixed independently of
   the residues.
2. **(major) A lower bound is read as an upper bound.** Prop. 4.1 and Cor. 6.4 are *lower* bounds for
   equal-weight profiles. They cannot show that "equal-weight granularity costs only a bounded
   factor" (author's summary; §6.3 item 1 "costs a flat factor"; the §6.2 paragraph after Cor. 6.4
   "discreteness costs an additive `~ M log X` ... a constant factor and not growth"). What is shown is
   only this: *the Plotkin-count obstruction* is a flat factor. The strictly equal-weight upper bounds
   found (`min r tau` of the structured SA, `best_ub_M*.jsonl`) rise, as medians, from 1.59 to 1.79 to
   1.96 times the Perron bound at `M = 12, 20, 24`.
3. **(medium) Example 6.2A is a statement about the relaxation `W_sc` only.** Every profile it
   constructs contains a Remark 4.2 tie. Each column is balanced inside each class, so for every pair of
   classes the sum of the `(k+1)^2` cross-pair margins vanishes identically in `w`. By the draft's own
   Remark 4.2, the profile is therefore never realisable by distinct primes (verified). Its use in the
   Verdict ("met at the prime floor", "a lower bound for (Q) cannot hold for all `D`") must say
   "`W_sc`". A realisable variant at cost `2 lambda_max + O(M log M)` is sketched in issue 3.

There are also eleven minor issues (§2, items 4-14). They include a miscount, a cross-class
comparison, wording overstatements and missing attributions.

## 1. Independent checks (all own code)

| claim | check | result |
|---|---|---|
| genuine unit-1 and any-unit moduli `m_xy` (§1, §2; inputs of all of §4) | `ref_demands.py`: direct Gaussian arithmetic `mod q^(A_q)` (no discrete logs), own Paley matrix and `construct` profile (construct.md §2), stored primes `pi_j = A_j + i B_j` | **all 120** stored `construct` instances (`M = 8,12,20,24,32`, `t = 0..11`, `X = 2M, 4M`) agree exactly for unit-1 **and** any-unit moduli; `D = 2 log m` to `1e-9`; every `p_j` prime, `= 1 mod 4`, `A^2+B^2 = p` |
| §2.1 table | `ref_sec2.py`: re-aggregation of `genuine_trials.jsonl`; own exact `E D = 2 sum_(q^a <= X) log q/((q -+ 1) q^(a-1))` | every entry reproduced (to the printed digits) |
| `2 lambda_max`, `W_cut` and the certified dual bound (§4, §6.2, Tables 4.1, 4.2a) | `ref_lp.py`: own column generation (HiGHS) with exhaustive max-cut pricing, certified bound `<nu,D>/surplus(nu)`; all 32 instances `M = 8, 12, 20, 24` | all 32 agree with `table42_rows.json` (`W_cut` to `1.2e-12` relative, `2 lambda_max` exactly); my certified gaps `<= 4.3e-10` relative |
| Prop. 4.1, `LB_eq` medians 1.12-1.22 | `ref_eqlb.py`: the **exact** infimum `min_r r max(w_0, tau_r)` (`tau_r` found by bisection over the breakpoints `D_xy/(n - r/2)`), all 120 instances | the author's grid value is `<=` mine in every instance and within 1.1% of it; medians 1.122-1.223, range [0.958, 1.278]; claim confirmed |
| 40 realised actual-prime profiles, `M = 12, 20, 24` (Table 4.2a, item 6) | `ref_profiles.py`: integer pair check `Pi_sep >= m^4 Pi_same`; primes distinct, prime, `= 1 mod 4`, `> X`; nonconstant columns, distinct rows; `W`; `stair(r)`; ternary support-4 and `2e_x - e_y - e_z` characters | all pass; `W` and `stair(r)` agree exactly; one prime per column throughout |
| character minimum over `|c_x| <= 12` (item 6) | `ref_charmilp.py`: own MILP (`M` subproblems `c_x0 >= 1`) | all 40 profiles: minimum 0.020-0.325, always attained at a pair, and equal to the least `mu` over `D = 0` pairs to `4e-8`; claim confirmed |
| `M = 8` actual-prime profiles with and without separation (Table 4.1) | `ref_m8sep.py`: exact pair check, and Lemma 1.1(2) placement with the stored rational positions | all 24 pass; `W` and `r` agree |
| `M = 8` `W_sc` optima ("7 of 12 proved optimal") | table: 6 stars; `ref_sc8.py`: own semi-continuous MILP over all 127 cuts | all 12 solved to optimality; the 6 starred values reproduced; the draft's value is also optimal at `X = 32`, `t = 3` (issue 4, §3) |
| Example 6.2A | `ref_ex62A.py`: construction for `(m, k+1) in {(2,2),(3,2),(5,2),(3,4),(5,4),(2,8),(6,8)}` | Gram matrix as claimed; **identical cross-class cancellation in every case**; LP `max_(w >= 1) min_cross mu = 0` |
| Remark 4.2 ties in stored solutions | `ref_ties.py`: LP `max_(w>=1) min_(D=0) mu` | every stored sign matrix at `M >= 12` is tie-free; ILS `M = 8`, `X = 16`, `t = 4, 5` have ties |
| §5.1 self-consistency | `ref_selfcons.py`: own genuine moduli for every unit `eta`, 400 random conjugation choices per profile, `k = 0` and all `k in (Z/4)^8` | issue 1 |
| Prop. 3.2 constants | direct | `(8/M) S_M = 22.77, 40.17, 76.25, 94.59, 131.6`; `4.515 M` at `M = 100`; `int = 0.59657` |

## 2. Issues

### Issue 1 (major): self-consistency is tested only at windings `k = 0`, but windings are free in P3-gen

*Where.* Verdict item 11 ("The fitted `M = 8` profiles are infeasible for their own genuine residues
(0 of 400 conjugation choices each)"), the last paragraph of §5.1 ("infeasible for their own
residues"; "adaptation makes things worse, not better"), and §6.3 item 2 ("At `M = 8` every fitted
profile is infeasible for its own genuine residues").

*Why the inference fails.* In P3-gen (`p3lower.md` Def. 1.1, inside `F^all_X` of `sharp.md` Def. 1.1)
the point residues are `rho_x(q) = i^(-k_x) prod_j r_j(q)^(a_xj) conj(r_j(q))^(1 - a_xj)`. Here the
windings `k` are determined by the angles through `<2a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x`.
The angles `phi` are free in P3-gen ("The angles `phi` are not required to be the arguments of the
`pi_j`"). For a profile whose rows `2a_x - e` are linearly independent, every `k in Z^M` is attained by some
`phi`. That holds for all 24 realised `M = 8` profiles, which have row rank 8 (checked).

The unit-1 collision of `zhat_x` and `zhat_y` is the collision `z_x = i^(k_x - k_y) z_y (mod Q)` of the
`k = 0` points. So the genuine demand of a member is
`D_xy(k) = D^(i^(k_x - k_y))_xy`, a different unit's demand for each pair. The draft (§1, §5.1,
`theorem_check.py` docstring: "random conjugates, `k = 0`") evaluates only `D(0)`.

*Test* (`ref_selfcons.py`, own moduli for all four units, 400 random conjugation choices per
profile). Results for the 24 realised `M = 8` profiles of `q_results_M8.jsonl`, of which the author's
16 are a subset:

| profile (`X`, `t`) | no-sep: feasible at `k = 0` | no-sep: feasible for some `k` | `+sep`: feasible at `k = 0` | `+sep`: feasible for some `k` |
|---|---|---|---|---|
| 16, 0 * | 0/400 | 8/400 | 0/400 | 86/400 |
| 16, 1 * | 0/400 | 3/400 | 0/400 | 35/400 |
| 16, 2 * | 0/400 | 8/400 | 0/400 | **230/400** |
| 16, 3 * | 0/400 | 8/400 | 0/400 | 70/400 |
| 16, 4 * | 0/400 | 12/400 | 0/400 | 36/400 |
| 16, 5 * | 0/400 | 10/400 | 0/400 | **359/400** |
| 32, 0 * | 0/400 | 6/400 | 0/400 | 97/400 |
| 32, 1 | 0/400 | 10/400 | 0/400 | **272/400** |
| 32, 2 | 0/400 | 12/400 | 0/400 | 65/400 |
| 32, 3 | **1/400** | **167/400** | **1/400** | **310/400** |
| 32, 4 | 0/400 | 2/400 | 0/400 | **219/400** |
| 32, 5 * | 0/400 | 0/400 | 0/400 | 1/400 |

`*` marks the eight instances (16 profiles) of the author's `theorem_check_out.txt`.
* On those 16 profiles at `k = 0` I reproduce the author's result exactly: 0/400 each. The mean
  numbers of violated pairs also agree (for example 9.1 against 9.0, and 6.4 against 6.6).
* With `k` free, 23 of the 24 profiles are pair-feasible for their own genuine residues for some
  conjugation choice. The separation profiles are feasible in a median of 23% of choices, and in up to
  90%.
* Even at `k = 0`, instance `X = 32`, `t = 3` (not in the author's set) has a feasible conjugation choice.
* All 24 sign matrices have row rank 8 (checked), so every `k` is realised by some `phi`.
* The check covers pair inequalities only, as the author's does. Characters, placement and `(e-all)`
  are not checked on either side.
* For an *actual* cluster, `k` is fixed by `phi_j = arg pi_j`. But then the points must really lie on
  one short arc, and nothing in §5.1 tests that either. The `k = 0` choice of the draft is not the
  actual-circle choice; it is simply one of the `4^(M-1)` choices.

*Consequences.*
* Verdict item 11 and the two sentences quoted above must be restricted to "with windings `k = 0`".
  With the windings free, which P3-gen permits, the fitted profiles are often pair-feasible for their
  own genuine residues. That is the opposite of the stated conclusion.
* The ranking "uniform beats Paley beats adapted" (Monte Carlo in `U_X`, §5.1) was obtained with `k`
  fixed. It is evidence only for that model.
* Theorem 5.1 and Cor. 5.2 are correct as stated: the model fixes `k` before the residues are drawn.
  But choosing `k` after seeing the residues is a selection of the kind the §6.3 remark discusses,
  with `4^(M-1)` options, and Thm 5.1 does not cover it.
  * The statement "the effect returns against independent residues" (Verdict bullet 3, Verdict
    item 10, §5 "What this means", §6.3 item 2) must say "with windings fixed independently of the
    residues".
  * *Referee's sketch (not checked in detail) that the order survives:* the proof of Thm 5.1 still
    works if the event `V_y` is replaced by `V'_y = {D^eta_(x_0 y) > mu_(x_0 y) + c_C for all four units eta}`.
    Feasibility for some `k` forces `#{y : V'_y} = 0`. Step 2 of the proof shows that the *full*
    collision patterns of `(x_0, y)` and `(x_0, z)` are independent, so the `V'_y` are pairwise
    independent. Then `G(t) := P(min_eta D^(eta, ind) > t)` replaces `Fbar_X(t)`.
  * In the window `(M^(1-eta), M]` each prime collides with at most one unit. Needing about `k'`
    collisions for each of the four units gives `log G ~ -4 k' log k'`, so the threshold becomes
    roughly `(1/2 - o(1)) M log^2 M/loglog M` instead of `(2 - o(1)) M log^2 M/loglog M`. If this is
    right, the order is unchanged and only the constant drops by a factor of about 4.

### Issue 2 (major): lower bounds for equal weights are read as statements of cost

*Where.*
* Author's summary: "so equal-weight granularity costs only a bounded factor".
* §6.3 item 1: "The granularity of equal weights costs a flat factor (exact, Prop. 4.1 ...)".
* §6.2, after Cor. 6.4: "So with equal weights, discreteness costs an additive `~ M log X`, ... i.e.
  a constant factor and not growth".

*Why it fails.*
* Prop. 4.1 gives `W >= LB_eq` for every equal-weight profile, and Cor. 6.4 gives
  `W >= 2 lambda_max + tau sum 4 u_x u_y phi_xy`. Both are **lower** bounds. A lower bound that is a
  flat multiple of the Perron bound says only that this obstruction (the Plotkin count) does not force
  growth. It does not say that equal-weight profiles of cost `O(lambda_max)` exist.
* The only equal-weight **upper** bounds in the draft come from the structured SA. Its strictly
  equal-weight values `r tau` (`per_tau`, before LP polish), minimised over the `tau` grid, have
  medians 1.59, 1.79 and 1.96 times `2 lambda_max` at `M = 12, 20, 24`. Against them,
  `LB_eq/2 lambda_max` is 1.16-1.20. After the polish, which makes the weights unequal, the values
  `W_sa` are 1.42, 1.64 and 1.83.
* The calibration says these are unreliable, so the question is open, which is what the
  hedged versions in the Verdict ("at the Plotkin level") and §6.3 item 3 say.

*Fix.* Replace the three quoted sentences with: "the Plotkin-count lower bound for equal weights is
a flat factor (1.12-1.22 times `2 lambda_max`, exact), so this obstruction alone cannot produce growth;
whether equal-weight profiles reach `O(lambda_max)` is open."

### Issue 3 (medium): Example 6.2A is valid for `W_sc` only; its profiles are never realisable by distinct primes

*Where.* Verdict item 5a ("met at the prime floor ... proved"); §6.1 "What this explains"; "It also
shows that a lower bound for (Q) cannot hold for all `D`"; and route (b) of the Verdict.

*The tie.* In the construction every column `(a, b)` takes the values `h_c(a) sigma_i(b)`. Within class
`c` these sum to `h_c(a) sum_i sigma_i(b) = 0` (the columns of `sigma` are balanced). Hence, for two
classes `c != c'` and every weight vector `w`,

```text
sum_(x in c, y in c') mu_xy(w) = -(1/2) sum_j w_j (sum_(x in c) s_xj)(sum_(y in c') s_yj) = 0.
```

All these pairs have `D = 0`. So Remark 4.2 applies with `c_xy = 1` on the cross pairs, and **no**
assignment of distinct primes satisfies the pair inequalities. For `k + 1 = 2` the rows `(c,0)` and
`(c,1)` are even complementary.

`ref_ex62A.py` confirms this for `(m, k+1) = (2,2), (3,2), (5,2), (3,4), (5,4), (2,8), (6,8)`. In each
case `max_(w >= 1) min_(cross pairs) mu = 0` exactly.

The example is still a correct theorem about the relaxation `W_sc`, which allows ties (§3). Since §3
defines (Q) through `W_prime` ("the discrete problem itself"), each use of the example must say
`W_sc`.

*A realisable variant (referee's sketch; it relies on primes `= 1 mod 4` in short intervals
`[x, x + x^theta]`, `theta > 7/12`, known since Huxley; not checked in detail).*
1. Add the Sylvester block of Prop. 6.2(3). It has `r' <= 2M - 3` columns, and every pair has margin
   `>= tau'/2` from it.
2. Take all main-block primes in `[P, P + P^0.6]` and all block primes in `[P', P' + P'^0.6]`, with
   `P, P' >= max(X, M^5, e^(2h/n_A))`.
3. The weight perturbations then change every margin by at most `(n_A k + r') P^(-0.4)/2 = O(M^(-1))`.
   So:
   * the cross-class margins become `>= tau'/2 - o(1) > 0`;
   * the in-class margins stay `>= h`.
4. This gives `W_prime(h(B - I)) <= 2 lambda_max + O(M log M)`.

The qualitative message of the example (a single layer is cheap) therefore survives for actual
primes, with `O(M log M)` in place of `2 M w_0`.

*Also.* "So at the floor the difficulty of genuine `D` is merging its `~ X/log X` layers into `O(M)`
shared columns" is an interpretation, not part of the proved example. Genuine layers have
Poisson-sized, unequal classes, which the example does not cover without padding. Verdict item 5a
labels the whole row "proved"; the second sentence should be labelled heuristic.

### Minor issues

4. **"7 of 12 proved optimal" at `M = 8`: only 6 are proved in the draft.** The claim appears in
   §4.1 observation 2, §4.2 reading 3, Verdict item 6 and §6.3 item 1. Table 4.1, `q_table_out.txt`
   and `table42_rows.json` (`milp*`) all mark exactly six optima. §3 and §4.2 (calibration) themselves
   speak of "the six proved MILP optima".
   * My own MILP solves all 12 instances to optimality (§3). It shows that a seventh value of the
     draft, `X = 32`, `t = 3`, is in fact optimal, so "7 of 12 values optimal" is true.
   * But it was not proved in the draft. The other five upper bounds are 0.6-6.4% above the optimum.
5. **The comparison with the Paley fakes is across classes** (Verdict item 6, §4.2 reading 5: "1.8-2.4
   times cheaper than the genuine Paley fakes"). The two sides are different kinds of object.
   * `W_req = r tau_req` is an *equal real-weight* scale. It cannot be realised by distinct primes:
     at `M = 12` it would need 100 distinct primes `= 1 mod 4` near `e^3.67 ~ 39`.
   * `W_prime` is a verified actual-prime profile, but it is fitted to a *foreign* `D`. The Paley fakes
     are pair-feasible for their own residues, at `k = 0`; the adapted profiles are not (issue 1).
   * The actual-prime floor of the `construct` profile at `X = 2M` is `stair(9M - 8)`: 615, 1165 and
     1454 at `M = 12, 20, 24`.

   So the ratio is neither like for like nor an upper bound on the advantage. Say "fixed-`D` cost
   comparison" and give the floor.
6. **Inaccurate range at `M = 8`.** §4.1 observation 3 says actual primes "add 30-75%". The verified
   values of `W_prime/W_sc(ub)` range over 1.11-1.76; the instance `X = 32`, `t = 5` adds only 11%
   (`ref_m8sep.py`).
7. **"Exact check of characters of support `<= 4`" is narrower than it sounds** (§3, §4, Verdict
   item 6). `char_small` checks ternary `c` of support 4 and `2e_x - e_y - e_z` only. Zero-sum `c` of
   support `<= 4` such as `3e_x - e_y - e_z - e_w` or `(2, 1, -1, -2)` are covered only by the
   floating-point MILP. My own MILP confirms the minima (§1), so nothing is wrong, but the wording
   should say "ternary support 4 and `(2,-1,-1)`".
8. **Verdict item 6 garbles its list.** "Verified actual-prime profiles: 1.7, 1.8 / 1.5, 2.3, 2.7 / 2.3"
   should read, at `X = 2M` / `X = 4M`: 1.74 / 1.73 (`M = 8`), 1.80 / 1.53 (12), 2.31 / 2.30 (20),
   2.70 / 2.25 (24) (Table 4.2b).
9. **"4 of the 10 instances have excess `< 0.01`" (§6.2) should be 5.** The excesses are 0.0027, 0.0058
   and 0.0049 at `M = 20`, and `1.3e-5` and `2.4e-5` at `M = 24` (`perron_check.json`; confirmed by
   `ref_lp.py`).
10. **Attribution.**
    * Theorem 5.1 has a direct precursor, `sharp.md` Prop. 6.5: independent fibres cannot reach
      `O(M log M)`. It uses the same Poisson tail at primes `q in (M^(1-eta_0), M]`, for single-prime
      monomials of one construction. Thm 5.1 is the extension to pair inequalities of *every* fixed
      profile. Prop. 6.5 should be cited.
    * The identity (6.1) of Lemma 6.3 is the exact residual of the Perron dual certificate
      `nu_xy = 4 u_x u_y` in the proof of `p3lower.md` Thm A(2), and should be credited as such.
11. **"At real weights the effect is absorbed by adaptivity (proved, `p3lower.md` Thm A)"** (§6.3
    item 1, Verdict). The proved part is `W_cut(D) <= pi lambda_max(D)`. That `lambda_max(D)` is of
    order `M log M` rather than `(M - 1) D_max` is a property of the genuine `D`, shown here by data
    (§2.1: `lambda_max/(M log M) = 1.5-2.0` for `M <= 84`). The sentence should say "proved for every
    `D` with `lambda_max(D) = O(M log M)`; the genuine `D` here have this property".
12. **Remark 4.2.**
    * It uses `c_C = 0`, i.e. `C = sqrt 2`. For `C > sqrt 2` the `D = 0` pairs need only
      `mu >= -c_C < 0`, so ties do no harm. The draft's general sentence "any other fixed `C` shifts
      every pair constraint by the same constant" does not carry this remark along.
    * "Equal-weight solutions often contain such ties" is not documented. The equal-weight SA matrices
      `S_sa` are not saved (`best_ub.py` stores only the post-ILS matrix). The only recorded ties are
      the ILS solutions at `M = 8`, `X = 16`, `t = 4, 5` (`ref_ties.py` and the author's `notes.md`).
13. **Prop. 6.2, "for actual primes ... a short-interval prime count, which is not proved here".** This
    is available unconditionally. The window is `[P, P(1 + O(1/M))]` with
    `P = e^w ~ e^2 (M - 1)^4`. Its length is `>> P^(3/4) > P^(7/12 + eps)`, and Huxley-type
    prime-number theorems in short intervals for primes `= 1 mod 4` give `>> P^(3/4)/log P >> M`
    primes there. So item 3 holds for actual primes as stated, without the detour through
    `construct.md`. This is an improvement, not an error.
14. **Calibration scope.**
    * The control targets of `test_unif.py` are *exactly* Plotkin-tight, with zero slack, which is the
      hardest case for local search.
    * The genuine targets have Plotkin slack: `LB_eq` is about `1.2 x 2 lambda_max`, and the best found
      is 1.6-1.8 times it.
    * Failure on the tight controls therefore shows that the solvers can miss structured optima. It
      does not measure how far the genuine upper bounds are from the optimum.
    * "Say nothing about the optimum" is too strong in the other direction: they are valid upper
      bounds. The conservative conclusion (no evidence of growth) is right.

## 3. Proofs re-derived, and claims confirmed

*Proofs* (each re-derived line by line; no gaps found beyond the issues above):
* **Lemma 1.1.**
  * (1) follows from `|sin(v/2)| <= |v|/2` and `2u/Delta = m e^(-mu/2)`.
  * (2): `arcsin u <= u (1 - u^2)^(-1/2)`, `e^(-2mu) = Pi_same/Pi_sep`, and
    `(1 - a/2)^2 >= 1 - a` with `a = 2u^2`. The hypothesis forces `m^2 < Pi_sep`.
* **Prop. 3.2.**
  * Gap sums: `sum_i (t_(i+k) - t_(i)) <= k`, then AM-GM.
  * Cut identity `sum mu <= M W/4`.
  * `S_M ~ (1/4 + (log 2)/2) M^2`.
  * The Fekete and equal-spacing comparisons (`(log 2) M^2`, `(3/4) M^2`) are correct.
* **Prop. 4.1.** `T_xy` is nonincreasing in `tau`, the Plotkin count is `c(M - c) <= floor(M/2) ceil(M/2)`,
  and the treatment of `tau >= a_n` through `T_inf` is right. No monotonicity in `r` is needed, since
  `r_P` is defined as a least element.
* **Remark 4.2** (at `c_C = 0`).
* **Example 6.2A** as a `W_sc` statement:
  * `(h_c . h_c')(sigma_i . sigma_i')` gives `-n_A` and `0`;
  * the columns are nonconstant;
  * `W = max(n_A k w_0, 2hk) <= 2 M w_0 + 2 lambda_max`.
* **Theorem 5.1.**
  * The pairwise independence argument is right: the `2 x 2` minor is `+-1`, so the map is onto
    `(Z/m)^2`.
  * Chebyshev step.
  * Markov step: at least `(M-1)/2` of the star margins are `<= W/(M-1) + 2c_C`.
* **Corollary 5.2.** The window `(M^(1-eta), M]`, `k <= (1 - eps/8) L/l + 1`, the Poisson-binomial lower
  bound of `p3lower.md` Lemma C.1, and `log P(c >= k) >= -(1 - eps/16) L` for large `M`.
* **Lemma 6.1.**
* **Prop. 6.2**, items 1-4:
  * row sums via `log m! >= m log m - m + 1`;
  * `n = 2^ceil(log2 M) <= 2M - 2`;
  * Sylvester margins `w(n - r')/2 >= w/2`;
  * the log form of Lemma 1.1(2), with `m^2/Pi_sep <= e^(-1)`.
* **Lemma 6.3 and identity (6.1)**, from `G = W I - 2 mu` off the diagonal.
* **Cor. 6.4.**
* **The remark on selection in §6.3**: rows of `S diag(w)^(1/2)` with pairwise nonpositive inner products
  give `M <= 2r`.

*Data confirmed by independent computation* (details in §1):
* every genuine modulus behind §§2 and 4 (120 instances, both unit-1 and any-unit);
* the §2.1 table;
* every `2 lambda_max` and `W_cut` with its certificate, at `M = 8, 12, 20` and `M = 24` (see below);
* every `LB_eq` value and the medians 1.12-1.22;
* every one of the 40 actual-prime profiles at `M >= 12` and the 24 at `M = 8`, with separation where
  claimed;
* every character minimum (all at pairs, 0.020-0.325);
* the ILS and SA percentages of §4.2 reading 3 (4.2-13.0% at `M = 20`, 17.9-24.5% at `M = 24`);
* `W_prime/stair(r)` (3.2-6.2% at `M = 20`, 2.2-4.1% for SA at `M = 24`, 73% for the ILS profile);
* the column ranges `1.5M-4.7M`.

*`W_cut` at `M = 20, 24`* (`ref_lp_2024.log`).
* All 12 instances agree with the draft to `1e-11` relative.
* My certified dual bounds close the gap to `<= 5e-12` relative.
* The excesses `W_cut - 2 lambda_max` at `M = 24` are `1.3e-5` and `2.4e-5` (`X = 48`), and 0.776 and
  `1.7e-5` (`X = 96`).
* So the "4 of the 10" of §6.2 is 5 (minor issue 9). Counting `X = 96` as well, it is 6 of the 12
  instances at `M = 20, 24`.

*`W_sc` at `M = 8`, proved optimal for all 12 instances* (`ref_sc8.py`: own semi-continuous MILP
over all 127 cuts, HiGHS, every instance solved to proven optimality).
* *Validity of the variable bound.* The bound `lambda_s <= U = 4 D_max + 4 w_0` loses nothing. For any
  feasible measure and any cut `s`, a pair on one side of `s` has `0 <= mu <= (W - 2 lambda_s)/2`, so
  `lambda_s <= W/2`, and `W/2 < U` at every optimum.
* *Agreement with the draft.* The six optima the draft marks `*` are reproduced to `1e-6`.

| `X`, `t` | draft `W_sc [lb, ub]` | referee optimum | optimum / `2 lambda_max` |
|---|---|---|---|
| 16, 0 | [51.0, 56.8] | 55.334 | 1.307 |
| 16, 2 | [48.0, 52.9] | 52.656 | 1.292 |
| 16, 3 | [57.6, 65.1] | 61.172 | 1.163 |
| 32, 1 | [75.6, 80.8] | 79.886 | 1.194 |
| 32, 3 | [94.9, 97.8] | 97.829 (= the draft's upper bound) | 1.110 |
| 32, 4 | [67.1, 73.4] | 70.780 | 1.236 |

* *Count of optimal values.* The draft's value is optimal in 7 of 12 instances: six were proved, and
  `X = 32`, `t = 3` is optimal but was not proved. So "7 of 12 proved optimal" is wrong as written, but
  "7 of 12 optimal" is true.
* *Corrected ratios.* With the true optima the `M = 8` medians of `W_sc/2 lambda_max` are 1.31
  (`X = 16`) and 1.21 (`X = 32`), against 1.32 and 1.22 in Table 4.2b. Nothing changes qualitatively.
* *Solver error at `M = 8`.* The draft's upper bounds exceeded the optimum by 0-6.4%.

## 4. Corrected statement

**`round7/p3data.md` (data).** There is no growth result and no new barrier for actual circles. The
discreteness question (Q) is undecided.

*Proved, with all quantifiers.*
* *Every member of `A_X(M, C)`, hence every actual cluster:* Lemma 1.1(1) and Prop. 3.2.
* *(Q) and its relaxations:*
  * Prop. 3.1, Lemma 6.1, Prop. 6.2, Lemma 6.3 with identity (6.1), Cor. 6.4.
  * Prop. 4.1, which is a **lower** bound for equal-weight profiles only.
  * Remark 4.2, at `C = sqrt 2`.
  * Example 6.2A, **for the relaxation `W_sc` only.** Its profiles always contain a Remark 4.2 tie,
    so they are not realisable by distinct primes as built. A tie-breaking block plus primes in short
    intervals should give `W_prime <= 2 lambda_max + O(M log M)` (referee's sketch).
* *The model `U_X`, with windings fixed independently of the residues:* Theorem 5.1 and Cor. 5.2.
  If `k` may be chosen after the residues, the star argument appears to survive with a constant about
  4 times smaller (referee's sketch, issue 1).

*Exact per instance; evidence about `M -> infinity`.* For genuine unit-1 demands of the fully flipped
Paley profiles with actual primes (`M <= 84` for the statistics, `M <= 24` for (Q)):
* The demands match the independent residue model.
* `W_cut/2 lambda_max = 1.00-1.05` for `M >= 12` (certified).
* The Plotkin obstruction for equal weights is a flat factor: `LB_eq/2 lambda_max` has median
  1.12-1.22 (`M <= 32`).
* Prime-floor upper bounds run 1.22-1.82 times `2 lambda_max` (medians). The calibration shows the
  solvers to be unreliable at `M >= 20`, so these say nothing about the asymptotics.
* 64 actual-prime profiles were verified exactly, with all character minima (`|c_x| <= 12`) attained at
  pairs.
* *Self-consistency at `M = 8`:*
  * at windings `k = 0`, the 16 fitted profiles of the author's check are infeasible for their own
    genuine residues;
  * with the winding freedom that P3-gen allows, 23 of 24 fitted profiles are pair-feasible for their
    own residues for some conjugation choice (separation profiles: median 23%, up to 90% of choices).
    So the data give **no** evidence that adapted profiles are self-inconsistent in P3-gen.

*Not shown, and to be removed:*
* "equal-weight granularity costs only a bounded factor";
* "the fitted profiles are infeasible for their own residues", unless qualified by `k = 0`;
* "7 of 12 proved optimal" (6 were proved; the referee's MILP proves all 12 and confirms that 7 of the draft's values are optimal, §3);
* "1.8-2.4 times cheaper than the Paley fakes", unless labelled as a cross-class, fixed-`D`
  comparison.

*Open, as the draft says:* whether adapted near-tight `+-1` designs with `O(M)` prime columns exist
for genuine `D`. Any growth proof through pair data needs an input about how actual residues relate to
the exponent profile; with free windings this now includes the choice of `k`.

## 5. Files (`round7/referee_p3data_checks/`)

| file | content |
|---|---|
| `ref_demands.py` | genuine moduli by direct Gaussian arithmetic; own Paley and `construct` profile; 120 instances against the stored data |
| `ref_sec2.py` | §2.1 aggregates; exact independent mean |
| `ref_lp.py` -> `ref_lp_out_8_12.json`, `ref_lp_2024.log` | own `W_cut` column generation with exhaustive pricing and certified bound |
| `ref_eqlb.py` | exact infimum of Prop. 4.1 for all 120 instances |
| `ref_profiles.py` -> `ref_profiles_out.json` | exact verification of the 40 realised profiles at `M >= 12` (pairs, primes, `stair`, ternary support-4 characters) |
| `ref_charmilp.py` -> `ref_charmilp_{12,20,24}.json/.log` | own MILP for `min s_c`, `|c_x| <= 12` |
| `ref_m8sep.py` | exact verification of the 24 `M = 8` profiles, including Lemma 1.1(2) placement |
| `ref_sc8.py` -> `ref_sc8.log`, `ref_sc8_out.json` | own semi-continuous MILP at `M = 8` |
| `ref_ex62A.py` | Example 6.2A: construction, and the tie that is always present |
| `ref_ties.py` | Remark 4.2 ties in stored sign matrices |
| `ref_selfcons.py` -> `ref_selfcons_400.log` | §5.1 self-consistency at `k = 0` and with free windings |
| `notes.md` | working log |
