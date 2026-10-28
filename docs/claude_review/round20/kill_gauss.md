# Kill report: gauss.md (round 20)

Target: `round20/gauss.md` (status claimed by its author: `reduced_to_cap`). Scripts and outputs:
`round20/kill_gauss_checks/` (own code; nothing imported from `gauss_checks/`). Labels: *proved*, *exact*,
*evidence*, *heuristic*, as in the overview. Nothing here is refereed.

## 0. Verdict

**KILLED as a route.** It gives no growth improvement and no uniform bound, unconditionally or from anything weaker
than Vojta. The needed statement `(L_5^circ)` / `A_5(δ<1/30)` is **capped**:

* Every method in the record is **exactly critical** on it. Lemma K1 (§6.1) shows that dividing a tuple by its
  Gaussian gcd is invisible to every forced-divisibility, auxiliary-form, pair or single-form argument: each such
  inequality for `z` (critical arc) is *literally the same inequality* for `ζ = z/g` (the "subcritical" arc).
  So the advertised move "strictly below the critical exponent 1/2" changes the normalisation of `R` and nothing
  else. For the cube profile on the `k` chosen rows, `construct.md` Cor. C3 still caps every auxiliary polynomial
  at ratio exactly `1`. Single-form Lang–Waldschmidt with `ε = 0` is also exactly tight (§6.2).
* The target `e` is a **rational** point of `T_5`. In approximation terms, `A_5(δ)` sits strictly between the
  Dirichlet exponent (`δ = 0`) and Liouville (`δ = 1/30` balanced, `δ = 1/10` general). That is the
  McKinnon–Roth regime for a rational point of a 4-fold, where nothing beyond Liouville is known
  unconditionally. The only thing that predicts it is a Vojta instance, `V_E(5)`, whose margin is
  `κ = 2δ ≈ 1/15`.
* Fakes violate `(L_5^circ)` (gauss §5 item 1, correct). So a proof must use the genuine Gaussian angles at every
  prime of `N''` at once, which is the "global, multi-place" input the record has never had.

Its position is the same as `(L_8)`'s after round 18, as the proposal itself concedes. Its proved content is
correct (with the minor corrections of §1.2). The model it works in was already set up in `threshold/L8.md` §6.1
(the class `γ̂`, and a "Gaussian analogue of Theorem D with 6 points … I did not do this").

**What survives** (all re-derived here; numbers recomputed with independent code; §7 gives the exact statements):

1. Theorem 3.4 and 3.4′: `A_k(δ)` ⇒ uniform. Prop 3.1: `A_k ⇒ (L_k^circ)`. Lemmas 1.2, 3.2, 3.3. Correct and
   elementary.
2. Theorem 4.1: `A_4(δ)` is false for every `δ < 1/14`, and `(L_4^circ)` is false. Exact Jacobian rank 3 was
   reproduced at 6 rational points. Corollary (new here): `A_3(δ)`, `δ < 1/6`, and `(L_3^circ)` are false.
3. Theorem 4.2 (conditional on `V_E(k,ε)`): a correct, non-vacuous Vojta instance at `k = 5`, `D = E`,
   `S = {∞}`. It improves `L8.md` §6.1's pointer from 6 points to 5.
4. Theorem 4.3(ii) is correct as stated. Theorem 4.3(i) needs a narrowing (§1.2 item 4).

## 1. What the proposal claims, and the issues found

### 1.1 Claims

* (P) Theorem 3.4: `A_k(δ)`, `δ ≤ 2^{-k-1}` ⇒ uniform. Theorem 3.4′: any `δ < 1/(2^k-2)`, and also via
  `(L_k^circ)` through `ptolemy` Steps 0–5 with Step 2 run in `T_k`.
* (R) Theorem 4.1: `k = 4` is refuted by a dense rational family at exponent `3/7`.
* (C) Theorem 4.2: `V_E(k,ε)` ⇒ `A_k(δ)` for `κ = (k-1)(1+2δ)/4 - 1 > 0`.
* (M) Theorem 4.3: the curve-family refutation mechanisms fail for circles at `k ≥ 5`.
* (S) Status: `(L_5^circ)` is open, Liouville is tight at `7/15`, and Roth/Schmidt do not apply.

### 1.2 Issues (severity; full list in the structured output)

1. **Major (interpretive): the "subcritical" gain is a change of units, not of difficulty.** Lemma K1 (§6.1)
   shows that every forced-divisibility, Liouville, pair or single-form inequality for the reduced tuple
   coincides with the one for the unreduced tuple. The only things that see the reduction are the heuristic count
   and Vojta's `D = E` term. The summary sentence "puts the reduced tuple … strictly below the critical exponent"
   must be read as "below `1/2` in the reduced radius, and still critical for every method in the record".
   The proposal half-says this in §5 ("the gain is exactly eaten by Liouville's tightness"). It does not say that
   the tightness is *forced by invariance*. That point is what makes this a cap and not an open lane.
2. **Minor: `t_ab` is not always an integer.** In the definition of `(L_k^circ)`, "Then `t_ab` is a nonzero
   integer" is false for tuples whose points lie in different Gaussian unit classes. In Theorem 4.1's own family,
   408 of 1866 tested pairs have `|c_ab|^2 ≡ 2 (mod 4)`, where `c_ab = (ζ_a-ζ_b)/gcd(ζ_a,ζ_b)`. So
   `|t_ab| = |c_ab|/2 ∈ √(odd/2)`, for example `1/√2` (exact, `k4_lattice2.py`). The cited integrality
   (`ptolemy` Thm B(i), item 532) holds inside **one unit class**, which ptolemy Step 0 enforces. Repair: add
   "in one unit class" (costing a factor 4 in `M` in Theorem 3.4′), or allow `t_ab ∈ 2^{-1/2}Z[√…]`. Either way
   `H_Σ` shifts by at most `C(k,2)·(log 2)/2`, and nothing downstream changes.
3. **Minor: one numerical range is overstated as general.** "Reduced diameter `/R''^{3/7}` between 8.3 and 43.4"
   (Thm 4.1 check (i)) depends on the sample. Over random `s ≤ 10^30` the asymptotic ratio is
   `≈ 840·Norm(g)^{-2/7}`, which ranges over `[8.31, 142.28]` as `Norm g` runs over its values `500 … 10368000`
   (exact, `k4_lattice2.out`). Still bounded, so harmless.
4. **Minor: Theorem 4.3(i) proves less than it states.** Its claim that the approached point `O` is rational or
   real quadratic is not proved for a **single** curve. The text excludes orbits of degree `≥ 3` "by the
   computation in (ii)", but (ii) is a nefness argument valid only for **covering** families. For one curve,
   Roth gives `μ_0 ≥ (1/4+δ/2)B.C/(2+ε)` and `E.C ≥ dμ_0`, with no contradiction. Correct (i) to: genus `0`,
   and `O` real algebraic with `μ_0 ≥ (1/4+δ/2)B.C/(1+ε)` (`O` rational) or `/(2+ε)` (otherwise). (ii) already
   treats all orbit degrees and is unaffected.
5. **Minor: §0(a) overclaims the density mechanism.** "[the] three mechanisms (density, one rational curve,
   covering genus-one family) provably do not transfer (Theorem 4.3)": Theorem 4.3 says nothing about the density
   mechanism (`L8.md` Prop C). For that, §2.3 itself offers only a **heuristic** count `e^{-w}`. Unlike round 18
   (theory Prop 2.3, a proved lattice-index argument at `k = 8`), there is no proof at `k = 5` for circles.
6. **Minor: the φ-test of §2.1 is not evidence.** It reads the families with the single choice `ρ = i`. By §2.2,
   every form-blind point has a 2-parameter family of circle models (any conjugate pair `ρ, ρ̄`). So `φ ≈ 0` at
   `ρ = i` is not evidence that the families "are not circle data"; only §2.2–2.3 carry content. The numbers
   themselves (`φ = 0` on C1, `≤ 0.027` on `E_t`) were not re-derived. They are tangential.
7. **Cosmetic: "equivalently".** "`(L_5^circ)`, or equivalently for this purpose `A_5(δ)`": Prop 3.1 proves only
   `A_5 ⇒ (L_5^circ)`. Both suffice, but they are not equivalent. The unrestricted `A_5` is strictly stronger
   (it covers unbalanced profiles, down to the Cilleruelo–Córdoba exponent `2/5`).
8. **Attribution.** The model (`γ̂` on `M̄_{0,k+2}`, isotropic points `P, Q`) and the observation that
   balanced `k`-subsets of large clusters converge to `γ̂` are in `threshold/L8.md` §6.1, together with the
   pointer "Gaussian analogue of Theorem D … would need 6 points … writing it out would need the archimedean
   `D_PQ`-height … I did not do this". The gcd-reduced balanced tuples are `ptolemy` Step 1 (item 532). The
   genuinely new content is: writing the reduction out, the archimedean divisor `D = E` (6 → 5), Theorem 4.1, and
   the circle version of round 18 item 7.

## 2. Re-derivation of the reduction (Thm 3.4 / 3.4′): correct

All quantifiers were re-checked by hand. The constraints of the project owner are respected: `Z`, `R_0`, `c`,
`η_0` and `w_0` are absolute; no gcd is discarded; there is no phase; and no fixed-`S` exceptional set is used.

* **Lemma 1.2.** `z_a - z_b = g(ζ_a - ζ_b)` with `ζ_a ≠ ζ_b` Gaussian integers gives `|g| ≤ L`, hence
  `R'' ≥ R/L`. Exponents:
  * `L/|g| ≤ e^{W/4 - W_∅/2}` and `R''^{1/2-δ} = e^{(W-W_∅)(1/4-δ/2)}`;
  * so the comparison is `-W_∅/4 ≤ -(W-W_∅)δ/2`, i.e. `δ ≤ W_∅/(2(W-W_∅))`;
  * and `W_∅ ≥ θW` gives `≥ θ/(2(1-θ)) ≥ θ/2`.

  Correct.
* **Lemma 3.2.** The gcd of the tuple is `∏ π^{min a}π̄^{e - max a}`, and `min a + (e - max a)` counts the
  layers on which the tuple lies on one side. Normalisation gives `1 ≤ S_{p,t} ≤ M-1`.
  * `φ(x) = C(x,k)` on `[k-1,∞)`, and `0` below, is convex: all roots of `C(x,k)` are `≤ k-1`, and the slope
    is `≥ 0` at `k-1`.
  * `∏_{j<k}(1 - j/(M-j)) ≥ 1 - k(k-1)/(2(M-k))`.
  * Reverse Markov: `E X ≤ θ + (1-θ)P(X ≥ θ)`.
  * For `M ≥ k + 2k(k-1)`: `ρ ≥ (3/2)2^{-k}`, so `P ≥ 2^{-k-1}`.

  Exact check over `k ≤ 8`, `M < 400`: no violation of either inequality (`gcdred_own.py`).
* **Lemma 3.3.** On the circle, `u = x+iy` identifies the points with a subset of `C^*`. `[z]` is the quotient
  `(C^*)^k → (C^*)^k/C^*`, which is surjective, so the preimage `Z'` is proper closed and independent of `N`.
  Grid points of a hypersurface of degree `d` number `≤ d M^{k-1}` (Schwartz–Zippel). Correct, and uniform in
  the circle.
* **Theorem 3.4.**
  * Choose `M ≥ M_k(Z)`. The tuple gets `W_∅ ≥ 2^{-k}W` and `[z] ∉ Z`.
  * Lemma 1.2 gives the reduced tuple at `R''^{1/2-δ}` with `R'' ≥ √R_1`.
  * `A_k` gives `R'' < R_0`, hence `R_1 < R_0^2` and `M ≤ R_0 + 1`.

  Correct.
* **Theorem 3.4′.**
  * Step 1 is item 532's balanced selection. Its balance comes from the cut identity
    `Σ log p (S/M - 1/2)^2 ≤ W/(4M) + O(log C)`, i.e. from `(UB)`. It holds for any profile, once restricted to
    one unit class (ptolemy Step 0).
  * Residue bound: the average of `log|c_ab|` is `≤ W/(4(M-1)) + log C`, and `|t_ab| = |c_ab|/2`. Markov at
    level `4E` gives `H_Σ ≤ k(k-1)W/(2(M-1))` on at least 3/4 of tuples. Together with balance (at least 1/2)
    and Lemma 3.3 (`O(M^{k-1})`), a positive proportion of tuples has all three properties.

  The proposal's constant `k(k-1)W/(2(M-1))` is exactly the level-`4E` Markov level, so it is consistent.
* **Prop 3.1.** With `m = k-1`, the exponent is `1 - 2^{m-1}(1-η_0)/((2^m-1)(1+η_0)) + 2c/((2^m-1)(1-η_0))`.
  At `η_0 = c = 0` this is `(2^{m-1}-1)/(2^m-1) = 1/2 - 1/(2^k-2)`, which is `7/15` at `k = 5`. Correct.
* **Monotonicity `A_5(δ) ⇒ A_6(δ)`.**
  * Take a 5-subtuple with gcd `g'`. Its diameter is `≤ R^{1/2-δ}/|g'| ≤ (R/|g'|)^{1/2-δ}`.
  * `R/|g'| ≥ R^{1/2+δ}`, since `|g'| ≤` diameter.
  * The pull-back of `Z_5` is proper in `T_6`.

  Correct.

**Conclusion:** (P) is a correct, elementary reduction. That the reduction is correct is not in question. What is
in question is whether the statement it reduces to is any easier than the target (§6).

## 3. Recomputation of numbers (own code, `kill_gauss_checks/`)

| claim in gauss.md | own script → output | result |
|---|---|---|
| `Φ` and the 4th-root formula (k4phi) | `phi_check.py` (sympy, symbolic) | `Im e_1 = Im e_2 = 0` identically; `Im e_3 = Φ/(v_1+v_2+v_3)`, with `Φ` exactly as printed. **Confirmed** |
| instance `X=-3, m=3`: roots and `Im P_x = -420, -180, -240` | `k4_instance.py` | all 7 roots reproduce; `Im P_x` is constant in `s` with these values. **Confirmed** |
| Jacobian rank 3 of `(X,m,s) ↦ T_4` | `k4_instance.py` (exact over `Q`, coordinates `Re P_x/Im P_x`) | rank 3, nonzero exact determinants at 6 points, including `(X,m) = (2,5), (1/2,3/2), (-5,1), (7,2)`. **Confirmed** |
| gcd norm `≤ 10368000`; `η → 0`; ratio in `[8.3, 43.4]`; `|t_ij| ≤ 45` | `k4_lattice2.py` → `.out` (311 values of `s` up to `10^30`) | max `Norm g = 10368000` ✓. `η = 0.033, 0.013, 0.0065` at `s = 10^6, 10^15, 10^30` (decays like `1/log s`) ✓. **Ratio range `[8.31, 142.28]`**, not `[8.3, 43.4]` (§1.2 item 3). Max `|t|` over all pairs, anchor edges included, is 240 (`≤ |Y_x|`) ✓. **408/1866 pairs have non-integral `t`** (§1.2 item 2) |
| `K.γ̂ = -3/8, -3/16, 1/32, 17/64` (`k = 4..7`); `(K+D_PQ).γ̂ = -1/4, -1/8, 1/16, 9/32, 33/64` (`k = 3..7`); count `= -(K+D_PQ).γ̂·2^{k-1}`; margins `1/15, 0.29, 0.52` | `kclass_own.py` (`K = Σ(s(n-s)/(n-1)-2)D_S`, `n = k+2`) | all reproduce exactly, `k = 3..9`. The direct count `W'(1-(k-1)(1/4+δ_k/2))` equals `2^{k-3}(5-k)-1`. `κ(δ_k) = 1/15, 9/31, 11/21`. **Confirmed** |
| `ρ_{10,5} = 0.0079`; `E[W_∅/W] = 0.0235`; `R'' ≥ √R`; reduced exponents `0.61–0.70` on the 10-point cluster | `gcdred_own.py` (all lattice points of `n = 1176852625`, best 10-window, chord diameter `8.2386√R`) | `0.007937`, `0.023455`, always, `0.606–0.702`. **Confirmed** |
| corollary: 3-subtuples of the Thm 4.1 family | `k3_projection.py` | each 3-subtuple, after its own gcd, has `diam/R_3^{1/3}` constant (60.5, 103.4, 105.8, 181.0). So `A_3(δ)` is false for `δ < 1/6` (dense, by projection of a dominant map) |

Not re-derived (tangential to the verdict): the C1/`E_t` φ-values of §2.1 (see §1.2 item 6 for why they carry no
weight).

## 4. Theorems 4.1–4.3 checked

**Theorem 4.1 (survives).**
* `z_x/z_0 = P_x/P̄_x`, and `|z_x - z_0| = |z_0|·2|Im P_x|/|P_x|`. With `Im P_x` constant, this is
  `≍ s^{7-4} = s^3` against `R'' ≍ s^7`. The diameter is `≤ 2 max_x |z_x - z_0|`. So the exponent is `3/7`.
* The common gcd is bounded, because the `z_x` are coprime in `Q(i)[s]`. Balance: all 7 cut norms are `≍ s^2`.
* Dominance was checked exactly (§3). The Zariski-closure argument (the set `B` and finitely many `s` per
  `(X,m)`) is correct.

One point needs care, and it is fine. `A_4(δ)` asks for all distances `≤ R''^{1/2-δ}` **for large `R''`**. With
`diam ≤ C R''^{3/7}` and `δ < 1/14`, this holds once `R''^{1/14-δ} > C`, i.e. for all large `s`. Correct.

**Dimension count, generalised (proved by parameter count; expected dimension only).** Consider a one-variable
polynomial block family at `k` points with block degrees `d_T`, `D = Σ d_T`, and exponent `e`. Its expected
moduli dimension is `≤ D(2 - (k-1)(1-e)) + k - 3`. Being dense in `T_k` needs `≥ k-2`, i.e.
`D(2-(k-1)(1-e)) ≥ 1`. At `k = 5` that needs `e > 1/2`, **for every profile, balanced or not**. So the
Theorem 4.1 mechanism cannot be rescued at `k = 5` by unbalanced degrees. Excess-dimensional solution sets are
excluded rigorously only in bounded degree, by Theorem 4.3(ii).

**Theorem 4.2 (survives, conditional).**
* `K_{Y_k} = -π^*B + (k-2)E` (blow-up of a point in dimension `k-1`).
* `h_E ≥ λ_{E,∞} - O(1)`, and `λ_{E,∞} ≥ (1/2+δ)log R'' = (1/4+δ/2)W'`. So
  `λ_{E,∞} + h_K ≥ ((k-1)(1/4+δ/2) - 1)W' = κW'`. At `k = 5`, `κ = 2δ`.
* Normalisation check over `Q`: at a split `p` (embedding via `π`) the valuation vector is the chain
  `Σ_t 1_{A_{p,t}}`, so `λ_{E_A+E_{A^c},p} = #{t : A_{p,t} ∈ {A, A^c}}·log p`. Summing over `p` gives
  `log n_A` (Lemma 1.1 holds over `Q`).
* Hand check of the class: at `k = 4`, `K.γ̂ = (1/4)(-2/5) + 2^{-4}·(-4.4) = -3/8`.

So `V_E(5,ε)` predicts `A_5(δ)` for **every** `δ > 0`, consistent with the heuristic count `e^{-2δW'}`, which is
profile-independent (§6.3). At `k = 4`, `κ > 0` iff `δ > 1/6`, consistent with Cilleruelo–Córdoba and with
Theorem 4.1. Correct as a conditional statement.

**Theorem 4.3 (survives, with (i) narrowed).**
* (ii) uses free-curve nefness. `F_ν ⊂ ν^*T_Y` has colength `μ_*` per fixed point. `dh` is surjective because
  `E` is smooth. Generic smoothness of the dominant evaluation map makes `F_ν` generically generated, hence nef.
  This gives `#fixed·μ_* + (k-2)#orbit·μ_* ≤ B.C`. At `k = 5`, the rational case gives `4(1/4+δ/2) > 1+ε` for
  `ε < 2δ`.
* The "dense union ⇒ one dominating component with a dense set of good members" step is correct. There are
  finitely many components at bounded `B.C`, and `μ_*` is uniform on a component.
* (i) for a single curve: see §1.2 item 4.
* Calibration at `k = 4` (`B.C = 14`, `μ_0 = 4`, `12 ≤ 14`): matches the family.

## 5. Barrier checklist, item by item (applied to the needed statement `(L_5^circ)` / `A_5(δ<1/30)`)

1. **Inputs that fakes satisfy.** The reduction uses only the gcd/layer identity, which fakes satisfy. So it
   proves nothing by itself, as the proposal says. The needed statement is **violated by fakes**: balanced
   5-row angle models with every monomial gap and diameter `≍ R''^{7/15}` exist, since the angles are continuous
   and `|t_ab| ≥ 1` is the pair gap. Any proof must therefore use the genuine Gaussian angles **at every prime
   of `N''` jointly**. That is the global, multi-place input that the record (overview §7.6) says no known
   method supplies. **Not a mechanism; a restatement of what is missing.**
2. **Averaged statements.** Not used. The statement is pointwise.
3. **Single-place auxiliary polynomials.** **Capped, exactly.** By Lemma K1 (§6.1) the forced-divisibility
   Liouville criterion for `ζ` is the same as for `z`. The `k`-row profile of a balanced cluster is the cube
   profile on `k` rows. `construct.md` Thm C2/C3 then gives ratio `≤ 1` for every form. Ru–Vojta at
   `S = {∞}` with one irreducible divisor: `β(L,E) ≤ τ(L,E)`, so it is no better than Liouville off a closed
   set. Both directions of the gauss §4.4 claim check out.
4. **Single linear form.** **Capped, exactly.** For an anchor edge at balanced `k = 5`:
   * the form is `Λ_x = Σ_{T∋x} 2 arg H_T - (π/2)k_x`;
   * it must satisfy `|Λ_x| ≲ R''^{-1/2-δ} = e^{-(1/4+δ/2)15w}`, which is `e^{-4w}` at `δ = 1/30`;
   * single-form LW (product form, naive heights `h(H_T/H̄_T) = (1/2)log n_T`, 8 blocks) gives only
     `|Λ_x| ≥ e^{-(1+ε)4w}`.

   That is tight at `ε = 0` and useless for `ε > 0`. SLW is false, and `SLW_w` is Vojta-strength. Not used.
5. **`S`-integer arguments with `|S|`-dependent thresholds.** Not used. Theorem 4.2 has `S = {∞}`. Any `S`-unit
   or Subspace attack on `(L_5^circ)` would have `S` equal to the primes of `N''` (unbounded): lit18 Props A/B
   apply, as `L8.md` §7 notes for `(L_8)`.
6. **Route 5 / CIRCLE(q).** Not used.
7. **`(L_5)/(L_6)`-type family mechanisms.** The proposal proves the circle analogue for curves (Theorem 4.3).
   * The density mechanism is excluded only heuristically (§1.2 item 5).
   * As in round 18 (correction M1), the following are **not excluded**: finitely many A-points on each member
     of a covering bounded-degree family; infinitely many non-covering curves; curves of unbounded degree;
     covering families of higher genus with finitely many A-points per member; sporadic points.

   Computational tests: a balanced 5-tuple needs 15 comparable cut norms. At one prime the threshold sets form a
   chain of proper nonempty subsets of `[5]`. Such a chain has length `≤ 4`, so one prime carries at most 4 of the
   15 cuts, and balanced data need at least 4 distinct split primes. The predicted count is `e^{-w}` per scale
   (**heuristic**). A census returning nothing would therefore say nothing. The `N ≤ 10^{12}` scan is evidence
   only.
8. **abc unit pencils.** Not used.
* **`hardness.md` (exact Dirichlet/Liouville exponent).** **Applies, in sharpened form** (§6.2):
  * the needed statement says Liouville is **not attained densely** at a **rational** point;
  * Theorem 4.1 shows that Liouville **is** attained densely at `k = 4`;
  * on `P^1`, Liouville is sharp at every rational point;
  * Roth/Schmidt beat Liouville only for irrational targets.

**Owner's constraints.** Respected by the reduction: no fixed-`S` exceptional set, no discarded gcd, no
phase-to-zero, no template-dependent threshold used as uniform.

## 6. Why the needed statement is capped

### 6.1 Lemma K1 (gcd rescaling is invisible to every method in the record; proved)

Let `z_1..z_k` be lattice points of `|z|^2 = N`, `g = gcd_{Z[i]}(z_j)`, and `ζ_j = z_j/g`. Let
`F ∈ Z[i][u_1..u_k, v_1..v_k]` be bihomogeneous of bidegree `(a,b)`, evaluated at `u_j = z_j`, `v_j = z̄_j`.
Every form on `X_k`, and every Laurent polynomial on `T_k` after homogenising, is a sum of such pieces, so the
argument applies piece by piece. Then:

* (i) **Values.** `F(z) = g^a ḡ^b F(ζ)`.
* (ii) **Forced divisors.** At every Gaussian prime `π`, `v_π(u_j(z)) = v_π(g) + v_π(ζ_j)` and
  `v_π(v_j(z)) = v_π(ḡ) + v_π(ζ̄_j)`. So the tropical forced divisor (`construct.md` Thm C1) satisfies
  `d_F(z) = g^a ḡ^b d_F(ζ)`.
* (iii) **Archimedean contact.** The angular width `Δ` of the tuple is the same for `z` and `ζ`, so
  `|F(z)| ≤ c|z|^{a+b}Δ^{ord}` holds iff `|F(ζ)| ≤ c|ζ|^{a+b}Δ^{ord}`.
* (iv) **Pair data, residues, block angles.** These are unchanged: `d_ab` is the separating-prime weight, and
  `t_ab = (ζ_a - ζ_b)/(2 gcd(ζ_a,ζ_b))` up to a unit.

*Proof.* (i) and (ii) are direct computations. (iii) holds because multiplication by `g` is a similarity.
(iv) holds because `gcd(z_a,z_b) = g·gcd(ζ_a,ζ_b)`. `□`

*Consequence.* For every auxiliary form, the Liouville contradiction "`F(z) ≠ 0` and `|F(z)| < |d_F(z)|`" holds
for `z` iff it holds for `ζ`. The same is true of every pair, cut, character and single-form inequality, since
these see only angles and separating primes. So the passage from "arc `√R`, critical" to "arc `R''^{7/15}`,
subcritical" is a change of the unit in which `R` is measured. For a 5-subtuple of a balanced cluster, the
`k`-row profile is the cube profile, and `construct.md` Cor C3 caps every form at ratio exactly `1`, before and
after the reduction alike. Only two things see the reduction:

* the **heuristic count** of classes, which becomes `e^{-w}` at `k = 5`;
* **Vojta's `D = E` term** (Theorem 4.2).

Neither is a proof method.

### 6.2 Where `A_5` sits in approximation terms

Approximation to `e ∈ T_5(Q)` (a rational point) with height `h_B = W'`: an A-point is
`dist(y,e) ≤ H_B(y)^{-(1/4+δ/2)}`.

* **Dirichlet / counting** (**heuristic**). Points of height `≤ X` within that distance number about
  `X·X^{-4(1/4+δ/2)} = X^{-2δ}`, uniformly in the profile. So `δ = 0` is the counting-critical exponent.
* **Liouville** (pair bound, equal to the forced-divisibility optimum). It excludes balanced points with
  `δ > 1/30` and arbitrary points with `δ > 1/10` (Cilleruelo–Córdoba). For `δ < 1/30` it is silent, even on
  balanced points.
* **Needed:** non-density at some `δ < 1/30` (balanced). In words: *the Liouville exponent at the rational point
  `e` is not attained by a Zariski-dense set of balanced points*.

Comparison:

* At `k = 4` exactly this fails. Theorem 4.1's dense family sits at the Liouville exponent `3/7` with bounded
  residues.
* On `P^1`, Liouville is attained at every rational point.
* The tools that beat Liouville (Roth, Schmidt, Ru–Vojta with several components) need an irrational, or
  Galois-split, target.
* For a rational point of a higher-dimensional variety, the only predictions beyond Liouville/Seshadri are
  conjectural: Vojta, and McKinnon's conjecture that best approximations lie on rational curves. As far as I
  know, no unconditional theorem improves Liouville at a rational point of such a variety for a Zariski-dense set
  (**not literature-checked this round**: no web search).

This is `round4/hardness.md`'s obstruction, restated for the reduced problem. It is not circumvented.

### 6.3 Is the needed statement at least as strong as the target, or as Vojta?

* `(L_5^circ) ⇒` uniform (Theorem 3.4′), so it is at least as strong as the target, as any sufficient condition
  is. The converse is not known: the uniform bound says nothing about the density of 5-point classes.
* It is implied by `V_E(5)`, not conversely. A refutation of `(L_5^circ)` would refute `V_E(5)`.
* The unrestricted `A_5(δ<1/30)` is strictly harder. On unbalanced profiles (for example, weight only on the ten
  2|3 cuts) Liouville allows `δ` up to `1/10`, so `A_5` must exclude points that Liouville permits by a wide
  margin.
* Against `(L_8)`, the heuristic margin is thinner: `e^{-w}` versus `e^{-10w}` per scale, and `κ = 2δ` versus
  `K.β = 10`. The proposal's "status: same as `(L_8)`" is accurate. Its "with increasing heuristic margin"
  remark (`k = 6, 7` better than `k = 5`) is also accurate.

No refutation was found. One-variable polynomial families cannot be dense at `k = 5` for any profile
(expected dimension, §4). Bounded-degree covering rational families are excluded (Theorem 4.3(ii)). Sporadic or
unbounded-degree supplies are exactly as unconstructible as for `(L_8)`.

## 7. What survives, stated exactly

1. **(proved; elementary)** For `k ≥ 3`: if `A_k(δ)` holds for some `0 < δ < 1/(2^k-2)`, then `M(C) < ∞` for
   every `C`. With `δ ≤ 2^{-k-1}` the proof needs only the convexity averaging of Lemma 3.2; otherwise item 532's
   balanced selection, restricted to one unit class.
   `A_k(δ)` reads: there are a proper Zariski-closed `Z ⊂ T_k` and an `R_0` such that every primitive `k`-tuple
   of distinct lattice points of a circle of radius `R ≥ R_0` with all distances `≤ R^{1/2-δ}` has `[ζ] ∈ Z`.
   Also `A_k(δ) ⇒ (L_k^circ)`, and `(L_k^circ) ⇒` uniform. Here `(L_k^circ)` must be stated **within one unit
   class**, or with `|t_ab| := |c_ab|/2` real and the bound read as `H_Σ ≥ cw - O(1)`. The reduction is
   invisible to every unconditional method (Lemma K1). It is a normalisation of `ptolemy` Steps 0–5 with the
   two isotropic marks of `L8.md` §6.1 restored.
2. **(proved, with exact rank input; reproduced)** `A_4(δ)` is false for every `δ < 1/14`, and `(L_4^circ)` is
   false for every `c, η_0, Z` (Theorem 4.1). Corollary (here, `k3_projection.py`): `A_3(δ)`, `δ < 1/6`, and
   `(L_3^circ)` are false, because the 3-subtuples sit at exponent `1/3` and the projection is dominant.
3. **(conditional on `V_E(k,ε)`, Vojta for `(Bl_e L̄_k^tw, D = E, S = {∞})`)** `A_k(δ)` whenever
   `(k-1)(1+2δ)/4 > 1`. At `k = 5` this holds for every `δ > 0`, so `V_E(5)` implies the uniform theorem. This is
   a new Vojta instance on a 4-fold. It does not improve on the existing conditional results as a route; it adds
   a fourth conditional form.
4. **(proved)** Theorem 4.3(ii), as stated. Theorem 4.3(i) narrowed to: genus `0`, `O` real algebraic, with the
   Liouville or Roth lower bound on `μ_0`.
5. **Open, capped, not a route:** `(L_5^circ)`, `(L_6^circ)`, `(L_7^circ)`. Each is Liouville-critical with a
   rational target, violated by fakes, predicted only by Vojta, and has no method in either direction.

**Net effect on the record.**
* Growth rate and uniform bound: unchanged.
* The barrier map gains one cap (Lemma K1 + construct C3: the gcd-reduced `k`-point statements are exactly
  Liouville-critical for every `k`), one refutation (`k ≤ 4`), and one conditional instance (`V_E(5)`).
* Recommended status line: "`(L_5^circ)`: open; Vojta-predicted (`κ = 2δ`); Liouville-tight by invariance; same
  standing as `(L_8)`".
