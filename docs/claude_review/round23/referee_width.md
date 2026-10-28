# Referee report on round23/width.md (adversarial)

## 0. Verdict

**The main result survives.** Theorem 2.5 (`t_max(m) = reg(S_m) = ⌊10m/3⌋` for every `m ≥ 1`, hence
`τ_e(-K_{LM_5}) = 10/3`, split and twisted) is correct. I re-derived every step with all quantifiers and found no gap.
The same holds for the `LM_4` / `dP_6` statements (`t_max(m) = 3m`), Corollary 3.1, the pair-box reduction, and the
consequence `inf_F δ_F = 1/10` (via `toric.md` Thm 3.3, whose last step I re-derived independently).

I found one false numerical statement, which is **minor**. It sits in a remark labelled *derived* (§3.2, echoed in
`width_notes.md`) and plays no part in any proof. Slicing `S_m` by `x_1`, bounding each slice with C_bal(4), and peeling
gives `7m/2`, not `15m/4`. The value `15m/4` comes from the refereed Theorem B(4) slice bound `(7/4)(a+b)`. The
qualitative conclusion (that the mean-width template cannot reach `10m/3`) still holds, because `7m/2 > 10m/3`. I also
note a few presentation points. None of them affects a claim.

Labels used below: *re-derived* means checked by hand with all quantifiers. *Exact/certified* means my own code (in
`referee_width_checks/`) with a different prime (`2^31 - 1`) from the author's (`2^61 - 1`), or exact integer arithmetic.
*Evidence* means a finite check that does not prove the statement.

## 1. Definitions and the geometric/elementary translation

* `P_{-K}(LM_5) = R_5 = {x ∈ M_R : Σ x_+ ≤ 1} = conv(e_a - e_b)`. The rays of the `A_4` fan are the classes of `1_A`,
  `∅ ≠ A ⊊ [5]`, and `-K` gives the polytope `{⟨x,1_A⟩ ≥ -1 ∀A} = {Σ x_- ≤ 1} = {Σ x_+ ≤ 1}` on `M`. So
  `S_m = mR_5 ∩ M = {Σ x = 0, Σ x_+ ≤ m} = T_5(m,m)`. Re-derived. The Ehrhart counts `21, 131, 471, 1251, 2751` were
  recomputed (`checks1.out`) and agree with `lit.md`.
* `toric.md` Lemma 3.1, `max ord_e f = reg(supp)`: re-derived. `ord_e f ≥ t` holds iff `c ⊥ Poly_{<t}|_S`, and such a
  `c ≠ 0` exists iff `t ≤ reg(S)`.
* `τ_e = lim t_max(m)/m = sup_m t_max(m)/m`: `t_max` is superadditive (orders add under products), and Fekete's lemma
  gives the equality. The pseudo-effective and big thresholds have the same supremum.
* Twisted form: a section of `mB` over `Q` is a `σ`-stable Laurent polynomial over `Q(i)`. The subspace of sections with
  `ord_e ≥ t` is `σ`-stable (`e` is rational), so by descent it has a nonzero `Q`-point whenever it is nonzero over
  `Q(i)`. `width.md` states this only as "the same `H^0(mB)` over `Q(i)`, `reg` field-independent", and the descent step
  is in `toric.md` §3. Minor clarity point, not a gap.
* The heights/distance/balance dictionary (`γ = 4/(1+2δ)`) is not used in Theorem 2.5. It enters only in §5, where
  `10/3 ↔ δ = 1/10` and `15/4 ↔ δ = 1/30` are arithmetically correct. Cilleruelo–Córdoba at `k = 5`:
  `1/(4⌊5/2⌋ + 2) = 1/10`. I checked this against the statement quoted from Cilleruelo–Granville (arXiv math/0608109,
  Thm 13) in the search snippet.

## 2. Proof steps, re-derived

**Lemma 2.1 (slices).** Correct.
* `c ≥ 0`: the slice is `T_{n-1}(a-c, b)`.
* `c < 0`: the slice is `T_{n-1}(a, b-|c|)`, using `Σ z_- ≤ b - |c|`.
* The slice is empty outside `[-b, a]`.
* Every `T_{n'}(a',b')` with `a', b' ≥ 0` is nonempty, so no slice is lost and none is double-counted.
* `T_2(a,b)` has `a + b + 1` collinear points, so `reg = a + b`.

Exhaustive enumeration for `n ≤ 5`, `a, b ≤ 4` and every `c ∈ [-b-2, a+2]` matches (`checks1.out`).

**Lemma 2.2 (peeling).** Correct. This is the standard Newton-type interpolation.
* The induction gives `reg(S') ≤ max_{i≥2}(reg S_i + i - 2) ≤ d - 1`.
* Then `P = p + (φ - c_1) q` has degree `≤ d` and interpolates `g`.
* `reg` of a slice may be computed inside its hyperplane: polynomials of degree `≤ d` on an affine subspace extend to
  polynomials of the same degree.
* Ordering the slices by decreasing upper bounds `ρ_i` gives `reg ≤ peel({ρ_i})`.

The criterion (2.2) is also correct: if `N_V(t) ≤ F + 1 - t` for all `t ∈ V`, then `i ≤ N_V(v_(i)) ≤ F + 1 - v_(i)`.

**Lemma 2.3 (recursion step).** Correct. `(R_n) ⇔ (R'_n)` by `F_{n-1}(a-p,b) - (a-p) - (b-q) = H_{n-1}(a-p,b) + q`, and
symmetrically for the other term. The three cases follow.
* Both sets nonempty: `(a*, b*) = (α, β)`, with `β ≤ b - 1 ≤ b`.
* Only the first set: `(α, b)`, i.e. `q = 0`. By monotonicity in the first argument the minimum is `F_{n-1}(α, b)`.
* Only the second set: `(a, β)`, i.e. `p = 0`.

In each case `N_V(t) ≤ F_n + 1 - t`. Edge cases `a = 0` or `b = 0` are fine: one block of `V` is empty and the ranges
for `p, q` degenerate.

**Lemma 2.4.** Correct.
* Symmetry: `(a,p) ↔ (b,q)` sends `X ↦ H_{n-1}(b-q, a) + p = Y`, because `H_{n-1}` is symmetric.
* `n = 3`: `min(X,Y) ≤ min(b, a)`.
* `n = 4`: `X + Y ≤ (a - p + q) + (b - q + p) = a + b`, `X ≤ 2b`, `Y ≤ 2a`.
* `n = 5`, `a ≤ b`:
  * `Y ≤ 2 min(a, b-q) + p ≤ 2a + p ≤ 3a`.
  * `X ≤ 2(a-p) + q` and `Y ≤ (a + b - q)/2 + p`, so `X/3 + 2Y/3 ≤ a + b/3`.
  * Integrality then gives `min(X,Y) ≤ a + ⌊b/3⌋ = min(a,b) + ⌊max(a,b)/3⌋`.

Each upper bound used for `H_4` is one of its two defining terms, so it is valid. Monotonicity of `F_n` holds because
`min`, `max` and `⌊(x+y)/2⌋` are each nondecreasing in `x`.

**Theorem 2.5.** Correct.
* Induction `n = 2 → 5` on all `(a, b)` simultaneously. The slice multiset of `T_n(a,b)` is exactly `V_n(a,b)`, with
  `F_{n-1}` replaced by the true `reg ≤ F_{n-1}`, and `peel` is monotone.
* `F_5(m,m) = 2m + min(m + ⌊m/3⌋, 3m) = 3m + ⌊m/3⌋ = ⌊10m/3⌋`.
* Lower bound: `V^s V_{123}^r`, `m = 3s + r`, has order `10s + 3r = ⌊10m/3⌋`.
* Hence `τ_e = sup ⌊10m/3⌋/m = 10/3`, attained at `m = 3`.
* `LM_4`: `F_4(m,m) = 3m`, and `V_4^s V_{123}^r` with `m = 2s + r` attains it. `V_4` has widths `|A||A^c| ≤ 4`, so
  `V_4 ∈ H^0(2B)`, with order 6.
* `dP_6`: `F_3(m,m) = 3m`, attained by `V_{123}^m`.

**Corollary 3.1 and the §3.2 reduction.** Correct.
* On `T_5(a,b)` every pair functional has range exactly `[-b, a]` (enumerated for `a, b ≤ 5`).
* `F_5 ≤ (5/3)(a+b)`, with equality iff `(a,b) = (3s,3s)`. Checked exhaustively on `[0,200)^2`.
* Pair-box reduction: the box contains `S` and has the same 10 widths. `reg` is monotone, and the 10 pair functionals
  span `M^∨ ⊗ Q` (e.g. `x_1 = (x_{12} + x_{13} - x_{23})/2`), so the box is finite.

**Remark 1.1 ((W) = C_bal at `M = 5`).** Correct against `threshold/construct.md` §2.8 and §2.7.
* `μ_bal` is uniform on the cuts with `|P| = ⌊M/2⌋`.
* `π_bal = 2·(2/5)(3/4) = 3/5`.
* `E w_P = (1/10)Σ_{|A|=2} w_A`, so `ord ≤ (5/3)(1/10)Σ = (1/6)Σ`.
* `M = 4`: `(3/2)(1/6)Σ = (1/4)Σ`.

## 3. Independent computations

All files are in `referee_width_checks/`. The code is independent of `width_checks/`, and the modular ranks use the
prime `2^31 - 1`. All runs were under `prlimit --as=3500000000`.

| script → output | what | result |
|---|---|---|
| `rec_check.py 45` → `rec_check.out` | monotonicity of `F_n`; `(R_n)` for all `0 ≤ p ≤ a`, `0 ≤ q ≤ b ≤ 45`; `peel(V_n(a,b)) ≤ F_n(a,b)` directly; exact recursion `ρ_n` | 0 failures; `max_{p,q} min(X,Y) = H_n` on every point (the bounds are tight); `ρ_n = F_n` on `[0,45]^2`; `F_5(m,m) = ⌊10m/3⌋`, `F_4 = F_3 = 3m` for `m < 10^4` |
| `checks1.py` → `checks1.out` | Lemma 2.1 by enumeration; `reg T_n(a,b)` against `F_n` (`n = 3: a ≤ b ≤ 6`; `n = 4: ≤ 5`; `n = 5: ≤ 3`); `reg T_5(2,3)`; `reg S_m`, `m ≤ 3` | slices match; full rank at `F_n` everywhere (certified upper bound); full already at `F_n - 1` only at `(2,3)`; `T_5(2,3)` is full at 7 and deficient at 6 (mod `p`, evidence); `S_1, S_2, S_3` are full at `3, 6, 10` and deficient one below (the lower bound is proved by Vandermonde in any case) |
| `checks2.py` → `checks2.out` | Vandermonde claims (exact, sympy, integer moments); slicing-template numbers | `V·(x_1⋯x_5)^{-2}` has max `Σ x_+ = 3`, `ord_e = 10`, widths 4 (`|A| = 1`) and 6 (`|A| = 2`); `V_{123}/(x_1x_2x_3)` has max `Σ x_+ = 1`, `ord_e = 3`; template peeling: see §5, issue 1 |
| `checks3.py` → `checks3.out` | Corollary 3.1 | confirmed (§2 above) |

All numbers in `width.md` §§2-4 that I recomputed agree, except the template value in §3.2 (§5, issue 1). The author's
own outputs (`treg5.out`: 14 cases `0 ≤ a ≤ b ≤ 4`, excluding `(0,0)`; `Rn_check.out`; `rho_verify.out`) say what
`width.md` reports.

## 4. Consequences, citations, constraints

**§5 item 1 (`inf_F δ_F = 1/10`).** Re-derived from `toric.md` Thm 3.3.
* Every element `f` of the semigroup generated by `F` is itself a Laurent polynomial over `Q(i)`.
* `f f^*` lies in `H^0(m'B)` with `m' = max_A width_f(1_A)`, because `Q_f - Q_f ⊂ m'R_5` and `R_5` is cut out by
  `⟨·, 1_A⟩ ≤ 1`. It has order `2 ord_e f`.
* So `2 ord_e f ≤ t_max(m') = ⌊10m'/3⌋ ≤ 10m'/3` (Theorem 2.5), i.e. `max_A width_f(1_A)/ord_e f ≥ 3/5`.
* This holds for every semigroup element, hence for the infimum, and `δ_F ≥ 3/5 - 1/2 = 1/10` for every finite `F`.
  Equality holds for `F = {V}` (`max_A width = 6`, `ord = 10`).

The quantifiers are stated correctly:
* the threshold `1/10` holds for every `F`;
* the profile witnessing it (the minimax `ℓ_F`) may depend on `F`, and `width.md` does not claim otherwise;
* profiles are probability vectors on the 15 cuts. Their realisability by actual circle points is not claimed: §5 item 6
  says the statement is a barrier for a method only.

**§5 item 3 (McKinnon–Roth Liouville at `e` gives exactly `10/3`).** Correct for `L = -K`. The section bound proves
`α ≥ r/m`, and `sup r/m = τ_e(-K) = 10/3`. Other toric `L`, compared through profiles, are covered by Thm 3.3.

**§5 item 4.** Correct: nothing about `α^ess > 15/4`, `A_5(δ < 1/30)` or `(L_5^∘)` is claimed. The balanced cap `15/4`
(construct.md Cor C3, attained by `V`) is unaffected.

**Side data (`ρ_6`, `ρ_7`, `LM_6`).** These are labelled evidence, and that label is right. The pull-back claim
`τ_e(-K_{LM_6}) ≥ 10/3` is correct: `V_5` has width `|A ∩ [5]|·|[5] \ A| ≤ 6`.

**Citations** (my own searches).
* McKinnon–Roth, *An analogue of Liouville's theorem and an application to cubic surfaces*, Eur. J. Math. 2 (2016)
  929–959, arXiv:1306.2977. Confirmed, including the asymptotic base-locus extension.
* Cilleruelo–Córdoba exponent `1/(4⌊k/2⌋+2)`. Confirmed via Cilleruelo–Granville arXiv math/0608109, Thm 13 (snippet).
* `width.md` cites nothing else externally. Internal references (`toric.md` Lemma 3.1 and Thms 3.3, 3.4(a);
  `construct.md` §2.4 and §2.8; `lit.md` Prop 1.3) were checked against the files and say what is attributed to them.
* Novelty: one search found no prior computation of the interpolation regularity of `(mΔ - mΔ) ∩ Z^5`. That is
  inconclusive (snippets only). `width.md` claims only "new to the record", which is appropriate.

**Owner constraints.**
* No Subspace exceptional sets.
* No gcd or denominator cost is dropped: Theorem 2.5 is pure linear algebra, and §5 inherits Thm 3.2, which keeps
  `N'^m`.
* No phase is involved.
* No template-dependent threshold is treated as uniform: the profile dependence on `F` is stated explicitly, and nothing
  depends on the radius.

The audit in `width.md` §6 is accurate.

## 5. Issues

1. **(minor, false derived remark.)** `width.md` §3.2 says: "slicing `S_m` by `x_1` with the C_bal(4) bound
   `(3/2)(2m - |c|)` on the slice `c` and peeling returns exactly `15m/4` (Theorem B's value)". This is false.
   * The second copy of slice `|c| = k` sits at peel index `2k + 1`. So the peel value is
     `max_k (3/2)(2m - k) + 2k = 3m + k/2`, which is `7m/2` at `k = m`.
   * With floors the values are `3, 7, 10, 14, 21, 28, 42, 84, 168` for `m = 1, 2, 3, 4, 6, 8, 12, 24, 48`, i.e. `7m/2`
     for even `m` (`checks2.out`).
   * The value `15m/4` (exactly `15, 30, 45, 90, 180` at `m = 4, 8, 12, 24, 48`) is what the refereed Theorem B(4)
     bound gives: `2^{-3}·14·(a+b) = (7/4)(a+b)`.
   * C_bal(4) is itself a conjecture, so it cannot serve as a slice bound in a proof anyway.
   * Conclusion unchanged: `7m/2 > 10m/3`, so the mean-width template still fails. The §2.4 Remark ("returns exactly
     Theorem B's `15m/4`" with mean-width slice bounds) is correct if it refers to Theorem B's own bound.
   * Fix: replace "C_bal(4) bound `(3/2)(2m-|c|)` … returns exactly `15m/4`" by "Theorem B(4) bound `(7/4)(2m-|c|)`
     returns exactly `15m/4`; even the conjectural C_bal(4) bound `(3/2)(2m-|c|)` returns `7m/2`".
2. **(minor, clarity.)** Theorem 2.5 part 1 for the twisted form should say that `σ`-descent (`toric.md` §3) makes the
   maximal order attained by a `Q`-rational section. The phrase "reg is field-independent" alone gives only the upper
   bound.
3. **(cosmetic.)** In §4, "`T_5(a,b)`, `0 ≤ a ≤ b ≤ 4` … all 14 cases" excludes `(0,0)` without saying so. In §2.2
   (Lemma 2.2's proof) the remark "`d ≥ reg(S_r) + r - 1 ≥ 1`" is unnecessary.
4. **(not an issue, recorded.)** The §3.2 *derived* coefficient-`1/2` analysis was not re-verified in full. It is
   labelled as a statement about one template and supports no claim.

No fatal or major issue was found.

## 6. Corrected statement

**Theorem (proved).** For all `a, b ≥ 0`: `reg T_n(a,b) ≤ F_n(a,b)` for `n = 2, …, 5`, where `F_n = a + b + H_n` with
`H_2 = 0`, `H_3 = min(a,b)`, `H_4 = min(⌊(a+b)/2⌋, 2min(a,b))` and
`H_5 = min(min(a,b) + ⌊max(a,b)/3⌋, 3min(a,b))`. Consequently:
* `t_max(m) = ⌊10m/3⌋` for every `m ≥ 1` on `LM_5` and on its twisted form `L̄_5^tw`, attained over `Q` (by descent),
  so `τ_e(-K) = 10/3`;
* `t_max(m) = 3m` on `LM_4` and `dP_6`;
* (W) holds on every `T_5(a,b)`, with equality exactly on `S_{3s}`.

With `toric.md` Thm 3.3: `δ_F ≥ 1/10` for every finite family `F` of auxiliary Laurent polynomials over `Q(i)` with toric
forced divisibility, and equality holds for `F = {V}`. The profile witnessing the cap may depend on `F`. General (W)
(C_bal at `M = 5`) remains open. Nothing changes for `A_5(δ < 1/30)`, `α^ess_e > 15/4` or `(L_5^∘)`.

§3.2 should be amended: with the conjectural C_bal(4) slice bound, slicing and peeling gives `7m/2`, not `15m/4`; the
`15m/4` comes from the Theorem B(4) bound.
