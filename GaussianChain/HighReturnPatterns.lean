import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Data.Fin.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

namespace GaussianChain
namespace HighReturnPatterns

open scoped BigOperators

/-!
# Abstract inequalities for high-return pattern descent

This file records only the elementary real and finite-sum inequalities used by
the high-return pattern argument.  In particular, it does not assert that a
Gaussian configuration satisfies their hypotheses.
-/

/--
Algebraic rearrangement of the max-cut inequality in the even-cardinality
case.  The coefficient `M ^ 2 / 4` in `hcut` is the exact balanced-cut
coefficient for even integral `M`.
-/
theorem even_maxCut_weight_cap
    {M H κ w : ℝ} (hM : 0 < M)
    (hcut :
      (M * (M - 1) / 2) * (H / 2 - κ) ≤
        (M ^ 2 / 4) * (H - w)) :
    M * w ≤ H + 2 * (M - 1) * κ := by
  have hfactor :
      0 ≤ (M / 4) * (H + 2 * (M - 1) * κ - M * w) := by
    calc
      0 ≤ (M ^ 2 / 4) * (H - w) -
          (M * (M - 1) / 2) * (H / 2 - κ) := sub_nonneg.mpr hcut
      _ = (M / 4) * (H + 2 * (M - 1) * κ - M * w) := by ring
  have hbracket : 0 ≤ H + 2 * (M - 1) * κ - M * w := by
    by_contra hneg
    have hbracket_neg : H + 2 * (M - 1) * κ - M * w < 0 := lt_of_not_ge hneg
    have hMfour : 0 < M / 4 := by positivity
    have := mul_neg_of_pos_of_neg hMfour hbracket_neg
    linarith
  linarith

/--
Algebraic rearrangement of the max-cut inequality in the odd-cardinality
case.  The coefficient `(M ^ 2 - 1) / 4` in `hcut` is the exact
balanced-cut coefficient for odd integral `M`.  The case `M = 1` is excluded:
there are then no pairs, so the max-cut inequality gives no weight bound.
-/
theorem odd_maxCut_weight_cap
    {M H κ w : ℝ} (hM : 1 < M)
    (hcut :
      (M * (M - 1) / 2) * (H / 2 - κ) ≤
        ((M ^ 2 - 1) / 4) * (H - w)) :
    (M + 1) * w ≤ H + 2 * M * κ := by
  have hfactor :
      0 ≤ ((M - 1) / 4) * (H + 2 * M * κ - (M + 1) * w) := by
    calc
      0 ≤ ((M ^ 2 - 1) / 4) * (H - w) -
          (M * (M - 1) / 2) * (H / 2 - κ) := sub_nonneg.mpr hcut
      _ = ((M - 1) / 4) * (H + 2 * M * κ - (M + 1) * w) := by ring
  have hbracket : 0 ≤ H + 2 * M * κ - (M + 1) * w := by
    by_contra hneg
    have hbracket_neg : H + 2 * M * κ - (M + 1) * w < 0 := lt_of_not_ge hneg
    have hMfour : 0 < (M - 1) / 4 := by positivity
    have := mul_neg_of_pos_of_neg hMfour hbracket_neg
    linarith
  linarith

/--
Summing the six pair lower bounds among four weights counts each weight three
times.  This is the abstract scalar four-body gain, independent of any
Gaussian interpretation.
-/
theorem four_weight_sum_lower
    {h₁ h₂ h₃ h₄ H L : ℝ}
    (h₁₂ : H / 2 + L ≤ h₁ + h₂)
    (h₁₃ : H / 2 + L ≤ h₁ + h₃)
    (h₁₄ : H / 2 + L ≤ h₁ + h₄)
    (h₂₃ : H / 2 + L ≤ h₂ + h₃)
    (h₂₄ : H / 2 + L ≤ h₂ + h₄)
    (h₃₄ : H / 2 + L ≤ h₃ + h₄) :
    H + 2 * L ≤ h₁ + h₂ + h₃ + h₄ := by
  have hA := add_le_add h₁₂ h₁₃
  have hB := add_le_add h₁₄ h₂₃
  have hC := add_le_add h₂₄ h₃₄
  have hAB := add_le_add hA hB
  have hABC := add_le_add hAB hC
  linarith [hABC]

/--
If the same four weights have total at most `H`, the six pair bounds force
the pair gain `L` to be nonpositive.
-/
theorem four_weight_pair_gain_nonpos
    {h₁ h₂ h₃ h₄ H L : ℝ}
    (htotal : h₁ + h₂ + h₃ + h₄ ≤ H)
    (h₁₂ : H / 2 + L ≤ h₁ + h₂)
    (h₁₃ : H / 2 + L ≤ h₁ + h₃)
    (h₁₄ : H / 2 + L ≤ h₁ + h₄)
    (h₂₃ : H / 2 + L ≤ h₂ + h₃)
    (h₂₄ : H / 2 + L ≤ h₂ + h₄)
    (h₃₄ : H / 2 + L ≤ h₃ + h₄) :
    4 * L ≤ 0 := by
  have hlower :=
    four_weight_sum_lower h₁₂ h₁₃ h₁₄ h₂₃ h₂₄ h₃₄
  linarith

/-- Total mass of indices whose return score meets the threshold `a * m`. -/
noncomputable def highReturnMass
    {ι : Type*} [Fintype ι]
    (m a : ℝ) (returns weight : ι → ℝ) : ℝ :=
  ∑ i, if a * m ≤ returns i then weight i else 0

/--
If every return score is at most `m / 2`, its weighted sum is bounded by
the threshold contribution from the full mass plus the excess coefficient
times the high-return mass.
-/
theorem weightedReturn_le_baseline_add_highReturnMass
    {ι : Type*} [Fintype ι]
    {m a H : ℝ} {returns weight : ι → ℝ}
    (hweight : ∀ i, 0 ≤ weight i)
    (hreturn : ∀ i, returns i ≤ m / 2)
    (htotal : ∑ i, weight i = H) :
    ∑ i, returns i * weight i ≤
      a * m * H + (m / 2 - a * m) * highReturnMass m a returns weight := by
  classical
  have hpoint : ∀ i,
      returns i * weight i ≤
        a * m * weight i +
          (m / 2 - a * m) * (if a * m ≤ returns i then weight i else 0) := by
    intro i
    by_cases hi : a * m ≤ returns i
    · simp only [hi, if_true]
      have hmul := mul_le_mul_of_nonneg_right (hreturn i) (hweight i)
      nlinarith
    · simp only [hi, if_false, mul_zero, add_zero]
      have hlow : returns i ≤ a * m := le_of_not_ge hi
      exact mul_le_mul_of_nonneg_right hlow (hweight i)
  calc
    ∑ i, returns i * weight i ≤
        ∑ i, (a * m * weight i +
          (m / 2 - a * m) * (if a * m ≤ returns i then weight i else 0)) :=
      Finset.sum_le_sum fun i _ ↦ hpoint i
    _ = a * m * H + (m / 2 - a * m) * highReturnMass m a returns weight := by
      rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, htotal]
      rfl

/--
Abstract threshold-mass lower bound.  It is the finite weighted form of
combining a lower bound for total return mass with `returns i ≤ m / 2`.
-/
theorem highReturnMass_lower_of_weightedReturn_lower
    {ι : Type*} [Fintype ι]
    {m a H Rlower : ℝ} {returns weight : ι → ℝ}
    (hm : 0 < m) (ha : a < 1 / 2)
    (hweight : ∀ i, 0 ≤ weight i)
    (hreturn : ∀ i, returns i ≤ m / 2)
    (htotal : ∑ i, weight i = H)
    (hlower : Rlower ≤ ∑ i, returns i * weight i) :
    (Rlower - a * m * H) / ((1 / 2 - a) * m) ≤
      highReturnMass m a returns weight := by
  have hupper :=
    weightedReturn_le_baseline_add_highReturnMass
      (m := m) (a := a) (H := H) hweight hreturn htotal
  have hdenom : 0 < (1 / 2 - a) * m := by positivity
  apply (div_le_iff₀ hdenom).2
  have hcoeff : m / 2 - a * m = (1 / 2 - a) * m := by ring
  rw [← hcoeff]
  linarith

/-- Number of entries at least `W / 8`. -/
noncomputable def denseEdgeCount
    {m : ℕ} (S : Fin m → ℝ) (W : ℝ) : ℕ :=
  (Finset.univ.filter fun k ↦ W / 8 ≤ S k).card

/--
If `m` entries are individually at most `W` and have total at
least `mW/4`, then at least `m/7` of them are at least `W/8` (with the
cardinality comparison made in `ℝ`).  Nonnegativity of the entries, present
in the intended application, is not needed for this upper-counting argument.
-/
theorem denseEdgeCount_lower
    {m : ℕ} {S : Fin m → ℝ} {W : ℝ}
    (hW : 0 < W)
    (hS_le : ∀ k, S k ≤ W)
    (hsum : (m : ℝ) * W / 4 ≤ ∑ k, S k) :
    (m : ℝ) / 7 ≤ (denseEdgeCount S W : ℝ) := by
  classical
  have hpoint : ∀ k,
      S k ≤ W / 8 + (7 * W / 8) * (if W / 8 ≤ S k then (1 : ℝ) else 0) := by
    intro k
    by_cases hk : W / 8 ≤ S k
    · simp only [hk, if_true]
      have := hS_le k
      linarith
    · simp only [hk, if_false, mul_zero, add_zero]
      exact le_of_not_ge hk
  have hupper :
      ∑ k, S k ≤
        ∑ k, (W / 8 + (7 * W / 8) *
          (if W / 8 ≤ S k then (1 : ℝ) else 0)) :=
    Finset.sum_le_sum fun k _ ↦ hpoint k
  have hindicator :
      ∑ k, (if W / 8 ≤ S k then (1 : ℝ) else 0) =
        (denseEdgeCount S W : ℝ) := by
    simp [denseEdgeCount]
  have hupper' :
      ∑ k, S k ≤
        (m : ℝ) * (W / 8) + (7 * W / 8) * (denseEdgeCount S W : ℝ) := by
    calc
      ∑ k, S k ≤
          ∑ k, (W / 8 + (7 * W / 8) *
            (if W / 8 ≤ S k then (1 : ℝ) else 0)) := hupper
      _ = (m : ℝ) * (W / 8) +
          (7 * W / 8) * (denseEdgeCount S W : ℝ) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, hindicator]
        simp
  have hmass :
      (m : ℝ) * W ≤ 7 * W * (denseEdgeCount S W : ℝ) := by
    nlinarith [hsum, hupper']
  have hcancel : (m : ℝ) ≤ 7 * (denseEdgeCount S W : ℝ) := by
    apply le_of_mul_le_mul_left _ hW
    calc
      W * (m : ℝ) = (m : ℝ) * W := by ring
      _ ≤ 7 * W * (denseEdgeCount S W : ℝ) := hmass
      _ = W * (7 * (denseEdgeCount S W : ℝ)) := by ring
  linarith

end HighReturnPatterns
end GaussianChain
