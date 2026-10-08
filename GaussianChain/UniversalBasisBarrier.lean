import Mathlib

namespace GaussianChain
namespace UniversalBasisBarrier

/--
Two distinct degree-`d` monomials in three variables have common exponent
weight at most `d - 1`.

This is the elementary integral gap behind the obstruction to replacing the
point-dependent adapted basis by one universal basis or a convex mixture of
universal bases.
-/
theorem min_exponent_sum_le_degree_sub_one
    (a₀ a₁ a₂ b₀ b₁ b₂ d : ℕ)
    (ha : a₀ + a₁ + a₂ = d)
    (hb : b₀ + b₁ + b₂ = d)
    (hne : a₀ ≠ b₀ ∨ a₁ ≠ b₁ ∨ a₂ ≠ b₂) :
    min a₀ b₀ + min a₁ b₁ + min a₂ b₂ ≤ d - 1 := by
  have hlt : b₀ < a₀ ∨ b₁ < a₁ ∨ b₂ < a₂ := by
    by_contra h
    push_neg at h
    rcases h with ⟨h₀, h₁, h₂⟩
    have heq₀ : a₀ = b₀ := by omega
    have heq₁ : a₁ = b₁ := by omega
    have heq₂ : a₂ = b₂ := by omega
    rcases hne with hne₀ | hne₁ | hne₂
    · exact hne₀ heq₀
    · exact hne₁ heq₁
    · exact hne₂ heq₂
  have hm₀l : min a₀ b₀ ≤ a₀ := min_le_left _ _
  have hm₀r : min a₀ b₀ ≤ b₀ := min_le_right _ _
  have hm₁l : min a₁ b₁ ≤ a₁ := min_le_left _ _
  have hm₁r : min a₁ b₁ ≤ b₁ := min_le_right _ _
  have hm₂l : min a₂ b₂ ≤ a₂ := min_le_left _ _
  have hm₂r : min a₂ b₂ ≤ b₂ := min_le_right _ _
  rcases hlt with h₀ | h₁ | h₂ <;> omega

/--
If three valuation types have total paired weight at most `432`, at least one
of them has weight at most `144`.
-/
theorem one_of_three_le_144
    (w₀ w₁ w₂ : ℚ)
    (hsum : w₀ + w₁ + w₂ ≤ 432) :
    w₀ ≤ 144 ∨ w₁ ≤ 144 ∨ w₂ ≤ 144 := by
  by_contra h
  push_neg at h
  linarith

/-- The endpoint threshold for a 27-dimensional sextic system is `297 / 2`. -/
theorem universal_weight_below_endpoint :
    (144 : ℚ) < 27 * (6 - 1 / 2) := by
  norm_num

/--
Consequently, three paired weights with total at most `432` cannot all exceed
`27 * (6 - 1/2)`, the weight required to retain the endpoint margin.
-/
theorem not_all_above_endpoint
    (w₀ w₁ w₂ : ℚ)
    (hsum : w₀ + w₁ + w₂ ≤ 432) :
    ¬ (27 * (6 - 1 / 2) < w₀ ∧
       27 * (6 - 1 / 2) < w₁ ∧
       27 * (6 - 1 / 2) < w₂) := by
  intro h
  obtain h₀ | h₁ | h₂ := one_of_three_le_144 w₀ w₁ w₂ hsum
  · have hg := universal_weight_below_endpoint
    linarith
  · have hg := universal_weight_below_endpoint
    linarith
  · have hg := universal_weight_below_endpoint
    linarith

end UniversalBasisBarrier
end GaussianChain
