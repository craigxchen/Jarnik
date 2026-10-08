import GaussianChain.Basic

namespace GaussianChain
namespace RootedPell

/-- If a rooted-circle coordinate has the squarefree-style factorization
`k = a*x^2`, `s = a*x*y`, then the rooted circle equation becomes a
sum-of-two-squares representation `M = a*(x^2+y^2)`. -/
theorem norm_representation_of_rooted_factorization
    {M k s a x y : ℤ}
    (hcircle : s ^ 2 = k * (M - k))
    (hk : k = a * x ^ 2)
    (hs : s = a * x * y)
    (hax : a * x ^ 2 ≠ 0) :
    M = a * (x ^ 2 + y ^ 2) := by
  rw [hk, hs] at hcircle
  have hmul :
      (a * x ^ 2) * M = (a * x ^ 2) * (a * (x ^ 2 + y ^ 2)) := by
    nlinarith
  exact mul_left_cancel₀ hax hmul

/-- Two rooted points represented by `M = a_i(x_i^2+y_i^2)` satisfy one
inhomogeneous Pell-type equation. -/
theorem pell_difference
    {M a₁ a₂ x₁ x₂ y₁ y₂ : ℤ}
    (h₁ : M = a₁ * (x₁ ^ 2 + y₁ ^ 2))
    (h₂ : M = a₂ * (x₂ ^ 2 + y₂ ^ 2)) :
    a₁ * y₁ ^ 2 - a₂ * y₂ ^ 2 =
      a₂ * x₂ ^ 2 - a₁ * x₁ ^ 2 := by
  nlinarith

/-- Three rooted points give a simultaneous Pell system with the first
complementary coordinate shared between the two equations. -/
theorem simultaneous_pell_system
    {M a₁ a₂ a₃ x₁ x₂ x₃ y₁ y₂ y₃ : ℤ}
    (h₁ : M = a₁ * (x₁ ^ 2 + y₁ ^ 2))
    (h₂ : M = a₂ * (x₂ ^ 2 + y₂ ^ 2))
    (h₃ : M = a₃ * (x₃ ^ 2 + y₃ ^ 2)) :
    (a₁ * y₁ ^ 2 - a₂ * y₂ ^ 2 =
        a₂ * x₂ ^ 2 - a₁ * x₁ ^ 2) ∧
      (a₁ * y₁ ^ 2 - a₃ * y₃ ^ 2 =
        a₃ * x₃ ^ 2 - a₁ * x₁ ^ 2) := by
  exact ⟨pell_difference h₁ h₂, pell_difference h₁ h₃⟩

/-- Equal squarefree coefficients reduce the Pell equation to a difference of
squares. -/
theorem equal_coefficient_difference
    {a x₁ x₂ y₁ y₂ : ℤ}
    (ha : a ≠ 0)
    (h : a * y₁ ^ 2 - a * y₂ ^ 2 =
      a * x₂ ^ 2 - a * x₁ ^ 2) :
    y₁ ^ 2 - y₂ ^ 2 = x₂ ^ 2 - x₁ ^ 2 := by
  have hmul : a * (y₁ ^ 2 - y₂ ^ 2) =
      a * (x₂ ^ 2 - x₁ ^ 2) := by
    nlinarith
  exact mul_left_cancel₀ ha hmul

end RootedPell
end GaussianChain
