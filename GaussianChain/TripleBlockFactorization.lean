import GaussianChain.Basic

namespace GaussianChain
namespace TripleBlockFactorization

/-- The first coordinate of the universal three-vector determinant relation. -/
theorem det2_linear_relation_re
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    det2 a₂ b₂ a₃ b₃ * a₁ -
        det2 a₁ b₁ a₃ b₃ * a₂ +
        det2 a₁ b₁ a₂ b₂ * a₃ = 0 := by
  unfold det2
  ring

/-- The second coordinate of the universal three-vector determinant relation. -/
theorem det2_linear_relation_im
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    det2 a₂ b₂ a₃ b₃ * b₁ -
        det2 a₁ b₁ a₃ b₃ * b₂ +
        det2 a₁ b₁ a₂ b₂ * b₃ = 0 := by
  unfold det2
  ring

/-- The pattern-product identity behind the product of three chord norms.
The variables are indexed by subsets of a three-point set. -/
theorem triple_pattern_product_identity
    (n₀ n₁ n₂ n₃ n₁₂ n₁₃ n₂₃ n₁₂₃ : ℤ) :
    (n₀ * n₃ * n₁₂ * n₁₂₃) *
        (n₀ * n₂ * n₁₃ * n₁₂₃) *
        (n₀ * n₁ * n₂₃ * n₁₂₃) =
      (n₀ * n₁ * n₂ * n₃ * n₁₂ * n₁₃ * n₂₃ * n₁₂₃) *
        (n₀ * n₁₂₃) ^ 2 := by
  ring

/-- Exact product formula obtained when the three squared chords have their
shared-prime pattern factorizations. -/
theorem triple_chord_product_factorization
    {h N n₀ n₁ n₂ n₃ n₁₂ n₁₃ n₂₃ n₁₂₃
      δ₁₂ δ₁₃ δ₂₃ d₁₂ d₁₃ d₂₃ : ℤ}
    (hN : N = n₀ * n₁ * n₂ * n₃ * n₁₂ * n₁₃ * n₂₃ * n₁₂₃)
    (hd₁₂ : d₁₂ = 4 * h ^ 2 * n₀ * n₃ * n₁₂ * n₁₂₃ * δ₁₂ ^ 2)
    (hd₁₃ : d₁₃ = 4 * h ^ 2 * n₀ * n₂ * n₁₃ * n₁₂₃ * δ₁₃ ^ 2)
    (hd₂₃ : d₂₃ = 4 * h ^ 2 * n₀ * n₁ * n₂₃ * n₁₂₃ * δ₂₃ ^ 2) :
    d₁₂ * d₁₃ * d₂₃ =
      64 * h ^ 6 * N * (n₀ * n₁₂₃) ^ 2 *
        (δ₁₂ * δ₁₃ * δ₂₃) ^ 2 := by
  rw [hd₁₂, hd₁₃, hd₂₃, hN]
  ring

/-- A convenient rearrangement of the product factorization. -/
theorem triple_chord_product_rearranged
    {h N U δprod dprod : ℤ}
    (hfactor : dprod = 64 * h ^ 6 * N * U ^ 2 * δprod ^ 2) :
    dprod = 64 * h ^ 6 * N * (U * δprod) ^ 2 := by
  rw [hfactor]
  ring

end TripleBlockFactorization
end GaussianChain
