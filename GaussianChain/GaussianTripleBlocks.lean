import GaussianChain.Basic

namespace GaussianChain
namespace GaussianTripleBlocks

/-- The ordinary dot product of two integer vectors in the plane. -/
def dot2 (a b c d : ℤ) : ℤ := a * c + b * d

/-- Lagrange's two-square identity in dot/determinant form. -/
theorem norm2_mul_eq_dot2_sq_add_det2_sq
    (a b c d : ℤ) :
    norm2 a b * norm2 c d =
      dot2 a b c d ^ 2 + det2 a b c d ^ 2 := by
  unfold norm2 dot2 det2
  ring

/-- The three planar vectors satisfy the exact Plücker syzygy in the first
coordinate. -/
theorem det2_syzygy_fst
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    det2 a₂ b₂ a₃ b₃ * a₁ -
        det2 a₁ b₁ a₃ b₃ * a₂ +
        det2 a₁ b₁ a₂ b₂ * a₃ = 0 := by
  unfold det2
  ring

/-- The three planar vectors satisfy the exact Plücker syzygy in the second
coordinate. -/
theorem det2_syzygy_snd
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    det2 a₂ b₂ a₃ b₃ * b₁ -
        det2 a₁ b₁ a₃ b₃ * b₂ +
        det2 a₁ b₁ a₂ b₂ * b₃ = 0 := by
  unfold det2
  ring

/-- Vanishing of the `3 x 3` Gram determinant for three vectors in a
two-dimensional space. -/
theorem gram_rank_two
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    norm2 a₁ b₁ * norm2 a₂ b₂ * norm2 a₃ b₃ +
        2 * dot2 a₁ b₁ a₂ b₂ * dot2 a₂ b₂ a₃ b₃ *
          dot2 a₃ b₃ a₁ b₁ =
      norm2 a₁ b₁ * dot2 a₂ b₂ a₃ b₃ ^ 2 +
        norm2 a₂ b₂ * dot2 a₃ b₃ a₁ b₁ ^ 2 +
        norm2 a₃ b₃ * dot2 a₁ b₁ a₂ b₂ ^ 2 := by
  unfold norm2 dot2
  ring

/-- Exact defect of the product of the three norms from the product of the
three pairwise dot products.  In a narrow sector the determinants are the
small residual chord coordinates. -/
theorem norm_product_sub_dot_product
    (a₁ b₁ a₂ b₂ a₃ b₃ : ℤ) :
    2 *
        (norm2 a₁ b₁ * norm2 a₂ b₂ * norm2 a₃ b₃ -
          dot2 a₁ b₁ a₂ b₂ * dot2 a₂ b₂ a₃ b₃ *
            dot2 a₃ b₃ a₁ b₁) =
      norm2 a₁ b₁ * det2 a₂ b₂ a₃ b₃ ^ 2 +
        norm2 a₂ b₂ * det2 a₃ b₃ a₁ b₁ ^ 2 +
        norm2 a₃ b₃ * det2 a₁ b₁ a₂ b₂ ^ 2 := by
  unfold norm2 dot2 det2
  ring

end GaussianTripleBlocks
end GaussianChain
