import GaussianChain.Basic

namespace GaussianChain

/-- The rank-two Plucker relation for four planar integer vectors.

For cyclically ordered vectors the determinants are positive, and this is the
algebraic form of Ptolemy after the vertex squareclass factors have been
removed. -/
theorem det2_plucker
    (a₁ b₁ a₂ b₂ a₃ b₃ a₄ b₄ : ℤ) :
    det2 a₁ b₁ a₃ b₃ * det2 a₂ b₂ a₄ b₄ =
      det2 a₁ b₁ a₂ b₂ * det2 a₃ b₃ a₄ b₄ +
      det2 a₁ b₁ a₄ b₄ * det2 a₂ b₂ a₃ b₃ := by
  unfold det2
  ring

/-- Rearranged form useful for solving for the crossing Plucker coordinate. -/
theorem det2_plucker_sub
    (a₁ b₁ a₂ b₂ a₃ b₃ a₄ b₄ : ℤ) :
    det2 a₁ b₁ a₃ b₃ * det2 a₂ b₂ a₄ b₄ -
        det2 a₁ b₁ a₄ b₄ * det2 a₂ b₂ a₃ b₃ =
      det2 a₁ b₁ a₂ b₂ * det2 a₃ b₃ a₄ b₄ := by
  rw [det2_plucker]
  ring

end GaussianChain
