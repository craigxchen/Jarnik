import GaussianChain.Basic

namespace GaussianChain

/-- Dot product of two integer vectors. -/
def dot2 (x₁ y₁ x₂ y₂ : ℤ) : ℤ := x₁ * x₂ + y₁ * y₂

/-- Oriented area of two integer vectors. -/
def cross2 (x₁ y₁ x₂ y₂ : ℤ) : ℤ := x₁ * y₂ - y₁ * x₂

/-- The dot-cross map with a fixed vector is a similarity: squared lengths are
multiplied by the squared length of the fixed vector. -/
theorem dot_cross_norm_identity
    (x y u v : ℤ) :
    cross2 x y u v ^ 2 + dot2 x y u v ^ 2 =
      norm2 x y * norm2 u v := by
  unfold cross2 dot2 norm2
  ring

/-- The dot-cross map also scales pairwise squared distances by the same factor. -/
theorem dot_cross_sqDist_identity
    (x y u₁ v₁ u₂ v₂ : ℤ) :
    (cross2 x y u₂ v₂ - cross2 x y u₁ v₁) ^ 2 +
        (dot2 x y u₂ v₂ - dot2 x y u₁ v₁) ^ 2 =
      norm2 x y * norm2 (u₂ - u₁) (v₂ - v₁) := by
  unfold cross2 dot2 norm2
  ring

/-- Algebraic relation making the certificate extension's second reduced
coordinate a dot product with the fixed third chord. -/
theorem reducedDelta_eq_dot
    (A B h r m n dm dn : ℤ)
    (hm : m = -2 * h * A + r * n)
    (hn : n = -2 * h * B - r * m) :
    2 * h * (-(A * dm + B * dn)) -
        r * cross2 m n dm dn =
      dot2 m n dm dn := by
  rw [hm, hn]
  unfold cross2 dot2
  ring

/-- Consequently the affine extension map, after subtracting the image of the
base point, is exactly the dot-cross similarity determined by the third chord. -/
theorem extensionDifference_is_similarity
    (A B h r m n dm dn : ℤ)
    (hm : m = -2 * h * A + r * n)
    (hn : n = -2 * h * B - r * m) :
    cross2 m n dm dn ^ 2 +
        (2 * h * (-(A * dm + B * dn)) -
          r * cross2 m n dm dn) ^ 2 =
      norm2 m n * norm2 dm dn := by
  rw [reducedDelta_eq_dot A B h r m n dm dn hm hn]
  exact dot_cross_norm_identity m n dm dn

end GaussianChain
