import Mathlib

namespace GaussianChain.ReciprocalCollision

/-- The exact cubic identity behind the signed reciprocal collision.
No endpoint or circle hypothesis is used here. -/
theorem reciprocal_cubic_identity {K : Type*} [Field K]
    (a b c d : K) (ha : a ≠ 0) (hb : b ≠ 0)
    (hc : c ≠ 0) (hd : d ≠ 0) :
    (a + b) * (a - c) * (a - d) =
      a ^ 2 * (a + b - c - d) +
        (a * b * c * d) * (a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹) := by
  field_simp
  ring

/-- Equal sums and reciprocal sums force the two pairs to coincide,
unless both pair sums vanish. -/
theorem pair_collision {K : Type*} [Field K]
    (a b c d : K) (ha : a ≠ 0) (hb : b ≠ 0)
    (hc : c ≠ 0) (hd : d ≠ 0)
    (hsum : a + b = c + d)
    (hinv : a⁻¹ + b⁻¹ = c⁻¹ + d⁻¹) :
    a + b = 0 ∨ (a = c ∧ b = d) ∨ (a = d ∧ b = c) := by
  have hu : a + b - c - d = 0 := by linear_combination hsum
  have hv : a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹ = 0 := by linear_combination hinv
  have hp : (a + b) * (a - c) * (a - d) = 0 := by
    rw [reciprocal_cubic_identity a b c d ha hb hc hd, hu, hv]
    ring
  rcases mul_eq_zero.mp hp with h | h
  · rcases mul_eq_zero.mp h with h | h
    · exact Or.inl h
    · have hac : a = c := sub_eq_zero.mp h
      have hbd : b = d := by linear_combination hsum - hac
      exact Or.inr (Or.inl ⟨hac, hbd⟩)
  · have had : a = d := sub_eq_zero.mp h
    have hbc : b = c := by linear_combination hsum - had
    exact Or.inr (Or.inr ⟨had, hbc⟩)

/-- Signs do not avoid a collision of absolute values. This is the
algebraic input to the four-coordinate row certificates. -/
theorem abs_collision
    (a b c d : ℝ) (ha : a ≠ 0) (hb : b ≠ 0)
    (hc : c ≠ 0) (hd : d ≠ 0)
    (hsum : a + b = c + d)
    (hinv : a⁻¹ + b⁻¹ = c⁻¹ + d⁻¹) :
    |a| = |b| ∨ |a| = |c| ∨ |a| = |d| := by
  rcases pair_collision a b c d ha hb hc hd hsum hinv with h | h | h
  · left
    have hab : a = -b := by linarith
    simp [hab]
  · exact Or.inr (Or.inl (congrArg abs h.1))
  · exact Or.inr (Or.inr (congrArg abs h.1))

/-- A quantitative version with explicitly exposed separation and height
hypotheses. The circle-to-residual reduction is not asserted by this lemma. -/
theorem stable_reciprocal_bound
    (a b c d δ H : ℝ) (ha : a ≠ 0) (hb : b ≠ 0)
    (hc : c ≠ 0) (hd : d ≠ 0)
    (hδ : 0 ≤ δ) (hH : 0 ≤ H)
    (hab : δ ≤ |a + b|) (hac : δ ≤ |a - c|) (had : δ ≤ |a - d|)
    (h_a : |a| ≤ H) (h_b : |b| ≤ H)
    (h_c : |c| ≤ H) (h_d : |d| ≤ H) :
    δ ^ 3 ≤ H ^ 2 * |a + b - c - d| +
      H ^ 4 * |a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹| := by
  have hprod : |a * b * c * d| ≤ H ^ 4 := by
    calc
      |a * b * c * d| = |a| * |b| * |c| * |d| := by simp only [abs_mul]
      _ ≤ H * H * H * H := by gcongr
      _ = H ^ 4 := by ring
  have hlo : δ ^ 3 ≤ |(a + b) * (a - c) * (a - d)| := by
    calc
      δ ^ 3 = δ * δ * δ := by ring
      _ ≤ |a + b| * |a - c| * |a - d| := by gcongr
      _ = |(a + b) * (a - c) * (a - d)| := by simp only [abs_mul]
  calc
    δ ^ 3 ≤ |(a + b) * (a - c) * (a - d)| := hlo
    _ = |a ^ 2 * (a + b - c - d) +
          (a * b * c * d) * (a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹)| :=
      congrArg abs (reciprocal_cubic_identity a b c d ha hb hc hd)
    _ ≤ |a ^ 2 * (a + b - c - d)| +
          |(a * b * c * d) * (a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹)| := abs_add_le _ _
    _ = |a| ^ 2 * |a + b - c - d| +
          |a * b * c * d| * |a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹| := by
      rw [abs_mul, abs_pow, abs_mul]
    _ ≤ H ^ 2 * |a + b - c - d| +
          H ^ 4 * |a⁻¹ + b⁻¹ - c⁻¹ - d⁻¹| := by gcongr

end GaussianChain.ReciprocalCollision
