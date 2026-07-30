import Mathlib

namespace GaussianChain
namespace MaxComplementIdentity

/--
For three additive valuations, the maximum is the total valuation minus the
minimum valuation of the three complementary products.

If `a`, `b`, and `c` are the valuations of three fixed product sections, then
`b+c`, `a+c`, and `a+b` are the valuations of the generators obtained by
omitting one product.  This identity packages the local pointwise maximum as a
fixed divisor minus the local height of one fixed ideal.
-/
theorem max_three_eq_sum_sub_min_complements (a b c : ℝ) :
    max a (max b c) =
      a + b + c - min (b + c) (min (a + c) (a + b)) := by
  rcases le_total a b with hab | hba
  · rcases le_total b c with hbc | hcb
    · have hac : a ≤ c := hab.trans hbc
      simp [max_eq_right hab, max_eq_right hbc, min_eq_left hbc,
        min_eq_left hac, min_eq_left hab]
      ring
    · rcases le_total a c with hac | hca
      · simp [max_eq_right hab, max_eq_left hcb, max_eq_right hac,
          min_eq_right hcb, min_eq_left hac, min_eq_left hab]
        ring
      · simp [max_eq_right hab, max_eq_left hcb, max_eq_left hca,
          min_eq_right hcb, min_eq_right hca, min_eq_left hab]
        ring
  · rcases le_total a c with hac | hca
    · rcases le_total b c with hbc | hcb
      · simp [max_eq_left hba, max_eq_right hbc, max_eq_right hac,
          min_eq_left hbc, min_eq_left hac, min_eq_right hba]
        ring
      · simp [max_eq_left hba, max_eq_left hcb, min_eq_right hcb,
          min_eq_left hac, min_eq_right hba]
        ring
    · have hcb : c ≤ b := hca.trans hba
      simp [max_eq_left hba, max_eq_left hcb, min_eq_right hcb,
        min_eq_right hca, min_eq_right hba]
      ring

/-- Multiplicative form of the same identity for positive real absolute values. -/
theorem min_three_eq_product_div_max_complements
    {x y z : ℝ} (hx : 0 < x) (hy : 0 < y) (hz : 0 < z) :
    min x (min y z) =
      x * y * z / max (y * z) (max (x * z) (x * y)) := by
  rcases le_total x y with hxy | hyx
  · rcases le_total y z with hyz | hzy
    · have hxz : x ≤ z := hxy.trans hyz
      have hyzmul : y * z ≤ x * z ∨ x * z ≤ y * z := le_total _ _
      simp [min_eq_left hxy, min_eq_left hyz, max_eq_right
        (mul_le_mul_of_nonneg_right hxy hz.le),
        max_eq_right (mul_le_mul_of_nonneg_left hyz hx.le)]
      field_simp
      ring
    · rcases le_total x z with hxz | hzx
      · simp [min_eq_left hxy, min_eq_right hzy, min_eq_left hxz]
        have h1 : x * y ≤ y * z := mul_le_mul_of_nonneg_left hxz hy.le
        have h2 : x * z ≤ y * z := mul_le_mul_of_nonneg_right hxy hz.le
        rw [max_eq_left h2, max_eq_left h1]
        field_simp
        ring
      · simp [min_eq_left hxy, min_eq_right hzy, min_eq_right hzx]
        have h1 : y * z ≤ x * y := mul_le_mul_of_nonneg_right hzx hy.le
        have h2 : x * z ≤ x * y := mul_le_mul_of_nonneg_left hzy hx.le
        rw [max_eq_right h1, max_eq_right h2]
        field_simp
        ring
  · rcases le_total x z with hxz | hzx
    · rcases le_total y z with hyz | hzy
      · simp [min_eq_right hyx, min_eq_left hyz]
        have h1 : x * y ≤ y * z := mul_le_mul_of_nonneg_left hxz hy.le
        have h2 : x * z ≤ y * z := mul_le_mul_of_nonneg_right hyx hz.le
        rw [max_eq_left h2, max_eq_left h1]
        field_simp
        ring
      · simp [min_eq_right hyx, min_eq_right hzy]
        have h1 : x * z ≤ x * y := mul_le_mul_of_nonneg_left hzy hx.le
        have h2 : y * z ≤ x * y := mul_le_mul_of_nonneg_right hyx hy.le
        rw [max_eq_right h1, max_eq_right h2]
        field_simp
        ring
    · have hzy : z ≤ y := hzx.trans hyx
      simp [min_eq_right hyx, min_eq_right hzy]
      have h1 : y * z ≤ x * y := mul_le_mul_of_nonneg_right hzx hy.le
      have h2 : x * z ≤ x * y := mul_le_mul_of_nonneg_left hzy hx.le
      rw [max_eq_right h2, max_eq_right h1]
      field_simp
      ring

end MaxComplementIdentity
end GaussianChain
