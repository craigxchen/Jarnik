import Mathlib

namespace GaussianChain
namespace ChordRatio

/--
Cross-multiplied chord-ratio identity.

For a norm-one involution, the conjugate of `(x - 1) / (y - 1)` is
`(x⁻¹ - 1) / (y⁻¹ - 1)`.  We take that transformed value as an explicit
argument here, keeping this lemma independent of the particular conjugation API.
-/
theorem chordRatio_cross_transformed
    {x y rbar : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hrbar : rbar = (x⁻¹ - 1) / (y⁻¹ - 1)) :
    y * ((x - 1) / (y - 1)) = x * rbar := by
  have hxi : x⁻¹ - 1 ≠ 0 := sub_ne_zero.mpr (inv_ne_one.mpr hx1)
  have hyi : y⁻¹ - 1 ≠ 0 := sub_ne_zero.mpr (inv_ne_one.mpr hy1)
  rw [hrbar]
  field_simp [hx0, hy0, hx1, hy1, hxi, hyi]
  ring

/-- Purely algebraic quotient form of `chordRatio_cross_transformed`. -/
theorem chordRatio_div_transformed
    {x y rbar : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hrbar : rbar = (x⁻¹ - 1) / (y⁻¹ - 1)) :
    ((x - 1) / (y - 1)) / rbar = x / y := by
  have hr : rbar ≠ 0 := by
    rw [hrbar]
    exact div_ne_zero (sub_ne_zero.mpr (inv_ne_one.mpr hx1))
      (sub_ne_zero.mpr (inv_ne_one.mpr hy1))
  apply (div_eq_div_iff hr hy0).2
  simpa [mul_comm] using chordRatio_cross_transformed hx0 hy0 hx1 hy1 hrbar

end ChordRatio
end GaussianChain
