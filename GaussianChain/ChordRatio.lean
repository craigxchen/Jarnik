import Mathlib

namespace GaussianChain
namespace ChordRatio

/--
Algebraic core of the chord-ratio identity.

For norm-one elements complex conjugation is inversion. Hence the ratio of two
chords has norm-one quotient equal to `x / y`.
-/
theorem chordRatio_div_conj
    {x y : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hx : Complex.conj x = x⁻¹) (hy : Complex.conj y = y⁻¹) :
    ((x - 1) / (y - 1)) /
        Complex.conj ((x - 1) / (y - 1)) = x / y := by
  rw [map_div, map_sub, map_one, map_sub, map_one, hx, hy]
  field_simp [hx0, hy0, hx1, hy1]
  ring

/-- Equivalent cross-multiplied form. -/
theorem chordRatio_cross
    {x y : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hx : Complex.conj x = x⁻¹) (hy : Complex.conj y = y⁻¹) :
    y * ((x - 1) / (y - 1)) =
      x * Complex.conj ((x - 1) / (y - 1)) := by
  have h := chordRatio_div_conj hx0 hy0 hx1 hy1 hx hy
  have hr : Complex.conj ((x - 1) / (y - 1)) ≠ 0 := by
    rw [map_div, map_sub, map_one, map_sub, map_one, hx, hy]
    exact div_ne_zero (sub_ne_zero.mpr (inv_ne_one.mpr hx1))
      (sub_ne_zero.mpr (inv_ne_one.mpr hy1))
  have hcross := (div_eq_div_iff hr hy0).mp h
  simpa [mul_comm] using hcross

end ChordRatio
end GaussianChain
