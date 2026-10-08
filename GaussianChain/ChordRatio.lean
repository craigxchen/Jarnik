import Mathlib

namespace GaussianChain
namespace ChordRatio

open scoped ComplexConjugate

/--
Algebraic core of the chord-ratio identity.

For norm-one elements the involution is inversion.  Hence the ratio of two
chords has norm-one quotient equal to `x / y`:

`r = (x - 1) / (y - 1)` implies `r / conj(r) = x / y`.

The theorem is stated with the inversion identities as hypotheses so that the
proof is purely algebraic and can later be reused for number fields with an
involution.
-/
theorem chordRatio_div_conj
    {x y : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hx : conj x = x⁻¹) (hy : conj y = y⁻¹) :
    ((x - 1) / (y - 1)) / conj ((x - 1) / (y - 1)) = x / y := by
  have hx' : x - 1 ≠ 0 := sub_ne_zero.mpr hx1
  have hy' : y - 1 ≠ 0 := sub_ne_zero.mpr hy1
  have hconj : conj ((x - 1) / (y - 1)) =
      y * (x - 1) / (x * (y - 1)) := by
    rw [map_div₀, map_sub, map_one, map_sub, map_one, hx, hy]
    apply (div_eq_div_iff
      (sub_ne_zero.mpr (inv_ne_one.mpr hy1)) (mul_ne_zero hx0 hy')).mpr
    field_simp [hx0, hy0]
    ring
  rw [hconj]
  field_simp [hx0, hy0, hx', hy']

/--
Equivalent cross-multiplied form, useful when avoiding division by the chord
ratio itself.
-/
theorem chordRatio_cross
    {x y : ℂ}
    (hx0 : x ≠ 0) (hy0 : y ≠ 0)
    (hx1 : x ≠ 1) (hy1 : y ≠ 1)
    (hx : conj x = x⁻¹) (hy : conj y = y⁻¹) :
    y * ((x - 1) / (y - 1)) =
      x * conj ((x - 1) / (y - 1)) := by
  have h := chordRatio_div_conj hx0 hy0 hx1 hy1 hx hy
  have hr : conj ((x - 1) / (y - 1)) ≠ 0 := by
    rw [map_div₀, map_sub, map_one, map_sub, map_one, hx, hy]
    exact div_ne_zero (sub_ne_zero.mpr (inv_ne_one.mpr hx1))
      (sub_ne_zero.mpr (inv_ne_one.mpr hy1))
  simpa only [mul_comm] using (div_eq_div_iff hr hy0).mp h

end ChordRatio
end GaussianChain
