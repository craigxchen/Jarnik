import Mathlib

namespace GaussianChain
namespace ChordRatio

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
  rw [map_div, map_sub, map_one, map_sub, map_one, hx, hy]
  field_simp [hx0, hy0, hx1, hy1]
  ring

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
    rw [map_div, map_sub, map_one, map_sub, map_one, hx, hy]
    exact div_ne_zero (sub_ne_zero.mpr (inv_ne_one.mpr hx1))
      (sub_ne_zero.mpr (inv_ne_one.mpr hy1))
  field_simp [hy0, hr] at h ⊢
  nlinarith

end ChordRatio
end GaussianChain
