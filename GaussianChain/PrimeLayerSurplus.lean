import Mathlib

namespace GaussianChain
namespace PrimeLayerSurplus

/--
For two exponents in the interval `[0,e]`, the common valuation at a split
prime plus the common valuation at the conjugate prime is exactly the conductor
exponent minus their valuation distance.

This is the scalar identity behind the conductor-prime contribution to the
Ramana determinant.
-/
theorem min_add_min_complement_add_abs_sub
    (e a b : ℤ)
    (ha0 : 0 ≤ a) (hae : a ≤ e)
    (hb0 : 0 ≤ b) (hbe : b ≤ e) :
    min a b + min (e - a) (e - b) + |a - b| = e := by
  rcases le_total a b with hab | hba
  · rw [min_eq_left hab]
    have hcomp : e - b ≤ e - a := by linarith
    rw [min_eq_right hcomp, abs_of_nonpos (by linarith : a - b ≤ 0)]
    ring
  · rw [min_eq_right hba]
    have hcomp : e - a ≤ e - b := by linarith
    rw [min_eq_left hcomp, abs_of_nonneg (by linarith : 0 ≤ a - b)]
    ring

/--
The deficit of a cut of an odd set from the maximal number `s(s+1)` of
separated pairs is a product of two consecutive integers.
-/
theorem odd_cut_defect_identity (s k : ℤ) :
    s * (s + 1) - k * (2 * s + 1 - k) =
      (k - s) * (k - s - 1) := by
  ring

/--
Consequently every cut of a set of cardinality `2s+1` separates at most
`s(s+1)` unordered pairs.  Equality is possible only at the two balanced cut
sizes `s` and `s+1`.
-/
theorem odd_cut_defect_nonneg
    (s k : ℤ)
    (hs : 0 ≤ s) (hk0 : 0 ≤ k) (hktop : k ≤ 2 * s + 1) :
    0 ≤ s * (s + 1) - k * (2 * s + 1 - k) := by
  rw [odd_cut_defect_identity]
  rcases (by omega : k ≤ s ∨ s + 1 ≤ k) with hlow | hhigh
  · exact mul_nonneg_of_nonpos_of_nonpos (by linarith) (by linarith)
  · exact mul_nonneg (by linarith) (by linarith)

/-- The odd cut defect vanishes at the lower balanced size. -/
theorem odd_cut_defect_eq_zero_at_lower (s : ℤ) :
    s * (s + 1) - s * (2 * s + 1 - s) = 0 := by
  ring

/-- The odd cut defect vanishes at the upper balanced size. -/
theorem odd_cut_defect_eq_zero_at_upper (s : ℤ) :
    s * (s + 1) - (s + 1) * (2 * s + 1 - (s + 1)) = 0 := by
  ring

/--
Away from the two balanced sizes, the defect is at least two.  This integral
gap is useful when a prime layer is known not to be perfectly balanced.
-/
theorem two_le_odd_cut_defect_of_unbalanced
    (s k : ℤ)
    (hs : 0 ≤ s) (hk0 : 0 ≤ k) (hktop : k ≤ 2 * s + 1)
    (hks : k ≠ s) (hks1 : k ≠ s + 1) :
    2 ≤ s * (s + 1) - k * (2 * s + 1 - k) := by
  rw [odd_cut_defect_identity]
  rcases (by omega : k ≤ s - 1 ∨ s + 2 ≤ k) with hlow | hhigh
  · have h1 : k - s ≤ -1 := by linarith
    have h2 : k - s - 1 ≤ -2 := by linarith
    nlinarith [mul_nonneg (show 0 ≤ -(k - s) by linarith)
      (show 0 ≤ -(k - s - 1) by linarith)]
  · have h1 : 2 ≤ k - s := by linarith
    have h2 : 1 ≤ k - s - 1 := by linarith
    nlinarith [mul_nonneg (show 0 ≤ k - s by linarith)
      (show 0 ≤ k - s - 1 by linarith)]

end PrimeLayerSurplus
end GaussianChain
