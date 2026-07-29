import Mathlib

namespace GaussianChain
namespace FiniteBeta

/--
The dimension of degree `N + j` plane forms vanishing to order at least `N`
at one point, written as its Hilbert-polynomial expression over `ℚ`.

Geometrically this is

`h⁰(Bl_P ℙ², (N + j) H - N E)`.
-/
def sectionCount (N : ℚ) (j : ℕ) : ℚ :=
  ((N + (j : ℚ) + 2) * (N + (j : ℚ) + 1) - (N + 1) * N) / 2

/-- The finite numerator appearing in the Ru--Vojta beta approximation. -/
def betaNumerator (N : ℚ) (k : ℕ) : ℚ :=
  Finset.sum (Finset.range k) (fun j ↦ sectionCount N j)

/--
Exact finite sum of the section dimensions.  Here `k = d - N`, where `d`
is the plane degree and `N` is the required multiplicity at the blown-up
point.
-/
theorem betaNumerator_closed (N : ℚ) (k : ℕ) :
    betaNumerator N k =
      (k : ℚ) * ((k : ℚ) + 1) * (3 * N + (k : ℚ) + 2) / 6 := by
  unfold betaNumerator
  induction k with
  | zero =>
      simp
  | succ k ih =>
      rw [Finset.sum_range_succ, ih]
      simp only [sectionCount, Nat.cast_succ]
      ring

/-- The finite beta ratio at level `N` and excess degree `k`. -/
def finiteBeta (N : ℚ) (k : ℕ) : ℚ :=
  betaNumerator N k / (N * sectionCount N k)

/-- Closed expression for the finite beta ratio. -/
theorem finiteBeta_formula (N : ℚ) (k : ℕ) :
    finiteBeta N k =
      ((k : ℚ) * ((k : ℚ) + 1) * (3 * N + (k : ℚ) + 2) / 6) /
        (N * sectionCount N k) := by
  rw [finiteBeta, betaNumerator_closed]

/-- The numerator factor in the limiting beta expression. -/
theorem betaLimit_numerator_factor (c : ℚ) :
    (c - 1) ^ 2 * (c + 2) = c ^ 3 - 3 * c + 2 := by
  ring

/--
An explicit finite level used in the uniform-GCD route:
`c = 347 / 100`, `N = 100`, `d = 347`, and `k = d - N = 247`.
-/
theorem sectionCount_100_247 : sectionCount 100 247 = 55676 := by
  norm_num [sectionCount]

/-- Exact numerator at the explicit finite level. -/
theorem betaNumerator_100_247 : betaNumerator 100 247 = 5604924 := by
  rw [betaNumerator_closed]
  norm_num

/-- Exact finite beta value at the explicit finite level. -/
theorem finiteBeta_100_247 : finiteBeta 100 247 = 45201 / 44900 := by
  rw [finiteBeta_formula]
  norm_num [sectionCount]

/-- In particular, this finite approximation already lies strictly above `1`. -/
theorem one_lt_finiteBeta_100_247 : 1 < finiteBeta 100 247 := by
  rw [finiteBeta_100_247]
  norm_num

end FiniteBeta
end GaussianChain
