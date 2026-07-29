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
A compact explicit finite level:
`c = 52 / 15`, `N = 15`, `d = 52`, and `k = d - N = 37`.
The section-space dimension is only `1311`.
-/
theorem sectionCount_15_37 : sectionCount 15 37 = 1311 := by
  norm_num [sectionCount]

/-- Exact finite numerator at the compact level. -/
theorem betaNumerator_15_37 : betaNumerator 15 37 = 19684 := by
  rw [betaNumerator_closed]
  norm_num

/-- Exact finite beta value at the compact level. -/
theorem finiteBeta_15_37 : finiteBeta 15 37 = 1036 / 1035 := by
  rw [finiteBeta_formula]
  norm_num [sectionCount]

/-- The compact finite approximation already lies strictly above `1`. -/
theorem one_lt_finiteBeta_15_37 : 1 < finiteBeta 15 37 := by
  rw [finiteBeta_15_37]
  norm_num

/--
The zero-error coefficient `c - 3β` at the compact level is exactly `32/69`.
-/
theorem coefficient_15_37 :
    (52 / 15 : ℚ) - 3 * finiteBeta 15 37 = 32 / 69 := by
  rw [finiteBeta_15_37]
  norm_num

/-- The compact finite-level coefficient has a strict margin below `1/2`. -/
theorem coefficient_15_37_lt_half :
    (52 / 15 : ℚ) - 3 * finiteBeta 15 37 < 1 / 2 := by
  rw [coefficient_15_37]
  norm_num

/--
The earlier, larger explicit level used in the first version of the route:
`c = 347 / 100`, `N = 100`, `d = 347`, and `k = d - N = 247`.
-/
theorem sectionCount_100_247 : sectionCount 100 247 = 55676 := by
  norm_num [sectionCount]

/-- Exact numerator at the earlier finite level. -/
theorem betaNumerator_100_247 : betaNumerator 100 247 = 5604924 := by
  rw [betaNumerator_closed]
  norm_num

/-- Exact finite beta value at the earlier finite level. -/
theorem finiteBeta_100_247 : finiteBeta 100 247 = 45201 / 44900 := by
  rw [finiteBeta_formula]
  norm_num [sectionCount]

/-- In particular, the earlier finite approximation also lies strictly above `1`. -/
theorem one_lt_finiteBeta_100_247 : 1 < finiteBeta 100 247 := by
  rw [finiteBeta_100_247]
  norm_num

end FiniteBeta
end GaussianChain
