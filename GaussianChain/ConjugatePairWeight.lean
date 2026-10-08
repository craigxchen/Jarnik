import Mathlib

namespace GaussianChain
namespace ConjugatePairWeight

/-- Coordinate local heights at one orientation of a split prime. -/
def forwardHeight (d₀ d₁ d₂ : ℝ) (i : Fin 3) : ℝ :=
  max d₀ (max d₁ d₂) - ![d₀, d₁, d₂] i

/-- Coordinate local heights at the conjugate orientation. -/
def backwardHeight (d₀ d₁ d₂ : ℝ) (i : Fin 3) : ℝ :=
  ![d₀, d₁, d₂] i - min d₀ (min d₁ d₂)

/-- The sum of the three forward local coordinate heights. -/
theorem sum_forwardHeight (d₀ d₁ d₂ : ℝ) :
    ∑ i : Fin 3, forwardHeight d₀ d₁ d₂ i =
      3 * max d₀ (max d₁ d₂) - (d₀ + d₁ + d₂) := by
  fin_cases i <;> simp [forwardHeight]
  ring

/-- The sum of the three conjugate local coordinate heights. -/
theorem sum_backwardHeight (d₀ d₁ d₂ : ℝ) :
    ∑ i : Fin 3, backwardHeight d₀ d₁ d₂ i =
      (d₀ + d₁ + d₂) - 3 * min d₀ (min d₁ d₂) := by
  fin_cases i <;> simp [backwardHeight]
  ring

/-- The largest forward height is the valuation range. -/
theorem max_forwardHeight (d₀ d₁ d₂ : ℝ) :
    max (forwardHeight d₀ d₁ d₂ 0)
        (max (forwardHeight d₀ d₁ d₂ 1)
             (forwardHeight d₀ d₁ d₂ 2)) =
      max d₀ (max d₁ d₂) - min d₀ (min d₁ d₂) := by
  simp [forwardHeight]
  rcases le_total d₀ d₁ with h01 | h10 <;>
  rcases le_total d₀ d₂ with h02 | h20 <;>
  rcases le_total d₁ d₂ with h12 | h21 <;>
  simp [max_eq_left, max_eq_right, min_eq_left, min_eq_right, *] <;>
  linarith

/-- The largest conjugate height is the same valuation range. -/
theorem max_backwardHeight (d₀ d₁ d₂ : ℝ) :
    max (backwardHeight d₀ d₁ d₂ 0)
        (max (backwardHeight d₀ d₁ d₂ 1)
             (backwardHeight d₀ d₁ d₂ 2)) =
      max d₀ (max d₁ d₂) - min d₀ (min d₁ d₂) := by
  simp [backwardHeight]
  rcases le_total d₀ d₁ with h01 | h10 <;>
  rcases le_total d₀ d₂ with h02 | h20 <;>
  rcases le_total d₁ d₂ with h12 | h21 <;>
  simp [max_eq_left, max_eq_right, min_eq_left, min_eq_right, *] <;>
  linarith

/--
Exact conductor-specific improvement.

At each orientation the explicit sextic star basis has weight

`56 * (sum coordinate heights) - 6 * (largest coordinate height)`.

Pairing a split prime with its conjugate gives exactly `156` times the
valuation range.
-/
theorem paired_sextic_weight_eq
    (d₀ d₁ d₂ : ℝ) :
    (56 * (∑ i : Fin 3, forwardHeight d₀ d₁ d₂ i) -
        6 * max (forwardHeight d₀ d₁ d₂ 0)
          (max (forwardHeight d₀ d₁ d₂ 1)
               (forwardHeight d₀ d₁ d₂ 2))) +
      (56 * (∑ i : Fin 3, backwardHeight d₀ d₁ d₂ i) -
        6 * max (backwardHeight d₀ d₁ d₂ 0)
          (max (backwardHeight d₀ d₁ d₂ 1)
               (backwardHeight d₀ d₁ d₂ 2))) =
      156 * (max d₀ (max d₁ d₂) - min d₀ (min d₁ d₂)) := by
  rw [sum_forwardHeight, sum_backwardHeight,
    max_forwardHeight, max_backwardHeight]
  ring

end ConjugatePairWeight
end GaussianChain
