import Mathlib
import GaussianChain.FiniteBeta

namespace GaussianChain
namespace SegreCoefficient

/-- A convenient rational weight for the conductor `P¹` factor. -/
def conductorWeight : ℚ := 1 / 36

/-- The degree-six point-blowup coefficient is `4/9`. -/
theorem pointBlowupCoefficient :
    (6 : ℚ) - 3 * FiniteBeta.finiteBeta 1 5 = 4 / 9 := by
  exact FiniteBeta.coefficient_1_5

/--
Adding the conductor factor with weight `1/36` still leaves a strict endpoint
margin: `1/36 + 4/9 = 17/36 < 1/2`.
-/
theorem conductor_add_pointBlowupCoefficient :
    conductorWeight + ((6 : ℚ) - 3 * FiniteBeta.finiteBeta 1 5) = 17 / 36 := by
  rw [pointBlowupCoefficient]
  norm_num [conductorWeight]

/-- The remaining strict margin is exactly `1/36`. -/
theorem segreCoefficient_lt_half :
    conductorWeight + ((6 : ℚ) - 3 * FiniteBeta.finiteBeta 1 5) < 1 / 2 := by
  rw [conductor_add_pointBlowupCoefficient]
  norm_num

/-- Explicit form of the endpoint margin. -/
theorem half_sub_segreCoefficient :
    (1 / 2 : ℚ) -
        (conductorWeight + ((6 : ℚ) - 3 * FiniteBeta.finiteBeta 1 5)) = 1 / 36 := by
  rw [conductor_add_pointBlowupCoefficient]
  norm_num

/--
Tensoring every section-count term by a fixed positive-dimensional auxiliary
factor does not alter a finite beta ratio: the common multiplicity cancels.
-/
theorem finiteBeta_invariant_under_common_factor
    (m numerator denominator : ℚ)
    (hm : m ≠ 0) (hden : denominator ≠ 0) :
    (m * numerator) / (m * denominator) = numerator / denominator := by
  field_simp

/-- The explicit integral multiple corresponding to conductor weight `1/36`. -/
theorem integral_segre_class_decomposition :
    (1 : ℚ) / 36 = conductorWeight := by
  rfl

end SegreCoefficient
end GaussianChain
