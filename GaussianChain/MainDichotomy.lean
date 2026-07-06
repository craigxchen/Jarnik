import GaussianChain.MissingPrimeIntervalBranch

namespace GaussianChain
namespace MainDichotomy

open SubcriticalBound

/-- The weighted-log condition that makes the missing-prime determinant branch fire on an
interval `m < p ≤ U`. -/
noncomputable def manyMissingWeightedLogCondition (s N m U : ℕ) (B : ℝ) : Prop :=
  ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
    ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
        2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
      ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)

/-- If a lower bound `μ` for the missing-prime mass is already enough to make the
weighted-log determinant inequality true, then the actual missing-prime mass also makes the
many-missing branch true. -/
theorem manyMissingWeightedLogCondition_of_missing_lower_bound
    {s N m U : ℕ} {B μ : ℝ}
    (hmissing : μ ≤ MertensLower.weightedMissingPrimeInterval N m U)
    (hineq :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * μ - 2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    manyMissingWeightedLogCondition s N m U B := by
  have hs_nonneg : 0 ≤ (s : ℝ) ^ 2 := sq_nonneg _
  have hgain :
      ((s : ℝ) ^ 2 * μ - 2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
        ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    nlinarith [mul_le_mul_of_nonneg_left hmissing hs_nonneg]
  exact hineq.trans_le hgain

end MainDichotomy
end GaussianChain
