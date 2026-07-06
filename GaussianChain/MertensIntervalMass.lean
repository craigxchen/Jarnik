/-
An explicit lower bound for the weighted prime mass of an interval, derived from the vendored
explicit Mertens first theorem (`GaussianChain.MertensFirst`): for `1 ≤ m ≤ n`,

  `log n - log m - mertensDeficit ≤ weightedPrimeInterval m n`,

where `mertensDeficit = (log 4 + 4) + (2 + Mertens.E₁)` is at most `10`.
-/
import GaussianChain.MertensFirst
import GaussianChain.MertensLower
import Mathlib.Analysis.Complex.ExponentialBounds

namespace GaussianChain
namespace MertensLower

/-- `Nat.primesLE n` is the `Ioc`-filtered prime finset appearing in the vendored Mertens
first theorem. -/
theorem primesLE_eq_filter_Ioc (n : ℕ) :
    Nat.primesLE n = (Finset.Ioc 0 n).filter Nat.Prime :=
  Nat.primesLE_eq_filter_Icc_one n

/-- Bridge to the Mertens normalization: the cumulative weighted prime sum up to `n` equals
`log n` plus the Mertens error term `E₁p n`. -/
theorem weightedPrimeSum_eq_log_add_E₁p (n : ℕ) :
    weightedPrimeSum n = Real.log (n : ℝ) + Mertens.E₁p (n : ℝ) := by
  have h := Mertens.sum_log_prime_div_eq (n : ℝ)
  rw [Nat.floor_natCast] at h
  rw [weightedPrimeSum, primesLE_eq_filter_Ioc]
  exact h

/-- The weighted prime interval sum is the difference of the cumulative weighted prime sums. -/
theorem weightedPrimeInterval_eq_sub {m n : ℕ} (hmn : m ≤ n) :
    weightedPrimeInterval m n = weightedPrimeSum n - weightedPrimeSum m := by
  unfold weightedPrimeInterval weightedPrimeSum
  exact Finset.sum_sdiff_eq_sub (Nat.primesLE_mono hmn)

/-- The explicit deficit constant in the Mertens interval lower bound: the upper Mertens error
`log 4 + 4` at the lower endpoint plus the lower Mertens error `2 + E₁` at the upper endpoint. -/
noncomputable def mertensDeficit : ℝ := (Real.log 4 + 4) + (2 + Mertens.E₁)

/-- Explicit Mertens lower bound for the weighted prime mass of the interval `m < p ≤ n`. -/
theorem weightedPrimeInterval_ge_log_sub_deficit {m n : ℕ} (hm : 1 ≤ m) (hmn : m ≤ n) :
    Real.log (n : ℝ) - Real.log (m : ℝ) - mertensDeficit ≤ weightedPrimeInterval m n := by
  have hm1 : (1 : ℝ) ≤ (m : ℝ) := by exact_mod_cast hm
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hm.trans hmn
  have hupper : Mertens.E₁p (m : ℝ) ≤ Real.log 4 + 4 := Mertens.E₁p.le hm1
  have hlower : Mertens.E₁p (n : ℝ) ≥ -2 - Mertens.E₁ := Mertens.E₁p.ge hn1
  rw [weightedPrimeInterval_eq_sub hmn, weightedPrimeSum_eq_log_add_E₁p,
    weightedPrimeSum_eq_log_add_E₁p]
  unfold mertensDeficit
  linarith

/-- Numeric bound on the deficit constant; its exact value is
`(13 * log 2 + 27) / 4 ≈ 9.003`. -/
theorem mertensDeficit_le_ten : mertensDeficit ≤ 10 := by
  have hE := Mertens.E₁.le
  have hlog2 : Real.log 2 < 0.6931471808 := Real.log_two_lt_d9
  unfold mertensDeficit
  rw [log_four_eq_two_mul_log_two]
  linarith

/-- The Mertens deficit constant is nonnegative. -/
theorem mertensDeficit_nonneg : 0 ≤ mertensDeficit := by
  have hE := Mertens.E₁.nonneg
  have hlog4 : 0 ≤ Real.log (4 : ℝ) := Real.log_nonneg (by norm_num)
  unfold mertensDeficit
  linarith

/-- If the missing-prime mass of `m < p ≤ n` is below the explicit Mertens lower bound
`log n - log m - mertensDeficit`, then some prime in the interval divides `N`. -/
theorem exists_prime_dvd_of_missing_lt_log_sub_deficit {N m n : ℕ}
    (hm : 1 ≤ m) (hmn : m ≤ n)
    (hmissing : weightedMissingPrimeInterval N m n <
      Real.log (n : ℝ) - Real.log (m : ℝ) - mertensDeficit) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n ∧ p ∣ N :=
  exists_prime_dvd_of_missing_lt_interval_lower
    (weightedPrimeInterval_ge_log_sub_deficit hm hmn) hmissing

end MertensLower
end GaussianChain
