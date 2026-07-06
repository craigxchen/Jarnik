import GaussianChain.MissingPrimeSubcritical
import GaussianChain.MertensLower

namespace GaussianChain

namespace MertensLower

/-- A member of a missing-prime interval is a prime. -/
theorem prime_of_mem_missingPrimeIntervalFinset {N m U p : ℕ}
    (hp : p ∈ missingPrimeIntervalFinset N m U) :
    Nat.Prime p := by
  classical
  rw [missingPrimeIntervalFinset, Finset.mem_filter, Finset.mem_sdiff] at hp
  exact Nat.prime_of_mem_primesLE hp.1.1

/-- A member of a missing-prime interval is at most the interval upper endpoint. -/
theorem le_upper_of_mem_missingPrimeIntervalFinset {N m U p : ℕ}
    (hp : p ∈ missingPrimeIntervalFinset N m U) :
    p ≤ U :=
  (bounds_of_mem_missingPrimeIntervalFinset (N := N) hp).2

/-- A member of a missing-prime interval does not divide the common norm, in integer form. -/
theorem int_not_dvd_of_mem_missingPrimeIntervalFinset {N m U p : ℕ}
    (hp : p ∈ missingPrimeIntervalFinset N m U) :
    ¬ (p : ℤ) ∣ (N : ℤ) := by
  classical
  rw [missingPrimeIntervalFinset, Finset.mem_filter] at hp
  exact fun hInt => hp.2 (by exact_mod_cast hInt)

end MertensLower

namespace MissingPrimeIntervalBranch

/-- The unweighted logarithmic mass of the missing primes in an interval. -/
noncomputable def missingPrimeLogSum (N m U : ℕ) : ℝ :=
  ∑ p ∈ MertensLower.missingPrimeIntervalFinset N m U, Real.log (p : ℝ)

/-- Crude upper bound for the unweighted missing-prime logarithmic mass. -/
theorem missingPrimeLogSum_le_succ_mul_log
    {N m U : ℕ} (hU : 1 ≤ U) :
    missingPrimeLogSum N m U ≤ ((U + 1 : ℕ) : ℝ) * Real.log (U : ℝ) := by
  classical
  let S := MertensLower.missingPrimeIntervalFinset N m U
  have hlogU_nonneg : 0 ≤ Real.log (U : ℝ) :=
    Real.log_nonneg (by exact_mod_cast hU)
  have hterm : ∀ p ∈ S, Real.log (p : ℝ) ≤ Real.log (U : ℝ) := by
    intro p hp
    have hp_prime := MertensLower.prime_of_mem_missingPrimeIntervalFinset
      (N := N) (m := m) (U := U) hp
    have hpU := MertensLower.le_upper_of_mem_missingPrimeIntervalFinset
      (N := N) (m := m) (U := U) hp
    exact Real.log_le_log (by exact_mod_cast hp_prime.pos) (by exact_mod_cast hpU)
  have hsum :
      missingPrimeLogSum N m U ≤ (S.card : ℝ) * Real.log (U : ℝ) := by
    calc
      missingPrimeLogSum N m U = ∑ p ∈ S, Real.log (p : ℝ) := by
        simp [missingPrimeLogSum, S]
      _ ≤ ∑ _p ∈ S, Real.log (U : ℝ) := by
        exact Finset.sum_le_sum hterm
      _ = (S.card : ℝ) * Real.log (U : ℝ) := by
        simp [nsmul_eq_mul]
  have hsub : S ⊆ Finset.range (U + 1) := by
    intro p hp
    have hpU := MertensLower.le_upper_of_mem_missingPrimeIntervalFinset
      (N := N) (m := m) (U := U) hp
    exact Finset.mem_range.mpr (Nat.lt_succ_of_le hpU)
  have hcard_nat : S.card ≤ U + 1 := by
    simpa using Finset.card_le_card hsub
  have hcard_real : (S.card : ℝ) ≤ ((U + 1 : ℕ) : ℝ) := by
    exact_mod_cast hcard_nat
  exact hsum.trans (mul_le_mul_of_nonneg_right hcard_real hlogU_nonneg)

end MissingPrimeIntervalBranch
end GaussianChain
