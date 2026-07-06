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

open MissingPrimeSubcritical
open SubcriticalBound

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

/-- Many missing primes in a fixed interval give the missing-prime determinant branch.

This packages `card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound` with
`P` specialized to the actual set of interval primes not dividing `N`. -/
theorem missing_prime_branch
    {M s N m U : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hsmallPrime : 4 * U ≤ s)
    (hN : 0 < N)
    (hweightedLog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
            2 * missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  classical
  let P := MertensLower.missingPrimeIntervalFinset N m U
  refine card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound
    (M := M) (s := s) (N := N) (a := a) (L := L) (B := B)
    (z := z) (t := t) (P := P)
    hk hB hfloor ?hprime ?hsmallPrime hN ?hpN ?hweightedLog hcircle hz hmono hparam hmem
  · intro p hp
    exact MertensLower.prime_of_mem_missingPrimeIntervalFinset (N := N) (m := m) (U := U) hp
  · intro p hp
    have hpU : p ≤ U :=
      MertensLower.le_upper_of_mem_missingPrimeIntervalFinset (N := N) (m := m) (U := U) hp
    exact (Nat.mul_le_mul_left 4 hpU).trans hsmallPrime
  · intro p hp
    exact MertensLower.int_not_dvd_of_mem_missingPrimeIntervalFinset
      (N := N) (m := m) (U := U) hp
  · simpa [P, MertensLower.weightedMissingPrimeInterval, missingPrimeLogSum] using hweightedLog

end MissingPrimeIntervalBranch
end GaussianChain
