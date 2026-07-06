import Mathlib.NumberTheory.Chebyshev
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics

namespace GaussianChain
namespace MertensLower

open Finset
open Filter Asymptotics
open scoped Nat.Prime

/-- The weighted prime reciprocal sum `∑_{p ≤ n} log p / p`. -/
noncomputable def weightedPrimeSum (n : ℕ) : ℝ :=
  ∑ p ∈ Nat.primesLE n, Real.log (p : ℝ) / (p : ℝ)

/-- The weighted prime reciprocal sum over the interval `m < p ≤ n`. -/
noncomputable def weightedPrimeInterval (m n : ℕ) : ℝ :=
  ∑ p ∈ Nat.primesLE n \ Nat.primesLE m, Real.log (p : ℝ) / (p : ℝ)

/-- The finite set of interval primes `m < p ≤ n` which do not divide `N`. -/
noncomputable def missingPrimeIntervalFinset (N m n : ℕ) : Finset ℕ :=
  (Nat.primesLE n \ Nat.primesLE m).filter fun p => ¬ p ∣ N

/-- The part of `weightedPrimeInterval m n` contributed by interval primes which do not divide
`N`. These are the "missing" primes in the determinant/descent dichotomy. -/
noncomputable def weightedMissingPrimeInterval (N m n : ℕ) : ℝ :=
  ∑ p ∈ missingPrimeIntervalFinset N m n, Real.log (p : ℝ) / (p : ℝ)

/-- If the weighted missing-prime contribution is strictly smaller than the full interval
contribution, then some interval prime divides `N`. -/
theorem exists_prime_dvd_of_missing_lt_total {N m n : ℕ}
    (hmissing : weightedMissingPrimeInterval N m n < weightedPrimeInterval m n) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n ∧ p ∣ N := by
  classical
  by_contra hnone
  have hfilter :
      ((Nat.primesLE n \ Nat.primesLE m).filter fun p => ¬ p ∣ N) =
        Nat.primesLE n \ Nat.primesLE m := by
    refine Finset.filter_true_of_mem ?_
    intro p hp
    have hp_mem_n : p ∈ Nat.primesLE n := (Finset.mem_sdiff.mp hp).1
    have hp_not_mem_m : p ∉ Nat.primesLE m := (Finset.mem_sdiff.mp hp).2
    have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp_mem_n
    have hp_le_n : p ≤ n := Nat.le_of_mem_primesLE hp_mem_n
    have hm_lt_p : m < p := by
      by_contra hnot
      have hp_le_m : p ≤ m := Nat.le_of_not_gt hnot
      exact hp_not_mem_m ((Nat.mem_primesLE).mpr ⟨hp_le_m, hp_prime⟩)
    by_contra hdiv
    exact hnone ⟨p, hp_prime, hm_lt_p, hp_le_n, hdiv⟩
  have heq : weightedMissingPrimeInterval N m n = weightedPrimeInterval m n := by
    rw [weightedMissingPrimeInterval, weightedPrimeInterval, missingPrimeIntervalFinset, hfilter]
  rw [heq] at hmissing
  exact (lt_irrefl _ hmissing).elim

/-- A lower-bound formulation of `exists_prime_dvd_of_missing_lt_total`. -/
theorem exists_prime_dvd_of_missing_lt_interval_lower {N m n : ℕ} {B : ℝ}
    (hlow : B ≤ weightedPrimeInterval m n)
    (hmissing : weightedMissingPrimeInterval N m n < B) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n ∧ p ∣ N :=
  exists_prime_dvd_of_missing_lt_total (lt_of_lt_of_le hmissing hlow)

/-- The elementary identity `log 4 = 2 log 2`, kept as a named helper for explicit
Chebyshev interval arithmetic. -/
theorem log_four_eq_two_mul_log_two : Real.log (4 : ℝ) = 2 * Real.log 2 := by
  rw [show (4 : ℝ) = 2 * 2 by norm_num, Real.log_mul (by norm_num) (by norm_num)]
  ring

/-- Membership in a missing-prime interval gives the corresponding strict interval bounds. -/
theorem bounds_of_mem_missingPrimeIntervalFinset {N m n p : ℕ}
    (hp : p ∈ missingPrimeIntervalFinset N m n) :
    m < p ∧ p ≤ n := by
  classical
  rw [missingPrimeIntervalFinset, Finset.mem_filter, Finset.mem_sdiff] at hp
  rcases hp with ⟨⟨hpn, hpm⟩, _hpN⟩
  have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hpn
  have hpn_le : p ≤ n := Nat.le_of_mem_primesLE hpn
  have hm_lt : m < p := by
    by_contra hnot
    have hp_le_m : p ≤ m := Nat.le_of_not_gt hnot
    exact hpm ((Nat.mem_primesLE).mpr ⟨hp_le_m, hp_prime⟩)
  exact ⟨hm_lt, hpn_le⟩

end MertensLower
end GaussianChain
