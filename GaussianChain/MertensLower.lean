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

/-- The summand `log p / p` is nonnegative on prime intervals. -/
theorem log_div_nonneg_of_mem_primesLE {p n : ℕ} (hp : p ∈ Nat.primesLE n) :
    0 ≤ Real.log (p : ℝ) / (p : ℝ) := by
  have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp
  exact div_nonneg
    (Real.log_nonneg (by exact_mod_cast hp_prime.one_le))
    (by exact_mod_cast hp_prime.pos.le)

/-- Missing-prime interval sums are nonnegative. -/
theorem weightedMissingPrimeInterval_nonneg (N m n : ℕ) :
    0 ≤ weightedMissingPrimeInterval N m n := by
  classical
  unfold weightedMissingPrimeInterval missingPrimeIntervalFinset
  exact Finset.sum_nonneg fun p hp =>
    log_div_nonneg_of_mem_primesLE (Finset.mem_sdiff.mp (Finset.mem_filter.mp hp).1).1

/-- Monotonicity of missing-prime interval sums under interval inclusion. -/
theorem weightedMissingPrimeInterval_mono {N m₁ m₂ n₁ n₂ : ℕ}
    (hm : m₁ ≤ m₂) (hn : n₂ ≤ n₁) :
    weightedMissingPrimeInterval N m₂ n₂ ≤ weightedMissingPrimeInterval N m₁ n₁ := by
  classical
  unfold weightedMissingPrimeInterval missingPrimeIntervalFinset
  refine Finset.sum_le_sum_of_subset_of_nonneg ?hsub ?hnonneg
  · intro p hp
    rw [Finset.mem_filter] at hp ⊢
    rcases hp with ⟨hpint, hpN⟩
    rw [Finset.mem_sdiff] at hpint ⊢
    rcases hpint with ⟨hpn₂, hpm₂⟩
    have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hpn₂
    have hp_le_n₁ : p ≤ n₁ := (Nat.le_of_mem_primesLE hpn₂).trans hn
    have hpn₁ : p ∈ Nat.primesLE n₁ := (Nat.mem_primesLE).mpr ⟨hp_le_n₁, hp_prime⟩
    have hpm₁ : p ∉ Nat.primesLE m₁ := by
      intro hp_m₁
      have hp_le_m₂ : p ≤ m₂ := (Nat.le_of_mem_primesLE hp_m₁).trans hm
      exact hpm₂ ((Nat.mem_primesLE).mpr ⟨hp_le_m₂, hp_prime⟩)
    exact ⟨⟨hpn₁, hpm₁⟩, hpN⟩
  · intro p hp _hnot
    exact log_div_nonneg_of_mem_primesLE
      (Finset.mem_sdiff.mp (Finset.mem_filter.mp hp).1).1

/-- A disjoint family of missing-prime intervals contributes at most the missing-prime sum over
any enclosing interval. -/
theorem sum_weightedMissingPrimeInterval_le_of_pairwiseDisjoint
    {ι : Type*} [Fintype ι]
    {N M U : ℕ} {m n : ι → ℕ}
    (hsub : ∀ i, missingPrimeIntervalFinset N (m i) (n i) ⊆ missingPrimeIntervalFinset N M U)
    (hdisj :
      (Set.univ : Set ι).PairwiseDisjoint fun i => missingPrimeIntervalFinset N (m i) (n i)) :
    (∑ i, weightedMissingPrimeInterval N (m i) (n i)) ≤
      weightedMissingPrimeInterval N M U := by
  classical
  let S : ι → Finset ℕ := fun i => missingPrimeIntervalFinset N (m i) (n i)
  let T : Finset ℕ := missingPrimeIntervalFinset N M U
  let w : ℕ → ℝ := fun p => Real.log (p : ℝ) / (p : ℝ)
  have hunion_sub : (Finset.univ.biUnion S) ⊆ T := by
    rw [Finset.biUnion_subset]
    intro i _hi
    exact hsub i
  have hsum_union :
      ∑ p ∈ Finset.univ.biUnion S, w p =
        ∑ i, ∑ p ∈ S i, w p := by
    have hdisjS : (((Finset.univ : Finset ι) : Set ι)).PairwiseDisjoint S := by
      simpa [S] using hdisj
    simpa [S, w] using (Finset.sum_biUnion (s := (Finset.univ : Finset ι))
      (t := S) (f := w) hdisjS)
  have hle_union :
      ∑ p ∈ Finset.univ.biUnion S, w p ≤ ∑ p ∈ T, w p := by
    refine Finset.sum_le_sum_of_subset_of_nonneg hunion_sub ?_
    intro p hp _hnot
    exact log_div_nonneg_of_mem_primesLE
      (Finset.mem_sdiff.mp (Finset.mem_filter.mp hp).1).1
  calc
    (∑ i, weightedMissingPrimeInterval N (m i) (n i))
        = ∑ i, ∑ p ∈ S i, w p := by
            simp [weightedMissingPrimeInterval, S, w]
    _ = ∑ p ∈ Finset.univ.biUnion S, w p := hsum_union.symm
    _ ≤ ∑ p ∈ T, w p := hle_union
    _ = weightedMissingPrimeInterval N M U := by
            simp [weightedMissingPrimeInterval, T, w]

/-- The Chebyshev `θ` increment over `m < p ≤ n` is the corresponding sum of `log p`. -/
theorem theta_sub_eq_sum_primesLE_sdiff {m n : ℕ} (hmn : m ≤ n) :
    Chebyshev.theta (n : ℝ) - Chebyshev.theta (m : ℝ) =
      ∑ p ∈ Nat.primesLE n \ Nat.primesLE m, Real.log (p : ℝ) := by
  rw [Chebyshev.theta_eq_sum_primesLE_log n, Chebyshev.theta_eq_sum_primesLE_log m]
  exact (Finset.sum_sdiff_eq_sub (Nat.primesLE_mono hmn)).symm

/-- On an interval `m < p ≤ n`, replacing every denominator `p` by the larger denominator `n`
can only decrease the weighted prime reciprocal sum. -/
theorem theta_sub_div_le_weightedPrimeInterval {m n : ℕ} (hmn : m ≤ n) :
    (Chebyshev.theta (n : ℝ) - Chebyshev.theta (m : ℝ)) / (n : ℝ) ≤
      weightedPrimeInterval m n := by
  rw [theta_sub_eq_sum_primesLE_sdiff hmn, weightedPrimeInterval, Finset.sum_div]
  refine sum_le_sum ?_
  intro p hp
  have hp_mem : p ∈ Nat.primesLE n := (Finset.mem_sdiff.mp hp).1
  have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp_mem
  have hp_le_n : p ≤ n := Nat.le_of_mem_primesLE hp_mem
  have hlog_nonneg : 0 ≤ Real.log (p : ℝ) :=
    Real.log_nonneg (by exact_mod_cast hp_prime.one_le)
  exact div_le_div_of_nonneg_left hlog_nonneg
    (by exact_mod_cast hp_prime.pos) (by exact_mod_cast hp_le_n)

/-- A reusable finite-interval Mertens lower-bound transfer: lower-bound `θ n`, upper-bound `θ m`,
then divide by `n` to get a weighted prime reciprocal interval lower bound. -/
theorem weightedPrimeInterval_lower_of_theta_bounds
    {m n : ℕ} (hmn : m ≤ n) (hn : 0 < n) {lowerN upperM : ℝ}
    (hlo : lowerN ≤ Chebyshev.theta (n : ℝ))
    (hhi : Chebyshev.theta (m : ℝ) ≤ upperM) :
    (lowerN - upperM) / (n : ℝ) ≤ weightedPrimeInterval m n := by
  have htheta : lowerN - upperM ≤
      Chebyshev.theta (n : ℝ) - Chebyshev.theta (m : ℝ) :=
    sub_le_sub hlo hhi
  exact (div_le_div_of_nonneg_right htheta (by exact_mod_cast Nat.le_of_lt hn)).trans
    (theta_sub_div_le_weightedPrimeInterval hmn)

/-- An explicit Chebyshev-based lower bound for the weighted prime reciprocal sum on `m < p ≤ n`.

This is not yet the full Mertens theorem; it is the finite interval estimate that the eventual
Mertens proof will sum over a geometric sequence of intervals. -/
theorem weightedPrimeInterval_lower_explicit (m n : ℕ) (hmn : m ≤ n) (hn : 0 < n) :
    (((n : ℝ) * Real.log 2 - Real.log ((n : ℝ) + 1) -
        2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ)) -
        Real.log 4 * (m : ℝ)) / (n : ℝ) ≤
      weightedPrimeInterval m n := by
  refine weightedPrimeInterval_lower_of_theta_bounds hmn hn ?_ ?_
  · simpa using Chebyshev.theta_ge n
  · simpa [mul_comm] using
      Chebyshev.theta_le_log4_mul_x (x := (m : ℝ)) (by positivity)

/-- A positive weighted prime interval contains an actual rational prime in that interval. -/
theorem exists_prime_of_weightedPrimeInterval_pos {m n : ℕ}
    (hpos : 0 < weightedPrimeInterval m n) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n := by
  classical
  rw [weightedPrimeInterval] at hpos
  have hnonneg :
      ∀ p ∈ Nat.primesLE n \ Nat.primesLE m,
        0 ≤ Real.log (p : ℝ) / (p : ℝ) := by
    intro p hp
    have hp_mem : p ∈ Nat.primesLE n := (Finset.mem_sdiff.mp hp).1
    have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp_mem
    exact div_nonneg
      (Real.log_nonneg (by exact_mod_cast hp_prime.one_le))
      (by exact_mod_cast hp_prime.pos.le)
  obtain ⟨p, hp, _hterm⟩ := (Finset.sum_pos_iff_of_nonneg hnonneg).mp hpos
  have hp_mem_n : p ∈ Nat.primesLE n := (Finset.mem_sdiff.mp hp).1
  have hp_not_mem_m : p ∉ Nat.primesLE m := (Finset.mem_sdiff.mp hp).2
  have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp_mem_n
  have hp_le_n : p ≤ n := Nat.le_of_mem_primesLE hp_mem_n
  have hm_lt_p : m < p := by
    by_contra hnot
    have hp_le_m : p ≤ m := Nat.le_of_not_gt hnot
    exact hp_not_mem_m ((Nat.mem_primesLE).mpr ⟨hp_le_m, hp_prime⟩)
  exact ⟨p, hp_prime, hm_lt_p, hp_le_n⟩

/-- Positivity corollary for the explicit Chebyshev lower bound. -/
theorem weightedPrimeInterval_pos_of_explicit_pos {m n : ℕ}
    (hmn : m ≤ n) (hn : 0 < n)
    (hpos : 0 < (((n : ℝ) * Real.log 2 - Real.log ((n : ℝ) + 1) -
        2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ)) -
        Real.log 4 * (m : ℝ)) / (n : ℝ)) :
    0 < weightedPrimeInterval m n :=
  hpos.trans_le (weightedPrimeInterval_lower_explicit m n hmn hn)

/-- Prime-existence corollary for a positive explicit Chebyshev interval lower bound. -/
theorem exists_prime_of_explicit_interval_lower_pos {m n : ℕ}
    (hmn : m ≤ n) (hn : 0 < n)
    (hpos : 0 < (((n : ℝ) * Real.log 2 - Real.log ((n : ℝ) + 1) -
        2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ)) -
        Real.log 4 * (m : ℝ)) / (n : ℝ)) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n :=
  exists_prime_of_weightedPrimeInterval_pos
    (weightedPrimeInterval_pos_of_explicit_pos hmn hn hpos)

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

/-- If the lower endpoint is at most one eighth of the upper endpoint, then the Chebyshev
upper-endpoint contribution `log 4 * m` is at most one quarter of `n log 2`. -/
theorem log_four_mul_le_quarter_log_two_mul_of_eight_mul_le {m n : ℕ}
    (hm : 8 * m ≤ n) :
    Real.log 4 * (m : ℝ) ≤ (Real.log 2 / 4) * (n : ℝ) := by
  have hlog2_nonneg : 0 ≤ Real.log (2 : ℝ) := le_of_lt (Real.log_pos (by norm_num))
  have hcoeff : 0 ≤ Real.log (2 : ℝ) / 4 := div_nonneg hlog2_nonneg (by norm_num)
  have hcast : (8 : ℝ) * (m : ℝ) ≤ (n : ℝ) := by exact_mod_cast hm
  calc
    Real.log 4 * (m : ℝ)
        = (Real.log 2 / 4) * ((8 : ℝ) * (m : ℝ)) := by
          rw [log_four_eq_two_mul_log_two]
          ring
    _ ≤ (Real.log 2 / 4) * (n : ℝ) :=
          mul_le_mul_of_nonneg_left hcast hcoeff

/-- A packaged finite-interval lower bound: once the explicit Chebyshev error terms consume at
most one quarter of the main term, and the lower endpoint consumes at most another quarter, the
weighted prime reciprocal sum on `m < p ≤ n` is at least `log 2 / 2`. -/
theorem weightedPrimeInterval_lower_log_two_half_of_error_bounds
    {m n : ℕ} (hmn : m ≤ n) (hn : 0 < n)
    (herr :
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ))
    (hm :
      Real.log 4 * (m : ℝ) ≤ (Real.log 2 / 4) * (n : ℝ)) :
    Real.log 2 / 2 ≤ weightedPrimeInterval m n := by
  have hbase := weightedPrimeInterval_lower_explicit m n hmn hn
  have hnpos : 0 < (n : ℝ) := by exact_mod_cast hn
  have hdiv :
      Real.log 2 / 2 ≤
        (((n : ℝ) * Real.log 2 - Real.log ((n : ℝ) + 1) -
            2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ)) -
            Real.log 4 * (m : ℝ)) / (n : ℝ) := by
    rw [le_div_iff₀ hnpos]
    nlinarith [herr, hm]
  exact hdiv.trans hbase

/-- One-eighth interval corollary of the explicit Chebyshev lower bound. -/
theorem weightedPrimeInterval_lower_log_two_half_of_eight_mul_le
    {m n : ℕ} (hmn : m ≤ n) (hn : 0 < n)
    (herr :
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ))
    (hm8 : 8 * m ≤ n) :
    Real.log 2 / 2 ≤ weightedPrimeInterval m n :=
  weightedPrimeInterval_lower_log_two_half_of_error_bounds hmn hn herr
    (log_four_mul_le_quarter_log_two_mul_of_eight_mul_le hm8)

/-- A prime-existence form of the one-eighth Chebyshev interval corollary. -/
theorem exists_prime_of_eighth_interval_error_bound
    {m n : ℕ} (hmn : m ≤ n) (hn : 0 < n)
    (herr :
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ))
    (hm8 : 8 * m ≤ n) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n := by
  have hlow := weightedPrimeInterval_lower_log_two_half_of_eight_mul_le hmn hn herr hm8
  have hposlog : 0 < Real.log (2 : ℝ) / 2 := half_pos (Real.log_pos (by norm_num))
  exact exists_prime_of_weightedPrimeInterval_pos (hposlog.trans_le hlow)

/-- If the missing-prime contribution is below `log 2 / 2`, then the Chebyshev one-eighth
interval contains a prime divisor of `N`. -/
theorem exists_prime_dvd_of_missing_lt_log_two_half_eighth_interval
    {N m n : ℕ} (hmn : m ≤ n) (hn : 0 < n)
    (herr :
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ))
    (hm8 : 8 * m ≤ n)
    (hmissing : weightedMissingPrimeInterval N m n < Real.log 2 / 2) :
    ∃ p, Nat.Prime p ∧ m < p ∧ p ≤ n ∧ p ∣ N :=
  exists_prime_dvd_of_missing_lt_interval_lower
    (weightedPrimeInterval_lower_log_two_half_of_eight_mul_le hmn hn herr hm8)
    hmissing

/-- A finite-average pigeonhole principle over real-valued data. -/
theorem exists_lt_of_sum_lt_card_mul
    {ι : Type*} [Fintype ι] [Nonempty ι] {f : ι → ℝ} {c : ℝ}
    (hsum : (∑ i, f i) < (Fintype.card ι : ℝ) * c) :
    ∃ i, f i < c := by
  classical
  by_contra hnone
  have hall : ∀ i, c ≤ f i := by
    intro i
    exact le_of_not_gt fun hi => hnone ⟨i, hi⟩
  have hle : (∑ _i : ι, c) ≤ ∑ i, f i := Finset.sum_le_sum fun i _ => hall i
  have hcard : (∑ _i : ι, c) = (Fintype.card ι : ℝ) * c := by
    simp [Finset.sum_const, nsmul_eq_mul]
  rw [hcard] at hle
  exact not_lt_of_ge hle hsum

/-- Finite-family form of the one-eighth interval prime-divisor extraction.

If the average missing-prime contribution over a finite family of one-eighth intervals is below
`log 2 / 2`, then one of those intervals contains a prime divisor of `N`. -/
theorem exists_prime_dvd_of_family_missing_average_lt_log_two_half_eighth_interval
    {ι : Type*} [Fintype ι] [Nonempty ι]
    {N : ℕ} {m n : ι → ℕ}
    (hmn : ∀ i, m i ≤ n i)
    (hn : ∀ i, 0 < n i)
    (herr : ∀ i,
      Real.log (((n i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((n i : ℕ) : ℝ) * Real.log ((n i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((n i : ℕ) : ℝ))
    (hm8 : ∀ i, 8 * m i ≤ n i)
    (havg :
      (∑ i, weightedMissingPrimeInterval N (m i) (n i)) <
        (Fintype.card ι : ℝ) * (Real.log 2 / 2)) :
    ∃ i p, Nat.Prime p ∧ m i < p ∧ p ≤ n i ∧ p ∣ N := by
  obtain ⟨i, hi⟩ := exists_lt_of_sum_lt_card_mul (f :=
    fun i => weightedMissingPrimeInterval N (m i) (n i)) (c := Real.log 2 / 2) havg
  obtain ⟨p, hp, hmp, hpn, hpdvd⟩ :=
    exists_prime_dvd_of_missing_lt_log_two_half_eighth_interval
      (N := N) (m := m i) (n := n i) (hmn i) (hn i) (herr i) (hm8 i) hi
  exact ⟨i, p, hp, hmp, hpn, hpdvd⟩

/-- Disjoint-family form of the one-eighth interval extraction controlled by a larger missing
prime interval. -/
theorem exists_prime_dvd_of_disjoint_family_total_missing_lt_log_two_half_eighth_interval
    {ι : Type*} [Fintype ι] [Nonempty ι]
    {N M U : ℕ} {m n : ι → ℕ}
    (hsub : ∀ i, missingPrimeIntervalFinset N (m i) (n i) ⊆ missingPrimeIntervalFinset N M U)
    (hdisj :
      (Set.univ : Set ι).PairwiseDisjoint fun i => missingPrimeIntervalFinset N (m i) (n i))
    (hmn : ∀ i, m i ≤ n i)
    (hn : ∀ i, 0 < n i)
    (herr : ∀ i,
      Real.log (((n i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((n i : ℕ) : ℝ) * Real.log ((n i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((n i : ℕ) : ℝ))
    (hm8 : ∀ i, 8 * m i ≤ n i)
    (htotal :
      weightedMissingPrimeInterval N M U <
        (Fintype.card ι : ℝ) * (Real.log 2 / 2)) :
    ∃ i p, Nat.Prime p ∧ m i < p ∧ p ≤ n i ∧ p ∣ N := by
  have hsum_le :=
    sum_weightedMissingPrimeInterval_le_of_pairwiseDisjoint
      (N := N) (M := M) (U := U) (m := m) (n := n) hsub hdisj
  exact exists_prime_dvd_of_family_missing_average_lt_log_two_half_eighth_interval
    hmn hn herr hm8 (lt_of_le_of_lt hsum_le htotal)

/-- The elementary analytic estimate used by the Chebyshev-interval Mertens step:
`sqrt x * log x` is `o(x)`. -/
theorem sqrt_mul_log_isLittleO_id :
    (fun x : ℝ => Real.sqrt x * Real.log x) =o[atTop] (fun x => x) := by
  have hsqrt :
      (fun x : ℝ => Real.sqrt x) =O[atTop] fun x => x ^ (1 / 2 : ℝ) := by
    simpa [Real.sqrt_eq_rpow] using
      (isBigO_refl (fun x : ℝ => x ^ (1 / 2 : ℝ)) atTop)
  have hlog :
      (fun x : ℝ => Real.log x) =o[atTop] fun x => x ^ (1 / 2 : ℝ) :=
    isLittleO_log_rpow_atTop one_half_pos
  have hmul := hsqrt.mul_isLittleO hlog
  refine hmul.congr' EventuallyEq.rfl ?_
  filter_upwards [eventually_ge_atTop (0 : ℝ)] with x hx
  rw [← Real.rpow_add_of_nonneg hx (by norm_num) (by norm_num),
    show (1 / 2 : ℝ) + 1 / 2 = 1 by norm_num, Real.rpow_one]

/-- A shifted endpoint logarithm is still sublinear. -/
theorem log_add_one_isLittleO_id :
    (fun x : ℝ => Real.log (x + 1)) =o[atTop] (fun x => x) := by
  have hshift : Tendsto (fun x : ℝ => x + 1) atTop atTop :=
    tendsto_atTop_add_const_right atTop 1 tendsto_id
  have hlog_shift :
      (fun x : ℝ => Real.log (x + 1)) =o[atTop] fun x => x + 1 :=
    Real.isLittleO_log_id_atTop.comp_tendsto hshift
  have hshiftO : (fun x : ℝ => x + 1) =O[atTop] (fun x => x) := by
    refine IsBigO.of_bound 2 ?_
    filter_upwards [eventually_ge_atTop (1 : ℝ)] with x hx
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg (by linarith),
      abs_of_nonneg (by linarith)]
    linarith
  exact hlog_shift.trans_isBigO hshiftO

/-- Natural-number form of `sqrt_mul_log_isLittleO_id`. -/
theorem sqrt_mul_log_nat_isLittleO_nat :
    (fun n : ℕ => Real.sqrt (n : ℝ) * Real.log (n : ℝ)) =o[atTop]
      (fun n => (n : ℝ)) :=
  sqrt_mul_log_isLittleO_id.natCast_atTop

/-- Natural-number form of `log_add_one_isLittleO_id`. -/
theorem log_nat_add_one_isLittleO_nat :
    (fun n : ℕ => Real.log ((n : ℝ) + 1)) =o[atTop] (fun n => (n : ℝ)) :=
  log_add_one_isLittleO_id.natCast_atTop

/-- The explicit Chebyshev endpoint error is eventually at most one quarter of the main
`n log 2` term. -/
theorem chebyshev_error_eventually_le_quarter_log_two :
    ∀ᶠ n : ℕ in atTop,
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ) := by
  have hlog2_pos : 0 < Real.log (2 : ℝ) := Real.log_pos (by norm_num)
  have hlog_coeff : 0 < Real.log (2 : ℝ) / 8 := div_pos hlog2_pos (by norm_num)
  have hprod_coeff : 0 < Real.log (2 : ℝ) / 16 := div_pos hlog2_pos (by norm_num)
  have hlog := log_nat_add_one_isLittleO_nat.bound hlog_coeff
  have hprod := sqrt_mul_log_nat_isLittleO_nat.bound hprod_coeff
  filter_upwards [hlog, hprod, eventually_ge_atTop 1] with n hlog hprod hn
  have hn_nonneg : 0 ≤ (n : ℝ) := by positivity
  have hn_one : 1 ≤ (n : ℝ) := by exact_mod_cast hn
  have hnp1_one : 1 ≤ (n : ℝ) + 1 := by linarith
  have hlog_nonneg : 0 ≤ Real.log ((n : ℝ) + 1) := Real.log_nonneg hnp1_one
  have hprod_nonneg : 0 ≤ Real.sqrt (n : ℝ) * Real.log (n : ℝ) :=
    mul_nonneg (Real.sqrt_nonneg _) (Real.log_nonneg hn_one)
  rw [Real.norm_eq_abs, abs_of_nonneg hlog_nonneg, Real.norm_eq_abs,
    abs_of_nonneg hn_nonneg] at hlog
  rw [Real.norm_eq_abs, abs_of_nonneg hprod_nonneg, Real.norm_eq_abs,
    abs_of_nonneg hn_nonneg] at hprod
  nlinarith [hlog, hprod]

/-- Eventual prime existence in one-eighth intervals, with the Chebyshev endpoint error supplied
by `chebyshev_error_eventually_le_quarter_log_two`. -/
theorem eventually_exists_prime_of_eighth_interval {m : ℕ → ℕ} :
    ∀ᶠ n : ℕ in atTop,
      m n ≤ n → 0 < n → 8 * m n ≤ n →
        ∃ p, Nat.Prime p ∧ m n < p ∧ p ≤ n := by
  filter_upwards [chebyshev_error_eventually_le_quarter_log_two] with n herr hmn hn hm8
  exact exists_prime_of_eighth_interval_error_bound hmn hn herr hm8

/-- Eventual one-eighth interval prime-divisor extraction from a `log 2 / 2` missing-prime
upper bound. -/
theorem eventually_exists_prime_dvd_of_missing_lt_log_two_half_eighth_interval
    {N m : ℕ → ℕ} :
    ∀ᶠ n : ℕ in atTop,
      m n ≤ n → 0 < n → 8 * m n ≤ n →
        weightedMissingPrimeInterval (N n) (m n) n < Real.log 2 / 2 →
          ∃ p, Nat.Prime p ∧ m n < p ∧ p ≤ n ∧ p ∣ N n := by
  filter_upwards [chebyshev_error_eventually_le_quarter_log_two] with n herr hmn hn hm8 hmissing
  exact exists_prime_dvd_of_missing_lt_log_two_half_eighth_interval hmn hn herr hm8 hmissing

end MertensLower
end GaussianChain
