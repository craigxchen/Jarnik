import GaussianChain.MertensLower

namespace GaussianChain

open Finset

namespace MertensLower

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

namespace IntervalStack

/-- Lower endpoint of the geometric stack `m₀, 8m₀, 8²m₀, ...`. -/
def geomLower (m₀ i : ℕ) : ℕ := m₀ * 8 ^ i

/-- Upper endpoint of the `i`-th one-eighth interval in the geometric stack. -/
def geomUpper (m₀ i : ℕ) : ℕ := m₀ * 8 ^ (i + 1)

theorem geomLower_le_geomUpper (m₀ i : ℕ) :
    geomLower m₀ i ≤ geomUpper m₀ i := by
  unfold geomLower geomUpper
  exact Nat.mul_le_mul_left _ (Nat.pow_le_pow_right (by norm_num) (Nat.le_succ i))

theorem eight_mul_geomLower_le_geomUpper (m₀ i : ℕ) :
    8 * geomLower m₀ i ≤ geomUpper m₀ i := by
  unfold geomLower geomUpper
  rw [pow_succ]
  have h : 8 * (m₀ * 8 ^ i) = m₀ * (8 ^ i * 8) := by ring
  rw [h]

theorem geomLower_mono_left {m₀ i j : ℕ} (hij : i ≤ j) :
    geomLower m₀ i ≤ geomLower m₀ j := by
  unfold geomLower
  exact Nat.mul_le_mul_left _ (Nat.pow_le_pow_right (by norm_num) hij)

theorem geomUpper_le_stackUpper {m₀ k : ℕ} (i : Fin k) :
    geomUpper m₀ i ≤ geomLower m₀ k := by
  unfold geomUpper geomLower
  exact Nat.mul_le_mul_left _ (Nat.pow_le_pow_right (by norm_num) i.isLt)

theorem geomLower_base_le (m₀ : ℕ) (i : ℕ) :
    m₀ ≤ geomLower m₀ i := by
  have hpow : 1 ≤ 8 ^ i := Nat.one_le_pow i 8 (by norm_num)
  simpa [geomLower] using Nat.mul_le_mul_left m₀ hpow

theorem geomUpper_le_geomLower_of_succ_le {m₀ i j : ℕ} (hij : i + 1 ≤ j) :
    geomUpper m₀ i ≤ geomLower m₀ j := by
  unfold geomUpper geomLower
  exact Nat.mul_le_mul_left _ (Nat.pow_le_pow_right (by norm_num) hij)

/-- The missing-prime intervals in a geometric one-eighth stack are pairwise disjoint. -/
theorem pairwiseDisjoint_missingPrimeIntervalFinset_geom
    (N m₀ k : ℕ) :
    (Set.univ : Set (Fin k)).PairwiseDisjoint fun i =>
      MertensLower.missingPrimeIntervalFinset N (geomLower m₀ i) (geomUpper m₀ i) := by
  classical
  rw [Finset.pairwiseDisjoint_iff]
  intro i _hi j _hj hnonempty
  rcases hnonempty with ⟨p, hp⟩
  rw [Finset.mem_inter] at hp
  rcases hp with ⟨hpi, hpj⟩
  rcases lt_trichotomy (i : ℕ) (j : ℕ) with hij | hij | hij
  · have hbounds_i :=
      MertensLower.bounds_of_mem_missingPrimeIntervalFinset (N := N) hpi
    have hbounds_j :=
      MertensLower.bounds_of_mem_missingPrimeIntervalFinset (N := N) hpj
    have hle : geomUpper m₀ i ≤ geomLower m₀ j :=
      geomUpper_le_geomLower_of_succ_le (m₀ := m₀) (i := i) (j := j) (by omega)
    have : p < p := lt_of_le_of_lt (hbounds_i.2.trans hle) hbounds_j.1
    exact (lt_irrefl p this).elim
  · exact Fin.ext hij
  · have hbounds_i :=
      MertensLower.bounds_of_mem_missingPrimeIntervalFinset (N := N) hpi
    have hbounds_j :=
      MertensLower.bounds_of_mem_missingPrimeIntervalFinset (N := N) hpj
    have hle : geomUpper m₀ j ≤ geomLower m₀ i :=
      geomUpper_le_geomLower_of_succ_le (m₀ := m₀) (i := j) (j := i) (by omega)
    have : p < p := lt_of_le_of_lt (hbounds_j.2.trans hle) hbounds_i.1
    exact (lt_irrefl p this).elim

/-- Every interval in the geometric stack is contained in the total stack interval. -/
theorem missingPrimeIntervalFinset_geom_subset_stack
    {N m₀ k : ℕ} (i : Fin k) :
    MertensLower.missingPrimeIntervalFinset N (geomLower m₀ i) (geomUpper m₀ i) ⊆
      MertensLower.missingPrimeIntervalFinset N m₀ (geomLower m₀ k) := by
  classical
  intro p hp
  rw [MertensLower.missingPrimeIntervalFinset, Finset.mem_filter] at hp ⊢
  rcases hp with ⟨hpint, hpN⟩
  rw [Finset.mem_sdiff] at hpint ⊢
  rcases hpint with ⟨hp_upper, hp_not_lower⟩
  have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hp_upper
  have hp_le_upper : p ≤ geomUpper m₀ i := Nat.le_of_mem_primesLE hp_upper
  have hp_le_stack : p ≤ geomLower m₀ k :=
    hp_le_upper.trans (geomUpper_le_stackUpper (m₀ := m₀) i)
  have hp_stack : p ∈ Nat.primesLE (geomLower m₀ k) :=
    (Nat.mem_primesLE).mpr ⟨hp_le_stack, hp_prime⟩
  have hp_not_base : p ∉ Nat.primesLE m₀ := by
    intro hbase
    have hp_le_lower : p ≤ geomLower m₀ i :=
      (Nat.le_of_mem_primesLE hbase).trans (geomLower_base_le m₀ i)
    exact hp_not_lower ((Nat.mem_primesLE).mpr ⟨hp_le_lower, hp_prime⟩)
  exact ⟨⟨hp_stack, hp_not_base⟩, hpN⟩

/-- If the Chebyshev endpoint error estimate holds for every integer above the base of the
stack, then it holds at every geometric upper endpoint in the stack. -/
theorem chebyshev_error_on_geomUpper_of_base
    {m₀ k : ℕ}
    (hbase : ∀ n, m₀ ≤ n →
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ)) :
    ∀ i : Fin k,
      Real.log (((geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((geomUpper m₀ i : ℕ) : ℝ) := by
  intro i
  exact hbase (geomUpper m₀ i)
    ((geomLower_base_le m₀ i).trans (geomLower_le_geomUpper m₀ i))

/-- A concrete geometric-stack form of the one-eighth interval prime-divisor extraction. -/
theorem exists_prime_dvd_of_geometric_stack_total_missing_lt
    {N m₀ k : ℕ} (hk : 0 < k) (hm₀ : 0 < m₀)
    (herr : ∀ i : Fin k,
      Real.log (((geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((geomUpper m₀ i : ℕ) : ℝ))
    (htotal :
      MertensLower.weightedMissingPrimeInterval N m₀ (geomLower m₀ k) <
        (k : ℝ) * (Real.log 2 / 2)) :
    ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ geomLower m₀ k ∧ p ∣ N := by
  classical
  letI : Nonempty (Fin k) := ⟨⟨0, hk⟩⟩
  obtain ⟨i, p, hp, hmp, hpn, hpdvd⟩ :=
    MertensLower.exists_prime_dvd_of_disjoint_family_total_missing_lt_log_two_half_eighth_interval
      (ι := Fin k) (N := N) (M := m₀) (U := geomLower m₀ k)
      (m := fun i => geomLower m₀ i) (n := fun i => geomUpper m₀ i)
      (fun i => missingPrimeIntervalFinset_geom_subset_stack (N := N) (m₀ := m₀) (k := k) i)
      (pairwiseDisjoint_missingPrimeIntervalFinset_geom N m₀ k)
      (fun i => geomLower_le_geomUpper m₀ i)
      (fun i => by
        unfold geomUpper
        exact Nat.mul_pos hm₀ (Nat.pow_pos (by norm_num)))
      herr
      (fun i => eight_mul_geomLower_le_geomUpper m₀ i)
      (by simpa [Fintype.card_fin] using htotal)
  exact ⟨p, hp, (geomLower_base_le m₀ i).trans_lt hmp,
    hpn.trans (geomUpper_le_stackUpper (m₀ := m₀) i), hpdvd⟩

end IntervalStack
end GaussianChain
