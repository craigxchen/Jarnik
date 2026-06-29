import GaussianChain.DescentSubcritical
import GaussianChain.IntervalStack
import GaussianChain.SubcriticalLog

namespace GaussianChain
namespace PrimeDescentDichotomy

open DescentSubcritical
open SubcriticalBound

/-- A single prime divisor of the common norm gives the appropriate scaled descent bound.

The result keeps the two residue classes as branches: inert primes divide the whole family by the
rational Gaussian integer `p`, while split primes descend an ordered half-subfamily by a Gaussian
factor above `p`. -/
theorem card_le_of_prime_divisor_descent_branch
    {M s N p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) (hpN : p ∣ N)
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hN : 0 < N)
    (hsmallSplit : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) * (N' : ℤ) →
      (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hsmallInert : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
      (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (p % 4 = 3 ∧
        (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
      (p % 4 = 1 ∧
        (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A)) := by
  rcases PrimeDescent.nat_prime_mod_four_eq_one_or_three_of_ne_two (p := p) hp2 with hp1 | hp3
  · obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_prime_dvd_norm hN hpN
    have hsmall := hsmallSplit N' hN' hfactor
    exact Or.inr ⟨hp1,
      card_le_two_mul_param_subcritical_bound_after_split_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp1 hM hlarge hA hN' hfactor hsmall hcircle hz hmono hparam hmem⟩
  · have hz0 : (z 0).norm = (N : ℤ) :=
      RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle 0)
    obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_inert_prime_dvd_norm
        (p := p) hp3 hN (z := z 0) hz0 hpN
    have hsmall := hsmallInert N' hN' hfactor
    have hk : 2 * s ≤ M := by omega
    exact Or.inl ⟨hp3,
      card_le_of_param_subcritical_windows_after_inert_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp3 hk hA hN' hfactor hsmall hcircle hz hmono hparam hmem⟩

/-- Logarithmic version of `card_le_of_prime_divisor_descent_branch`.

The split and inert subcritical thresholds are supplied as real logarithmic inequalities for
the quotient norms `N / p` and `N / p²`. -/
theorem card_le_of_prime_divisor_descent_branch_log
    {M s N p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) (hpN : p ∣ N)
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hN : 0 < N)
    (hlogSplit : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (p % 4 = 3 ∧
        (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
      (p % 4 = 1 ∧
        (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A)) := by
  exact card_le_of_prime_divisor_descent_branch
    (M := M) (s := s) (N := N) (p := p)
    hp2 hpN hM hlarge hA hN
    (fun N' hN' hfactor =>
      SubcriticalLog.split_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := A)
        hB (Fact.out : Nat.Prime p).pos hN' hfactor (hlogSplit N' hN' hfactor))
    (fun N' hN' hfactor =>
      SubcriticalLog.inert_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := A)
        hB (Fact.out : Nat.Prime p).pos hN' hfactor (hlogInert N' hN' hfactor))
    hcircle hz hmono hparam hmem

/-- Few missing primes on a geometric stack produce a prime divisor and then a scaled descent
bound for that divisor. -/
theorem exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (htotal :
      MertensLower.weightedMissingPrimeInterval N m₀ (IntervalStack.geomLower m₀ k) <
        (k : ℝ) * (Real.log 2 / 2))
    (hsmallSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hsmallInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          (Nat.floor (A ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ IntervalStack.geomLower m₀ k ∧ p ∣ N ∧
      ((p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A))) := by
  obtain ⟨p, hp, hmp, hpU, hpN⟩ :=
    IntervalStack.exists_prime_dvd_of_geometric_stack_total_missing_lt
      (N := N) (m₀ := m₀) (k := k) hk hm₀ herr htotal
  have hp2 : p ≠ 2 := by omega
  letI : Fact p.Prime := ⟨hp⟩
  have hbranch :=
    card_le_of_prime_divisor_descent_branch
      (M := M) (s := s) (N := N) (p := p)
      hp2 hpN hM hlarge hA hN
      (fun N' hN' hfactor => hsmallSplit p N' hp hmp hpU hpN hN' hfactor)
      (fun N' hN' hfactor => hsmallInert p N' hp hmp hpU hpN hN' hfactor)
      hcircle hz hmono hparam hmem
  exact ⟨p, hp, hmp, hpU, hpN, hbranch⟩

/-- Logarithmic version of
`exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt`. -/
theorem exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log
    {M s N m₀ k : ℕ}
    {a L A : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hA : 0 < A)
    (hB : 0 < Nat.floor (A ^ 2))
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (htotal :
      MertensLower.weightedMissingPrimeInterval N m₀ (IntervalStack.geomLower m₀ k) <
        (k : ℝ) * (Real.log 2 / 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ IntervalStack.geomLower m₀ k ∧ p ∣ N ∧
      ((p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A))) := by
  obtain ⟨p, hp, hmp, hpU, hpN⟩ :=
    IntervalStack.exists_prime_dvd_of_geometric_stack_total_missing_lt
      (N := N) (m₀ := m₀) (k := k) hk hm₀ herr htotal
  have hp2 : p ≠ 2 := by omega
  letI : Fact p.Prime := ⟨hp⟩
  have hbranch :=
    card_le_of_prime_divisor_descent_branch_log
      (M := M) (s := s) (N := N) (p := p)
      hp2 hpN hM hlarge hA hB hN
      (fun N' hN' hfactor => hlogSplit p N' hp hmp hpU hpN hN' hfactor)
      (fun N' hN' hfactor => hlogInert p N' hp hmp hpU hpN hN' hfactor)
      hcircle hz hmono hparam hmem
  exact ⟨p, hp, hmp, hpU, hpN, hbranch⟩

/-- Two-scale logarithmic version of `card_le_of_prime_divisor_descent_branch`.

The split branch uses `Asplit`; the inert branch uses `Ainert`. This matches the final geometric
scales `C * sqrt R / sqrt p` and `C * sqrt R / p`. -/
theorem card_le_of_prime_divisor_descent_branch_log_twoScale
    {M s N p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) (hpN : p ∣ N)
    {a L Asplit Ainert : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hAsplit : 0 < Asplit)
    (hBs : 0 < Nat.floor (Asplit ^ 2))
    (hAinert : 0 < Ainert)
    (hBi : 0 < Nat.floor (Ainert ^ 2))
    (hN : 0 < N)
    (hlogSplit : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (Asplit ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (Ainert ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (p % 4 = 3 ∧
        (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert) ∨
      (p % 4 = 1 ∧
        (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit)) := by
  rcases PrimeDescent.nat_prime_mod_four_eq_one_or_three_of_ne_two (p := p) hp2 with hp1 | hp3
  · obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_prime_dvd_norm hN hpN
    have hsmall :=
      SubcriticalLog.split_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := Asplit)
        hBs (Fact.out : Nat.Prime p).pos hN' hfactor (hlogSplit N' hN' hfactor)
    exact Or.inr ⟨hp1,
      card_le_two_mul_param_subcritical_bound_after_split_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp1 hM hlarge hAsplit hN' hfactor hsmall hcircle hz hmono hparam hmem⟩
  · have hz0 : (z 0).norm = (N : ℤ) :=
      RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle 0)
    obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_inert_prime_dvd_norm
        (p := p) hp3 hN (z := z 0) hz0 hpN
    have hsmall :=
      SubcriticalLog.inert_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (A := Ainert)
        hBi (Fact.out : Nat.Prime p).pos hN' hfactor (hlogInert N' hN' hfactor)
    have hk : 2 * s ≤ M := by omega
    exact Or.inl ⟨hp3,
      card_le_of_param_subcritical_windows_after_inert_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp3 hk hAinert hN' hfactor hsmall hcircle hz hmono hparam hmem⟩

/-- Two-scale prime-dependent logarithmic version of the geometric-stack descent branch. -/
theorem exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log_twoScale
    {M s N m₀ k : ℕ}
    {a L : ℝ} {Asplit Ainert : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (htotal :
      MertensLower.weightedMissingPrimeInterval N m₀ (IntervalStack.geomLower m₀ k) <
        (k : ℝ) * (Real.log 2 / 2))
    (hAsplit : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Asplit p)
    (hBs : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Asplit p) ^ 2))
    (hAinert : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Ainert p)
    (hBi : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((Ainert p) ^ 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Asplit p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Ainert p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ IntervalStack.geomLower m₀ k ∧ p ∣ N ∧
      ((p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert p) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit p))) := by
  obtain ⟨p, hp, hmp, hpU, hpN⟩ :=
    IntervalStack.exists_prime_dvd_of_geometric_stack_total_missing_lt
      (N := N) (m₀ := m₀) (k := k) hk hm₀ herr htotal
  have hp2 : p ≠ 2 := by omega
  letI : Fact p.Prime := ⟨hp⟩
  have hbranch :=
    card_le_of_prime_divisor_descent_branch_log_twoScale
      (M := M) (s := s) (N := N) (p := p)
      (Asplit := Asplit p) (Ainert := Ainert p)
      hp2 hpN hM hlarge
      (hAsplit p hp hmp hpU hpN) (hBs p hp hmp hpU hpN)
      (hAinert p hp hmp hpU hpN) (hBi p hp hmp hpU hpN) hN
      (fun N' hN' hfactor => hlogSplit p N' hp hmp hpU hpN hN' hfactor)
      (fun N' hN' hfactor => hlogInert p N' hp hmp hpU hpN hN' hfactor)
      hcircle hz hmono hparam hmem
  exact ⟨p, hp, hmp, hpU, hpN, hbranch⟩

/-- Prime-dependent logarithmic version of the geometric-stack descent branch.

This is the sharper final-use form: after the stack produces a divisor prime `p`, the
subcritical scale may be chosen as `A p`, so the descended length can be measured at the actual
prime rather than at the lower endpoint of the stack. -/
theorem exists_prime_dvd_with_descent_bound_of_geometric_stack_total_missing_lt_log_primeScale
    {M s N m₀ k : ℕ}
    {a L : ℝ} {A : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 0 < k)
    (hm₀ : 0 < m₀)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (htotal :
      MertensLower.weightedMissingPrimeInterval N m₀ (IntervalStack.geomLower m₀ k) <
        (k : ℝ) * (Real.log 2 / 2))
    (hA : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < A p)
    (hB : ∀ p, Nat.Prime p → m₀ < p → p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
      0 < Nat.floor ((A p) ^ 2))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((A p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((A p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : ∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt)))
    (hz : Function.Injective z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    ∃ p, Nat.Prime p ∧ m₀ < p ∧ p ≤ IntervalStack.geomLower m₀ k ∧ p ∣ N ∧
      ((p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A p) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A p))) := by
  obtain ⟨p, hp, hmp, hpU, hpN⟩ :=
    IntervalStack.exists_prime_dvd_of_geometric_stack_total_missing_lt
      (N := N) (m₀ := m₀) (k := k) hk hm₀ herr htotal
  have hp2 : p ≠ 2 := by omega
  letI : Fact p.Prime := ⟨hp⟩
  have hbranch :=
    card_le_of_prime_divisor_descent_branch_log
      (M := M) (s := s) (N := N) (p := p) (A := A p)
      hp2 hpN hM hlarge (hA p hp hmp hpU hpN) (hB p hp hmp hpU hpN) hN
      (fun N' hN' hfactor => hlogSplit p N' hp hmp hpU hpN hN' hfactor)
      (fun N' hN' hfactor => hlogInert p N' hp hmp hpU hpN hN' hfactor)
      hcircle hz hmono hparam hmem
  exact ⟨p, hp, hmp, hpU, hpN, hbranch⟩

/-- A split/inert descent alternative is bounded by a uniform branch estimate once the two
prime-scaled length terms are bounded by the same multiple of `S = 2s`. -/
theorem descent_branch_le_of_prime_scale_terms
    {M s p q : ℕ} {L A : ℝ}
    (htermInert :
      ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A ≤
        (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermSplit :
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A ≤
        (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (hbranch :
      (p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A))) :
    (M : ℝ) ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  rcases hbranch with ⟨_hp3, hM⟩ | ⟨_hp1, hM⟩
  · calc
      (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A := hM
      _ ≤ ((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ) := by
        nlinarith [htermInert]
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        nlinarith [hS_nonneg]
  · calc
      (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A) := hM
      _ ≤ 2 * (((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        nlinarith [htermSplit]
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        ring_nf
        exact le_rfl

/-- Two-scale version of `descent_branch_le_of_prime_scale_terms`. -/
theorem descent_branch_le_of_prime_twoScale_terms
    {M s p q : ℕ} {L Asplit Ainert : ℝ}
    (htermInert :
      ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert ≤
        (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (htermSplit :
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit ≤
        (q : ℝ) * ((2 * s : ℕ) : ℝ))
    (hbranch :
      (p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit))) :
    (M : ℝ) ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  rcases hbranch with ⟨_hp3, hM⟩ | ⟨_hp1, hM⟩
  · calc
      (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Ainert := hM
      _ ≤ ((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ) := by
        nlinarith [htermInert]
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        nlinarith [hS_nonneg]
  · calc
      (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Asplit) := hM
      _ ≤ 2 * (((2 * s : ℕ) : ℝ) + (q : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        nlinarith [htermSplit]
      _ ≤ 2 * (((q + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
        norm_num
        ring_nf
        exact le_rfl

/-- The split/inert descent alternatives are both bounded by the split-scale bound at the lower
endpoint `m₀`.

This is the common estimate used after a divisor prime has been found in a stack interval above
`m₀`: the inert branch is even smaller, while the split branch is controlled by
`sqrt m₀ ≤ sqrt p`. -/
theorem descent_branch_le_uniform_lower_endpoint
    {M s m₀ p : ℕ} {L A : ℝ}
    (hm₀ : 1 ≤ m₀) (hmp : m₀ < p)
    (hL : 0 ≤ L) (hA : 0 < A)
    (hbranch :
      (p % 4 = 3 ∧
          (M : ℝ) ≤
            ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / A) ∨
        (p % 4 = 1 ∧
          (M : ℝ) ≤
            2 * (((2 * s : ℕ) : ℝ) +
              ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / A))) :
    (M : ℝ) ≤
      2 * (((2 * s : ℕ) : ℝ) +
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) / A) := by
  let S : ℝ := ((2 * s : ℕ) : ℝ)
  have hS_nonneg : 0 ≤ S := by positivity
  have hm₀_pos_real : 0 < (m₀ : ℝ) := by exact_mod_cast (lt_of_lt_of_le zero_lt_one hm₀)
  have hp_pos_real : 0 < (p : ℝ) := by
    have hp_pos_nat : 0 < p := lt_trans (by omega : 0 < m₀) hmp
    exact_mod_cast hp_pos_nat
  have hsqrt_m₀_pos : 0 < Real.sqrt (m₀ : ℝ) := Real.sqrt_pos_of_pos hm₀_pos_real
  have hsqrt_p_pos : 0 < Real.sqrt (p : ℝ) := Real.sqrt_pos_of_pos hp_pos_real
  have hsqrt_m₀_le_sqrt_p : Real.sqrt (m₀ : ℝ) ≤ Real.sqrt (p : ℝ) := by
    exact Real.sqrt_le_sqrt (by exact_mod_cast hmp.le)
  have hsqrt_m₀_le_p : Real.sqrt (m₀ : ℝ) ≤ (p : ℝ) := by
    have hm₀_le_p : (m₀ : ℝ) ≤ (p : ℝ) := by exact_mod_cast hmp.le
    have hp_one : (1 : ℝ) ≤ (p : ℝ) := by
      exact_mod_cast (le_trans hm₀ hmp.le)
    rw [Real.sqrt_le_left hp_pos_real.le]
    nlinarith
  have hsplit_scale :
      L / Real.sqrt (p : ℝ) ≤ L / Real.sqrt (m₀ : ℝ) :=
    div_le_div_of_nonneg_left hL hsqrt_m₀_pos hsqrt_m₀_le_sqrt_p
  have hinert_scale :
      L / (p : ℝ) ≤ L / Real.sqrt (m₀ : ℝ) :=
    div_le_div_of_nonneg_left hL hsqrt_m₀_pos hsqrt_m₀_le_p
  have hsplit_term :
      S * (L / Real.sqrt (p : ℝ)) / A ≤
        S * (L / Real.sqrt (m₀ : ℝ)) / A := by
    exact div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_left hsplit_scale hS_nonneg) hA.le
  have hinert_term :
      S * (L / (p : ℝ)) / A ≤
        S * (L / Real.sqrt (m₀ : ℝ)) / A := by
    exact div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_left hinert_scale hS_nonneg) hA.le
  have huniform_nonneg :
      0 ≤ S + S * (L / Real.sqrt (m₀ : ℝ)) / A := by
    have hscale_nonneg : 0 ≤ L / Real.sqrt (m₀ : ℝ) :=
      div_nonneg hL hsqrt_m₀_pos.le
    have hterm_nonneg : 0 ≤ S * (L / Real.sqrt (m₀ : ℝ)) / A :=
      div_nonneg (mul_nonneg hS_nonneg hscale_nonneg) hA.le
    exact add_nonneg hS_nonneg hterm_nonneg
  rcases hbranch with ⟨_hp3, hM⟩ | ⟨_hp1, hM⟩
  · have hbase :
        (M : ℝ) ≤ S + S * (L / Real.sqrt (m₀ : ℝ)) / A := by
      calc
        (M : ℝ) ≤ S + S * (L / (p : ℝ)) / A := by simpa [S] using hM
        _ ≤ S + S * (L / Real.sqrt (m₀ : ℝ)) / A := by
          exact add_le_add (le_refl S) hinert_term
    exact hbase.trans (by nlinarith [huniform_nonneg])
  · have hsplit_add :
        S + S * (L / Real.sqrt (p : ℝ)) / A ≤
          S + S * (L / Real.sqrt (m₀ : ℝ)) / A :=
      add_le_add (le_refl S) hsplit_term
    exact hM.trans (by
      simpa [S] using mul_le_mul_of_nonneg_left hsplit_add (by norm_num : (0 : ℝ) ≤ 2))

end PrimeDescentDichotomy
end GaussianChain
