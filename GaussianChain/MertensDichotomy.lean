import GaussianChain.MertensIntervalMass
import GaussianChain.SaturatedWindows
import GaussianChain.MainSublogBound

/-!
Mertens-interval replacement for the geometric-stack dichotomy.

The Chebyshev stack certified `k * log 2 / 2` of weighted prime mass on a stack of
ratio-8 intervals.  Mertens' first theorem certifies the full mass
`log m₁ - log m₀ - mertensDeficit` of a single interval `(m₀, m₁]`, so the stack,
its height `k`, and the per-interval Chebyshev error hypotheses all disappear.
Combined with the saturated window bounds (scales equal to the descended arc
lengths), the branch bounds become `M ≤ 2s` (missing, inert) and `M ≤ 4s` (split).
-/

namespace GaussianChain
namespace MertensDichotomy

open SubcriticalBound

/-- The "few missing primes" condition on the Mertens interval `(m₀, m₁]`: the missing
primes account for less than the mass guaranteed by Mertens' first theorem. -/
noncomputable def fewMissingMertensCondition (N m₀ m₁ : ℕ) : Prop :=
  MertensLower.weightedMissingPrimeInterval N m₀ m₁ <
    Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit

/-- Mertens-interval dichotomy from a threshold inequality: either the missing mass reaches
the Mertens floor and the weighted-log determinant condition fires, or few primes are
missing. -/
theorem missing_or_few_mertens_dichotomy_of_threshold
    {s N m₀ m₁ : ℕ} {B : ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 *
            (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀ m₁) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    MainDichotomy.manyMissingWeightedLogCondition s N m₀ m₁ B ∨
      fewMissingMertensCondition N m₀ m₁ := by
  classical
  by_cases hfew : fewMissingMertensCondition N m₀ m₁
  · exact Or.inr hfew
  · exact Or.inl
      (MainDichotomy.manyMissingWeightedLogCondition_of_missing_lower_bound
        (s := s) (N := N) (m := m₀) (U := m₁) (B := B)
        (μ := Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit)
        (le_of_not_gt hfew) hthreshold)

/-- Version of `missing_or_few_mertens_dichotomy_of_threshold` using the crude bound
`missingPrimeLogSum ≤ (m₁ + 1) log m₁`. -/
theorem missing_or_few_mertens_dichotomy_of_crude_threshold
    {s N m₀ m₁ : ℕ} {B : ℝ}
    (hm₁ : 1 ≤ m₁)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 *
            (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
            2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ)) :
    MainDichotomy.manyMissingWeightedLogCondition s N m₀ m₁ B ∨
      fewMissingMertensCondition N m₀ m₁ := by
  have hlogsum :=
    MissingPrimeIntervalBranch.missingPrimeLogSum_le_succ_mul_log
      (N := N) (m := m₀) (U := m₁) hm₁
  refine missing_or_few_mertens_dichotomy_of_threshold ?_
  have hle :
      ((s : ℝ) ^ 2 *
          (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
          2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
        ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
      ((s : ℝ) ^ 2 *
          (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
          2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀ m₁) +
        ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    nlinarith [hlogsum]
  exact hthreshold.trans_le hle

set_option linter.style.longLine false in
/-- Saturated missing-prime branch: when the scale `B` is at least the arc length, the
weighted-log condition on `(m, U]` bounds the family by a single window, `M ≤ 2s`. -/
theorem missing_prime_branch_saturated
    {M s N m U : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hLB : L ≤ B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hsmallPrime : 4 * U ≤ s)
    (hN : 0 < N)
    (hweightedLog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * MertensLower.weightedMissingPrimeInterval N m U -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m U) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 2 * s := by
  classical
  let P := MertensLower.missingPrimeIntervalFinset N m U
  refine
    MissingPrimeSubcritical.card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound_saturated
      (M := M) (s := s) (N := N) (a := a) (L := L) (B := B)
      (z := z) (t := t) (P := P)
      hk hB hLB hfloor ?hprime ?hsmallPrime hN ?hpN ?hweightedLog hcircle hz hmono hparam hmem
  · intro p hp
    exact MertensLower.prime_of_mem_missingPrimeIntervalFinset (N := N) (m := m) (U := U) hp
  · intro p hp
    have hpU : p ≤ U :=
      MertensLower.le_upper_of_mem_missingPrimeIntervalFinset (N := N) (m := m) (U := U) hp
    exact (Nat.mul_le_mul_left 4 hpU).trans hsmallPrime
  · intro p hp
    exact MertensLower.int_not_dvd_of_mem_missingPrimeIntervalFinset
      (N := N) (m := m) (U := U) hp
  · simpa [P, MertensLower.weightedMissingPrimeInterval,
      MissingPrimeIntervalBranch.missingPrimeLogSum] using hweightedLog

set_option linter.style.longLine false in
/-- Saturated two-scale descent bound for a single prime divisor: with the split and inert
scales at least the corresponding descended arc lengths, the inert branch gives `M ≤ 2s` and
the split branch `M ≤ 4s`. -/
theorem card_le_of_prime_divisor_descent_branch_log_twoScale_saturated
    {M s N p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) (hpN : p ∣ N)
    {a L Bsplit Binert : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hBsplit : 0 < Bsplit)
    (hLBsplit : L / Real.sqrt (p : ℝ) ≤ Bsplit)
    (hfloorS : 0 < Nat.floor (Bsplit ^ 2))
    (hBinert : 0 < Binert)
    (hLBinert : L / (p : ℝ) ≤ Binert)
    (hfloorI : 0 < Nat.floor (Binert ^ 2))
    (hN : 0 < N)
    (hlogSplit : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (Bsplit ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ N' : ℕ, 0 < N' → (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (Binert ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 4 * s := by
  rcases PrimeDescent.nat_prime_mod_four_eq_one_or_three_of_ne_two (p := p) hp2 with hp1 | hp3
  · obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_prime_dvd_norm hN hpN
    have hsmall :=
      SubcriticalLog.split_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (B := Bsplit)
        hfloorS (Fact.out : Nat.Prime p).pos hN' hfactor (hlogSplit N' hN' hfactor)
    exact
      DescentSubcritical.card_le_two_mul_param_subcritical_bound_after_split_descent_scaled_saturated
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp1 hM hlarge hBsplit hLBsplit hN' hfactor hsmall hcircle hz hmono hparam hmem
  · have hz0 : (z 0).norm = (N : ℤ) :=
      RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle 0 hM)
    obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_inert_prime_dvd_norm
        (p := p) hp3 hN (z := z 0) hz0 hpN
    have hsmall :=
      SubcriticalLog.inert_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (B := Binert)
        hfloorI (Fact.out : Nat.Prime p).pos hN' hfactor (hlogInert N' hN' hfactor)
    have hk : 2 * s ≤ M := by omega
    have h2s :=
      DescentSubcritical.card_le_of_param_subcritical_windows_after_inert_descent_scaled_saturated
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp3 hk hBinert hLBinert hN' hfactor hsmall hcircle hz hmono hparam hmem
    omega

/-- Mertens-interval dichotomy with saturated branch bounds: `M ≤ 4s`.

Either the missing mass on `(m₀, m₁]` reaches the Mertens floor and the saturated missing
branch gives `M ≤ 2s`, or some prime of the interval divides `N` and the saturated descent
branch gives `M ≤ 4s`. -/
theorem card_le_four_mul_s_of_mertens_interval_log
    {M s N m₀ m₁ : ℕ}
    {a L B : ℝ} {Bsplit Binert : ℕ → ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (h2m₀ : 2 ≤ m₀)
    (hm₀m₁ : m₀ ≤ m₁)
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hB : 0 < B)
    (hLB : L ≤ B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hsmallPrime : 4 * m₁ ≤ s)
    (hN : 0 < N)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 *
            (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
            2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hBsplitPos : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N → 0 < Bsplit p)
    (hLBsplit : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      L / Real.sqrt (p : ℝ) ≤ Bsplit p)
    (hfloorS : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      0 < Nat.floor ((Bsplit p) ^ 2))
    (hBinertPos : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N → 0 < Binert p)
    (hLBinert : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      L / (p : ℝ) ≤ Binert p)
    (hfloorI : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      0 < Nat.floor ((Binert p) ^ 2))
    (hlogSplit : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Bsplit p) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p, Nat.Prime p → m₀ < p → p ≤ m₁ → p ∣ N →
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor ((Binert p) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 4 * s := by
  have hm₁ : 1 ≤ m₁ := le_trans (by omega) hm₀m₁
  rcases missing_or_few_mertens_dichotomy_of_crude_threshold (B := B) hm₁ hthreshold with
    hmany | hfew
  · have hk : 2 * s ≤ M := by omega
    have h2s :=
      missing_prime_branch_saturated (M := M) (s := s) (N := N) (m := m₀) (U := m₁)
        hk hB hLB hfloor hsmallPrime hN hmany hcircle hz hmono hparam hmem
    omega
  · obtain ⟨p, hp, hmp, hpm₁, hpN⟩ :=
      MertensLower.exists_prime_dvd_of_missing_lt_log_sub_deficit
        (N := N) (m := m₀) (n := m₁) (by omega) hm₀m₁ hfew
    haveI : Fact p.Prime := ⟨hp⟩
    have hp2 : p ≠ 2 := by omega
    exact
      card_le_of_prime_divisor_descent_branch_log_twoScale_saturated
        (M := M) (s := s) (N := N) (p := p)
        hp2 hpN hM hlarge
        (hBsplitPos p hp hmp hpm₁ hpN)
        (hLBsplit p hp hmp hpm₁ hpN)
        (hfloorS p hp hmp hpm₁ hpN)
        (hBinertPos p hp hmp hpm₁ hpN)
        (hLBinert p hp hmp hpm₁ hpN)
        (hfloorI p hp hmp hpm₁ hpN)
        hN
        (fun N' _ _ => hlogSplit p hp hmp hpm₁ hpN)
        (fun N' _ _ => hlogInert p hp hmp hpm₁ hpN)
        hcircle hz hmono hparam hmem

open MainSublogBound in
/-- Radius-scale form of the Mertens-interval bound.

The missing branch runs at the full scale `C√R`, the descent branches at the descended
lengths `C√R/√p` and `C√R/p`; all three saturate, so `M ≤ 4s ≤ 4Q log R / log log R`. -/
theorem card_le_four_mul_log_div_loglog_of_mertens_radius_sq_endpoint_bounds
    {M s N m₀ m₁ : ℕ}
    {C R Q a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (honeFull : 1 ≤ C ^ 2 * R)
    (honeSplitUpper : 1 ≤ C ^ 2 * R / (m₁ : ℝ))
    (honeInertUpper : 1 ≤ C ^ 2 * R / (m₁ : ℝ) ^ 2)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 *
            (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
            2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hsplitLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)))
    (hinertLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)))
    (h2m₀ : 2 ≤ m₀)
    (hm₀m₁ : m₀ ≤ m₁)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * m₁ ≤ s)
    (hN : 0 < N)
    (hs : (s : ℝ) ≤ Q * (Real.log R / Real.log (Real.log R)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 4 * Q * (Real.log R / Real.log (Real.log R)) := by
  have hm₀ : 0 < m₀ := by omega
  have hCR_pos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hCR_nonneg : 0 ≤ C ^ 2 * R := hCR_pos.le
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  have hlogN := log_nat_eq_two_log_radius_of_sq (N := N) hR hR2
  have h4s : M ≤ 4 * s := by
    by_cases hlarge : 4 * s ≤ M
    · -- the dichotomy applies
      have hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2) :=
        fullRadiusScale_floor_pos_of_one_le hR.le honeFull
      have hlogFull :
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤
            Real.log (C ^ 2 * R) :=
        fullRadiusScale_log_floor_le_log_of_sq_le hR.le hfloorFull le_rfl
      have hthresholdFloor :
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
            ((s : ℝ) ^ 2 *
                (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - MertensLower.mertensDeficit) -
                2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
              ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
        have hstep :=
          (mul_le_mul_of_nonneg_left hlogFull hcoeff_nonneg).trans_lt hthreshold
        simpa [hlogN] using hstep
      have hLB : L ≤ fullRadiusScale C R := hLbound
      refine
        card_le_four_mul_s_of_mertens_interval_log
          (M := M) (s := s) (N := N) (m₀ := m₀) (m₁ := m₁)
          (a := a) (L := L) (B := fullRadiusScale C R)
          (Bsplit := fun p => radiusScale C R p)
          (Binert := fun p => inertRadiusScale C R p) (z := z) (t := t)
          h2m₀ hm₀m₁ hM hlarge (fullRadiusScale_pos hC hR) hLB hfloorFull
          hsmallPrime hN hthresholdFloor
          (fun p hp _hmp _hpm₁ _hpN => radiusScale_pos hC hR hp.pos)
          (fun p hp _hmp _hpm₁ _hpN => by
            change L / Real.sqrt (p : ℝ) ≤ C * Real.sqrt R / Real.sqrt (p : ℝ)
            gcongr)
          (fun p hp _hmp hpm₁ _hpN => by
            have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
            have hdiv : C ^ 2 * R / (m₁ : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
              div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpm₁)
            exact radiusScale_floor_pos_of_one_le hR.le hp.pos
              (honeSplitUpper.trans hdiv))
          (fun p hp _hmp _hpm₁ _hpN => inertRadiusScale_pos hC hR hp.pos)
          (fun p hp _hmp _hpm₁ _hpN => by
            change L / (p : ℝ) ≤ C * Real.sqrt R / (p : ℝ)
            gcongr)
          (fun p hp _hmp hpm₁ _hpN => by
            have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
            have hpm₁real : (p : ℝ) ≤ (m₁ : ℝ) := by exact_mod_cast hpm₁
            have hpm₁_sq : (p : ℝ) ^ 2 ≤ (m₁ : ℝ) ^ 2 :=
              pow_le_pow_left₀ hp_pos_real.le hpm₁real 2
            have hdiv : C ^ 2 * R / (m₁ : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
              div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpm₁_sq
            exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
              (honeInertUpper.trans hdiv))
          (fun p hp hmp hpm₁ hpN => by
            have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
            have hfloorP : 0 < Nat.floor ((radiusScale C R p) ^ 2) := by
              have hdiv : C ^ 2 * R / (m₁ : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
                div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpm₁)
              exact radiusScale_floor_pos_of_one_le hR.le hp.pos
                (honeSplitUpper.trans hdiv)
            have hlogP :
                Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) ≤
                  Real.log (C ^ 2 * R / (p : ℝ)) :=
              radiusScale_log_floor_le_log_of_sq_le hR.le hp.pos hfloorP le_rfl
            have hsplit :=
              split_real_scale_log_condition_of_lower_endpoint
                (s := s) (m₀ := m₀) (p := p) hC hR hm₀ hp hmp hsplitLower
            have hsplitN :
                ((s * (2 * s + 1) : ℕ) : ℝ) *
                    Real.log (C ^ 2 * R / (p : ℝ)) <
                  ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)) := by
              simpa [hlogN] using hsplit
            exact (mul_le_mul_of_nonneg_left hlogP hcoeff_nonneg).trans_lt hsplitN)
          (fun p hp hmp hpm₁ hpN => by
            have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
            have hfloorP : 0 < Nat.floor ((inertRadiusScale C R p) ^ 2) := by
              have hpm₁real : (p : ℝ) ≤ (m₁ : ℝ) := by exact_mod_cast hpm₁
              have hpm₁_sq : (p : ℝ) ^ 2 ≤ (m₁ : ℝ) ^ 2 :=
                pow_le_pow_left₀ hp_pos_real.le hpm₁real 2
              have hdiv : C ^ 2 * R / (m₁ : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
                div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpm₁_sq
              exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
                (honeInertUpper.trans hdiv)
            have hlogP :
                Real.log (Nat.floor ((inertRadiusScale C R p) ^ 2) : ℝ) ≤
                  Real.log (C ^ 2 * R / (p : ℝ) ^ 2) :=
              inertRadiusScale_log_floor_le_log_of_sq_le hR.le hp.pos hfloorP le_rfl
            have hinert :=
              inert_sq_real_scale_log_condition_of_lower_endpoint
                (s := s) (m₀ := m₀) (p := p) hC hR hm₀ hp hmp hinertLower
            have hinertN :
                ((s * (2 * s + 1) : ℕ) : ℝ) *
                    Real.log (C ^ 2 * R / (p : ℝ) ^ 2) <
                  ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)) := by
              simpa [hlogN] using hinert
            exact (mul_le_mul_of_nonneg_left hlogP hcoeff_nonneg).trans_lt hinertN)
          hcircle hz hmono hparam hmem
    · omega
  have h4s_real : (M : ℝ) ≤ 4 * (s : ℝ) := by exact_mod_cast h4s
  calc
    (M : ℝ) ≤ 4 * (s : ℝ) := h4s_real
    _ ≤ 4 * (Q * (Real.log R / Real.log (Real.log R))) := by nlinarith [hs]
    _ = 4 * Q * (Real.log R / Real.log (Real.log R)) := by ring

end MertensDichotomy
end GaussianChain
