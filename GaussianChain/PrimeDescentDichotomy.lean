import GaussianChain.DescentSubcritical
import GaussianChain.SubcriticalLog

namespace GaussianChain
namespace PrimeDescentDichotomy

open DescentSubcritical
open SubcriticalBound

/-- Two-scale logarithmic version of the single-prime descent branch bound.

The split branch uses `Bsplit`; the inert branch uses `Binert`. This matches the final geometric
scales `C * sqrt R / sqrt p` and `C * sqrt R / p`. -/
theorem card_le_of_prime_divisor_descent_branch_log_twoScale
    {M s N p : ℕ} [Fact p.Prime]
    (hp2 : p ≠ 2) (hpN : p ∣ N)
    {a L Bsplit Binert : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hlarge : 4 * s ≤ M)
    (hBsplit : 0 < Bsplit)
    (hfloorS : 0 < Nat.floor (Bsplit ^ 2))
    (hBinert : 0 < Binert)
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
    (p % 4 = 3 ∧
        (M : ℝ) ≤
          ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / Binert) ∨
      (p % 4 = 1 ∧
        (M : ℝ) ≤
          2 * (((2 * s : ℕ) : ℝ) +
            ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / Bsplit)) := by
  rcases PrimeDescent.nat_prime_mod_four_eq_one_or_three_of_ne_two (p := p) hp2 with hp1 | hp3
  · obtain ⟨N', hN', hfactor⟩ :=
      PrimeDescent.exists_pos_nat_factor_of_prime_dvd_norm hN hpN
    have hsmall :=
      SubcriticalLog.split_natFloor_sq_pow_lt_of_log_sub
        (s := s) (N := N) (p := p) (N' := N') (B := Bsplit)
        hfloorS (Fact.out : Nat.Prime p).pos hN' hfactor (hlogSplit N' hN' hfactor)
    exact Or.inr ⟨hp1,
      card_le_two_mul_param_subcritical_bound_after_split_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp1 hM hlarge hBsplit hN' hfactor hsmall hcircle hz hmono hparam hmem⟩
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
    exact Or.inl ⟨hp3,
      card_le_of_param_subcritical_windows_after_inert_descent_scaled
        (M := M) (s := s) (N := N) (N' := N') (p := p)
        hp3 hk hBinert hN' hfactor hsmall hcircle hz hmono hparam hmem⟩

end PrimeDescentDichotomy
end GaussianChain
