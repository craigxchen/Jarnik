import GaussianChain.MertensDichotomy
import GaussianChain.MertensSchedule

/-!
Final asymptotic theorems of the Mertens-interval architecture.

The schedule `s = ⌈(9/4) log R / log log R⌉`, `m₀ = ⌊√(log R)⌋`, `m₁ = ⌊(log R)^{49/50}⌋`
feeds the saturated Mertens dichotomy, giving `M ≤ 4Q log R / log log R` for every
`Q > 9/4`, i.e. every constant `K = 4Q > 9`; the rounded corollary uses `K = 10`.
-/

namespace GaussianChain
namespace MertensMain

open Filter
open SubcriticalBound
open MertensSchedule
open AsymptoticParameters

/-- Mertens-architecture asymptotic bound: every constant strictly larger than the formal
threshold `4 * (9/4) = 9` is admissible. -/
theorem mertens_asymptotic_parametrized_theorem
    (C D : ℝ) (hC : 0 < C) (hD : (9 : ℝ) < D) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        OnCircleUpTo M N z →
        InjectiveUpTo M z →
        (∀ i j, i ≤ j → j < M → t i ≤ t j) →
        (∀ i j, i < M → j < M → SubcriticalBound.gaussianSqDist (z i) (z j) ≤
          (t j - t i) ^ 2) →
        (∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) →
        (M : ℝ) ≤ D * (Real.log R / Real.log (Real.log R)) := by
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hscale := eventually_mertens_scale_hypotheses C hC
  have hthreshold :=
    eventually_mertens_threshold_condition C MertensLower.mertensDeficit hC
      MertensLower.mertensDeficit_le_ten
  have hlower := eventually_mertens_lower_endpoint_log_conditions C hC
  have hm₀ := eventually_two_le_sublogM₀
  have hm₀m₁ := eventually_sublogM₀_le_mertensM₁
  have hsmall := eventually_four_mul_mertensM₁_le_mertensS
  have hs := eventually_mertensS_upper_of_gt (Q := D / 4) (by nlinarith)
  filter_upwards [hR_gt_one, hscale, hthreshold, hlower, hm₀, hm₀m₁, hsmall, hs] with
    R hR_gt_one hscale hthreshold hlower hm₀ hm₀m₁ hsmall hs
  intro M N a L z t hR2 hM hN hLbound hcircle hz hmono hparam hmem
  have hR : 0 < R := lt_trans zero_lt_one hR_gt_one
  have hbound :=
    MertensDichotomy.card_le_four_mul_log_div_loglog_of_mertens_radius_sq_endpoint_bounds
      (M := M) (s := mertensS R) (N := N) (m₀ := sublogM₀ R) (m₁ := mertensM₁ R)
      (C := C) (R := R) (Q := D / 4) (a := a) (L := L) (z := z) (t := t)
      hR2 hscale.1 hscale.2.1 hscale.2.2
      (by simpa using hthreshold)
      (by simpa using hlower.1)
      (by simpa using hlower.2)
      hm₀ hm₀m₁ hM hC hR hLbound hsmall hN hs
      hcircle hz hmono hparam hmem
  calc
    (M : ℝ) ≤ 4 * (D / 4) * (Real.log R / Real.log (Real.log R)) := hbound
    _ = D * (Real.log R / Real.log (Real.log R)) := by ring

/-- Rounded corollary of `mertens_asymptotic_parametrized_theorem`. -/
theorem eventually_jarnik_arc_sublog_mertens (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        OnCircleUpTo M N z →
        InjectiveUpTo M z →
        (∀ i j, i ≤ j → j < M → t i ≤ t j) →
        (∀ i j, i < M → j < M → SubcriticalBound.gaussianSqDist (z i) (z j) ≤
          (t j - t i) ^ 2) →
        (∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) →
        (M : ℝ) ≤ 10 * (Real.log R / Real.log (Real.log R)) :=
  mertens_asymptotic_parametrized_theorem C 10 hC (by norm_num)

/-- Any constant strictly above `10` in the rounded form. -/
theorem eventually_jarnik_arc_sublog_mertens_of_constant_gt_ten
    (C D : ℝ) (hC : 0 < C) (hD : 10 < D) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        OnCircleUpTo M N z →
        InjectiveUpTo M z →
        (∀ i j, i ≤ j → j < M → t i ≤ t j) →
        (∀ i j, i < M → j < M → SubcriticalBound.gaussianSqDist (z i) (z j) ≤
          (t j - t i) ^ 2) →
        (∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) →
        (M : ℝ) ≤ D * (Real.log R / Real.log (Real.log R)) :=
  mertens_asymptotic_parametrized_theorem C D hC (by linarith)

/-- Natural-norm form of the Mertens-architecture bound with an arbitrary constant above the
formal threshold `9`. -/
theorem eventually_jarnik_arc_sublog_mertens_nat_norm_of_constant_gt_nine
    (C D : ℝ) (hC : 0 < C) (hD : (9 : ℝ) < D) :
    ∀ᶠ N : ℕ in atTop,
      ∀ {M : ℕ} {a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ},
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt (Real.sqrt (N : ℝ)) →
        OnCircleUpTo M N z →
        InjectiveUpTo M z →
        (∀ i j, i ≤ j → j < M → t i ≤ t j) →
        (∀ i j, i < M → j < M → SubcriticalBound.gaussianSqDist (z i) (z j) ≤
          (t j - t i) ^ 2) →
        (∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) →
        (M : ℝ) ≤
          D * (Real.log (Real.sqrt (N : ℝ)) /
            Real.log (Real.log (Real.sqrt (N : ℝ)))) := by
  have hR := mertens_asymptotic_parametrized_theorem C D hC hD
  have hN_event := tendsto_sqrt_natCast_atTop.eventually hR
  filter_upwards [hN_event] with N hbound
  intro M a L z t hM hNpos hLbound hcircle hz hmono hparam hmem
  have hR2 : (Real.sqrt (N : ℝ)) ^ 2 = (N : ℝ) := Real.sq_sqrt (by positivity)
  exact hbound hR2 hM hNpos hLbound hcircle hz hmono hparam hmem

end MertensMain
end GaussianChain
