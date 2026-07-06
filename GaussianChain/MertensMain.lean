import GaussianChain.MertensDichotomy
import GaussianChain.MertensSchedule

/-!
Final asymptotic theorems of the Mertens-interval architecture.

Given `D > 8`, put `Q = D/4 > 2`, fix `q = (2 + Q)/2 ∈ (2, Q)` and
`θ₁ = 3/4 + 1/(2q) ∈ (1/2 + 1/q, 1)`.  The schedule `s = ⌈q log R / log log R⌉`,
`m₀ = ⌊√(log R)⌋`, `m₁ = ⌊(log R)^θ₁⌋` feeds the saturated Mertens dichotomy, giving
`M ≤ 4Q log R / log log R = D log R / log log R` for every constant `D > 8`;
the rounded corollary uses `D = 10`.
-/

namespace GaussianChain
namespace MertensMain

open Filter
open SubcriticalBound
open MertensSchedule
open AsymptoticParameters

/-- Mertens-architecture asymptotic bound: every constant strictly larger than the formal
threshold `4 * 2 = 8` is admissible. -/
theorem mertens_asymptotic_parametrized_theorem
    (C D : ℝ) (hC : 0 < C) (hD : (8 : ℝ) < D) :
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
  have hQ2 : 2 < D / 4 := by linarith
  set q : ℝ := (2 + D / 4) / 2 with hq_def
  have hq2 : 2 < q := by rw [hq_def]; linarith
  have hq0 : 0 < q := by linarith
  have hqQ : q < D / 4 := by rw [hq_def]; linarith
  have h2q_pos : (0 : ℝ) < 2 * q := by linarith
  have hx : 1 / (2 * q) < 1 / 4 := by
    rw [div_lt_div_iff₀ h2q_pos (by norm_num : (0 : ℝ) < 4)]
    linarith
  have hx0 : (0 : ℝ) < 1 / (2 * q) := by positivity
  have h1q : (1 : ℝ) / q = 2 * (1 / (2 * q)) := by
    rw [mul_one_div, div_eq_div_iff (ne_of_gt hq0) (ne_of_gt h2q_pos)]
    ring
  set θ₁ : ℝ := 3 / 4 + 1 / (2 * q) with hθ₁_def
  have hθ0 : 0 < θ₁ := by rw [hθ₁_def]; linarith
  have hθhalf : 1 / 2 < θ₁ := by rw [hθ₁_def]; linarith
  have hθ1 : θ₁ < 1 := by rw [hθ₁_def]; linarith
  have hθgap : 1 / 2 + 1 / q < θ₁ := by rw [hθ₁_def]; linarith
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hscale := eventually_mertens_scale_hypotheses C θ₁ hC hθ0 hθ1
  have hthreshold :=
    eventually_mertens_threshold_condition C MertensLower.mertensDeficit q θ₁ hC
      hq2 hθgap hθ1
  have hlower := eventually_mertens_lower_endpoint_log_conditions C q hC hq2
  have hm₀ := eventually_two_le_sublogM₀
  have hm₀m₁ := eventually_sublogM₀_le_mertensM₁ hθhalf
  have hsmall := eventually_four_mul_mertensM₁_le_mertensS hq0 hθ1
  have hs := eventually_mertensS_upper_of_gt hq0 hqQ
  filter_upwards [hR_gt_one, hscale, hthreshold, hlower, hm₀, hm₀m₁, hsmall, hs] with
    R hR_gt_one hscale hthreshold hlower hm₀ hm₀m₁ hsmall hs
  intro M N a L z t hR2 hM hN hLbound hcircle hz hmono hparam hmem
  have hR : 0 < R := lt_trans zero_lt_one hR_gt_one
  have hbound :=
    MertensDichotomy.card_le_four_mul_log_div_loglog_of_mertens_radius_sq_endpoint_bounds
      (M := M) (s := mertensS q R) (N := N) (m₀ := sublogM₀ R) (m₁ := mertensM₁ θ₁ R)
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
formal threshold `8`. -/
theorem eventually_jarnik_arc_sublog_mertens_nat_norm_of_constant_gt_eight
    (C D : ℝ) (hC : 0 < C) (hD : (8 : ℝ) < D) :
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
