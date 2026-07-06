import GaussianChain.MainSublogBound
import Mathlib.Analysis.SpecialFunctions.Log.Base

namespace GaussianChain
namespace AsymptoticParameters

open Filter Asymptotics
open SubcriticalBound

/-- The lower endpoint of the prime stack: a square-root of `log R`. -/
noncomputable def sublogM₀ (R : ℝ) : ℕ :=
  Nat.floor (Real.sqrt (Real.log R))

theorem sublogM₀_le_sqrt_log {R : ℝ} :
    (sublogM₀ R : ℝ) ≤ Real.sqrt (Real.log R) := by
  unfold sublogM₀
  exact Nat.floor_le (Real.sqrt_nonneg _)

theorem two_le_sublogM₀_of {R : ℝ} (h : (2 : ℝ) ≤ Real.sqrt (Real.log R)) :
    2 ≤ sublogM₀ R := by
  unfold sublogM₀
  exact Nat.le_floor h

theorem log_log_isLittleO_log_rpow {a : ℝ} (ha : 0 < a) :
    (fun R : ℝ => Real.log (Real.log R)) =o[atTop]
      (fun R : ℝ => (Real.log R) ^ a) :=
  (isLittleO_log_rpow_atTop ha).comp_tendsto Real.tendsto_log_atTop

theorem eventually_log_log_le_log_rpow_mul {a ε : ℝ}
    (ha : 0 < a) (hε : 0 < ε) :
    ∀ᶠ R : ℝ in atTop,
      Real.log (Real.log R) ≤ ε * (Real.log R) ^ a := by
  have hbound := (log_log_isLittleO_log_rpow ha).bound hε
  filter_upwards [hbound, Real.tendsto_log_atTop.eventually_ge_atTop 1] with R hbound hW
  have hLL_nonneg : 0 ≤ Real.log (Real.log R) := Real.log_nonneg hW
  have hpow_nonneg : 0 ≤ (Real.log R) ^ a :=
    Real.rpow_nonneg (le_trans zero_le_one hW) a
  rw [Real.norm_eq_abs, abs_of_nonneg hLL_nonneg, Real.norm_eq_abs,
    abs_of_nonneg hpow_nonneg] at hbound
  simpa [mul_comm] using hbound

theorem eventually_log_sq_le_const_mul_self {c : ℝ} (hc : 0 < c) :
    ∀ᶠ R : ℝ in atTop, (Real.log R) ^ 2 ≤ c * R := by
  have hbound := (Real.isLittleO_pow_log_id_atTop (n := 2)).bound hc
  filter_upwards [hbound, eventually_ge_atTop (1 : ℝ)] with R hbound hR
  have hlog_nonneg : 0 ≤ Real.log R := Real.log_nonneg hR
  have hpow_nonneg : 0 ≤ Real.log R ^ 2 := sq_nonneg _
  have hR_nonneg : 0 ≤ R := le_trans zero_le_one hR
  rw [Real.norm_eq_abs, abs_of_nonneg hpow_nonneg] at hbound
  simpa [id, abs_of_nonneg hR_nonneg, mul_comm] using hbound

theorem eventually_two_le_sublogM₀ :
    ∀ᶠ R : ℝ in atTop, 2 ≤ sublogM₀ R := by
  have hsqrt_atTop :
      Tendsto (fun R : ℝ => Real.sqrt (Real.log R)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp Real.tendsto_log_atTop
  filter_upwards [hsqrt_atTop.eventually_ge_atTop (2 : ℝ)] with R hR
  exact two_le_sublogM₀_of hR

/-- The square-root radius associated to an integer norm tends to infinity with the norm. -/
theorem tendsto_sqrt_natCast_atTop :
    Tendsto (fun N : ℕ => Real.sqrt (N : ℝ)) atTop atTop :=
  Real.tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop

end AsymptoticParameters
end GaussianChain
