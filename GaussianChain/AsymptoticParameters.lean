import GaussianChain.MainSublogBound
import Mathlib.Analysis.SpecialFunctions.Log.Base

namespace GaussianChain
namespace AsymptoticParameters

open Filter Asymptotics

/-- The window parameter used in the final sublogarithmic bound. -/
noncomputable def sublogS (R : ℝ) : ℕ :=
  Nat.ceil (1024 * (Real.log R / Real.log (Real.log R)))

/-- The lower endpoint of the prime stack: a square-root of `log R`. -/
noncomputable def sublogM₀ (R : ℝ) : ℕ :=
  Nat.floor (Real.sqrt (Real.log R))

/-- The height of the geometric prime stack. -/
noncomputable def sublogK (R : ℝ) : ℕ :=
  Nat.floor (Real.log (Real.log R) / (8 * Real.log 8))

theorem sublogS_lower (R : ℝ) :
    1024 * (Real.log R / Real.log (Real.log R)) ≤ (sublogS R : ℝ) := by
  unfold sublogS
  exact Nat.le_ceil _

theorem sublogS_upper_of_nonneg
    {R : ℝ}
    (hX_nonneg : 0 ≤ 1024 * (Real.log R / Real.log (Real.log R)))
    (hX_one : 1 ≤ 1024 * (Real.log R / Real.log (Real.log R))) :
    (sublogS R : ℝ) ≤ 2048 * (Real.log R / Real.log (Real.log R)) := by
  have hceil_lt :
      (sublogS R : ℝ) <
        1024 * (Real.log R / Real.log (Real.log R)) + 1 := by
    unfold sublogS
    exact Nat.ceil_lt_add_one hX_nonneg
  linarith

theorem sublogM₀_le_sqrt_log {R : ℝ} :
    (sublogM₀ R : ℝ) ≤ Real.sqrt (Real.log R) := by
  unfold sublogM₀
  exact Nat.floor_le (Real.sqrt_nonneg _)

theorem two_le_sublogM₀_of {R : ℝ} (h : (2 : ℝ) ≤ Real.sqrt (Real.log R)) :
    2 ≤ sublogM₀ R := by
  unfold sublogM₀
  exact Nat.le_floor h

theorem sublogK_lower_of {R : ℝ}
    (h :
      2 ≤ Real.log (Real.log R) / (8 * Real.log 8)) :
    Real.log (Real.log R) / (16 * Real.log 8) ≤ (sublogK R : ℝ) := by
  unfold sublogK
  have hfloor : (1 : ℕ) ≤ Nat.floor (Real.log (Real.log R) / (8 * Real.log 8)) :=
    Nat.le_floor (by linarith)
  have hfloor_real :
      (1 : ℝ) ≤ (Nat.floor (Real.log (Real.log R) / (8 * Real.log 8)) : ℝ) := by
    exact_mod_cast hfloor
  have hfloor_add :
      (Nat.floor (Real.log (Real.log R) / (8 * Real.log 8)) : ℝ) + 1 ≤
        2 * (Nat.floor (Real.log (Real.log R) / (8 * Real.log 8)) : ℝ) := by
    linarith
  have hlt :
      Real.log (Real.log R) / (8 * Real.log 8) <
        (Nat.floor (Real.log (Real.log R) / (8 * Real.log 8)) : ℝ) + 1 :=
    Nat.lt_floor_add_one _
  have hhalf :
      Real.log (Real.log R) / (16 * Real.log 8) =
        (Real.log (Real.log R) / (8 * Real.log 8)) / 2 := by ring
  rw [hhalf]
  linarith

theorem sublogK_pow_le_exp {R : ℝ}
    (hy_nonneg : 0 ≤ Real.log (Real.log R) / (8 * Real.log 8)) :
    (((8 ^ sublogK R : ℕ) : ℝ)) ≤
      Real.exp (Real.log (Real.log R) / 8) := by
  have hbase : (1 : ℝ) ≤ 8 := by norm_num
  have hk_le :
      (sublogK R : ℝ) ≤ Real.log (Real.log R) / (8 * Real.log 8) := by
    unfold sublogK
    exact Nat.floor_le hy_nonneg
  calc
    (((8 ^ sublogK R : ℕ) : ℝ)) = (8 : ℝ) ^ sublogK R := by norm_num
    _ = (8 : ℝ) ^ ((sublogK R : ℕ) : ℝ) := by
      exact (Real.rpow_natCast (8 : ℝ) (sublogK R)).symm
    _ ≤ (8 : ℝ) ^ (Real.log (Real.log R) / (8 * Real.log 8)) :=
      Real.rpow_le_rpow_of_exponent_le hbase hk_le
    _ = Real.exp (Real.log (Real.log R) / 8) := by
      rw [Real.rpow_def_of_pos (by norm_num : (0 : ℝ) < 8)]
      congr 1
      have hlog8 : Real.log (8 : ℝ) ≠ 0 := by
        exact (Real.log_pos (by norm_num : (1 : ℝ) < 8)).ne'
      field_simp [hlog8]

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

theorem eventually_sublogS_bounds :
    ∀ᶠ R : ℝ in atTop,
      1024 * (Real.log R / Real.log (Real.log R)) ≤ (sublogS R : ℝ) ∧
        (sublogS R : ℝ) ≤ 2048 * (Real.log R / Real.log (Real.log R)) := by
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (1 / 2048) * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1)
      (by norm_num : (0 : ℝ) < 1 / 2048)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hloglog_pos, hLL_small, hlog_gt_one] with R hLL hsmall hW_gt_one
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  rw [Real.rpow_one] at hsmall
  have hratio_large : 2048 ≤ Real.log R / Real.log (Real.log R) := by
    rw [le_div_iff₀ hLL]
    nlinarith
  have hX : 1 ≤ 1024 * (Real.log R / Real.log (Real.log R)) := by
    nlinarith
  have hX_nonneg : 0 ≤ 1024 * (Real.log R / Real.log (Real.log R)) := by positivity
  exact ⟨sublogS_lower R, sublogS_upper_of_nonneg hX_nonneg hX⟩

theorem eventually_two_le_sublogM₀ :
    ∀ᶠ R : ℝ in atTop, 2 ≤ sublogM₀ R := by
  have hsqrt_atTop :
      Tendsto (fun R : ℝ => Real.sqrt (Real.log R)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp Real.tendsto_log_atTop
  filter_upwards [hsqrt_atTop.eventually_ge_atTop (2 : ℝ)] with R hR
  exact two_le_sublogM₀_of hR

theorem eventually_sublogK_positive_and_lower :
    ∀ᶠ R : ℝ in atTop,
      0 < sublogK R ∧
        Real.log (Real.log R) / (16 * Real.log 8) ≤ (sublogK R : ℝ) := by
  have hlog8_pos : 0 < Real.log (8 : ℝ) := Real.log_pos (by norm_num)
  have htarget :
      ∀ᶠ R : ℝ in atTop,
        2 * (8 * Real.log 8) ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (2 * (8 * Real.log 8))
  filter_upwards [htarget] with R hR
  have hy : 2 ≤ Real.log (Real.log R) / (8 * Real.log 8) := by
    rw [le_div_iff₀ (by positivity : (0 : ℝ) < 8 * Real.log 8)]
    simpa [mul_assoc] using hR
  have hk_one : 1 ≤ sublogK R := by
    unfold sublogK
    exact Nat.le_floor (by linarith)
  exact ⟨lt_of_lt_of_le zero_lt_one hk_one, sublogK_lower_of hy⟩

theorem eventually_sublogK_pow_le_log_rpow :
    ∀ᶠ R : ℝ in atTop,
      (((8 ^ sublogK R : ℕ) : ℝ)) ≤ (Real.log R) ^ (1 / 8 : ℝ) := by
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hlog_gt_one] with R hW_gt_one
  have hLL_nonneg : 0 ≤ Real.log (Real.log R) :=
    Real.log_nonneg hW_gt_one.le
  have hy_nonneg : 0 ≤ Real.log (Real.log R) / (8 * Real.log 8) := by positivity
  calc
    (((8 ^ sublogK R : ℕ) : ℝ)) ≤
        Real.exp (Real.log (Real.log R) / 8) :=
      sublogK_pow_le_exp hy_nonneg
    _ = (Real.log R) ^ (1 / 8 : ℝ) := by
      rw [Real.rpow_def_of_pos (lt_trans zero_lt_one hW_gt_one)]
      congr 1
      ring

theorem eventually_smallPrime_sublog_parameters :
    ∀ᶠ R : ℝ in atTop,
      4 * IntervalStack.geomLower (sublogM₀ R) (sublogK R) ≤ sublogS R := by
  have hs_bounds := eventually_sublogS_bounds
  have hpow_bound := eventually_sublogK_pow_le_log_rpow
  have hLL_bound :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ 256 * (Real.log R) ^ (3 / 8 : ℝ) :=
    eventually_log_log_le_log_rpow_mul
      (by norm_num : (0 : ℝ) < 3 / 8)
      (by norm_num : (0 : ℝ) < 256)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  filter_upwards [hs_bounds, hpow_bound, hLL_bound, hlog_gt_one, hloglog_pos] with
    R hs hpow hLL hW_gt_one hLL_pos
  let W := Real.log R
  let U := IntervalStack.geomLower (sublogM₀ R) (sublogK R)
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hW_nonneg : 0 ≤ W := hW_pos.le
  have hm₀_le : (sublogM₀ R : ℝ) ≤ Real.sqrt W := by
    simpa [W] using (sublogM₀_le_sqrt_log (R := R))
  have hU_le :
      (U : ℝ) ≤ Real.sqrt W * W ^ (1 / 8 : ℝ) := by
    dsimp [U]
    unfold IntervalStack.geomLower
    rw [Nat.cast_mul]
    exact mul_le_mul hm₀_le (by simpa [W] using hpow)
      (by positivity) (by positivity)
  have hpowprod : W ^ (1 / 8 : ℝ) * W ^ (3 / 8 : ℝ) = Real.sqrt W := by
    rw [← Real.rpow_add hW_pos]
    norm_num
    rw [← Real.sqrt_eq_rpow]
  have hscale :
      4 * (Real.sqrt W * W ^ (1 / 8 : ℝ)) ≤
        1024 * (W / Real.log W) := by
    rw [show 1024 * (W / Real.log W) = (1024 * W) / Real.log W by ring]
    rw [le_div_iff₀ (by simpa [W] using hLL_pos)]
    calc
      (4 * (Real.sqrt W * W ^ (1 / 8 : ℝ))) * Real.log W
          ≤ (4 * (Real.sqrt W * W ^ (1 / 8 : ℝ))) *
              (256 * W ^ (3 / 8 : ℝ)) := by
            exact mul_le_mul_of_nonneg_left (by simpa [W] using hLL) (by positivity)
      _ = 1024 * W := by
            calc
              4 * (Real.sqrt W * W ^ (1 / 8 : ℝ)) *
                    (256 * W ^ (3 / 8 : ℝ)) =
                  1024 * (Real.sqrt W * (W ^ (1 / 8 : ℝ) * W ^ (3 / 8 : ℝ))) := by
                    ring
              _ = 1024 * (Real.sqrt W * Real.sqrt W) := by rw [hpowprod]
              _ = 1024 * W := by
                    rw [← pow_two, Real.sq_sqrt hW_nonneg]
  have hfourU :
      ((4 * U : ℕ) : ℝ) ≤ (sublogS R : ℝ) := by
    have hfourU_real : (4 : ℝ) * (U : ℝ) ≤ (sublogS R : ℝ) := by
      calc
        (4 : ℝ) * (U : ℝ) ≤ 4 * (Real.sqrt W * W ^ (1 / 8 : ℝ)) := by
          exact mul_le_mul_of_nonneg_left hU_le (by norm_num)
        _ ≤ 1024 * (W / Real.log W) := hscale
        _ = 1024 * (Real.log R / Real.log (Real.log R)) := by rfl
        _ ≤ (sublogS R : ℝ) := hs.1
    simpa [Nat.cast_mul] using hfourU_real
  exact_mod_cast hfourU

theorem eventually_scale_hypotheses (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      let U := IntervalStack.geomLower (sublogM₀ R) (sublogK R)
      1 ≤ C ^ 2 * R ∧
        1 ≤ C ^ 2 * R / (U : ℝ) ∧
          1 ≤ C ^ 2 * R / (U : ℝ) ^ 2 := by
  have hs_bounds := eventually_sublogS_bounds
  have hsmall := eventually_smallPrime_sublog_parameters
  have hm₀_two := eventually_two_le_sublogM₀
  have hloglog_ge_one :
      ∀ᶠ R : ℝ in atTop, 1 ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop 1
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_large : ∀ᶠ R : ℝ in atTop, 1 ≤ C ^ 2 * R := by
    have htendsto : Tendsto (fun R : ℝ => C ^ 2 * R) atTop atTop :=
      Tendsto.const_mul_atTop (sq_pos_of_pos hC) tendsto_id
    exact htendsto.eventually_ge_atTop 1
  have hlog_sq_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log R ^ 2 ≤ (C ^ 2 / (2048 : ℝ) ^ 2) * R :=
    eventually_log_sq_le_const_mul_self
      (div_pos (sq_pos_of_pos hC) (sq_pos_of_pos (by norm_num : (0 : ℝ) < 2048)))
  filter_upwards [hs_bounds, hsmall, hm₀_two, hloglog_ge_one, hlog_gt_one, hR_large,
    hlog_sq_small] with R hs hsmall hm₀_two hLL_ge_one hW_gt_one hfull hlog_sq
  let U := IntervalStack.geomLower (sublogM₀ R) (sublogK R)
  let s := sublogS R
  let W := Real.log R
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hW_nonneg : 0 ≤ W := hW_pos.le
  have hLL_pos : 0 < Real.log W := lt_of_lt_of_le zero_lt_one (by simpa [W] using hLL_ge_one)
  have hX_le_W : W / Real.log W ≤ W := by
    rw [div_le_iff₀ hLL_pos]
    nlinarith
  have hs_le_W : (s : ℝ) ≤ 2048 * W := by
    calc
      (s : ℝ) ≤ 2048 * (Real.log R / Real.log (Real.log R)) := hs.2
      _ = 2048 * (W / Real.log W) := by rfl
      _ ≤ 2048 * W := by nlinarith
  have hU_le_s_nat : U ≤ s := by
    have : 4 * U ≤ s := by simpa [U, s] using hsmall
    omega
  have hU_le_s : (U : ℝ) ≤ (s : ℝ) := by exact_mod_cast hU_le_s_nat
  have hU_pos_nat : 0 < U := by
    dsimp [U]
    unfold IntervalStack.geomLower
    exact Nat.mul_pos (by omega) (Nat.pow_pos (by norm_num))
  have hU_pos : 0 < (U : ℝ) := by exact_mod_cast hU_pos_nat
  have hbig : (2048 * W) ^ 2 ≤ C ^ 2 * R := by
    have h2048sq_pos : (0 : ℝ) < (2048 : ℝ) ^ 2 := sq_pos_of_pos (by norm_num)
    have hmul :
        (2048 : ℝ) ^ 2 * W ^ 2 ≤ (2048 : ℝ) ^ 2 * ((C ^ 2 / (2048 : ℝ) ^ 2) * R) :=
      mul_le_mul_of_nonneg_left (by simpa [W] using hlog_sq) h2048sq_pos.le
    have hcancel :
        (2048 : ℝ) ^ 2 * ((C ^ 2 / (2048 : ℝ) ^ 2) * R) = C ^ 2 * R := by
      field_simp [(ne_of_gt h2048sq_pos)]
    calc
      (2048 * W) ^ 2 = (2048 : ℝ) ^ 2 * W ^ 2 := by ring
      _ ≤ (2048 : ℝ) ^ 2 * ((C ^ 2 / (2048 : ℝ) ^ 2) * R) := hmul
      _ = C ^ 2 * R := hcancel
  have hU_le_CR : (U : ℝ) ≤ C ^ 2 * R := by
    have hU_le_2048W : (U : ℝ) ≤ 2048 * W := hU_le_s.trans hs_le_W
    have h2048W_ge_one : 1 ≤ 2048 * W := by nlinarith
    calc
      (U : ℝ) ≤ 2048 * W := hU_le_2048W
      _ ≤ (2048 * W) ^ 2 := by nlinarith [sq_nonneg (2048 * W)]
      _ ≤ C ^ 2 * R := hbig
  have hU_sq_le_CR : (U : ℝ) ^ 2 ≤ C ^ 2 * R := by
    have hU_le_2048W : (U : ℝ) ≤ 2048 * W := hU_le_s.trans hs_le_W
    calc
      (U : ℝ) ^ 2 ≤ (2048 * W) ^ 2 := by
        nlinarith [hU_pos.le, hU_le_2048W]
      _ ≤ C ^ 2 * R := hbig
  exact ⟨hfull, (one_le_div hU_pos).mpr hU_le_CR,
    (one_le_div (sq_pos_of_pos hU_pos)).mpr hU_sq_le_CR⟩

theorem eventually_log_sublogM₀_lower :
    ∀ᶠ R : ℝ in atTop,
      (1 / 4 : ℝ) * Real.log (Real.log R) ≤ Real.log (sublogM₀ R : ℝ) := by
  have hW14_atTop :
      Tendsto (fun R : ℝ => (Real.log R) ^ (1 / 4 : ℝ)) atTop atTop :=
    (tendsto_rpow_atTop (by norm_num : (0 : ℝ) < 1 / 4)).comp Real.tendsto_log_atTop
  have hW14_large :
      ∀ᶠ R : ℝ in atTop, 2 ≤ (Real.log R) ^ (1 / 4 : ℝ) :=
    hW14_atTop.eventually_ge_atTop 2
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hW14_large, hlog_gt_one] with R hW14_large hW_gt_one
  let W := Real.log R
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hW_nonneg : 0 ≤ W := hW_pos.le
  have hW14_pos : 0 < W ^ (1 / 4 : ℝ) :=
    Real.rpow_pos_of_pos hW_pos _
  have hW14_sq : (W ^ (1 / 4 : ℝ)) ^ 2 = Real.sqrt W := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul hW_nonneg]
    norm_num
    rw [← Real.sqrt_eq_rpow]
  have hW14_add_le_sqrt : W ^ (1 / 4 : ℝ) + 1 ≤ Real.sqrt W := by
    calc
      W ^ (1 / 4 : ℝ) + 1 ≤
          W ^ (1 / 4 : ℝ) + W ^ (1 / 4 : ℝ) := by linarith
      _ ≤ (W ^ (1 / 4 : ℝ)) * (W ^ (1 / 4 : ℝ)) := by
            nlinarith
      _ = Real.sqrt W := by
            rw [← pow_two, hW14_sq]
  have hfloor_lt :
      Real.sqrt W < (sublogM₀ R : ℝ) + 1 := by
    unfold sublogM₀
    simpa [W] using Nat.lt_floor_add_one (Real.sqrt (Real.log R))
  have hW14_le_m₀ : W ^ (1 / 4 : ℝ) ≤ (sublogM₀ R : ℝ) := by
    linarith
  calc
    (1 / 4 : ℝ) * Real.log (Real.log R)
        = Real.log (W ^ (1 / 4 : ℝ)) := by
          rw [Real.log_rpow hW_pos]
    _ ≤ Real.log (sublogM₀ R : ℝ) :=
        Real.log_le_log hW14_pos hW14_le_m₀

theorem split_lower_endpoint_algebra
    {S W LL logm c : ℝ}
    (hS : 1 ≤ S)
    (hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : 1024 * W ≤ S * LL)
    (hlogm : LL / 4 ≤ logm)
    (hc : 2 * c ≤ LL / 16) :
    (S * (2 * S + 1)) * (2 * c + W - logm) <
      S ^ 2 * (2 * W - logm) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hconst :
      (S * (2 * S + 1)) * (2 * c) ≤ (3 / 16) * S ^ 2 * LL := by
    have hcoeff_nonneg : 0 ≤ S * (2 * S + 1) := by nlinarith
    have hcoeff_le : S * (2 * S + 1) ≤ 3 * S ^ 2 := by nlinarith
    calc
      (S * (2 * S + 1)) * (2 * c) ≤
          (S * (2 * S + 1)) * (LL / 16) :=
        mul_le_mul_of_nonneg_left hc hcoeff_nonneg
      _ ≤ (3 * S ^ 2) * (LL / 16) := by
        exact mul_le_mul_of_nonneg_right hcoeff_le (by nlinarith)
      _ = (3 / 16) * S ^ 2 * LL := by ring
  have hloggain : S ^ 2 * (LL / 4) ≤ S ^ 2 * logm := by
    exact mul_le_mul_of_nonneg_left hlogm (sq_nonneg S)
  have hSW_small : S * W < (1 / 16) * S ^ 2 * LL := by
    have hmul : 1024 * S * W ≤ S ^ 2 * LL := by
      nlinarith
    nlinarith [hSpos, hW, hLL]
  nlinarith

theorem inert_lower_endpoint_algebra
    {S W LL logm c : ℝ}
    (hS : 1 ≤ S)
    (hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : 1024 * W ≤ S * LL)
    (hlogm : LL / 4 ≤ logm)
    (hc : 2 * c ≤ LL / 16) :
    (S * (2 * S + 1)) * (2 * c + W - 2 * logm) <
      S ^ 2 * (2 * W - 2 * logm) := by
  have hsplit :=
    split_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := c)
      hS hW hLL hSLL hlogm hc
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hextra_pos :
      0 < (S * (2 * S + 1)) * logm - S ^ 2 * logm := by
    have hlogm_pos : 0 < logm := by nlinarith
    have hcoeff_gap : 0 < S ^ 2 + S := by nlinarith
    nlinarith
  nlinarith

theorem crude_threshold_algebra
    {S W LL K U logU c : ℝ}
    (hS : 1 ≤ S)
    (hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : 1024 * W ≤ S * LL)
    (hKgain : LL / 96 ≤ K * (Real.log 2 / 2))
    (hc : 2 * c ≤ LL / 1024)
    (hUterm : 2 * ((U + 1) * logU) ≤ 4 * S * W) :
    (S * (2 * S + 1)) * (2 * c + W) <
      (S ^ 2 * (K * (Real.log 2 / 2)) - 2 * ((U + 1) * logU)) +
        S ^ 2 * (2 * W) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hstack : S ^ 2 * (LL / 96) ≤ S ^ 2 * (K * (Real.log 2 / 2)) :=
    mul_le_mul_of_nonneg_left hKgain (sq_nonneg S)
  have hconst :
      (S * (2 * S + 1)) * (2 * c) ≤ (3 / 1024) * S ^ 2 * LL := by
    have hcoeff_nonneg : 0 ≤ S * (2 * S + 1) := by nlinarith
    have hcoeff_le : S * (2 * S + 1) ≤ 3 * S ^ 2 := by nlinarith
    calc
      (S * (2 * S + 1)) * (2 * c) ≤
          (S * (2 * S + 1)) * (LL / 1024) :=
        mul_le_mul_of_nonneg_left hc hcoeff_nonneg
      _ ≤ (3 * S ^ 2) * (LL / 1024) := by
        exact mul_le_mul_of_nonneg_right hcoeff_le (by nlinarith)
      _ = (3 / 1024) * S ^ 2 * LL := by ring
  have hSW_small : S * W ≤ (1 / 1024) * S ^ 2 * LL := by
    have hmul : 1024 * S * W ≤ S ^ 2 * LL := by nlinarith
    nlinarith
  have herr_small :
      S * W + 2 * ((U + 1) * logU) + (S * (2 * S + 1)) * (2 * c) ≤
        (8 / 1024) * S ^ 2 * LL := by
    nlinarith
  have hgap : (8 / 1024 : ℝ) * S ^ 2 * LL < S ^ 2 * (LL / 96) := by
    have hSsq_pos : 0 < S ^ 2 := sq_pos_of_pos hSpos
    nlinarith
  nlinarith

theorem eventually_lower_endpoint_log_conditions (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      let s := sublogS R
      let m₀ := sublogM₀ R
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)) ∧
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)) := by
  have hs_bounds := eventually_sublogS_bounds
  have hlogm_lower := eventually_log_sublogM₀_lower
  have hm₀_two := eventually_two_le_sublogM₀
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (1 / 2048) * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1)
      (by norm_num : (0 : ℝ) < 1 / 2048)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hconst :
      ∀ᶠ R : ℝ in atTop, 32 * Real.log C ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (32 * Real.log C)
  filter_upwards [hs_bounds, hlogm_lower, hm₀_two, hLL_small, hlog_gt_one, hR_gt_one,
    hconst] with R hs hlogm_lower hm₀_two hLL_small hW_gt_one hR_gt_one hconst
  let s := sublogS R
  let m₀ := sublogM₀ R
  let S : ℝ := s
  let W := Real.log R
  let LL := Real.log W
  let logm := Real.log (m₀ : ℝ)
  have hR_pos : 0 < R := lt_trans zero_lt_one hR_gt_one
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hLL_pos : 0 < LL := by
    dsimp [LL, W]
    exact Real.log_pos hW_gt_one
  have hLL_small' : LL ≤ (1 / 2048) * W := by
    simpa [LL, W, Real.rpow_one] using hLL_small
  have hratio_large : 2048 ≤ W / LL := by
    have hLL_mul : 2048 * LL ≤ W := by
      calc
        2048 * LL ≤ 2048 * ((1 / 2048 : ℝ) * W) :=
          mul_le_mul_of_nonneg_left hLL_small' (by norm_num)
        _ = W := by ring
    rw [le_div_iff₀ hLL_pos]
    exact hLL_mul
  have hS_one : 1 ≤ S := by
    have hs_lower : 1024 * (W / LL) ≤ S := by
      simpa [s, S, W, LL] using hs.1
    nlinarith
  have hSLL : 1024 * W ≤ S * LL := by
    have hs_lower : 1024 * (W / LL) ≤ S := by
      simpa [s, S, W, LL] using hs.1
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : 1024 * (W / LL) * LL = 1024 * W := by
      field_simp [hLL_pos.ne']
    nlinarith
  have hlogm : LL / 4 ≤ logm := by
    have hraw : (1 / 4 : ℝ) * LL ≤ logm := by
      simpa [m₀, LL, W, logm] using hlogm_lower
    nlinarith
  have hc : 2 * Real.log C ≤ LL / 16 := by
    dsimp [LL, W] at hconst ⊢
    nlinarith
  have hm_pos_nat : 0 < m₀ := by omega
  have hm_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm_pos_nat
  have hCR_pos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR_pos
  have hlog_split_scale :
      Real.log (C ^ 2 * R / (m₀ : ℝ)) =
        2 * Real.log C + W - logm := by
    rw [Real.log_div hCR_pos.ne' hm_pos.ne']
    rw [Real.log_mul (pow_ne_zero 2 hC.ne') hR_pos.ne']
    rw [Real.log_pow, show (2 : ℕ) * Real.log C = 2 * Real.log C by norm_num]
  have hlog_inert_scale :
      Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) =
        2 * Real.log C + W - 2 * logm := by
    rw [Real.log_div hCR_pos.ne' (pow_ne_zero 2 hm_pos.ne')]
    rw [Real.log_mul (pow_ne_zero 2 hC.ne') hR_pos.ne']
    rw [Real.log_pow, Real.log_pow,
      show (2 : ℕ) * Real.log C = 2 * Real.log C by norm_num,
      show (2 : ℕ) * Real.log (m₀ : ℝ) = 2 * Real.log (m₀ : ℝ) by norm_num]
  have hcoeff :
      (((s * (2 * s + 1) : ℕ) : ℝ)) = S * (2 * S + 1) := by
    simp [S]
  have hs_sq : (((s * s : ℕ) : ℝ)) = S ^ 2 := by
    simp [S]
    ring
  have hsplit :=
    split_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C)
      hS_one hW_pos hLL_pos hSLL hlogm hc
  have hinert :=
    inert_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C)
      hS_one hW_pos hLL_pos hSLL hlogm hc
  constructor
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_split_scale] using hsplit
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_inert_scale] using hinert

theorem eventually_crude_threshold_condition (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      let s := sublogS R
      let m₀ := sublogM₀ R
      let k := sublogK R
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R) := by
  have hs_bounds := eventually_sublogS_bounds
  have hsmall := eventually_smallPrime_sublog_parameters
  have hk_lower := eventually_sublogK_positive_and_lower
  have hm₀_two := eventually_two_le_sublogM₀
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (1 / 2048) * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1)
      (by norm_num : (0 : ℝ) < 1 / 2048)
  have hLL_ge_one :
      ∀ᶠ R : ℝ in atTop, 1 ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop 1
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hconst :
      ∀ᶠ R : ℝ in atTop, 2048 * Real.log C ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (2048 * Real.log C)
  have hlog_sq_small :
      ∀ᶠ R : ℝ in atTop, Real.log R ^ 2 ≤ (1 / 2048 : ℝ) * R :=
    eventually_log_sq_le_const_mul_self (by norm_num : (0 : ℝ) < 1 / 2048)
  filter_upwards [hs_bounds, hsmall, hk_lower, hm₀_two, hLL_small, hLL_ge_one,
    hlog_gt_one, hR_gt_one, hconst, hlog_sq_small] with
    R hs hsmall hk hm₀_two hLL_small hLL_ge_one hW_gt_one hR_gt_one hconst hlog_sq
  let s := sublogS R
  let m₀ := sublogM₀ R
  let k := sublogK R
  let U := IntervalStack.geomLower m₀ k
  let S : ℝ := s
  let W := Real.log R
  let LL := Real.log W
  let K : ℝ := k
  let logU := Real.log (U : ℝ)
  have hR_pos : 0 < R := lt_trans zero_lt_one hR_gt_one
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hLL_pos : 0 < LL := by
    dsimp [LL, W]
    exact Real.log_pos hW_gt_one
  have hLL_small' : LL ≤ (1 / 2048) * W := by
    simpa [LL, W, Real.rpow_one] using hLL_small
  have hratio_large : 2048 ≤ W / LL := by
    have hLL_mul : 2048 * LL ≤ W := by
      calc
        2048 * LL ≤ 2048 * ((1 / 2048 : ℝ) * W) :=
          mul_le_mul_of_nonneg_left hLL_small' (by norm_num)
        _ = W := by ring
    rw [le_div_iff₀ hLL_pos]
    exact hLL_mul
  have hS_one : 1 ≤ S := by
    have hs_lower : 1024 * (W / LL) ≤ S := by
      simpa [s, S, W, LL] using hs.1
    nlinarith
  have hSLL : 1024 * W ≤ S * LL := by
    have hs_lower : 1024 * (W / LL) ≤ S := by
      simpa [s, S, W, LL] using hs.1
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : 1024 * (W / LL) * LL = 1024 * W := by
      field_simp [hLL_pos.ne']
    nlinarith
  have hc : 2 * Real.log C ≤ LL / 1024 := by
    dsimp [LL, W] at hconst ⊢
    nlinarith
  have hlog8 : Real.log (8 : ℝ) = 3 * Real.log (2 : ℝ) := by
    rw [show (8 : ℝ) = 2 * 2 * 2 by norm_num]
    rw [Real.log_mul (by norm_num : (2 : ℝ) * 2 ≠ 0) (by norm_num : (2 : ℝ) ≠ 0)]
    rw [Real.log_mul (by norm_num : (2 : ℝ) ≠ 0) (by norm_num : (2 : ℝ) ≠ 0)]
    ring
  have hlog2_pos : 0 < Real.log (2 : ℝ) := Real.log_pos (by norm_num)
  have hKgain : LL / 96 ≤ K * (Real.log 2 / 2) := by
    have hk_lower' : LL / (16 * Real.log 8) ≤ K := by
      simpa [k, K, LL, W] using hk.2
    calc
      LL / 96 = (LL / (16 * Real.log 8)) * (Real.log 2 / 2) := by
        rw [hlog8]
        field_simp [hlog2_pos.ne']
        ring
      _ ≤ K * (Real.log 2 / 2) :=
        mul_le_mul_of_nonneg_right hk_lower' (by positivity)
  have hU_le_s_nat : U ≤ s := by
    have : 4 * U ≤ s := by simpa [U, m₀, k, s] using hsmall
    omega
  have hU_le_s : (U : ℝ) ≤ S := by
    have hcast : (U : ℝ) ≤ (s : ℝ) := by exact_mod_cast hU_le_s_nat
    simpa [S] using hcast
  have hs_le_W : S ≤ 2048 * W := by
    calc
      S ≤ 2048 * (Real.log R / Real.log (Real.log R)) := by
        simpa [s, S] using hs.2
      _ = 2048 * (W / LL) := by rfl
      _ ≤ 2048 * W := by
        have hX_le_W : W / LL ≤ W := by
          have hLL_ge_one' : 1 ≤ LL := by simpa [LL, W] using hLL_ge_one
          rw [div_le_iff₀ hLL_pos]
          calc
            W = W * 1 := by ring
            _ ≤ W * LL := mul_le_mul_of_nonneg_left hLL_ge_one' hW_pos.le
        nlinarith
  have hU_pos_nat : 0 < U := by
    dsimp [U]
    unfold IntervalStack.geomLower
    exact Nat.mul_pos (by omega) (Nat.pow_pos (by norm_num))
  have hU_pos : 0 < (U : ℝ) := by exact_mod_cast hU_pos_nat
  have hU_one : 1 ≤ (U : ℝ) := by exact_mod_cast Nat.succ_le_of_lt hU_pos_nat
  have h2048W_le_R : 2048 * W ≤ R := by
    have hW_sq : W ^ 2 ≤ (1 / 2048 : ℝ) * R := by
      simpa [W] using hlog_sq
    have hW_le_sq : W ≤ W ^ 2 := by
      have hW_one : 1 ≤ W := le_of_lt hW_gt_one
      calc
        W = W * 1 := by ring
        _ ≤ W * W := mul_le_mul_of_nonneg_left hW_one hW_pos.le
        _ = W ^ 2 := by ring
    calc
      2048 * W ≤ 2048 * W ^ 2 := by nlinarith
      _ ≤ 2048 * ((1 / 2048 : ℝ) * R) :=
        mul_le_mul_of_nonneg_left hW_sq (by norm_num)
      _ = R := by ring
  have hU_le_R : (U : ℝ) ≤ R := (hU_le_s.trans hs_le_W).trans h2048W_le_R
  have hlogU_nonneg : 0 ≤ logU := by
    dsimp [logU]
    exact Real.log_nonneg hU_one
  have hlogU_le_W : logU ≤ W := by
    dsimp [logU, W]
    exact Real.log_le_log hU_pos hU_le_R
  have hUplus_le : (U : ℝ) + 1 ≤ 2 * S := by
    have hS_nonneg : 0 ≤ S := le_trans (by norm_num) hS_one
    have hOne_le_S : 1 ≤ S := hS_one
    nlinarith
  have hUterm : 2 * (((U : ℝ) + 1) * logU) ≤ 4 * S * W := by
    calc
      2 * (((U : ℝ) + 1) * logU) ≤ 2 * ((2 * S) * W) := by
        have hmul := mul_le_mul hUplus_le hlogU_le_W hlogU_nonneg (by nlinarith)
        nlinarith
      _ = 4 * S * W := by ring
  have hthreshold :=
    crude_threshold_algebra
      (S := S) (W := W) (LL := LL) (K := K) (U := (U : ℝ)) (logU := logU)
      (c := Real.log C) hS_one hW_pos hLL_pos hSLL hKgain hc hUterm
  have hlog_full :
      Real.log (C ^ 2 * R) = 2 * Real.log C + W := by
    rw [Real.log_mul (pow_ne_zero 2 hC.ne') hR_pos.ne']
    rw [Real.log_pow, show (2 : ℕ) * Real.log C = 2 * Real.log C by norm_num]
  have hcoeff :
      (((s * (2 * s + 1) : ℕ) : ℝ)) = S * (2 * S + 1) := by
    simp [S]
  have hs_sq : (((s * s : ℕ) : ℝ)) = S ^ 2 := by
    simp [S]
    ring
  simpa [s, m₀, k, U, S, W, LL, K, logU, hcoeff, hs_sq, hlog_full] using hthreshold

theorem tendsto_sublogM₀_atTop : Tendsto sublogM₀ atTop atTop := by
  rw [tendsto_atTop]
  intro n
  have hsqrt_atTop :
      Tendsto (fun R : ℝ => Real.sqrt (Real.log R)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp Real.tendsto_log_atTop
  filter_upwards [hsqrt_atTop.eventually_ge_atTop (n : ℝ)] with R hR
  unfold sublogM₀
  exact Nat.le_floor hR

theorem eventually_chebyshev_errors_sublog :
    ∀ᶠ R : ℝ in atTop,
      ∀ i : Fin (sublogK R),
        Real.log (((IntervalStack.geomUpper (sublogM₀ R) i : ℕ) : ℝ) + 1) +
            2 * Real.sqrt ((IntervalStack.geomUpper (sublogM₀ R) i : ℕ) : ℝ) *
              Real.log ((IntervalStack.geomUpper (sublogM₀ R) i : ℕ) : ℝ) ≤
          (Real.log 2 / 4) * ((IntervalStack.geomUpper (sublogM₀ R) i : ℕ) : ℝ) := by
  have hcheb := MertensLower.chebyshev_error_eventually_le_quarter_log_two
  rw [eventually_atTop] at hcheb
  rcases hcheb with ⟨N₀, hN₀⟩
  have hm₀_event : ∀ᶠ R : ℝ in atTop, N₀ ≤ sublogM₀ R :=
    tendsto_sublogM₀_atTop.eventually_ge_atTop N₀
  filter_upwards [hm₀_event] with R hm₀
  exact IntervalStack.chebyshev_error_on_geomUpper_of_base
    (m₀ := sublogM₀ R) (k := sublogK R)
    (fun n hn => hN₀ n (hm₀.trans hn))

/-- Final asymptotic Jarnik bound in the parametrized arc form used by the finite theorem.

For each fixed `C > 0`, once `R` is sufficiently large, every injective ordered family of
Gaussian-integer points on `x^2 + y^2 = N = R^2` whose arclength parameters lie in an interval of
length at most `C * sqrt R` has cardinality at most a fixed multiple of
`log R / log log R`. -/
theorem eventually_jarnik_arc_sublog (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      ∀ {M N : ℕ} {a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ},
        R ^ 2 = (N : ℝ) →
        0 < M →
        0 < N →
        L ≤ C * Real.sqrt R →
        (∀ n, z n * star (z n) = (((N : ℤ) : GaussianInt))) →
        Function.Injective z →
        (∀ i j, i ≤ j → j < M → t i ≤ t j) →
        (∀ i j, i < M → j < M → SubcriticalBound.gaussianSqDist (z i) (z j) ≤
          (t j - t i) ^ 2) →
        (∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) →
        (M : ℝ) ≤ 16384 * (Real.log R / Real.log (Real.log R)) := by
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hscale := eventually_scale_hypotheses C hC
  have hthreshold := eventually_crude_threshold_condition C hC
  have hlower := eventually_lower_endpoint_log_conditions C hC
  have hk := eventually_sublogK_positive_and_lower
  have hm₀ := eventually_two_le_sublogM₀
  have hsmall := eventually_smallPrime_sublog_parameters
  have hs := eventually_sublogS_bounds
  have herr := eventually_chebyshev_errors_sublog
  filter_upwards [hR_gt_one, hscale, hthreshold, hlower, hk, hm₀, hsmall, hs, herr] with
    R hR_gt_one hscale hthreshold hlower hk hm₀ hsmall hs herr
  intro M N a L z t hR2 hM hN hLbound hcircle hz hmono hparam hmem
  let s := sublogS R
  let m₀ := sublogM₀ R
  let k := sublogK R
  have hR : 0 < R := lt_trans zero_lt_one hR_gt_one
  have hfinite :=
    MainSublogBound.card_le_eight_mul_log_div_loglog_of_radius_sq_twoScale_endpoint_bounds
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (C := C) (R := R) (B := (2048 : ℝ)) (a := a) (L := L) (z := z) (t := t)
      hR2
      (by simpa [m₀, k] using hscale.1)
      (by simpa [m₀, k] using hscale.2.1)
      (by simpa [m₀, k] using hscale.2.2)
      (by simpa [s, m₀, k] using hthreshold)
      (by simpa [s, m₀] using hlower.1)
      (by simpa [s, m₀] using hlower.2)
      (by simpa [k] using hk.1)
      (by simpa [m₀] using hm₀)
      hM hC hR hLbound
      (by simpa [s, m₀, k] using hsmall)
      hN
      (by simpa [s] using hs.2)
      (by simpa [m₀, k] using herr)
      hcircle hz hmono hparam hmem
  simpa [show (8 : ℝ) * 2048 = 16384 by norm_num, mul_assoc] using hfinite

end AsymptoticParameters
end GaussianChain
