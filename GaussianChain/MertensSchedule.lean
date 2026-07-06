import GaussianChain.AsymptoticParameters

/-!
# Parameter schedule for the Mertens-based variant of the main theorem

This file introduces the window parameter `mertensS` and the upper interval endpoint
`mertensM₁` for a Mertens-based prime interval `[sublogM₀ R, mertensM₁ R]`, together
with the `∀ᶠ R in atTop` facts consumed by the downstream wiring.  The statements
mirror their counterparts in `GaussianChain.AsymptoticParameters`, with the schedule
constants `9/4` (window) and `49/50` (upper endpoint exponent).
-/

namespace GaussianChain
namespace MertensSchedule

open Filter Asymptotics
open AsymptoticParameters

/-- The window parameter of the Mertens schedule. -/
noncomputable def mertensS (R : ℝ) : ℕ :=
  Nat.ceil ((9 / 4 : ℝ) * (Real.log R / Real.log (Real.log R)))

/-- The upper endpoint of the Mertens prime interval: the power `(log R) ^ (49/50)`. -/
noncomputable def mertensM₁ (R : ℝ) : ℕ :=
  Nat.floor ((Real.log R) ^ ((49 : ℝ) / 50))

theorem mertensS_lower (R : ℝ) :
    (9 / 4 : ℝ) * (Real.log R / Real.log (Real.log R)) ≤ (mertensS R : ℝ) := by
  unfold mertensS
  exact Nat.le_ceil _

theorem mertensM₁_le_rpow {R : ℝ} (hW : 0 ≤ Real.log R) :
    (mertensM₁ R : ℝ) ≤ (Real.log R) ^ ((49 : ℝ) / 50) := by
  unfold mertensM₁
  exact Nat.floor_le (Real.rpow_nonneg hW _)

theorem mertensM₁_le_log {R : ℝ} (hW : 1 ≤ Real.log R) :
    (mertensM₁ R : ℝ) ≤ Real.log R := by
  have h₂ : (Real.log R) ^ ((49 : ℝ) / 50) ≤ (Real.log R) ^ (1 : ℝ) :=
    Real.rpow_le_rpow_of_exponent_le hW (by norm_num)
  rw [Real.rpow_one] at h₂
  exact (mertensM₁_le_rpow (le_trans zero_le_one hW)).trans h₂

/-- Passing to the floor loses at most a factor `2` in the logarithm once the value is `≥ 2`. -/
theorem log_nat_floor_ge {x : ℝ} (hx : 2 ≤ x) :
    Real.log x - Real.log 2 ≤ Real.log (Nat.floor x : ℝ) := by
  have hx_pos : 0 < x := by linarith
  have hhalf_pos : 0 < x / 2 := by linarith
  have hfloor : x / 2 ≤ (Nat.floor x : ℝ) := by
    have h1 : x < (Nat.floor x : ℝ) + 1 := Nat.lt_floor_add_one x
    linarith
  calc
    Real.log x - Real.log 2 = Real.log (x / 2) :=
      (Real.log_div hx_pos.ne' (by norm_num)).symm
    _ ≤ Real.log (Nat.floor x : ℝ) := Real.log_le_log hhalf_pos hfloor

theorem eventually_mertensS_lower :
    ∀ᶠ R : ℝ in atTop,
      (9 / 4 : ℝ) * (Real.log R / Real.log (Real.log R)) ≤ (mertensS R : ℝ) :=
  Eventually.of_forall mertensS_lower

theorem eventually_mertensS_pos : ∀ᶠ R : ℝ in atTop, 0 < mertensS R := by
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  filter_upwards [hlog_gt_one, hloglog_pos] with R hW_gt_one hLL_pos
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  have hpos : 0 < (9 / 4 : ℝ) * (Real.log R / Real.log (Real.log R)) := by positivity
  unfold mertensS
  exact Nat.ceil_pos.mpr hpos

theorem eventually_mertensS_upper_of_gt {Q : ℝ} (hQ : (9 / 4 : ℝ) < Q) :
    ∀ᶠ R : ℝ in atTop,
      (mertensS R : ℝ) ≤ Q * (Real.log R / Real.log (Real.log R)) := by
  have hδ : 0 < Q - (9 / 4 : ℝ) := by linarith
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (Q - (9 / 4 : ℝ)) * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1) hδ
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hloglog_pos, hLL_small, hlog_gt_one] with R hLL hsmall hW_gt_one
  let W := Real.log R
  let LL := Real.log W
  let X := W / LL
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hsmall' : LL ≤ (Q - (9 / 4 : ℝ)) * W := by
    simpa [LL, W, Real.rpow_one] using hsmall
  have hone : 1 ≤ (Q - (9 / 4 : ℝ)) * X := by
    rw [show (Q - (9 / 4 : ℝ)) * X =
        ((Q - (9 / 4 : ℝ)) * W) / LL by ring]
    rw [le_div_iff₀ (by simpa [LL, W] using hLL)]
    simpa using hsmall'
  have hX_nonneg : 0 ≤ (9 / 4 : ℝ) * X := by positivity
  have hceil_lt : (mertensS R : ℝ) < (9 / 4 : ℝ) * X + 1 := by
    unfold mertensS
    exact Nat.ceil_lt_add_one hX_nonneg
  have htarget : (9 / 4 : ℝ) * X + 1 ≤ Q * X := by
    nlinarith
  exact le_of_lt (hceil_lt.trans_le htarget)

theorem eventually_two_le_mertensM₁ :
    ∀ᶠ R : ℝ in atTop, 2 ≤ mertensM₁ R := by
  have hpow_atTop :
      Tendsto (fun R : ℝ => (Real.log R) ^ ((49 : ℝ) / 50)) atTop atTop :=
    (tendsto_rpow_atTop (by norm_num : (0 : ℝ) < 49 / 50)).comp Real.tendsto_log_atTop
  filter_upwards [hpow_atTop.eventually_ge_atTop (2 : ℝ)] with R hR
  unfold mertensM₁
  exact Nat.le_floor hR

theorem eventually_four_mul_mertensM₁_le_mertensS :
    ∀ᶠ R : ℝ in atTop, 4 * mertensM₁ R ≤ mertensS R := by
  have hLL_bound :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (9 / 16 : ℝ) * (Real.log R) ^ (1 / 50 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1 / 50)
      (by norm_num : (0 : ℝ) < 9 / 16)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  filter_upwards [hLL_bound, hlog_gt_one, hloglog_pos] with R hLL hW_gt_one hLL_pos
  let W := Real.log R
  have hW_pos : 0 < W := by
    dsimp [W]
    exact lt_trans zero_lt_one hW_gt_one
  have hm₁_le : (mertensM₁ R : ℝ) ≤ W ^ ((49 : ℝ) / 50) :=
    mertensM₁_le_rpow hW_pos.le
  have hpowprod : W ^ ((49 : ℝ) / 50) * W ^ (1 / 50 : ℝ) = W := by
    have h1 : (49 : ℝ) / 50 + 1 / 50 = 1 := by norm_num
    rw [← Real.rpow_add hW_pos, h1, Real.rpow_one]
  have hscale :
      4 * (W ^ ((49 : ℝ) / 50)) ≤ (9 / 4 : ℝ) * (W / Real.log W) := by
    rw [show (9 / 4 : ℝ) * (W / Real.log W) =
        ((9 / 4 : ℝ) * W) / Real.log W by ring]
    rw [le_div_iff₀ (by simpa [W] using hLL_pos)]
    calc
      4 * W ^ ((49 : ℝ) / 50) * Real.log W
          ≤ 4 * W ^ ((49 : ℝ) / 50) * ((9 / 16 : ℝ) * W ^ (1 / 50 : ℝ)) := by
            exact mul_le_mul_of_nonneg_left (by simpa [W] using hLL) (by positivity)
      _ = (9 / 4 : ℝ) * (W ^ ((49 : ℝ) / 50) * W ^ (1 / 50 : ℝ)) := by ring
      _ = (9 / 4 : ℝ) * W := by rw [hpowprod]
  have hfour : ((4 * mertensM₁ R : ℕ) : ℝ) ≤ (mertensS R : ℝ) := by
    have hreal : (4 : ℝ) * (mertensM₁ R : ℝ) ≤ (mertensS R : ℝ) := by
      calc
        (4 : ℝ) * (mertensM₁ R : ℝ) ≤ 4 * W ^ ((49 : ℝ) / 50) :=
          mul_le_mul_of_nonneg_left hm₁_le (by norm_num)
        _ ≤ (9 / 4 : ℝ) * (W / Real.log W) := hscale
        _ = (9 / 4 : ℝ) * (Real.log R / Real.log (Real.log R)) := by rfl
        _ ≤ (mertensS R : ℝ) := mertensS_lower R
    simpa [Nat.cast_mul] using hreal
  exact_mod_cast hfour

theorem eventually_sublogM₀_le_mertensM₁ :
    ∀ᶠ R : ℝ in atTop, sublogM₀ R ≤ mertensM₁ R := by
  have hlog_ge_one : ∀ᶠ R : ℝ in atTop, 1 ≤ Real.log R :=
    Real.tendsto_log_atTop.eventually_ge_atTop 1
  filter_upwards [hlog_ge_one] with R hW_one
  have hsqrt_le : Real.sqrt (Real.log R) ≤ (Real.log R) ^ ((49 : ℝ) / 50) := by
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le hW_one (by norm_num)
  unfold AsymptoticParameters.sublogM₀ mertensM₁
  exact Nat.floor_mono hsqrt_le

theorem eventually_log_mertensM₁_ge :
    ∀ᶠ R : ℝ in atTop,
      (49 / 50 : ℝ) * Real.log (Real.log R) - Real.log 2 ≤
        Real.log (mertensM₁ R : ℝ) := by
  have hpow_atTop :
      Tendsto (fun R : ℝ => (Real.log R) ^ ((49 : ℝ) / 50)) atTop atTop :=
    (tendsto_rpow_atTop (by norm_num : (0 : ℝ) < 49 / 50)).comp Real.tendsto_log_atTop
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hpow_atTop.eventually_ge_atTop (2 : ℝ), hlog_gt_one] with
    R hpow_two hW_gt_one
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  have hfloor := log_nat_floor_ge hpow_two
  rw [Real.log_rpow hW_pos] at hfloor
  unfold mertensM₁
  linarith

theorem eventually_log_sublogM₀_le_half :
    ∀ᶠ R : ℝ in atTop,
      Real.log (sublogM₀ R : ℝ) ≤ (1 / 2 : ℝ) * Real.log (Real.log R) := by
  have hm₀_two := eventually_two_le_sublogM₀
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hm₀_two, hlog_gt_one] with R hm₀ hW_gt_one
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  have hm₀_pos : 0 < (sublogM₀ R : ℝ) := by
    have h : 0 < sublogM₀ R := by omega
    exact_mod_cast h
  have hle : Real.log (sublogM₀ R : ℝ) ≤ Real.log (Real.sqrt (Real.log R)) :=
    Real.log_le_log hm₀_pos sublogM₀_le_sqrt_log
  rw [Real.log_sqrt hW_pos.le] at hle
  linarith

theorem eventually_log_sublogM₀_ge :
    ∀ᶠ R : ℝ in atTop,
      (1 / 2 : ℝ) * Real.log (Real.log R) - Real.log 2 ≤
        Real.log (sublogM₀ R : ℝ) := by
  have hsqrt_atTop :
      Tendsto (fun R : ℝ => Real.sqrt (Real.log R)) atTop atTop :=
    Real.tendsto_sqrt_atTop.comp Real.tendsto_log_atTop
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  filter_upwards [hsqrt_atTop.eventually_ge_atTop (2 : ℝ), hlog_gt_one] with
    R hsqrt_two hW_gt_one
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  have hfloor := log_nat_floor_ge hsqrt_two
  rw [Real.log_sqrt hW_pos.le] at hfloor
  unfold AsymptoticParameters.sublogM₀
  linarith

theorem eventually_mertens_scale_hypotheses (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      1 ≤ C ^ 2 * R ∧
        1 ≤ C ^ 2 * R / (mertensM₁ R : ℝ) ∧
          1 ≤ C ^ 2 * R / (mertensM₁ R : ℝ) ^ 2 := by
  have hm₁_two := eventually_two_le_mertensM₁
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_large : ∀ᶠ R : ℝ in atTop, 1 ≤ C ^ 2 * R := by
    have htendsto : Tendsto (fun R : ℝ => C ^ 2 * R) atTop atTop :=
      Tendsto.const_mul_atTop (sq_pos_of_pos hC) tendsto_id
    exact htendsto.eventually_ge_atTop 1
  have hlog_sq_small :
      ∀ᶠ R : ℝ in atTop, Real.log R ^ 2 ≤ C ^ 2 * R :=
    eventually_log_sq_le_const_mul_self (sq_pos_of_pos hC)
  filter_upwards [hm₁_two, hlog_gt_one, hR_large, hlog_sq_small] with
    R hm₁_two hW_gt_one hfull hlog_sq
  have hW_pos : 0 < Real.log R := lt_trans zero_lt_one hW_gt_one
  have hm₁_pos : 0 < (mertensM₁ R : ℝ) := by
    have h : 0 < mertensM₁ R := by omega
    exact_mod_cast h
  have hm₁_le_W : (mertensM₁ R : ℝ) ≤ Real.log R := mertensM₁_le_log hW_gt_one.le
  have hW_le_sq : Real.log R ≤ Real.log R ^ 2 := by nlinarith
  have hm₁_le_CR : (mertensM₁ R : ℝ) ≤ C ^ 2 * R :=
    hm₁_le_W.trans (hW_le_sq.trans hlog_sq)
  have hm₁_sq_le_CR : (mertensM₁ R : ℝ) ^ 2 ≤ C ^ 2 * R := by
    have h : (mertensM₁ R : ℝ) ^ 2 ≤ Real.log R ^ 2 := by nlinarith
    exact h.trans hlog_sq
  exact ⟨hfull, (one_le_div hm₁_pos).mpr hm₁_le_CR,
    (one_le_div (sq_pos_of_pos hm₁_pos)).mpr hm₁_sq_le_CR⟩

theorem mertens_split_lower_endpoint_algebra
    {S W LL logm c : ℝ}
    (hS : 1 ≤ S)
    (_hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : (9 / 4 : ℝ) * W ≤ S * LL)
    (hlogm : LL / 2 - Real.log 2 ≤ logm)
    (hlog2 : Real.log 2 ≤ LL / 72)
    (hc : 2 * c ≤ LL / 1152) :
    (S * (2 * S + 1)) * (2 * c + W - logm) <
      S ^ 2 * (2 * W - logm) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hlogm' : (35 / 72 : ℝ) * LL ≤ logm := by linarith
  have hlogm_pos : 0 < logm := by nlinarith
  have hconst :
      (S * (2 * S + 1)) * (2 * c) ≤ (3 / 1152) * S ^ 2 * LL := by
    have hcoeff_nonneg : 0 ≤ S * (2 * S + 1) := by nlinarith
    have hcoeff_le : S * (2 * S + 1) ≤ 3 * S ^ 2 := by nlinarith
    calc
      (S * (2 * S + 1)) * (2 * c) ≤
          (S * (2 * S + 1)) * (LL / 1152) :=
        mul_le_mul_of_nonneg_left hc hcoeff_nonneg
      _ ≤ (3 * S ^ 2) * (LL / 1152) := by
        exact mul_le_mul_of_nonneg_right hcoeff_le (by nlinarith)
      _ = (3 / 1152) * S ^ 2 * LL := by ring
  have hloggain : S ^ 2 * ((35 / 72 : ℝ) * LL) ≤ S ^ 2 * logm :=
    mul_le_mul_of_nonneg_left hlogm' (sq_nonneg S)
  have hSW_small : S * W ≤ (4 / 9 : ℝ) * S ^ 2 * LL := by
    have hmul : (9 / 4 : ℝ) * S * W ≤ S ^ 2 * LL := by nlinarith
    nlinarith
  nlinarith [mul_pos (mul_pos hSpos hSpos) hLL, mul_pos hSpos hlogm_pos]

theorem mertens_inert_lower_endpoint_algebra
    {S W LL logm c : ℝ}
    (hS : 1 ≤ S)
    (hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : (9 / 4 : ℝ) * W ≤ S * LL)
    (hlogm : LL / 2 - Real.log 2 ≤ logm)
    (hlog2 : Real.log 2 ≤ LL / 72)
    (hc : 2 * c ≤ LL / 1152) :
    (S * (2 * S + 1)) * (2 * c + W - 2 * logm) <
      S ^ 2 * (2 * W - 2 * logm) := by
  have hsplit :=
    mertens_split_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := c)
      hS hW hLL hSLL hlogm hlog2 hc
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hlogm_pos : 0 < logm := by nlinarith
  have hextra_pos :
      0 < (S * (2 * S + 1)) * logm - S ^ 2 * logm := by
    have hcoeff_gap : 0 < S ^ 2 + S := by nlinarith
    nlinarith
  nlinarith

theorem mertens_threshold_algebra
    {S W LL G U logU c : ℝ}
    (hS : 1 ≤ S)
    (hW : 0 < W)
    (hLL : 0 < LL)
    (hSLL : (9 / 4 : ℝ) * W ≤ S * LL)
    (hGgain : (5 / 11 : ℝ) * LL ≤ G)
    (hc : 2 * c ≤ LL / 1152)
    (hUterm : 2 * ((U + 1) * logU) ≤ (1 / 1000) * S ^ 2 * LL) :
    (S * (2 * S + 1)) * (2 * c + W) <
      (S ^ 2 * G - 2 * ((U + 1) * logU)) + S ^ 2 * (2 * W) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hstack : S ^ 2 * ((5 / 11 : ℝ) * LL) ≤ S ^ 2 * G :=
    mul_le_mul_of_nonneg_left hGgain (sq_nonneg S)
  have hconst :
      (S * (2 * S + 1)) * (2 * c) ≤ (3 / 1152) * S ^ 2 * LL := by
    have hcoeff_nonneg : 0 ≤ S * (2 * S + 1) := by nlinarith
    have hcoeff_le : S * (2 * S + 1) ≤ 3 * S ^ 2 := by nlinarith
    calc
      (S * (2 * S + 1)) * (2 * c) ≤
          (S * (2 * S + 1)) * (LL / 1152) :=
        mul_le_mul_of_nonneg_left hc hcoeff_nonneg
      _ ≤ (3 * S ^ 2) * (LL / 1152) := by
        exact mul_le_mul_of_nonneg_right hcoeff_le (by nlinarith)
      _ = (3 / 1152) * S ^ 2 * LL := by ring
  have hSW_small : S * W ≤ (4 / 9 : ℝ) * S ^ 2 * LL := by
    have hmul : (9 / 4 : ℝ) * S * W ≤ S ^ 2 * LL := by nlinarith
    nlinarith
  have herr_small :
      S * W + 2 * ((U + 1) * logU) + (S * (2 * S + 1)) * (2 * c) ≤
        (4 / 9 + 1 / 1000 + 3 / 1152) * S ^ 2 * LL := by
    nlinarith
  have hgap : (4 / 9 + 1 / 1000 + 3 / 1152 : ℝ) * S ^ 2 * LL <
      S ^ 2 * ((5 / 11 : ℝ) * LL) := by
    have hSsq_pos : 0 < S ^ 2 := sq_pos_of_pos hSpos
    nlinarith
  nlinarith

theorem eventually_mertens_lower_endpoint_log_conditions (C : ℝ) (hC : 0 < C) :
    ∀ᶠ R : ℝ in atTop,
      let s := mertensS R
      let m₀ := sublogM₀ R
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)) ∧
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)) := by
  have hlogm_lower := eventually_log_sublogM₀_ge
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
      ∀ᶠ R : ℝ in atTop, 2304 * Real.log C ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (2304 * Real.log C)
  have hlog2_small :
      ∀ᶠ R : ℝ in atTop, 72 * Real.log 2 ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (72 * Real.log 2)
  filter_upwards [hlogm_lower, hm₀_two, hLL_small, hlog_gt_one, hR_gt_one, hconst,
    hlog2_small] with R hlogm_lower hm₀_two hLL_small hW_gt_one hR_gt_one hconst hlog2_small
  let s := mertensS R
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
  have hs_lower : (9 / 4 : ℝ) * (W / LL) ≤ S := by
    simpa [s, S, W, LL] using mertensS_lower R
  have hS_one : 1 ≤ S := by nlinarith
  have hSLL : (9 / 4 : ℝ) * W ≤ S * LL := by
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : (9 / 4 : ℝ) * (W / LL) * LL = (9 / 4 : ℝ) * W := by
      field_simp [hLL_pos.ne']
    nlinarith
  have hlogm : LL / 2 - Real.log 2 ≤ logm := by
    have hraw : (1 / 2 : ℝ) * LL - Real.log 2 ≤ logm := by
      simpa [m₀, LL, W, logm] using hlogm_lower
    linarith
  have hlog2 : Real.log 2 ≤ LL / 72 := by
    dsimp [LL, W]
    linarith
  have hc : 2 * Real.log C ≤ LL / 1152 := by
    dsimp [LL, W]
    linarith
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
    mertens_split_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C)
      hS_one hW_pos hLL_pos hSLL hlogm hlog2 hc
  have hinert :=
    mertens_inert_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C)
      hS_one hW_pos hLL_pos hSLL hlogm hlog2 hc
  constructor
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_split_scale] using hsplit
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_inert_scale] using hinert

set_option maxHeartbeats 400000 in
-- The chain of `nlinarith` calls and the final cast-normalizing `simpa` exceed the
-- default heartbeat budget, as in `eventually_crude_threshold_condition`.
theorem eventually_mertens_threshold_condition
    (C c : ℝ) (hC : 0 < C) (hc : c ≤ 10) :
    ∀ᶠ R : ℝ in atTop,
      let s := mertensS R
      let m₀ := sublogM₀ R
      let m₁ := mertensM₁ R
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c) -
            2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R) := by
  have hfour := eventually_four_mul_mertensM₁_le_mertensS
  have hm₁_two := eventually_two_le_mertensM₁
  have hlogm₁ := eventually_log_mertensM₁_ge
  have hlogm₀ := eventually_log_sublogM₀_le_half
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (1 / 2048) * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul (by norm_num : (0 : ℝ) < 1)
      (by norm_num : (0 : ℝ) < 1 / 2048)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hconst :
      ∀ᶠ R : ℝ in atTop, 2304 * Real.log C ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop
      (2304 * Real.log C)
  have hLL_large :
      ∀ᶠ R : ℝ in atTop, (433 : ℝ) ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop 433
  filter_upwards [hfour, hm₁_two, hlogm₁, hlogm₀, hLL_small, hlog_gt_one, hR_gt_one,
    hconst, hLL_large] with
    R hfour hm₁_two hlogm₁ hlogm₀ hLL_small hW_gt_one hR_gt_one hconst hLL_large
  let s := mertensS R
  let m₀ := sublogM₀ R
  let m₁ := mertensM₁ R
  let S : ℝ := s
  let W := Real.log R
  let LL := Real.log W
  let logU := Real.log (m₁ : ℝ)
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
  have hs_lower : (9 / 4 : ℝ) * (W / LL) ≤ S := by
    simpa [s, S, W, LL] using mertensS_lower R
  have hS_large : 4608 ≤ S := by nlinarith
  have hS_one : 1 ≤ S := by linarith
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS_one
  have hSLL : (9 / 4 : ℝ) * W ≤ S * LL := by
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : (9 / 4 : ℝ) * (W / LL) * LL = (9 / 4 : ℝ) * W := by
      field_simp [hLL_pos.ne']
    nlinarith
  have hcLL : 2 * Real.log C ≤ LL / 1152 := by
    dsimp [LL, W]
    linarith
  have hGgain :
      (5 / 11 : ℝ) * LL ≤ Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c := by
    have hm₁_log : (49 / 50 : ℝ) * LL - Real.log 2 ≤ Real.log (m₁ : ℝ) := by
      simpa [m₁, LL, W] using hlogm₁
    have hm₀_log : Real.log (m₀ : ℝ) ≤ (1 / 2 : ℝ) * LL := by
      simpa [m₀, LL, W] using hlogm₀
    have hLL_large' : (433 : ℝ) ≤ LL := by
      dsimp [LL, W]
      exact hLL_large
    have hlog2_le_one : Real.log 2 ≤ 1 := by
      have h := Real.log_le_sub_one_of_pos (by norm_num : (0 : ℝ) < 2)
      linarith
    linarith
  have hU_two : (2 : ℝ) ≤ (m₁ : ℝ) := by
    have h : 2 ≤ m₁ := hm₁_two
    exact_mod_cast h
  have hU_pos : 0 < (m₁ : ℝ) := by linarith
  have hlogU_nonneg : 0 ≤ logU := Real.log_nonneg (by linarith)
  have hm₁_le_W : (m₁ : ℝ) ≤ W := mertensM₁_le_log hW_gt_one.le
  have hlogU_le_LL : logU ≤ LL := Real.log_le_log hU_pos hm₁_le_W
  have hfourU_le_S : (4 : ℝ) * (m₁ : ℝ) ≤ S := by
    have hcast : ((4 * m₁ : ℕ) : ℝ) ≤ (s : ℝ) := by
      exact_mod_cast (hfour : 4 * mertensM₁ R ≤ mertensS R)
    simpa [Nat.cast_mul, S] using hcast
  have hUplus_le : (m₁ : ℝ) + 1 ≤ S / 2 := by linarith
  have hUterm :
      2 * (((m₁ : ℝ) + 1) * logU) ≤ (1 / 1000) * S ^ 2 * LL := by
    have h1 : 2 * (((m₁ : ℝ) + 1) * logU) ≤ 2 * ((S / 2) * LL) := by
      have hhalf_nonneg : 0 ≤ S / 2 := by linarith
      have hmul := mul_le_mul hUplus_le hlogU_le_LL hlogU_nonneg hhalf_nonneg
      linarith
    have h2 : 2 * ((S / 2) * LL) = S * LL := by ring
    have h3 : S * LL ≤ (1 / 1000) * S ^ 2 * LL := by
      nlinarith [mul_nonneg (mul_nonneg
        (by linarith : (0 : ℝ) ≤ S - 1000) hSpos.le) hLL_pos.le]
    linarith
  have hthreshold :=
    mertens_threshold_algebra
      (S := S) (W := W) (LL := LL)
      (G := Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c)
      (U := (m₁ : ℝ)) (logU := logU) (c := Real.log C)
      hS_one hW_pos hLL_pos hSLL hGgain hcLL hUterm
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
  simpa [s, m₀, m₁, S, W, LL, logU, hcoeff, hs_sq, hlog_full] using hthreshold

end MertensSchedule
end GaussianChain
