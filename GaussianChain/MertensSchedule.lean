import GaussianChain.AsymptoticParameters

/-!
# Parameter schedule for the Mertens-based variant of the main theorem

This file introduces the window parameter `mertensS q` and the upper interval endpoint
`mertensM₁ θ₁` for a Mertens-based prime interval `[sublogM₀ R, mertensM₁ θ₁ R]`, together
with the `∀ᶠ R in atTop` facts consumed by the downstream wiring.  The schedule is fully
parametric, matching the paper: the window constant `q` ranges over `(2, ∞)` and the upper
endpoint exponent `θ₁` over `(1/2 + 1/q, 1)`.  Every nuisance term (rounding losses, the
fixed constants `log C` and the Mertens deficit, and the `(m₁ + 1) log m₁` toll) is
`o(log log R)` against the positive gaps `1/2 - 1/q` and `θ₁ - 1/2 - 1/q`, so each
hypothesis of the finite dichotomy holds for all sufficiently large radii.
-/

namespace GaussianChain
namespace MertensSchedule

open Filter Asymptotics
open AsymptoticParameters

/-- The window parameter of the Mertens schedule with window constant `q`. -/
noncomputable def mertensS (q R : ℝ) : ℕ :=
  Nat.ceil (q * (Real.log R / Real.log (Real.log R)))

/-- The upper endpoint of the Mertens prime interval: the power `(log R) ^ θ₁`. -/
noncomputable def mertensM₁ (θ₁ R : ℝ) : ℕ :=
  Nat.floor ((Real.log R) ^ θ₁)

theorem mertensS_lower (q R : ℝ) :
    q * (Real.log R / Real.log (Real.log R)) ≤ (mertensS q R : ℝ) := by
  unfold mertensS
  exact Nat.le_ceil _

theorem mertensM₁_le_rpow {θ₁ R : ℝ} (hW : 0 ≤ Real.log R) :
    (mertensM₁ θ₁ R : ℝ) ≤ (Real.log R) ^ θ₁ := by
  unfold mertensM₁
  exact Nat.floor_le (Real.rpow_nonneg hW _)

theorem mertensM₁_le_log {θ₁ R : ℝ} (hθ : θ₁ ≤ 1) (hW : 1 ≤ Real.log R) :
    (mertensM₁ θ₁ R : ℝ) ≤ Real.log R := by
  have h₂ : (Real.log R) ^ θ₁ ≤ (Real.log R) ^ (1 : ℝ) :=
    Real.rpow_le_rpow_of_exponent_le hW hθ
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

theorem eventually_mertensS_lower {q : ℝ} :
    ∀ᶠ R : ℝ in atTop,
      q * (Real.log R / Real.log (Real.log R)) ≤ (mertensS q R : ℝ) :=
  Eventually.of_forall (mertensS_lower q)

/-- The ratio `log R / log log R` eventually dominates every fixed constant. -/
theorem eventually_le_log_div_loglog (K : ℝ) :
    ∀ᶠ R : ℝ in atTop, K ≤ Real.log R / Real.log (Real.log R) := by
  have hKpos : 0 < max K 1 := lt_of_lt_of_le zero_lt_one (le_max_right K 1)
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (max K 1)⁻¹ * (Real.log R) ^ (1 : ℝ) :=
    eventually_log_log_le_log_rpow_mul one_pos (inv_pos.mpr hKpos)
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  filter_upwards [hLL_small, hloglog_pos] with R hsmall hLL
  rw [Real.rpow_one] at hsmall
  have hmax : max K 1 ≤ Real.log R / Real.log (Real.log R) := by
    rw [le_div_iff₀ hLL]
    calc
      max K 1 * Real.log (Real.log R) ≤ max K 1 * ((max K 1)⁻¹ * Real.log R) :=
        mul_le_mul_of_nonneg_left hsmall hKpos.le
      _ = Real.log R := by
        rw [← mul_assoc, mul_inv_cancel₀ (ne_of_gt hKpos), one_mul]
  exact (le_max_left K 1).trans hmax

/-- Any fixed constant is eventually at most any fixed positive multiple of `log log R`. -/
theorem eventually_const_le_mul_loglog {ε : ℝ} (hε : 0 < ε) (x : ℝ) :
    ∀ᶠ R : ℝ in atTop, x ≤ ε * Real.log (Real.log R) := by
  have hεne : ε ≠ 0 := ne_of_gt hε
  have hlarge :
      ∀ᶠ R : ℝ in atTop, x / ε ≤ Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_ge_atTop (x / ε)
  filter_upwards [hlarge] with R hR
  have hmul : ε * (x / ε) ≤ ε * Real.log (Real.log R) :=
    mul_le_mul_of_nonneg_left hR hε.le
  have hcancel : ε * (x / ε) = x := by field_simp
  linarith

theorem eventually_mertensS_upper_of_gt {q Q : ℝ} (hq : 0 < q) (hQ : q < Q) :
    ∀ᶠ R : ℝ in atTop,
      (mertensS q R : ℝ) ≤ Q * (Real.log R / Real.log (Real.log R)) := by
  have hδ : 0 < Q - q := by linarith
  have hloglog_pos :
      ∀ᶠ R : ℝ in atTop, 0 < Real.log (Real.log R) :=
    (Real.tendsto_log_atTop.comp Real.tendsto_log_atTop).eventually_gt_atTop 0
  have hLL_small :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (Q - q) * (Real.log R) ^ (1 : ℝ) :=
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
  have hLL_pos : 0 < LL := by
    dsimp [LL, W]
    exact hLL
  have hsmall' : LL ≤ (Q - q) * W := by
    simpa [LL, W, Real.rpow_one] using hsmall
  have hone : 1 ≤ (Q - q) * X := by
    rw [show (Q - q) * X = ((Q - q) * W) / LL by ring]
    rw [le_div_iff₀ hLL_pos]
    simpa using hsmall'
  have hX_nonneg : 0 ≤ q * X :=
    mul_nonneg hq.le (le_of_lt (div_pos hW_pos hLL_pos))
  have hceil_lt : (mertensS q R : ℝ) < q * X + 1 := by
    unfold mertensS
    exact Nat.ceil_lt_add_one hX_nonneg
  have htarget : q * X + 1 ≤ Q * X := by nlinarith
  exact le_of_lt (hceil_lt.trans_le htarget)

theorem eventually_two_le_mertensM₁ {θ₁ : ℝ} (hθ : 0 < θ₁) :
    ∀ᶠ R : ℝ in atTop, 2 ≤ mertensM₁ θ₁ R := by
  have hpow_atTop :
      Tendsto (fun R : ℝ => (Real.log R) ^ θ₁) atTop atTop :=
    (tendsto_rpow_atTop hθ).comp Real.tendsto_log_atTop
  filter_upwards [hpow_atTop.eventually_ge_atTop (2 : ℝ)] with R hR
  unfold mertensM₁
  exact Nat.le_floor hR

theorem eventually_four_mul_mertensM₁_le_mertensS {q θ₁ : ℝ}
    (hq : 0 < q) (hθ1 : θ₁ < 1) :
    ∀ᶠ R : ℝ in atTop, 4 * mertensM₁ θ₁ R ≤ mertensS q R := by
  have hexp : 0 < 1 - θ₁ := by linarith
  have hLL_bound :
      ∀ᶠ R : ℝ in atTop,
        Real.log (Real.log R) ≤ (q / 16) * (Real.log R) ^ (1 - θ₁) :=
    eventually_log_log_le_log_rpow_mul hexp (by linarith)
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
  have hm₁_le : (mertensM₁ θ₁ R : ℝ) ≤ W ^ θ₁ :=
    mertensM₁_le_rpow hW_pos.le
  have hpowprod : W ^ θ₁ * W ^ (1 - θ₁) = W := by
    rw [← Real.rpow_add hW_pos, show θ₁ + (1 - θ₁) = 1 by ring, Real.rpow_one]
  have hpow_nonneg : 0 ≤ W ^ θ₁ := Real.rpow_nonneg hW_pos.le _
  have hscale : 4 * (W ^ θ₁) ≤ q * (W / Real.log W) := by
    rw [show q * (W / Real.log W) = (q * W) / Real.log W by ring]
    rw [le_div_iff₀ (by simpa [W] using hLL_pos)]
    calc
      4 * W ^ θ₁ * Real.log W
          ≤ 4 * W ^ θ₁ * ((q / 16) * W ^ (1 - θ₁)) := by
            exact mul_le_mul_of_nonneg_left (by simpa [W] using hLL)
              (by linarith [hpow_nonneg])
      _ = (q / 4) * (W ^ θ₁ * W ^ (1 - θ₁)) := by ring
      _ = (q / 4) * W := by rw [hpowprod]
      _ ≤ q * W := by nlinarith [mul_pos hq hW_pos]
  have hfour : ((4 * mertensM₁ θ₁ R : ℕ) : ℝ) ≤ (mertensS q R : ℝ) := by
    have hreal : (4 : ℝ) * (mertensM₁ θ₁ R : ℝ) ≤ (mertensS q R : ℝ) := by
      calc
        (4 : ℝ) * (mertensM₁ θ₁ R : ℝ) ≤ 4 * W ^ θ₁ :=
          mul_le_mul_of_nonneg_left hm₁_le (by norm_num)
        _ ≤ q * (W / Real.log W) := hscale
        _ = q * (Real.log R / Real.log (Real.log R)) := by rfl
        _ ≤ (mertensS q R : ℝ) := mertensS_lower q R
    simpa [Nat.cast_mul] using hreal
  exact_mod_cast hfour

theorem eventually_sublogM₀_le_mertensM₁ {θ₁ : ℝ} (hθ : 1 / 2 < θ₁) :
    ∀ᶠ R : ℝ in atTop, sublogM₀ R ≤ mertensM₁ θ₁ R := by
  have hlog_ge_one : ∀ᶠ R : ℝ in atTop, 1 ≤ Real.log R :=
    Real.tendsto_log_atTop.eventually_ge_atTop 1
  filter_upwards [hlog_ge_one] with R hW_one
  have hsqrt_le : Real.sqrt (Real.log R) ≤ (Real.log R) ^ θ₁ := by
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le hW_one hθ.le
  unfold AsymptoticParameters.sublogM₀ mertensM₁
  exact Nat.floor_mono hsqrt_le

theorem eventually_log_mertensM₁_ge {θ₁ : ℝ} (hθ : 0 < θ₁) :
    ∀ᶠ R : ℝ in atTop,
      θ₁ * Real.log (Real.log R) - Real.log 2 ≤
        Real.log (mertensM₁ θ₁ R : ℝ) := by
  have hpow_atTop :
      Tendsto (fun R : ℝ => (Real.log R) ^ θ₁) atTop atTop :=
    (tendsto_rpow_atTop hθ).comp Real.tendsto_log_atTop
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

theorem eventually_mertens_scale_hypotheses (C θ₁ : ℝ)
    (hC : 0 < C) (hθ0 : 0 < θ₁) (hθ1 : θ₁ < 1) :
    ∀ᶠ R : ℝ in atTop,
      1 ≤ C ^ 2 * R ∧
        1 ≤ C ^ 2 * R / (mertensM₁ θ₁ R : ℝ) ∧
          1 ≤ C ^ 2 * R / (mertensM₁ θ₁ R : ℝ) ^ 2 := by
  have hm₁_two := eventually_two_le_mertensM₁ hθ0
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
  have hm₁_pos : 0 < (mertensM₁ θ₁ R : ℝ) := by
    have h : 0 < mertensM₁ θ₁ R := by omega
    exact_mod_cast h
  have hm₁_le_W : (mertensM₁ θ₁ R : ℝ) ≤ Real.log R :=
    mertensM₁_le_log hθ1.le hW_gt_one.le
  have hW_le_sq : Real.log R ≤ Real.log R ^ 2 := by nlinarith
  have hm₁_le_CR : (mertensM₁ θ₁ R : ℝ) ≤ C ^ 2 * R :=
    hm₁_le_W.trans (hW_le_sq.trans hlog_sq)
  have hm₁_sq_le_CR : (mertensM₁ θ₁ R : ℝ) ^ 2 ≤ C ^ 2 * R := by
    have h : (mertensM₁ θ₁ R : ℝ) ^ 2 ≤ Real.log R ^ 2 := by nlinarith
    exact h.trans hlog_sq
  exact ⟨hfull, (one_le_div hm₁_pos).mpr hm₁_le_CR,
    (one_le_div (sq_pos_of_pos hm₁_pos)).mpr hm₁_sq_le_CR⟩

theorem mertens_split_lower_endpoint_algebra
    {S W LL logm c δ : ℝ}
    (hS : 1 ≤ S)
    (hLL : 0 < LL)
    (hδ : 0 < δ)
    (hWS : W ≤ S * ((1 / 2 - δ) * LL))
    (hlogm : LL / 2 - Real.log 2 ≤ logm)
    (hlog2 : Real.log 2 ≤ δ / 4 * LL)
    (hc : 2 * c ≤ δ / 8 * LL) :
    (S * (2 * S + 1)) * (2 * c + W - logm) <
      S ^ 2 * (2 * W - logm) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have h1 : (S + 1) * (LL / 2 - Real.log 2) ≤ (S + 1) * logm :=
    mul_le_mul_of_nonneg_left hlogm (by linarith)
  have h2 : (S + 1) * Real.log 2 ≤ (S + 1) * (δ / 4 * LL) :=
    mul_le_mul_of_nonneg_left hlog2 (by linarith)
  have h3 : (2 * S + 1) * (2 * c) ≤ (2 * S + 1) * (δ / 8 * LL) :=
    mul_le_mul_of_nonneg_left hc (by linarith)
  have hgap : 0 < (S + 1) * logm - W - (2 * S + 1) * (2 * c) := by
    nlinarith [h1, h2, h3, hWS, mul_pos hδ hLL,
      mul_nonneg (sub_nonneg.mpr hS) (mul_pos hδ hLL).le]
  nlinarith [mul_pos hSpos hgap]

theorem mertens_inert_lower_endpoint_algebra
    {S W LL logm c δ : ℝ}
    (hS : 1 ≤ S)
    (hLL : 0 < LL)
    (hδ : 0 < δ)
    (hδhalf : δ ≤ 1 / 2)
    (hWS : W ≤ S * ((1 / 2 - δ) * LL))
    (hlogm : LL / 2 - Real.log 2 ≤ logm)
    (hlog2 : Real.log 2 ≤ δ / 4 * LL)
    (hc : 2 * c ≤ δ / 8 * LL) :
    (S * (2 * S + 1)) * (2 * c + W - 2 * logm) <
      S ^ 2 * (2 * W - 2 * logm) := by
  have hsplit :=
    mertens_split_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := c) (δ := δ)
      hS hLL hδ hWS hlogm hlog2 hc
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have hlogm_pos : 0 < logm := by
    nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ 1 / 2 - δ) hLL.le]
  have hextra_pos :
      0 < (S * (2 * S + 1)) * logm - S ^ 2 * logm := by
    have hcoeff_gap : 0 < S ^ 2 + S := by nlinarith
    nlinarith [hcoeff_gap, hlogm_pos]
  nlinarith [hsplit, hextra_pos]

theorem mertens_threshold_algebra
    {S W LL G U logU c δ β : ℝ}
    (hS : 1 ≤ S)
    (hLL : 0 < LL)
    (hδ : 0 < δ)
    (hWS : W ≤ S * (β * LL))
    (hGgain : (β + 3 / 4 * δ) * LL ≤ G)
    (hc : 2 * c ≤ δ / 8 * LL)
    (hUterm : 2 * ((U + 1) * logU) ≤ δ / 8 * S ^ 2 * LL) :
    (S * (2 * S + 1)) * (2 * c + W) <
      (S ^ 2 * G - 2 * ((U + 1) * logU)) + S ^ 2 * (2 * W) := by
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS
  have h1 : S * W ≤ S * (S * (β * LL)) :=
    mul_le_mul_of_nonneg_left hWS hSpos.le
  have h2 : S ^ 2 * ((β + 3 / 4 * δ) * LL) ≤ S ^ 2 * G :=
    mul_le_mul_of_nonneg_left hGgain (sq_nonneg S)
  have h3 : (2 * S ^ 2 + S) * (2 * c) ≤ (2 * S ^ 2 + S) * (δ / 8 * LL) :=
    mul_le_mul_of_nonneg_left hc (by nlinarith [sq_nonneg S])
  nlinarith [h1, h2, h3, hUterm, mul_pos (mul_pos hδ hLL) hSpos,
    mul_nonneg (mul_nonneg hδ.le hLL.le) (mul_nonneg hSpos.le (sub_nonneg.mpr hS))]

theorem eventually_mertens_lower_endpoint_log_conditions
    (C q : ℝ) (hC : 0 < C) (hq : 2 < q) :
    ∀ᶠ R : ℝ in atTop,
      let s := mertensS q R
      let m₀ := sublogM₀ R
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)) ∧
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)) := by
  have hq0 : 0 < q := by linarith
  have hqne : q ≠ 0 := ne_of_gt hq0
  have h2q : (1 : ℝ) / q < 1 / 2 := by
    rw [div_lt_div_iff₀ hq0 (by norm_num : (0 : ℝ) < 2)]
    linarith
  have h1q_pos : (0 : ℝ) < 1 / q := by positivity
  set δ : ℝ := 1 / 2 - 1 / q with hδ_def
  have hδ : 0 < δ := by rw [hδ_def]; linarith
  have hδhalf : δ ≤ 1 / 2 := by rw [hδ_def]; linarith
  have hlogm_lower := eventually_log_sublogM₀_ge
  have hm₀_two := eventually_two_le_sublogM₀
  have hratio := eventually_le_log_div_loglog (1 / q)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hconst :
      ∀ᶠ R : ℝ in atTop, 2 * Real.log C ≤ δ / 8 * Real.log (Real.log R) :=
    eventually_const_le_mul_loglog (by linarith) _
  have hlog2_small :
      ∀ᶠ R : ℝ in atTop, Real.log 2 ≤ δ / 4 * Real.log (Real.log R) :=
    eventually_const_le_mul_loglog (by linarith) _
  filter_upwards [hlogm_lower, hm₀_two, hratio, hlog_gt_one, hR_gt_one, hconst,
    hlog2_small] with R hlogm_lower hm₀_two hratio hW_gt_one hR_gt_one hconst hlog2_small
  let s := mertensS q R
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
  have hLLne : LL ≠ 0 := ne_of_gt hLL_pos
  have hs_lower : q * (W / LL) ≤ S := by
    simpa [s, S, W, LL] using mertensS_lower q R
  have hratio' : 1 / q ≤ W / LL := by
    simpa [W, LL] using hratio
  have hS_one : 1 ≤ S := by
    have hmul : q * (1 / q) ≤ q * (W / LL) :=
      mul_le_mul_of_nonneg_left hratio' hq0.le
    have hcancel : q * (1 / q) = 1 := by field_simp
    linarith
  have hWS : W ≤ S * ((1 / 2 - δ) * LL) := by
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : q * (W / LL) * LL = q * W := by field_simp
    have hqW : q * W ≤ S * LL := by linarith
    have hhalf : (1 / 2 - δ) = 1 / q := by rw [hδ_def]; ring
    rw [hhalf, show S * (1 / q * LL) = (S * LL) / q by ring, le_div_iff₀ hq0]
    linarith
  have hlogm : LL / 2 - Real.log 2 ≤ logm := by
    have hraw : (1 / 2 : ℝ) * LL - Real.log 2 ≤ logm := by
      simpa [m₀, LL, W, logm] using hlogm_lower
    linarith
  have hlog2 : Real.log 2 ≤ δ / 4 * LL := by
    simpa [LL, W] using hlog2_small
  have hc : 2 * Real.log C ≤ δ / 8 * LL := by
    simpa [LL, W] using hconst
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
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C) (δ := δ)
      hS_one hLL_pos hδ hWS hlogm hlog2 hc
  have hinert :=
    mertens_inert_lower_endpoint_algebra
      (S := S) (W := W) (LL := LL) (logm := logm) (c := Real.log C) (δ := δ)
      hS_one hLL_pos hδ hδhalf hWS hlogm hlog2 hc
  constructor
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_split_scale] using hsplit
  · simpa [s, m₀, S, W, LL, logm, hcoeff, hs_sq, hlog_inert_scale] using hinert

set_option maxHeartbeats 400000 in
-- The chain of `nlinarith` calls and the final cast-normalizing `simpa` exceed the
-- default heartbeat budget, as in `eventually_crude_threshold_condition`.
theorem eventually_mertens_threshold_condition
    (C c q θ₁ : ℝ) (hC : 0 < C) (hq : 2 < q)
    (hgap : 1 / 2 + 1 / q < θ₁) (hθ1 : θ₁ < 1) :
    ∀ᶠ R : ℝ in atTop,
      let s := mertensS q R
      let m₀ := sublogM₀ R
      let m₁ := mertensM₁ θ₁ R
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * (Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c) -
            2 * (((m₁ + 1 : ℕ) : ℝ) * Real.log (m₁ : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R) := by
  have hq0 : 0 < q := by linarith
  have h1q_pos : (0 : ℝ) < 1 / q := by positivity
  have hθ0 : 0 < θ₁ := by linarith
  set δ : ℝ := θ₁ - 1 / 2 - 1 / q with hδ_def
  have hδ : 0 < δ := by rw [hδ_def]; linarith
  have hδne : δ ≠ 0 := ne_of_gt hδ
  have hδ_half : δ < 1 / 2 := by rw [hδ_def]; linarith
  have hfour := eventually_four_mul_mertensM₁_le_mertensS hq0 hθ1
  have hm₁_two := eventually_two_le_mertensM₁ hθ0
  have hlogm₁ := eventually_log_mertensM₁_ge hθ0
  have hlogm₀ := eventually_log_sublogM₀_le_half
  have hratio := eventually_le_log_div_loglog (4 / δ)
  have hlog_gt_one : ∀ᶠ R : ℝ in atTop, 1 < Real.log R :=
    Real.tendsto_log_atTop.eventually_gt_atTop 1
  have hR_gt_one : ∀ᶠ R : ℝ in atTop, 1 < R := eventually_gt_atTop 1
  have hconst :
      ∀ᶠ R : ℝ in atTop, 2 * Real.log C ≤ δ / 8 * Real.log (Real.log R) :=
    eventually_const_le_mul_loglog (by linarith) _
  have hmass :
      ∀ᶠ R : ℝ in atTop,
        Real.log 2 + c ≤ δ / 4 * Real.log (Real.log R) :=
    eventually_const_le_mul_loglog (by linarith) _
  filter_upwards [hfour, hm₁_two, hlogm₁, hlogm₀, hratio, hlog_gt_one,
    hR_gt_one, hconst, hmass] with
    R hfour hm₁_two hlogm₁ hlogm₀ hratio hW_gt_one hR_gt_one hconst hmass
  let s := mertensS q R
  let m₀ := sublogM₀ R
  let m₁ := mertensM₁ θ₁ R
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
  have hLLne : LL ≠ 0 := ne_of_gt hLL_pos
  have hs_lower : q * (W / LL) ≤ S := by
    simpa [s, S, W, LL] using mertensS_lower q R
  have hratio' : 4 / δ ≤ W / LL := by
    simpa [W, LL] using hratio
  have hS_eight : 8 / δ ≤ S := by
    have hmul : q * (4 / δ) ≤ q * (W / LL) :=
      mul_le_mul_of_nonneg_left hratio' hq0.le
    have h2mul : 2 * (4 / δ) ≤ q * (4 / δ) :=
      mul_le_mul_of_nonneg_right (by linarith)
        (le_of_lt (div_pos (by norm_num) hδ))
    have h8 : 8 / δ = 2 * (4 / δ) := by ring
    linarith
  have hS_one : 1 ≤ S := by
    have h8 : (1 : ℝ) ≤ 8 / δ := by
      rw [le_div_iff₀ hδ]
      linarith
    linarith
  have hSpos : 0 < S := lt_of_lt_of_le zero_lt_one hS_one
  have hWS : W ≤ S * (1 / q * LL) := by
    have hmul := mul_le_mul_of_nonneg_right hs_lower hLL_pos.le
    have hcancel : q * (W / LL) * LL = q * W := by
      field_simp
    have hqW : q * W ≤ S * LL := by linarith
    rw [show S * (1 / q * LL) = (S * LL) / q by ring, le_div_iff₀ hq0]
    linarith
  have hcLL : 2 * Real.log C ≤ δ / 8 * LL := by
    simpa [LL, W] using hconst
  have hmassLL : Real.log 2 + c ≤ δ / 4 * LL := by
    simpa [LL, W] using hmass
  have hθ_split : θ₁ = 1 / 2 + 1 / q + δ := by rw [hδ_def]; ring
  have hGgain :
      (1 / q + 3 / 4 * δ) * LL ≤ Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c := by
    have hm₁_log : θ₁ * LL - Real.log 2 ≤ Real.log (m₁ : ℝ) := by
      simpa [m₁, LL, W] using hlogm₁
    have hm₀_log : Real.log (m₀ : ℝ) ≤ (1 / 2 : ℝ) * LL := by
      simpa [m₀, LL, W] using hlogm₀
    rw [hθ_split] at hm₁_log
    nlinarith [hm₁_log, hm₀_log, hmassLL]
  have hU_two : (2 : ℝ) ≤ (m₁ : ℝ) := by
    have h : 2 ≤ m₁ := hm₁_two
    exact_mod_cast h
  have hU_pos : 0 < (m₁ : ℝ) := by linarith
  have hlogU_nonneg : 0 ≤ logU := Real.log_nonneg (by linarith)
  have hm₁_le_W : (m₁ : ℝ) ≤ W := mertensM₁_le_log hθ1.le hW_gt_one.le
  have hlogU_le_LL : logU ≤ LL := Real.log_le_log hU_pos hm₁_le_W
  have hfourU_le_S : (4 : ℝ) * (m₁ : ℝ) ≤ S := by
    have hcast : ((4 * m₁ : ℕ) : ℝ) ≤ (s : ℝ) := by
      exact_mod_cast (hfour : 4 * mertensM₁ θ₁ R ≤ mertensS q R)
    simpa [Nat.cast_mul, S] using hcast
  have hUplus_le : (m₁ : ℝ) + 1 ≤ S / 2 := by linarith
  have hUterm :
      2 * (((m₁ : ℝ) + 1) * logU) ≤ δ / 8 * S ^ 2 * LL := by
    have h1 : 2 * (((m₁ : ℝ) + 1) * logU) ≤ 2 * ((S / 2) * LL) := by
      have hhalf_nonneg : 0 ≤ S / 2 := by linarith
      have hmul := mul_le_mul hUplus_le hlogU_le_LL hlogU_nonneg hhalf_nonneg
      linarith
    have h2 : 2 * ((S / 2) * LL) = S * LL := by ring
    have h3 : S * LL ≤ δ / 8 * S ^ 2 * LL := by
      have hSS : S ≤ δ / 8 * S ^ 2 := by
        have hmul := mul_le_mul_of_nonneg_right hS_eight hSpos.le
        have hcancel : δ / 8 * (8 / δ * S) = S := by field_simp
        have hstep : δ / 8 * (8 / δ * S) ≤ δ / 8 * (S * S) :=
          mul_le_mul_of_nonneg_left hmul (by linarith)
        calc
          S = δ / 8 * (8 / δ * S) := hcancel.symm
          _ ≤ δ / 8 * (S * S) := hstep
          _ = δ / 8 * S ^ 2 := by ring
      exact mul_le_mul_of_nonneg_right hSS hLL_pos.le
    linarith
  have hthreshold :=
    mertens_threshold_algebra
      (S := S) (W := W) (LL := LL)
      (G := Real.log (m₁ : ℝ) - Real.log (m₀ : ℝ) - c)
      (U := (m₁ : ℝ)) (logU := logU) (c := Real.log C) (δ := δ) (β := 1 / q)
      hS_one hLL_pos hδ hWS hGgain hcLL hUterm
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
