import GaussianChain.IntervalStack
import Mathlib.Analysis.SpecificLimits.Basic

namespace GaussianChain
namespace MertensInterval

open Filter
open scoped Nat.Prime

/-- The weighted prime reciprocal sum over a real interval `lo ≤ p ≤ hi`. -/
noncomputable def weightedPrimeRealInterval (lo hi : ℝ) : ℝ :=
  ∑ p ∈ (Nat.primesLE (Nat.floor hi)).filter (fun p : ℕ => lo ≤ (p : ℝ)),
    Real.log (p : ℝ) / (p : ℝ)

/-- Lower natural endpoint used to approximate the interval `s ^ η ≤ p`. -/
noncomputable def lowerEndpoint (η : ℝ) (s : ℕ) : ℕ :=
  Nat.ceil ((s : ℝ) ^ η)

/-- Upper natural endpoint used to approximate the interval `p ≤ s / 4`. -/
noncomputable def upperEndpoint (s : ℕ) : ℕ :=
  Nat.floor ((s : ℝ) / 4)

/-- The coefficient in the logarithmic height of the geometric stack. -/
noncomputable def stackCoeff (η : ℝ) : ℝ :=
  (1 - η) / (8 * Real.log 8)

/-- A stack of this many eighth-scale intervals fits between `s ^ η` and `s / 4`
for all sufficiently large `s`, when `η < 1`. -/
noncomputable def stackHeight (η : ℝ) (s : ℕ) : ℕ :=
  Nat.floor (stackCoeff η * Real.log (s : ℝ))

theorem stackCoeff_pos {η : ℝ} (hη : η < 1) : 0 < stackCoeff η := by
  unfold stackCoeff
  exact div_pos (sub_pos.mpr hη)
    (mul_pos (by norm_num) (Real.log_pos (by norm_num : (1 : ℝ) < 8)))

/-- Ordinary weighted prime intervals are monotone in the upper endpoint. -/
theorem weightedPrimeInterval_mono_upper {m n₁ n₂ : ℕ} (hn : n₁ ≤ n₂) :
    MertensLower.weightedPrimeInterval m n₁ ≤
      MertensLower.weightedPrimeInterval m n₂ := by
  classical
  unfold MertensLower.weightedPrimeInterval
  refine Finset.sum_le_sum_of_subset_of_nonneg ?hsub ?hnonneg
  · intro p hp
    rw [Finset.mem_sdiff] at hp ⊢
    rcases hp with ⟨hpn₁, hpm⟩
    have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hpn₁
    have hp_le_n₂ : p ≤ n₂ := (Nat.le_of_mem_primesLE hpn₁).trans hn
    exact ⟨(Nat.mem_primesLE).mpr ⟨hp_le_n₂, hp_prime⟩, hpm⟩
  · intro p hp _hnot
    exact MertensLower.log_div_nonneg_of_mem_primesLE (Finset.mem_sdiff.mp hp).1

/-- The rounded natural interval is contained in the corresponding real interval. -/
theorem weightedPrimeInterval_le_weightedPrimeRealInterval
    {m n : ℕ} {lo hi : ℝ}
    (hlo : lo ≤ (m : ℝ)) (hhi : (n : ℝ) ≤ hi) :
    MertensLower.weightedPrimeInterval m n ≤ weightedPrimeRealInterval lo hi := by
  classical
  unfold MertensLower.weightedPrimeInterval weightedPrimeRealInterval
  refine Finset.sum_le_sum_of_subset_of_nonneg ?hsub ?hnonneg
  · intro p hp
    have hbounds := MertensLower.bounds_of_mem_primeIntervalFinset hp
    rw [Finset.mem_sdiff] at hp
    rcases hp with ⟨hpn, _hpm⟩
    have hp_prime : Nat.Prime p := Nat.prime_of_mem_primesLE hpn
    have hp_le_n_real : (p : ℝ) ≤ (n : ℝ) := by exact_mod_cast hbounds.2
    have hp_le_hi : (p : ℝ) ≤ hi := hp_le_n_real.trans hhi
    have hp_floor : p ≤ Nat.floor hi := Nat.le_floor hp_le_hi
    have hlo_p : lo ≤ (p : ℝ) :=
      hlo.trans (by exact_mod_cast hbounds.1.le)
    exact Finset.mem_filter.mpr
      ⟨(Nat.mem_primesLE).mpr ⟨hp_floor, hp_prime⟩, hlo_p⟩
  · intro p hp _hnot
    exact MertensLower.log_div_nonneg_of_mem_primesLE (Finset.mem_filter.mp hp).1

theorem tendsto_lowerEndpoint_atTop {η : ℝ} (hη : 0 < η) :
    Tendsto (lowerEndpoint η) atTop atTop := by
  unfold lowerEndpoint
  exact tendsto_nat_ceil_atTop.comp
    ((tendsto_rpow_atTop hη).comp tendsto_natCast_atTop_atTop)

theorem floor_half_le_of_two_le {x : ℝ} (hx : 2 ≤ x) :
    x / 2 ≤ (Nat.floor x : ℝ) := by
  have hfloor_gt : x - 1 < (Nat.floor x : ℝ) := by
    have hceil_like : x < (Nat.floor x : ℝ) + 1 := by
      simpa [add_comm] using Nat.lt_floor_add_one x
    linarith
  exact le_of_lt (lt_of_le_of_lt (by linarith) hfloor_gt)

theorem eventually_stackHeight_log_lower {η : ℝ} (hη : η < 1) :
    ∀ᶠ s : ℕ in atTop,
      (stackCoeff η * Real.log (s : ℝ)) / 2 ≤ (stackHeight η s : ℝ) := by
  have hA : 0 < stackCoeff η := stackCoeff_pos hη
  have hlog :
      Tendsto (fun s : ℕ => Real.log (s : ℝ)) atTop atTop :=
    Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop
  have hmul :
      Tendsto (fun s : ℕ => stackCoeff η * Real.log (s : ℝ)) atTop atTop :=
    hlog.const_mul_atTop hA
  filter_upwards [hmul.eventually_ge_atTop 2] with s hs
  exact floor_half_le_of_two_le hs

theorem eventually_geomLower_le_upperEndpoint {η : ℝ} (hη0 : 0 < η) (hη1 : η < 1) :
    ∀ᶠ s : ℕ in atTop,
      IntervalStack.geomLower (lowerEndpoint η s) (stackHeight η s) ≤ upperEndpoint s := by
  let stackExp : ℝ := (1 - η) / 8
  let alpha : ℝ := η + stackExp
  let gap : ℝ := 1 - alpha
  have hA : 0 < stackCoeff η := stackCoeff_pos hη1
  have hstackExp_nonneg : 0 ≤ stackExp := by
    dsimp [stackExp]
    exact div_nonneg (sub_nonneg.mpr hη1.le) (by norm_num)
  have hgap_pos : 0 < gap := by
    dsimp [gap, alpha, stackExp]
    linarith
  have hpowη :
      Tendsto (fun s : ℕ => (s : ℝ) ^ η) atTop atTop :=
    (tendsto_rpow_atTop hη0).comp tendsto_natCast_atTop_atTop
  have hpow_gap :
      Tendsto (fun s : ℕ => (s : ℝ) ^ gap) atTop atTop :=
    (tendsto_rpow_atTop hgap_pos).comp tendsto_natCast_atTop_atTop
  have hlog_nonneg : ∀ᶠ s : ℕ in atTop, 0 ≤ Real.log (s : ℝ) := by
    filter_upwards [eventually_ge_atTop (1 : ℕ)] with s hs
    exact Real.log_nonneg (by exact_mod_cast hs)
  filter_upwards
    [eventually_ge_atTop (1 : ℕ),
      hpowη.eventually_ge_atTop 1,
      hpow_gap.eventually_ge_atTop 8,
      hlog_nonneg] with s hs_pos hs_powη_one hs_gap hlog_s_nonneg
  have hs_pos_real : 0 < (s : ℝ) := by exact_mod_cast hs_pos
  have hs_nonneg_real : 0 ≤ (s : ℝ) := hs_pos_real.le
  have hs_one_real : 1 ≤ (s : ℝ) := by exact_mod_cast hs_pos
  have hceil_le :
      (lowerEndpoint η s : ℝ) ≤ 2 * (s : ℝ) ^ η := by
    unfold lowerEndpoint
    have hceil_lt :
        (Nat.ceil ((s : ℝ) ^ η) : ℝ) < (s : ℝ) ^ η + 1 :=
      Nat.ceil_lt_add_one (Real.rpow_nonneg hs_nonneg_real η)
    linarith
  have hheight_le :
      (stackHeight η s : ℝ) ≤ stackCoeff η * Real.log (s : ℝ) := by
    unfold stackHeight
    exact Nat.floor_le (mul_nonneg hA.le hlog_s_nonneg)
  have hpow_stack :
      (((8 ^ stackHeight η s : ℕ) : ℝ)) ≤ (s : ℝ) ^ stackExp := by
    calc
      (((8 ^ stackHeight η s : ℕ) : ℝ))
          = (8 : ℝ) ^ stackHeight η s := by norm_num
      _ = (8 : ℝ) ^ ((stackHeight η s : ℕ) : ℝ) := by
          exact (Real.rpow_natCast (8 : ℝ) (stackHeight η s)).symm
      _ ≤ (8 : ℝ) ^ (stackCoeff η * Real.log (s : ℝ)) :=
          Real.rpow_le_rpow_of_exponent_le (by norm_num) hheight_le
      _ = (s : ℝ) ^ stackExp := by
          rw [Real.rpow_def_of_pos (by norm_num : (0 : ℝ) < 8),
            Real.rpow_def_of_pos hs_pos_real]
          congr 1
          dsimp [stackCoeff, stackExp]
          have hlog8_ne : Real.log 8 ≠ 0 :=
            (Real.log_pos (by norm_num : (1 : ℝ) < 8)).ne'
          field_simp [hlog8_ne]
  have hgeom_real :
      ((IntervalStack.geomLower (lowerEndpoint η s) (stackHeight η s) : ℕ) : ℝ) ≤
        2 * (s : ℝ) ^ alpha := by
    unfold IntervalStack.geomLower
    rw [Nat.cast_mul]
    calc
      ((lowerEndpoint η s : ℝ) * (((8 ^ stackHeight η s : ℕ) : ℝ)))
          ≤ (2 * (s : ℝ) ^ η) * (s : ℝ) ^ stackExp := by
            exact mul_le_mul hceil_le hpow_stack
              (by exact_mod_cast (Nat.zero_le (8 ^ stackHeight η s)))
              (mul_nonneg (by norm_num) (Real.rpow_nonneg hs_nonneg_real η))
      _ = 2 * ((s : ℝ) ^ η * (s : ℝ) ^ stackExp) := by ring
      _ = 2 * (s : ℝ) ^ alpha := by
          rw [← Real.rpow_add hs_pos_real]
  have hquarter :
      2 * (s : ℝ) ^ alpha ≤ (s : ℝ) / 4 := by
    have halpha_nonneg : 0 ≤ alpha := by
      dsimp [alpha, stackExp]
      nlinarith [hη0.le, sub_nonneg.mpr hη1.le]
    have hpow_alpha_nonneg : 0 ≤ (s : ℝ) ^ alpha :=
      Real.rpow_nonneg hs_nonneg_real alpha
    have hmul :
        8 * (s : ℝ) ^ alpha ≤ (s : ℝ) ^ gap * (s : ℝ) ^ alpha :=
      mul_le_mul_of_nonneg_right hs_gap hpow_alpha_nonneg
    have hpow_mul :
        (s : ℝ) ^ gap * (s : ℝ) ^ alpha = (s : ℝ) := by
      rw [← Real.rpow_add hs_pos_real]
      have hgap_alpha : gap + alpha = 1 := by dsimp [gap]; ring
      rw [hgap_alpha, Real.rpow_one]
    have : 8 * (s : ℝ) ^ alpha ≤ (s : ℝ) := by
      simpa [hpow_mul] using hmul
    linarith
  have hgeom_quarter :
      ((IntervalStack.geomLower (lowerEndpoint η s) (stackHeight η s) : ℕ) : ℝ) ≤
        (s : ℝ) / 4 :=
    hgeom_real.trans hquarter
  unfold upperEndpoint
  exact Nat.le_floor hgeom_quarter

/-- A logarithmic Mertens lower bound on the rounded natural interval
`ceil(s ^ η) < p ≤ floor(s / 4)`. -/
theorem eventually_weightedPrimeInterval_lowerEndpoint_upperEndpoint
    {η : ℝ} (hη0 : 0 < η) (hη1 : η < 1) :
    ∃ c > 0, ∀ᶠ s : ℕ in atTop,
      c * Real.log (s : ℝ) ≤
        MertensLower.weightedPrimeInterval (lowerEndpoint η s) (upperEndpoint s) := by
  classical
  let c : ℝ := stackCoeff η * Real.log 2 / 4
  have hc_pos : 0 < c := by
    dsimp [c]
    exact div_pos (mul_pos (stackCoeff_pos hη1)
      (Real.log_pos (by norm_num : (1 : ℝ) < 2))) (by norm_num)
  refine ⟨c, hc_pos, ?_⟩
  obtain ⟨Ncheb, hcheb⟩ :=
    Filter.eventually_atTop.mp MertensLower.chebyshev_error_eventually_le_quarter_log_two
  have hm0_atTop := tendsto_lowerEndpoint_atTop hη0
  have hm0_ge_cheb : ∀ᶠ s : ℕ in atTop, Ncheb ≤ lowerEndpoint η s :=
    hm0_atTop.eventually_ge_atTop Ncheb
  have hpowη :
      Tendsto (fun s : ℕ => (s : ℝ) ^ η) atTop atTop :=
    (tendsto_rpow_atTop hη0).comp tendsto_natCast_atTop_atTop
  filter_upwards
    [eventually_stackHeight_log_lower hη1,
      eventually_geomLower_le_upperEndpoint hη0 hη1,
      hm0_ge_cheb,
      hpowη.eventually_ge_atTop 1] with s hheight hgeom_upper hm0_ge hpow_one
  let m₀ := lowerEndpoint η s
  let k := stackHeight η s
  have hm₀_pos : 0 < m₀ := by
    have hone : 1 ≤ m₀ := by
      dsimp [m₀, lowerEndpoint]
      have hone_real : (1 : ℝ) ≤ (Nat.ceil ((s : ℝ) ^ η) : ℝ) :=
        hpow_one.trans (Nat.le_ceil ((s : ℝ) ^ η))
      exact_mod_cast hone_real
    omega
  have herr_base : ∀ n, m₀ ≤ n →
      Real.log ((n : ℝ) + 1) + 2 * Real.sqrt (n : ℝ) * Real.log (n : ℝ) ≤
        (Real.log 2 / 4) * (n : ℝ) := by
    intro n hn
    exact hcheb n (hm0_ge.trans hn)
  have hstack :
      (k : ℝ) * (Real.log 2 / 2) ≤
        MertensLower.weightedPrimeInterval m₀ (IntervalStack.geomLower m₀ k) :=
    IntervalStack.weightedPrimeInterval_geom_stack_lower hm₀_pos
      (IntervalStack.chebyshev_error_on_geomUpper_of_base herr_base)
  have hmono :
      MertensLower.weightedPrimeInterval m₀ (IntervalStack.geomLower m₀ k) ≤
        MertensLower.weightedPrimeInterval m₀ (upperEndpoint s) :=
    weightedPrimeInterval_mono_upper hgeom_upper
  have hc_height :
      c * Real.log (s : ℝ) ≤ (k : ℝ) * (Real.log 2 / 2) := by
    have hlog2_half_nonneg : 0 ≤ Real.log 2 / 2 := by positivity
    calc
      c * Real.log (s : ℝ)
          = ((stackCoeff η * Real.log (s : ℝ)) / 2) * (Real.log 2 / 2) := by
              dsimp [c]
              ring
      _ ≤ (k : ℝ) * (Real.log 2 / 2) :=
          mul_le_mul_of_nonneg_right hheight hlog2_half_nonneg
  exact hc_height.trans (hstack.trans hmono)

/-- Ordinary-prime Mertens lower bound on the exact real interval
`s ^ η ≤ p ≤ s / 4`. This is the formal analytic dependency needed by the
logarithmic-width limiting argument. -/
theorem eventually_weightedPrimeRealInterval_rpow_quarter_lower
    {η : ℝ} (hη0 : 0 < η) (hη1 : η < 1) :
    ∃ c > 0, ∀ᶠ s : ℕ in atTop,
      c * Real.log (s : ℝ) ≤
        weightedPrimeRealInterval ((s : ℝ) ^ η) ((s : ℝ) / 4) := by
  obtain ⟨c, hc_pos, hc_eventual⟩ :=
    eventually_weightedPrimeInterval_lowerEndpoint_upperEndpoint hη0 hη1
  refine ⟨c, hc_pos, ?_⟩
  filter_upwards [hc_eventual, eventually_ge_atTop (1 : ℕ)] with s hs_lower hs_pos
  have hlo : (s : ℝ) ^ η ≤ (lowerEndpoint η s : ℝ) := by
    unfold lowerEndpoint
    exact Nat.le_ceil ((s : ℝ) ^ η)
  have hhi : (upperEndpoint s : ℝ) ≤ (s : ℝ) / 4 := by
    unfold upperEndpoint
    exact Nat.floor_le (by positivity)
  exact hs_lower.trans
    (weightedPrimeInterval_le_weightedPrimeRealInterval (m := lowerEndpoint η s)
      (n := upperEndpoint s) (lo := (s : ℝ) ^ η) (hi := (s : ℝ) / 4) hlo hhi)

end MertensInterval
end GaussianChain
