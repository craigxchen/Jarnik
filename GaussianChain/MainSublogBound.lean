import GaussianChain.MainDichotomy

namespace GaussianChain
namespace MainSublogBound

open SubcriticalBound

/-- The scale used after the descent step: the original permitted length `C * sqrt R`,
reduced by the fixed lower prime endpoint `sqrt m₀`. -/
noncomputable def radiusScale (C R : ℝ) (m₀ : ℕ) : ℝ :=
  C * Real.sqrt R / Real.sqrt (m₀ : ℝ)

theorem radiusScale_pos {C R : ℝ} {m₀ : ℕ}
    (hC : 0 < C) (hR : 0 < R) (hm₀ : 0 < m₀) :
    0 < radiusScale C R m₀ := by
  have hsqrtR : 0 < Real.sqrt R := Real.sqrt_pos_of_pos hR
  have hm₀_real : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀
  have hsqrtm₀ : 0 < Real.sqrt (m₀ : ℝ) := Real.sqrt_pos_of_pos hm₀_real
  exact div_pos (mul_pos hC hsqrtR) hsqrtm₀

/-- The sharper scale for inert descent by a rational prime `p`: the descended length is
`C * sqrt R / p`. -/
noncomputable def inertRadiusScale (C R : ℝ) (p : ℕ) : ℝ :=
  C * Real.sqrt R / (p : ℝ)

theorem inertRadiusScale_pos {C R : ℝ} {p : ℕ}
    (hC : 0 < C) (hR : 0 < R) (hp : 0 < p) :
    0 < inertRadiusScale C R p := by
  unfold inertRadiusScale
  exact div_pos (mul_pos hC (Real.sqrt_pos_of_pos hR)) (by exact_mod_cast hp)

/-- The full original arc-length scale `C * sqrt R`. -/
noncomputable def fullRadiusScale (C R : ℝ) : ℝ :=
  C * Real.sqrt R

theorem fullRadiusScale_pos {C R : ℝ} (hC : 0 < C) (hR : 0 < R) :
    0 < fullRadiusScale C R := by
  unfold fullRadiusScale
  exact mul_pos hC (Real.sqrt_pos_of_pos hR)

theorem fullRadiusScale_sq {C R : ℝ} (hR : 0 ≤ R) :
    (fullRadiusScale C R) ^ 2 = C ^ 2 * R := by
  unfold fullRadiusScale
  rw [mul_pow, Real.sq_sqrt hR]

theorem fullRadiusScale_floor_pos_of_one_le
    {C R : ℝ}
    (hR : 0 ≤ R)
    (hone : 1 ≤ C ^ 2 * R) :
    0 < Nat.floor ((fullRadiusScale C R) ^ 2) := by
  apply Nat.floor_pos.mpr
  rwa [fullRadiusScale_sq hR]

theorem log_natFloor_le_log_of_le
    {x y : ℝ}
    (hx0 : 0 ≤ x)
    (hfloor : 0 < Nat.floor x)
    (hxy : x ≤ y) :
    Real.log (Nat.floor x : ℝ) ≤ Real.log y := by
  have hfloor_pos_real : 0 < (Nat.floor x : ℝ) := by exact_mod_cast hfloor
  have hfloor_le_x : (Nat.floor x : ℝ) ≤ x := Nat.floor_le hx0
  exact Real.log_le_log hfloor_pos_real (hfloor_le_x.trans hxy)

theorem fullRadiusScale_log_floor_le_of_sq_le
    {C R Y : ℝ}
    (hfloor : 0 < Nat.floor ((fullRadiusScale C R) ^ 2))
    (hY : (fullRadiusScale C R) ^ 2 ≤ Y) :
    Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤ Real.log Y :=
  log_natFloor_le_log_of_le (sq_nonneg _) hfloor hY

theorem fullRadiusScale_log_floor_le_log_of_sq_le
    {C R Y : ℝ}
    (hR : 0 ≤ R)
    (hfloor : 0 < Nat.floor ((fullRadiusScale C R) ^ 2))
    (hY : C ^ 2 * R ≤ Y) :
    Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤ Real.log Y := by
  exact fullRadiusScale_log_floor_le_of_sq_le hfloor
    ((fullRadiusScale_sq hR).trans_le hY)

theorem radiusScale_log_floor_le_of_sq_le
    {C R Y : ℝ} {p : ℕ}
    (hfloor : 0 < Nat.floor ((radiusScale C R p) ^ 2))
    (hY : (radiusScale C R p) ^ 2 ≤ Y) :
    Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) ≤ Real.log Y :=
  log_natFloor_le_log_of_le (sq_nonneg _) hfloor hY

theorem split_log_condition_mono_from_lower
    {s : ℕ} {X Y u v : ℝ}
    (huv : u ≤ v)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) * (X - u) <
        ((s * s : ℕ) : ℝ) * (Y - u)) :
    ((s * (2 * s + 1) : ℕ) : ℝ) * (X - v) <
      ((s * s : ℕ) : ℝ) * (Y - v) := by
  have hgap : 0 ≤ (((s * (2 * s + 1) : ℕ) : ℝ) - ((s * s : ℕ) : ℝ)) := by
    norm_num
    nlinarith [sq_nonneg (s : ℝ)]
  nlinarith

theorem inert_sq_log_condition_mono_from_lower
    {s : ℕ} {X Y u v : ℝ}
    (huv : u ≤ v)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) * (X - 2 * u) <
        ((s * s : ℕ) : ℝ) * (Y - 2 * u)) :
    ((s * (2 * s + 1) : ℕ) : ℝ) * (X - 2 * v) <
      ((s * s : ℕ) : ℝ) * (Y - 2 * v) := by
  have hgap :
      0 ≤ 2 * (((s * (2 * s + 1) : ℕ) : ℝ) - ((s * s : ℕ) : ℝ)) := by
    norm_num
    nlinarith [sq_nonneg (s : ℝ)]
  nlinarith

theorem split_real_scale_log_condition_of_lower_endpoint
    {s m₀ p : ℕ} {C R : ℝ}
    (hC : 0 < C) (hR : 0 < R)
    (hm₀ : 0 < m₀) (hp : Nat.Prime p) (hmp : m₀ < p)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ))) :
    ((s * (2 * s + 1) : ℕ) : ℝ) *
        Real.log (C ^ 2 * R / (p : ℝ)) <
      ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (p : ℝ)) := by
  have hCRpos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hm₀_real_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀
  have hp_real_pos : 0 < (p : ℝ) := by exact_mod_cast hp.pos
  have hlog_mono : Real.log (m₀ : ℝ) ≤ Real.log (p : ℝ) :=
    Real.log_le_log hm₀_real_pos (by exact_mod_cast hmp.le)
  have hlower' :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          (Real.log (C ^ 2 * R) - Real.log (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)) := by
    rwa [Real.log_div hCRpos.ne' hm₀_real_pos.ne'] at hlower
  have htarget :=
    split_log_condition_mono_from_lower
      (s := s) (X := Real.log (C ^ 2 * R)) (Y := 2 * Real.log R)
      (u := Real.log (m₀ : ℝ)) (v := Real.log (p : ℝ)) hlog_mono hlower'
  rwa [← Real.log_div hCRpos.ne' hp_real_pos.ne'] at htarget

theorem inert_sq_real_scale_log_condition_of_lower_endpoint
    {s m₀ p : ℕ} {C R : ℝ}
    (hC : 0 < C) (hR : 0 < R)
    (hm₀ : 0 < m₀) (hp : Nat.Prime p) (hmp : m₀ < p)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ))) :
    ((s * (2 * s + 1) : ℕ) : ℝ) *
        Real.log (C ^ 2 * R / (p : ℝ) ^ 2) <
      ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (p : ℝ)) := by
  have hCRpos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hm₀_real_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀
  have hp_real_pos : 0 < (p : ℝ) := by exact_mod_cast hp.pos
  have hm₀_sq_ne : (m₀ : ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 hm₀_real_pos.ne'
  have hp_sq_ne : (p : ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 hp_real_pos.ne'
  have hlog_mono : Real.log (m₀ : ℝ) ≤ Real.log (p : ℝ) :=
    Real.log_le_log hm₀_real_pos (by exact_mod_cast hmp.le)
  have hlog_m₀_sq : Real.log ((m₀ : ℝ) ^ 2) = 2 * Real.log (m₀ : ℝ) := by
    rw [pow_two, Real.log_mul hm₀_real_pos.ne' hm₀_real_pos.ne']
    ring
  have hlog_p_sq : Real.log ((p : ℝ) ^ 2) = 2 * Real.log (p : ℝ) := by
    rw [pow_two, Real.log_mul hp_real_pos.ne' hp_real_pos.ne']
    ring
  have hlower' :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          (Real.log (C ^ 2 * R) - 2 * Real.log (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)) := by
    rw [Real.log_div hCRpos.ne' hm₀_sq_ne, hlog_m₀_sq] at hlower
    exact hlower
  have htarget :=
    inert_sq_log_condition_mono_from_lower
      (s := s) (X := Real.log (C ^ 2 * R)) (Y := 2 * Real.log R)
      (u := Real.log (m₀ : ℝ)) (v := Real.log (p : ℝ)) hlog_mono hlower'
  rwa [Real.log_div hCRpos.ne' hp_sq_ne, hlog_p_sq]

theorem radiusScale_sq {C R : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (hp : 0 < p) :
    (radiusScale C R p) ^ 2 = C ^ 2 * R / (p : ℝ) := by
  unfold radiusScale
  rw [div_pow, mul_pow, Real.sq_sqrt hR, Real.sq_sqrt (by positivity : 0 ≤ (p : ℝ))]

theorem inertRadiusScale_sq {C R : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (_hp : 0 < p) :
    (inertRadiusScale C R p) ^ 2 = C ^ 2 * R / (p : ℝ) ^ 2 := by
  unfold inertRadiusScale
  rw [div_pow, mul_pow, Real.sq_sqrt hR]

theorem radiusScale_floor_pos_of_one_le
    {C R : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (hp : 0 < p)
    (hone : 1 ≤ C ^ 2 * R / (p : ℝ)) :
    0 < Nat.floor ((radiusScale C R p) ^ 2) := by
  apply Nat.floor_pos.mpr
  rwa [radiusScale_sq hR hp]

theorem radiusScale_log_floor_le_log_of_sq_le
    {C R Y : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (hp : 0 < p)
    (hfloor : 0 < Nat.floor ((radiusScale C R p) ^ 2))
    (hY : C ^ 2 * R / (p : ℝ) ≤ Y) :
    Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) ≤ Real.log Y := by
  exact radiusScale_log_floor_le_of_sq_le hfloor
    ((radiusScale_sq hR hp).trans_le hY)

theorem inertRadiusScale_floor_pos_of_one_le
    {C R : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (hp : 0 < p)
    (hone : 1 ≤ C ^ 2 * R / (p : ℝ) ^ 2) :
    0 < Nat.floor ((inertRadiusScale C R p) ^ 2) := by
  apply Nat.floor_pos.mpr
  rwa [inertRadiusScale_sq hR hp]

theorem inertRadiusScale_log_floor_le_log_of_sq_le
    {C R Y : ℝ} {p : ℕ}
    (hR : 0 ≤ R) (hp : 0 < p)
    (hfloor : 0 < Nat.floor ((inertRadiusScale C R p) ^ 2))
    (hY : C ^ 2 * R / (p : ℝ) ^ 2 ≤ Y) :
    Real.log (Nat.floor ((inertRadiusScale C R p) ^ 2) : ℝ) ≤ Real.log Y := by
  exact log_natFloor_le_log_of_le (sq_nonneg _) hfloor
    ((inertRadiusScale_sq hR hp).trans_le hY)

/-- If the real radius satisfies `R^2 = N`, then `log N = 2 log R`. -/
theorem log_nat_eq_two_log_radius_of_sq {N : ℕ} {R : ℝ}
    (hR : 0 < R) (hR2 : R ^ 2 = (N : ℝ)) :
    Real.log (N : ℝ) = 2 * Real.log R := by
  rw [← hR2, pow_two, Real.log_mul hR.ne' hR.ne']
  ring

end MainSublogBound
end GaussianChain
