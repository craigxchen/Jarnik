import GaussianChain.MainDichotomy

namespace GaussianChain
namespace MainSublogBound

open SubcriticalBound

/-- The scale used after the descent step: the original permitted length `C * sqrt R`,
reduced by the fixed lower prime endpoint `sqrt m₀`. -/
noncomputable def radiusScale (C R : ℝ) (m₀ : ℕ) : ℝ :=
  C * Real.sqrt R / Real.sqrt (m₀ : ℝ)

theorem sqrt_nat_le_self_of_one_le {m : ℕ} (hm : 1 ≤ m) :
    Real.sqrt (m : ℝ) ≤ (m : ℝ) := by
  have hm_nonneg : 0 ≤ (m : ℝ) := by positivity
  rw [Real.sqrt_le_left hm_nonneg]
  have hm_real : 1 ≤ (m : ℝ) := by exact_mod_cast hm
  nlinarith [sq_nonneg ((m : ℝ) - 1)]

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

theorem radiusScale_many_term_le
    {s m₀ : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hm₀ : 1 ≤ m₀)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * L / radiusScale C R m₀ ≤
      (m₀ : ℝ) * ((2 * s : ℕ) : ℝ) := by
  have hm₀_pos : 0 < m₀ := Nat.pos_of_ne_zero (by omega)
  have hsqrtR : 0 < Real.sqrt R := Real.sqrt_pos_of_pos hR
  have hm₀_real_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀_pos
  have hsqrtm₀_pos : 0 < Real.sqrt (m₀ : ℝ) := Real.sqrt_pos_of_pos hm₀_real_pos
  have hCR_pos : 0 < C * Real.sqrt R := mul_pos hC hsqrtR
  have hdiv_le_sqrt :
      L / radiusScale C R m₀ ≤ Real.sqrt (m₀ : ℝ) := by
    calc
      L / radiusScale C R m₀ =
          (L / (C * Real.sqrt R)) * Real.sqrt (m₀ : ℝ) := by
            unfold radiusScale
            field_simp [hCR_pos.ne', hsqrtm₀_pos.ne']
      _ ≤ 1 * Real.sqrt (m₀ : ℝ) := by
            exact mul_le_mul_of_nonneg_right
              (div_le_one_of_le₀ hLbound hCR_pos.le) (Real.sqrt_nonneg _)
      _ = Real.sqrt (m₀ : ℝ) := by ring
  have hdiv_le_m₀ : L / radiusScale C R m₀ ≤ (m₀ : ℝ) :=
    hdiv_le_sqrt.trans (sqrt_nat_le_self_of_one_le hm₀)
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  calc
    ((2 * s : ℕ) : ℝ) * L / radiusScale C R m₀ =
        ((2 * s : ℕ) : ℝ) * (L / radiusScale C R m₀) := by ring
    _ ≤ ((2 * s : ℕ) : ℝ) * (m₀ : ℝ) :=
        mul_le_mul_of_nonneg_left hdiv_le_m₀ hS_nonneg
    _ = (m₀ : ℝ) * ((2 * s : ℕ) : ℝ) := by ring

theorem radiusScale_descent_term_le
    {s m₀ : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hm₀ : 0 < m₀)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) /
        radiusScale C R m₀ ≤
      ((2 * s : ℕ) : ℝ) := by
  have hsqrtR : 0 < Real.sqrt R := Real.sqrt_pos_of_pos hR
  have hm₀_real_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀
  have hsqrtm₀_pos : 0 < Real.sqrt (m₀ : ℝ) := Real.sqrt_pos_of_pos hm₀_real_pos
  have hCR_pos : 0 < C * Real.sqrt R := mul_pos hC hsqrtR
  have hquot_le_one :
      (L / Real.sqrt (m₀ : ℝ)) / radiusScale C R m₀ ≤ 1 := by
    calc
      (L / Real.sqrt (m₀ : ℝ)) / radiusScale C R m₀ =
          L / (C * Real.sqrt R) := by
            unfold radiusScale
            field_simp [hCR_pos.ne', hsqrtm₀_pos.ne']
      _ ≤ 1 := div_le_one_of_le₀ hLbound hCR_pos.le
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  calc
    ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) /
        radiusScale C R m₀ =
        ((2 * s : ℕ) : ℝ) *
          ((L / Real.sqrt (m₀ : ℝ)) / radiusScale C R m₀) := by ring
    _ ≤ ((2 * s : ℕ) : ℝ) * 1 :=
        mul_le_mul_of_nonneg_left hquot_le_one hS_nonneg
    _ = ((2 * s : ℕ) : ℝ) := by ring

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

theorem fullRadiusScale_term_le
    {s : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * L / fullRadiusScale C R ≤
      ((2 * s : ℕ) : ℝ) := by
  have hscale_pos : 0 < fullRadiusScale C R := fullRadiusScale_pos hC hR
  have hquot : L / fullRadiusScale C R ≤ 1 := by
    unfold fullRadiusScale
    exact div_le_one_of_le₀ hLbound hscale_pos.le
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  calc
    ((2 * s : ℕ) : ℝ) * L / fullRadiusScale C R =
        ((2 * s : ℕ) : ℝ) * (L / fullRadiusScale C R) := by ring
    _ ≤ ((2 * s : ℕ) : ℝ) * 1 :=
        mul_le_mul_of_nonneg_left hquot hS_nonneg
    _ = ((2 * s : ℕ) : ℝ) := by ring

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

theorem logSplit_of_log_floor_le
    {s N p : ℕ} {A X : ℝ}
    (hlogA : Real.log (Nat.floor (A ^ 2) : ℝ) ≤ X)
    (hineq :
      ((s * (2 * s + 1) : ℕ) : ℝ) * X <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ))) :
    ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
      ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)) := by
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  exact (mul_le_mul_of_nonneg_left hlogA hcoeff_nonneg).trans_lt hineq

theorem logInert_of_log_floor_le
    {s N p : ℕ} {A X : ℝ}
    (hlogA : Real.log (Nat.floor (A ^ 2) : ℝ) ≤ X)
    (hineq :
      ((s * (2 * s + 1) : ℕ) : ℝ) * X <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ))) :
    ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (A ^ 2) : ℝ) <
      ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)) := by
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  exact (mul_le_mul_of_nonneg_left hlogA hcoeff_nonneg).trans_lt hineq

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

theorem inert_log_condition_mono_from_lower
    {s : ℕ} {X Y u v : ℝ}
    (huv : u ≤ v)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) * (X - u) <
        ((s * s : ℕ) : ℝ) * (Y - 2 * u)) :
    ((s * (2 * s + 1) : ℕ) : ℝ) * (X - v) <
      ((s * s : ℕ) : ℝ) * (Y - 2 * v) := by
  have hgap :
      0 ≤ (((s * (2 * s + 1) : ℕ) : ℝ) -
        2 * ((s * s : ℕ) : ℝ)) := by
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

theorem inert_real_scale_log_condition_of_lower_endpoint
    {s m₀ p : ℕ} {C R : ℝ}
    (hC : 0 < C) (hR : 0 < R)
    (hm₀ : 0 < m₀) (hp : Nat.Prime p) (hmp : m₀ < p)
    (hlower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ))) :
    ((s * (2 * s + 1) : ℕ) : ℝ) *
        Real.log (C ^ 2 * R / (p : ℝ)) <
      ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (p : ℝ)) := by
  have hCRpos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hm₀_real_pos : 0 < (m₀ : ℝ) := by exact_mod_cast hm₀
  have hp_real_pos : 0 < (p : ℝ) := by exact_mod_cast hp.pos
  have hlog_mono : Real.log (m₀ : ℝ) ≤ Real.log (p : ℝ) :=
    Real.log_le_log hm₀_real_pos (by exact_mod_cast hmp.le)
  have hlower' :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          (Real.log (C ^ 2 * R) - Real.log (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)) := by
    rwa [Real.log_div hCRpos.ne' hm₀_real_pos.ne'] at hlower
  have htarget :=
    inert_log_condition_mono_from_lower
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

theorem primeRadiusScale_split_term_le
    {s p : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hp : 0 < p)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) /
        radiusScale C R p ≤
    ((2 * s : ℕ) : ℝ) :=
  radiusScale_descent_term_le (s := s) (m₀ := p) hC hR hp hLbound

theorem primeRadiusScale_split_term_le_div_const
    {s p : ℕ} {C R L K : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hp : 0 < p)
    (hK : 0 < K)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) /
        (K * radiusScale C R p) ≤
      (1 / K) * ((2 * s : ℕ) : ℝ) := by
  have hbase :=
    primeRadiusScale_split_term_le (s := s) (p := p) hC hR hp hLbound
  have hscale_pos : 0 < radiusScale C R p := radiusScale_pos hC hR hp
  calc
    ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) /
        (K * radiusScale C R p) =
        (1 / K) *
          (((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) /
            radiusScale C R p) := by
          field_simp [hK.ne', hscale_pos.ne']
    _ ≤ (1 / K) * ((2 * s : ℕ) : ℝ) :=
      mul_le_mul_of_nonneg_left hbase (by positivity)

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

theorem primeRadiusScale_inert_term_le
    {s p : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hp : 1 ≤ p)
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) /
        radiusScale C R p ≤
      ((2 * s : ℕ) : ℝ) := by
  have hp_pos : 0 < p := by omega
  have hp_real_pos : 0 < (p : ℝ) := by exact_mod_cast hp_pos
  have hsqrtp_pos : 0 < Real.sqrt (p : ℝ) := Real.sqrt_pos_of_pos hp_real_pos
  have hscale_pos : 0 < radiusScale C R p := radiusScale_pos hC hR hp_pos
  have hsqrtp_le_p : Real.sqrt (p : ℝ) ≤ (p : ℝ) := sqrt_nat_le_self_of_one_le hp
  have hdiv :
      L / (p : ℝ) ≤ L / Real.sqrt (p : ℝ) :=
    div_le_div_of_nonneg_left hL hsqrtp_pos hsqrtp_le_p
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  have hterm :
      ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / radiusScale C R p ≤
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / radiusScale C R p := by
    exact div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_left hdiv hS_nonneg) hscale_pos.le
  exact hterm.trans
    (primeRadiusScale_split_term_le (s := s) (p := p) hC hR hp_pos hLbound)

theorem inertRadiusScale_inert_term_le
    {s p : ℕ} {C R L : ℝ}
    (hC : 0 < C) (hR : 0 < R) (hp : 0 < p)
    (hLbound : L ≤ C * Real.sqrt R) :
    ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) /
        inertRadiusScale C R p ≤
      ((2 * s : ℕ) : ℝ) := by
  have hp_real_pos : 0 < (p : ℝ) := by exact_mod_cast hp
  have hscale_pos : 0 < C * Real.sqrt R := mul_pos hC (Real.sqrt_pos_of_pos hR)
  have hquot : (L / (p : ℝ)) / inertRadiusScale C R p ≤ 1 := by
    calc
      (L / (p : ℝ)) / inertRadiusScale C R p = L / (C * Real.sqrt R) := by
        unfold inertRadiusScale
        field_simp [hp_real_pos.ne', hscale_pos.ne']
      _ ≤ 1 := div_le_one_of_le₀ hLbound hscale_pos.le
  have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
  calc
    ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / inertRadiusScale C R p =
        ((2 * s : ℕ) : ℝ) * ((L / (p : ℝ)) / inertRadiusScale C R p) := by ring
    _ ≤ ((2 * s : ℕ) : ℝ) * 1 :=
        mul_le_mul_of_nonneg_left hquot hS_nonneg
    _ = ((2 * s : ℕ) : ℝ) := by ring

/-- If the real radius satisfies `R^2 = N`, then `log N = 2 log R`. -/
theorem log_nat_eq_two_log_radius_of_sq {N : ℕ} {R : ℝ}
    (hR : 0 < R) (hR2 : R ^ 2 = (N : ℝ)) :
    Real.log (N : ℝ) = 2 * Real.log R := by
  rw [← hR2, pow_two, Real.log_mul hR.ne' hR.ne']
  ring

/-- Radius-scale finite bound with an explicit factor depending on the fixed lower endpoint
`m₀`. This is the finite form closest to the final `O_C(log R / log log R)` theorem: once the
analytic choices make `s` comparable to `log R / log log R`, this theorem supplies the
cardinality comparison. -/
theorem card_le_radiusScale_const_mul_s_of_interval_stack_log
    {M s N m₀ k : ℕ}
    {C R a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      MainDichotomy.manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k)
          (radiusScale C R m₀) ∨
        MainDichotomy.fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  have hm₀ : 0 < m₀ := by omega
  have hA : 0 < radiusScale C R m₀ := radiusScale_pos hC hR hm₀
  have htermMany :
      ((2 * s : ℕ) : ℝ) * L / radiusScale C R m₀ ≤
        (m₀ : ℝ) * ((2 * s : ℕ) : ℝ) :=
    radiusScale_many_term_le hC hR (by omega) hLbound
  have htermDesc :
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) /
          radiusScale C R m₀ ≤
        (m₀ : ℝ) * ((2 * s : ℕ) : ℝ) := by
    have hbase :=
      radiusScale_descent_term_le (s := s) (m₀ := m₀) hC hR hm₀ hLbound
    have hS_nonneg : 0 ≤ ((2 * s : ℕ) : ℝ) := by positivity
    have hm₀_one : 1 ≤ (m₀ : ℝ) := by exact_mod_cast (by omega : 1 ≤ m₀)
    calc
      ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (m₀ : ℝ)) /
          radiusScale C R m₀ ≤ ((2 * s : ℕ) : ℝ) := hbase
      _ ≤ (m₀ : ℝ) * ((2 * s : ℕ) : ℝ) := by
        nlinarith [hS_nonneg, hm₀_one]
  exact MainDichotomy.card_le_two_mul_succ_q_mul_two_s_of_interval_stack_log
    (M := M) (s := s) (q := m₀) (N := N) (m₀ := m₀) (k := k)
    (a := a) (L := L) (A := radiusScale C R m₀) (z := z) (t := t)
    hdichotomy hk hm₀ h2m₀ hM hA hfloor hL hsmallPrime hN
    htermMany htermDesc herr hlogSplit hlogInert hcircle hz hmono hparam hmem

/-- Big-O-shaped corollary of the radius-scale finite bound. The remaining analytic work in
the full Jarnik theorem is precisely to instantiate the hypotheses so that `s` has this
`log R / log log R` upper bound. -/
theorem card_le_radiusScale_const_mul_log_div_loglog_of_interval_stack_log
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hdichotomy :
      MainDichotomy.manyMissingWeightedLogCondition s N m₀ (IntervalStack.geomLower m₀ k)
          (radiusScale C R m₀) ∨
        MainDichotomy.fewMissingStackCondition N m₀ k)
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
        (Real.log R / Real.log (Real.log R)) := by
  have hfinite :=
    card_le_radiusScale_const_mul_s_of_interval_stack_log
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
      hdichotomy hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN
      herr hlogSplit hlogInert hcircle hz hmono hparam hmem
  have hcoeff_nonneg : 0 ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) := by positivity
  calc
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := hfinite
    _ = 4 * ((m₀ + 1 : ℕ) : ℝ) * (s : ℝ) := by norm_num; ring
    _ ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) *
          (B * (Real.log R / Real.log (Real.log R))) := by
        exact mul_le_mul_of_nonneg_left hs hcoeff_nonneg
    _ = (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
          (Real.log R / Real.log (Real.log R)) := by ring

/-- Radius-scale finite bound from the explicit stack threshold inequality, rather than an
externally supplied disjunction. -/
theorem card_le_radiusScale_const_mul_s_of_stack_threshold
    {M s N m₀ k : ℕ}
    {C R a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀
              (IntervalStack.geomLower m₀ k)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  exact card_le_radiusScale_const_mul_s_of_interval_stack_log
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
    (MainDichotomy.missing_or_few_stack_dichotomy_of_threshold
      (s := s) (N := N) (m₀ := m₀) (k := k) (A := radiusScale C R m₀)
      hthreshold)
    hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN
    herr hlogSplit hlogInert hcircle hz hmono hparam hmem

/-- Big-O-shaped radius-scale corollary from the explicit stack threshold inequality. -/
theorem card_le_radiusScale_const_mul_log_div_loglog_of_stack_threshold
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * MissingPrimeIntervalBranch.missingPrimeLogSum N m₀
              (IntervalStack.geomLower m₀ k)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
        (Real.log R / Real.log (Real.log R)) := by
  have hfinite :=
    card_le_radiusScale_const_mul_s_of_stack_threshold
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
      hthreshold hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN
      herr hlogSplit hlogInert hcircle hz hmono hparam hmem
  have hcoeff_nonneg : 0 ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) := by positivity
  calc
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := hfinite
    _ = 4 * ((m₀ + 1 : ℕ) : ℝ) * (s : ℝ) := by norm_num; ring
    _ ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) *
          (B * (Real.log R / Real.log (Real.log R))) := by
        exact mul_le_mul_of_nonneg_left hs hcoeff_nonneg
    _ = (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
          (Real.log R / Real.log (Real.log R)) := by ring

/-- Radius-scale finite bound from a crude stack threshold that contains no unweighted
missing-prime log sum. -/
theorem card_le_radiusScale_const_mul_s_of_crude_stack_threshold
    {M s N m₀ k : ℕ}
    {C R a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := by
  have hU : 1 ≤ IntervalStack.geomLower m₀ k := by
    unfold IntervalStack.geomLower
    have hm₀ : 1 ≤ m₀ := by omega
    have hpow : 1 ≤ 8 ^ k := Nat.one_le_pow k 8 (by norm_num)
    simpa using Nat.mul_le_mul hm₀ hpow
  exact card_le_radiusScale_const_mul_s_of_interval_stack_log
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
    (MainDichotomy.missing_or_few_stack_dichotomy_of_crude_threshold
      (s := s) (N := N) (m₀ := m₀) (k := k) (A := radiusScale C R m₀)
      hU hthreshold)
    hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN
    herr hlogSplit hlogInert hcircle hz hmono hparam hmem

/-- Big-O-shaped radius-scale corollary from the crude stack threshold. -/
theorem card_le_radiusScale_const_mul_log_div_loglog_of_crude_stack_threshold
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
        (Real.log R / Real.log (Real.log R)) := by
  have hfinite :=
    card_le_radiusScale_const_mul_s_of_crude_stack_threshold
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
      hthreshold hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN
      herr hlogSplit hlogInert hcircle hz hmono hparam hmem
  have hcoeff_nonneg : 0 ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) := by positivity
  calc
    (M : ℝ) ≤ 2 * (((m₀ + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := hfinite
    _ = 4 * ((m₀ + 1 : ℕ) : ℝ) * (s : ℝ) := by norm_num; ring
    _ ≤ 4 * ((m₀ + 1 : ℕ) : ℝ) *
          (B * (Real.log R / Real.log (Real.log R))) := by
        exact mul_le_mul_of_nonneg_left hs hcoeff_nonneg
    _ = (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
          (Real.log R / Real.log (Real.log R)) := by ring

/-- Big-O-shaped radius-scale corollary with the split/inert descent checks reduced to the
single stack upper endpoint. -/
theorem card_le_radiusScale_const_mul_log_div_loglog_of_crude_threshold_and_endpoint_checks
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hsplitUpper :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) *
          (Real.log (N : ℝ) - Real.log (IntervalStack.geomLower m₀ k : ℝ)))
    (hinertUpper :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) *
          (Real.log (N : ℝ) - 2 * Real.log (IntervalStack.geomLower m₀ k : ℝ)))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
        (Real.log R / Real.log (Real.log R)) := by
  exact card_le_radiusScale_const_mul_log_div_loglog_of_crude_stack_threshold
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (B := B) (a := a) (L := L) (z := z) (t := t)
    hthreshold hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN hs herr
    (MainDichotomy.logSplit_condition_of_upper_endpoint
      (s := s) (N := N) (m₀ := m₀) (U := IntervalStack.geomLower m₀ k)
      (A := radiusScale C R m₀) hsplitUpper)
    (MainDichotomy.logInert_condition_of_upper_endpoint
      (s := s) (N := N) (m₀ := m₀) (U := IntervalStack.geomLower m₀ k)
      (A := radiusScale C R m₀) hinertUpper)
    hcircle hz hmono hparam hmem

/-- Radius-squared version of
`card_le_radiusScale_const_mul_log_div_loglog_of_crude_threshold_and_endpoint_checks`.
The remaining inequalities are expressed using `2 * log R`, matching the geometric hypothesis
`R^2 = n`. -/
theorem card_le_radiusScale_const_mul_log_div_loglog_of_radius_sq_checks
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hsplitUpper :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) *
          (2 * Real.log R - Real.log (IntervalStack.geomLower m₀ k : ℝ)))
    (hinertUpper :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((radiusScale C R m₀) ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) *
          (2 * Real.log R - 2 * Real.log (IntervalStack.geomLower m₀ k : ℝ)))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloor : 0 < Nat.floor ((radiusScale C R m₀) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      (4 * ((m₀ + 1 : ℕ) : ℝ) * B) *
        (Real.log R / Real.log (Real.log R)) := by
  have hlogN := log_nat_eq_two_log_radius_of_sq (N := N) hR hR2
  exact card_le_radiusScale_const_mul_log_div_loglog_of_crude_threshold_and_endpoint_checks
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (B := B) (a := a) (L := L) (z := z) (t := t)
    (by simpa [hlogN] using hthreshold)
    (by simpa [hlogN] using hsplitUpper)
    (by simpa [hlogN] using hinertUpper)
    hk h2m₀ hM hC hR hfloor hL hLbound hsmallPrime hN hs herr
    hcircle hz hmono hparam hmem

/-- Sharper finite radius-scale bound with a prime-dependent descent scale.

The many-missing branch uses the original scale `C * sqrt R`; after a prime divisor `p` is found,
the descent branch uses the actual scale `C * sqrt R / sqrt p`. This removes the extra factor
depending on the stack lower endpoint. -/
theorem card_le_eight_mul_s_of_full_and_prime_radius_scales
    {M s N m₀ k : ℕ}
    {C R a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2))
    (hfloorPrime : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        0 < Nat.floor ((radiusScale C R p) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ (8 * s : ℕ) := by
  have hm₀ : 0 < m₀ := by omega
  have hU : 1 ≤ IntervalStack.geomLower m₀ k := by
    unfold IntervalStack.geomLower
    have hm₀_one : 1 ≤ m₀ := by omega
    have hpow : 1 ≤ 8 ^ k := Nat.one_le_pow k 8 (by norm_num)
    simpa using Nat.mul_le_mul hm₀_one hpow
  have hdichotomy :=
    MainDichotomy.missing_or_few_stack_dichotomy_of_crude_threshold
      (s := s) (N := N) (m₀ := m₀) (k := k) (A := fullRadiusScale C R)
      hU hthreshold
  have hfinite :=
    MainDichotomy.card_le_two_mul_succ_q_mul_two_s_of_interval_stack_primeScale_log
      (M := M) (s := s) (q := 1) (N := N) (m₀ := m₀) (k := k)
      (a := a) (L := L) (A₀ := fullRadiusScale C R)
      (A := fun p => radiusScale C R p) (z := z) (t := t)
      hdichotomy hk hm₀ h2m₀ hM (fullRadiusScale_pos hC hR) hfloorFull
      hsmallPrime hN
      (by simpa using fullRadiusScale_term_le (s := s) hC hR hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using
          primeRadiusScale_inert_term_le (s := s) (p := p) hC hR hp.one_le hL hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using
          primeRadiusScale_split_term_le (s := s) (p := p) hC hR hp.pos hLbound)
      herr
      (fun p hp _hmp _hpU _hpN => radiusScale_pos hC hR hp.pos)
      hfloorPrime hlogSplit hlogInert hcircle hz hmono hparam hmem
  calc
    (M : ℝ) ≤ 2 * (((1 + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := hfinite
    _ ≤ (8 * s : ℕ) := by
      norm_num
      ring_nf
      exact le_rfl

/-- Big-O-shaped corollary of `card_le_eight_mul_s_of_full_and_prime_radius_scales`. -/
theorem card_le_eight_mul_log_div_loglog_of_full_and_prime_radius_scales
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2))
    (hfloorPrime : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        0 < Nat.floor ((radiusScale C R p) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 8 * B * (Real.log R / Real.log (Real.log R)) := by
  have hfinite :=
    card_le_eight_mul_s_of_full_and_prime_radius_scales
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (C := C) (R := R) (a := a) (L := L) (z := z) (t := t)
      hthreshold hk h2m₀ hM hC hR hfloorFull hfloorPrime hL hLbound
      hsmallPrime hN herr hlogSplit hlogInert hcircle hz hmono hparam hmem
  calc
    (M : ℝ) ≤ (8 * s : ℕ) := hfinite
    _ = 8 * (s : ℝ) := by norm_num
    _ ≤ 8 * (B * (Real.log R / Real.log (Real.log R))) := by
      nlinarith [hs]
    _ = 8 * B * (Real.log R / Real.log (Real.log R)) := by ring

/-- Radius-squared version of the sharper prime-scaled bound. -/
theorem card_le_eight_mul_log_div_loglog_of_radius_sq_prime_scales
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2))
    (hfloorPrime : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        0 < Nat.floor ((radiusScale C R p) ^ 2))
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) <
            ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 8 * B * (Real.log R / Real.log (Real.log R)) := by
  have hlogN := log_nat_eq_two_log_radius_of_sq (N := N) hR hR2
  exact card_le_eight_mul_log_div_loglog_of_full_and_prime_radius_scales
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (B := B) (a := a) (L := L) (z := z) (t := t)
    (by simpa [hlogN] using hthreshold)
    hk h2m₀ hM hC hR hfloorFull hfloorPrime hL hLbound
    hsmallPrime hN hs herr
    (fun p N' hp hmp hpU hpN hN' hfactor => by
      simpa [hlogN] using hlogSplit p N' hp hmp hpU hpN hN' hfactor)
    (fun p N' hp hmp hpU hpN hN' hfactor => by
      simpa [hlogN] using hlogInert p N' hp hmp hpU hpN hN' hfactor)
    hcircle hz hmono hparam hmem

/-- Final-facing prime-scaled bound with `Nat.floor` hidden behind elementary real estimates.

The caller proves positivity of the full and prime scales through `1 ≤ C^2 R` and
`1 ≤ C^2 R / p`, and proves the logarithmic inequalities using the simpler real logarithms
`log (C^2 R)` and `log (C^2 R / p)`. -/
theorem card_le_eight_mul_log_div_loglog_of_radius_sq_real_scale_bounds
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (honeFull : 1 ≤ C ^ 2 * R)
    (honePrime : ∀ p, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N →
        1 ≤ C ^ 2 * R / (p : ℝ))
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hlogSplit : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (C ^ 2 * R / (p : ℝ)) <
            ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (p : ℝ)))
    (hlogInert : ∀ p N' : ℕ, Nat.Prime p → m₀ < p →
      p ≤ IntervalStack.geomLower m₀ k → p ∣ N → 0 < N' →
        (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ) →
          ((s * (2 * s + 1) : ℕ) : ℝ) *
              Real.log (C ^ 2 * R / (p : ℝ)) <
            ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (p : ℝ)))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 8 * B * (Real.log R / Real.log (Real.log R)) := by
  have hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2) :=
    fullRadiusScale_floor_pos_of_one_le hR.le honeFull
  have hlogFull :
      Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤
        Real.log (C ^ 2 * R) :=
    fullRadiusScale_log_floor_le_log_of_sq_le hR.le hfloorFull le_rfl
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  have hthresholdFloor :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R) :=
    (mul_le_mul_of_nonneg_left hlogFull hcoeff_nonneg).trans_lt hthreshold
  exact card_le_eight_mul_log_div_loglog_of_radius_sq_prime_scales
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (B := B) (a := a) (L := L) (z := z) (t := t)
    hR2 hthresholdFloor hk h2m₀ hM hC hR hfloorFull
    (fun p hp hmp hpU hpN =>
      radiusScale_floor_pos_of_one_le hR.le hp.pos (honePrime p hp hmp hpU hpN))
    hL hLbound hsmallPrime hN hs herr
    (fun p N' hp hmp hpU hpN hN' hfactor => by
      have hfloorP :=
        radiusScale_floor_pos_of_one_le hR.le hp.pos (honePrime p hp hmp hpU hpN)
      have hlogP :
          Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) ≤
            Real.log (C ^ 2 * R / (p : ℝ)) :=
        radiusScale_log_floor_le_log_of_sq_le hR.le hp.pos hfloorP le_rfl
      exact (mul_le_mul_of_nonneg_left hlogP hcoeff_nonneg).trans_lt
        (hlogSplit p N' hp hmp hpU hpN hN' hfactor))
    (fun p N' hp hmp hpU hpN hN' hfactor => by
      have hfloorP :=
        radiusScale_floor_pos_of_one_le hR.le hp.pos (honePrime p hp hmp hpU hpN)
      have hlogP :
          Real.log (Nat.floor ((radiusScale C R p) ^ 2) : ℝ) ≤
            Real.log (C ^ 2 * R / (p : ℝ)) :=
        radiusScale_log_floor_le_log_of_sq_le hR.le hp.pos hfloorP le_rfl
      exact (mul_le_mul_of_nonneg_left hlogP hcoeff_nonneg).trans_lt
        (hlogInert p N' hp hmp hpU hpN hN' hfactor))
    hcircle hz hmono hparam hmem

/-- Prime-scaled bound with the remaining prime-dependent hypotheses reduced to stack endpoint
checks: positivity at the upper endpoint and log-subcriticality at the lower endpoint. -/
theorem card_le_eight_mul_log_div_loglog_of_radius_sq_endpoint_bounds
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (honeFull : 1 ≤ C ^ 2 * R)
    (honeUpper :
      1 ≤ C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ))
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hsplitLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)))
    (hinertLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hL : 0 ≤ L)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 8 * B * (Real.log R / Real.log (Real.log R)) := by
  have hm₀ : 0 < m₀ := by omega
  have hCR_nonneg : 0 ≤ C ^ 2 * R := by
    exact le_of_lt (mul_pos (sq_pos_of_pos hC) hR)
  exact card_le_eight_mul_log_div_loglog_of_radius_sq_real_scale_bounds
    (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
    (C := C) (R := R) (B := B) (a := a) (L := L) (z := z) (t := t)
    hR2 honeFull
    (fun p hp _hmp hpU _hpN => by
      have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
      have hdiv :
          C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ) ≤
            C ^ 2 * R / (p : ℝ) :=
        div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpU)
      exact honeUpper.trans hdiv)
    hthreshold hk h2m₀ hM hC hR hL hLbound hsmallPrime hN hs herr
    (fun p N' hp hmp _hpU _hpN _hN' _hfactor =>
      split_real_scale_log_condition_of_lower_endpoint
        (s := s) (m₀ := m₀) (p := p) hC hR hm₀ hp hmp hsplitLower)
    (fun p N' hp hmp _hpU _hpN _hN' _hfactor =>
      inert_real_scale_log_condition_of_lower_endpoint
        (s := s) (m₀ := m₀) (p := p) hC hR hm₀ hp hmp hinertLower)
    hcircle hz hmono hparam hmem

/-- Correct two-scale endpoint version: split descent uses `C * sqrt R / sqrt p`, while inert
descent uses `C * sqrt R / p`. -/
theorem card_le_eight_mul_log_div_loglog_of_radius_sq_twoScale_endpoint_bounds
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (honeFull : 1 ≤ C ^ 2 * R)
    (honeSplitUpper :
      1 ≤ C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ))
    (honeInertUpper :
      1 ≤ C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ) ^ 2)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hsplitLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)))
    (hinertLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ 8 * B * (Real.log R / Real.log (Real.log R)) := by
  let U := IntervalStack.geomLower m₀ k
  have hm₀ : 0 < m₀ := by omega
  have hU_one : 1 ≤ U := by
    dsimp [U]
    unfold IntervalStack.geomLower
    have hm₀_one : 1 ≤ m₀ := by omega
    have hpow : 1 ≤ 8 ^ k := Nat.one_le_pow k 8 (by norm_num)
    simpa using Nat.mul_le_mul hm₀_one hpow
  have hU_pos : 0 < U := Nat.pos_of_ne_zero (by omega)
  have hU_pos_real : 0 < (U : ℝ) := by exact_mod_cast hU_pos
  have hCR_pos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hCR_nonneg : 0 ≤ C ^ 2 * R := hCR_pos.le
  have hlogN := log_nat_eq_two_log_radius_of_sq (N := N) hR hR2
  have hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2) :=
    fullRadiusScale_floor_pos_of_one_le hR.le honeFull
  have hlogFull :
      Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤
        Real.log (C ^ 2 * R) :=
    fullRadiusScale_log_floor_le_log_of_sq_le hR.le hfloorFull le_rfl
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  have hthresholdFloor :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((U + 1 : ℕ) : ℝ) * Real.log (U : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    have hthresholdFloorR :
        ((s * (2 * s + 1) : ℕ) : ℝ) *
            Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
          ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
              2 * (((U + 1 : ℕ) : ℝ) * Real.log (U : ℝ))) +
            ((s * s : ℕ) : ℝ) * (2 * Real.log R) :=
      (mul_le_mul_of_nonneg_left hlogFull hcoeff_nonneg).trans_lt (by simpa [U] using hthreshold)
    simpa [hlogN] using hthresholdFloorR
  have hdichotomy :=
    MainDichotomy.missing_or_few_stack_dichotomy_of_crude_threshold
      (s := s) (N := N) (m₀ := m₀) (k := k) (A := fullRadiusScale C R)
      (by simpa [U] using hU_one) (by simpa [U] using hthresholdFloor)
  have hfinite :=
    MainDichotomy.card_le_two_mul_succ_q_mul_two_s_of_interval_stack_twoScale_log
      (M := M) (s := s) (q := 1) (N := N) (m₀ := m₀) (k := k)
      (a := a) (L := L) (A₀ := fullRadiusScale C R)
      (Asplit := fun p => radiusScale C R p)
      (Ainert := fun p => inertRadiusScale C R p) (z := z) (t := t)
      hdichotomy hk hm₀ h2m₀ hM (fullRadiusScale_pos hC hR) hfloorFull
      hsmallPrime hN
      (by simpa using fullRadiusScale_term_le (s := s) hC hR hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using inertRadiusScale_inert_term_le (s := s) (p := p) hC hR hp.pos hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using primeRadiusScale_split_term_le (s := s) (p := p) hC hR hp.pos hLbound)
      herr
      (fun p hp _hmp _hpU _hpN => radiusScale_pos hC hR hp.pos)
      (fun p hp _hmp hpU _hpN => by
        have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
        have hdiv :
            C ^ 2 * R / (U : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
          div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpU)
        have honeSplitU : 1 ≤ C ^ 2 * R / (U : ℝ) := by
          simpa [U] using honeSplitUpper
        exact radiusScale_floor_pos_of_one_le hR.le hp.pos
          (honeSplitU.trans hdiv))
      (fun p hp _hmp _hpU _hpN => inertRadiusScale_pos hC hR hp.pos)
      (fun p hp _hmp hpU _hpN => by
        have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
        have hpUreal : (p : ℝ) ≤ (U : ℝ) := by exact_mod_cast hpU
        have hpU_sq : (p : ℝ) ^ 2 ≤ (U : ℝ) ^ 2 :=
          pow_le_pow_left₀ hp_pos_real.le hpUreal 2
        have hdiv :
            C ^ 2 * R / (U : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
          div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpU_sq
        have honeInertU : 1 ≤ C ^ 2 * R / (U : ℝ) ^ 2 := by
          simpa [U] using honeInertUpper
        exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
          (honeInertU.trans hdiv))
      (fun p N' hp hmp hpU hpN hN' hfactor => by
        have hfloorP : 0 < Nat.floor ((radiusScale C R p) ^ 2) := by
          have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
          have hdiv :
              C ^ 2 * R / (U : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
            div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpU)
          have honeSplitU : 1 ≤ C ^ 2 * R / (U : ℝ) := by
            simpa [U] using honeSplitUpper
          exact radiusScale_floor_pos_of_one_le hR.le hp.pos
            (honeSplitU.trans hdiv)
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
      (fun p N' hp hmp hpU hpN hN' hfactor => by
        have hfloorP : 0 < Nat.floor ((inertRadiusScale C R p) ^ 2) := by
          have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
          have hpUreal : (p : ℝ) ≤ (U : ℝ) := by exact_mod_cast hpU
          have hpU_sq : (p : ℝ) ^ 2 ≤ (U : ℝ) ^ 2 :=
            pow_le_pow_left₀ hp_pos_real.le hpUreal 2
          have hdiv :
              C ^ 2 * R / (U : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
            div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpU_sq
          have honeInertU : 1 ≤ C ^ 2 * R / (U : ℝ) ^ 2 := by
            simpa [U] using honeInertUpper
          exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
            (honeInertU.trans hdiv)
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
  calc
    (M : ℝ) ≤ 2 * (((1 + 1 : ℕ) : ℝ) * ((2 * s : ℕ) : ℝ)) := hfinite
    _ ≤ (8 * s : ℕ) := by
      norm_num
      ring_nf
      exact le_rfl
    _ = 8 * (s : ℝ) := by norm_num
    _ ≤ 8 * (B * (Real.log R / Real.log (Real.log R))) := by
      nlinarith [hs]
    _ = 8 * B * (Real.log R / Real.log (Real.log R)) := by ring

/-- Sharper two-scale endpoint version with a constant-enlarged split scale.

The split descent branch uses `48 * C * sqrt R / sqrt p` as its subcritical scale.
This keeps the determinant endpoint inequality asymptotically valid while making the split
sliding-window length contribution only `1/48` of the base window size.  The finite branch
factor is therefore `4 + 4/48 = 49/12` instead of `8`. -/
theorem card_le_49_div_12_mul_log_div_loglog_of_radius_sq_twoScale_endpoint_bounds
    {M s N m₀ k : ℕ}
    {C R B a L : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hR2 : R ^ 2 = (N : ℝ))
    (honeFull : 1 ≤ C ^ 2 * R)
    (honeSplitUpper :
      1 ≤ C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ))
    (honeInertUpper :
      1 ≤ C ^ 2 * R / (IntervalStack.geomLower m₀ k : ℝ) ^ 2)
    (hthreshold :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (C ^ 2 * R) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((IntervalStack.geomLower m₀ k + 1 : ℕ) : ℝ) *
              Real.log (IntervalStack.geomLower m₀ k : ℝ))) +
          ((s * s : ℕ) : ℝ) * (2 * Real.log R))
    (hsplitScaledLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log ((48 * C) ^ 2 * R / (m₀ : ℝ)) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - Real.log (m₀ : ℝ)))
    (hinertLower :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (C ^ 2 * R / (m₀ : ℝ) ^ 2) <
        ((s * s : ℕ) : ℝ) * (2 * Real.log R - 2 * Real.log (m₀ : ℝ)))
    (hk : 0 < k)
    (h2m₀ : 2 ≤ m₀)
    (hM : 0 < M)
    (hC : 0 < C)
    (hR : 0 < R)
    (hLbound : L ≤ C * Real.sqrt R)
    (hsmallPrime : 4 * IntervalStack.geomLower m₀ k ≤ s)
    (hN : 0 < N)
    (hs :
      (s : ℝ) ≤ B * (Real.log R / Real.log (Real.log R)))
    (herr : ∀ i : Fin k,
      Real.log (((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) + 1) +
          2 * Real.sqrt ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) *
            Real.log ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ) ≤
        (Real.log 2 / 4) * ((IntervalStack.geomUpper m₀ i : ℕ) : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ (49 / 12 : ℝ) * B * (Real.log R / Real.log (Real.log R)) := by
  let U := IntervalStack.geomLower m₀ k
  have hm₀ : 0 < m₀ := by omega
  have hU_one : 1 ≤ U := by
    dsimp [U]
    unfold IntervalStack.geomLower
    have hm₀_one : 1 ≤ m₀ := by omega
    have hpow : 1 ≤ 8 ^ k := Nat.one_le_pow k 8 (by norm_num)
    simpa using Nat.mul_le_mul hm₀_one hpow
  have hU_pos : 0 < U := Nat.pos_of_ne_zero (by omega)
  have hCR_pos : 0 < C ^ 2 * R := mul_pos (sq_pos_of_pos hC) hR
  have hCR_nonneg : 0 ≤ C ^ 2 * R := hCR_pos.le
  have hlogN := log_nat_eq_two_log_radius_of_sq (N := N) hR hR2
  have hfloorFull : 0 < Nat.floor ((fullRadiusScale C R) ^ 2) :=
    fullRadiusScale_floor_pos_of_one_le hR.le honeFull
  have hlogFull :
      Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) ≤
        Real.log (C ^ 2 * R) :=
    fullRadiusScale_log_floor_le_log_of_sq_le hR.le hfloorFull le_rfl
  have hcoeff_nonneg : 0 ≤ ((s * (2 * s + 1) : ℕ) : ℝ) := by positivity
  have hthresholdFloor :
      ((s * (2 * s + 1) : ℕ) : ℝ) *
          Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
            2 * (((U + 1 : ℕ) : ℝ) * Real.log (U : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    have hthresholdFloorR :
        ((s * (2 * s + 1) : ℕ) : ℝ) *
            Real.log (Nat.floor ((fullRadiusScale C R) ^ 2) : ℝ) <
          ((s : ℝ) ^ 2 * ((k : ℝ) * (Real.log 2 / 2)) -
              2 * (((U + 1 : ℕ) : ℝ) * Real.log (U : ℝ))) +
            ((s * s : ℕ) : ℝ) * (2 * Real.log R) :=
      (mul_le_mul_of_nonneg_left hlogFull hcoeff_nonneg).trans_lt
        (by simpa [U] using hthreshold)
    simpa [hlogN] using hthresholdFloorR
  have hdichotomy :=
    MainDichotomy.missing_or_few_stack_dichotomy_of_crude_threshold
      (s := s) (N := N) (m₀ := m₀) (k := k) (A := fullRadiusScale C R)
      (by simpa [U] using hU_one) (by simpa [U] using hthresholdFloor)
  have hfinite :=
    MainDichotomy.card_le_four_add_split_eps_mul_s_of_interval_stack_twoScale_log
      (M := M) (s := s) (N := N) (m₀ := m₀) (k := k)
      (a := a) (L := L) (A₀ := fullRadiusScale C R) (eps := (1 / 48 : ℝ))
      (Asplit := fun p => 48 * radiusScale C R p)
      (Ainert := fun p => inertRadiusScale C R p) (z := z) (t := t)
      hdichotomy hk hm₀ h2m₀ hM (fullRadiusScale_pos hC hR) hfloorFull
      hsmallPrime hN (by norm_num)
      (by simpa using fullRadiusScale_term_le (s := s) hC hR hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using inertRadiusScale_inert_term_le (s := s) (p := p) hC hR hp.pos hLbound)
      (fun p hp _hmp _hpU _hpN => by
        simpa using primeRadiusScale_split_term_le_div_const
          (s := s) (p := p) (K := (48 : ℝ)) hC hR hp.pos (by norm_num) hLbound)
      herr
      (fun p hp _hmp _hpU _hpN => by
        exact mul_pos (by norm_num) (radiusScale_pos hC hR hp.pos))
      (fun p hp _hmp hpU _hpN => by
        have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
        have hdiv :
            C ^ 2 * R / (U : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
          div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpU)
        have honeSplitU : 1 ≤ C ^ 2 * R / (U : ℝ) := by
          simpa [U] using honeSplitUpper
        have hscale_sq :
            (48 * radiusScale C R p) ^ 2 = (48 * C) ^ 2 * R / (p : ℝ) := by
          rw [mul_pow, radiusScale_sq hR.le hp.pos]
          ring
        have hscale_ge_one : 1 ≤ (48 * C) ^ 2 * R / (p : ℝ) := by
          have honeP : 1 ≤ C ^ 2 * R / (p : ℝ) := honeSplitU.trans hdiv
          have heq :
              (48 * C) ^ 2 * R / (p : ℝ) =
                (48 : ℝ) ^ 2 * (C ^ 2 * R / (p : ℝ)) := by
            ring
          rw [heq]
          nlinarith
        apply Nat.floor_pos.mpr
        rw [hscale_sq]
        exact hscale_ge_one)
      (fun p hp _hmp _hpU _hpN => inertRadiusScale_pos hC hR hp.pos)
      (fun p hp _hmp hpU _hpN => by
        have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
        have hpUreal : (p : ℝ) ≤ (U : ℝ) := by exact_mod_cast hpU
        have hpU_sq : (p : ℝ) ^ 2 ≤ (U : ℝ) ^ 2 :=
          pow_le_pow_left₀ hp_pos_real.le hpUreal 2
        have hdiv :
            C ^ 2 * R / (U : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
          div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpU_sq
        have honeInertU : 1 ≤ C ^ 2 * R / (U : ℝ) ^ 2 := by
          simpa [U] using honeInertUpper
        exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
          (honeInertU.trans hdiv))
      (fun p N' hp hmp hpU _hpN _hN' _hfactor => by
        have hfloorP : 0 < Nat.floor ((48 * radiusScale C R p) ^ 2) := by
          have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
          have hdiv :
              C ^ 2 * R / (U : ℝ) ≤ C ^ 2 * R / (p : ℝ) :=
            div_le_div_of_nonneg_left hCR_nonneg hp_pos_real (by exact_mod_cast hpU)
          have honeSplitU : 1 ≤ C ^ 2 * R / (U : ℝ) := by
            simpa [U] using honeSplitUpper
          have hscale_sq :
              (48 * radiusScale C R p) ^ 2 = (48 * C) ^ 2 * R / (p : ℝ) := by
            rw [mul_pow, radiusScale_sq hR.le hp.pos]
            ring
          have hscale_ge_one : 1 ≤ (48 * C) ^ 2 * R / (p : ℝ) := by
            have honeP : 1 ≤ C ^ 2 * R / (p : ℝ) := honeSplitU.trans hdiv
            have heq :
                (48 * C) ^ 2 * R / (p : ℝ) =
                  (48 : ℝ) ^ 2 * (C ^ 2 * R / (p : ℝ)) := by
              ring
            rw [heq]
            nlinarith
          apply Nat.floor_pos.mpr
          rw [hscale_sq]
          exact hscale_ge_one
        have hscale_sq :
            (48 * radiusScale C R p) ^ 2 = (48 * C) ^ 2 * R / (p : ℝ) := by
          rw [mul_pow, radiusScale_sq hR.le hp.pos]
          ring
        have hlogP :
            Real.log (Nat.floor ((48 * radiusScale C R p) ^ 2) : ℝ) ≤
              Real.log ((48 * C) ^ 2 * R / (p : ℝ)) :=
          log_natFloor_le_log_of_le (sq_nonneg _) hfloorP (by rw [hscale_sq])
        have hC48 : 0 < 48 * C := mul_pos (by norm_num) hC
        have hsplit :=
          split_real_scale_log_condition_of_lower_endpoint
            (s := s) (m₀ := m₀) (p := p) (C := 48 * C)
            hC48 hR hm₀ hp hmp hsplitScaledLower
        have hsplitN :
            ((s * (2 * s + 1) : ℕ) : ℝ) *
                Real.log ((48 * C) ^ 2 * R / (p : ℝ)) <
              ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ)) := by
          simpa [hlogN] using hsplit
        exact (mul_le_mul_of_nonneg_left hlogP hcoeff_nonneg).trans_lt hsplitN)
      (fun p N' hp hmp hpU _hpN _hN' _hfactor => by
        have hfloorP : 0 < Nat.floor ((inertRadiusScale C R p) ^ 2) := by
          have hp_pos_real : 0 < (p : ℝ) := by exact_mod_cast hp.pos
          have hpUreal : (p : ℝ) ≤ (U : ℝ) := by exact_mod_cast hpU
          have hpU_sq : (p : ℝ) ^ 2 ≤ (U : ℝ) ^ 2 :=
            pow_le_pow_left₀ hp_pos_real.le hpUreal 2
          have hdiv :
              C ^ 2 * R / (U : ℝ) ^ 2 ≤ C ^ 2 * R / (p : ℝ) ^ 2 :=
            div_le_div_of_nonneg_left hCR_nonneg (sq_pos_of_pos hp_pos_real) hpU_sq
          have honeInertU : 1 ≤ C ^ 2 * R / (U : ℝ) ^ 2 := by
            simpa [U] using honeInertUpper
          exact inertRadiusScale_floor_pos_of_one_le hR.le hp.pos
            (honeInertU.trans hdiv)
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
  calc
    (M : ℝ) ≤ (4 + 4 * (1 / 48 : ℝ)) * (s : ℝ) := hfinite
    _ = (49 / 12 : ℝ) * (s : ℝ) := by norm_num
    _ ≤ (49 / 12 : ℝ) * (B * (Real.log R / Real.log (Real.log R))) := by
      exact mul_le_mul_of_nonneg_left hs (by norm_num)
    _ = (49 / 12 : ℝ) * B * (Real.log R / Real.log (Real.log R)) := by ring

end MainSublogBound
end GaussianChain
