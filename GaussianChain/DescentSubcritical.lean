import GaussianChain.DescentGeometry

namespace GaussianChain
namespace DescentSubcritical

open SubcriticalBound
open DescentGeometry
open scoped ComplexConjugate

/-- Inert-prime descent followed by the ordinary subcritical window bound.

Dividing by an inert rational prime decreases squared chord lengths by `p^2`; hence any
arclength-parameter chord bound valid before descent remains valid after descent. -/
theorem card_le_of_param_subcritical_windows_after_inert_descent
    {M s N N' p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ))
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
  let zFin : Fin M → GaussianInt := fun i => z i
  have hzNorm : ∀ i : Fin M, (zFin i).norm = (N : ℤ) := by
    intro i
    exact RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle i i.isLt)
  obtain ⟨wFin, hwFin⟩ :=
    PrimeDescent.exists_inert_quotient_family_norm_eq_of_common_norm_factor
      (ι := Fin M) hp3 hfactor zFin hzNorm
  let w : ℕ → GaussianInt := fun i => if hi : i < M then wFin ⟨i, hi⟩ else 0
  have hwfac : ∀ i, i < M → z i = (p : GaussianInt) * w i := by
    intro i hi
    simpa [w, zFin, hi] using (hwFin ⟨i, hi⟩).1
  have hwnorm : ∀ i, i < M → (w i).norm = (N' : ℤ) := by
    intro i hi
    simpa [w, hi] using (hwFin ⟨i, hi⟩).2
  have hcircle_w : OnCircleUpTo M N' w := by
    intro i hi
    rw [← Zsqrtd.norm_eq_mul_conj, hwnorm i hi]
  have hw_inj : InjectiveUpTo M w := by
    intro i j hi hj hij
    apply hz hi hj
    rw [hwfac i hi, hwfac j hj, hij]
  have hp_sq_one : (1 : ℝ) ≤ (p : ℝ) ^ 2 := by
    have hp_one : (1 : ℝ) ≤ (p : ℝ) := by
      exact_mod_cast (Fact.out : Nat.Prime p).one_le
    nlinarith [sq_nonneg ((p : ℝ) - 1)]
  have hparam_w :
      ∀ i j, i < M → j < M → gaussianSqDist (w i) (w j) ≤ (t j - t i) ^ 2 := by
    intro i j hi hj
    rw [gaussianSqDist_quotient_eq_div_of_natCast_mul_left (hwfac i hi) (hwfac j hj)]
    exact (div_le_self (sq_nonneg _) hp_sq_one).trans'
      (div_le_div_of_nonneg_right (hparam i j hi hj) (sq_nonneg _))
  exact card_le_of_param_subcritical_windows (M := M) (s := s) (N := N')
    (a := a) (L := L) (B := B) (z := w) (t := t)
    hk hB hN' hsmall hcircle_w hw_inj hmono hparam_w hmem

/-- Inert-prime descent with the arclength parameter scaled by `1 / p`.

This is the sharper form used in the final descent branch: after dividing every point by the
rational Gaussian integer `p`, the descended arc length is `L / p`. -/
theorem card_le_of_param_subcritical_windows_after_inert_descent_scaled
    {M s N N' p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ))
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * (L / (p : ℝ)) / B := by
  let zFin : Fin M → GaussianInt := fun i => z i
  have hzNorm : ∀ i : Fin M, (zFin i).norm = (N : ℤ) := by
    intro i
    exact RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle i i.isLt)
  obtain ⟨wFin, hwFin⟩ :=
    PrimeDescent.exists_inert_quotient_family_norm_eq_of_common_norm_factor
      (ι := Fin M) hp3 hfactor zFin hzNorm
  let w : ℕ → GaussianInt := fun i => if hi : i < M then wFin ⟨i, hi⟩ else 0
  have hwfac : ∀ i, i < M → z i = (p : GaussianInt) * w i := by
    intro i hi
    simpa [w, zFin, hi] using (hwFin ⟨i, hi⟩).1
  have hwnorm : ∀ i, i < M → (w i).norm = (N' : ℤ) := by
    intro i hi
    simpa [w, hi] using (hwFin ⟨i, hi⟩).2
  have hcircle_w : OnCircleUpTo M N' w := by
    intro i hi
    rw [← Zsqrtd.norm_eq_mul_conj, hwnorm i hi]
  have hw_inj : InjectiveUpTo M w := by
    intro i j hi hj hij
    apply hz hi hj
    rw [hwfac i hi, hwfac j hj, hij]
  let tdesc : ℕ → ℝ := fun i => t i / (p : ℝ)
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).pos
  have hp_sq_nonneg : 0 ≤ (p : ℝ) ^ 2 := sq_nonneg _
  have hmono_desc : ∀ i j, i ≤ j → j < M → tdesc i ≤ tdesc j := by
    intro i j hij hj
    exact div_le_div_of_nonneg_right (hmono i j hij hj) hp_pos.le
  have hparam_desc :
      ∀ i j, i < M → j < M → gaussianSqDist (w i) (w j) ≤ (tdesc j - tdesc i) ^ 2 := by
    intro i j hi hj
    have hquot :
        gaussianSqDist (w i) (w j) =
          gaussianSqDist (z i) (z j) / (p : ℝ) ^ 2 :=
      gaussianSqDist_quotient_eq_div_of_natCast_mul_left (hwfac i hi) (hwfac j hj)
    have hsq :
        (tdesc j - tdesc i) ^ 2 = (t j - t i) ^ 2 / (p : ℝ) ^ 2 := by
      dsimp [tdesc]
      field_simp [hp_pos.ne']
    rw [hquot, hsq]
    exact div_le_div_of_nonneg_right (hparam i j hi hj) hp_sq_nonneg
  have hmem_desc : ∀ i, i < M → a / (p : ℝ) ≤ tdesc i ∧
      tdesc i ≤ a / (p : ℝ) + L / (p : ℝ) := by
    intro i hi
    rcases hmem i hi with ⟨hlo, hhi⟩
    constructor
    · exact div_le_div_of_nonneg_right hlo hp_pos.le
    · have h := div_le_div_of_nonneg_right hhi hp_pos.le
      simpa [tdesc, add_div] using h
  exact card_le_of_param_subcritical_windows (M := M) (s := s) (N := N')
    (a := a / (p : ℝ)) (L := L / (p : ℝ)) (B := B) (z := w) (t := tdesc)
    hk hB hN' hsmall hcircle_w hw_inj hmono_desc hparam_desc hmem_desc

/-- Split-prime descent followed by the ordinary subcritical window bound, for an already
ordered descended subfamily.

This is the analytic/combinatorial bridge needed after the split-prime pigeonhole step: if at
least half of the original `M` points are listed in an ordered subfamily of length `Md`, and that
subfamily descends through a split Gaussian factor to a subcritical circle, then the original
cardinality is at most twice the corresponding subcritical bound. -/
theorem card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence
    {M Md s N' p : ℕ} [Fact p.Prime] {ρ : GaussianInt}
    (hρ : ρ * star ρ = (p : GaussianInt))
    {a L B : ℝ} {z w : ℕ → GaussianInt} {t : ℕ → ℝ} {idx : ℕ → Fin M}
    (hmany : M ≤ 2 * Md)
    (hk : 2 * s ≤ Md)
    (hB : 0 < B)
    (hN' : 0 < N')
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hfac : ∀ i, i < Md → z (idx i) = ρ * w i)
    (hwnorm : ∀ i, (w i).norm = (N' : ℤ))
    (hwinj : ∀ j, j < Md - 2 * s →
      Function.Injective fun i : Fin (2 * s + 1) => w (j + (i : ℕ)))
    (hidx_mono : ∀ i j, i ≤ j → j < Md → (idx i : ℕ) ≤ (idx j : ℕ))
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      2 * (((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B) := by
  let tsub : ℕ → ℝ := fun i => t (idx i)
  have hcircle_w : OnCircleUpTo Md N' w := by
    intro n _hn
    rw [← Zsqrtd.norm_eq_mul_conj, hwnorm n]
  have hmono_sub : ∀ i j, i ≤ j → j < Md → tsub i ≤ tsub j := by
    intro i j hij hj
    exact hmono (idx i) (idx j) (hidx_mono i j hij hj) (idx j).isLt
  have hp_one : (1 : ℝ) ≤ (p : ℝ) := by
    exact_mod_cast (Fact.out : Nat.Prime p).one_le
  have hparam_sub :
      ∀ i j, i < Md → j < Md → gaussianSqDist (w i) (w j) ≤ (tsub j - tsub i) ^ 2 := by
    intro i j hi hj
    have hquot :
        gaussianSqDist (w i) (w j) =
          gaussianSqDist (z (idx i)) (z (idx j)) / (p : ℝ) :=
      gaussianSqDist_quotient_eq_div_of_mul_left hρ (hfac i hi) (hfac j hj)
    rw [hquot]
    exact (div_le_self (sq_nonneg _) hp_one).trans'
      (div_le_div_of_nonneg_right (hparam (idx i) (idx j) (idx i).isLt (idx j).isLt)
        (by exact_mod_cast (Fact.out : Nat.Prime p).pos.le))
  have hmem_sub : ∀ i, i < Md → a ≤ tsub i ∧ tsub i ≤ a + L := by
    intro i _hi
    exact hmem (idx i) (idx i).isLt
  have hMd_bound :
      (Md : ℝ) ≤ ((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) * L / B := by
    refine card_le_of_subcritical_windows (M := Md) (s := s) (N := N')
      (a := a) (L := L) (K := B ^ 2) (B := B) (z := w) (t := tsub)
      hk hB hN' hsmall hcircle_w hwinj ?_ hmem_sub
    intro j hj hspan i k
    exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono_sub hparam_sub hj hspan i k
  have hMle : (M : ℝ) ≤ 2 * (Md : ℝ) := by
    exact_mod_cast hmany
  exact hMle.trans (mul_le_mul_of_nonneg_left hMd_bound (by norm_num))

/-- Split-prime descent with the arclength parameter scaled by `1 / sqrt p`, for an already
ordered descended subfamily. -/
theorem card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence_scaled
    {M Md s N' p : ℕ} [Fact p.Prime] {ρ : GaussianInt}
    (hρ : ρ * star ρ = (p : GaussianInt))
    {a L B : ℝ} {z w : ℕ → GaussianInt} {t : ℕ → ℝ} {idx : ℕ → Fin M}
    (hmany : M ≤ 2 * Md)
    (hk : 2 * s ≤ Md)
    (hB : 0 < B)
    (hN' : 0 < N')
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hfac : ∀ i, i < Md → z (idx i) = ρ * w i)
    (hwnorm : ∀ i, (w i).norm = (N' : ℤ))
    (hwinj : ∀ j, j < Md - 2 * s →
      Function.Injective fun i : Fin (2 * s + 1) => w (j + (i : ℕ)))
    (hidx_mono : ∀ i j, i ≤ j → j < Md → (idx i : ℕ) ≤ (idx j : ℕ))
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      2 * (((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) *
        (L / Real.sqrt (p : ℝ)) / B) := by
  let tsub : ℕ → ℝ := fun i => t (idx i) / Real.sqrt (p : ℝ)
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).pos
  have hsqrt_pos : 0 < Real.sqrt (p : ℝ) := Real.sqrt_pos_of_pos hp_pos
  have hsqrt_ne : Real.sqrt (p : ℝ) ≠ 0 := ne_of_gt hsqrt_pos
  have hcircle_w : OnCircleUpTo Md N' w := by
    intro n _hn
    rw [← Zsqrtd.norm_eq_mul_conj, hwnorm n]
  have hmono_sub : ∀ i j, i ≤ j → j < Md → tsub i ≤ tsub j := by
    intro i j hij hj
    exact div_le_div_of_nonneg_right
      (hmono (idx i) (idx j) (hidx_mono i j hij hj) (idx j).isLt) hsqrt_pos.le
  have hparam_sub :
      ∀ i j, i < Md → j < Md → gaussianSqDist (w i) (w j) ≤ (tsub j - tsub i) ^ 2 := by
    intro i j hi hj
    have hquot :
        gaussianSqDist (w i) (w j) =
          gaussianSqDist (z (idx i)) (z (idx j)) / (p : ℝ) :=
      gaussianSqDist_quotient_eq_div_of_mul_left hρ (hfac i hi) (hfac j hj)
    have hsq :
        (tsub j - tsub i) ^ 2 = (t (idx j) - t (idx i)) ^ 2 / (p : ℝ) := by
      dsimp [tsub]
      field_simp [hsqrt_ne, hp_pos.ne']
      rw [Real.sq_sqrt hp_pos.le]
      ring
    rw [hquot, hsq]
    exact div_le_div_of_nonneg_right
      (hparam (idx i) (idx j) (idx i).isLt (idx j).isLt) hp_pos.le
  have hmem_sub : ∀ i, i < Md → a / Real.sqrt (p : ℝ) ≤ tsub i ∧
      tsub i ≤ a / Real.sqrt (p : ℝ) + L / Real.sqrt (p : ℝ) := by
    intro i _hi
    rcases hmem (idx i) (idx i).isLt with ⟨hlo, hhi⟩
    constructor
    · exact div_le_div_of_nonneg_right hlo hsqrt_pos.le
    · have h := div_le_div_of_nonneg_right hhi hsqrt_pos.le
      simpa [tsub, add_div] using h
  have hMd_bound :
      (Md : ℝ) ≤ ((2 * s : ℕ) : ℝ) +
        ((2 * s : ℕ) : ℝ) * (L / Real.sqrt (p : ℝ)) / B := by
    refine card_le_of_subcritical_windows (M := Md) (s := s) (N := N')
      (a := a / Real.sqrt (p : ℝ)) (L := L / Real.sqrt (p : ℝ))
      (K := B ^ 2) (B := B) (z := w) (t := tsub)
      hk hB hN' hsmall hcircle_w hwinj ?_ hmem_sub
    intro j hj hspan i k
    exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono_sub hparam_sub hj hspan i k
  have hMle : (M : ℝ) ≤ 2 * (Md : ℝ) := by
    exact_mod_cast hmany
  exact hMle.trans (mul_le_mul_of_nonneg_left hMd_bound (by norm_num))

/-- Split-prime descent with automatic choice of the ordered half-subfamily. -/
theorem card_le_two_mul_param_subcritical_bound_after_split_descent_scaled
    {M s N N' p : ℕ} [Fact p.Prime] (hp1 : p % 4 = 1)
    {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hmany_for_window : 4 * s ≤ M)
    (hB : 0 < B)
    (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) * (N' : ℤ))
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (M : ℝ) ≤
      2 * (((2 * s : ℕ) : ℝ) + ((2 * s : ℕ) : ℝ) *
        (L / Real.sqrt (p : ℝ)) / B) := by
  let zFin : Fin M → GaussianInt := fun i => z i
  have hzNorm : ∀ i : Fin M, (zFin i).norm = (N : ℤ) := by
    intro i
    exact RamanaDeterminant.norm_eq_int_of_mul_star_eq (hcircle i i.isLt)
  obtain ⟨Md, ρ, idx, w, hρ, hmany, hfac, hwnorm, hidx_mono, hwinj_all⟩ :=
    PrimeDescent.exists_split_descended_ordered_subsequence_norm_eq
      hp1 hfactor zFin hzNorm hM
  have hk : 2 * s ≤ Md := by omega
  have hzFin_inj : Function.Injective zFin := by
    intro i j hij
    apply Fin.ext
    exact hz i.isLt j.isLt hij
  have hfac' : ∀ i, i < Md → z (idx i) = ρ * w i := by
    intro i hi
    simpa [zFin] using hfac i hi
  exact card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence_scaled
    (M := M) (Md := Md) (s := s) (N' := N') (p := p) (ρ := ρ)
    hρ hmany hk hB hN' hsmall hfac' hwnorm (hwinj_all hzFin_inj s) hidx_mono
    hmono hparam hmem

end DescentSubcritical
end GaussianChain
