import GaussianChain.MissingPrimeSubcritical
import GaussianChain.DescentSubcritical

/-!
Saturated versions of the sliding-window cardinality bounds.

When the window span threshold `B` is at least the interval length `L`, no complete window of
`2 * s + 1` points fits inside the interval at all, so the point count is at most `2 * s`, with
no `2 * s * L / B` length term.
-/

namespace GaussianChain

namespace ArcCombinatorics

/-- Saturated window count: if the span threshold `B` is at least the interval length `L`, then
no complete `k + 1`-point window fits inside the interval, so there are at most `k` points. -/
theorem card_le_of_window_span_saturated
    {M k : ℕ} {a L B : ℝ} (t : ℕ → ℝ)
    (hLB : L ≤ B)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L)
    (hspan : ∀ j, j < M - k → B < t (j + k) - t j) :
    M ≤ k := by
  by_contra hMk
  have hk_lt : k < M := by omega
  have hspan0 : B < t k - t 0 := by simpa using hspan 0 (by omega)
  have hmemk := hmem k hk_lt
  have hmem0 := hmem 0 (by omega)
  linarith [hmemk.2, hmem0.1]

end ArcCombinatorics

namespace SubcriticalBound

set_option linter.unusedVariables false in
/-- Saturated form of `card_le_of_param_subcritical_windows`: when the span threshold `B` is at
least the arc length `L`, no complete `2 * s + 1`-point window fits inside the arc, so at most
`2 * s` points occur. The hypothesis `hk` is kept for call-site uniformity. -/
theorem card_le_of_param_subcritical_windows_saturated
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hLB : L ≤ B)
    (hN : 0 < N)
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 2 * s := by
  refine ArcCombinatorics.card_le_of_window_span_saturated (M := M) (k := 2 * s) t hLB hmem ?_
  refine window_span_gt_of_subcritical hN hsmall hcircle
    (fun _j hj => nat_window_injective_of_injectiveUpTo hz hj) ?_
  intro j hj hspan i k
  exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono hparam hj hspan i k

end SubcriticalBound

namespace MissingPrimeSubcritical

open SubcriticalBound
open LogProductBounds

set_option linter.unusedVariables false in
/-- Saturated form of
`card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound`: when the span
threshold `B` is at least the arc length `L`, at most `2 * s` points occur. The hypothesis `hk`
is kept for call-site uniformity. -/
theorem card_le_of_param_missing_prime_subcritical_windows_of_weighted_log_bound_saturated
    {M s N : ℕ} {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ} {P : Finset ℕ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hLB : L ≤ B)
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hprime : ∀ p ∈ P, Nat.Prime p)
    (hsmallPrime : ∀ p ∈ P, 4 * p ≤ s)
    (hN : 0 < N)
    (hpN : ∀ p ∈ P, ¬ (p : ℤ) ∣ (N : ℤ))
    (hweightedLog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
            2 * (∑ p ∈ P, Real.log (p : ℝ))) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ p q, p ≤ q → q < M → t p ≤ t q)
    (hparam : ∀ p q, p < M → q < M → gaussianSqDist (z p) (z q) ≤ (t q - t p) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 2 * s := by
  have hlog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        (∑ p ∈ P,
          ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
            Real.log (p : ℝ)) +
          ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
    have hfloorLower :=
      sum_two_floor_half_s_sq_div_mul_log_lower
        (P := P) (s := s)
        (fun p hp => (hprime p hp).pos)
        (fun p hp => (hprime p hp).one_le)
    have hfloorWithNorm :
        ((s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
              2 * (∑ p ∈ P, Real.log (p : ℝ))) +
            ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) ≤
          (∑ p ∈ P,
            ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
              Real.log (p : ℝ)) +
            ((s * s : ℕ) : ℝ) * Real.log (N : ℝ) := by
      simpa [add_comm, add_left_comm, add_assoc] using
        add_le_add_right hfloorLower (((s * s : ℕ) : ℝ) * Real.log (N : ℝ))
    exact hweightedLog.trans_le hfloorWithNorm
  have hsmallE :
      (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) <
        (∏ p ∈ P, p ^ (2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))))) *
          N ^ (s * s) :=
    nat_pow_lt_prod_pow_mul_pow_of_log_lt hfloor hN (fun p hp => (hprime p hp).pos) hlog
  have hsmall :=
    window_natFloor_pow_lt_missing_prime_product_of_exponent_le
      (M := M) (s := s) (N := N) (K := B ^ 2) (z := z)
      (fun p hp => (hprime p hp).one_le)
      (window_floor_half_s_sq_div_le_pairCollisionCount hprime hsmallPrime hcircle)
      hsmallE
  refine ArcCombinatorics.card_le_of_window_span_saturated (M := M) (k := 2 * s) t hLB hmem ?_
  refine window_span_gt_of_missing_prime_subcritical hprime hN hpN hsmall hcircle
    (fun _j hj => nat_window_injective_of_injectiveUpTo hz hj) ?_
  intro j hj hspan i k
  exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono hparam hj hspan i k

end MissingPrimeSubcritical

namespace DescentSubcritical

open SubcriticalBound
open DescentGeometry
open scoped ComplexConjugate

/-- Saturated form of `card_le_of_param_subcritical_windows_after_inert_descent_scaled`: when
the span threshold `B` is at least the descended arc length `L / p`, at most `2 * s` points
occur. -/
theorem card_le_of_param_subcritical_windows_after_inert_descent_scaled_saturated
    {M s N N' p : ℕ} [Fact p.Prime] (hp3 : p % 4 = 3)
    {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hk : 2 * s ≤ M)
    (hB : 0 < B)
    (hLB : L / (p : ℝ) ≤ B)
    (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ))
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 2 * s := by
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
  exact card_le_of_param_subcritical_windows_saturated (M := M) (s := s) (N := N')
    (a := a / (p : ℝ)) (L := L / (p : ℝ)) (B := B) (z := w) (t := tdesc)
    hk hB hLB hN' hsmall hcircle_w hw_inj hmono_desc hparam_desc hmem_desc

set_option linter.unusedVariables false in
/-- Saturated form of
`card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence_scaled`: when the span
threshold `B` is at least the descended arc length `L / √p`, the descended half-subfamily has at
most `2 * s` points, so the original family has at most `4 * s`. The hypothesis `hk` is kept for
call-site uniformity. -/
theorem card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence_scaled_saturated
    {M Md s N' p : ℕ} [Fact p.Prime] {ρ : GaussianInt}
    (hρ : ρ * star ρ = (p : GaussianInt))
    {a L B : ℝ} {z w : ℕ → GaussianInt} {t : ℕ → ℝ} {idx : ℕ → Fin M}
    (hmany : M ≤ 2 * Md)
    (hk : 2 * s ≤ Md)
    (hB : 0 < B)
    (hLB : L / Real.sqrt (p : ℝ) ≤ B)
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
    M ≤ 4 * s := by
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
  have hMd_bound : Md ≤ 2 * s := by
    refine ArcCombinatorics.card_le_of_window_span_saturated (M := Md) (k := 2 * s)
      tsub hLB hmem_sub ?_
    refine window_span_gt_of_subcritical hN' hsmall hcircle_w hwinj ?_
    intro j hj hspan i k
    exact sqDist_le_sq_span_of_param_sq_bound hB.le hmono_sub hparam_sub hj hspan i k
  omega

/-- Saturated split-prime descent with automatic choice of the ordered half-subfamily: when the
span threshold `B` is at least the descended arc length `L / √p`, the original family has at
most `4 * s` points. -/
theorem card_le_two_mul_param_subcritical_bound_after_split_descent_scaled_saturated
    {M s N N' p : ℕ} [Fact p.Prime] (hp1 : p % 4 = 1)
    {a L B : ℝ} {z : ℕ → GaussianInt} {t : ℕ → ℝ}
    (hM : 0 < M)
    (hmany_for_window : 4 * s ≤ M)
    (hB : 0 < B)
    (hLB : L / Real.sqrt (p : ℝ) ≤ B)
    (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) * (N' : ℤ))
    (hsmall : (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s))
    (hcircle : OnCircleUpTo M N z)
    (hz : InjectiveUpTo M z)
    (hmono : ∀ i j, i ≤ j → j < M → t i ≤ t j)
    (hparam : ∀ i j, i < M → j < M → gaussianSqDist (z i) (z j) ≤ (t j - t i) ^ 2)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    M ≤ 4 * s := by
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
  exact
    card_le_two_mul_param_subcritical_bound_after_split_descent_subsequence_scaled_saturated
      (M := M) (Md := Md) (s := s) (N' := N') (p := p) (ρ := ρ)
      hρ hmany hk hB hLB hN' hsmall hfac' hwnorm (hwinj_all hzFin_inj s) hidx_mono
      hmono hparam hmem

end DescentSubcritical

end GaussianChain
