import Mathlib

namespace GaussianChain
namespace LogProductBounds

open Finset

/-- A natural-number inequality follows from a strict logarithmic inequality, with positivity
made explicit. -/
theorem nat_lt_of_log_lt_log {m n : ℕ}
    (hm : 0 < m) (hn : 0 < n)
    (hlog : Real.log (m : ℝ) < Real.log (n : ℝ)) :
    m < n := by
  exact_mod_cast (Real.log_lt_log_iff (by exact_mod_cast hm) (by exact_mod_cast hn)).mp hlog

/-- Logarithm of a finite product of natural-prime-style powers, stated for any positive natural
bases. -/
theorem log_prod_pow_natCast_eq_sum {P : Finset ℕ} {E : ℕ → ℕ}
    (hpos : ∀ p ∈ P, 0 < p) :
    Real.log (∏ p ∈ P, (p : ℝ) ^ E p) =
      ∑ p ∈ P, (E p : ℝ) * Real.log (p : ℝ) := by
  classical
  revert hpos
  refine Finset.induction_on P ?base ?step
  · intro _
    simp
  · intro a s has ih hpos
    have ha_pos_nat : 0 < a := hpos a (Finset.mem_insert_self a s)
    have ha_pos : 0 < (a : ℝ) := by exact_mod_cast ha_pos_nat
    have hspos : ∀ p ∈ s, 0 < p := fun p hp => hpos p (Finset.mem_insert_of_mem hp)
    have hprod_ne : (∏ p ∈ s, (p : ℝ) ^ E p) ≠ 0 := by
      exact ne_of_gt (Finset.prod_pos (fun p hp => pow_pos (by exact_mod_cast hspos p hp) _))
    rw [Finset.prod_insert has, Finset.sum_insert has]
    rw [Real.log_mul (pow_ne_zero _ ha_pos.ne') hprod_ne, Real.log_pow, ih hspos]

/-- Natural-cast version of `log_prod_pow_natCast_eq_sum`. -/
theorem log_nat_prod_pow_eq_sum {P : Finset ℕ} {E : ℕ → ℕ}
    (hpos : ∀ p ∈ P, 0 < p) :
    Real.log ((∏ p ∈ P, p ^ E p : ℕ) : ℝ) =
      ∑ p ∈ P, (E p : ℝ) * Real.log (p : ℝ) := by
  simpa using log_prod_pow_natCast_eq_sum (P := P) (E := E) hpos

/-- Logarithm of a product of prime-power factors times a common natural power. -/
theorem log_nat_prod_pow_mul_pow_eq_sum
    {P : Finset ℕ} {E : ℕ → ℕ} {N y : ℕ}
    (hpos : ∀ p ∈ P, 0 < p) (hN : 0 < N) :
    Real.log (((∏ p ∈ P, p ^ E p) * N ^ y : ℕ) : ℝ) =
      (∑ p ∈ P, (E p : ℝ) * Real.log (p : ℝ)) +
        (y : ℝ) * Real.log (N : ℝ) := by
  have hprod_pos_nat : 0 < ∏ p ∈ P, p ^ E p := by
    exact Finset.prod_pos (fun p hp => Nat.pow_pos (hpos p hp))
  have hprod_ne : ((∏ p ∈ P, p ^ E p : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast ne_of_gt hprod_pos_nat
  have hNpow_ne : ((N ^ y : ℕ) : ℝ) ≠ 0 := by
    exact_mod_cast ne_of_gt (Nat.pow_pos hN)
  rw [Nat.cast_mul, Real.log_mul hprod_ne hNpow_ne]
  rw [show ((N ^ y : ℕ) : ℝ) = (N : ℝ) ^ y by norm_num, Real.log_pow]
  rw [show ((∏ p ∈ P, p ^ E p : ℕ) : ℝ) = ∏ p ∈ P, (p : ℝ) ^ E p by simp]
  rw [log_prod_pow_natCast_eq_sum hpos]

/-- A logarithmic lower bound for a product gives the corresponding natural-number product
threshold. -/
theorem nat_lt_prod_pow_mul_pow_of_log_lt
    {m N y : ℕ} {P : Finset ℕ} {E : ℕ → ℕ}
    (hm : 0 < m) (hN : 0 < N) (hpos : ∀ p ∈ P, 0 < p)
    (hlog :
      Real.log (m : ℝ) <
        (∑ p ∈ P, (E p : ℝ) * Real.log (p : ℝ)) +
          (y : ℝ) * Real.log (N : ℝ)) :
    m < (∏ p ∈ P, p ^ E p) * N ^ y := by
  have hright_pos : 0 < (∏ p ∈ P, p ^ E p) * N ^ y := by
    exact Nat.mul_pos (Finset.prod_pos (fun p hp => Nat.pow_pos (hpos p hp)))
      (Nat.pow_pos hN)
  refine nat_lt_of_log_lt_log hm hright_pos ?_
  rwa [log_nat_prod_pow_mul_pow_eq_sum hpos hN]

/-- A convenient power-left version of `nat_lt_prod_pow_mul_pow_of_log_lt`. -/
theorem nat_pow_lt_prod_pow_mul_pow_of_log_lt
    {B x N y : ℕ} {P : Finset ℕ} {E : ℕ → ℕ}
    (hB : 0 < B) (hN : 0 < N) (hpos : ∀ p ∈ P, 0 < p)
    (hlog :
      (x : ℝ) * Real.log (B : ℝ) <
        (∑ p ∈ P, (E p : ℝ) * Real.log (p : ℝ)) +
          (y : ℝ) * Real.log (N : ℝ)) :
    B ^ x < (∏ p ∈ P, p ^ E p) * N ^ y := by
  have hleft_pos : 0 < B ^ x := Nat.pow_pos hB
  refine nat_lt_prod_pow_mul_pow_of_log_lt hleft_pos hN hpos ?_
  rwa [show ((B ^ x : ℕ) : ℝ) = (B : ℝ) ^ x by norm_num, Real.log_pow]

/-- The floor exponent `⌊s^2/(2p)⌋` loses at most one from `s^2/(2p)`, after multiplication by
the nonnegative factor `log p`. -/
theorem floor_half_s_sq_div_log_term_lower {s p : ℕ}
    (hp : 0 < p) (hp1 : 1 ≤ p) :
    ((s : ℝ) ^ 2 / (p : ℝ) - 2) * Real.log (p : ℝ) ≤
      ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
        Real.log (p : ℝ) := by
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp
  have hlog_nonneg : 0 ≤ Real.log (p : ℝ) := Real.log_nonneg (by exact_mod_cast hp1)
  have hfloor_lt := Nat.sub_one_lt_floor (((s : ℝ) ^ 2) / (2 * (p : ℝ)))
  have hfloor_le :
      ((s : ℝ) ^ 2 / (p : ℝ) - 2) ≤
        ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) := by
    norm_num at hfloor_lt ⊢
    field_simp [hpR.ne'] at hfloor_lt ⊢
    nlinarith
  exact mul_le_mul_of_nonneg_right hfloor_le hlog_nonneg

/-- Sum form of `floor_half_s_sq_div_log_term_lower`. This is the bridge from weighted
prime-reciprocal sums to the floor-exponent logarithmic product. -/
theorem sum_two_floor_half_s_sq_div_mul_log_lower {P : Finset ℕ} {s : ℕ}
    (hpos : ∀ p ∈ P, 0 < p) (hone : ∀ p ∈ P, 1 ≤ p) :
    (s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
        2 * (∑ p ∈ P, Real.log (p : ℝ)) ≤
      ∑ p ∈ P,
        ((2 * Nat.floor (((s : ℝ) ^ 2) / (2 * (p : ℝ))) : ℕ) : ℝ) *
          Real.log (p : ℝ) := by
  have hrewrite :
      (s : ℝ) ^ 2 * (∑ p ∈ P, Real.log (p : ℝ) / (p : ℝ)) -
          2 * (∑ p ∈ P, Real.log (p : ℝ)) =
        ∑ p ∈ P, (((s : ℝ) ^ 2 / (p : ℝ) - 2) * Real.log (p : ℝ)) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl ?_
    intro p hp
    have hpR : (p : ℝ) ≠ 0 := by exact_mod_cast ne_of_gt (hpos p hp)
    field_simp [hpR]
  rw [hrewrite]
  exact Finset.sum_le_sum fun p hp =>
    floor_half_s_sq_div_log_term_lower (hpos p hp) (hone p hp)

end LogProductBounds
end GaussianChain
