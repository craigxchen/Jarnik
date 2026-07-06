import GaussianChain.LogProductBounds

namespace GaussianChain
namespace SubcriticalLog

/-- Logarithmic sufficient condition for the ordinary subcritical Ramana threshold. -/
theorem natFloor_sq_pow_lt_of_log_lt
    {s N' : ℕ} {B : ℝ}
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hN' : 0 < N')
    (hlog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * Real.log (N' : ℝ)) :
    (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s) := by
  simpa using
    LogProductBounds.nat_pow_lt_prod_pow_mul_pow_of_log_lt
      (B := Nat.floor (B ^ 2)) (x := s * (2 * s + 1))
      (N := N') (y := s * s) (P := (∅ : Finset ℕ)) (E := fun _ => 0)
      hfloor hN' (by intro p hp; simp at hp) (by simpa using hlog)

/-- Natural-number factorization extracted from the integer split quotient equation. -/
theorem nat_eq_mul_of_int_eq_mul
    {N p N' : ℕ} (hfactor : (N : ℤ) = (p : ℤ) * (N' : ℤ)) :
    N = p * N' := by
  exact_mod_cast hfactor

/-- Natural-number factorization extracted from the integer inert quotient equation. -/
theorem nat_eq_sq_mul_of_int_eq_sq_mul
    {N p N' : ℕ} (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ)) :
    N = p ^ 2 * N' := by
  exact_mod_cast hfactor

/-- Logarithm of a positive split quotient norm. -/
theorem log_eq_log_mul_split_factor
    {N p N' : ℕ} (hp : 0 < p) (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) * (N' : ℤ)) :
    Real.log (N : ℝ) = Real.log (p : ℝ) + Real.log (N' : ℝ) := by
  have hfac : N = p * N' := nat_eq_mul_of_int_eq_mul hfactor
  rw [hfac, Nat.cast_mul]
  exact Real.log_mul (by exact_mod_cast hp.ne') (by exact_mod_cast hN'.ne')

/-- Logarithm of a positive inert quotient norm. -/
theorem log_eq_log_sq_mul_inert_factor
    {N p N' : ℕ} (hp : 0 < p) (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ)) :
    Real.log (N : ℝ) = 2 * Real.log (p : ℝ) + Real.log (N' : ℝ) := by
  have hfac : N = p ^ 2 * N' := nat_eq_sq_mul_of_int_eq_sq_mul hfactor
  rw [hfac, Nat.cast_mul, Nat.cast_pow]
  rw [Real.log_mul (by
      exact_mod_cast (pow_pos hp 2).ne') (by exact_mod_cast hN'.ne')]
  rw [Real.log_pow]
  ring

/-- Split-quotient subcritical threshold from a logarithmic inequality written against
`log N - log p`. -/
theorem split_natFloor_sq_pow_lt_of_log_sub
    {s N p N' : ℕ} {B : ℝ}
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hp : 0 < p) (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) * (N' : ℤ))
    (hlog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - Real.log (p : ℝ))) :
    (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s) := by
  refine natFloor_sq_pow_lt_of_log_lt hfloor hN' ?_
  have hlogfac := log_eq_log_mul_split_factor hp hN' hfactor
  have hrewrite : Real.log (N : ℝ) - Real.log (p : ℝ) = Real.log (N' : ℝ) := by
    rw [hlogfac]
    ring
  simpa [hrewrite] using hlog

/-- Inert-quotient subcritical threshold from a logarithmic inequality written against
`log N - 2 log p`. -/
theorem inert_natFloor_sq_pow_lt_of_log_sub
    {s N p N' : ℕ} {B : ℝ}
    (hfloor : 0 < Nat.floor (B ^ 2))
    (hp : 0 < p) (hN' : 0 < N')
    (hfactor : (N : ℤ) = (p : ℤ) ^ 2 * (N' : ℤ))
    (hlog :
      ((s * (2 * s + 1) : ℕ) : ℝ) * Real.log (Nat.floor (B ^ 2) : ℝ) <
        ((s * s : ℕ) : ℝ) * (Real.log (N : ℝ) - 2 * Real.log (p : ℝ))) :
    (Nat.floor (B ^ 2)) ^ (s * (2 * s + 1)) < N' ^ (s * s) := by
  refine natFloor_sq_pow_lt_of_log_lt hfloor hN' ?_
  have hlogfac := log_eq_log_sq_mul_inert_factor hp hN' hfactor
  have hrewrite : Real.log (N : ℝ) - 2 * Real.log (p : ℝ) = Real.log (N' : ℝ) := by
    rw [hlogfac]
    ring
  simpa [hrewrite] using hlog

end SubcriticalLog
end GaussianChain
