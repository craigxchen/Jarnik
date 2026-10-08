import Mathlib

/-!
# Analytic cutoff for multiplicatively dependent endpoint ratios

This file formalizes two ingredients of docs/uniformity_continuation.md:
Gaussian-integer separation, and the uniform cutoff of the primitive relation
exponent. The factorization of circle ratios, the Bezout construction, removal
of units, and counting of the surviving square relations remain in prose.

No endpoint cardinality theorem is asserted here.
-/

namespace GaussianChain
namespace DependentEndpoint

/-- A nonzero Gaussian integer has complex absolute value at least one. -/
theorem one_le_complex_norm_of_ne_zero
    {z : GaussianInt} (hz : z ≠ 0) : 1 ≤ ‖(z : ℂ)‖ := by
  have hint : (1 : ℤ) ≤ z.norm := GaussianInt.norm_pos.mpr hz
  have hreal : (1 : ℝ) ≤ (z.norm : ℝ) := by exact_mod_cast hint
  rw [GaussianInt.intCast_real_norm, Complex.normSq_eq_norm_sq] at hreal
  nlinarith [norm_nonneg (z : ℂ)]

/-- Equal Gaussian norms make the squared separation an even integer. -/
theorem equal_norm_sub_norm_identity {a b : GaussianInt}
    (h : a.norm = b.norm) :
    (a - b).norm = 2 * (a.norm - (a.re * b.re + a.im * b.im)) := by
  simp only [Zsqrtd.norm_def, Zsqrtd.re_sub, Zsqrtd.im_sub] at *
  nlinarith

/-- Distinct lattice points of equal norm have squared separation at least two. -/
theorem two_le_norm_sub_of_norm_eq {a b : GaussianInt}
    (h : a.norm = b.norm) (hne : a ≠ b) : 2 ≤ (a - b).norm := by
  have hpos := GaussianInt.norm_pos.mpr (sub_ne_zero.mpr hne)
  have heven := equal_norm_sub_norm_identity h
  omega

/-- The sharp universal equal-radius lattice separation used by the finite
rectangle, angular-lift, and bounded-prefactor arguments. -/
theorem sqrt_two_le_complex_dist_of_norm_eq {a b : GaussianInt}
    (h : a.norm = b.norm) (hne : a ≠ b) :
    Real.sqrt 2 ≤ ‖(a : ℂ) - (b : ℂ)‖ := by
  have hint := two_le_norm_sub_of_norm_eq h hne
  have hreal : (2 : ℝ) ≤ ((a - b).norm : ℝ) := by exact_mod_cast hint
  rw [GaussianInt.intCast_real_norm, Complex.normSq_eq_norm_sq,
    GaussianInt.toComplex_sub] at hreal
  have hs := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)
  nlinarith [norm_nonneg ((a : ℂ) - (b : ℂ)), Real.sqrt_nonneg (2 : ℝ)]

/-- Equal-height Gaussian prefactors give a sharp relative separation. -/
theorem sqrt_two_div_le_quotient_sub_one {a b : GaussianInt}
    (h : a.norm = b.norm) (hne : a ≠ b) (hb : b ≠ 0) :
    Real.sqrt 2 / ‖(b : ℂ)‖ ≤ ‖(a : ℂ) / (b : ℂ) - 1‖ := by
  have hb' : (b : ℂ) ≠ 0 := GaussianInt.toComplex_eq_zero.not.mpr hb
  calc
    Real.sqrt 2 / ‖(b : ℂ)‖ ≤
        ‖(a : ℂ) - (b : ℂ)‖ / ‖(b : ℂ)‖ :=
      div_le_div_of_nonneg_right
        (sqrt_two_le_complex_dist_of_norm_eq h hne) (norm_nonneg _)
    _ = ‖(a : ℂ) / (b : ℂ) - 1‖ := by
      rw [← norm_div, sub_div, div_self hb']

/-- Gaussian numerator separation, without any assumption on the denominator.
For a zero denominator both sides are zero by the field convention. -/
theorem gaussian_fraction_separation
    {a : GaussianInt} (ha : a ≠ 0) (b : GaussianInt) :
    1 / ‖(b : ℂ)‖ ≤ ‖(a : ℂ) / (b : ℂ)‖ := by
  rw [norm_div]
  exact div_le_div_of_nonneg_right
    (one_le_complex_norm_of_ne_zero ha) (norm_nonneg _)

/-- The separation inequality for the actual unit-conjugate quotient.
Only nonvanishing of the Gaussian numerator is assumed; primitivity is the
number-theoretic input used to supply it in the accompanying note. -/
theorem unit_conjugate_separation
    (η W : GaussianInt) (hW : W ≠ 0) (hnum : η * W - star W ≠ 0) :
    1 / ‖(W : ℂ)‖ ≤ ‖(η : ℂ) * (W : ℂ) / star (W : ℂ) - 1‖ := by
  have hstar : star (W : ℂ) ≠ 0 := by
    simpa using (GaussianInt.toComplex_eq_zero.not.mpr hW)
  have hconj : (starRingEnd ℂ) (W : ℂ) ≠ 0 := hstar
  have h := gaussian_fraction_separation hnum (star W)
  simpa only [GaussianInt.toComplex_sub, GaussianInt.toComplex_mul,
    GaussianInt.toComplex_star, Complex.norm_conj, sub_div, div_self hconj] using h

/-- A convenient exact identity used in the exponent cutoff. -/
theorem twelfth_power_sq (R : ℝ) (hR : 0 ≤ R) :
    (R ^ (1 / 12 : ℝ)) ^ 2 = R ^ (1 / 6 : ℝ) := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul hR]
  norm_num

/-- Uniform analytic cutoff.

The capacity hypothesis comes from |W|^M ≤ R and |W| ≥ sqrt 2.
The separation hypothesis comes from the nonzero Gaussian numerator and the
Bezout estimate for the short arc. The large-radius condition is explicit and
independent of the conductor support.

The threshold used here is slightly weaker than the prose threshold:
the proof uses log R ≤ 12 R^(1/12) rather than its sharpened version
with the factor 1 / exp 1. -/
theorem primitive_exponent_le_two
    {C R : ℝ} (M : ℕ)
    (hC : 0 ≤ C) (hR : 1 ≤ R)
    (hcapacity : (M : ℝ) * Real.log (Real.sqrt 2) ≤ Real.log R)
    (hseparation : R ^ ((1 / 2 : ℝ) - 1 / (M : ℝ)) ≤ 2 * C * M)
    (hlarge : 24 * C <
      Real.log (Real.sqrt 2) * R ^ (1 / 12 : ℝ)) :
    M ≤ 2 := by
  by_contra hM
  have hM3 : (3 : ℝ) ≤ M := by exact_mod_cast (show 3 ≤ M by omega)
  have hR0 : 0 ≤ R := le_trans zero_le_one hR
  have hRpos : 0 < R := lt_of_lt_of_le zero_lt_one hR
  have hlpos : 0 < Real.log (Real.sqrt 2) := by
    apply Real.log_pos
    have hs := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)
    nlinarith [Real.sqrt_nonneg (2 : ℝ)]
  have hinv : 1 / (M : ℝ) ≤ (1 / 3 : ℝ) :=
    one_div_le_one_div_of_le (by norm_num) hM3
  have hexp : (1 / 6 : ℝ) ≤ 1 / 2 - 1 / (M : ℝ) := by linarith
  have hsmall : R ^ (1 / 6 : ℝ) ≤ 2 * C * M :=
    (Real.rpow_le_rpow_of_exponent_le hR hexp).trans hseparation
  have hsmall' := mul_le_mul_of_nonneg_right hsmall hlpos.le
  have hcapacity' := mul_le_mul_of_nonneg_left hcapacity (by positivity : 0 ≤ 2 * C)
  have hlog := Real.log_le_rpow_div hR0 (by norm_num : 0 < (1 / 12 : ℝ))
  have hlog' := mul_le_mul_of_nonneg_left hlog (by positivity : 0 ≤ 2 * C)
  have hbound : R ^ (1 / 6 : ℝ) * Real.log (Real.sqrt 2) ≤
      24 * C * R ^ (1 / 12 : ℝ) := by
    norm_num at hlog'
    nlinarith
  have hspos : 0 < R ^ (1 / 12 : ℝ) := Real.rpow_pos_of_pos hRpos _
  have hstrict := mul_lt_mul_of_pos_right hlarge hspos
  rw [← twelfth_power_sq R hR0] at hbound
  nlinarith

/-- The explicit sufficient threshold in an ordinary integer power. -/
theorem primitive_exponent_le_two_of_radius_bound
    {C R : ℝ} (M : ℕ)
    (hC : 0 ≤ C) (hR : 1 ≤ R)
    (hcapacity : (M : ℝ) * Real.log (Real.sqrt 2) ≤ Real.log R)
    (hseparation : R ^ ((1 / 2 : ℝ) - 1 / (M : ℝ)) ≤ 2 * C * M)
    (hlarge : (24 * C / Real.log (Real.sqrt 2)) ^ 12 < R) :
    M ≤ 2 := by
  apply primitive_exponent_le_two M hC hR hcapacity hseparation
  have hlpos : 0 < Real.log (Real.sqrt 2) := by
    apply Real.log_pos
    have hs := Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)
    nlinarith [Real.sqrt_nonneg (2 : ℝ)]
  have hR0 : 0 ≤ R := le_trans zero_le_one hR
  have hroot : (R ^ (1 / 12 : ℝ)) ^ 12 = R := by
    simpa using Real.rpow_inv_natCast_pow hR0 (by norm_num : (12 : ℕ) ≠ 0)
  have hlt : 24 * C / Real.log (Real.sqrt 2) < R ^ (1 / 12 : ℝ) := by
    by_contra hnot
    have hpow := pow_le_pow_left₀ (Real.rpow_nonneg hR0 _) (le_of_not_gt hnot) 12
    rw [hroot] at hpow
    linarith
  simpa [mul_comm] using (div_lt_iff₀ hlpos).mp hlt

/-- A bound strictly below one cannot contain a nonzero Gaussian numerator.
This is the integrality step behind the C < 1 no-dependent-pairs corollary. -/
theorem not_small_gaussian_numerator
    {z : GaussianInt} {C : ℝ} (hz : z ≠ 0) (hC : C < 1) :
    ¬ ‖(z : ℂ)‖ ≤ C := by
  have h := one_le_complex_norm_of_ne_zero hz
  linarith

end DependentEndpoint
end GaussianChain
