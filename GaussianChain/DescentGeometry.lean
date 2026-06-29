import GaussianChain.PrimeDescent
import GaussianChain.SubcriticalBound

namespace GaussianChain
namespace DescentGeometry

open SubcriticalBound
open scoped ComplexConjugate

/-- Multiplication by a Gaussian integer scales squared Euclidean distances by its Gaussian norm. -/
theorem gaussianSqDist_mul_left (π z w : GaussianInt) :
    gaussianSqDist (π * z) (π * w) = (π.norm : ℝ) * gaussianSqDist z w := by
  unfold gaussianSqDist
  have hsub :
      ((π * w : GaussianInt) : ℂ) - ((π * z : GaussianInt) : ℂ) =
        (π : ℂ) * (((w : GaussianInt) : ℂ) - ((z : GaussianInt) : ℂ)) := by
    rw [map_mul, map_mul]
    ring
  calc
    Complex.normSq (((π * w : GaussianInt) : ℂ) - ((π * z : GaussianInt) : ℂ))
        = Complex.normSq ((π : ℂ) * (((w : GaussianInt) : ℂ) - ((z : GaussianInt) : ℂ))) := by
          rw [hsub]
    _ = Complex.normSq (π : ℂ) *
          Complex.normSq (((w : GaussianInt) : ℂ) - ((z : GaussianInt) : ℂ)) :=
        Complex.normSq_mul _ _
    _ = (π.norm : ℝ) * gaussianSqDist z w := by
        rw [← GaussianInt.intCast_real_norm]
        rfl

/-- Multiplication by a split Gaussian factor above `p` scales squared distances by `p`. -/
theorem gaussianSqDist_mul_left_of_mul_star_eq_nat
    {π z w : GaussianInt} {p : ℕ}
    (hπ : π * star π = (p : GaussianInt)) :
    gaussianSqDist (π * z) (π * w) = (p : ℝ) * gaussianSqDist z w := by
  rw [gaussianSqDist_mul_left, PrimeDescent.norm_eq_int_of_mul_star_eq_nat hπ]
  simp

/-- If two points have been divided by a split Gaussian factor above `p`, their squared distance
is the original squared distance divided by `p`. -/
theorem gaussianSqDist_quotient_eq_div_of_mul_left
    {π z w z' w' : GaussianInt} {p : ℕ} [Fact p.Prime]
    (hπ : π * star π = (p : GaussianInt))
    (hz : z = π * z') (hw : w = π * w') :
    gaussianSqDist z' w' = gaussianSqDist z w / (p : ℝ) := by
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).pos
  rw [hz, hw, gaussianSqDist_mul_left_of_mul_star_eq_nat hπ]
  field_simp [hp_pos.ne']

/-- Inequality form of split-factor descent for squared distances. -/
theorem gaussianSqDist_quotient_le_div_of_mul_left
    {π z w z' w' : GaussianInt} {p : ℕ} [Fact p.Prime]
    (hπ : π * star π = (p : GaussianInt))
    (hz : z = π * z') (hw : w = π * w')
    {B : ℝ} (hB : gaussianSqDist z w ≤ B) :
    gaussianSqDist z' w' ≤ B / (p : ℝ) := by
  rw [gaussianSqDist_quotient_eq_div_of_mul_left hπ hz hw]
  exact div_le_div_of_nonneg_right hB (by exact_mod_cast (Fact.out : Nat.Prime p).pos.le)

/-- Family form of split-factor distance scaling. -/
theorem gaussianSqDist_quotient_family_eq_div_of_mul_left
    {ι : Type*} {π : GaussianInt} {p : ℕ} [Fact p.Prime]
    {z w : ι → GaussianInt}
    (hπ : π * star π = (p : GaussianInt))
    (hzfac : ∀ i, z i = π * w i) (i j : ι) :
    gaussianSqDist (w i) (w j) = gaussianSqDist (z i) (z j) / (p : ℝ) :=
  gaussianSqDist_quotient_eq_div_of_mul_left hπ (hzfac i) (hzfac j)

/-- Family inequality form of split-factor distance scaling. -/
theorem gaussianSqDist_quotient_family_le_div_of_mul_left
    {ι : Type*} {π : GaussianInt} {p : ℕ} [Fact p.Prime]
    {z w : ι → GaussianInt}
    (hπ : π * star π = (p : GaussianInt))
    (hzfac : ∀ i, z i = π * w i) {i j : ι} {B : ℝ}
    (hB : gaussianSqDist (z i) (z j) ≤ B) :
    gaussianSqDist (w i) (w j) ≤ B / (p : ℝ) :=
  gaussianSqDist_quotient_le_div_of_mul_left hπ (hzfac i) (hzfac j) hB

/-- Multiplication by the rational Gaussian integer `p` scales squared distances by `p^2`. -/
theorem gaussianSqDist_natCast_mul_left (p : ℕ) (z w : GaussianInt) :
    gaussianSqDist ((p : GaussianInt) * z) ((p : GaussianInt) * w) =
      (p : ℝ) ^ 2 * gaussianSqDist z w := by
  rw [gaussianSqDist_mul_left]
  simp [pow_two]

/-- In inert descent, after dividing by the rational Gaussian integer `p`, squared distances are
divided by `p^2`. -/
theorem gaussianSqDist_quotient_eq_div_of_natCast_mul_left
    {p : ℕ} [Fact p.Prime] {z w z' w' : GaussianInt}
    (hz : z = (p : GaussianInt) * z') (hw : w = (p : GaussianInt) * w') :
    gaussianSqDist z' w' = gaussianSqDist z w / (p : ℝ) ^ 2 := by
  have hp_pos : 0 < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).pos
  rw [hz, hw, gaussianSqDist_natCast_mul_left]
  field_simp [hp_pos.ne']

/-- Inequality form of inert rational-prime descent for squared distances. -/
theorem gaussianSqDist_quotient_le_div_of_natCast_mul_left
    {p : ℕ} [Fact p.Prime] {z w z' w' : GaussianInt}
    (hz : z = (p : GaussianInt) * z') (hw : w = (p : GaussianInt) * w')
    {B : ℝ} (hB : gaussianSqDist z w ≤ B) :
    gaussianSqDist z' w' ≤ B / (p : ℝ) ^ 2 := by
  rw [gaussianSqDist_quotient_eq_div_of_natCast_mul_left hz hw]
  exact div_le_div_of_nonneg_right hB (sq_nonneg _)

/-- Family form of inert rational-prime distance scaling. -/
theorem gaussianSqDist_quotient_family_eq_div_of_natCast_mul_left
    {ι : Type*} {p : ℕ} [Fact p.Prime] {z w : ι → GaussianInt}
    (hzfac : ∀ i, z i = (p : GaussianInt) * w i) (i j : ι) :
    gaussianSqDist (w i) (w j) = gaussianSqDist (z i) (z j) / (p : ℝ) ^ 2 :=
  gaussianSqDist_quotient_eq_div_of_natCast_mul_left (hzfac i) (hzfac j)

/-- Family inequality form of inert rational-prime distance scaling. -/
theorem gaussianSqDist_quotient_family_le_div_of_natCast_mul_left
    {ι : Type*} {p : ℕ} [Fact p.Prime] {z w : ι → GaussianInt}
    (hzfac : ∀ i, z i = (p : GaussianInt) * w i) {i j : ι} {B : ℝ}
    (hB : gaussianSqDist (z i) (z j) ≤ B) :
    gaussianSqDist (w i) (w j) ≤ B / (p : ℝ) ^ 2 :=
  gaussianSqDist_quotient_le_div_of_natCast_mul_left (hzfac i) (hzfac j) hB

end DescentGeometry
end GaussianChain
