import Mathlib

namespace GaussianChain
namespace RationalPhaseObstruction

open scoped BigOperators

/-- A nonzero integral transverse coordinate has square at least one. -/
theorem one_le_sq_of_int_ne_zero (d : ℤ) (hd : d ≠ 0) :
    (1 : ℝ) ≤ (d : ℝ) ^ 2 := by
  have hz : 1 ≤ d ^ 2 := by
    rcases lt_or_gt_of_ne hd with h | h
    · have : d ≤ -1 := by omega
      nlinarith [sq_nonneg (d + 1)]
    · have : 1 ≤ d := by omega
      nlinarith [sq_nonneg (d - 1)]
  exact_mod_cast hz

/--
The product step in the rational-phase obstruction. There are `m + 1`
complementary factors with product of squared moduli `R ^ m`. If an
integral transverse coordinate forces `R ≤ K * n i` for every factor,
then the radius itself is bounded by `K ^ (m + 1)`.
-/
theorem radius_le_of_product
    (m : ℕ) (R K : ℝ) (n : Fin (m + 1) → ℝ)
    (hR : 0 < R)
    (hprod : (∏ i, n i) = R ^ m)
    (hgap : ∀ i, R ≤ K * n i) :
    R ≤ K ^ (m + 1) := by
  have h : (∏ _i : Fin (m + 1), R) ≤ ∏ i, K * n i :=
    Finset.prod_le_prod (fun _ _ ↦ hR.le) (fun i _ ↦ hgap i)
  simp only [Finset.prod_const, Finset.card_univ, Fintype.card_fin,
    Finset.prod_mul_distrib] at h
  rw [hprod, pow_succ] at h
  have hp : 0 < R ^ m := pow_pos hR m
  nlinarith

/--
Integer coordinates close to rational lines give the radius bound when
the product of their squared moduli has the balanced complementary exponent.
The hypotheses expose the two genuine inputs: nonzero integral transverse
coordinates, and the product identity. No endpoint cardinality theorem is
assumed or asserted.
-/
theorem radius_le_of_integral_transverse_bound
    (m : ℕ) (R K : ℝ)
    (x y a b : Fin (m + 1) → ℤ)
    (hR : 0 < R)
    (htrans : ∀ i, a i * y i - b i * x i ≠ 0)
    (hprod : (∏ i, ((x i : ℝ) ^ 2 + (y i : ℝ) ^ 2)) = R ^ m)
    (hnear : ∀ i,
      ((a i * y i - b i * x i : ℤ) : ℝ) ^ 2 * R ≤
        K * ((x i : ℝ) ^ 2 + (y i : ℝ) ^ 2)) :
    R ≤ K ^ (m + 1) := by
  apply radius_le_of_product m R K
    (fun i ↦ (x i : ℝ) ^ 2 + (y i : ℝ) ^ 2) hR hprod
  intro i
  have hi := one_le_sq_of_int_ne_zero (a i * y i - b i * x i) (htrans i)
  exact (le_mul_of_one_le_left hR.le hi).trans (hnear i)

/-- The seven complementary factors in the `k = 4` rational-ray branches. -/
theorem seven_factor_radius_bound
    (R C : ℝ) (x y a b : Fin 7 → ℤ)
    (hR : 0 < R)
    (htrans : ∀ i, a i * y i - b i * x i ≠ 0)
    (hprod : (∏ i, ((x i : ℝ) ^ 2 + (y i : ℝ) ^ 2)) = R ^ 6)
    (hnear : ∀ i,
      ((a i * y i - b i * x i : ℤ) : ℝ) ^ 2 * R ≤
        (9 * C ^ 2 / 8) * ((x i : ℝ) ^ 2 + (y i : ℝ) ^ 2)) :
    R ≤ (9 * C ^ 2 / 8) ^ 7 :=
  radius_le_of_integral_transverse_bound 6 R (9 * C ^ 2 / 8)
    x y a b hR htrans hprod hnear

end RationalPhaseObstruction
end GaussianChain
