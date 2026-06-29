import Mathlib.LinearAlgebra.Vandermonde
import Mathlib.NumberTheory.Zsqrtd.GaussianInt

namespace GaussianChain
namespace RamanaDeterminant

open Matrix
open scoped ComplexConjugate

/-- The Laurent monomial row used in Ramana's determinant.

For columns `0, ..., s` this is `conj z` to descending powers `s, ..., 0`; for
columns `s, ..., 2s` it is `z` to ascending powers `0, ..., s`. -/
def ramanaEntry (s : ℕ) (z : GaussianInt) (j : Fin (2 * s + 1)) : GaussianInt :=
  if (j : ℕ) ≤ s then star z ^ (s - (j : ℕ)) else z ^ ((j : ℕ) - s)

/-- The square Ramana matrix for `2s+1` Gaussian integers. -/
def ramanaMatrix (s : ℕ) (z : Fin (2 * s + 1) → GaussianInt) :
    Matrix (Fin (2 * s + 1)) (Fin (2 * s + 1)) GaussianInt :=
  fun i j => ramanaEntry s (z i) j

/-- Column scaling that turns the row-scaled Ramana matrix into an ordinary Vandermonde matrix
on a common norm circle. -/
def ramanaColumnScale (s : ℕ) (N : ℤ) (j : Fin (2 * s + 1)) : GaussianInt :=
  if (j : ℕ) ≤ s then (N : GaussianInt) ^ (s - (j : ℕ)) else 1

private theorem sum_range_sub_id (s : ℕ) :
    (∑ k ∈ Finset.range (s + 1), (s - k)) = s * (s + 1) / 2 := by
  have h := Finset.sum_range_reflect (fun k : ℕ => k) (s + 1)
  rw [Finset.sum_range_id] at h
  have hleft : (∑ j ∈ Finset.range (s + 1), (s + 1 - 1 - j)) =
      ∑ j ∈ Finset.range (s + 1), (s - j) := by
    simp
  have hright : (s + 1) * (s + 1 - 1) / 2 = s * (s + 1) / 2 := by
    rw [Nat.add_sub_cancel_right]
    rw [Nat.mul_comm]
  rw [hleft, hright] at h
  exact h

/-- The product of Ramana's column scalars is the expected triangular power of the common norm. -/
theorem prod_ramanaColumnScale (s : ℕ) (N : ℤ) :
    (∏ j : Fin (2 * s + 1), ramanaColumnScale s N j) =
      (N : GaussianInt) ^ (s * (s + 1) / 2) := by
  unfold ramanaColumnScale
  rw [Fin.prod_univ_eq_prod_range
    (fun k : ℕ => if k ≤ s then (N : GaussianInt) ^ (s - k) else 1) (2 * s + 1)]
  have hsubset : Finset.range (s + 1) ⊆ Finset.range (2 * s + 1) := by
    intro k hk
    rw [Finset.mem_range] at hk ⊢
    omega
  have hprodSubset := Finset.prod_subset (s₁ := Finset.range (s + 1))
    (s₂ := Finset.range (2 * s + 1))
    (f := fun k : ℕ => if k ≤ s then (N : GaussianInt) ^ (s - k) else 1)
    hsubset (by
      intro k _hk hnot
      have hnotle : ¬ k ≤ s := by
        intro hks
        exact hnot (Finset.mem_range.mpr (Nat.lt_succ_of_le hks))
      simp [hnotle])
  rw [← hprodSubset]
  have hdropIf : (∏ k ∈ Finset.range (s + 1),
      (if k ≤ s then (N : GaussianInt) ^ (s - k) else 1)) =
      ∏ k ∈ Finset.range (s + 1), (N : GaussianInt) ^ (s - k) := by
    refine Finset.prod_congr rfl ?_
    intro k hk
    have hk_le : k ≤ s := Nat.lt_succ_iff.mp (Finset.mem_range.mp hk)
    simp [hk_le]
  rw [hdropIf, Finset.prod_pow_eq_pow_sum, sum_range_sub_id]

/-- On a circle `z * conj z = N`, multiplying a Ramana row by `z^s` gives the corresponding
Vandermonde row, up to the fixed column scaling `N^(s-j)` in the first half of the columns. -/
theorem row_mul_ramanaEntry_eq_columnScale_mul_pow
    {s : ℕ} {N : ℤ} {z : GaussianInt}
    (hcircle : z * star z = (N : GaussianInt)) (j : Fin (2 * s + 1)) :
    z ^ s * ramanaEntry s z j = ramanaColumnScale s N j * z ^ (j : ℕ) := by
  unfold ramanaEntry ramanaColumnScale
  by_cases hjs : (j : ℕ) ≤ s
  · simp only [hjs, ↓reduceIte]
    calc
      z ^ s * star z ^ (s - (j : ℕ))
          = (z ^ (s - (j : ℕ)) * z ^ (j : ℕ)) * star z ^ (s - (j : ℕ)) := by
              rw [← pow_add, Nat.sub_add_cancel hjs]
      _ = (z ^ (s - (j : ℕ)) * star z ^ (s - (j : ℕ))) * z ^ (j : ℕ) := by
              ac_rfl
      _ = (z * star z) ^ (s - (j : ℕ)) * z ^ (j : ℕ) := by
              rw [mul_pow]
      _ = (N : GaussianInt) ^ (s - (j : ℕ)) * z ^ (j : ℕ) := by
              rw [hcircle]
  · have hsj : s ≤ (j : ℕ) := Nat.le_of_not_ge hjs
    simp only [hjs, ↓reduceIte, one_mul]
    rw [← pow_add, Nat.add_sub_of_le hsj]

/-- Matrix form of `row_mul_ramanaEntry_eq_columnScale_mul_pow`. -/
theorem rowScale_ramana_eq_vandermonde_colScale
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    Matrix.diagonal (fun i => z i ^ s) * ramanaMatrix s z =
      Matrix.vandermonde z * Matrix.diagonal (ramanaColumnScale s N) := by
  apply Matrix.ext
  intro i j
  calc
    (Matrix.diagonal (fun i => z i ^ s) * ramanaMatrix s z) i j
        = z i ^ s * ramanaEntry s (z i) j := by
            simp [ramanaMatrix, Matrix.mul_apply, Matrix.diagonal]
    _ = ramanaColumnScale s N j * z i ^ (j : ℕ) :=
            row_mul_ramanaEntry_eq_columnScale_mul_pow (hcircle i) j
    _ = (Matrix.vandermonde z * Matrix.diagonal (ramanaColumnScale s N)) i j := by
            simp [Matrix.mul_apply, Matrix.diagonal, mul_comm]

/-- Determinant form of the row/column scaling identity. -/
theorem prod_pow_mul_det_ramana_eq_det_vandermonde_mul_prod_colScale
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (∏ i, z i ^ s) * (ramanaMatrix s z).det =
      (Matrix.vandermonde z).det * ∏ j, ramanaColumnScale s N j := by
  have h := congrArg Matrix.det (rowScale_ramana_eq_vandermonde_colScale (s := s) hcircle)
  simpa [Matrix.det_mul, Matrix.det_diagonal] using h

/-- Ramana determinant identity, with the common norm contribution retained as an explicit
finite product of column scalars. -/
theorem prod_pow_mul_det_ramana_eq_colScale_mul_vandermondeProduct
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (∏ i, z i ^ s) * (ramanaMatrix s z).det =
      (∏ j, ramanaColumnScale s N j) *
        ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i) := by
  rw [prod_pow_mul_det_ramana_eq_det_vandermonde_mul_prod_colScale hcircle,
    Matrix.det_vandermonde]
  ac_rfl

/-- Ramana determinant identity on a common norm circle. -/
theorem prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeProduct
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (∏ i, z i ^ s) * (ramanaMatrix s z).det =
      (N : GaussianInt) ^ (s * (s + 1) / 2) *
        ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i) := by
  rw [prod_pow_mul_det_ramana_eq_colScale_mul_vandermondeProduct hcircle,
    prod_ramanaColumnScale]

/-- Norm of the row-scaling product in Ramana's determinant identity. -/
theorem norm_prod_pow {s : ℕ} {z : Fin (2 * s + 1) → GaussianInt} :
    (∏ i, z i ^ s).norm = ∏ i, (z i).norm ^ s := by
  calc
    (∏ i, z i ^ s).norm = ∏ i, (z i ^ s).norm := by
      exact (map_prod (Zsqrtd.normMonoidHom (d := -1)) (fun i => z i ^ s) Finset.univ)
    _ = ∏ i, (z i).norm ^ s := by
      refine Finset.prod_congr rfl ?_
      intro i _
      exact (map_pow (Zsqrtd.normMonoidHom (d := -1)) (z i) s)

/-- Norm of a rational-integer power, viewed inside `ℤ[i]`. -/
theorem norm_intCast_pow {k : ℕ} {N : ℤ} :
    (((N : GaussianInt) ^ k).norm) = (N * N) ^ k := by
  change (Zsqrtd.normMonoidHom (d := -1)) ((N : GaussianInt) ^ k) = (N * N) ^ k
  rw [map_pow]
  change ((N : GaussianInt).norm) ^ k = (N * N) ^ k
  rw [Zsqrtd.norm_intCast]

/-- Norm of the Vandermonde pair product is the product of the pair norms. -/
theorem norm_vandermondeProduct {s : ℕ} {z : Fin (2 * s + 1) → GaussianInt} :
    (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i)).norm =
      ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm := by
  calc
    (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i)).norm
        = ∏ i : Fin (2 * s + 1), (∏ j ∈ Finset.Ioi i, (z j - z i)).norm := by
            exact (map_prod (Zsqrtd.normMonoidHom (d := -1))
              (fun i : Fin (2 * s + 1) => ∏ j ∈ Finset.Ioi i, (z j - z i)) Finset.univ)
    _ = ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm := by
            refine Finset.prod_congr rfl ?_
            intro i _
            exact (map_prod (Zsqrtd.normMonoidHom (d := -1))
              (fun j : Fin (2 * s + 1) => z j - z i) (Finset.Ioi i))

/-- Natural absolute value distributes across a finite integer product. -/
theorem int_natAbs_prod {ι : Type*} (s : Finset ι) (f : ι → ℤ) :
    (∏ i ∈ s, f i).natAbs = ∏ i ∈ s, (f i).natAbs := by
  classical
  refine Finset.induction_on s ?base ?step
  · simp
  · intro a s ha hs
    simp [ha, hs, Int.natAbs_mul]

/-- Natural absolute value of the normed Vandermonde product, expanded pairwise. -/
theorem natAbs_vandermondeNormProduct {s : ℕ} {z : Fin (2 * s + 1) → GaussianInt} :
    (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs =
      ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm.natAbs := by
  rw [int_natAbs_prod Finset.univ]
  refine Finset.prod_congr rfl ?_
  intro i _
  rw [int_natAbs_prod (Finset.Ioi i)]

/-- The number of strict ordered-by-index pairs in `Fin (2s+1)` is `s(2s+1)`. -/
theorem sum_card_Ioi_fin_two_mul_add_one (s : ℕ) :
    (∑ i : Fin (2 * s + 1), (Finset.Ioi i).card) = s * (2 * s + 1) := by
  simp_rw [Fin.card_Ioi]
  rw [Fin.sum_univ_eq_sum_range (fun k => 2 * s + 1 - 1 - k) (2 * s + 1)]
  calc
    (∑ k ∈ Finset.range (2 * s + 1), (2 * s + 1 - 1 - k))
        = (2 * s) * (2 * s + 1) / 2 := by
            simpa using sum_range_sub_id (2 * s)
    _ = s * (2 * s + 1) := by
            rw [show (2 * s) * (2 * s + 1) = 2 * (s * (2 * s + 1)) by ring]
            exact Nat.mul_div_right (s * (2 * s + 1)) (by decide)

/-- If every pair norm in the Vandermonde product is at most `B`, then the whole product is
bounded by `B` to the number of pairs. -/
theorem natAbs_vandermondeNormProduct_le_pow
    {s : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B) :
    (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs
      ≤ B ^ (s * (2 * s + 1)) := by
  rw [natAbs_vandermondeNormProduct]
  calc
    (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm.natAbs)
        ≤ ∏ i : Fin (2 * s + 1), B ^ (Finset.Ioi i).card := by
            refine Finset.prod_le_prod' ?_
            intro i _
            calc
              (∏ j ∈ Finset.Ioi i, (z j - z i).norm.natAbs)
                  ≤ ∏ _j ∈ Finset.Ioi i, B := by
                      exact Finset.prod_le_prod' (fun j hj => hB i j hj)
              _ = B ^ (Finset.Ioi i).card := by
                      rw [Finset.prod_const]
    _ = B ^ (s * (2 * s + 1)) := by
            rw [Finset.prod_pow_eq_pow_sum, sum_card_Ioi_fin_two_mul_add_one]

/-- Extract the integer norm from an equation `z * conj z = N`. -/
theorem norm_eq_int_of_mul_star_eq {z : GaussianInt} {N : ℤ}
    (h : z * star z = (N : GaussianInt)) :
    z.norm = N := by
  have hcast : ((z.norm : ℤ) : GaussianInt) = (N : GaussianInt) := by
    rw [Zsqrtd.norm_eq_mul_conj]
    exact h
  simpa using congrArg Zsqrtd.re hcast

/-- On a common norm circle, the norm of the row-scaling product is a pure power of the
common norm. -/
theorem prod_norm_pow_eq_common
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (∏ i : Fin (2 * s + 1), (z i).norm ^ s) = N ^ (s * (2 * s + 1)) := by
  calc
    (∏ i : Fin (2 * s + 1), (z i).norm ^ s)
        = ∏ _i : Fin (2 * s + 1), N ^ s := by
            refine Finset.prod_congr rfl ?_
            intro i _
            rw [norm_eq_int_of_mul_star_eq (hcircle i)]
    _ = (N ^ s) ^ Fintype.card (Fin (2 * s + 1)) := by simp
    _ = N ^ (s * (2 * s + 1)) := by
          rw [Fintype.card_fin, pow_mul]

/-- Norm form of Ramana's determinant identity. This is the algebraic input for the eventual
subcritical chord-product lower bound. -/
theorem norm_prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeNorm
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    N ^ (s * (2 * s + 1)) * (ramanaMatrix s z).det.norm =
      (N * N) ^ (s * (s + 1) / 2) *
        (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm) := by
  have h := congrArg (Zsqrtd.normMonoidHom (d := -1))
    (prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeProduct hcircle)
  change (((∏ i, z i ^ s) * (ramanaMatrix s z).det).norm) =
    (((N : GaussianInt) ^ (s * (s + 1) / 2) *
      ∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i)).norm) at h
  rw [Zsqrtd.norm_mul, norm_prod_pow, prod_norm_pow_eq_common hcircle,
    Zsqrtd.norm_mul, norm_intCast_pow, norm_vandermondeProduct] at h
  exact h

/-- Natural-absolute-value form of the normed Ramana determinant identity. -/
theorem natAbs_norm_prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeNorm
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt)) :
    (N ^ (s * (2 * s + 1))).natAbs * (ramanaMatrix s z).det.norm.natAbs =
      ((N * N) ^ (s * (s + 1) / 2)).natAbs *
        (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs := by
  have h := congrArg Int.natAbs
    (norm_prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeNorm hcircle)
  simpa [Int.natAbs_mul] using h

/-- If the points are distinct and the common circle norm is nonzero, then Ramana's determinant
is nonzero. -/
theorem det_ramana_ne_zero_of_injective
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hN : (N : GaussianInt) ≠ 0) (hz : Function.Injective z) :
    (ramanaMatrix s z).det ≠ 0 := by
  have hvand : (Matrix.vandermonde z).det ≠ 0 :=
    Matrix.det_vandermonde_ne_zero_iff.mpr hz
  have hcolEntry : ∀ j : Fin (2 * s + 1), ramanaColumnScale s N j ≠ 0 := by
    intro j
    unfold ramanaColumnScale
    split
    · exact pow_ne_zero _ hN
    · exact one_ne_zero
  have hcol : (Matrix.diagonal (ramanaColumnScale s N)).det ≠ 0 := by
    rw [Matrix.det_diagonal]
    exact Finset.prod_ne_zero_iff.mpr fun j _ => hcolEntry j
  have h := congrArg Matrix.det (rowScale_ramana_eq_vandermonde_colScale (s := s) hcircle)
  rw [Matrix.det_mul, Matrix.det_mul] at h
  intro hramana
  have hright : (Matrix.vandermonde z).det *
      (Matrix.diagonal (ramanaColumnScale s N)).det ≠ 0 :=
    mul_ne_zero hvand hcol
  rw [hramana, mul_zero] at h
  exact hright h.symm

/-- Nonvanishing of Ramana's determinant as a nonzero Gaussian norm. -/
theorem det_ramana_norm_ne_zero_of_injective
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hN : (N : GaussianInt) ≠ 0) (hz : Function.Injective z) :
    (ramanaMatrix s z).det.norm ≠ 0 := by
  have hdet := det_ramana_ne_zero_of_injective hcircle hN hz
  exact fun hnorm => hdet (GaussianInt.norm_eq_zero.mp hnorm)

/-- A nonzero Ramana determinant has Gaussian norm of natural absolute value at least one. -/
theorem one_le_det_ramana_norm_natAbs_of_injective
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hN : (N : GaussianInt) ≠ 0) (hz : Function.Injective z) :
    1 ≤ (ramanaMatrix s z).det.norm.natAbs := by
  exact Nat.succ_le_of_lt
    (Int.natAbs_pos.mpr (det_ramana_norm_ne_zero_of_injective hcircle hN hz))

/-- Algebraic lower-bound hook: after taking norms in Ramana's identity, a nonzero determinant
contributes at least one to the natural absolute value. -/
theorem common_norm_pow_natAbs_le_vandermonde_side_of_injective
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hN : (N : GaussianInt) ≠ 0) (hz : Function.Injective z) :
    (N ^ (s * (2 * s + 1))).natAbs ≤
      ((N * N) ^ (s * (s + 1) / 2)).natAbs *
        (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs := by
  have hdet := one_le_det_ramana_norm_natAbs_of_injective hcircle hN hz
  have hident := natAbs_norm_prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeNorm hcircle
  calc
    (N ^ (s * (2 * s + 1))).natAbs
        = (N ^ (s * (2 * s + 1))).natAbs * 1 := by rw [mul_one]
    _ ≤ (N ^ (s * (2 * s + 1))).natAbs * (ramanaMatrix s z).det.norm.natAbs := by
        exact Nat.mul_le_mul_left _ hdet
    _ = ((N * N) ^ (s * (s + 1) / 2)).natAbs *
        (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs := hident

/-- Subcritical product bridge: a uniform upper bound for all pair norms converts Ramana's
determinant lower hook into a pure finite inequality. -/
theorem common_norm_pow_natAbs_le_of_pair_norm_le
    {s : ℕ} {N : ℤ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = (N : GaussianInt))
    (hN : (N : GaussianInt) ≠ 0) (hz : Function.Injective z)
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B) :
    (N ^ (s * (2 * s + 1))).natAbs ≤
      ((N * N) ^ (s * (s + 1) / 2)).natAbs * B ^ (s * (2 * s + 1)) := by
  exact (common_norm_pow_natAbs_le_vandermonde_side_of_injective hcircle hN hz).trans
    (Nat.mul_le_mul_left _ (natAbs_vandermondeNormProduct_le_pow hB))

/-- Natural-number version of `common_norm_pow_natAbs_le_of_pair_norm_le`, for circles
`z * conj z = N` with `0 < N`. -/
theorem common_nat_norm_pow_le_of_pair_norm_le
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N) (hz : Function.Injective z)
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B) :
    N ^ (s * (2 * s + 1)) ≤
      (N * N) ^ (s * (s + 1) / 2) * B ^ (s * (2 * s + 1)) := by
  have hNgi : ((N : ℤ) : GaussianInt) ≠ 0 := by
    intro hzero
    have hNzero_int : (N : ℤ) = 0 := by
      simpa using congrArg Zsqrtd.re hzero
    have hNzero : N = 0 := by
      exact_mod_cast hNzero_int
    omega
  have h := common_norm_pow_natAbs_le_of_pair_norm_le
    (N := (N : ℤ)) hcircle hNgi hz hB
  simpa using h

/-- The triangular norm factor in Ramana's inequality, simplified over natural numbers. -/
theorem norm_sq_pow_triangular (N s : ℕ) :
    (N * N) ^ (s * (s + 1) / 2) = N ^ (s * (s + 1)) := by
  rw [mul_pow, ← pow_add]
  congr 1
  rw [← Nat.two_mul]
  exact Nat.two_mul_div_two_of_even (Nat.even_mul_succ_self s)

/-- Cancel the common circle-norm factor in the natural-number Ramana product bridge. -/
theorem common_nat_norm_pow_cancelled_le_of_pair_norm_le
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N) (hz : Function.Injective z)
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B) :
    N ^ (s * s) ≤ B ^ (s * (2 * s + 1)) := by
  have hle := common_nat_norm_pow_le_of_pair_norm_le hcircle hN hz hB
  have hleft :
      N ^ (s * (2 * s + 1)) = N ^ (s * (s + 1)) * N ^ (s * s) := by
    rw [← pow_add]
    congr 1
    ring
  rw [hleft, norm_sq_pow_triangular] at hle
  exact Nat.le_of_mul_le_mul_left hle (Nat.pow_pos hN)

/-- Determinant-norm upper bound after cancelling the common circle-norm factor.

This is the strengthened version of `common_nat_norm_pow_cancelled_le_of_pair_norm_le` used when
extra divisibility gives a nontrivial lower bound for the Ramana determinant. -/
theorem det_norm_common_nat_norm_pow_le_of_pair_norm_le
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B) :
    N ^ (s * s) * (ramanaMatrix s z).det.norm.natAbs ≤
      B ^ (s * (2 * s + 1)) := by
  have hident := natAbs_norm_prod_pow_mul_det_ramana_eq_normPow_mul_vandermondeNorm hcircle
  have hV := natAbs_vandermondeNormProduct_le_pow (s := s) (z := z) (B := B) hB
  have hidentNat :
      N ^ (s * (2 * s + 1)) * (ramanaMatrix s z).det.norm.natAbs =
        (N * N) ^ (s * (s + 1) / 2) *
          (∏ i : Fin (2 * s + 1), ∏ j ∈ Finset.Ioi i, (z j - z i).norm).natAbs := by
    simpa using hident
  have hle :
      N ^ (s * (2 * s + 1)) * (ramanaMatrix s z).det.norm.natAbs ≤
        (N * N) ^ (s * (s + 1) / 2) * B ^ (s * (2 * s + 1)) := by
    rw [hidentNat]
    exact Nat.mul_le_mul_left _ hV
  rw [norm_sq_pow_triangular] at hle
  have hleft :
      N ^ (s * (2 * s + 1)) * (ramanaMatrix s z).det.norm.natAbs =
        N ^ (s * (s + 1)) *
          (N ^ (s * s) * (ramanaMatrix s z).det.norm.natAbs) := by
    rw [← mul_assoc, ← pow_add]
    congr 1
    ring
  rw [hleft] at hle
  exact Nat.le_of_mul_le_mul_left hle (Nat.pow_pos hN)

/-- If the uniform pair-norm bound would make the cancelled Ramana inequality false, then the
points cannot all be distinct. -/
theorem not_injective_of_pair_norm_pow_lt
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N)
    (hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B)
    (hsmall : B ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ¬ Function.Injective z := by
  intro hz
  have hle := common_nat_norm_pow_cancelled_le_of_pair_norm_le hcircle hN hz hB
  exact (not_lt_of_ge hle) hsmall

/-- Contrapositive subcritical chord product bound: if the Ramana threshold is subcritical for
`B`, then some pair has squared chord norm strictly larger than `B`. -/
theorem exists_pair_norm_gt_of_pow_lt
    {s N : ℕ} {z : Fin (2 * s + 1) → GaussianInt} {B : ℕ}
    (hcircle : ∀ i, z i * star (z i) = ((N : ℤ) : GaussianInt))
    (hN : 0 < N) (hz : Function.Injective z)
    (hsmall : B ^ (s * (2 * s + 1)) < N ^ (s * s)) :
    ∃ i j, j ∈ Finset.Ioi i ∧ B < (z j - z i).norm.natAbs := by
  classical
  by_contra hnone
  have hB : ∀ i j, j ∈ Finset.Ioi i → (z j - z i).norm.natAbs ≤ B := by
    intro i j hij
    exact le_of_not_gt fun hgt => hnone ⟨i, j, hij, hgt⟩
  exact (not_injective_of_pair_norm_pow_lt hcircle hN hB hsmall) hz

end RamanaDeterminant
end GaussianChain
