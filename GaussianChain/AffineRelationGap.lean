import Mathlib

namespace GaussianChain
namespace AffineRelationGap

open scoped BigOperators

/-- Squared Euclidean norm of an integer planar vector. -/
def normSq (x y : ℤ) : ℤ := x ^ 2 + y ^ 2

/-- Squared distance between two integer planar vectors. -/
def distSq (x y x₀ y₀ : ℤ) : ℤ := (x - x₀) ^ 2 + (y - y₀) ^ 2

/--
Global affine-relation identity for points on one centered circle.

If all points `(x i,y i)` have squared norm `N`, and the integer coefficient
vector `a` gives a linear relation in both coordinates, then its coefficient
sum is measured exactly by the weighted sum of the squared chords from any
chosen root point:

`sum_i a_i * |z_i-z₀|² = 2 N * sum_i a_i`.

This uses the whole relation at once and is stronger than a pairwise chord
estimate.
-/
theorem weighted_distSq_sum
    {ι : Type*} [Fintype ι]
    (a x y : ι → ℤ) (x₀ y₀ N : ℤ)
    (hnorm : ∀ i, normSq (x i) (y i) = N)
    (hroot : normSq x₀ y₀ = N)
    (hx : ∑ i, a i * x i = 0)
    (hy : ∑ i, a i * y i = 0) :
    ∑ i, a i * distSq (x i) (y i) x₀ y₀ =
      2 * N * ∑ i, a i := by
  classical
  simp only [distSq, normSq] at hnorm hroot ⊢
  calc
    ∑ i, a i * ((x i - x₀) ^ 2 + (y i - y₀) ^ 2)
        = ∑ i, (2 * N * a i - 2 * x₀ * (a i * x i) -
            2 * y₀ * (a i * y i)) := by
              apply Finset.sum_congr rfl
              intro i hi
              have hi := hnorm i
              nlinarith [hroot]
    _ = 2 * N * ∑ i, a i - 2 * x₀ * (∑ i, a i * x i) -
          2 * y₀ * (∑ i, a i * y i) := by
            simp only [Finset.mul_sum]
            ring
    _ = 2 * N * ∑ i, a i := by rw [hx, hy]; ring

/--
The same identity with the root chosen from the indexed family.
-/
theorem weighted_distSq_sum_root
    {ι : Type*} [Fintype ι]
    (a x y : ι → ℤ) (r : ι) (N : ℤ)
    (hnorm : ∀ i, normSq (x i) (y i) = N)
    (hx : ∑ i, a i * x i = 0)
    (hy : ∑ i, a i * y i = 0) :
    ∑ i, a i * distSq (x i) (y i) (x r) (y r) =
      2 * N * ∑ i, a i := by
  exact weighted_distSq_sum a x y (x r) (y r) N hnorm (hnorm r) hx hy

end AffineRelationGap
end GaussianChain
