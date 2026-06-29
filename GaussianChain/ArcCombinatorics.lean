import Mathlib

namespace GaussianChain
namespace ArcCombinatorics

open Finset

/-- A sliding-window sum over points in a real interval is at most the window
width times the interval length.

The monotonicity hypothesis that is natural for ordered arc parameters is not
needed for this estimate; the interval bounds alone suffice after telescoping. -/
theorem windowSum_le
    {M k : ℕ} {a L : ℝ} (t : ℕ → ℝ)
    (hk : k ≤ M)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L) :
    (∑ j ∈ range (M - k), (t (j + k) - t j)) ≤ k * L := by
  let S : ℕ → ℝ := fun n => ∑ i ∈ range n, t i
  have hshift :
      (∑ j ∈ range (M - k), t (j + k)) = S M - S k := by
    have hkm : k + (M - k) = M := Nat.add_sub_of_le hk
    calc
      (∑ j ∈ range (M - k), t (j + k))
          = ∑ j ∈ range (M - k), t (k + j) := by
              refine sum_congr rfl ?_
              intro j _
              rw [Nat.add_comm]
      _ = S M - S k := by
              rw [← Finset.sum_range_add_sub_sum_range t k (M - k)]
              simp [S, hkm]
  have hlast_eq :
      S M - S (M - k) = ∑ j ∈ range k, t ((M - k) + j) := by
    have hmk : M - k + k = M := Nat.sub_add_cancel hk
    rw [← Finset.sum_range_add_sub_sum_range t (M - k) k]
    simp [S, hmk]
  have hlast_le : S M - S (M - k) ≤ k * (a + L) := by
    rw [hlast_eq]
    calc
      (∑ j ∈ range k, t ((M - k) + j))
          ≤ ∑ _j ∈ range k, (a + L) := by
              refine sum_le_sum ?_
              intro j hj
              exact (hmem ((M - k) + j) (by
                have hjlt : j < k := mem_range.mp hj
                omega)).2
      _ = k * (a + L) := by simp [nsmul_eq_mul]; ring
  have hfirst_le : k * a ≤ S k := by
    calc
      k * a = ∑ _i ∈ range k, a := by simp [nsmul_eq_mul]
      _ ≤ ∑ i ∈ range k, t i := by
              refine sum_le_sum ?_
              intro i hi
              exact (hmem i (lt_of_lt_of_le (mem_range.mp hi) hk)).1
      _ = S k := rfl
  rw [sum_sub_distrib, hshift]
  nlinarith

/-- If every consecutive `k + 1` block has span at least `B`, then the number
of points is bounded by the sliding-window estimate. -/
theorem card_le_of_window_span
    {M k : ℕ} {a L B : ℝ} (t : ℕ → ℝ)
    (hk : k ≤ M)
    (hB : 0 < B)
    (hmem : ∀ i, i < M → a ≤ t i ∧ t i ≤ a + L)
    (hspan : ∀ j, j < M - k → B ≤ t (j + k) - t j) :
    (M : ℝ) ≤ k + k * L / B := by
  have hsum_lower :
      (M - k : ℕ) * B ≤
        ∑ j ∈ range (M - k), (t (j + k) - t j) := by
    calc
      (M - k : ℕ) * B = ∑ _j ∈ range (M - k), B := by simp [nsmul_eq_mul]
      _ ≤ ∑ j ∈ range (M - k), (t (j + k) - t j) := by
              refine sum_le_sum ?_
              intro j hj
              exact hspan j (mem_range.mp hj)
  have hsum_upper :
      (∑ j ∈ range (M - k), (t (j + k) - t j)) ≤ k * L :=
    windowSum_le t hk hmem
  have hsub_le : ((M - k : ℕ) : ℝ) ≤ k * L / B := by
    rw [le_div_iff₀ hB]
    exact hsum_lower.trans hsum_upper
  have hM : (M : ℝ) = ((M - k : ℕ) : ℝ) + k := by
    exact_mod_cast (Nat.sub_add_cancel hk).symm
  calc
    (M : ℝ) = ((M - k : ℕ) : ℝ) + k := hM
    _ ≤ k * L / B + k := by gcongr
    _ = k + k * L / B := by ring

end ArcCombinatorics
end GaussianChain
