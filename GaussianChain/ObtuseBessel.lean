import Mathlib

/-!
# An obtuse-family Bessel inequality

This is the Hilbert-space ingredient of `docs/uniform_profile_extraction.md`.
The Gaussian construction and the probabilistic profile extraction remain
prose proofs. No uniform endpoint cardinality theorem is asserted here.
-/

namespace GaussianChain.ObtuseBessel

open scoped BigOperators

variable {ι E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

theorem nonneg_combination_norm_sq_le
    (s : Finset ι) (f : ι → E) (a : ι → ℝ)
    (ha : ∀ i ∈ s, 0 ≤ a i)
    (hf : ∀ i ∈ s, ‖f i‖ = 1)
    (hob : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → inner ℝ (f i) (f j) ≤ 0) :
    ‖∑ i ∈ s, a i • f i‖ ^ 2 ≤ ∑ i ∈ s, (a i) ^ 2 := by
  classical
  rw [← real_inner_self_eq_norm_sq]
  simp only [sum_inner, inner_sum, real_inner_smul_left, real_inner_smul_right]
  apply Finset.sum_le_sum
  intro i hi
  calc
    ∑ j ∈ s, a i * (a j * inner ℝ (f j) (f i)) ≤
        ∑ j ∈ s, if j = i then (a i) ^ 2 else 0 := by
      apply Finset.sum_le_sum
      intro j hj
      by_cases hji : j = i
      · subst j
        simp [hf i hi, sq]
      · simp only [if_neg hji]
        exact mul_nonpos_of_nonneg_of_nonpos (ha i hi)
          (mul_nonpos_of_nonneg_of_nonpos (ha j hj) (hob j hj i hi hji))
    _ = (a i) ^ 2 := by simp [hi]

theorem positive_inner_sq_sum_le
    (s : Finset ι) (f : ι → E) (h : E)
    (hf : ∀ i ∈ s, ‖f i‖ = 1)
    (hob : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → inner ℝ (f i) (f j) ≤ 0) :
    (∑ i ∈ s, (max (inner ℝ (f i) h) 0) ^ 2) ≤ ‖h‖ ^ 2 := by
  classical
  let a : ι → ℝ := fun i => max (inner ℝ (f i) h) 0
  let v : E := ∑ i ∈ s, a i • f i
  have hv : ‖v‖ ^ 2 ≤ ∑ i ∈ s, (a i) ^ 2 :=
    nonneg_combination_norm_sq_le s f a (fun i _ => le_max_right _ _) hf hob
  have hav : inner ℝ v h = ∑ i ∈ s, (a i) ^ 2 := by
    dsimp [v]
    rw [sum_inner]
    apply Finset.sum_congr rfl
    intro i hi
    rw [real_inner_smul_left]
    dsimp [a]
    by_cases ht : 0 ≤ inner ℝ (f i) h
    · rw [max_eq_left ht]
      ring
    · rw [max_eq_right (le_of_not_ge ht)]
      ring
  have hp := sq_nonneg ‖h - v‖
  have hav' : inner ℝ h v = ∑ i ∈ s, (a i) ^ 2 := by
    simpa only [real_inner_comm] using hav
  rw [norm_sub_sq_real, hav'] at hp
  change (∑ i ∈ s, (a i) ^ 2) ≤ ‖h‖ ^ 2
  linarith

/-- Pairwise nonpositive inner products give Bessel's inequality with
constant two, with no dependence on the number of vectors. -/
theorem inner_sq_sum_le_twice
    (s : Finset ι) (f : ι → E) (h : E)
    (hf : ∀ i ∈ s, ‖f i‖ = 1)
    (hob : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → inner ℝ (f i) (f j) ≤ 0) :
    (∑ i ∈ s, (inner ℝ (f i) h) ^ 2) ≤ 2 * ‖h‖ ^ 2 := by
  have hp := positive_inner_sq_sum_le s f h hf hob
  have hn := positive_inner_sq_sum_le s f (-h) hf hob
  simp only [inner_neg_right, norm_neg] at hn
  have hid (x : ℝ) : x ^ 2 = (max x 0) ^ 2 + (max (-x) 0) ^ 2 := by
    by_cases hx : 0 ≤ x
    · rw [max_eq_left hx, max_eq_right (neg_nonpos.mpr hx)]
      ring
    · have hx' : x ≤ 0 := le_of_not_ge hx
      rw [max_eq_right hx', max_eq_left (neg_nonneg.mpr hx')]
      ring
  simp_rw [hid (inner ℝ (f _) h)]
  rw [Finset.sum_add_distrib]
  linarith

end GaussianChain.ObtuseBessel
