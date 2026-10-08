import GaussianChain.PrimeLayerSurplus

namespace GaussianChain
namespace PrimeLayerCentered

open scoped BigOperators

def centeredCut (s k : ℤ) : ℤ := 2 * k - (2 * s + 1)

theorem abs_centeredCut_le_one_add_defect
    (s k : ℤ)
    (hs : 0 ≤ s) (hk0 : 0 ≤ k) (hktop : k ≤ 2 * s + 1) :
    |centeredCut s k| ≤
      1 + (s * (s + 1) - k * (2 * s + 1 - k)) := by
  rw [PrimeLayerSurplus.odd_cut_defect_identity]
  rcases (by omega : k ≤ s ∨ s + 1 ≤ k) with hlow | hhigh
  · have hc : centeredCut s k ≤ 0 := by
      unfold centeredCut
      linarith
    rw [abs_of_nonpos hc]
    by_cases hks : k = s
    · subst k
      norm_num [centeredCut]
    · have hmul : 0 ≤ (s - k) * (s - k - 1) :=
        mul_nonneg (by omega) (by omega)
      unfold centeredCut
      nlinarith
  · have hc : 0 ≤ centeredCut s k := by
      unfold centeredCut
      linarith
    rw [abs_of_nonneg hc]
    by_cases hks : k = s + 1
    · subst k
      norm_num [centeredCut]
    · have hmul : 0 ≤ (k - s - 1) * (k - s - 2) :=
        mul_nonneg (by omega) (by omega)
      unfold centeredCut
      nlinarith

theorem defect_eq_zero_iff_balanced
    (s k : ℤ) :
    s * (s + 1) - k * (2 * s + 1 - k) = 0 ↔
      k = s ∨ k = s + 1 := by
  rw [PrimeLayerSurplus.odd_cut_defect_identity]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h₁ | h₂
    · exact Or.inl (sub_eq_zero.mp h₁)
    · right
      linarith
  · rintro (rfl | rfl) <;> ring

private theorem abs_finset_sum_le_sum_abs
    {ι : Type*} (t : Finset ι) (f : ι → ℤ) :
    |(∑ i in t, f i)| ≤ ∑ i in t, |f i| := by
  classical
  induction t using Finset.induction_on with
  | empty => simp
  | @insert a t ha ih =>
      rw [Finset.sum_insert ha, Finset.sum_insert ha]
      exact (abs_add (f a) (∑ i in t, f i)).trans
        (add_le_add_left ih |f a|)

theorem abs_sum_centeredCut_le_card_add_sum_defect
    {ι : Type*} [Fintype ι]
    (s : ℤ) (k : ι → ℤ)
    (hs : 0 ≤ s)
    (hk0 : ∀ i, 0 ≤ k i)
    (hktop : ∀ i, k i ≤ 2 * s + 1) :
    |∑ i, centeredCut s (k i)| ≤
      (Fintype.card ι : ℤ) +
        ∑ i, (s * (s + 1) - k i * (2 * s + 1 - k i)) := by
  calc
    |∑ i, centeredCut s (k i)|
        ≤ ∑ i, |centeredCut s (k i)| := by
          simpa using
            (abs_finset_sum_le_sum_abs (Finset.univ) (fun i => centeredCut s (k i)))
    _ ≤ ∑ i, (1 + (s * (s + 1) - k i * (2 * s + 1 - k i))) := by
      exact Finset.sum_le_sum fun i _ =>
        abs_centeredCut_le_one_add_defect s (k i) hs (hk0 i) (hktop i)
    _ = (Fintype.card ι : ℤ) +
        ∑ i, (s * (s + 1) - k i * (2 * s + 1 - k i)) := by
      simp [Finset.sum_add_distrib]

end PrimeLayerCentered
end GaussianChain
