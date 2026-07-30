import GaussianChain.PrimeLayerSurplus

namespace GaussianChain
namespace PrimeLayerCentered

open scoped BigOperators

/-- The centered signed size of a cut of an odd set of cardinality `2s+1`. -/
def centeredCut (s k : ℤ) : ℤ := 2 * k - (2 * s + 1)

/--
The absolute centered cut size is controlled by one plus the Ramana cut defect.
The two zero-defect cuts `k=s` and `k=s+1` have centered size exactly one;
every step away from them is paid for by the quadratic defect.
-/
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
      simp [centeredCut]
    · have hd1 : 1 ≤ s - k := by omega
      have hmul : 0 ≤ (s - k) * (s - k - 1) :=
        mul_nonneg (by linarith) (by linarith)
      unfold centeredCut
      nlinarith
  · have hc : 0 ≤ centeredCut s k := by
      unfold centeredCut
      linarith
    rw [abs_of_nonneg hc]
    by_cases hks : k = s + 1
    · subst k
      simp [centeredCut]
    · have hd1 : 1 ≤ k - (s + 1) := by omega
      have hmul : 0 ≤ (k - s - 1) * (k - s - 2) :=
        mul_nonneg (by linarith) (by linarith)
      unfold centeredCut
      nlinarith

/--
The cut defect vanishes exactly at the two balanced cut sizes.
-/
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

/--
Summed version of `abs_centeredCut_le_one_add_defect`.

For a chain of prime-power threshold cuts, the absolute centered exponent sum
is at most one unit per layer plus the total Ramana surplus.  This is the exact
bridge from determinant near-equality to the height of the product of all
selected norm-one ratios.
-/
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
        ≤ ∑ i, |centeredCut s (k i)| := abs_sum_le_sum_abs _
    _ ≤ ∑ i, (1 + (s * (s + 1) - k i * (2 * s + 1 - k i))) := by
      exact Finset.sum_le_sum fun i _ =>
        abs_centeredCut_le_one_add_defect s (k i) hs (hk0 i) (hktop i)
    _ = (Fintype.card ι : ℤ) +
        ∑ i, (s * (s + 1) - k i * (2 * s + 1 - k i)) := by
      simp [Finset.sum_add_distrib]

end PrimeLayerCentered
end GaussianChain
