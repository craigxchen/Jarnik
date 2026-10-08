import Mathlib

namespace GaussianChain
namespace PrimePowerMajority

/-- The half exponent used to split a prime-power orientation block. -/
def halfExponent (e : ℕ) : ℕ := (e + 1) / 2

/-- Every exponent in `[0,e]` lies in one of the two half-power orientation classes. -/
theorem halfExponent_le_or_halfExponent_le_complement
    {e a : ℕ} (ha : a ≤ e) :
    halfExponent e ≤ a ∨ halfExponent e ≤ e - a := by
  unfold halfExponent
  omega

/--
For a finite family of exponents in `[0,e]`, at least half the family lies in
one of the two half-power classes.
-/
theorem many_high_or_many_low
    {n e : ℕ} (a : Fin n → ℕ) (ha : ∀ i, a i ≤ e) :
    n ≤ 2 * (Finset.univ.filter fun i => halfExponent e ≤ a i).card ∨
      n ≤ 2 * (Finset.univ.filter fun i => halfExponent e ≤ e - a i).card := by
  classical
  let A : Finset (Fin n) :=
    Finset.univ.filter fun i => halfExponent e ≤ a i
  let B : Finset (Fin n) :=
    Finset.univ.filter fun i => halfExponent e ≤ e - a i
  have hunion : Finset.univ ⊆ A ∪ B := by
    intro i hi
    rcases halfExponent_le_or_halfExponent_le_complement (ha i) with h | h
    · exact Finset.mem_union.mpr (Or.inl (by simp [A, h]))
    · exact Finset.mem_union.mpr (Or.inr (by simp [B, h]))
  have hcard : n ≤ A.card + B.card := by
    have h1 : n ≤ (A ∪ B).card := by
      simpa using Finset.card_le_card hunion
    exact h1.trans (Finset.card_union_le A B)
  have hhalf : n ≤ 2 * A.card ∨ n ≤ 2 * B.card := by omega
  simpa [A, B] using hhalf

/-- The chosen half-power removes at least half of the full prime-power exponent. -/
theorem two_mul_halfExponent_ge (e : ℕ) : e ≤ 2 * halfExponent e := by
  unfold halfExponent
  omega

/-- For positive exponent, the half exponent is positive. -/
theorem halfExponent_pos {e : ℕ} (he : 0 < e) : 0 < halfExponent e := by
  unfold halfExponent
  omega

end PrimePowerMajority
end GaussianChain
