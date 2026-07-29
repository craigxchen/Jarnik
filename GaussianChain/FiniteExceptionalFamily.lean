import Mathlib

namespace GaussianChain
namespace FiniteExceptionalFamily

/--
A finite-family grid lemma.

For each exceptional relation `R k`, at most `D` first coordinates are bad and,
outside those, each fiber has at most `D` points.  If every off-diagonal pair
from `A` lies in at least one relation from the fixed family `K`, then
`A.card ≤ K.card * D + 1`.

This is the exact combinatorial endgame needed when pair-dependent twisted
weights produce exceptional hypersurfaces from one fixed finite candidate
family, rather than one single hypersurface.
-/
theorem card_le_of_finite_exceptional_family
    {α κ : Type*} [DecidableEq α]
    (A : Finset α)
    (K : Finset κ)
    (R : κ → α → α → Prop)
    [∀ k, DecidableRel (R k)]
    (bad : κ → Finset α)
    (D : ℕ)
    (hbad : ∀ k ∈ K, (bad k).card ≤ D)
    (hfiber : ∀ k ∈ K, ∀ x ∈ A, x ∉ bad k →
      (A.filter (R k x)).card ≤ D)
    (hoff : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → ∃ k ∈ K, R k x y) :
    A.card ≤ K.card * D + 1 := by
  classical
  let B : Finset α := K.biUnion bad
  have hBcard : B.card ≤ K.card * D := by
    calc
      B.card ≤ ∑ k ∈ K, (bad k).card := Finset.card_biUnion_le
      _ ≤ ∑ _k ∈ K, D := by
        gcongr with k hk
        exact hbad k hk
      _ = K.card * D := by simp [Nat.mul_comm]
  by_cases hex : ∃ x ∈ A, x ∉ B
  · obtain ⟨x, hxA, hxB⟩ := hex
    have hxnotbad : ∀ k ∈ K, x ∉ bad k := by
      intro k hk hx
      exact hxB (Finset.mem_biUnion.mpr ⟨k, hk, hx⟩)
    have hsub : A.erase x ⊆ K.biUnion (fun k => A.filter (R k x)) := by
      intro y hy
      have hy' := Finset.mem_erase.mp hy
      obtain ⟨k, hk, hR⟩ := hoff x hxA y hy'.2 (Ne.symm hy'.1)
      exact Finset.mem_biUnion.mpr ⟨k, hk, Finset.mem_filter.mpr ⟨hy'.2, hR⟩⟩
    have hunion : (K.biUnion (fun k => A.filter (R k x))).card ≤ K.card * D := by
      calc
        (K.biUnion (fun k => A.filter (R k x))).card
            ≤ ∑ k ∈ K, (A.filter (R k x)).card := Finset.card_biUnion_le
        _ ≤ ∑ _k ∈ K, D := by
          gcongr with k hk
          exact hfiber k hk x hxA (hxnotbad k hk)
        _ = K.card * D := by simp [Nat.mul_comm]
    have herase := Finset.card_le_card hsub
    rw [Finset.card_erase_of_mem hxA] at herase
    omega
  · have hsub : A ⊆ B := by
      intro x hxA
      by_contra hxB
      exact hex ⟨x, hxA, hxB⟩
    exact (Finset.card_le_card hsub).trans (hBcard.trans (Nat.le_add_right _ _))

end FiniteExceptionalFamily
end GaussianChain
