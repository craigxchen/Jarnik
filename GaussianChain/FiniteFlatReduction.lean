import Mathlib

namespace GaussianChain
namespace FiniteFlatReduction

/--
If an exceptional object is determined by a subset of a fixed finite index set,
then only finitely many exceptional objects can occur.  This is the abstract
combinatorial core of the proposed canonical-flat reduction: a flat is encoded
by the adapted linear forms that vanish on it.
-/
theorem exceptional_range_card_le_two_pow
    {ι α β : Type*}
    [Fintype ι] [DecidableEq ι] [DecidableEq β]
    (A : Finset α)
    (signature : α → Finset ι)
    (canonical : Finset ι → β) :
    (A.image (fun x => canonical (signature x))).card ≤ 2 ^ Fintype.card ι := by
  let P : Finset (Finset ι) := Finset.univ.powerset
  have hsub : A.image (fun x => canonical (signature x)) ⊆ P.image canonical := by
    intro y hy
    obtain ⟨x, hxA, rfl⟩ := Finset.mem_image.mp hy
    apply Finset.mem_image.mpr
    refine ⟨signature x, ?_, rfl⟩
    exact Finset.mem_powerset.mpr (Finset.subset_univ _)
  calc
    (A.image (fun x => canonical (signature x))).card
        ≤ (P.image canonical).card := Finset.card_le_card hsub
    _ ≤ P.card := Finset.card_image_le
    _ = 2 ^ Fintype.card ι := by
      simp [P]

/--
A finite-family version with an arbitrary explicit signature family `Σ`.
This is useful when geometric compatibility excludes most subsets of forms.
-/
theorem exceptional_range_card_le_signature_family
    {ι α β : Type*}
    [DecidableEq ι] [DecidableEq β]
    (A : Finset α)
    (Σ : Finset (Finset ι))
    (signature : α → Finset ι)
    (canonical : Finset ι → β)
    (hsignature : ∀ x ∈ A, signature x ∈ Σ) :
    (A.image (fun x => canonical (signature x))).card ≤ Σ.card := by
  have hsub : A.image (fun x => canonical (signature x)) ⊆ Σ.image canonical := by
    intro y hy
    obtain ⟨x, hxA, rfl⟩ := Finset.mem_image.mp hy
    exact Finset.mem_image.mpr ⟨signature x, hsignature x hxA, rfl⟩
  exact (Finset.card_le_card hsub).trans Finset.card_image_le

end FiniteFlatReduction
end GaussianChain
