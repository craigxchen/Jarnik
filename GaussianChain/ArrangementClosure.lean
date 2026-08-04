import Mathlib

namespace GaussianChain
namespace ArrangementClosure

/--
The closure of a set with respect to a fixed family of vanishing predicates.
It consists of all points on every arrangement hyperplane that contains `S`.
-/
def closure
    {α ι : Type*}
    (zero : ι → α → Prop) (S : Set α) : Set α :=
  {x | ∀ i, (∀ y ∈ S, zero i y) → zero i x}

/-- Every set is contained in its arrangement closure. -/
theorem subset_closure
    {α ι : Type*}
    (zero : ι → α → Prop) (S : Set α) :
    S ⊆ closure zero S := by
  intro x hx i hi
  exact hi x hx

/-- A form vanishes on the closure exactly when it vanishes on the original set. -/
theorem vanishes_on_closure_iff
    {α ι : Type*}
    (zero : ι → α → Prop) (S : Set α) (i : ι) :
    (∀ x ∈ closure zero S, zero i x) ↔ (∀ x ∈ S, zero i x) := by
  constructor
  · intro h x hx
    exact h x (subset_closure zero S hx)
  · intro h x hx
    exact hx i h

/-- Arrangement closure is idempotent. -/
theorem closure_closure
    {α ι : Type*}
    (zero : ι → α → Prop) (S : Set α) :
    closure zero (closure zero S) = closure zero S := by
  apply Set.Subset.antisymm
  · intro x hx i hi
    apply hx i
    exact (vanishes_on_closure_iff zero S i).2 hi
  · exact subset_closure zero (closure zero S)

/--
Abstract canonical-maximizer closure theorem.

Assume `score` cannot decrease after arrangement closure, `T` has globally
maximal score, and `T` contains every other maximizer.  Then `T` is closed
under the fixed arrangement.
-/
theorem canonical_eq_closure
    {α ι : Type*}
    (zero : ι → α → Prop)
    (score : Set α → ℝ)
    (T : Set α)
    (hscore : ∀ S, score S ≤ score (closure zero S))
    (hmax : ∀ S, score S ≤ score T)
    (hgreatest : ∀ S, score S = score T → S ⊆ T) :
    T = closure zero T := by
  have hle : score T ≤ score (closure zero T) := hscore T
  have hge : score (closure zero T) ≤ score T := hmax (closure zero T)
  have heq : score (closure zero T) = score T := le_antisymm hge hle
  apply Set.Subset.antisymm
  · exact subset_closure zero T
  · exact hgreatest (closure zero T) heq

/--
A convenient variant where maximality is stated by equality to a prescribed
maximum value `μ`.
-/
theorem canonical_eq_closure_of_value
    {α ι : Type*}
    (zero : ι → α → Prop)
    (score : Set α → ℝ)
    (T : Set α) (μ : ℝ)
    (hscore : ∀ S, score S ≤ score (closure zero S))
    (hupper : ∀ S, score S ≤ μ)
    (hT : score T = μ)
    (hgreatest : ∀ S, score S = μ → S ⊆ T) :
    T = closure zero T := by
  apply canonical_eq_closure zero score T hscore
  · intro S
    simpa [hT] using hupper S
  · intro S hS
    apply hgreatest S
    simpa [hT] using hS

end ArrangementClosure
end GaussianChain
