import Mathlib

namespace GaussianChain
namespace UniformExceptionalReduction

/--
The numerical gap at the heart of the endpoint argument.

If every off-diagonal pair has generalized-GCD proximity at least
`W / 2 - logC`, every joint height is at most `W`, and the uniform
exceptional theorem gives coefficient `η < 1 / 2`, then every pair must be
exceptional once `logC < (1 / 2 - η) W`.
-/
theorem offDiagonal_exceptional_of_gap
    {α : Type*} [DecidableEq α]
    (A : Finset α)
    (exceptional : α → α → Prop)
    (gcdPlus height : α → α → ℝ)
    (η W logC : ℝ)
    (hη : 0 ≤ η)
    (hgap : logC < ((1 / 2 : ℝ) - η) * W)
    (hgcd : ∀ x ∈ A, ∀ y ∈ A, x ≠ y →
      W / 2 - logC ≤ gcdPlus x y)
    (hheight : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → height x y ≤ W)
    (hOutside : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → ¬ exceptional x y →
      gcdPlus x y < η * height x y) :
    ∀ x ∈ A, ∀ y ∈ A, x ≠ y → exceptional x y := by
  intro x hx y hy hxy
  by_contra hnot
  have hg := hgcd x hx y hy hxy
  have hh := hheight x hx y hy hxy
  have ho := hOutside x hx y hy hxy hnot
  have hηheight : η * height x y ≤ η * W :=
    mul_le_mul_of_nonneg_left hh hη
  have hupper : gcdPlus x y < η * W := lt_of_lt_of_le ho hηheight
  nlinarith

/--
Abstract grid-zero lemma.

For a polynomial of total degree at most `D`, `bad` models the at most `D`
first-coordinate values at which the specialization in the second variable is
identically zero.  At every other first coordinate, the fiber has at most `D`
zeros.  If every ordered off-diagonal pair from `A` is a zero, then
`A.card ≤ D + 1`.
-/
theorem card_le_of_offDiagonal_and_fiber_bound
    {α : Type*} [DecidableEq α]
    (A bad : Finset α)
    (relation : α → α → Prop)
    [DecidableRel relation]
    (D : ℕ)
    (hbad : bad.card ≤ D)
    (hOffDiagonal : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → relation x y)
    (hFiber : ∀ x ∈ A, x ∉ bad → (A.filter (relation x)).card ≤ D) :
    A.card ≤ D + 1 := by
  by_cases hex : ∃ x ∈ A, x ∉ bad
  · obtain ⟨x, hxA, hxnotbad⟩ := hex
    have hsub : A.erase x ⊆ A.filter (relation x) := by
      intro y hy
      have hydata := Finset.mem_erase.mp hy
      have hyne : y ≠ x := hydata.1
      have hyA : y ∈ A := hydata.2
      exact Finset.mem_filter.mpr
        ⟨hyA, hOffDiagonal x hxA y hyA (Ne.symm hyne)⟩
    have hcard := Finset.card_le_card hsub
    rw [Finset.card_erase_of_mem hxA] at hcard
    have hfiber := hFiber x hxA hxnotbad
    omega
  · have hsub : A ⊆ bad := by
      intro x hxA
      by_contra hxnotbad
      exact hex ⟨x, hxA, hxnotbad⟩
    have hcard := Finset.card_le_card hsub
    omega

/--
The complete formalized endgame, conditional on the deep exceptional relation
and its degree/fiber bounds.

The root point is not included in `A`; consequently the corresponding cluster
has cardinality `A.card + 1`, bounded here by `D + 2`.
-/
theorem cluster_card_le_of_uniform_exceptional_relation
    {α : Type*} [DecidableEq α]
    (A bad : Finset α)
    (exceptional : α → α → Prop)
    [DecidableRel exceptional]
    (gcdPlus height : α → α → ℝ)
    (η W logC : ℝ)
    (D : ℕ)
    (hη : 0 ≤ η)
    (hgap : logC < ((1 / 2 : ℝ) - η) * W)
    (hgcd : ∀ x ∈ A, ∀ y ∈ A, x ≠ y →
      W / 2 - logC ≤ gcdPlus x y)
    (hheight : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → height x y ≤ W)
    (hOutside : ∀ x ∈ A, ∀ y ∈ A, x ≠ y → ¬ exceptional x y →
      gcdPlus x y < η * height x y)
    (hbad : bad.card ≤ D)
    (hFiber : ∀ x ∈ A, x ∉ bad → (A.filter (exceptional x)).card ≤ D) :
    A.card + 1 ≤ D + 2 := by
  have hOffDiagonal :=
    offDiagonal_exceptional_of_gap A exceptional gcdPlus height η W logC
      hη hgap hgcd hheight hOutside
  have hcard :=
    card_le_of_offDiagonal_and_fiber_bound A bad exceptional D hbad
      hOffDiagonal hFiber
  omega

end UniformExceptionalReduction
end GaussianChain
