import Mathlib

namespace GaussianChain
namespace NearMaximalReduction

/-- The candidates whose score lies within `ε` of a prescribed maximal value `μ`. -/
noncomputable def nearMaximal
    {α : Type*} [DecidableEq α]
    (F : Finset α) (score : α → ℝ) (μ ε : ℝ) : Finset α :=
  F.filter (fun x ↦ μ - ε ≤ score x)

/-- Every candidate outside the near-maximal set has a uniform strict gap `ε`. -/
theorem score_lt_of_mem_not_mem_nearMaximal
    {α : Type*} [DecidableEq α]
    {F : Finset α} {score : α → ℝ} {μ ε : ℝ} {x : α}
    (hxF : x ∈ F)
    (hx : x ∉ nearMaximal F score μ ε) :
    score x < μ - ε := by
  have h := hx
  simp [nearMaximal, hxF] at h
  linarith

/-- A genuine maximizer belongs to every nonnegative near-maximal set. -/
theorem maximizer_mem_nearMaximal
    {α : Type*} [DecidableEq α]
    {F : Finset α} {score : α → ℝ} {μ ε : ℝ} {x : α}
    (hε : 0 ≤ ε)
    (hxF : x ∈ F)
    (hx : score x = μ) :
    x ∈ nearMaximal F score μ ε := by
  unfold nearMaximal
  refine Finset.mem_filter.mpr ⟨hxF, ?_⟩
  rw [hx]
  linarith

/--
As the score data vary, the near-maximal signature is still only a subset of
one fixed finite candidate family.  Hence there are at most `2^(card F)`
signatures, with no dependence on the number of places producing the score.
-/
theorem nearMaximal_signature_range_card_le
    {α Ω : Type*} [DecidableEq α]
    (F : Finset α) (W : Finset Ω)
    (score : Ω → α → ℝ) (μ ε : Ω → ℝ) :
    (W.image (fun ω ↦ nearMaximal F (score ω) (μ ω) (ε ω))).card ≤
      2 ^ F.card := by
  let P : Finset (Finset α) := F.powerset
  have hsub :
      W.image (fun ω ↦ nearMaximal F (score ω) (μ ω) (ε ω)) ⊆ P := by
    intro S hS
    obtain ⟨ω, hω, rfl⟩ := Finset.mem_image.mp hS
    exact Finset.mem_powerset.mpr (Finset.filter_subset _ _)
  calc
    (W.image (fun ω ↦ nearMaximal F (score ω) (μ ω) (ε ω))).card
        ≤ P.card := Finset.card_le_card hsub
    _ = 2 ^ F.card := by simp [P]

end NearMaximalReduction
end GaussianChain
