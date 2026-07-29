import Mathlib

namespace GaussianChain
namespace ConductorWeights

open scoped BigOperators

/-- Aggregate local conductor weights according to a fixed finite type map. -/
def aggregateWeight
    {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq κ]
    (τ : ι → κ) (w : ι → ℝ) (t : κ) : ℝ :=
  ∑ i, if τ i = t then w i else 0

/-- Aggregation preserves the total conductor mass. -/
theorem sum_aggregateWeight
    {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq κ]
    (τ : ι → κ) (w : ι → ℝ) :
    ∑ t, aggregateWeight τ w t = ∑ i, w i := by
  classical
  unfold aggregateWeight
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i hi
  simp

/-- Nonnegative local weights give nonnegative aggregate weights. -/
theorem aggregateWeight_nonneg
    {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq κ]
    (τ : ι → κ) (w : ι → ℝ)
    (hw : ∀ i, 0 ≤ w i) (t : κ) :
    0 ≤ aggregateWeight τ w t := by
  classical
  unfold aggregateWeight
  positivity

/--
Every linear slope score factors through the aggregate conductor weights.
Thus the optimization depends on the fixed finite set of filtration types,
not on the number of active places.
-/
theorem weightedScore_eq_aggregateScore
    {ι κ α : Type*} [Fintype ι] [Fintype κ] [DecidableEq κ]
    (τ : ι → κ) (w : ι → ℝ) (score : κ → α → ℝ) (a : α) :
    ∑ i, w i * score (τ i) a =
      ∑ t, aggregateWeight τ w t * score t a := by
  classical
  unfold aggregateWeight
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i hi
  simp [ite_mul]

/-- Two local systems with the same aggregate weights give identical scores. -/
theorem weightedScore_eq_of_aggregateWeight_eq
    {ι ι' κ α : Type*}
    [Fintype ι] [Fintype ι'] [Fintype κ] [DecidableEq κ]
    (τ : ι → κ) (τ' : ι' → κ)
    (w : ι → ℝ) (w' : ι' → ℝ)
    (score : κ → α → ℝ)
    (hagg : ∀ t, aggregateWeight τ w t = aggregateWeight τ' w' t)
    (a : α) :
    ∑ i, w i * score (τ i) a =
      ∑ i, w' i * score (τ' i) a := by
  rw [weightedScore_eq_aggregateScore, weightedScore_eq_aggregateScore]
  apply Finset.sum_congr rfl
  intro t ht
  rw [hagg t]

end ConductorWeights
end GaussianChain
