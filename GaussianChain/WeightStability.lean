import Mathlib

namespace GaussianChain
namespace WeightStability

open scoped BigOperators

/-- The L1 distance between two weight vectors on a fixed finite type set. -/
def l1Dist
    {κ : Type*} [Fintype κ]
    (w w' : κ → ℝ) : ℝ :=
  ∑ t, |w t - w' t|

/-- A weighted score attached to a finite family of local types. -/
def weightedScore
    {κ α : Type*} [Fintype κ]
    (w : κ → ℝ) (score : κ → α → ℝ) (x : α) : ℝ :=
  ∑ t, w t * score t x

/--
Weighted scores are Lipschitz in the L1 distance of the weights.  The bound is
scaled by a common point-height bound on every type score.
-/
theorem abs_weightedScore_sub_le
    {κ α : Type*} [Fintype κ]
    (w w' : κ → ℝ) (score : κ → α → ℝ) (x : α) (H : ℝ)
    (hscore : ∀ t, |score t x| ≤ H) :
    |weightedScore w score x - weightedScore w' score x| ≤
      H * l1Dist w w' := by
  classical
  have hrewrite :
      weightedScore w score x - weightedScore w' score x =
        ∑ t, (w t - w' t) * score t x := by
    unfold weightedScore
    rw [← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro t ht
    ring
  rw [hrewrite]
  calc
    |∑ t, (w t - w' t) * score t x|
        ≤ ∑ t, |(w t - w' t) * score t x| :=
      Finset.abs_sum_le_sum_abs _ _
    _ = ∑ t, |w t - w' t| * |score t x| := by
      apply Finset.sum_congr rfl
      intro t ht
      rw [abs_mul]
    _ ≤ ∑ t, |w t - w' t| * H := by
      apply Finset.sum_le_sum
      intro t ht
      exact mul_le_mul_of_nonneg_left (hscore t) (abs_nonneg _)
    _ = H * l1Dist w w' := by
      unfold l1Dist
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro t ht
      ring

/--
A strict coefficient gap is stable under a sufficiently small L1 perturbation
of the aggregate conductor weights.
-/
theorem weightedScore_lt_of_close
    {κ α : Type*} [Fintype κ]
    (w w₀ : κ → ℝ) (score : κ → α → ℝ) (x : α)
    (H η₀ η : ℝ)
    (hH : 0 ≤ H)
    (hscore : ∀ t, |score t x| ≤ H)
    (hbase : weightedScore w₀ score x < η₀ * H)
    (hclose : l1Dist w w₀ ≤ η - η₀) :
    weightedScore w score x < η * H := by
  have habs := abs_weightedScore_sub_le w w₀ score x H hscore
  have hdiff :
      weightedScore w score x - weightedScore w₀ score x ≤
        H * l1Dist w w₀ :=
    (le_abs_self _).trans habs
  have hmul : H * l1Dist w w₀ ≤ H * (η - η₀) :=
    mul_le_mul_of_nonneg_left hclose hH
  linarith

/--
Finite-net uniformization.

For each realized weight system `ω`, choose a nearby center `center ω` from a
fixed finite net.  If the fixed-weight theorem holds at every center with
coefficient `η₀`, then every point having score at least `η` at an arbitrary
weight is exceptional for its chosen center.
-/
theorem exceptional_of_finite_weight_net
    {κ Ω α : Type*} [Fintype κ]
    (centers : Finset Ω)
    (center : Ω → Ω)
    (weights : Ω → κ → ℝ)
    (score : κ → α → ℝ)
    (height : α → ℝ)
    (exceptional : Ω → α → Prop)
    (η₀ η : ℝ)
    (hcenter : ∀ ω, center ω ∈ centers)
    (hheight : ∀ x, 0 ≤ height x)
    (hscore : ∀ t x, |score t x| ≤ height x)
    (hclose : ∀ ω, l1Dist (weights ω) (weights (center ω)) ≤ η - η₀)
    (hfixed : ∀ c ∈ centers, ∀ x, ¬ exceptional c x →
      weightedScore (weights c) score x < η₀ * height x)
    {ω : Ω} {x : α}
    (hlarge : η * height x ≤ weightedScore (weights ω) score x) :
    exceptional (center ω) x := by
  by_contra hnot
  have hbase := hfixed (center ω) (hcenter ω) x hnot
  have hstable := weightedScore_lt_of_close
    (weights ω) (weights (center ω)) score x (height x) η₀ η
    (hheight x) (fun t => hscore t x) hbase (hclose ω)
  linarith

end WeightStability
end GaussianChain
