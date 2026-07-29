import Mathlib

namespace GaussianChain
namespace SlopeClosure

/--
Abstract Harder--Narasimhan closure lemma.

Suppose `rk` is modular, `deg` is supermodular, and every object has slope at
most `μ`.  If `x` and `y` attain the maximal slope `μ`, then both their meet
and their join attain the same slope.

For subspaces, `rk` is dimension and `deg` is the sum of weighted filtration
intersection dimensions.  This is the algebraic reason the maximal-slope
subspaces have a unique largest member.
-/
theorem maximalSlope_inf_sup
    {α : Type*} [Lattice α]
    (rk deg : α → ℝ) (μ : ℝ)
    (hrank : ∀ x y,
      rk (x ⊓ y) + rk (x ⊔ y) = rk x + rk y)
    (hdeg : ∀ x y,
      deg x + deg y ≤ deg (x ⊓ y) + deg (x ⊔ y))
    (hupper : ∀ z, deg z ≤ μ * rk z)
    {x y : α}
    (hx : deg x = μ * rk x)
    (hy : deg y = μ * rk y) :
    deg (x ⊓ y) = μ * rk (x ⊓ y) ∧
      deg (x ⊔ y) = μ * rk (x ⊔ y) := by
  have hsuper := hdeg x y
  have hinter := hupper (x ⊓ y)
  have hjoin := hupper (x ⊔ y)
  have hmod := hrank x y
  constructor <;> nlinarith

/-- The join part of `maximalSlope_inf_sup`. -/
theorem maximalSlope_sup
    {α : Type*} [Lattice α]
    (rk deg : α → ℝ) (μ : ℝ)
    (hrank : ∀ x y,
      rk (x ⊓ y) + rk (x ⊔ y) = rk x + rk y)
    (hdeg : ∀ x y,
      deg x + deg y ≤ deg (x ⊓ y) + deg (x ⊔ y))
    (hupper : ∀ z, deg z ≤ μ * rk z)
    {x y : α}
    (hx : deg x = μ * rk x)
    (hy : deg y = μ * rk y) :
    deg (x ⊔ y) = μ * rk (x ⊔ y) :=
  (maximalSlope_inf_sup rk deg μ hrank hdeg hupper hx hy).2

/-- The meet part of `maximalSlope_inf_sup`. -/
theorem maximalSlope_inf
    {α : Type*} [Lattice α]
    (rk deg : α → ℝ) (μ : ℝ)
    (hrank : ∀ x y,
      rk (x ⊓ y) + rk (x ⊔ y) = rk x + rk y)
    (hdeg : ∀ x y,
      deg x + deg y ≤ deg (x ⊓ y) + deg (x ⊔ y))
    (hupper : ∀ z, deg z ≤ μ * rk z)
    {x y : α}
    (hx : deg x = μ * rk x)
    (hy : deg y = μ * rk y) :
    deg (x ⊓ y) = μ * rk (x ⊓ y) :=
  (maximalSlope_inf_sup rk deg μ hrank hdeg hupper hx hy).1

end SlopeClosure
end GaussianChain
