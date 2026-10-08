import Mathlib

namespace GaussianChain
namespace SexticTropicalWeight

open scoped BigOperators

/-- A three-coordinate weighted degree. -/
def weight3
    {ι : Type*}
    (e₀ e₁ e₂ : ι → ℝ) (ℓ₀ ℓ₁ ℓ₂ : ℝ) (a : ι) : ℝ :=
  e₀ a * ℓ₀ + e₁ a * ℓ₁ + e₂ a * ℓ₂

/--
Abstract star-basis weight calculation.

Assume a finite monomial family has total exponent `S` in each of three
coordinates.  Remove a root monomial of pure degree `d` in coordinate zero.
Then the remaining total weighted degree is

`S*(ℓ₀+ℓ₁+ℓ₂) - d*ℓ₀`.
-/
theorem sum_erase_pure_root_zero
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (e₀ e₁ e₂ : ι → ℝ) (root : ι)
    (S d ℓ₀ ℓ₁ ℓ₂ : ℝ)
    (h₀ : ∑ a, e₀ a = S)
    (h₁ : ∑ a, e₁ a = S)
    (h₂ : ∑ a, e₂ a = S)
    (hr₀ : e₀ root = d)
    (hr₁ : e₁ root = 0)
    (hr₂ : e₂ root = 0) :
    ∑ a in Finset.univ.erase root, weight3 e₀ e₁ e₂ ℓ₀ ℓ₁ ℓ₂ a =
      S * (ℓ₀ + ℓ₁ + ℓ₂) - d * ℓ₀ := by
  classical
  have hroot :
      weight3 e₀ e₁ e₂ ℓ₀ ℓ₁ ℓ₂ root = d * ℓ₀ := by
    simp [weight3, hr₀, hr₁, hr₂]
  have htotal :
      ∑ a, weight3 e₀ e₁ e₂ ℓ₀ ℓ₁ ℓ₂ a =
        S * (ℓ₀ + ℓ₁ + ℓ₂) := by
    unfold weight3
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
    rw [← Finset.sum_mul, ← Finset.sum_mul, ← Finset.sum_mul]
    rw [h₀, h₁, h₂]
    ring
  rw [← Finset.sum_erase_add _ _ (Finset.mem_univ root), hroot] at htotal
  linarith

/--
If coordinate zero has maximal nonnegative local height, deleting its pure
root loses at most `d` times the total local height.
-/
theorem sum_erase_pure_root_zero_lower
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (e₀ e₁ e₂ : ι → ℝ) (root : ι)
    (S d ℓ₀ ℓ₁ ℓ₂ : ℝ)
    (h₀ : ∑ a, e₀ a = S)
    (h₁ : ∑ a, e₁ a = S)
    (h₂ : ∑ a, e₂ a = S)
    (hr₀ : e₀ root = d)
    (hr₁ : e₁ root = 0)
    (hr₂ : e₂ root = 0)
    (hd : 0 ≤ d)
    (hℓ₀ : 0 ≤ ℓ₀) (hℓ₁ : 0 ≤ ℓ₁) (hℓ₂ : 0 ≤ ℓ₂) :
    (S - d) * (ℓ₀ + ℓ₁ + ℓ₂) ≤
      ∑ a in Finset.univ.erase root, weight3 e₀ e₁ e₂ ℓ₀ ℓ₁ ℓ₂ a := by
  rw [sum_erase_pure_root_zero e₀ e₁ e₂ root S d ℓ₀ ℓ₁ ℓ₂
    h₀ h₁ h₂ hr₀ hr₁ hr₂]
  have hsum : ℓ₀ ≤ ℓ₀ + ℓ₁ + ℓ₂ := by linarith
  nlinarith

/-- The numerical sextic specialization: `56 - 6 = 50`. -/
theorem sextic_weight_fifty
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    (e₀ e₁ e₂ : ι → ℝ) (root : ι)
    (ℓ₀ ℓ₁ ℓ₂ : ℝ)
    (h₀ : ∑ a, e₀ a = 56)
    (h₁ : ∑ a, e₁ a = 56)
    (h₂ : ∑ a, e₂ a = 56)
    (hr₀ : e₀ root = 6)
    (hr₁ : e₁ root = 0)
    (hr₂ : e₂ root = 0)
    (hℓ₀ : 0 ≤ ℓ₀) (hℓ₁ : 0 ≤ ℓ₁) (hℓ₂ : 0 ≤ ℓ₂) :
    50 * (ℓ₀ + ℓ₁ + ℓ₂) ≤
      ∑ a in Finset.univ.erase root, weight3 e₀ e₁ e₂ ℓ₀ ℓ₁ ℓ₂ a := by
  have h := sum_erase_pure_root_zero_lower
    e₀ e₁ e₂ root 56 6 ℓ₀ ℓ₁ ℓ₂
    h₀ h₁ h₂ hr₀ hr₁ hr₂ (by norm_num) hℓ₀ hℓ₁ hℓ₂
  norm_num at h ⊢
  exact h

end SexticTropicalWeight
end GaussianChain
