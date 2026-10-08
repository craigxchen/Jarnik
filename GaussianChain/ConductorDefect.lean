import Mathlib

namespace GaussianChain
namespace ConductorDefect

open scoped BigOperators

/--
The contribution of one split-prime valuation layer to the cancelled
Ramana determinant for `2 * s + 1` points.  Here `r` is the number of
points on one side of the layer cut.
-/
def layerDefect (s r : ℤ) : ℤ :=
  s * (s + 1) - r * (2 * s + 1 - r)

/-- The layer defect is a product of two consecutive integers. -/
theorem layerDefect_eq_consecutive (s r : ℤ) :
    layerDefect s r = (r - s) * (r - s - 1) := by
  unfold layerDefect
  ring

/-- Every integral layer defect is nonnegative. -/
theorem layerDefect_nonneg (s r : ℤ) :
    0 ≤ layerDefect s r := by
  rw [layerDefect_eq_consecutive]
  by_cases h : r ≤ s
  · exact mul_nonneg_of_nonpos_of_nonpos (by omega) (by omega)
  · exact mul_nonneg (by omega) (by omega)

/-- A layer has zero defect exactly when its two sides have sizes `s` and `s + 1`. -/
theorem layerDefect_eq_zero_iff (s r : ℤ) :
    layerDefect s r = 0 ↔ r = s ∨ r = s + 1 := by
  rw [layerDefect_eq_consecutive, mul_eq_zero]
  constructor
  · rintro (h | h)
    · left
      omega
    · right
      omega
  · rintro (h | h)
    · left
      omega
    · right
      omega

/-- Every layer defect is even. -/
theorem layerDefect_even (s r : ℤ) :
    Even (layerDefect s r) := by
  rw [layerDefect_eq_consecutive]
  exact Int.even_mul_pred_self (r - s)

/--
A sum of layer defects vanishes exactly when every layer is maximally
balanced.  This is the equality case needed when auditing conductor-only
determinant arguments.
-/
theorem sum_layerDefect_eq_zero_iff
    {ι : Type*} [Fintype ι] (s : ℤ) (r : ι → ℤ) :
    (∑ i, layerDefect s (r i)) = 0 ↔
      ∀ i, r i = s ∨ r i = s + 1 := by
  rw [Fintype.sum_eq_zero_iff_of_nonneg (fun i ↦ layerDefect_nonneg s (r i))]
  constructor
  · intro h i
    exact (layerDefect_eq_zero_iff s (r i)).mp (congrFun h i)
  · intro h
    funext i
    exact (layerDefect_eq_zero_iff s (r i)).mpr (h i)

/-- The nonzero vectors of `𝔽₂³`, indexed by the naturals `1,…,7`. -/
def simplexBit (x a : Fin 7) : Bool :=
  let x' := x.val + 1
  let a' := a.val + 1
  (((x' % 2) * (a' % 2) +
      ((x' / 2) % 2) * ((a' / 2) % 2) +
      ((x' / 4) % 2) * ((a' / 4) % 2)) % 2) == 1

/-- Each nonzero linear functional on `𝔽₂³` cuts the seven nonzero vectors `4+3`. -/
theorem simplexBit_column_card (a : Fin 7) :
    (Finset.univ.filter fun x : Fin 7 ↦ simplexBit x a).card = 4 := by
  fin_cases a <;> decide

/-- Any two rows of the seven-column simplex pattern differ in exactly four columns. -/
theorem simplexBit_pair_distance
    (x y : Fin 7) (hxy : x ≠ y) :
    (Finset.univ.filter fun a : Fin 7 ↦ simplexBit x a != simplexBit y a).card = 4 := by
  revert hxy
  fin_cases x <;> fin_cases y <;> decide

/-- The sign matrix attached to the seven-point simplex incidence matrix. -/
def simplexSign (x a : Fin 7) : ℤ :=
  if simplexBit x a then -1 else 1

/--
The integer matrix identity `B * (J - 2B) = -4I` for the punctured
three-dimensional Walsh matrix.  It is the exact inverse used in the
seven-row angular argument.
-/
theorem simplexBit_mul_simplexSign :
    ∀ x y : Fin 7,
      (∑ a : Fin 7, if simplexBit x a then simplexSign a y else 0) =
        if x = y then -4 else 0 := by
  decide

/-- The row sum of a vector over the simplex incidence matrix, modulo four. -/
def simplexIncidenceSumModFour (j : Fin 7 → ZMod 4) (x : Fin 7) : ZMod 4 :=
  ∑ a : Fin 7, if simplexBit x a then j a else 0

/--
Every vector in the mod-four kernel of the seven-point simplex incidence
matrix has constant parity.  This is the branch reduction in the angular
argument.
-/
theorem simplex_modFour_kernel_constant_parity :
    ∀ j : Fin 7 → ZMod 4,
      (∀ x, simplexIncidenceSumModFour j x = 0) →
        ∀ a b, (j a).val % 2 = (j b).val % 2 := by
  native_decide

/-- The parity of one incidence-row sum for a Boolean support. -/
def simplexIncidenceParity (e : Fin 7 → Bool) (x : Fin 7) : Nat :=
  ((Finset.univ.filter fun a : Fin 7 ↦ simplexBit x a && e a).card) % 2

/--
A zero-sum subset of the seven nonzero vectors of `𝔽₂³` has size
`0`, `3`, `4`, or `7`.  Hence its two phase branches always contain three
disjoint same-branch pairs.
-/
theorem simplex_kernel_support_card :
    ∀ e : Fin 7 → Bool,
      (∀ x, simplexIncidenceParity e x = 0) →
        let n := (Finset.univ.filter fun a : Fin 7 ↦ e a).card
        n = 0 ∨ n = 3 ∨ n = 4 ∨ n = 7 := by
  native_decide

end ConductorDefect
end GaussianChain
