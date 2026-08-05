import GaussianChain.PrimeRevealFaces

namespace GaussianChain.PrimeReveal

/-- The additive mixed difference on a two-prime reveal face:
`C₊₊ - C₊₋ - C₋₊ + C₋₋`.
Unlike the multiplicative cross-product, this is linear in the completion
count and therefore compatible with Fourier and convolution arguments. -/
def faceAdditiveCurvature {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i j : Fin n) : ℤ :=
  (faceCount A a i j true true : ℤ) -
    (faceCount A a i j true false : ℤ) -
    (faceCount A a i j false true : ℤ) +
    (faceCount A a i j false false : ℤ)

/-- Additive face curvature is symmetric in the two coordinates. -/
theorem faceAdditiveCurvature_swap {n : ℕ} (A : Finset (SignVector n))
    (a : PartialSignVector n) (i j : Fin n) (h : i ≠ j) :
    faceAdditiveCurvature A a i j = faceAdditiveCurvature A a j i := by
  simp only [faceAdditiveCurvature]
  rw [faceCount_swap_coordinates A a i j true true h,
    faceCount_swap_coordinates A a i j true false h,
    faceCount_swap_coordinates A a i j false true h,
    faceCount_swap_coordinates A a i j false false h]
  ring

/-- The additive curvature is the difference of the `i`-biases in the two
`j`-children. -/
theorem faceAdditiveCurvature_eq_bias_difference {n : ℕ}
    (A : Finset (SignVector n)) (a : PartialSignVector n) (i j : Fin n) :
    faceAdditiveCurvature A a i j =
      ((faceCount A a i j true true : ℤ) -
        faceCount A a i j false true) -
      ((faceCount A a i j true false : ℤ) -
        faceCount A a i j false false) := by
  simp [faceAdditiveCurvature]
  ring

/-- Equivalently, it is the difference of the `j`-biases in the two
`i`-children. -/
theorem faceAdditiveCurvature_eq_transverse_bias_difference {n : ℕ}
    (A : Finset (SignVector n)) (a : PartialSignVector n) (i j : Fin n) :
    faceAdditiveCurvature A a i j =
      ((faceCount A a i j true true : ℤ) -
        faceCount A a i j true false) -
      ((faceCount A a i j false true : ℤ) -
        faceCount A a i j false false) := by
  simp [faceAdditiveCurvature]
  ring

/-- Vanishing additive curvature means that revealing one coordinate does
not change the signed completion-count bias of the other coordinate. -/
theorem faceAdditiveCurvature_eq_zero_iff {n : ℕ}
    (A : Finset (SignVector n)) (a : PartialSignVector n) (i j : Fin n) :
    faceAdditiveCurvature A a i j = 0 ↔
      (faceCount A a i j true true : ℤ) -
          faceCount A a i j false true =
        (faceCount A a i j true false : ℤ) -
          faceCount A a i j false false := by
  rw [faceAdditiveCurvature_eq_bias_difference]
  exact sub_eq_zero

end GaussianChain.PrimeReveal
