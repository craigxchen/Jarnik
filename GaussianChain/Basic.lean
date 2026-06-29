import Mathlib

namespace GaussianChain

def norm2 (a b : ℤ) : ℤ := a ^ 2 + b ^ 2

def det2 (m1 n1 m2 n2 : ℤ) : ℤ := m1 * n2 - m2 * n1

def kCoord (A B m n : ℤ) : ℤ := -(A * m + B * n)

def sCoord (A B m n : ℤ) : ℤ := A * n - B * m

theorem norm_k_s (A B m n : ℤ) :
    kCoord A B m n ^ 2 + sCoord A B m n ^ 2 =
      norm2 A B * norm2 m n := by
  unfold kCoord sCoord norm2
  ring

theorem sameCircle_norm2
    {A B m n h : ℤ}
    (hsame : norm2 (h * A + m) (h * B + n) = norm2 (h * A) (h * B)) :
    norm2 m n = 2 * h * kCoord A B m n := by
  unfold norm2 kCoord at *
  ring_nf at hsame ⊢
  omega

theorem rooted_circle_equation
    {A B m n h : ℤ}
    (N : ℤ)
    (hN : N = norm2 A B)
    (hsame : norm2 (h * A + m) (h * B + n) = norm2 (h * A) (h * B)) :
    kCoord A B m n ^ 2 + sCoord A B m n ^ 2 =
      2 * h * N * kCoord A B m n := by
  rw [norm_k_s, hN, sameCircle_norm2 hsame]
  ring

theorem determinant_identity (A B m1 n1 m2 n2 : ℤ) :
    kCoord A B m2 n2 * sCoord A B m1 n1 -
        kCoord A B m1 n1 * sCoord A B m2 n2 =
      norm2 A B * det2 m1 n1 m2 n2 := by
  unfold kCoord sCoord norm2 det2
  ring

end GaussianChain
