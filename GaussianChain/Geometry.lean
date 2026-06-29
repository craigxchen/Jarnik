import GaussianChain.Basic

namespace GaussianChain

structure RootedPoint where
  A : ℤ
  B : ℤ
  h : ℤ
  m : ℤ
  n : ℤ

namespace RootedPoint

def N (P : RootedPoint) : ℤ := norm2 P.A P.B

def k (P : RootedPoint) : ℤ := kCoord P.A P.B P.m P.n

def s (P : RootedPoint) : ℤ := sCoord P.A P.B P.m P.n

def sameCircle (P : RootedPoint) : Prop :=
  norm2 (P.h * P.A + P.m) (P.h * P.B + P.n) =
    norm2 (P.h * P.A) (P.h * P.B)

theorem square_relation (P : RootedPoint) (hcircle : P.sameCircle) :
    P.k ^ 2 + P.s ^ 2 = 2 * P.h * P.N * P.k := by
  unfold k s N sameCircle at *
  exact rooted_circle_equation (norm2 P.A P.B) rfl hcircle

end RootedPoint

end GaussianChain
