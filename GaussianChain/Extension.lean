import GaussianChain.Certificate

namespace GaussianChain

def extA (C : Certificate) : ℤ := 2 * C.h * C.y * C.z * C.V

def extB (C : Certificate) (W : ℤ) : ℤ :=
  C.N * C.V ^ 2 * W - 2 * C.h * C.y * C.z * (C.q * W + C.z * C.V)

def extC (C : Certificate) (W : ℤ) : ℤ := C.N * C.V * C.y * W ^ 2

def extPoly (C : Certificate) (W t : ℤ) : ℤ :=
  extA C * t ^ 2 + extB C W * t + extC C W

def extDisc (C : Certificate) (W : ℤ) : ℤ :=
  extB C W ^ 2 - 4 * extA C * extC C W

structure Extension (C : Certificate) where
  t : ℤ
  W : ℤ
  wp : ℤ
  cp : ℤ
  Jp : ℤ
  hwp : C.z * wp = t * C.V + C.y * W
  hcp : C.z * cp = C.q * W - (t - C.z) * C.V
  hJp : C.h * Jp = C.V * W * wp
  hNp : C.N * Jp = 2 * C.y * C.z * t * cp

theorem extension_quadratic_from_reconstruction
    (C : Certificate) (E : Extension C) :
    extPoly C E.W E.t = 0 := by
  have hmain : C.N * C.V * E.W * E.wp = 2 * C.h * C.y * C.z * E.t * E.cp := by
    calc
      C.N * C.V * E.W * E.wp = C.N * (C.h * E.Jp) := by
        rw [E.hJp]
        ring
      _ = C.h * (C.N * E.Jp) := by ring
      _ = C.h * (2 * C.y * C.z * E.t * E.cp) := by
        rw [E.hNp]
      _ = 2 * C.h * C.y * C.z * E.t * E.cp := by ring
  have hzmain :
      C.N * C.V * E.W * (C.z * E.wp) =
        2 * C.h * C.y * C.z * E.t * (C.z * E.cp) := by
    calc
      C.N * C.V * E.W * (C.z * E.wp) =
          C.z * (C.N * C.V * E.W * E.wp) := by ring
      _ = C.z * (2 * C.h * C.y * C.z * E.t * E.cp) := by
        rw [hmain]
      _ = 2 * C.h * C.y * C.z * E.t * (C.z * E.cp) := by ring
  have hzero :
      C.N * C.V * E.W * (E.t * C.V + C.y * E.W) -
          2 * C.h * C.y * C.z * E.t *
            (C.q * E.W - (E.t - C.z) * C.V) = 0 := by
    rw [← E.hwp, ← E.hcp]
    linarith [hzmain]
  have hid :
      extPoly C E.W E.t =
        C.N * C.V * E.W * (E.t * C.V + C.y * E.W) -
          2 * C.h * C.y * C.z * E.t *
            (C.q * E.W - (E.t - C.z) * C.V) := by
    unfold extPoly extA extB extC
    ring
  rw [hid]
  exact hzero

namespace Extension

theorem quadratic {C : Certificate} (E : Extension C) :
    extPoly C E.W E.t = 0 :=
  extension_quadratic_from_reconstruction C E

theorem discriminant_square {C : Certificate} (E : Extension C) :
    extDisc C E.W = (2 * extA C * E.t + extB C E.W) ^ 2 := by
  have hid :
      extDisc C E.W - (2 * extA C * E.t + extB C E.W) ^ 2 =
        -4 * extA C * extPoly C E.W E.t := by
    unfold extDisc extPoly
    ring
  rw [E.quadratic] at hid
  nlinarith

end Extension

def ExtensionAdmissible (C : Certificate) (E : Extension C) : Prop :=
  C.z < E.t ∧ 4 * E.t ^ 2 ≤ C.N ∧
    0 < E.wp ∧ 0 < E.cp ∧ 0 < E.Jp ∧
    (E.t - C.z) * C.V * E.wp > C.y * E.cp * E.W

end GaussianChain
