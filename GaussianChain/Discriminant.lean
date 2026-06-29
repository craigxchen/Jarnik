import GaussianChain.Extension

namespace GaussianChain

def discLeading (C : Certificate) : ℤ :=
  (C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q) ^ 2 -
    8 * C.h * C.N * C.y ^ 2 * C.z * C.V ^ 2

def discLinear (C : Certificate) : ℤ :=
  -4 * C.h * C.y * C.z ^ 2 * C.V *
    (C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q)

def discConst (C : Certificate) : ℤ :=
  (2 * C.h * C.y * C.z ^ 2 * C.V) ^ 2

theorem extDisc_expand (C : Certificate) (W : ℤ) :
    extDisc C W = discLeading C * W ^ 2 + discLinear C * W + discConst C := by
  unfold extDisc extB extA extC discLeading discLinear discConst
  ring

theorem GeomCertificate.discSlope_mulJ
    (C : GeomCertificate) :
    C.J * (C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q) =
      -2 * C.y * C.z * C.V * C.n3num := by
  calc
    C.J * (C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q)
        = (C.N * C.J) * C.V ^ 2 -
            2 * (C.h * C.J) * C.y * C.z * C.q := by ring
    _ = (2 * C.x * C.y * C.z * C.c) * C.V ^ 2 -
            2 * (C.U * C.V * C.w) * C.y * C.z * C.q := by
          rw [C.hN, C.hJ]
    _ = -2 * C.y * C.z * C.V *
            (C.q * C.U * C.w - C.x * C.c * C.V) := by ring
    _ = -2 * C.y * C.z * C.V * C.n3num := by
          rw [C.hn3]
          unfold Certificate.n3
          ring

theorem GeomCertificate.discLeading_mulJ2
    (C : GeomCertificate) :
    C.J ^ 2 * discLeading C.toCertificate =
      C.J ^ 2 * (-4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2) := by
  let D : ℤ := C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q
  have hD : C.J * D = -2 * C.y * C.z * C.V * C.n3num := by
    dsimp [D]
    exact C.discSlope_mulJ
  have hsquare :
      C.n3num ^ 2 + C.J ^ 2 - 4 * C.x * C.y * C.c * C.U * C.V * C.w = 0 := by
    linarith [C.hsquare3]
  have htail :
      (-2 * C.y * C.z * C.V * C.n3num) ^ 2 -
          8 * (C.U * C.V * C.w) * (2 * C.x * C.y * C.z * C.c) *
            C.y ^ 2 * C.z * C.V ^ 2 =
        -4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2 * C.J ^ 2 := by
    have hid :
        ((-2 * C.y * C.z * C.V * C.n3num) ^ 2 -
            8 * (C.U * C.V * C.w) * (2 * C.x * C.y * C.z * C.c) *
              C.y ^ 2 * C.z * C.V ^ 2) -
          (-4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2 * C.J ^ 2) =
            4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2 *
              (C.n3num ^ 2 + C.J ^ 2 -
                4 * C.x * C.y * C.c * C.U * C.V * C.w) := by
      ring
    rw [hsquare] at hid
    nlinarith
  calc
    C.J ^ 2 * discLeading C.toCertificate
        = (C.J * D) ^ 2 -
            8 * (C.h * C.J) * (C.N * C.J) *
              C.y ^ 2 * C.z * C.V ^ 2 := by
          dsimp [D]
          unfold discLeading
          ring
    _ = (-2 * C.y * C.z * C.V * C.n3num) ^ 2 -
            8 * (C.U * C.V * C.w) * (2 * C.x * C.y * C.z * C.c) *
              C.y ^ 2 * C.z * C.V ^ 2 := by
          rw [hD, C.hJ, C.hN]
    _ = -4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2 * C.J ^ 2 := htail
    _ = C.J ^ 2 * (-4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2) := by ring

theorem GeomCertificate.discLeading_eq
    (C : GeomCertificate) :
    discLeading C.toCertificate = -4 * C.y ^ 2 * C.z ^ 2 * C.V ^ 2 := by
  have hmul := C.discLeading_mulJ2
  have hJ2 : C.J ^ 2 ≠ 0 := pow_ne_zero 2 C.hJ_ne_zero
  exact mul_left_cancel₀ hJ2 hmul

end GaussianChain
