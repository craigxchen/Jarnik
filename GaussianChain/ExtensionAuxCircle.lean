import GaussianChain.Discriminant

namespace GaussianChain

/-- Completing the square for a quadratic with negative-square leading coefficient.
If `Y² = -A² W² + L W + C`, then `(W,Y)` maps to an integral point on a
fixed circle. -/
theorem negativeSquareQuadratic_auxCircle
    (A L C W Y : ℤ)
    (hY : Y ^ 2 = -(A ^ 2) * W ^ 2 + L * W + C) :
    (2 * A * Y) ^ 2 + (2 * A ^ 2 * W - L) ^ 2 =
      L ^ 2 + 4 * A ^ 2 * C := by
  nlinarith [hY]

/-- The square root of the negative leading coefficient in the geometric
extension discriminant. -/
def extensionCircleScale (C : Certificate) : ℤ :=
  2 * C.y * C.z * C.V

/-- The affine `X`-coordinate in the auxiliary circle attached to an extension
parameter `W`. -/
def extensionCircleX (C : Certificate) (W : ℤ) : ℤ :=
  2 * (extensionCircleScale C) ^ 2 * W - discLinear C

/-- The affine `Y`-coordinate in the auxiliary circle attached to a square root
of the extension discriminant. -/
def extensionCircleY (C : Certificate) (Y : ℤ) : ℤ :=
  2 * extensionCircleScale C * Y

/-- Every integral square value of the extension discriminant lies on one fixed
integral auxiliary circle.  The geometric certificate is used exactly to turn
the quadratic leading coefficient into a negative square. -/
theorem GeomCertificate.extensionDiscriminant_auxCircle
    (C : GeomCertificate) (W Y : ℤ)
    (hY : extDisc C.toCertificate W = Y ^ 2) :
    extensionCircleX C.toCertificate W ^ 2 +
        extensionCircleY C.toCertificate Y ^ 2 =
      discLinear C.toCertificate ^ 2 +
        4 * (extensionCircleScale C.toCertificate) ^ 2 *
          discConst C.toCertificate := by
  have hexpand := extDisc_expand C.toCertificate W
  have hlead := C.discLeading_eq
  have hquad :
      Y ^ 2 =
        -((extensionCircleScale C.toCertificate) ^ 2) * W ^ 2 +
          discLinear C.toCertificate * W + discConst C.toCertificate := by
    rw [hexpand, hlead] at hY
    simpa [extensionCircleScale] using hY.symm
  simpa [extensionCircleX, extensionCircleY] using
    (negativeSquareQuadratic_auxCircle
      (extensionCircleScale C.toCertificate)
      (discLinear C.toCertificate)
      (discConst C.toCertificate) W Y hquad)

/-- An actual extension supplies the square root required by the auxiliary-circle
identity. -/
theorem GeomCertificate.extension_auxCircle
    (C : GeomCertificate) (E : Extension C.toCertificate) :
    extensionCircleX C.toCertificate E.W ^ 2 +
        extensionCircleY C.toCertificate
          (2 * extA C.toCertificate * E.t + extB C.toCertificate E.W) ^ 2 =
      discLinear C.toCertificate ^ 2 +
        4 * (extensionCircleScale C.toCertificate) ^ 2 *
          discConst C.toCertificate := by
  apply C.extensionDiscriminant_auxCircle
  exact E.discriminant_square.symm

/-- If the third geometric numerator is an integral multiple `J*r`, then the
slope term in the extension quadratic factors by the same integer `r`. -/
theorem GeomCertificate.discSlope_eq_of_n3_eq_J_mul
    (C : GeomCertificate) (r : ℤ) (hr : C.n3num = C.J * r) :
    C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q =
      -2 * C.y * C.z * C.V * r := by
  have hslope := C.discSlope_mulJ
  rw [hr] at hslope
  apply mul_left_cancel₀ C.hJ_ne_zero
  calc
    C.J * (C.N * C.V ^ 2 - 2 * C.h * C.y * C.z * C.q)
        = -2 * C.y * C.z * C.V * (C.J * r) := hslope
    _ = C.J * (-2 * C.y * C.z * C.V * r) := by ring

/-- The same integral slope `r` lies on the reduced rooted circle
`z (r²+1) = 2 h N`. -/
theorem GeomCertificate.slope_norm_eq_of_n3_eq_J_mul
    (C : GeomCertificate) (r : ℤ) (hr : C.n3num = C.J * r) :
    C.z * (r ^ 2 + 1) = 2 * C.h * C.N := by
  have hsquare := C.hsquare3
  rw [hr] at hsquare
  have hJrank : C.J * (r ^ 2 + 1) = 4 * C.h * C.x * C.y * C.c := by
    apply mul_left_cancel₀ C.hJ_ne_zero
    calc
      C.J * (C.J * (r ^ 2 + 1))
          = (C.J * r) ^ 2 + C.J ^ 2 := by ring
      _ = 4 * C.x * C.y * C.c * C.U * C.V * C.w := hsquare
      _ = C.J * (4 * C.h * C.x * C.y * C.c) := by
        rw [← C.hJ]
        ring
  apply mul_left_cancel₀ C.hJ_ne_zero
  calc
    C.J * (C.z * (r ^ 2 + 1))
        = C.z * (C.J * (r ^ 2 + 1)) := by ring
    _ = C.z * (4 * C.h * C.x * C.y * C.c) := by rw [hJrank]
    _ = 2 * C.h * (C.N * C.J) := by rw [C.hN]; ring
    _ = C.J * (2 * C.h * C.N) := by ring

/-- Reduced integral coordinate associated with an extension parameter, once the
geometric slope `n3/J` is integral. -/
def reducedExtensionX (C : Certificate) (r W : ℤ) : ℤ :=
  W - C.h * C.z * r

/-- Reduced integral coordinate associated with an extension root. -/
def reducedExtensionY (C : Certificate) (r W t : ℤ) : ℤ :=
  C.h * (2 * t - C.z) - r * W

/-- Critical renormalization identity.

When `n3 = J*r`, every actual extension maps to an integral point on the much
smaller circle

`X² + Y² = 2 h³ N z`.

The nonvanishing hypothesis is exactly what is needed to remove the common
factor `y*V` from the extension quadratic. -/
theorem GeomCertificate.extension_reducedCircle
    (C : GeomCertificate) (r : ℤ) (hr : C.n3num = C.J * r)
    (hyV : C.y * C.V ≠ 0) (E : Extension C.toCertificate) :
    reducedExtensionX C.toCertificate r E.W ^ 2 +
        reducedExtensionY C.toCertificate r E.W E.t ^ 2 =
      2 * C.h ^ 3 * C.N * C.z := by
  have hD := C.discSlope_eq_of_n3_eq_J_mul r hr
  have hrnorm := C.slope_norm_eq_of_n3_eq_J_mul r hr
  have hquad := E.quadratic
  have hid :
      extPoly C.toCertificate E.W E.t =
        (C.y * C.V) *
          (2 * C.h * C.z * E.t ^ 2 -
            2 * C.z * (r * E.W + C.h * C.z) * E.t +
            C.N * E.W ^ 2) := by
    unfold extPoly extA extB extC
    rw [hD]
    ring
  have hmul :
      (C.y * C.V) *
          (2 * C.h * C.z * E.t ^ 2 -
            2 * C.z * (r * E.W + C.h * C.z) * E.t +
            C.N * E.W ^ 2) = 0 := by
    rw [← hid, hquad]
  have hreduced :
      2 * C.h * C.z * E.t ^ 2 -
          2 * C.z * (r * E.W + C.h * C.z) * E.t +
          C.N * E.W ^ 2 = 0 :=
    (mul_eq_zero.mp hmul).resolve_left hyV
  unfold reducedExtensionX reducedExtensionY
  nlinarith [hreduced, hrnorm]

end GaussianChain
