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

end GaussianChain
