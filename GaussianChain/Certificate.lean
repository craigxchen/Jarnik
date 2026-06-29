import GaussianChain.Basic

namespace GaussianChain

structure Certificate where
  x : ℤ
  y : ℤ
  z : ℤ
  p : ℤ
  q : ℤ
  r : ℤ
  U : ℤ
  V : ℤ
  w : ℤ
  c : ℤ
  J : ℤ
  h : ℤ
  N : ℤ
  hp : p = y - x
  hq : q = z - y
  hr : r = z - x
  hw : y * w = z * U + x * V
  hc : y * c = p * V - q * U
  hJ : h * J = U * V * w
  hN : N * J = 2 * x * y * z * c

namespace Certificate

def n2Left (C : Certificate) : ℤ := C.p * C.V * C.w - C.z * C.c * C.U

def n2Right (C : Certificate) : ℤ := C.q * C.U * C.w + C.x * C.c * C.V

def n3 (C : Certificate) : ℤ := C.q * C.U * C.w - C.x * C.c * C.V

def endpointValid (C : Certificate) : Prop := 4 * C.z ^ 2 ≤ C.N

end Certificate

structure GeomCertificate extends Certificate where
  n2 : ℤ
  n3num : ℤ
  hn2_left : n2 = toCertificate.n2Left
  hn2_right : n2 = toCertificate.n2Right
  hn3 : n3num = toCertificate.n3
  hsquare3 : n3num ^ 2 + J ^ 2 = 4 * x * y * c * U * V * w
  hJ_ne_zero : J ≠ 0

end GaussianChain
