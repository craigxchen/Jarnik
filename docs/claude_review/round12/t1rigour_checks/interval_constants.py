"""Interval-arithmetic certification of the constants of Lemma 2.2 / Cor 2.3 (mpmath.iv, 40 digits)."""
from mpmath import iv
iv.dps = 40
c = iv.mpf('0.6417')
g = (1 + c) * iv.log(1 + c) - c * iv.log(c)
print("g(0.6417) - ln 3 in", g - iv.log(3), "(must be < 0)")
kappa = iv.log(8) / iv.log(1 + 1 / c)
print("kappa = ln8/ln(1+1/c*) in", kappa, "(< 2.2137)")
phi = c * 316 - 4 * iv.sqrt(iv.mpf(8 * 316))
print("phi(316) = c*316 - 4 sqrt(2528) in", phi, "; ln4*phi(316) in", iv.log(4) * phi, "(must exceed kappa)")
print("24 ln 3 =", 24 * iv.log(3), " 6 ln3/sqrt8 =", 6 * iv.log(3) / iv.sqrt(8))
