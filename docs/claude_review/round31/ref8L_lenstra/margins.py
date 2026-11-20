#!/usr/bin/env python3
"""Exact check of the numerical margins of VIII.1 / VIII.3 with c=6/69, eta0=1/381."""
from fractions import Fraction as Fr
from math import ceil
c = Fr(6, 69); eta0 = Fr(1, 381)
a = Fr(2, 10**5); lam = Fr(1, 10**16)
print("a < eta0:", a < eta0)
print("28 a/4 < c/256:", 28 * a / 4 < c / 256, float(28*a/4), float(c/256))
print("28(a/4+lam) < c/256:", 28 * (a / 4 + lam) < c / 256)
# 127 sqrt(8 lam) < a/4  <=>  127^2 * 8 lam < a^2/16
print("127 sqrt(8lam) < a/4:", 127**2 * 8 * lam < a * a / 16, (127**2*8*lam/(a*a/16)))
# admissible M: 127^2 (2/M + 8 lam) < a^2/16
Mmin = 2 / (a * a / 16 / 127**2 - 8 * lam)
print("need M >", float(Mmin))
# largest lambda allowed for this a, and largest over a
print("lam max for a=2e-5:", float(a*a/16/127**2/8))
amax = c/256*4/28
print("a max:", float(amax), " lam max overall ~", float(amax*amax/16/127**2/8))
# VIII.2: eps/2 <= lam ; VIII.3 second assertion: kappa <= 2 lam/3 suffices; Astra kappa < lam/2
eps = Fr(1, 10**16); kappa = Fr(1, 10**17)
print("eps/2 <= lam:", eps/2 <= lam, " kappa < lam/2:", kappa < lam/2, " kappa<=2lam/3:", kappa <= 2*lam/3)
# Step 4 of ptolemy: (7-381 eta)w <= 69 H + C'; eta<=1/381 gives H >= (6w-C')/69
print("7-381*eta0 =", 7 - 381*eta0)
# second assertion: partition count. intervals [n^{(k-1)kappa}, n^{k kappa}], k=1..ceil(1/kappa)
for kap in [Fr(1,3), Fr(1,7), Fr(2,9), Fr(1,10)]:
    print("kappa", kap, "intervals", ceil(1/kap), "<= 1+ceil(1/kappa)")
# monotonicity of f(rho)=sqrt(rho)-1/sqrt(rho) on rho>=1 and bound (B-A)/sqrt(AB) <= n^{kappa/2}
import random
rng = random.Random(1)
for _ in range(10000):
    A = rng.uniform(1, 1e6); B = A * rng.uniform(1, 50)
    di = rng.uniform(A, B); dj = rng.uniform(A, B)
    lo, hi = min(di, dj), max(di, dj)
    assert (hi - lo) / (lo*hi)**0.5 <= (B - A)/(A*B)**0.5 * (1 + 1e-12)
    assert (B - A)/(A*B)**0.5 <= (B/A)**0.5
print("monotonicity sampled ok")
