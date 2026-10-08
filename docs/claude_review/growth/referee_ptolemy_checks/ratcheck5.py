import numpy as np, sys
from fractions import Fraction
from rchart5b import *
rng = np.random.default_rng(11)
c = lambda n: (rng.normal(size=n)) * 1.0
rho = dict(zip(triples, np.array([-5,-3,-2,-1,1,2,3,5,7,-7], dtype=float)*0.5))
a1 = 0.0
sols=[]
import rchart5b
for tr in range(400):
    z0 = np.concatenate([c(4) * 2, c(4), c(12)]).astype(complex)
    z, nf = newton(z0, a1, rho)
    if nf > 1e-11: continue
    info, special = analyse(z, a1, rho)
    if info is None or min(info) < 1e-6: continue
    if not any(np.linalg.norm(z - s) < 1e-6 for s in sols):
        sols.append(z)
print(len(sols), "distinct nondegenerate real-start solutions")
for s in sols:
    fr = [Fraction(float(x.real)).limit_denominator(2000) for x in s[:8]]
    err = max(abs(float(f) - x.real) for f, x in zip(fr, s[:8]))
    print([str(f) for f in fr[:4]], "approx err", err)
