import numpy as np, sys
from rchart5b import *
rng = np.random.default_rng(int(sys.argv[1]))
ntr = int(sys.argv[2])
c = lambda n: (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)
# fixed rational-looking real parameters
rho = dict(zip(triples, np.array([-5,-3,-2,-1,1,2,3,5,7,-7], dtype=complex)*0.5))
a1 = 0.0+0j
sols = []; degen = 0; fail = 0
for tr in range(ntr):
    z0 = np.concatenate([c(4) * 2, c(4), c(12)])
    z, nf = newton(z0, a1, rho)
    if nf > 1e-10: fail += 1; continue
    info, special = analyse(z, a1, rho)
    if info is None or min(info) < 1e-6: degen += 1; continue
    if not any(np.linalg.norm(z - s) < 1e-6 for s in sols):
        sols.append(z)
print("distinct nondegenerate solutions:", len(sols), "degenerate", degen, "fail", fail)
for s in sols:
    print(np.round(s[:4], 6), " max|Im| =", np.max(np.abs(s.imag)))
