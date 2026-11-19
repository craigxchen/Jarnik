"""Real (not rational) strict solutions of the Im-moment system on a real line:
Im sum_{e ni i} r_e^k independent of i, for k=1..K (K=4 or 5), roots r_e = x_e + i y_e,
strict = 20 numbers r_e, conj(r_e) pairwise distinct (no real root).
K=5 is the strict linear-block system (constant cofactors) -- open in the record.
K=4 corresponds to all ten cofactor polynomials of degree <= 1, i.e. nu = 1 (Delta ~ P^-5),
exactly Astra's VI.5 boundary.  Gauge: x_0 = 0, y_0 = 1 (real affine maps of T)."""
import numpy as np, sys
from itertools import combinations
from scipy.optimize import least_squares
edges = list(combinations(range(5), 2))
Inc = np.array([[int(i in e) for e in edges] for i in range(5)], dtype=float)
E = Inc[:4] - Inc[4]
K = int(sys.argv[1]); trials = int(sys.argv[2])
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 0)
def unpack(p):
    x = np.concatenate([[0.0], p[:9]]); y = np.concatenate([[1.0], p[9:]])
    return x + 1j*y
def F(p):
    r = unpack(p)
    return np.concatenate([E @ np.imag(r**k) for k in range(1, K+1)])
def strictness(r):
    pts = np.concatenate([r, np.conj(r)])
    d = min(abs(a-b) for a, b in combinations(pts, 2))
    return d / max(1.0, np.max(np.abs(pts)))
found = []
for t in range(trials):
    p0 = rng.normal(size=18) * 2
    s = least_squares(F, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    r = unpack(s.x)
    res = np.max(np.abs(F(s.x))) / max(1, np.max(np.abs(r))**K)
    st = strictness(r)
    if res < 1e-11:
        J = s.jac
        rank = np.linalg.matrix_rank(J, tol=1e-8)
        found.append((st, res, rank, r))
found.sort(key=lambda z: -z[0])
print(f'K={K}: converged {len(found)}/{trials}')
for st, res, rank, r in found[:6]:
    print(f'  strictness={st:.3e} residual={res:.1e} jac rank={rank}  roots=', np.round(r, 4))
