"""M8 (Syl(8), L=2, m=8) and M8' (L=3, m=11): exact basis of V(A,m), then numerical search for real points of
TC(A,m) = {Phi_c = 0 for all c in V} and classification of the rational relations they satisfy.
Question (barrier.md §6): does TC contain collision-free real points (necessarily on light 4-term hyperplanes, by J4)?
Exact part: V(A,m) over Q (flint fmpz_mat nullspace).  Search part: EVIDENCE (floating point, mpmath refinement).
Usage: python3 m8_tc.py L m nstarts
"""
import sys, time, itertools
import numpy as np
import flint
from flib import sylvester_patterns, box, reduced_monomials

L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
m = int(sys.argv[2]) if len(sys.argv) > 2 else 8
NST = int(sys.argv[3]) if len(sys.argv) > 3 else 400
t0 = time.time()
M, pats = sylvester_patterns(3)
A = box(M, pats, L)
d = M - 1
print(f"Syl(8) L={L} m={m}: |A|={len(A)} density={len(A)/3**d:.4f}", flush=True)
mons = [mm for deg in range(m) for mm in reduced_monomials(d, deg)]
E = flint.fmpz_mat([[int(np.prod([k[i] ** e[i] for i in range(d)])) for k in A] for e in mons])  # N x |A|
X, nul = E.nullspace()
basis = [[int(X[i, j]) for i in range(X.nrows())] for j in range(nul)]
print(f"dim V(A,{m}) = {nul}  (monomials {len(mons)}), {time.time()-t0:.1f}s", flush=True)
Amat = np.array(A, dtype=float)
C = np.array(basis, dtype=float)
C = C / np.abs(C).sum(axis=1, keepdims=True)


def resid(u):
    v = Amat @ u
    return C @ (v ** m) / d ** m


def jac(u):
    v = Amat @ u
    return (C * (m * v ** (m - 1))) @ Amat / d ** m


from scipy.optimize import least_squares
rng = np.random.default_rng(1)
found = []
for s in range(NST):
    u0 = rng.uniform(-1, 1, d)
    f = lambda y: np.concatenate([resid(y) * 1e6, [np.dot(y, y) - 1.0]])
    J = lambda y: np.vstack([jac(y) * 1e6, 2 * y[None, :]])
    r = least_squares(f, u0 / np.linalg.norm(u0), jac=J, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    u = r.x / np.abs(r.x).max()
    kap = np.abs(resid(u)).max()
    z = np.concatenate([[0.0], u])
    gaps = [abs(z[a] - z[b]) for a in range(M) for b in range(a + 1, M)]
    found.append((kap, min(gaps), u))
found.sort(key=lambda t: t[0])
print("smallest kappa_inf values found:", [f"{t[0]:.2e}" for t in found[:10]])
zeros = [t for t in found if t[0] < 1e-14]
print(f"{len(zeros)} numerical zeros (kappa<1e-14) of {NST} starts")
# classify relations: small integer relations lambda.z = 0 with ||lambda||_1 <= 4 among z = (0,u)
quads = []
for x, y, zz, w in itertools.permutations(range(M), 4):
    if x < y and zz < w and x < zz and (x ^ y) == (zz ^ w):
        quads.append((x, y, zz, w))
cf = 0
for kap, g, u in zeros:
    z = np.concatenate([[0.0], u])
    coll = [(a, b) for a in range(M) for b in range(a + 1, M) if abs(z[a] - z[b]) < 1e-7]
    qrel = [q for q in quads if abs(z[q[0]] + z[q[1]] - z[q[2]] - z[q[3]]) < 1e-7]
    if not coll:
        cf += 1
        print(f"  COLLISION-FREE zero: kappa={kap:.1e} mingap={g:.3e} quad relations={qrel} u={np.round(u,6)}")
print(f"collision-free numerical zeros: {cf}")
classes = {}
for kap, g, u in zeros:
    z = np.concatenate([[0.0], u])
    vals = np.round(z, 6)
    key = len(set(vals))
    classes[key] = classes.get(key, 0) + 1
print("number of distinct values among the 8 positions, over zeros:", dict(sorted(classes.items())))
print(f"total {time.time()-t0:.1f}s")
