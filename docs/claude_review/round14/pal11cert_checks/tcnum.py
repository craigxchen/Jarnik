"""EVIDENCE ONLY (floating point).  Real zeros of the 436 tangent octics of pal11 (L=2, m=8) and the rational
subspaces containing them.

1. A real orthonormal basis C (12093 x 436) of V(A,8) = ker E_{<8} from the Z_5-isotypic blocks (complex SVD per
   block; real and imaginary parts of the lifted kernel vectors of blocks 1, 2; block 0 is real).
   Check: ||E_{<8} C|| small, rank 436.
2. F(u) = ||C^T y(u)||^2 / ||y(u)||^2, y(u)_k = <k',u>^8: F(u) = 0 iff u in TC(A,8).  Minimise F on the unit sphere
   from random starts (L-BFGS on u/|u|).
3. For each numerical zero (F < 1e-24), compute the rational relation lattice {z in Z^11 : <z,u> ~ 0} by LLL
   (flint) on [I | round(N u)], N = 1e9, keeping vectors with |<z,u>| < 1e-7 and small norm; report the dimension of
   the smallest rational subspace containing u (= 11 - rank of relations), and which known cheap spans it lies in.
usage: python3 tcnum.py nstarts seed"""
import sys, time, pickle
import numpy as np
from scipy.optimize import minimize
from flint import fmpz_mat
from fplib import pal_patterns, box, reduced_monomials

nstarts, seed = int(sys.argv[1]), int(sys.argv[2])
M_, pats = pal_patterns(11); A = box(M_, pats, 2); d = 11
Aarr = np.array(A, dtype=np.float64)
Aidx = {k: i for i, k in enumerate(A)}
perm = [(3 * i) % 11 for i in range(11)]
def act(k):
    o = [0] * d
    for i in range(d): o[perm[i]] = k[i]
    return tuple(o)
def orbits(items):
    seen = set(); reps = []; fixed = []
    for it in items:
        if it in seen: continue
        orb = [it]
        for _ in range(4): orb.append(act(orb[-1]))
        seen.update(orb)
        if len(set(orb)) == 1: fixed.append(it)
        else: reps.append(orb)
    return reps, fixed
t0 = time.time()
Pre, Pfix = orbits(A)
Pidx = np.array([[Aidx[o] for o in orb] for orb in Pre])
Fidx = np.array([Aidx[k] for k in Pfix])
low = [a for e in range(8) for a in reduced_monomials(d, e)]
Mre, Mfix = orbits(low)
def monovals(a):
    col = np.ones(len(A))
    for x, ax in enumerate(a):
        if ax == 1: col = col * Aarr[:, x]
        elif ax == 2: col = col * Aarr[:, x] ** 2
    return col
w = np.exp(2j * np.pi / 5)
blocks = [[] for _ in range(5)]
for orb in Mre:
    v = monovals(orb[0]); V = v[Pidx]
    for j in range(5):
        r = V @ (w ** (-j * np.arange(5)))
        if j == 0: r = np.concatenate([r.real, v[Fidx]])
        blocks[j].append(r)
for a in Mfix:
    v = monovals(a); V = v[Pidx]
    blocks[0].append(np.concatenate([5 * V[:, 0], v[Fidx]]))
cols = []
for j in (0, 1, 2):
    B = np.array(blocks[j])
    U_, S, Vh = np.linalg.svd(B, full_matrices=True)
    ns = Vh[np.sum(S > 1e-8 * S[0]):].conj()    # null space rows
    print('block', j, 'shape', B.shape, 'nullity', ns.shape[0], 'smallest kept sv', S[np.sum(S > 1e-8 * S[0]) - 1], flush=True)
    for kv in ns:
        c = np.zeros(len(A), dtype=complex)
        for o in range(len(Pre)):
            for t in range(5):
                c[Pidx[o, t]] += kv[o] * w ** (-j * t)
        if j == 0:
            for f in range(len(Pfix)): c[Fidx[f]] += kv[len(Pre) + f]
            cols.append(c.real)
        else:
            cols.append(c.real); cols.append(c.imag)
C = np.array(cols).T
Q, Rr = np.linalg.qr(C)
print('basis vectors', C.shape[1], 'min |diag R|', np.abs(np.diag(Rr)).min(), '(%.0fs)' % (time.time() - t0), flush=True)
# check kernel property on a random subset of low-degree monomials
idx = np.random.default_rng(0).choice(len(low), 400, replace=False)
Echk = np.array([monovals(low[i]) for i in idx])
print('max |E_<8 Q| on 400 random monomials:', np.abs(Echk @ Q).max(), flush=True)
C = Q
def F(u):
    n = np.linalg.norm(u); u = u / n
    s = Aarr @ u; y = s ** 8
    r = C.T @ y
    val = (r @ r) / (y @ y)
    # gradient wrt u (then project for normalisation)
    gy = 2 * (C @ r) / (y @ y) - 2 * val * y / (y @ y)
    g = Aarr.T @ (gy * 8 * s ** 7)
    g = (g - (g @ u) * u) / n
    return val, g
rng = np.random.default_rng(seed)
zeros = []
for it in range(nstarts):
    u0 = rng.standard_normal(d)
    res = minimize(F, u0, jac=True, method='L-BFGS-B', options={'maxiter': 3000, 'ftol': 1e-30, 'gtol': 1e-16})
    u = res.x / np.linalg.norm(res.x)
    zeros.append((res.fun, u))
pickle.dump(zeros, open('tcnum_zeros_%d.pkl' % seed, 'wb'))
def relations(u, N=10**9):
    n = len(u)
    rows = []
    for i in range(n):
        r = [0] * (n + 1); r[i] = 1; r[n] = int(round(N * u[i])); rows.append(r)
    L = fmpz_mat(rows).lll()
    rel = []
    for i in range(n):
        z = np.array([int(L[i, j]) for j in range(n)])
        if np.abs(z).max() <= 50 and abs(z @ u) < 1e-6 * np.abs(z).sum():
            rel.append(z)
    if not rel: return 0, []
    Mr = np.array(rel)
    return np.linalg.matrix_rank(Mr), rel
vals = sorted(zeros, key=lambda x: x[0])
print('F values (sorted):', ['%.1e' % v for v, _ in vals], flush=True)
for v, u in vals:
    if v > 1e-20: continue
    u = u / np.abs(u).max()
    rk, rel = relations(u)
    un = u / np.linalg.norm(u); sv = Aarr @ un
    G = 8 * (C.T * (sv ** 7)) @ Aarr / np.linalg.norm(sv ** 8)
    S = np.linalg.svd(G, compute_uv=False)
    jr = int((S > 1e-6 * S[0]).sum())
    print('zero F=%.1e  rational-span dim = %d  Jacobian rank = %d (sv: %s)  u=%s' % (v, 11 - rk, jr,
          np.array2string(S, precision=1, max_line_width=300), np.array2string(u, precision=4, max_line_width=200)), flush=True)
    for z in rel[:12]: print('      relation', z.tolist())
