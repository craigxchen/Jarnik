"""sdc11 model (n = 12, T = 1, L = 2): exact V(A,8) via the Z_11-invariant subspace, and the 5 tangent octics.

dim V(A,8) = H(8) = 5 mod p (cyc_sdc11_L2_out.txt, two primes).  Over Q, the Z_11-isotypic components of V(A,8)
for the 10 nontrivial characters are Galois conjugate and have equal dimension h_1, so 5 = h_0 + 10 h_1 forces
h_1 = 0 (if the Q-dimension is 5).  Here we compute the invariant part exactly: c constant on shift-orbits,
one linear condition per monomial orbit of degree < 8.  If its dimension is 5 = dim_p V(A,8) >= dim_Q V(A,8),
then V(A,8) is exactly this space.
Then: numerical search (evidence) for common real zeros maximising the separation, saved for certification.
"""
import sys, json, time, itertools
import numpy as np
import flint
from scipy.optimize import minimize
from flib import *
from cyc import orbit_reps, shift

M, pats = sdc_patterns(11)
q = 11; d = 11; m = 8
A = box(M, pats, 2)
Aset = {k: i for i, k in enumerate(A)}
t0 = time.time()
preps = orbit_reps(A, q)
orbit_of = {}
orbits = []
for oi, (k, sz) in enumerate(preps):
    orb = sorted({shift(k, g, q) for g in range(q)})
    orbits.append(orb)
    for kk in orb:
        orbit_of[kk] = oi
mons = []
for j in range(m):
    mons += reduced_monomials(d, j)
mreps = orbit_reps(mons, q)
print('|A|=%d point orbits=%d monomial orbits (deg<8)=%d (%.1fs)' % (len(A), len(preps), len(mreps), time.time() - t0), flush=True)
Aa = np.array(A, dtype=np.int64)
rows = []
for a, _ in mreps:
    v = np.ones(len(A), dtype=np.int64)
    for x, e in enumerate(a):
        if e:
            v = v * Aa[:, x] ** e
    row = [0] * len(preps)
    for i, k in enumerate(A):
        row[orbit_of[k]] += int(v[i])
    rows.append(row)
Mt = flint.fmpz_mat(rows)
X, nul = Mt.nullspace()
print('invariant V(A,8): dim = %d (%.1fs)' % (nul, time.time() - t0), flush=True)
B = flint.fmpz_mat([[X[i, j] for i in range(len(preps))] for j in range(nul)]).lll()
basis_orb = [[int(B[i, j]) for j in range(len(preps))] for i in range(B.nrows())]
basis = [[bo[orbit_of[k]] for k in A] for bo in basis_orb]
# exact verification against ALL monomials of degree < 8 (not only orbit reps)
Ev = None
cnt = 0
for a in mons:
    v = np.ones(len(A), dtype=np.int64)
    for x, e in enumerate(a):
        if e:
            v = v * Aa[:, x] ** e
    for c in basis:
        assert int(np.dot(np.array(c, dtype=object), v.astype(object))) == 0
    cnt += 1
print('verified exactly: all %d basis vectors annihilate all %d monomials of degree < 8' % (len(basis), cnt))
idx = {k: i for i, k in enumerate(A)}
for c in basis:
    ev = all(c[idx[tuple(-x for x in k)]] == c[i] for i, k in enumerate(A))
    print('  basis vector: |c|_inf=%d |c|_1=%d support=%d even=%s' % (max(abs(x) for x in c), sum(abs(x) for x in c), sum(1 for x in c if x), ev))
json.dump(dict(A=A, basis=basis), open('sdc11_m8_basis.json', 'w'))
