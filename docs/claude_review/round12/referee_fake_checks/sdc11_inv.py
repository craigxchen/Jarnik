"""sdc11 (T=1, L=2, base inf): exact Z_11-invariant part of V(A,8) (orbit-reduced kernel, flint), and check of
the cyclic symmetry of the box.  Writes sdc11_inv_basis.json (relative coordinates k', 11 entries)."""
import json, time, itertools, flint, numpy as np
from rlib import *
n = 12; pats = sdc(11); L = 2; d = 11
A, others = box(n, pats, L, 0)
Aset = set(map(tuple, A.tolist()))
shift = lambda k: (k[-1],) + tuple(k[:-1])      # field element i -> i+1 on coordinates 1..11
print('box invariant under the cyclic shift:', all(shift(k) in Aset for k in Aset), flush=True)
# orbits
orb = {}; reps = []
for k in sorted(Aset):
    if k in orb: continue
    o = []; x = k
    for _ in range(d):
        if x not in orb: orb[x] = len(reps); o.append(x)
        x = shift(x)
    reps.append(o)
print('|A| =', len(Aset), ' orbits =', len(reps), flush=True)
# monomial orbit sums of degree <= 7 (reduced exponents 0..2)
monos = reduced_monos(d, 7)
mo = {}; mreps = []
for e in monos:
    if e in mo: continue
    x = e; o = []
    for _ in range(d):
        if x not in mo: mo[x] = len(mreps); o.append(x)
        x = shift(x)
    mreps.append(o)
print('monomial orbits (deg<=7) =', len(mreps), flush=True)
pts = np.array([k for o in reps for k in o], dtype=np.int64)
col = np.array([orb[tuple(k)] for k in pts.tolist()])
rows = []
for o in mreps:
    e = o[0]                       # invariant c annihilates orbit sum iff it annihilates e (c invariant)
    v = np.ones(len(pts), dtype=np.int64)
    for j, ej in enumerate(e):
        if ej: v *= pts[:, j]**ej
    r = np.zeros(len(reps), dtype=np.int64); np.add.at(r, col, v)
    rows.append(r.tolist())
t = time.time()
M = flint.fmpz_mat(rows)
X, nul = M.nullspace()
print('dim of invariant V(A,8) (exact) =', nul, ' (%.0fs)' % (time.time()-t), flush=True)
basis = []
for j in range(nul):
    v = [int(X[i, j]) for i in range(X.nrows())]
    from math import gcd
    g = 0
    for x in v: g = gcd(g, abs(x))
    v = [x//g for x in v]
    c = {}
    for i, o in enumerate(reps):
        for k in o: c[k] = v[i]
    basis.append(c)
# exact verification on ALL reduced monomials of degree <= 7 (not just orbit reps)
Kl = list(Aset); Ka = np.array(Kl, dtype=np.int64)
ok = True
for c in basis:
    cv = np.array([c[k] for k in Kl], dtype=object)
    for e in monos:
        v = np.ones(len(Kl), dtype=np.int64)
        for j, ej in enumerate(e):
            if ej: v *= Ka[:, j]**ej
        if int(np.dot(cv, v.astype(object))) != 0: ok = False; break
    print('  basis vector |c|_1 =', sum(abs(x) for x in c.values()), ' |c|_inf =', max(abs(x) for x in c.values()), flush=True)
print('all invariant basis vectors annihilate all', len(monos), 'reduced monomials of degree <= 7:', ok)
json.dump({'K': Kl, 'C': [[c[k] for k in Kl] for c in basis]}, open('sdc11_inv_basis.json', 'w'))
