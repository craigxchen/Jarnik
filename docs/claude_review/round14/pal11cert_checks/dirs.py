"""Exact widths tau_A(v) = #{<v,k'> : k' in A} - 1 for test directions v; keep those with tau < m = 8.
Directions: all of {-1,0,1}^11 and all of {-2..2}^11 with support <= 4 (projectively normalised: first nonzero > 0,
primitive).  Output: dirs_pal11.npz (D, tau)."""
import itertools, time
from math import gcd
import numpy as np
from fplib import pal_patterns, box
M, pats = pal_patterns(11); A = box(M, pats, 2); d = 11
Aarr = np.array(A, dtype=np.int64)
t0 = time.time()
dirs = set()
for v in itertools.product((-1, 0, 1), repeat=d):
    if any(v):
        f = next(x for x in v if x)
        if f > 0: dirs.add(v)
for s in range(1, 5):
    for S in itertools.combinations(range(d), s):
        for vals in itertools.product((-2, -1, 1, 2), repeat=s):
            if vals[0] < 0: continue
            g = 0
            for x in vals: g = gcd(g, abs(x))
            if g != 1: continue
            v = [0]*d
            for x, a in zip(S, vals): v[x] = a
            dirs.add(tuple(v))
D = np.array(sorted(dirs), dtype=np.int64)
tau = np.zeros(len(D), dtype=np.int64)
for s0 in range(0, len(D), 4000):
    vals = Aarr @ D[s0:s0+4000].T
    vs = np.sort(vals, axis=0)
    tau[s0:s0+4000] = (np.diff(vs, axis=0) != 0).sum(axis=0)
print('directions tested', len(D), 'time %.1f' % (time.time()-t0))
for t in range(0, 12):
    print(' tau =', t, ':', int((tau == t).sum()))
keep = tau < 8
np.savez('dirs_pal11.npz', D=D[keep], tau=tau[keep])
print('kept (tau < 8):', int(keep.sum()))
# base direction (1,...,1): <v,k'> = -k_0
one = np.ones(d, dtype=np.int64)
print('tau(all-ones) =', len(set((Aarr @ one).tolist())) - 1)
