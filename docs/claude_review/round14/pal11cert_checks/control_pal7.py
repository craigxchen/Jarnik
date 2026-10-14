"""Positive control for the count-search programs: pal7 (n = 8, L = 2), m = 7, pool = balanced z with support <= 4 and
entries +-1, test directions = all of {-2..2}^7 (as in fatpoint prodsearch.py for d <= 8).  fatpoint §5.2 reports 57
count-feasible multisets (all of them zonotope certificates).  We write the files for csearch4 / csearch9 with NO
symmetry reduction (every generator flagged as representative) and weight vectors adapted to pal7 (core columns,
coordinates + all-ones, all directions)."""
import itertools
from math import gcd
import numpy as np
from fplib import pal_patterns, box
M, pats = pal_patterns(7); A = box(M, pats, 2); d = M - 1
print('pal7 |A| =', len(A))
Aarr = np.array(A, dtype=np.int64)
dirs = []
for v in itertools.product(range(-2, 3), repeat=d):
    if any(v):
        f = next(x for x in v if x)
        if f > 0: dirs.append(v)
D = np.array(dirs, dtype=np.int64)
vals = Aarr @ D.T
vs = np.sort(vals, axis=0)
tau = (np.diff(vs, axis=0) != 0).sum(axis=0)
keep = tau < 7
D = D[keep]; tau = tau[keep]
print('test dirs with tau < 7:', len(D))
pool = []
for s in range(2, 5):
    for S in itertools.combinations(range(M), s):
        for vv in itertools.product([-1, 1], repeat=s):
            if sum(vv) != 0 or vv[0] < 0: continue
            z = [0] * M
            for x, a in zip(S, vv): z[x] = a
            pool.append(tuple(z))
P = np.array(pool, dtype=np.int64)[:, 1:]
NZ = (P @ D.T) != 0
with open('ctl_pool.txt', 'w') as f:
    f.write(f'{len(pool)} {len(D)}\n' + ' '.join(map(str, tau.tolist())) + '\n')
    for i, z in enumerate(pool):
        idx = np.nonzero(NZ[i])[0]
        # csearch reads exactly 12 z-entries: pad with zeros (unused by the search logic)
        zz = list(z) + [0] * (12 - len(z))
        f.write('1 ' + ' '.join(map(str, zz)) + f' {len(idx)} ' + ' '.join(map(str, idx)) + '\n')
with open('ctl_dirs.txt', 'w') as f:
    f.write(f'{len(D)}\n')
    for v in D: f.write(' '.join(map(str, list(v) + [0] * (11 - d))) + '\n')
core = []
for P_ in pats:
    v = [0] * d
    for x in P_: v[x - 1] = 1
    core.append(tuple(v))
dirindex = {tuple(v): j for j, v in enumerate(D.tolist())}
Y = []
y = [0] * len(D)
for v in core:
    if v in dirindex: y[dirindex[v]] = 1
Y.append(y)
y = [1 if (np.count_nonzero(D[j]) == 1 and tau[j] == 2) or tuple(D[j]) == tuple([1] * d) else 0 for j in range(len(D))]
Y.append(y)
Y.append([1] * len(D))
with open('ctl_weights.txt', 'w') as f:
    f.write(f'{len(Y)}\n')
    for y in Y: f.write(' '.join(map(str, y)) + '\n')
open('ctl_pool.txt.cls', 'w').write(' '.join(['0'] * len(pool)) + '\n')
print('pool', len(pool))
