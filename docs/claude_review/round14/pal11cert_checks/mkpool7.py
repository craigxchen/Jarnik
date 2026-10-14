"""Pool P7 for the universal search: ALL nonzero balanced z (up to scale) of the types (2,6), (4,6), (6,2), (6,4),
(12,1).  Each of these types consists of flats with 1-dimensional complement (typeflats/excflats), so the explicit
vectors are all such generators.  Written in csearch format, ordered by t (exceptional first)."""
import pickle
from math import gcd
import numpy as np
from fplib import pal_patterns
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
TYPES = {(2, 6), (4, 6), (6, 2), (6, 4), (12, 1)}
fl = pickle.load(open('typeflats_2_3_4_5_6_7_8.pkl', 'rb')) + pickle.load(open('excflats.pkl', 'rb'))
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
Cv = np.zeros((11, 11), dtype=np.int64)
for c, P in enumerate(core):
    for x in P: Cv[c, x - 1] = 1
def norm(z):
    g = 0
    for x in z: g = gcd(g, abs(int(x)))
    z = [int(x) // g for x in z]
    f = next(x for x in z if x)
    return tuple(z) if f > 0 else tuple(-x for x in z)
S = set()
for (s, t, R, Z, dim, basis) in fl:
    if (s, t) in TYPES:
        assert dim == 1
        z = [0] * 12
        for x, v in zip(R, basis[0]): z[x] = v
        S.add(norm(z))
def tc(z): return int(((Cv @ np.array(z[1:])) != 0).sum())
pool = sorted(S, key=lambda z: (tc(z), sum(1 for x in z if x), z))
from collections import Counter
print('P7 size', len(pool), Counter((sum(1 for x in z if x), tc(z)) for z in pool))
print('max |entry|', max(max(abs(x) for x in z) for z in pool))
index = {z: i for i, z in enumerate(pool)}
n = 12
def act(a, z):
    w = [0]*n; w[0] = z[0]
    for i in range(11): w[1 + (a*i) % 11] = z[1+i]
    return norm(w)
rep = [1 if i == min(index[act(a, z)] for a in (1, 3, 4, 5, 9)) else 0 for i, z in enumerate(pool)]
P = np.array(pool, dtype=np.int64)[:, 1:]
NZ = (P @ D.T) != 0
with open('pool_P7.txt', 'w') as f:
    f.write(f'{len(pool)} {len(D)}\n')
    f.write(' '.join(str(int(t)) for t in tau) + '\n')
    for i, z in enumerate(pool):
        idx = np.nonzero(NZ[i])[0]
        f.write(f'{rep[i]} ' + ' '.join(map(str, z)) + f' {len(idx)} ' + ' '.join(map(str, idx)) + '\n')
print('orbit reps', sum(rep))
