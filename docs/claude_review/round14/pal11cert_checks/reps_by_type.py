"""Universal representative pool U for the width-count problem of pal11 (L=2, m=8).

Input: the exact flat lists typeflats_2_3_4_5_6_7_8.pkl (supports <= 8, all t) and excflats.pkl (support 12, t <= 5).
Allowed (sigma,t) types: those occurring in some feasible profile (profiles2.py):
  (2,6) (3,9) (4,6) (4,8) (4,9) (4,10) (6,2) (6,4) (6,6) (7,5) (7,6) (8,4) (8,5) (8,6) (12,1) (12,2).
For every flat of an allowed type with complement Pi (basis B, dim delta):
  delta = 1: the vector itself;
  delta = 2: for each test direction v with g_v = (<b_i', v>)_i != 0, the line of Pi orthogonal to g_v;
  delta = 3: for each pair of independent g_v, g_w, the line g_v x g_w.
Every nonzero point of Pi lies in a minimal flat of the arrangement induced on Pi by the 1519 test directions; minimal
flats are lines (a non-base coordinate direction always cuts a flat of dim >= 2 properly), and a point of the
minimal flat inside the flat of z touches a subset of the directions touched by z (dominance).  Representatives are
made primitive, sign-normalised, deduplicated, and kept only if their own (sigma,t) type is allowed (a dominated
count-feasible multiset is count-feasible, so its members have allowed types).
Output: pool file in the csearch format (mkpool2 conventions) and univpool_stats.txt."""
import pickle, sys, itertools
from math import gcd
from collections import Counter
import numpy as np
from fplib import pal_patterns

dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
Cv = np.zeros((11, 11), dtype=np.int64)
for c, P in enumerate(core):
    for x in P: Cv[c, x - 1] = 1
ALLOWED_ALL = {(2, 6), (3, 9), (4, 6), (4, 8), (4, 9), (4, 10), (6, 2), (6, 4), (6, 6), (7, 5), (7, 6), (8, 4), (8, 5),
           (8, 6), (12, 1), (12, 2)}
ALLOWED = set(eval(sys.argv[1]))
flats = pickle.load(open('typeflats_2_3_4_5_6_7_8.pkl', 'rb')) + \
    [f for f in pickle.load(open('excflats.pkl', 'rb')) if f[0] == 12]
def norm(z):
    g = 0
    for x in z: g = gcd(g, abs(int(x)))
    z = [int(x) // g for x in z]
    f = next(x for x in z if x)
    return tuple(z) if f > 0 else tuple(-x for x in z)
def typ(z):
    zz = np.array(z, dtype=np.int64)
    s = int((zz != 0).sum())
    t = int(((Cv @ zz[1:]) != 0).sum())
    return (s, t)
U = set()
stat = Counter()
for (s, t, R, Z, dim, basis) in flats:
    if (s, t) not in ALLOWED: continue
    stat[(s, t, dim)] += 1
    Bf = np.zeros((dim, 12), dtype=np.int64)
    for i, b in enumerate(basis):
        for x, val in zip(R, b): Bf[i, x] = val
    assert (Bf.sum(axis=1) == 0).all()
    if dim == 1:
        U.add(norm(Bf[0])); continue
    G = Bf[:, 1:] @ D.T          # dim x ND
    gs = set()
    for j in range(G.shape[1]):
        g = G[:, j]
        if g.any(): gs.add(norm(g))
    gs = [np.array(g, dtype=np.int64) for g in gs]
    if dim == 2:
        for g in gs:
            lam = np.array([g[1], -g[0]])
            z = lam @ Bf
            if z.any(): U.add(norm(z))
    elif dim == 3:
        for g, h in itertools.combinations(gs, 2):
            lam = np.cross(g, h)
            if lam.any():
                z = lam @ Bf
                if z.any(): U.add(norm(z))
    else:
        raise ValueError('unexpected dimension %d' % dim)
print('flats used by (sigma,t,dim):', sorted(stat.items()), flush=True)
U = sorted(U)
types = Counter(typ(z) for z in U)
keep = [z for z in U if typ(z) in ALLOWED_ALL]
print('representatives', len(U), ' kept (allowed type)', len(keep))
print('types of kept representatives:', sorted(Counter(typ(z) for z in keep).items()))
pickle.dump(keep, open(sys.argv[2], 'wb'))
