"""LP relaxation of the width-count design problem: max sum_i x_i s.t. sum_i x_i [<z_i',v_j> != 0] <= tau_j, x >= 0.
If the LP optimum is < m, no multiset of m pool generators satisfies the width counts (rigorous via the dual)."""
import sys
import numpy as np
from scipy.optimize import linprog
from math import gcd
import itertools
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
def pool(maxsupp, maxentry):
    out = []
    for s in range(2, maxsupp + 1):
        for S in itertools.combinations(range(12), s):
            for vals in itertools.product([e for e in range(-maxentry, maxentry + 1) if e], repeat=s):
                if sum(vals) != 0 or vals[0] < 0: continue
                g = 0
                for v in vals: g = gcd(g, abs(v))
                if g != 1: continue
                z = [0]*12
                for x, v in zip(S, vals): z[x] = v
                out.append(z)
    return np.array(out, dtype=np.int64)
for ms, me in [(2, 1), (4, 1), (3, 2), (4, 2), (5, 1), (6, 1)]:
    P = pool(ms, me)
    NZ = ((P[:, 1:] @ D.T) != 0).astype(float)
    res = linprog(-np.ones(len(P)), A_ub=NZ.T, b_ub=tau.astype(float), bounds=(0, None), method='highs')
    print(f'pool supp<={ms} entries<={me}: {len(P)} generators; LP max sum x = {-res.fun:.6f}', flush=True)
