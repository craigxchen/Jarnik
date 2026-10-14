"""LP relaxation with width counts and 2-direction family counts (fatpoint Prop 2.6(b)), by constraint generation:
start from the width counts, add the violated family constraints #{gens touching v_i or v_j} <= rho(B_F), resolve.
The final LP value is the value of the full LP (all 28444 informative family constraints are checked at the end).
usage: python3 lpfam2.py pool.txt"""
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack

dat = np.load('dirs_pal11.npz')
D = dat['D']
tau = dat['tau']
fam = np.load('fam2_pal11.npy')
keep = fam[:, 2] < np.minimum(tau[fam[:, 0]] + tau[fam[:, 1]], 8)
fam = fam[keep]


def readpool(fn):
    Z = []
    with open(fn) as f:
        f.readline()
        f.readline()
        for line in f:
            t = line.split()
            Z.append([int(x) for x in t[1:13]])
    return np.array(Z, dtype=np.int64)


P = readpool(sys.argv[1])
NZ = (P[:, 1:] @ D.T) != 0             # pool x dirs (bool)
A0 = csr_matrix(NZ.T.astype(np.float64))
active = []
for it in range(100):
    if active:
        fa = fam[active]
        NF = NZ[:, fa[:, 0]] | NZ[:, fa[:, 1]]
        Aub = vstack([A0, csr_matrix(NF.T.astype(np.float64))]).tocsr()
        bub = np.concatenate([tau, fa[:, 2]]).astype(float)
    else:
        Aub, bub = A0, tau.astype(float)
    res = linprog(-np.ones(len(P)), A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs')
    x = res.x
    # loads of all families
    load = np.zeros(len(fam))
    for s0 in range(0, len(fam), 2000):
        f = fam[s0:s0 + 2000]
        NFc = NZ[:, f[:, 0]] | NZ[:, f[:, 1]]
        load[s0:s0 + 2000] = x @ NFc
    viol = np.nonzero(load > fam[:, 2] + 1e-7)[0]
    print(f'iter {it}: LP value {-res.fun:.6f}, active families {len(active)}, violated {len(viol)}', flush=True)
    if len(viol) == 0:
        break
    order = viol[np.argsort(-(load[viol] - fam[viol, 2]))][:500]
    active = sorted(set(active) | set(order.tolist()))
print(sys.argv[1], len(P), 'final LP value with all family constraints:', -res.fun)
np.save(sys.argv[1] + '.lpfam_x.npy', x)
