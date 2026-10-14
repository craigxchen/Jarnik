"""LP relaxation with width counts AND 2-direction family counts (fatpoint Prop 2.6(b)): for F = (v_i, v_j),
#{gens touching v_i or v_j} <= rho(B_F).  usage: python3 lpfam.py MODE (as mkpool2)"""
import sys, numpy as np, itertools
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
fam = np.load('fam2_pal11.npy')
keep = fam[:, 2] < np.minimum(tau[fam[:, 0]] + tau[fam[:, 1]], 8)
fam = fam[keep]
print('informative families', len(fam))
def readpool(fn):
    Z = []
    with open(fn) as f:
        f.readline(); f.readline()
        for line in f:
            t = line.split(); Z.append([int(x) for x in t[1:13]])
    return np.array(Z, dtype=np.int64)
P = readpool(sys.argv[1])
NZ = (P[:, 1:] @ D.T) != 0
NF = NZ[:, fam[:, 0]] | NZ[:, fam[:, 1]]
Aub = vstack([csr_matrix(NZ.T, dtype=np.float64), csr_matrix(NF.T, dtype=np.float64)]).tocsr()
del NF
bub = np.concatenate([tau, fam[:, 2]]).astype(float)
res = linprog(-np.ones(len(P)), A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs')
print(sys.argv[1], len(P), 'LP with families:', -res.fun)
np.save(sys.argv[1] + '.lpx.npy', res.x)
