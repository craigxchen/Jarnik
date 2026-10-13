import sys, json, math
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from rcore import gauss, Lprime, tight_list
from farkas import load
ps = load(sys.argv[1]); M = int(sys.argv[2])
gp = [gauss(p) for p in ps]; r = len(gp); L = Lprime(2 * M, set(ps))
rels = tight_list(gp, L)
n = len(rels); A = lil_matrix((n, r + 1)); b = np.zeros(n)
for i, (v, s, sg, s2) in enumerate(rels):
    for col, co in v: A[i, col] += -sg * co
    A[i, r] = 1.0
    b[i] = -(2 * math.asin(math.sqrt(float(s2))) - math.pi / 2)
c = np.zeros(r + 1); c[r] = -1
res = linprog(c, A_ub=A.tocsr(), b_ub=b, bounds=[(None, None)] * r + [(None, 1)], method='highs')
print(res.status, res.x[r])
du = res.ineqlin.marginals
for i in np.argsort(du)[:12]:
    if du[i] < -1e-9:
        v, s, sg, s2 = rels[i]
        print(round(du[i], 4), v, s, sg, round(math.asin(math.sqrt(float(s2))), 4))
