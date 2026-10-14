import numpy as np, itertools
from math import gcd
from scipy.optimize import linprog
from fractions import Fraction
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
exec(open('poolgen.py').read())
P = pool(5, 1)
NZ = ((P[:, 1:] @ D.T) != 0).astype(float)
# dual: min sum tau_j y_j s.t. NZ y >= 1, y >= 0
res = linprog(tau.astype(float), A_ub=-NZ, b_ub=-np.ones(len(P)), bounds=(0, None), method='highs')
print('dual value', res.fun)
y = res.x
sup = np.nonzero(y > 1e-9)[0]
print('support size', len(sup))
for j in sup:
    print('  y=%.6f tau=%d v=%s' % (y[j], tau[j], D[j].tolist()))
