import numpy as np, itertools
from math import gcd
from scipy.optimize import linprog
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
exec(open('poolgen.py').read())
P = pool(6, 1)
NZ = ((P[:, 1:] @ D.T) != 0)
core = [j for j in range(len(D)) if tau[j] == 4 and sorted(set(D[j].tolist())) == [0, 1] and D[j].sum() == 6]
print('core dirs found', len(core))
tc = NZ[:, core].sum(axis=1)
supp = (P != 0).sum(axis=1)
for s in (2, 4, 6):
    print('support', s, 'min #core touched', tc[supp == s].min(), ' #gens touching <=5 core:', int(((supp == s) & (tc <= 5)).sum()))
res = linprog(-np.ones(len(P)), A_ub=NZ.T.astype(float), b_ub=tau.astype(float), bounds=(0, None), method='highs')
print('LP', -res.fun)
for i in np.nonzero(res.x > 1e-9)[0]:
    print('  x=%.4f z=%s supp=%d coretouch=%d ntouch=%d' % (res.x[i], P[i].tolist(), supp[i], tc[i], NZ[i].sum()))
