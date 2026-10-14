"""Is there a semidefinite combination of the 5 sdc11 octics?  LP on samples (evidence)."""
import json, numpy as np
from scipy.optimize import linprog
D = json.load(open('sdc11_m8_basis.json'))
A = np.array(D['A'], float); C = np.array(D['basis'], float); C /= np.abs(C).sum(1, keepdims=True)
rng = np.random.default_rng(5)
U = rng.normal(size=(20000, 11))
P = ((U @ A.T) ** 8) @ C.T
P /= np.linalg.norm(P, axis=1, keepdims=True)
# maximise t s.t. mu.P_j >= t for all j, -1<=mu<=1
n = 5
c = np.zeros(n + 1); c[-1] = -1
Aub = np.hstack([-P, np.ones((len(P), 1))]); bub = np.zeros(len(P))
r = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(-1, 1)] * n + [(None, None)], method='highs')
mu = r.x[:n]; print('LP: max_t min_j mu.P_j/|P_j| = %.4e, mu=%s' % (r.x[-1], np.round(mu, 5)))
# test on fresh samples incl. separated points
for trial in range(3):
    U2 = rng.normal(size=(50000, 11)); P2 = ((U2 @ A.T) ** 8) @ C.T
    v = P2 @ mu
    print('fresh: min mu.P = %.3e, max = %.3e, frac negative %.4f' % (v.min(), v.max(), (v < 0).mean()))
np.save('sdc11_mu.npy', mu)
