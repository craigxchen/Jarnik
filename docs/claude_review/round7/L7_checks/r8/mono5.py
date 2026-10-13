import numpy as np, sys
from mono import *
idx, k = sys.argv[1], sys.argv[2]
V = np.load('m12_type%s_%s.npy' % (idx, k)); E = eval(open('m12_type%s_%s_E.txt' % (idx, k)).read())
rng = np.random.default_rng(11)
X, y = member_to_X(V, E, t0=-0.4+0.2j)
sols = [X]; fails = 0; tot = 0
for L in range(30):
    sc = rng.choice([0.3, 1, 3])
    y1 = y + sc*(rng.normal(size=7) + 1j*rng.normal(size=7)); y2 = y + sc*(rng.normal(size=7) + 1j*rng.normal(size=7))
    for S in list(sols):
        tot += 1
        Z = track(S, y, y1, E)
        if Z is not None: Z = track(Z, y1, y2, E)
        if Z is not None: Z = track(Z, y2, y, E)
        if Z is None: fails += 1; continue
        if not any(np.linalg.norm(Z - T) < 1e-6*(1+np.linalg.norm(T)) for T in sols): sols.append(Z)
print('type', idx, 'member', k, 'sols', len(sols), 'track fails', fails, 'of', tot)
