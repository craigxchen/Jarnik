# kill_search: vertical-line locus q_k = i*y_k, c=0, D=10 (L=4): look for rank-4 clusters (lane claims r<=3).
import numpy as np, itertools
from scipy.optimize import least_squares
D=10; rng=np.random.default_rng(7)
subs = np.array(list(itertools.product([0,1], repeat=D)), float)
found = []
for trial in range(4000):
    Fs = []
    while len(Fs) < 4:
        F = (rng.random(D) < 0.5).astype(float)
        if 5 <= F.sum() <= 8: Fs.append(F)
    Fs = np.array(Fs)
    if np.linalg.matrix_rank(Fs) < 4: continue
    def res(y):
        return np.concatenate([Fs @ y, Fs @ y**3, [np.sum(y**2) - D]])
    y0 = rng.normal(size=D)
    s = least_squares(res, y0, method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    y = s.x
    if np.linalg.norm(res(y)) > 1e-11: continue
    ay = np.sort(np.abs(y))
    if ay[0] < 1e-3 or np.min(np.diff(ay)) < 1e-3: continue
    R = np.sqrt((subs @ y)**2 + (subs @ y**3)**2)
    zs = subs[R < 1e-9]
    rk = np.linalg.matrix_rank(zs)
    found.append((len(zs), rk, ay[0], np.min(np.diff(ay)), y))
print("nondegenerate solutions:", len(found))
from collections import Counter
print(Counter((f[0], f[1]) for f in found))
best = sorted(found, key=lambda f: -min(f[2], f[3]))[:3]
for f in best: print(f[0], f[1], "min|y|", f[2], "min gap |y|", f[3], repr(f[4]))
np.save('vert_sols.npy', np.array([f[4] for f in found]))
