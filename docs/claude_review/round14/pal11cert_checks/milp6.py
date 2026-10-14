import numpy as np, itertools, time
from math import gcd
from scipy.optimize import milp, LinearConstraint, Bounds
dat = np.load('dirs_pal11.npz'); D = dat['D']; tau = dat['tau']
def pool(maxsupp, maxentry, minsupp=2):
    out = []
    for s in range(minsupp, maxsupp + 1):
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
P = pool(6, 1)
NZ = ((P[:, 1:] @ D.T) != 0).astype(float)
t0 = time.time()
res = milp(-np.ones(len(P)), constraints=LinearConstraint(NZ.T, -np.inf, tau.astype(float)), integrality=np.ones(len(P)),
           bounds=Bounds(0, 8), options={'time_limit': 1500})
print(res.status, res.message, -res.fun if res.x is not None else None, '%.0fs' % (time.time() - t0))
if res.x is not None:
    for i in np.nonzero(res.x > 0.5)[0]:
        print(int(round(res.x[i])), P[i].tolist())
