# validation of the numeric pipeline: (i) one section octic alone (known to have nondegenerate real zeros);
# (ii) ten random sections together; (iii) all 30 classes x 11 shifts (the full run is paley_numeric.py)
import numpy as np, itertools, sys
from scipy.optimize import least_squares
exec(open('paley_numeric.py').read().split("rng = np.random.default_rng")[0])
rng = np.random.default_rng(5)
def run(secs, label, ntr=30):
    def res(v):
        u = v/np.linalg.norm(v)
        return np.array([((K @ u[orig])**8) @ c for orig, K, c in secs])
    out = []
    for tr in range(ntr):
        sol = least_squares(res, rng.standard_normal(q), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        u = sol.x/np.linalg.norm(sol.x)
        # degeneracy restricted to the coordinates involved
        out.append((np.linalg.norm(res(u)), degeneracy(u)))
    out.sort()
    good = [o for o in out if o[0] < 1e-12]
    print(label, 'zeros found (resid<1e-12):', len(good), ' max degeneracy among them: %.3e' % (max([g[1] for g in good]) if good else -1))
run(sections[:1], 'one section')
run([sections[i] for i in rng.choice(len(sections), 3, replace=False)], 'three sections')
run([sections[i] for i in rng.choice(len(sections), 10, replace=False)], 'ten sections')
