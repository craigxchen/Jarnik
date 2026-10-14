import numpy as np, itertools, flint
from tclib import box_points, cube_monomials, mon_eval
def kernel(pts, m, d):
    mons = [a for a in cube_monomials(d) if sum(a) < m]
    M = flint.fmpz_mat(len(mons), len(pts), [mon_eval(k, a) for a in mons for k in pts])
    X, nullity = M.nullspace()
    return [[int(X[i, j]) for i in range(len(pts))] for j in range(nullity)]
q = 11; D = [0,1,3,4,5,9]; L = 2
cols = [sorted((i + x) % q for x in D) for i in range(q)]
rng = np.random.default_rng(0)
for Y in [(0,1,2,3,4,5,6), (0,1,2,3,4,5,7), (0,1,2,3,4,5,9)]:
    idx = {y: i for i, y in enumerate(Y)}
    rc = [[idx[x] for x in c if x in idx] for c in cols]; rc = [c for c in rc if c]
    pts = box_points(7, [(c, L) for c in rc]); K = np.array(pts, dtype=float)
    c = np.array(kernel(pts, 8, 7)[0], dtype=float)
    U = rng.standard_normal((200000, 7))
    vals = ((U @ K.T) ** 8) @ c
    print(Y, 'min %.3e max %.3e  frac positive %.3f' % (vals.min(), vals.max(), (vals > 0).mean()))
