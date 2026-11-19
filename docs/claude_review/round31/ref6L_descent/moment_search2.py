"""Vectorized version of moment_search.py: T real, c_e=x_e+i y_e; find strict (x,y) with
Im star power sums equal for k=1,2 (then test 3,4,5 exactly). x normalized by x_0 = 0."""
import itertools, sys
import numpy as np
from itertools import combinations
edges = list(combinations(range(5), 2))
Inc = np.array([[int(i in e) for e in edges] for i in range(5)])
E = Inc[:4] - Inc[4]
Ymax = int(sys.argv[1]); Xmax = int(sys.argv[2])
rng = [v for v in range(-Ymax, Ymax+1) if v != 0]
def pivots(Mat, exclude=()):
    for cols in combinations([c for c in range(10) if c not in exclude], 4):
        if abs(np.linalg.det(Mat[:, cols].astype(float))) > 1e-9:
            return list(cols)
piv = pivots(E); free = [c for c in range(10) if c not in piv]
Einv = np.linalg.inv(E[:, piv].astype(float))
YF = np.array(list(itertools.product(rng, repeat=6)))
YP = -(Einv @ (E[:, free] @ YF.T)).T
YPr = np.rint(YP)
mask = (np.abs(YP-YPr).max(1) < 1e-9) & (YPr != 0).all(1) & (np.abs(YPr).max(1) <= Ymax)
Ys = np.zeros((mask.sum(), 10), dtype=np.int64); Ys[:, free] = YF[mask]; Ys[:, piv] = YPr[mask]
# remove global sign duplicates
Ys = np.array(sorted({tuple(y) if tuple(y) > tuple(-y) else tuple(-y) for y in Ys}))
print('y count', len(Ys))
def impow(x, y, k):
    z = (x.astype(object) + 1j*0)
    out = []
    for a, b in zip(x, y):
        a = int(a); b = int(b); s = 0
        from math import comb
        for j in range(1, k+1, 2):
            s += comb(k, j) * a**(k-j) * b**j * (-1)**((j-1)//2)
        out.append(s)
    return np.array(out, dtype=object)
def ok(x, y, k):
    v = impow(x, y, k); st = [sum(v[e] for e in range(10) if i in edges[e]) for i in range(5)]
    return len(set(st)) == 1
cnt = {2: 0, 3: 0, 4: 0, 5: 0}; ex = {}
xr = np.arange(-Xmax, Xmax+1)
for y in Ys:
    M = E * y[None, :]
    p = pivots(M, exclude=(0,))
    if p is None: continue
    f = [c for c in range(10) if c not in p]  # includes 0
    f_nz = [c for c in f if c != 0]
    Minv = np.linalg.inv(M[:, p].astype(float))
    XF = np.array(list(itertools.product(xr, repeat=len(f_nz))))
    XP = -(Minv @ (M[:, f_nz] @ XF.T)).T
    XPr = np.rint(XP)
    mk = (np.abs(XP-XPr).max(1) < 1e-9) & (np.abs(XPr).max(1) <= Xmax)
    for xf, xp in zip(XF[mk], XPr[mk]):
        x = np.zeros(10, dtype=np.int64); x[f_nz] = xf; x[p] = xp
        c = set(zip(x, y)) | set(zip(x, -y))
        if len(c) != 20: continue
        cnt[2] += 1
        for k in (3, 4, 5):
            if ok(x, y, k):
                cnt[k] += 1; ex.setdefault(k, (x.tolist(), y.tolist()))
            else:
                break
print('counts', cnt)
for k, v in ex.items(): print('example k>=%d' % k, v)
