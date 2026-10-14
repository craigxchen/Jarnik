"""Syl(16,2): find ALL primitive v in {-1,0,1}^15 whose width power has degree <= 6 and is nonzero in R
(|<v,A>| <= 6, support >= 3).  A = -A and 0 in A, so |<v,A>| <= 6 iff max <v,A> <= 2.
Filter: points of A with support <= 3 (exact subset of A); survivors are checked on all of A."""
import itertools, numpy as np, time, models
from syl16 import syl16_singer
t0 = time.time()
d, cols, Ls, colA = syl16_singer(2)
A = models.box(d, cols, Ls).astype(np.int8)
small = A[(A != 0).sum(axis=1) <= 3].astype(np.int16)
print("|A| =", len(A), " small points:", len(small), flush=True)
surv = []
V = np.array(list(itertools.product((-1, 0, 1), repeat=15)), dtype=np.int8)
nz = (V != 0).sum(axis=1)
first = np.array([v[np.nonzero(v)[0][0]] if n else 0 for v, n in zip(V, nz)])
V = V[(nz >= 3) & (first > 0)]
print("candidate directions (support >= 3, first nonzero > 0):", len(V), flush=True)
for s in range(0, len(V), 20000):
    blk = V[s:s + 20000].astype(np.int16)
    mx = (small @ blk.T).max(axis=0)
    surv.extend(blk[mx <= 2].tolist())
print("survivors after small-point filter:", len(surv), f"[{time.time()-t0:.0f}s]", flush=True)
res = []
Ai = A.astype(np.int16)
for s in range(0, len(surv), 2000):
    blk = np.array(surv[s:s + 2000], dtype=np.int16)
    mx = (Ai @ blk.T).max(axis=0)
    for v, m_ in zip(blk, mx):
        if m_ <= 2: res.append(tuple(int(z) for z in v))
print("directions with |<v,A>| <= 5 and support >= 3:", len(res))
cs = {tuple(c) for c in cols}
print("  of which columns:", sum(1 for v in res if v in cs), "; others:", [v for v in res if v not in cs][:10])
