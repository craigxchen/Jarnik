# Image-size scaling test for the dimension-4 variant over F_{2^k} (k=5,6), numpy, log tables.
import itertools, sys
import numpy as np
from gf import GF2k
k = int(sys.argv[1]); F = GF2k(k); q = F.q; N = q - 1
g = next(c for c in range(2, q) if len({F.pow(c, i) for i in range(N)}) == N)
exp = [F.pow(g, i) for i in range(N)]; logt = np.zeros(q, dtype=np.int64)
for i, e in enumerate(exp): logt[e] = i
sel = [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,1)] if (len(sys.argv) < 3 or sys.argv[2] == "a") else \
      [(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1)]
# coefficient vectors c in F_2^4 on (1,x2,x3,x4) with l(c)=0, c != 0 for each H
Hc = []
for lv in sel:
    Hc.append([c for c in itertools.product([0,1], repeat=4) if any(c) and
               (c[0] ^ (c[1] & lv[0]) ^ (c[2] & lv[1]) ^ (c[3] & lv[2])) == 0])
x3, x4, Y = np.meshgrid(np.arange(q), np.arange(q), np.arange(q), indexing="ij")
x3 = x3.ravel(); x4 = x4.ravel(); Y = Y.ravel()
keys = []; dom = 0
for x2 in range(q):
    # all 16 elements of V; independence: all nonzero combos nonzero; Y not in V
    ok = np.ones(x3.shape, dtype=bool)
    for c in itertools.product([0,1], repeat=4):
        v = (c[0] * 1) ^ (c[1] * x2) ^ (c[2] * x3) ^ (c[3] * x4)
        if any(c): ok &= (v != 0)
        ok &= (Y != v)
    if not ok.any(): continue
    a3, a4, aY = x3[ok], x4[ok], Y[ok]; dom += len(aY)
    logs = []
    for H in Hc:
        s = np.zeros(aY.shape, dtype=np.int64)
        for c in H:
            v = (c[0] * 1) ^ (c[1] * x2) ^ (c[2] * a3) ^ (c[3] * a4)
            s += logt[aY ^ v]
        logs.append(s % N)
    key = np.zeros(aY.shape, dtype=np.int64)
    for s in logs[1:]:
        key = key * N + ((s - logs[0]) % N)
    keys.append(key)
keys = np.concatenate(keys)
u = np.unique(keys)
print(f"F_{q} sel={sys.argv[2] if len(sys.argv)>2 else 'a'}: domain {dom}, distinct images {len(u)},"
      f" ratio images/q^3 = {len(u)/q**3:.2f}, images/q^4 = {len(u)/q**4:.3f}")
