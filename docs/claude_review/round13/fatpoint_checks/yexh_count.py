import itertools, numpy as np, sys
import yblock as yb
M, Y = yb.build_Y(3, 2)
Zs = list(range(8)); heavy = []; seen = set()
for h in range(1, 8):
    cos = sorted({tuple(sorted({x ^ y for y in [0, h]})) for x in Zs})
    for C1, C2 in itertools.combinations(cos, 2):
        c = np.zeros(8, dtype=np.int64)
        for x in C1: c[x] += 1
        for x in C2: c[x] -= 1
        a = tuple(M @ c)
        if a not in seen and tuple(-x for x in a) not in seen: seen.add(a); heavy.append(np.array(a))
lines = [(a, c, a ^ c) for a in range(1, 8) for c in range(a + 1, 8) if (a ^ c) > c]
nz = [np.nonzero(g)[0] for g in heavy]
tl = [[li for li, L in enumerate(lines) if any(heavy[i][q - 1] for q in L)] for i in range(len(heavy))]
for e, nmax in [(7, 5), (7, 6)]:
    cnt = np.zeros(7, dtype=np.int64); lc = np.zeros(7, dtype=np.int64); tot = [0]
    def rec(s, k):
        if k == e: tot[0] += 1; return
        for i in range(s, len(heavy)):
            cnt[nz[i]] += 1; lc[tl[i]] += 1
            if (cnt <= 4).all() and (lc <= nmax).all(): rec(i, k + 1)
            cnt[nz[i]] -= 1; lc[tl[i]] -= 1
    rec(0, 0); print(e, nmax, tot[0], flush=True)
