"""Restricted DFS for zonotopes in the sdc11 / pal11 T=1 L=2 boxes (base inf): generators with |v|_1 <= 2,
entries in {-1,0,1}, all lengths 1.  Finds a zonotope of total length >= target if one of this type exists."""
import sys, itertools, numpy as np, time
from rlib import *
name, target = sys.argv[1], int(sys.argv[2])
n = 12; d = 11
pats = sdc(11) if name == 'sdc11' else pal(11)
A, _ = box(n, pats, 2, 0)
pts = [tuple(x) for x in A.tolist()]; idx = {p: i for i, p in enumerate(pts)}; N = len(pts); Aarr = np.array(pts)
cands = []
for i in range(d):
    v = [0]*d; v[i] = 1; cands.append(tuple(v))
for i, j in itertools.combinations(range(d), 2):
    for s in (1, -1):
        v = [0]*d; v[i] = 1; v[j] = s; cands.append(tuple(v))
shiftidx = {}
for v in cands:
    arr = np.full(N, -1)
    for a in range(N):
        arr[a] = idx.get(tuple(Aarr[a] + np.array(v)), -1)
    shiftidx[v] = arr
print(name, '|A| =', N, ' candidates', len(cands), flush=True)
t0 = time.time(); best = [0, None]
def dfs(mask, gens, start):
    if len(gens) > best[0]:
        best[0] = len(gens); best[1] = list(gens)
        print('  length %d: %s (%.0fs)' % (len(gens), gens, time.time()-t0), flush=True)
        if len(gens) >= target: raise StopIteration
    for ci in range(start, len(cands)):
        v = cands[ci]; sh = shiftidx[v]
        ok = sh >= 0
        m1 = np.zeros_like(mask); m1[ok] = mask[ok] & mask[sh[ok]]
        if m1.any():
            gens.append(v); dfs(m1, gens, ci + 1); gens.pop()
try:
    dfs(np.ones(N, dtype=bool), [], 0)
    print('%s: exhaustive over this generator class: max length %d' % (name, best[0]))
except StopIteration:
    print('%s: zonotope of length %d found: %s' % (name, best[0], best[1]))
