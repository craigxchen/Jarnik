"""Formal g for q = 243, k = 3 (n = 243^3), m = 7, 8: pattern counting + all signs; a few form subsets."""
import numpy as np, itertools, random
from f9lib import Fq
from subsets import forms_list
F = Fq(5); q = 243
addt = np.array(F.addt, dtype=np.int32); neg = np.array(F.neg, dtype=np.int32); chi = np.array(F.chi, dtype=np.int8)
x0, x1, x2 = [a.ravel() for a in np.meshgrid(np.arange(q, dtype=np.int32), np.arange(q, dtype=np.int32), np.arange(q, dtype=np.int32), indexing='ij')]
forms = forms_list(3)
C = []
for lam in forms:
    s = np.zeros(len(x0), np.int32)
    for c, col in zip(lam, (x0, x1, x2)):
        if c == 1: s = addt[s, col]
        elif c == 2: s = addt[s, neg[col]]
    C.append(chi[s])
del x0, x1, x2
n = q**3; h = q**2
def best_g(S):
    m = len(S)
    code = np.zeros(n, np.int64)
    for j in S: code = code*3 + (C[j].astype(np.int64) + 1)
    cnt = np.bincount(code, minlength=3**m)
    pats = np.array(list(itertools.product((-1, 0, 1), repeat=m)))
    best = (-1, None)
    for eps in itertools.product((1, -1), repeat=m-1):
        e = np.array((1,) + eps)
        g = int((np.min(1 + e*pats, axis=1) * cnt).sum())
        if g > best[0]: best = (g, (1,) + eps)
    return best
random.seed(1)
for m in (7, 8):
    cands = [tuple(range(m))] + [tuple(sorted(random.sample(range(13), m))) for _ in range(12)]
    res = max((best_g(S), S) for S in cands)
    print('q=243 k=3 m=%d: n=%d h=%d best g found=%d (g/h=%.2f), forms %s, signs %s' % (m, n, h, res[0][0], res[0][0]/h, [forms[j] for j in res[1]], res[0][1]), flush=True)
