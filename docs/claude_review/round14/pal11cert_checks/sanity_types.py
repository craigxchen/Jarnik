"""Sanity check of the exact (sigma,t) tables: random balanced integer vectors (entries in [-3,3], all supports) and
all vectors of the pools bal:12:1 must have t in the attainable set of typeflats (sigma <= 8) and t >= f(sigma)."""
import numpy as np, itertools
from fplib import pal_patterns
M_, pats = pal_patterns(11)
core = [P for P in pats if len(P) == 6]
Cv = np.zeros((11, 11), dtype=np.int64)
for c, P in enumerate(core):
    for x in P: Cv[c, x - 1] = 1
att = {2: [6], 3: [9], 4: [6, 8, 9, 10, 11], 5: [7, 8, 9, 10, 11], 6: [2, 4, 6, 7, 8, 9, 10, 11],
       7: [5, 6, 7, 8, 9, 10, 11], 8: [4, 5, 6, 7, 8, 9, 10, 11]}
f = {2: 6, 3: 9, 4: 6, 5: 7, 6: 2, 7: 5, 8: 4, 9: 3, 10: 4, 11: 4, 12: 1}
rng = np.random.default_rng(5)
seen = set(); bad = 0; n = 0
for trial in range(300000):
    s = rng.integers(2, 13)
    R = rng.choice(12, s, replace=False)
    z = np.zeros(12, dtype=np.int64)
    z[R] = rng.integers(-3, 4, s)
    z[R[-1]] -= z.sum()
    if (z != 0).sum() < 2: continue
    sg = int((z != 0).sum()); t = int(((Cv @ z[1:]) != 0).sum())
    n += 1; seen.add((sg, t))
    if t < f[sg] or (sg <= 8 and t not in att[sg]): bad += 1
print('random vectors tested', n, 'violations', bad)
print('observed (sigma,t) with t<=6:', sorted(x for x in seen if x[1] <= 6))
