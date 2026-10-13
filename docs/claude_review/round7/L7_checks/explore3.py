import numpy as np
import cub13f as F
from collections import Counter
rng = np.random.default_rng(1)
tau, s34, s35 = -2.0+0j, 3.0+0j, 2.5+0j
sols = []
for k in range(300):
    v0 = rng.normal(size=2) * 3 + 1j * rng.normal(size=2) * 3
    v = F.newton(v0, tau, s34, s35)
    if v is None: continue
    if all(np.linalg.norm(v - w) > 1e-6 for w in sols):
        sols.append(v)
print(len(sols), 'distinct (s45, sigma) solutions')
for v in sols:
    d = F.build(tau, s34, s35, v[0])
    ev = F.events(d['c'])
    cnt = Counter()
    for p, big in ev:
        for g in big: cnt[g] += 1
    trip = {g for g in cnt if len(g) == 3}
    pairs = {g: cnt[g] for g in cnt if len(g) == 2}
    big4 = {g for g in cnt if len(g) >= 4}
    ok = (trip == F.TRIPLES13 and all(cnt[g] == 1 for g in trip) and pairs == {(1,4):1,(1,5):1,(2,3):1,(1,6):3} and not big4)
    print('s45=%s sigma=%s  nevents=%d  ok=%s' % (np.round(v[0], 5), np.round(v[1], 5), len(ev), ok))
    if not ok: print('   ', dict(cnt))
