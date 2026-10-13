import numpy as np
from q4b import *
rng = np.random.default_rng(11)
for trial in range(150):
    u0 = rng.normal(size=16) * 1.5 + 1j * rng.normal(size=16) * 1.5
    fixed = list(rng.choice(16, size=4, replace=False))
    u = newton(u0, fixed)
    if u is None: continue
    cs, d = coords(u)
    vals = {k: d[k] for k in NAMES if k[0] in 'rt'}; vals['r34'] = 0; vals['r35'] = 1
    coll = [(a, b) for a, b in itertools.combinations(vals, 2) if abs(vals[a] - vals[b]) < 1e-6]
    if coll: continue
    ok, cnt, bad, md = genuine(u)
    diff = {tuple(sorted(k)): v for k, v in cnt.items() if EXPECT.get(k, 0) != v}
    miss = {tuple(sorted(k)): v for k, v in EXPECT.items() if cnt.get(k, 0) != v}
    print('ok', ok, 'md', md, 'bad', [(np.round(r, 4), [tuple(sorted(x)) for x in b]) for r, b in bad][:4])
    print('   extra', diff, '\n   missing', miss)
    print('   lam1, lam2 =', d['lam1'], d['lam2'])
