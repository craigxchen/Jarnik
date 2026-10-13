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
    for S, tn in [((1,3,4),'t134'), ((1,3,5),'t135'), ((2,3,4),'t234')]:
        t = d[tn]
        print(tn, np.round(t, 5), [np.round(P.polyval(t, cs[i-1]), 6) for i in range(1, 6)])
    print('params', {k: np.round(v, 4) for k, v in d.items()})
    break
