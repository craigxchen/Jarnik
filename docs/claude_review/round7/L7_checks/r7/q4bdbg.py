import numpy as np
from q4b import *
from collections import Counter
rng = np.random.default_rng(11)
pat = Counter(); shown = 0
for trial in range(150):
    u0 = rng.normal(size=16) * 1.5 + 1j * rng.normal(size=16) * 1.5
    fixed = list(rng.choice(16, size=4, replace=False))
    u = newton(u0, fixed)
    if u is None: continue
    cs, d = coords(u)
    # collisions among parameters
    vals = {k: d[k] for k in NAMES if k[0] in 'rt'}
    vals['r34'] = 0; vals['r35'] = 1
    coll = [(a, b) for a, b in itertools.combinations(vals, 2) if abs(vals[a] - vals[b]) < 1e-6]
    small = [k for k in ('lam1', 'lam2') if abs(d[k]) < 1e-6 or abs(d[k] - 1) < 1e-6]
    key = tuple(sorted(coll)) + tuple(small)
    pat[key] += 1
for k, v in pat.most_common(12): print(v, k)
