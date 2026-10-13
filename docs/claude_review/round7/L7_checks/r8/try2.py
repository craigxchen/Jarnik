import numpy as np, itertools, sys
from fwd import *
seed = int(sys.argv[1]) if len(sys.argv)>1 else 0
rng = np.random.default_rng(seed)
E = eval(sys.argv[2]) if len(sys.argv)>2 else None
if E is None:
    from try1 import random_type
E = [tuple(e) for e in E]
stats = {}
for trial in range(400):
    R = setup(E, rng)
    V0 = (rng.normal(size=2*len(E)+28) + 1j*rng.normal(size=2*len(E)+28))
    V, ok = gn(V0, E, R)
    if not ok: stats['fail'] = stats.get('fail',0)+1; continue
    res, flag = check(V, E)
    stats[flag] = stats.get(flag,0)+1
    if flag == 'ok':
        g, times = res
        pair = sum(1 for t,S in times if len(S)==2); trip = sum(1 for t,S in times if len(S)==3)
        print('FOUND trial', trial, 'pair', pair, 'trip', trip); np.save('member_%d.npy'%seed, V); break
print(stats)
