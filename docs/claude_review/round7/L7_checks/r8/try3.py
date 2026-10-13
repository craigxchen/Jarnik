import numpy as np, itertools, sys
from fwd import *
rng = np.random.default_rng(int(sys.argv[1]))
def random_type(m):
    while True:
        T = list(itertools.combinations(range(7),3)); rng.shuffle(T)
        E = []; cnt = {}
        for t in T:
            if all(cnt.get(pr,0) < 2 for pr in itertools.combinations(t,2)):
                E.append(t)
                for pr in itertools.combinations(t,2): cnt[pr] = cnt.get(pr,0)+1
            if len(E) == m: return E
for m in [1, 4, 6, 8, 9, 10, 11, 12]:
    E = random_type(m); stats = {}
    for trial in range(200):
        R = setup(E, rng)
        V0 = (rng.normal(size=2*len(E)+28) + 1j*rng.normal(size=2*len(E)+28))
        V, ok = gn(V0, E, R)
        if not ok: stats['fail'] = stats.get('fail',0)+1; continue
        res, flag = check(V, E)
        stats[flag] = stats.get(flag,0)+1
    print(m, E, stats, flush=True)
