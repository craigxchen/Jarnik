import numpy as np, itertools, sys
from mob import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def random_type(m=12):
    while True:
        T = list(itertools.combinations(range(7),3)); rng.shuffle(T)
        E = []; cnt = {}
        for t in T:
            if all(cnt.get(pr,0) < 2 for pr in itertools.combinations(t,2)):
                E.append(t)
                for pr in itertools.combinations(t,2): cnt[pr] = cnt.get(pr,0)+1
            if len(E) == m: return E
E = random_type()
print('type', E)
y = rng.normal(size=7) + 1j*rng.normal(size=7)
sols = []; stats = {}
for trial in range(3000):
    X0 = (rng.normal(size=24) + 1j*rng.normal(size=24))*rng.choice([0.5,1,2,4])
    X, ok = newton(X0, y, E)
    if not ok: stats['fail'] = stats.get('fail',0)+1; continue
    times, flag = contacts(X, y, E)
    stats[flag] = stats.get(flag,0)+1
    if flag != 'ok': continue
    if not any(np.linalg.norm(X - S) < 1e-6*(1+np.linalg.norm(S)) for S in sols):
        sols.append(X); print('new sol', len(sols), 'trial', trial, classvec(times), flush=True)
print(stats, 'distinct ok sols:', len(sols))
