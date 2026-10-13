import numpy as np, itertools, sys, pickle
from graph import *
rng = np.random.default_rng(int(sys.argv[1]))
E = list(itertools.combinations(range(7), 3))
r = grow([3]*7, E, rng, verbose=True)
cur, V, P, full = r
print('reached', len(cur))
remaining = [e for e in E if e not in cur]
for e in remaining[:4]:
    Pn = GraphProblem([3]*7, cur + [e], rng); Pn.R = P.R
    flags = {}
    for (a, b) in itertools.combinations(e, 2):
        pp = P.pairpoly(V, a, b)
        for t in np.roots(pp[::-1]):
            if any(abs(t - V[i]) < 1e-6 for i in range(P.m)): continue
            z0 = P.val(V, e[0], t)
            V0 = np.concatenate([np.append(V[:P.m], t), np.append(V[P.m:2*P.m], z0), V[2*P.m:]])
            Vn, ok = Pn.gn(V0)
            fl = Pn.contacts(Vn)[1] if ok else 'noconv'
            flags[fl] = flags.get(fl, 0) + 1
    # rank of Jacobian at current point: dimension check
    print('event', e, flags)
Jm = P.J(V); sv = np.linalg.svd(Jm, compute_uv=False)
print('J shape', Jm.shape, 'rank', np.sum(sv > 1e-8*sv[0]), 'expected dim (with gauge 6):', P.N - np.sum(sv > 1e-8*sv[0]))
pickle.dump((cur, V, P.R), open('beta_partial_%s.pkl' % sys.argv[1], 'wb'))
