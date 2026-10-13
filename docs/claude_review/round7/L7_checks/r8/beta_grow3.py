# Resumable class-beta growth (7 cubic graphs, all 35 triple events), saves state after each addition.
import numpy as np, itertools, sys, pickle, os
from graph import *
seed = int(sys.argv[1]); rng = np.random.default_rng(seed)
E = list(itertools.combinations(range(7), 3)); degs = [3]*7
st = 'beta_state_%d.pkl' % seed
if os.path.exists(st):
    cur, V, R = pickle.load(open(st, 'rb')); P = GraphProblem(degs, cur, rng); P.R = R
else:
    cur = []; P = GraphProblem(degs, cur, rng)
    for _ in range(50):
        V, ok = P.gn(rng.normal(size=P.N) + 1j*rng.normal(size=P.N))
        if ok and P.contacts(V)[1] == 'ok': break
rem = [e for e in E if e not in cur]
fails = {}
while rem:
    order = list(rem); rng.shuffle(order); added = False
    for e in order:
        r = try_add(P, V, cur, e, degs, rng, tries=60)
        if r is not None:
            P, V = r; cur = cur + [e]; rem.remove(e); added = True
            pickle.dump((cur, V, P.R), open(st, 'wb'))
            print('added', len(cur), e, flush=True); break
    if not added:
        print('STUCK', len(cur), rem, flush=True); break
if not rem:
    print('FULL', flush=True)
    Jm = P.J(V); sv = np.linalg.svd(Jm, compute_uv=False)
    rank = int(np.sum(sv > 1e-9*sv[0]))
    print('unknowns', P.N, 'rank', rank, 'solution-set dim (incl 6-dim gauge + 7 scalings fixed)', P.N - rank, flush=True)
    print('contacts', P.contacts(V)[1], 'residual', np.linalg.norm(P.F(V)), flush=True)
