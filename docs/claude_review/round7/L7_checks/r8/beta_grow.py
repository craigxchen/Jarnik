import numpy as np, itertools, sys, pickle
from graph import *
rng = np.random.default_rng(int(sys.argv[1]))
E = list(itertools.combinations(range(7), 3))
best = 0
for att in range(int(sys.argv[2])):
    r = grow([3]*7, E, rng, verbose=False)
    if r is None: continue
    cur, V, P, full = r
    best = max(best, len(cur))
    print('attempt', att, 'reached', len(cur), 'events', flush=True)
    if full:
        pickle.dump((cur, V, P.R), open('beta_member_%s.pkl' % sys.argv[1], 'wb')); print('FULL'); break
print('best', best)
