import numpy as np, itertools, sys, pickle
from graph import *
rng = np.random.default_rng(int(sys.argv[1]))
E = list(itertools.combinations(range(7), 3))
r = grow2([3]*7, E, rng, verbose=True)
cur, V, P, full = r
print('reached', len(cur), 'full', full)
pickle.dump((cur, V, P.R, full), open('beta_g2_%s.pkl' % sys.argv[1], 'wb'))
