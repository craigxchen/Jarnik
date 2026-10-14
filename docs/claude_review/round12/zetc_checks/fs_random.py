# random search for T=1 models with rho >= max(2 h_max, zeta_ub + 1)  (zeta_ub: LP with coords+columns+pair functionals)
import random, sys, time
from faithful_search import *
d = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
random.seed(seed)
mons = cube_monomials(d)
best = []
t0 = time.time()
for tr in range(trials):
    r = random.randint(max(2, d//2), 2*d)
    cols = []
    for _ in range(r):
        sz = random.randint(2, d)
        cols.append(sorted(random.sample(range(d), sz)))
    Ls = [random.choice([1, 1, 2]) for _ in cols]
    res = analyze(d, cols, Ls, mons=mons)
    margin = res['rho'] - res['mlo']
    if margin >= -1 or res['zeta_ub'] < res['rho'] - 1e-6:
        print('trial', tr, 'margin', margin, 'cols', cols, 'Ls', Ls, res, flush=True)
print('done %d trials in %.0fs' % (trials, time.time()-t0))
