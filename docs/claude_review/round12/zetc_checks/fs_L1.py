# targeted search: all slabs L=1 (width 2, so LP-supercritical needs m >= 8), want zeta_LP < rho_p and rho_p >= 8
import random, sys, time
from faithful_search import analyze
from tclib import cube_monomials
d = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]); smax = int(sys.argv[4])
random.seed(seed); mons = cube_monomials(d); t0 = time.time()
for tr in range(trials):
    r = random.randint(d - 2, 2 * d)
    cols = set()
    while len(cols) < r:
        sz = random.randint(2, smax)
        cols.add(tuple(sorted(random.sample(range(d), sz))))
    cols = [list(c) for c in cols]
    # twin-free check (points incl. base 'd' with patterns)
    pats = [tuple(1 if x in c else 0 for c in cols) for x in range(d)] + [tuple(0 for c in cols)]
    if len(set(pats)) < d + 1:
        continue
    res = analyze(d, cols, [1] * len(cols), mons=mons)
    if res['rho'] >= 8 or res['zeta_ub'] < res['rho'] - 1e-6:
        print('trial', tr, 'r=%d' % r, 'cols', cols, {k: res[k] for k in ('N', 'rho', 'zeta_ub', 'mlo', 'ok')}, flush=True)
print('done %.0fs' % (time.time() - t0))
