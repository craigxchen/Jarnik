# clean random survey: every model logged; summary of zonotope-freeness (zeta_LP < rho_p) and LP-supercriticality
import random, sys, time
from faithful_search import analyze
from tclib import cube_monomials
d = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3])
random.seed(seed); mons = cube_monomials(d)
stats = dict(models=0, twinfree=0, zfree=0, zfree_twinfree=0, supercrit_and_zfree=0)
for tr in range(trials):
    r = random.randint(2, 2*d)
    cols = set()
    while len(cols) < r:
        cols.add(tuple(sorted(random.sample(range(d), random.randint(2, d)))))
    cols = [list(c) for c in cols]
    Ls = [random.choice([1, 1, 2]) for _ in cols]
    res = analyze(d, cols, Ls, mons=mons)
    pats = [tuple(1 if x in c else 0 for c in cols) for x in range(d)] + [tuple(0 for c in cols)]
    tf = len(set(pats)) == d + 1
    zf = res['zeta_ub'] < res['rho'] - 1e-9
    stats['models'] += 1; stats['twinfree'] += tf; stats['zfree'] += zf; stats['zfree_twinfree'] += (zf and tf)
    stats['supercrit_and_zfree'] += res['ok']
    print('r=%d twinfree=%d rho_p=%d zeta_LP=%.3f 2hmax=%d ok=%d cols=%s Ls=%s' % (r, tf, res['rho'], res['zeta_ub'], 2*2*max(res['widths']), res['ok'], cols, Ls), flush=True)
print('SUMMARY d=%d' % d, stats)
