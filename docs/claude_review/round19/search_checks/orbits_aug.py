# S_7-orbit representatives of D-edge multisets of K_7 by augmentation: reps(D) = {canon(r + e)}.
# Complete: removing any edge of M gives sigma(r) for a rep r, so sigma^{-1}(M) = r + sigma^{-1}(e0).
# Writes reps_D<D>.json in the same canonical form as common.canon (so it agrees with orbit_reps for D <= 5).
import sys, json, os, time
from math import comb
import common
W = os.path.dirname(os.path.abspath(__file__))
Dmax = int(sys.argv[1])
reps = [()]
for D in range(1, Dmax + 1):
    t0 = time.time()
    fn = os.path.join(W, 'reps_D%d.json' % D)
    if os.path.exists(fn):
        new = [tuple(tuple(e) for e in r) for r in json.load(open(fn))]
    else:
        new = sorted(set(common.canon(tuple(sorted(r + (e,)))) for r in reps for e in common.edges))
        json.dump(new, open(fn, 'w'))
    reps = new
    pairs = sum(comb(21 - len(set(A)) + D - 1, D) for A in reps)
    print('D', D, 'orbits', len(reps), 'pairs', pairs, 'time', round(time.time() - t0, 1)); sys.stdout.flush()
