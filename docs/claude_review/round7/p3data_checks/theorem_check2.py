"""Adapted vs uniform vs Paley profiles at EQUAL total weight W in the uniform residue model U_X
(Monte Carlo, 2000 draws each), with the Theorem 5.1 bound.  Adapted profiles: the realised prime
profiles (no separation) of q_run, weights scaled to the common W.  Evidence.
Usage: python3 theorem_check2.py M Xf [max_profiles]"""
import sys, json, math
import numpy as np
from theorem_check import fbar_tables, fbar_lo, check_profile
from gauss7 import paley_normalised, profile_construct
from qcore import mask_to_s

M = int(sys.argv[1]); Xf = int(sys.argv[2]); kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
X = Xf * M
rng = np.random.default_rng(11)
tabs = fbar_tables(X)
recs = [json.loads(l) for l in open('q_results_M%d.jsonl' % M)]
recs = [r for r in recs if r['X'] == X and 'prime_cols' in r][:kmax]
out = []
for rec in recs:
    cols = rec['prime_cols']
    S_a = np.array([mask_to_s(m, M) for m, pl in cols for p in pl]).T
    w_a = np.array([math.log(p) for m, pl in cols for p in pl])
    W0 = w_a.sum()
    grid = [W0 * f for f in (1, 2, 3, 4)]
    print('M=%d X=%d adapted profile of trial %d (r=%d, W0=%.1f):' % (M, X, rec['t'], len(w_a), W0))
    out += check_profile('adapted', S_a, w_a, X, tabs, 2000, rng, [g / W0 for g in grid])
    Hm = paley_normalised(M - 1)
    out += check_profile('uniform', Hm[:, 1:], np.ones(M - 1), X, tabs, 2000, rng, [g / (M - 1) for g in grid])
    S_c = profile_construct(M - 1, 8)
    out += check_profile('paley-b8', S_c, np.ones(S_c.shape[1]), X, tabs, 2000, rng, [g / S_c.shape[1] for g in grid])
json.dump(out, open('theorem_check2_M%d_X%d.json' % (M, X), 'w'))
