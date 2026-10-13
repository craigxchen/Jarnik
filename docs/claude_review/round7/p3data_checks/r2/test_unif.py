"""SA calibration on targets known to be feasible: k stacked copies of the Paley Hadamard design
(M - 1 balanced columns, pairwise distance M/2 each) give r = k(M-1), d_xy = kM/2 = r/2 + k/2 exactly.
Target: d_xy >= r/2 + k/2 for all pairs (Plotkin-tight).  Reports best violation V over 3 SA runs."""
import sys, numpy as np
sys.path.insert(0, '..')
from eqw import *
from gauss7 import paley_normalised
rng = np.random.default_rng(3)
for M in (8, 12, 20, 24):
    P = paley_normalised(M - 1)[:, 1:]
    for k in (2, 3):
        r = k * (M - 1)
        T = np.full((M, M), k * M // 2, dtype=np.int32); np.fill_diagonal(T, 0)
        S = np.hstack([P] * k).astype(np.int8)
        assert exact_ok(S, T)                       # the stacked design meets the target exactly
        vs = [sa_run(T, random_balanced(M, r, rng), 20_000_000, seed=j + 1)[0] for j in range(3)]
        print('M=%d k=%d r=%d target d >= %d: Paley stack feasible; SA best V over 3 runs %s' % (
            M, k, r, k * M // 2, vs), flush=True)
