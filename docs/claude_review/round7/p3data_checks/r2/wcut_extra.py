"""Certified W_cut (column generation with exact max-cut pricing, dual bound) for instances not in
../perron_check.json: M = 24, X = 96, t = 0, 1.  Output: wcut_extra.json."""
import sys, os, json
import numpy as np
sys.path.insert(0, '..')
from qsolve import lp_colgen
out = []
for M, t, Xf in ((24, 0, 4), (24, 1, 4)):
    D = np.load('../data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    W, x, masks, nu, lb = lp_colgen(D)
    lam2 = 2 * float(np.linalg.eigvalsh(D)[-1])
    out.append(dict(M=M, t=t, X=M * Xf, W=W, lb=lb, lam2=lam2))
    print(M, t, M * Xf, W, lb, lam2, W / lam2, flush=True)
json.dump(out, open('wcut_extra.json', 'w'), indent=1)
