"""Section 7 evidence: base coset H_0 = R b through 0 carrying a non-character slack seed
(b = (1,1,1,1,1,1,1,1), seed theta_0 = arg(2+i): phi_i = arg (2+i)^{b_i}, norms 5 < p_i), and its translates
C_eta = eta + H_0 (eta transverse, random direction).  Characters (n = <b,v> = 0) are IGNORED here (they are the
profile's margin problem); we test only whether NON-character windows cover the line.  Primes 13..73, T = 22."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *
from line_cosets import line_survivors
ps = [13, 17, 29, 37, 41, 53, 61, 73]
V, W = monomials([math.log(p) for p in ps], 22.0)
b = np.ones(8, dtype=np.int64)
rng = np.random.default_rng(5)
res = []
for trial in range(3):
    u = rng.normal(size=8); u -= u.mean(); u /= np.abs(u).max()      # transverse to b, sup-norm 1
    for eta in (0.0, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        r = line_survivors(V, W, b, eta * u)
        res.append({'trial': trial, 'eta': eta, 'surv_len': r['surv_len'], 'n_gaps': r['n_gaps']})
        print(res[-1], flush=True)
json.dump(res, open('seeded_base_out.json', 'w'), indent=1)
