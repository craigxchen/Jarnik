"""Theorem 2.1 check: transplant vectors (beta_j over the preceding split prime) are density points of G(p).
Full torus T^r, primes 13..73 (no 5).  Evidence for the proved lemma; exact check of gap-admissibility up to T."""
import sys, math, json
import numpy as np
sys.path.insert(0, '.')
from common import *

allsp = split_primes(200)
ps = [13, 17, 29, 37, 41, 53, 61, 73]
qs = [allsp[allsp.index(p) - 1] for p in ps]
xb = np.array([math.atan2(*gaussian_prime(q)[::-1]) / (math.pi / 2) for q in qs])  # arg(a+bi)/(pi/2)
logs = [math.log(p) for p in ps]
T = 24.0
V, W = monomials(logs, T)
a0, a1 = halfwidths(W)
res = {'primes': ps, 'q': qs, 'T': T, 'n_monomials': int(len(V))}
t = V @ xb
m0 = dist_Z(t) - a0
m1 = dist_Z(t - 0.5) - a1
res['transplant_min_margin_0'] = float(m0.min()); res['transplant_min_margin_half'] = float(m1.min())
rng = np.random.default_rng(1)
out = []
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
    n = 4000
    X = xb + rng.uniform(-eps, eps, size=(n, len(ps)))
    ok = np.ones(n, bool)
    for i0 in range(0, len(V), 5000):
        Vb = V[i0:i0+5000]; tt = X @ Vb.T
        ok &= np.all((dist_Z(tt) >= a0[i0:i0+5000]) & (dist_Z(tt - 0.5) >= a1[i0:i0+5000]), axis=1)
    out.append((eps, float(ok.mean())))
    print(eps, ok.mean(), flush=True)
res['cube_survival'] = out
# uniform points on the torus for comparison
X = rng.uniform(0, 1, size=(4000, len(ps)))
ok = np.ones(len(X), bool)
for i0 in range(0, len(V), 5000):
    Vb = V[i0:i0+5000]; tt = X @ Vb.T
    ok &= np.all((dist_Z(tt) >= a0[i0:i0+5000]) & (dist_Z(tt - 0.5) >= a1[i0:i0+5000]), axis=1)
res['uniform_survival'] = float(ok.mean())
print(json.dumps(res))
json.dump(res, open('density_check_out.json', 'w'), indent=1)
