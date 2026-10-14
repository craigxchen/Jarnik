"""Remark 2 check: the Theorem-2 vector (2.1) at q > q_T (random sigma_q) has, for norm 2h+2 > K, monomials v with
<v,a> = 0 as an INTEGER at a positive fraction of moduli; such v collide at every level q^a, so their demand
is unbounded as X grows (the system violates norm-(2h+2) budgets).  K = 2 (h = 1, beta_j = j), r = 4 and 8."""
import random, math
import numpy as np
from rlib import *
for r in (4, 8):
    K = 2; h = 1; B = r
    Lam = K * (K + 1) * B + K
    AK = K * (2 * (K + 1) ** 2 * B + 2 * K + 1); qT = 4 * AK + 5
    ps = split_primes_upto(200)[:r]
    V4 = np.array([v for v in VK(r, 4) if sum(abs(x) for x in v) == 4])
    V2 = np.array(VK(r, 2))
    rng = random.Random(11)
    qs = [q for q in primes_upto(20000) if q > qT and all(q != p for p in ps)]
    hit = 0; ex = None
    for q in qs:
        sig = list(range(r)); rng.shuffle(sig)
        nv = nu(ps, q)
        a = np.array([2 * (Lam + (K + 1) * (sig[j] + 1)) + nv[j] for j in range(r)])
        assert np.all(V2 @ a != 0)
        z = (V4 @ a == 0)
        if z.any():
            hit += 1
            if ex is None: ex = (q, tuple(V4[np.argmax(z)]))
    print(f"r={r}: q_T={qT}; primes q in (q_T, 20000]: {len(qs)}; with some norm-4 v having <v,a>=0 in Z: {hit} ({hit/len(qs):.3f}); example {ex}")
