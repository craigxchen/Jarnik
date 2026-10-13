# Exact forced/free status for K = 4, 5 at r = 4 (first, last sets), q <= 1200.
import time
from lib import *
from free import exists_free, verify
allp = split_primes(5, 200)
for K in (4, 5):
    for tag, P in (('first', allp[:4]), ('last', allp[-4:])):
        t0 = time.time(); nonfree = []
        for q in primerange(3, 1201):
            if q in P: continue
            m = mq(q); n = m // 4; need = (n % 2 == 0); par = nu(P, q)
            res, nodes = exists_free(n, par, K, need)
            if res is None: nonfree.append(q)
            else: assert verify(res, n, K, par, need)
        print(f"K={K} r=4 {tag} P={P}: forced primes q<=1200: {nonfree}; max={max(nonfree)} "
              f"(|S_2|={2*16+8+1}; pigeonhole forced if n<{2*16+8+1}) [{time.time()-t0:.0f}s]", flush=True)
