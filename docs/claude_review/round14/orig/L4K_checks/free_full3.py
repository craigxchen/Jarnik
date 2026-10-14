# Complete the K=3 scan of free.py from q=3001 up to Q_3 (first-found DFS; any non-free q is reported).
import time
from lib import *
from free import exists_free, verify
import random
rng = random.Random(7)
allp = split_primes(5, 200)
sets = []
for r in range(4, 9):
    sets.append(('first', allp[:r])); sets.append(('last', allp[-r:])); sets.append(('rand', sorted(rng.sample(allp, r))))
K = 3
for tag, P in sets:
    r = len(P); QK = 8 * K * (2 * r + 1) ** K + 1
    nonfree = []; t0 = time.time()
    for q in primerange(3001, QK + 1):
        m = mq(q); n = m // 4; need = (n % 2 == 0); par = nu(P, q)
        res, nodes = exists_free(n, par, K, need)
        if res is None: nonfree.append(q)
        else: assert verify(res, n, K, par, need)
    print(f"K=3 r={r} {tag}: q in (3000, {QK}]: non-free = {nonfree} ({time.time()-t0:.0f}s)", flush=True)
