# Global load accounting on the forced primes F (level 1 only, a lower bound for the true demand):
#   forced mass  = sum_q  min_c #{v in V_K/+-: <v,c>=0 mod n_q} * log q      (exact minimum over parity vectors)
#   genuine mass = sum_q #{genuine collisions at q (level 1)} * log q
#   capacity     = sum_v sum_s B(v,s)   (each collision of v uses log q of ONE unit's budget)
#   genuine total demand (all levels, on F) = sum_(v,s) log(F-smooth part of D_s(Z_v))
import sys, math, itertools
import numpy as np
from lib import *
from forced_coh import load_set

def vq(x, q):
    x = abs(x); k = 0
    while x and x % q == 0: x //= q; k += 1
    return k

for K, r, tag in [(2, 4, 'first'), (2, 4, 'last'), (3, 4, 'first'), (3, 4, 'last'), (4, 4, 'first'), (4, 4, 'last')]:
    P, F = load_set(K, r, tag); V = monos(r, K); Va = np.array(V, dtype=np.int64)
    pis = [gauss_prime(p) for p in P]
    forced = 0.0; gen = 0.0
    for q in F:
        m = mq(q); n = m // 4; par = nu(P, q)
        ranges = [range(par[j], n, 2) if n % 2 == 0 else range(n) for j in range(r)]
        best = None
        it = itertools.product(*ranges)
        while True:
            chunk = list(itertools.islice(it, 300000))
            if not chunk: break
            C = np.array(chunk, dtype=np.int64)
            cnt = ((C @ Va.T) % n == 0).sum(axis=1).min()
            best = cnt if best is None else min(best, cnt)
        forced += best * math.log(q)
        T = Torus(q); L = np.array([T.L(pi) % n for pi in pis])
        gen += int(((Va @ L) % n == 0).sum()) * math.log(q)
    cap = sum(budget(P, v, s) for v in V for s in range(4))
    capmax = sum(max(budget(P, v, s) for s in range(4)) for v in V)
    gdem = sum(sum(vq(Ds(Zv(pis, v), s), q) * math.log(q) for q in F) for v in V for s in range(4))
    print(f"K={K} r={r} {tag} P={P} |F|={len(F)} #V/2={len(V)}: forced mass (level 1) {forced:.1f}; genuine level-1 mass {gen:.1f}; "
          f"genuine all-level demand on F {gdem:.1f}; capacity sum_(v,s) B {cap:.1f}; ratios forced/cap {forced/cap:.3f}, gen-demand/cap {gdem/cap:.3f}", flush=True)
