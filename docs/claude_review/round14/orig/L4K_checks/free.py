# Exact: for which primes q does a parity-respecting label vector with NO collision of any
# v (1<=||v||_1<=K) with any unit exist?  (Level 1; collision-free at q implies at every q^a.)
# Work mod n = m/4 with c = ell mod n. q in Q_1 (n even): need c_j = nu_j mod 2.
# q in Q_3 (n odd): parity free (lift ell = c + n t).
# Exhaustive DFS; symmetry: scale by units of Z/n (parity-preserving when n even since units odd).
import sys, itertools, math, time
from lib import *

def sums_update(Sig, x, n, K):
    # Sig[t] = set of <v,c> mod n over ||v||_1 <= t (prefix); add column value x
    new = []
    for t in range(K + 1):
        S = set()
        for k in range(-t, t + 1):
            kx = (k * x) % n
            for s in Sig[t - abs(k)]:
                S.add((s + kx) % n)
        new.append(S)
    return new

def forbidden(Sig, n, K):
    F = set()
    for k in range(1, K + 1):
        T = Sig[K - k]
        for x in range(n):
            if (k * x) % n in T:
                F.add(x)
    return F

def exists_free(n, par, K, need_par, count_nodes=None):
    r = len(par)
    # column 1: representatives of unit orbits: divisors d of n (d in 1..n-1, d|n) -- orbit of x is determined by gcd(x,n)
    reps = [d for d in range(1, n) if n % d == 0] if n > 1 else []
    Sig0 = [{0} for _ in range(K + 1)]
    nodes = [0]
    def dfs(j, Sig, cur):
        nodes[0] += 1
        if j == r:
            return list(cur)
        F = forbidden(Sig, n, K)
        cand = reps if j == 0 else range(n)
        for x in cand:
            if x in F: continue
            if need_par and x % 2 != par[j]: continue
            if j == 0 and need_par and x % 2 != par[0]:
                continue
            res = dfs(j + 1, sums_update(Sig, x, n, K), cur + [x])
            if res is not None:
                return res
        return None
    # caution: when need_par (n even) units are odd so gcd parity class preserved; orbit reps of the right
    # parity: x with gcd(x,n)=d has parity of d when n even. fine.
    res = dfs(0, Sig0, [])
    return res, nodes[0]

def verify(c, n, K, par, need_par):
    r = len(c)
    if need_par and any(c[j] % 2 != par[j] for j in range(r)): return False
    for v in monos(r, K):
        if sum(a * b for a, b in zip(v, c)) % n == 0: return False
    return True

if __name__ == '__main__':
    import random
    rng = random.Random(7)
    allp = split_primes(5, 200)
    sets = []
    for r in range(4, 9):
        sets.append(('first', allp[:r]))
        sets.append(('last', allp[-r:]))
        sets.append(('rand', sorted(rng.sample(allp, r))))
    for K in (2, 3):
        for tag, P in sets:
            r = len(P)
            QK = 8 * K * (2 * r + 1) ** K + 1
            hi = min(QK, 3000 if K == 3 else QK)
            nonfree = []
            t0 = time.time(); maxnodes = 0
            for q in primerange(3, hi + 1):
                if q in P: continue
                m = mq(q); n = m // 4
                need_par = (n % 2 == 0)
                par = nu(P, q)
                res, nodes = exists_free(n, par, K, need_par)
                maxnodes = max(maxnodes, nodes)
                if res is None:
                    nonfree.append(q)
                else:
                    assert verify(res, n, K, par, need_par)
            # pigeonhole bound: |S_{floor K/2}| > n forces collision
            h = K // 2
            S = 1 + 2 * sum(math.comb(r, i) * math.comb(h, i) * 2 ** (i - 1) for i in range(1, h + 1))  # not used
            print(f"K={K} r={r} {tag} P={P} Q_K={QK} checked q<={hi}: non-free q = {nonfree}; "
                  f"max non-free = {max(nonfree) if nonfree else None}; nodes<= {maxnodes}; {time.time()-t0:.1f}s", flush=True)
